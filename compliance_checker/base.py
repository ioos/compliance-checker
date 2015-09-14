#!/usr/bin/env python

"""
Compliance Checker
"""

from functools import wraps
import collections
import pprint

from wicken.netcdf_dogma import NetCDFDogma
from wicken.xml_dogma import MultipleXmlDogma
from wicken.exceptions import DogmaGetterSetterException
from netCDF4 import Dataset
from owslib.swe.observation.sos100 import SensorObservationService_1_0_0
from owslib.swe.sensor.sml import SensorML
from owslib.namespaces import Namespaces
from petulantbear.netcdf_etree import namespaces as pb_namespaces

def get_namespaces():
    n = Namespaces()
    ns = n.get_namespaces(["ogc","sml","gml","sos","swe","xlink"])
    ns["ows"] = n.get_namespace("ows110")
    return ns

class DSPair(object):
    """
    Structure to hold a dataset/SOS instance and dogma pairing.

    Passed to each check method.
    """
    dataset = None
    dogma   = None
    def __init__(self, ds, dogma):
        self.dataset = ds
        self.dogma = dogma

class BaseCheck(object):
    HIGH   = 3
    MEDIUM = 2
    LOW    = 1

    supported_ds = []

    @classmethod
    def beliefs(cls):
        raise NotImplementedError("Define this in your derived Checker class")

    def setup(self, ds):
        """
        Common setup method for a Checker.

        Automatically run when running a CheckSuite. Define this method in your Checker class.
        """
        pass

    def load_datapair(self, ds):
        """
        Returns a DSPair object with the passed ds as one side and the proper Dogma object on the other.

        Override this in your derived class.
        """
        raise NotImplementedError("Define this in your derived checker class")

class BaseNCCheck(object):
    """
    Base Class for NetCDF Dataset supporting Check Suites.
    """
    supported_ds = [Dataset]

    def load_datapair(self, ds):
        # allow ncml as well as nc prefixes
        namespaces         = pb_namespaces.copy()
        namespaces['nc']   = namespaces['ncml']
        namespaces['ncml'] = namespaces['ncml']

        data_object = NetCDFDogma('ds', self.beliefs(), ds, namespaces=namespaces)
        return DSPair(ds, data_object)

class BaseSOSGCCheck(object):
    """
    Base class for SOS-GetCapabilities supporting Check Suites.
    """
    supported_ds = [SensorObservationService_1_0_0]

    def load_datapair(self, ds):
        data_object = MultipleXmlDogma('sos-gc', self.beliefs(), ds._capabilities, namespaces=get_namespaces())
        return DSPair(ds, data_object)

class BaseSOSDSCheck(object):
    """
    Base class for SOS-DescribeSensor supporting Check Suites.
    """
    supported_ds = [SensorML]

    def load_datapair(self, ds):
        data_object = MultipleXmlDogma('sos-ds', self.beliefs(), ds._root, namespaces=get_namespaces())
        return DSPair(ds, data_object)

class Result(object):
    """
    Holds the result of a check method.

    Stores such information as the check's value (True, False, a 2-tuple of (pass, total) or None for a skip),
    weight of the check, any granular messages, or a hierarchy of results.

    Stores the checker instance and the check method that produced this result.
    """
    def __init__(self, weight=BaseCheck.MEDIUM, value=None, name=None, msgs=None, children=None, checker=None, check_method=None):

        self.weight = weight
        self.value  = value
        self.name   = name
        self.msgs   = msgs or []

        self.children = children or []

        self.checker = checker
        self.check_method = check_method

    def __repr__(self):
        ret = "%s (*%s): %s" % (self.name, self.weight, self.value)

        if len(self.msgs):
            if len(self.msgs) == 1:
                ret += " (%s)" % self.msgs[0]
            else:
                ret += " (%d msgs)" % len(self.msgs)

        if len(self.children):
            ret += " (%d children)" % len(self.children)
            ret += "\n" + pprint.pformat(self.children)
        return ret

    def serialize(self):
        '''
        Returns a serializable dictionary that represents the result object
        '''
        return {
            'name' : self.name,
            'weight' : self.weight,
            'value' : self.value,
            'msgs' : self.msgs,
            'children' : [i.serialize() for i in self.children]
        }


def std_check_in(dataset_dogma, name, allowed_vals):
    """
    Returns 0 if attr not present, 1 if present but not in correct value, 2 if good
    """
    if not hasattr(dataset_dogma, name):
        return 0

    ret_val = 1
    try:
        if getattr(dataset_dogma, name) in allowed_vals:
            ret_val += 1
    except DogmaGetterSetterException:
        pass

    return ret_val

def std_check(dataset_dogma, name):
    if hasattr(dataset_dogma, name):
        getattr(dataset_dogma, name)
        return True

    return False

def check_has(priority=BaseCheck.HIGH):

    def _inner(func):
        def _dec(s, ds):
            list_vars = func(s, ds)

            ret_val = []
            for l in list_vars:
                msgs = []
                if isinstance(l, tuple):
                    name, allowed = l
                    res = std_check_in(ds.dogma, name, allowed)
                    if res == 0:
                        msgs.append("Attr %s not present" % name)
                    elif res == 1:
                        msgs.append("Attr %s present but not in expected value list (%s)" % (name, allowed))

                    ret_val.append(Result(priority, (res, 2), name, msgs))

                else:
                    res = std_check(ds.dogma, l)
                    if not res:
                        msgs = ["Attr %s not present" % l]
                    ret_val.append(Result(priority, res, l, msgs))

            return ret_val

        return wraps(func)(_dec)

    return _inner

def fix_return_value(v, method_name, method=None, checker=None):
    """
    Transforms scalar return values into Result.
    """
    method_name = (method_name or method.im_func.func_name).replace("check_", "")     # remove common check prefix

    if v is None or not isinstance(v, Result):
        v = Result(value=v, name=method_name)

    v.name         = v.name or method_name
    v.checker      = checker
    v.check_method = method

    return v

def score_group(group_name=None):
    def _inner(func):
        def _dec(s, ds):
            ret_val = func(s, ds)
            """
            if group_name != None and not isinstance(ret_val[0], tuple):
                return tuple([(group_name, ret_val[0])] + list(ret_val[1:]))
            """
            # multiple returns
            if not isinstance(ret_val, list):
                ret_val = [ret_val]

            def dogroup(r):
                cur_grouping = r.name

                if isinstance(cur_grouping, tuple):
                    cur_grouping = list(cur_grouping)
                elif not isinstance(cur_grouping, list):
                    cur_grouping = [cur_grouping]

                cur_grouping.insert(0, group_name)



                return Result(r.weight, r.value, tuple(cur_grouping), r.msgs)


            ret_val = map(lambda x: fix_return_value(x, func.func_name, func, s), ret_val)
            ret_val = map(dogroup, ret_val)

            return ret_val
        return wraps(func)(_dec)
    return _inner

