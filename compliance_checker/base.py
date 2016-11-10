#!/usr/bin/env python

"""
Compliance Checker
"""
from __future__ import unicode_literals

from functools import wraps
import pprint

from netCDF4 import Dataset
from owslib.swe.observation.sos100 import SensorObservationService_1_0_0
from owslib.swe.sensor.sml import SensorML
from owslib.namespaces import Namespaces
from compliance_checker import __version__
from distutils.version import StrictVersion as V
from lxml import etree
import sys


def get_namespaces():
    n = Namespaces()
    ns = n.get_namespaces(["ogc", "sml", "gml", "sos", "swe", "xlink"])
    ns["ows"] = n.get_namespace("ows110")
    return ns


class BaseCheck(object):
    HIGH   = 3
    MEDIUM = 2
    LOW    = 1

    _cc_checker_version = __version__

    supported_ds = []

    def setup(self, ds):
        """
        Common setup method for a Checker.

        Automatically run when running a CheckSuite. Define this method in your Checker class.
        """
        pass


class BaseNCCheck(object):
    """
    Base Class for NetCDF Dataset supporting Check Suites.
    """
    supported_ds = [Dataset]

    @classmethod
    def std_check_in(cls, dataset, name, allowed_vals):
        """
        Returns 0 if attr not present, 1 if present but not in correct value, 2 if good
        """
        if name not in dataset.ncattrs():
            return 0

        ret_val = 1
        if dataset.getncattr(name) in allowed_vals:
            ret_val += 1

        return ret_val

    @classmethod
    def std_check(cls, dataset, name):
        return name in dataset.ncattrs()


class BaseSOSGCCheck(object):
    """
    Base class for SOS-GetCapabilities supporting Check Suites.
    """
    supported_ds = [SensorObservationService_1_0_0]


class BaseSOSDSCheck(object):
    """
    Base class for SOS-DescribeSensor supporting Check Suites.
    """
    supported_ds = [SensorML]


class Result(object):
    """
    Holds the result of a check method.

    Stores such information as the check's value (True, False, a 2-tuple of (pass, total) or None for a skip),
    weight of the check, any granular messages, or a hierarchy of results. If given value is not a tuple, it
    is cast as a boolean using the bool() function.

    Stores the checker instance and the check method that produced this result.
    """

    def __init__(self,
                 weight=BaseCheck.MEDIUM,
                 value=None,
                 name=None,
                 msgs=None,
                 children=None,
                 checker=None,
                 check_method=None,
                 variable_name=None):

        self.weight = weight

        if value is None:
            self.value = None
        elif isinstance(value, tuple):
            assert len(value)==2, 'Result value must be 2-tuple or boolean!'
            self.value = value
        else:
            self.value = bool(value)
        self.name   = name
        self.msgs   = msgs or []

        self.children = children or []

        self.checker = checker
        self.check_method = check_method
        self.variable_name = variable_name

    def __repr__(self):
        ret = '{} (*{}): {}'.format(self.name, self.weight, self.value)

        if len(self.msgs):
            if len(self.msgs) == 1:
                ret += ' ({})'.format(self.msgs[0])
            else:
                ret += ' ({!s} msgs)'.format(len(self.msgs))

        if len(self.children):
            ret += ' ({!s} children)'.format(len(self.children))
            ret += '\n' + pprint.pformat(self.children)

        # python 2 requires repr to be an ASCII string
        # python 3 requires repr to be a unicode string
        if sys.version_info[0] < 3:
            return ret.encode("utf-8")

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

    def __eq__(self, other):
        return self.serialize() == other.serialize()


class TestCtx(object):
    '''
    Simple struct object that holds score values and messages to compile into a result
    '''
    def __init__(self, category=None, description='', out_of=0, score=0,
                 messages=None, variable=None):
        self.category = category or BaseCheck.LOW
        self.out_of = out_of
        self.score = score
        self.messages = messages or []
        self.description = description or ''
        self.variable = variable

    def to_result(self):
        return Result(self.category, (self.score, self.out_of), self.description, self.messages, variable_name=self.variable)

    def assert_true(self, test, message):
        '''
        Increments score if test is true otherwise appends a message
        '''
        self.out_of += 1

        if test:
            self.score += 1
        else:
            self.messages.append(message)


def std_check_in(dataset, name, allowed_vals):
    """
    Returns 0 if attr not present, 1 if present but not in correct value, 2 if good
    """
    if not hasattr(dataset, name):
        return 0

    ret_val = 1
    if getattr(dataset, name) in allowed_vals:
        ret_val += 1

    return ret_val


def std_check(dataset, name):
    if hasattr(dataset, name):
        getattr(dataset, name)
        return True

    return False

def xpath_check(tree, xpath):
    """Checks whether tree contains one or more elements matching xpath"""
    return len(xpath(tree)) > 0

def attr_check(l, ds, priority, ret_val):
    """
    Handles attribute checks for simple presence of an attribute, presence of
    one of several attributes, and passing a validation function.  Returns a
    status along with an error message in the event of a failure.  Mutates
    ret_val parameter
    """
    msgs = []
    if isinstance(l, tuple):
        name, other = l
        if hasattr(other, '__iter__'):
            # redundant, we could easily do this with a hasattr
            # check instead
            res = std_check_in(ds, name, other)
            if res == 0:
                msgs.append("Attr %s not present" % name)
            elif res == 1:
                msgs.append("Attr %s present, but not in expected value list (%s)" % (name, other))

            ret_val.append(Result(priority, (res, 2), name, msgs))
        # if we have an XPath expression, call it on the document
        elif type(other) is etree.XPath:
            # TODO: store tree instead of creating it each time?
            res = xpath_check(ds._root, other)
            if not res:
                msgs = ["XPath for {} not found".format(name)]
            ret_val.append(Result(priority, res, name, msgs))
        # if the attribute is a function, call it
        # right now only supports single attribute
        # important note: current magic approach uses all functions
        # starting with "check".  Avoid naming check functions
        # starting with check if you want to pass them in with
        # a tuple to avoid them being checked more than once
        elif hasattr(other, '__call__'):
            # check that the attribute is actually present.
            # This reduces boilerplate in functions by not needing
            # to check whether the attribute is present every time
            # and instead focuses on the core functionality of the
            # test
            res = std_check(ds, name)
            if not res:
                msgs = ["Attr %s not present" % name]
                ret_val.append(Result(priority, res, name, msgs))
            else:
                ret_val.append(other(ds)(priority))
        # unsupported second type in second
        else:
            raise TypeError("Second arg in tuple has unsupported type: {}".format(type(other)))

    else:
        res = std_check(ds, l)
        if not res:
            msgs = ["Attr %s not present" % l]
        else:
            try:
                # see if this attribute is a string, try stripping
                # whitespace, and return an error if empty
                att_strip = getattr(ds, l).strip()
                if not att_strip:
                    res = False
                    msgs = ["Attr %s is empty or completely whitespace" % l]
            # if not a string/has no strip method we should be OK
            except AttributeError:
                pass

        ret_val.append(Result(priority, res, l, msgs))

    return ret_val


def check_has(priority=BaseCheck.HIGH):

    def _inner(func):
        def _dec(s, ds):
            list_vars = func(s, ds)

            ret_val = []
            # could potentially run tests in parallel if we eliminated side
            # effects on `ret_val`
            for l in list_vars:
                # function mutates ret_val
                attr_check(l, ds, priority, ret_val)
            return ret_val

        return wraps(func)(_dec)

    return _inner


def fix_return_value(v, method_name, method=None, checker=None):
    """
    Transforms scalar return values into Result.
    """
    method_name = (method_name or method.__func__.__name__).replace("check_", "")     # remove common check prefix

    if v is None or not isinstance(v, Result):
        v = Result(value=v, name=method_name)

    v.name         = v.name or method_name
    v.checker      = checker
    v.check_method = method

    return v

def ratable_result(value, name, msgs):
    """Returns a partial function with a Result that has not been weighted."""
    return lambda w: Result(w, value, name, msgs)

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

            ret_val = [fix_return_value(x, func.__name__, func, s) for x in
                       ret_val]
            ret_val = list(map(dogroup, ret_val))

            return ret_val
        return wraps(func)(_dec)
    return _inner
