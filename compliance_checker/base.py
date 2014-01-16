#!/usr/bin/env python

"""
Compliance Checker
"""

from functools import wraps
import collections
import os
import os.path
import itertools
import pprint

from lxml import etree
from wicken.exceptions import DogmaGetterSetterException
from udunitspy import Unit, UdunitsError, Converter


class BaseCheck(object):
    HIGH   = 3
    MEDIUM = 2
    LOW    = 1

    @classmethod
    def beliefs(cls):
        raise NotImplementedError("Define this in your derived Checker class")

    def setup(self, ds):
        """
        Common setup method for a Checker.

        Automatically run when running a CheckSuite. Define this method in your Checker class.
        """
        pass

class Result(object):
    """
    Holds the result of a check method.

    Stores such information as the check's value (True, False, a 2-tuple of (pass, total) or None for a skip),
    weight of the check, any granular messages, or a hierarchy of results.
    """
    def __init__(self, weight=BaseCheck.MEDIUM, value=None, name=None, msgs=None, children=None):
        
        self.weight = weight
        self.value  = value
        self.name   = name
        self.msgs   = msgs or []

        self.children = children or []



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

class StandardNameTable(object):

    class NameEntry(object):
        def __init__(self, entrynode):
            self.canonical_units = self._get(entrynode, 'canonical_units', True)
            self.grib            = self._get(entrynode, 'grib')
            self.amip            = self._get(entrynode, 'amip')
            self.description     = self._get(entrynode, 'description')

        def _get(self, entrynode, attrname, required=False):
            vals = entrynode.xpath(attrname)
            if len(vals) > 1:
                raise StandardError("Multiple attrs (%s) found" % attrname)
            elif required and len(vals) == 0:
                raise StandardError("Required attr (%s) not found" % attrname)

            return vals[0].text

    def __init__(self, filename):
        if not os.path.isfile(filename):
            raise StandardError('File not found')

        parser = etree.XMLParser(remove_blank_text=True)
        self._tree = etree.parse(filename, parser)
        self._root = self._tree.getroot()

        # generate and save a list of all standard names in file
        self._names = [node.get('id') for node in self._root.iter('entry')]
        self._aliases = [node.get('id') for node in self._root.iter('alias')]

    def __len__(self):
        return len(self._names) + len(self._aliases)

    def __getitem__(self, key):
        if not (key in self._names or key in self._aliases):
            raise KeyError("%s not found in standard name table" % key)

        if key in self._aliases:
            idx = self._aliases.index(key)
            entryids = self._root.xpath('alias')[idx].xpath('entry_id')

            if len(entryids) != 1:
                raise StandardError("Inconsistency in standard name table, could not lookup alias for %s" % key)

            key = entryids[0].text

        if not key in self._names:
            raise KeyError("%s not found in standard name table" % key)

        idx = self._names.index(key)
        entry = self.NameEntry(self._root.xpath('entry')[idx])
        return entry

    def __contains__(self, key):
        return key in self._names or key in self._aliases

    def __iter__(self):
        return iter(itertools.chain(self._names, self._aliases))

def std_check_in(dataset_dogma, name, allowed_vals):
    #return name in dataset_dogma.variables and dataset_dogma.variables[name] in allowed_vals
    try:
        return hasattr(dataset_dogma, name) and getattr(dataset_dogma, name) in allowed_vals
    except DogmaGetterSetterException:
        pass

    return False

def std_check(dataset_dogma, name):
    #return name in dataset_dogma.variables
    return hasattr(dataset_dogma, name)

def check_has(priority=BaseCheck.HIGH):

    def _inner(func):
        def _dec(s, ds):
            list_vars = func(s, ds)

            ret_val = []
            for l in list_vars:
                if isinstance(l, tuple):
                    name, allowed = l
                    ret_val.append(Result(priority, std_check_in(ds.dogma, name, allowed), name))

                else:
                    ret_val.append(Result(priority, std_check(ds.dogma, l), l))


            return ret_val

        return wraps(func)(_dec)

    return _inner

def fix_return_value(v, method_name=None):
    """
    Transforms scalar return values into Result.
    """
    method_name = method_name.replace("check_", "")     # remove common check prefix

    if v is None or not isinstance(v, Result):

        v = Result(value=v, name=method_name)
    v.name = v.name or method_name


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


            ret_val = map(lambda x: fix_return_value(x, func.func_name), ret_val)
            ret_val = map(dogroup, ret_val)

            return ret_val
        return wraps(func)(_dec)
    return _inner

def units_known(units):
    try:
        Unit(str(units))
    except UdunitsError:
        return False
    return True

def units_convertible(units1, units2, reftimeistime=True):
    try:
        Converter(str(units1), str(units2))
    except UdunitsError:
        return False

    return True

def units_temporal(units):
    r = False
    try:
        u = Unit('seconds since 1900-01-01')
        r = u.are_convertible(str(units))
    except UdunitsError:
        return False
    return r
