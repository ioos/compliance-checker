#!/usr/bin/env python

"""
Compliance Checker
"""

from functools import wraps
import collections
from wicken.exceptions import DogmaGetterSetterException

class BaseCheck(object):
    HIGH   = 3
    MEDIUM = 2
    LOW    = 1

    @classmethod
    def beliefs(cls):
        raise NotImplementedError("Define this in your derived Checker class")

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
            ret += " (%d msgs)" % len(self.msgs)

        if len(self.children):
            ret += " (%d children)" % len(self.children)

        return ret

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

