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
                    ret_val.append((priority, std_check_in(ds.dogma, l[0], l[1]), l[0]))
                else:
                    ret_val.append((priority, std_check(ds.dogma, l), l))

            return ret_val

        return wraps(func)(_dec)

    return _inner

def fix_return_value(v, method_name=None):
    """
    Fixes up the return value of any check method.

    Full return format is (weight, value, identifier).
    Check method authors only have to specify value, or weight/value as a 2-tuple.

    This method ensures it is a 3-tuple.
    """
    method_name = method_name.replace("check_", "")     # remove common check prefix

    if v is None or not isinstance(v, tuple):
        v = (BaseCheck.MEDIUM, v, method_name)
    elif isinstance(v, tuple) and len(v) == 2:
        v = (v[0], v[1], method_name)
    elif isinstance(v, tuple):
        pass
    else:
        raise StandardError("Unknown return value from check method %s: %s, expected either value or 2/3-tuple of (weight, value, optional name)" % (check_method.im_func.func_name, type(v)))

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
                cur_grouping = r[2] if len(r) >= 3 else []

                if isinstance(cur_grouping, tuple):
                    cur_grouping = list(cur_grouping)
                elif not isinstance(cur_grouping, list):
                    cur_grouping = [cur_grouping]

                cur_grouping.insert(0, group_name)

                return tuple(list(r[0:2]) + [tuple(cur_grouping)])

            ret_val = map(lambda x: fix_return_value(x, func.func_name), ret_val)
            ret_val = map(dogroup, ret_val)

            return ret_val
        return wraps(func)(_dec)
    return _inner

