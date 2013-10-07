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

    def beliefs(self):
        raise NotImplementedError("Define this in your derived Checker class")

def std_check_in(dataset, name, allowed_vals):
    #return name in dataset.variables and dataset.variables[name] in allowed_vals
    try:
        return hasattr(dataset, name) and getattr(dataset, name) in allowed_vals
    except DogmaGetterSetterException:
        pass

    return False

def std_check(dataset, name):
    #return name in dataset.variables
    return hasattr(dataset, name)

def check_has(priority=BaseCheck.HIGH):

    def _inner(func):
        def _dec(s, ds):
            list_vars = func(s, ds)

            ret_val = []
            for l in list_vars:
                if isinstance(l, tuple):
                    ret_val.append((priority, std_check_in(ds, l[0], l[1]), l[0]))
                else:
                    ret_val.append((priority, std_check(ds, l), l))

            return ret_val

        return wraps(func)(_dec)

    return _inner

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

            ret_val = map(dogroup, ret_val)

            return ret_val
        return wraps(func)(_dec)
    return _inner

