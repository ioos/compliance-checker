#!/usr/bin/env python

"""
Compliance Checker
"""
import csv
import itertools
import pprint
import re
import warnings
from collections import defaultdict
from functools import wraps
from io import StringIO

import validators
from lxml import etree
from netCDF4 import Dataset
from owslib.namespaces import Namespaces
from owslib.swe.observation.sos100 import SensorObservationService_1_0_0
from owslib.swe.sensor.sml import SensorML

from compliance_checker import __version__
from compliance_checker.util import kvp_convert

# Python 3.5+ should work, also have a fallback
try:
    from re import Pattern

    re_pattern_type = Pattern
except ImportError:
    re_pattern_type = type(re.compile(""))


def get_namespaces():
    n = Namespaces()
    ns = n.get_namespaces(["ogc", "sml", "gml", "sos", "swe", "xlink"])
    ns["ows"] = n.get_namespace("ows110")
    return ns


def csv_splitter(input_string):
    """
    csv_splitter(input_string)

    Splits a string in CSV format and returns a flattened list

    Parameters:
    -----------
    input_string: str
        The string to be processed

    Returns:
    --------
    list of str
        A flattened list from the CSV processing contents
    """
    csv_contents = csv.reader(StringIO(input_string))
    return list(itertools.chain.from_iterable(csv_contents))


class ValidationObject:
    validator_fail_msg = ""
    expected_type = None

    def __init__(self, split_func=None):
        if split_func is None:
            self.split_func = lambda x: [x]
        else:
            self.split_func = split_func

    def validator_func(self, input_value):
        """
        validator_func(self, input_value)

        Function that should validate the result of a given input value
        """
        raise NotImplementedError

    def validate(self, input_name, input_value):
        if self.expected_type is not None:
            type_result = self.validate_type(input_name, input_value)
            if not type_result[0]:
                return type_result
        for processed_value in self.split_func(input_value):
            validator_result = self.validator_func(processed_value)
            if not validator_result:
                return False, [self.validator_fail_msg.format(input_name)]
        # if all pass, then we're good.
        return True, None

    def validate_type(self, input_name, input_value):
        if not isinstance(input_value, self.expected_type):
            expected_type_fmt = "Attribute {} should be instance of type {}"
            return (
                False,
                [expected_type_fmt.format(input_name, self.expected_type.__name__)],
            )
        else:
            return True, None


class EmailValidator(ValidationObject):
    validator_fail_msg = "{} must be a valid email address"
    expected_type = str

    def validator_func(self, input_value):
        return validators.email(input_value)


class RegexValidator(ValidationObject):
    expected_type = str
    validator_regex = r"^.+$"
    validator_fail_msg = "{} must not be an empty string"

    def validator_func(self, input_value):
        return bool(re.search(self.validator_regex, input_value))


class UrlValidator(ValidationObject):
    validator_fail_msg = "{} must be a valid URL"
    expected_type = str

    def validator_func(self, input_value):
        return bool(validators.url(input_value))


# Simple class for Generic File type (default to this if file not recognised)
class GenericFile:
    """
    Simple class for any file. Has same path lookup as netCDF4.Dataset.
    """

    def __init__(self, fpath):
        self.fpath = fpath

    def filepath(self):
        return self.fpath


class BaseCheck:
    HIGH = 3
    MEDIUM = 2
    LOW = 1

    _cc_checker_version = __version__
    _cc_display_headers = {3: "High Priority", 2: "Medium Priority", 1: "Low Priority"}

    supported_ds = []

    def setup(self, ds):
        """
        Common setup method for a Checker.

        Automatically run when running a CheckSuite. Define this method in your Checker class.
        """

    def __init__(self, options=None):
        self._defined_results = defaultdict(lambda: defaultdict(dict))
        if options is None:
            self.options = {}
        else:
            self.options = options

    def get_test_ctx(self, severity, name, variable=None):
        """
        Creates an existing TestCtx object in _defined_results dict if it does
        not exist for the current checker instance, or an returns the existing
        TestCtx for modification. Takes a severity level and name and uses the
        two element tuple formed by the arguments as a key into the dict.

        :param int severity: A BaseCheck severity level
        :param str name: The name of the check
        :rtype compliance_checker.base.TestCtx:
        :returns: A new or or existing `TestCtx` instance taken from this
                  instance's _defined_results dict
        """
        # Is it necessary to key out by severity?  Is severity level unique
        # per check?  If so, it could be eliminated from key hierarchy
        if severity not in self._defined_results[name][variable]:
            self._defined_results[name][variable][severity] = TestCtx(
                severity,
                name,
                variable=variable,
            )
        return self._defined_results[name][variable][severity]


class BaseNCCheck:
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


class BaseSOSGCCheck:
    """
    Base class for SOS-GetCapabilities supporting Check Suites.
    """

    supported_ds = [SensorObservationService_1_0_0]


class BaseSOSDSCheck:
    """
    Base class for SOS-DescribeSensor supporting Check Suites.
    """

    supported_ds = [SensorML]


class Result:
    """
    Holds the result of a check method.

    Stores such information as the check's value (True, False, a 2-tuple of (pass, total) or None for a skip),
    weight of the check, any granular messages, or a hierarchy of results. If given value is not a tuple, it
    is cast as a boolean using the bool() function.

    Stores the checker instance and the check method that produced this result.
    """

    def __init__(
        self,
        weight=BaseCheck.MEDIUM,
        value=None,
        name=None,
        msgs=None,
        children=None,
        checker=None,
        check_method=None,
        variable_name=None,
    ):
        self.weight = weight

        if value is None:
            self.value = None
        elif isinstance(value, tuple):
            if len(value) != 2:
                raise ValueError(
                    f"Result value must be 2-tuple or boolean! Got {value}",
                )
            self.value = value
        else:
            self.value = bool(value)
        self.name = name
        self.msgs = msgs or []

        self.children = children or []

        self.checker = checker
        self.check_method = check_method
        self.variable_name = variable_name

    def __repr__(self):
        ret = f"{self.name} (*{self.weight}): {self.value}"

        if len(self.msgs):
            if len(self.msgs) == 1:
                ret += f" ({self.msgs[0]})"
            else:
                ret += f" ({len(self.msgs)!s} msgs)"

        if len(self.children):
            ret += f" ({len(self.children)!s} children)"
            ret += "\n" + pprint.pformat(self.children)

        return ret

    def serialize(self):
        """
        Returns a serializable dictionary that represents the result object
        """
        return {
            "name": self.name,
            "weight": self.weight,
            "value": self.value,
            "msgs": self.msgs,
            "children": [i.serialize() for i in self.children],
        }

    def __eq__(self, other):
        return self.serialize() == other.serialize()


class TestCtx:
    """
    Simple struct object that holds score values and messages to compile into a result
    """

    def __init__(
        self,
        category=None,
        description="",
        out_of=0,
        score=0,
        messages=None,
        variable=None,
    ):
        self.category = category or BaseCheck.LOW
        self.out_of = out_of
        self.score = score
        self.messages = messages or []
        self.description = description or ""
        self.variable = variable

    def to_result(self):
        return Result(
            self.category,
            (self.score, self.out_of),
            self.description,
            self.messages,
            variable_name=self.variable,
        )

    def assert_true(self, test, message):
        """
        Increments score if test is true otherwise appends a message
        :rtype: bool
        :return: Boolean indicating whether test condition passed or not
        """
        self.out_of += 1

        if test:
            self.score += 1
        else:
            self.messages.append(message)

        return test

    def add_failure(self, message):
        """
        Adds a failure along with a message
        :rtype: None
        """
        self.assert_true(False, message)

    def add_pass(self):
        """
        Adds a pass condition
        :rtype: None
        """
        self.assert_true(True, None)


def std_check_in(base_context, name, allowed_vals):
    """
    Check that a value is contained within an iterable

    Parameters:
    -----------
    base_context: netCDF4.Dataset or netCDF4.variable
       The context in which to look for the attribute, either a
       netCDF4.Dataset or netCDF4.Variable.  If a netCDF dataset,
       the attribute is searched for in the global attributes.
       If a variable, the attributes are limited to those contained
       in the corresponding variable.
    name: str
       The name of the attribute to search for.
    allowed_vals: iterable
       An iterable, usually a set, which provides the possible valid values for
       the attribute.
    Returns:
    --------
    int
        Returns 0 if attr not present, 1 if present but not in correct value, 2
        if good.
    """
    if not hasattr(base_context, name):
        return 0

    ret_val = 1
    if base_context.getncattr(name) in allowed_vals:
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


def maybe_get_global_attr(attr_name, ds):
    if attr_name in ds.ncattrs():
        return True, ds.getncattr(attr_name)
    else:
        err_msg = "{} not present"
        return False, [err_msg.format(attr_name)]


def attr_check(kvp, ds, priority, ret_val, gname=None, var_name=None):
    """
    Handles attribute checks for simple presence of an attribute, presence of
    one of several attributes, and passing a validation function.  Returns a
    status along with an error message in the event of a failure.  Mutates
    ret_val parameter

    :param tuple(str, func) or str l: the attribute being checked
    :param netCDF4 dataset ds       : dataset being checked
    :param int priority             : priority level of check
    :param list ret_val             : result to be returned
    :param str or None gname        : group name assigned to a group of attribute Results
    :param str or None var_name     : name of the variable which contains this attribute
    """

    msgs = []
    name, other = kvp
    if var_name is not None:
        display_name = f"attribute {name} in variable {var_name}"
        base_context = ds.variables[var_name]
    else:
        display_name = name
        base_context = ds
    if other is None:
        res = std_check(ds, name)
        if not res:
            msgs = [f"{display_name} not present"]
        else:
            try:
                # see if this attribute is a string, try stripping
                # whitespace, and return an error if empty
                att_strip = base_context.getncattr(name).strip()
                if not att_strip:
                    res = False
                    msgs = [f"{display_name} is empty or completely whitespace"]
            # if not a string/has no strip method we should be OK
            except AttributeError:
                pass

        # gname arg allows the global attrs to be grouped together
        ret_val.append(
            Result(
                priority,
                value=res,
                name=gname if gname else name,
                msgs=msgs,
                variable_name=var_name,
            ),
        )
    elif hasattr(other, "__iter__"):
        # redundant, we could easily do this with a hasattr
        # check instead
        res = std_check_in(base_context, name, other)
        if res == 0:
            msgs.append(f"{display_name} not present")
        elif res == 1:
            msgs.append(
                f"{display_name} present, but not in expected value list ({sorted(other)})",
            )

        ret_val.append(
            Result(
                priority,
                (res, 2),
                gname if gname else name,  # groups Globals if supplied
                msgs,
                variable_name=var_name,
            ),
        )
    # if we have an XPath expression, call it on the document
    elif type(other) is etree.XPath:
        # TODO: store tree instead of creating it each time?
        # no execution path for variable
        res = xpath_check(ds._root, other)
        if not res:
            msgs = [f"XPath for {display_name} not found"]
        ret_val.append(
            Result(
                priority,
                res,
                gname if gname else name,
                msgs,
                variable_name=var_name,
            ),
        )
    # check if this is a subclass of ValidationObject
    elif isinstance(other, ValidationObject):
        attr_result = maybe_get_global_attr(name, ds)
        if not attr_result[0]:
            res_tup = attr_result
        else:
            check_val = attr_result[1]
            res_tup = other.validate(name, check_val)

        msgs = [] if res_tup[1] is None else res_tup[1]

        ret_val.append(Result(priority, res_tup[0], name, msgs))
    elif isinstance(other, re_pattern_type):
        attr_result = maybe_get_global_attr(name, ds)
        if not attr_result[0]:
            return attr_result
        else:
            check_val = attr_result[1]
        if not isinstance(check_val, str):
            res = False
            msgs = [f"{name} must be a string"]
        elif not other.search(check_val):
            res = False
            msgs = [f"{name} must match regular expression {other}"]
        else:
            res = True
            msgs = []

        ret_val.append(
            Result(priority, value=res, name=gname if gname else name, msgs=msgs),
        )

    # if the attribute is a function, call it
    # right now only supports single attribute
    # important note: current magic approach uses all functions
    # starting with "check".  Avoid naming check functions
    # starting with check if you want to pass them in with
    # a tuple to avoid them being checked more than once
    elif callable(other):
        # check that the attribute is actually present.
        # This reduces boilerplate in functions by not needing
        # to check whether the attribute is present every time
        # and instead focuses on the core functionality of the
        # test

        res = other(base_context)  # call the method on the dataset
        if not res:
            msgs = [f"{display_name} not present"]
            ret_val.append(
                Result(
                    priority,
                    res,
                    gname if gname else name,
                    msgs,
                    variable_name=var_name,
                ),
            )
        else:
            ret_val.append(res(priority))
    # unsupported second type in second
    else:
        raise TypeError(
            f"Second arg in tuple has unsupported type: {type(other)}",
        )

    return ret_val


def check_has(priority=BaseCheck.HIGH, gname=None):
    """Decorator to wrap a function to check if a dataset has given attributes.
    :param function func: function to wrap"""

    def _inner(func):
        def _dec(s, ds):
            attr_process = kvp_convert(func(s, ds))

            ret_val = []
            # could potentially run tests in parallel if we eliminated side
            # effects on `ret_val`
            for kvp in attr_process.items():
                # function mutates ret_val
                attr_check(kvp, ds, priority, ret_val, gname)
            return ret_val

        return wraps(func)(_dec)

    return _inner


def fix_return_value(v, method_name, method=None, checker=None):
    """
    Transforms scalar return values into Result.
    """
    # remove common check prefix
    method_name = (method_name or method.__func__.__name__).replace("check_", "")
    if v is None or not isinstance(v, Result):
        v = Result(value=v, name=method_name)

    v.name = v.name or method_name

    v.checker = checker
    v.check_method = method

    return v


def ratable_result(value, name, msgs, variable_name=None):
    """Returns a partial function with a Result that has not been weighted."""
    return lambda w: Result(w, value, name, msgs, variable_name=variable_name)


def score_group(group_name=None):
    """
    Warning this is deprecated as of Compliance Checker v3.2!

    Please do not using scoring groups and update your plugins
    if necessary
    """
    warnings.warn(
        "Score_group is deprecated as of Compliance Checker v3.2.",
        stacklevel=2,
    )

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

            ret_val = [fix_return_value(x, func.__name__, func, s) for x in ret_val]
            ret_val = list(map(dogroup, ret_val))

            return ret_val

        return wraps(func)(_dec)

    return _inner
