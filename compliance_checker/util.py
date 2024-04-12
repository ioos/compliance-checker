"""
General purpose utility functions to aid in compliance checking tasks
"""

from collections import OrderedDict

import isodate
import pendulum


def datetime_is_iso(date_str):
    """Attempts to parse a date formatted in ISO 8601 format"""
    try:
        if len(date_str) > 10:
            isodate.parse_datetime(date_str)
        else:
            isodate.parse_date(date_str)
        return True, []
    # The following errors qualify as non ISO8601 format
    except (TypeError, AttributeError, ValueError, isodate.ISO8601Error):
        return False, ["Datetime provided is not in a valid ISO 8601 format"]


def dateparse(date_str):
    """
    Returns a naive datetime. parsed from an ISO-8601 input string

    :param str date_str: An ISO-8601 string
    """

    return pendulum.parse(date_str)


def kvp_convert(input_coll):
    """
    Converts a list of string attributes and/or tuples into an OrderedDict.
    If passed in an OrderedDict, function is idempotent.
    Key/value pairs map to `first_tuple_element` -> `second_tuple_element` if
    a tuple, or `scalar_value` -> None if not a tuple.

    :param input_coll: An iterable with string and/or 2-tuple elements
    :returns: collections.OrderedDict
    """
    if isinstance(input_coll, OrderedDict):
        return input_coll
    else:
        return OrderedDict(
            (thing, None) if not isinstance(thing, tuple) else (thing[0], thing[1])
            for thing in input_coll
        )
