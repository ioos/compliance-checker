"""
General purpose utility functions to aid in compliance checking tasks
"""
import isodate
import arrow
from datetime import datetime

def isstring(obj):
    try:
        return isinstance(obj, basestring)
    except NameError:
        return isinstance(obj, str)


def datetime_is_iso(dt):
    """Attempts to parse a date formatted in ISO 8601 format"""
    try:
        dt = dateparse_iso(dt)
        return True, []
    except isodate.ISO8601Error:
        return False, ['Datetime provided is not in a valid ISO 8601 format']


def dateparse_iso(date_str):
    '''
    Returns a naive datetime. parsed from an ISO-8601 input string

    :param str date_str: An ISO-8601 string
    '''
    if len(date_str) > 10:
        return isodate.parse_datetime(date_str).replace(tzinfo=None)
    # Must be just a date. Parse as a python date, then convert to naive datetime
    date_iso = isodate.parse_date(date_str)
    return datetime.combine(date_iso, datetime.min.time()).replace(tzinfo=None)


def dateparse(date_str):
    '''
    Returns a datetime string parsed from an ISO-8601 input

    :param str date_str: An ISO-8601 string
    '''
    if isstring(date_str):
        if date_str.endswith('+00'):
            date_str = date_str.replace('+00', 'Z')
    arrow_obj = arrow.get(date_str)
    return arrow_obj.to('utc').naive
