"""
General purpose utility functions to aid in compliance checking tasks
"""
import isodate
import pendulum


def isstring(obj):
    try:
        return isinstance(obj, basestring)
    except NameError:
        return isinstance(obj, str)


def datetime_is_iso(date_str):
    """Attempts to parse a date formatted in ISO 8601 format"""
    try:
        if len(date_str) > 10:
            dt = isodate.parse_datetime(date_str)
        else:
            dt = isodate.parse_date(date_str)
        return True, []
    except:  # Any error qualifies as not ISO format
        return False, ['Datetime provided is not in a valid ISO 8601 format']


def dateparse(date_str):
    '''
    Returns a naive datetime. parsed from an ISO-8601 input string

    :param str date_str: An ISO-8601 string
    '''

    return pendulum.parse(date_str)

