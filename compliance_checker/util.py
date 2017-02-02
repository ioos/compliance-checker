"""
General purpose utility functions to aid in compliance checking tasks
"""
import isodate
import arrow


def isstring(obj):
    try:
        return isinstance(obj, basestring)
    except NameError:
        return isinstance(obj, str)


def datetime_is_iso(dt):
    """Attempts to parse a date formatted in ISO 8601 format"""
    try:
        if len(dt) > 10:
            isodate.parse_datetime(dt)
        else:
            isodate.parse_date(dt)
        return True, []
    except isodate.ISO8601Error:
        return False, ['Datetime provided is not in a valid ISO 8601 format']


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
