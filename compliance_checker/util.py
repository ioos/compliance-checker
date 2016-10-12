"""
General purpose utility functions to aid in compliance checking tasks
"""
import isodate
import arrow


def datetime_is_iso(dt):
    """Attempts to parse a date formatted in ISO 8601 format"""
    try:
        isodate.parse_datetime(dt)
        return True, []
    except isodate.ISO8601Error:
        return False, ['Datetime provided is not in a valid ISO 8602 format']


def dateparse(date_str):
    '''
    Returns a datetime string parsed from an ISO-8601 input

    :param str date_str: An ISO-8601 string
    '''
    if isinstance(date_str, basestring):
        if date_str.endswith('+00'):
            date_str = date_str.replace('+00', 'Z')
    arrow_obj = arrow.get(date_str)
    return arrow_obj.to('utc').naive
