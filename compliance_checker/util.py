import isodate
import re
"""General purpose utility functions to aid in compliance checking tasks"""

def is_readable(word):
    """Checks if a word is "readable" by applying a set of arbitrary regex
    heuristics
    """
    re_format = re.compile("^['\"0-9A-Za-z_\/()-. ]+$")
    return re_format.match(word)

def datetime_is_iso(dt):
    """Attempts to parse a date formatted in ISO 8601 format"""
    try:
        isodate.parse_datetime(dt)
        return True, []
    except isodate.ISO8601Error:
        return False, ['Datetime provided is not in a valid ISO 8601 format']
