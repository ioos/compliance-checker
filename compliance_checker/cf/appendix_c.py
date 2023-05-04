"""Appendix C - Standard Name modifiers"""

# Dict of standard name modifiers with values for units.  "u" indicates to use
# the same units as the standard name canonical units, "1" is unitless for
# observation counts, and None is used for status_flag, which expects units
# not to be present
valid_modifiers = {
    "detection_minimum": "u",
    "number_of_observations": "1",
    "standard_error": "u",
    "status_flag": None,
}
