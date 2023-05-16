#!/usr/bin/python

"""
Appendix E: Cell Methods
To be imported into cf.py upon initialization of a CF Checker class.
"""

cell_methods16 = {
    "point",
    "sum",
    "mean",
    "maximum",
    "minimum",
    "mid_range",
    "standard_deviation",
    "variance",
    "mode",
    "median",
    "sum_of_squares",
}

cell_methods17 = cell_methods16.union(
    {  # returns new set with elements from both
        "maximum_absolute_value",
        "minimum_absolute_value",
        "mean_absolute_value",
        "mean_of_upper_decile",
        "range",
        "root_mean_square",
    },
)
