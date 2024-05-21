#!/usr/bin/env python
"""
Appendix D compliance support for CF 1.6 and CF 1.7

The definitions given here allow an application to compute dimensional
coordinate values from the dimensionless ones and associated variables. The
formulas are expressed for a gridpoint (n,k,j,i) where i and j are the
horizontal indices, k is the vertical index and n is the time index. A
coordinate variable is associated with its definition by the value of the
standard_name attribute. The terms in the definition are associated with file
variables by the formula_terms attribute. The formula_terms attribute takes a
string value, the string being comprised of blank-separated elements of the form
"term: variable", where term is a keyword that represents one of the terms in
the definition, and variable is the name of the variable in a netCDF file that
contains the values for that term. The order of elements is not significant.

The gridpoint indices are not formally part of the definitions, but are included
to illustrate the indices that might be present in the file variables. For
example, a vertical coordinate whose definition contains a time index is not
necessarily time dependent in all netCDF files. Also, the definitions are given
in general forms that may be simplified by omitting certain terms. A term that
is omitted from the formula_terms attribute should be assumed to be zero.
"""

# Contains the standard name followed by a 2-tuple:
# (the set of expected formula terms, set of computed_standard_name(s)). Most
# vertical coordinates only have one computed_standard_name, but some have
# multiple acceptable values.
ocean_computed_standard_names = {
    "altitude": {
        "zlev": "altitude",
        "eta": "sea_surface_height_above_geoid",
        "depth": "sea_floor_depth_below_geoid",
    },
    "height_above_geopotential_datum": {
        "zlev": "height_above_geopotential_datum",
        "eta": "sea_surface_height_above_geopotential_datum",
        "depth": "sea_floor_depth_below_geopotential_datum",
    },
    "height_above_reference_ellipsoid": {
        "zlev": "height_above_reference_ellipsoid",
        "eta": "sea_surface_height_above_reference_ellipsoid",
        "depth": "sea_floor_depth_below_reference_ellipsoid",
    },
    "height_above_mean_sea_level": {
        "zlev": "height_above_mean_sea_level",
        "eta": "sea_surface_height_above_mean_sea_level",
        "depth": "sea_floor_depth_below_mean_sea_level",
    },
}

dimless_vertical_coordinates_1_6 = {  # only for CF-1.6
    "atmosphere_ln_pressure_coordinate": ({"p0", "lev"}, {"air_pressure"}),
    "atmosphere_sigma_coordinate": ({"sigma", "ps", "ptop"}, {"air_pressure"}),
    "atmosphere_hybrid_sigma_pressure_coordinate": (
        ({"a", "b", "ps"}, {"ap", "b", "ps"}),
        {"air_pressure"},
    ),
    "atmosphere_hybrid_height_coordinate": (
        {"a", "b", "orog"},
        {"altitude", "height_above_geopotential_datum"},
    ),
    "atmosphere_sleve_coordinate": (
        {"a", "b1", "b2", "ztop", "zsurf1", "zsurf2"},
        {"altitude", "height_above_geopotential_datum"},
    ),
    "ocean_sigma_coordinate": (
        {"sigma", "eta", "depth"},
        ocean_computed_standard_names,
    ),
    "ocean_s_coordinate": (
        {"s", "eta", "depth", "a", "b", "depth_c"},
        ocean_computed_standard_names,
    ),
    "ocean_sigma_z_coordinate": (
        {"sigma", "eta", "depth", "depth_c", "zlev"},
        ocean_computed_standard_names,
    ),
    "ocean_double_sigma_coordinate": (
        {"sigma", "depth", "z1", "z2", "a", "href", "k_c"},
        ocean_computed_standard_names,
    ),
}

dimless_vertical_coordinates_1_7 = (
    dimless_vertical_coordinates_1_6.copy()
)  # shallow copy
dimless_vertical_coordinates_1_7.update(
    {  # extends 1.6
        "ocean_s_coordinate_g1": (
            {"s", "C", "eta", "depth", "depth_c"},
            ocean_computed_standard_names,
        ),
        "ocean_s_coordinate_g2": (
            {"s", "C", "eta", "depth", "depth_c"},
            ocean_computed_standard_names,
        ),
    },
)


def no_missing_terms(formula_name, term_set, dimless_vertical_coordinates):
    """
    Returns true if the set is not missing terms corresponding to the
    entries in Appendix D, False otherwise.  The set of terms should be exactly
    equal, and not contain more or less terms than expected.
    """
    reqd_terms = dimless_vertical_coordinates[formula_name][0]

    def has_all_terms(reqd_termset):
        return len(reqd_termset ^ term_set) == 0

    if isinstance(reqd_terms, set):
        return has_all_terms(reqd_terms)
    # if it's not a set, it's likely some other form of iterable with multiple
    # possible definitions i.e. a/ap are interchangeable in
    else:
        return any(has_all_terms(req) for req in reqd_terms)
