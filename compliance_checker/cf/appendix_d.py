#!/usr/bin/env python
'''
Appendix D compliance support for CF 1.6

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
'''

# Contains the standard name followed by the set of expected formula terms
dimless_vertical_coordinates = {
    "atmosphere_ln_pressure_coordinate": {'p0', 'lev'},
    "atmosphere_sigma_coordinate": {'sigma', 'ps', 'ptop'},
    "atmosphere_hybrid_sigma_pressure_coordinate":
       ({'a', 'b', 'ps'}, {'ap', 'b', 'ps'}),
    "atmosphere_hybrid_height_coordinate": {'a', 'b', 'orog'},
    "atmosphere_sleve_coordinate":
       {'a', 'b1', 'b2', 'ztop', 'zsurf1', 'zsurf2'},
    "ocean_sigma_coordinate": {'sigma', 'eta', 'depth'},
    "ocean_s_coordinate": {'s', 'eta', 'depth', 'a', 'b', 'depth_c'},
    "ocean_sigma_z_coordinate":
       {'sigma', 'eta', 'depth', 'depth_c', 'nsigma', 'zlev'},
    "ocean_double_sigma_coordinate":
       {'sigma', 'depth', 'z1', 'z2', 'a', 'href', 'k_c'},
    # This comes from CF 1.7 but is used in circulation so we include it
    # TODO (badams): include these *only* with the CF 1.7 checker once we get
    # around to writing it
    "ocean_s_coordinate_g1": {'s', 'C', 'eta', 'depth', 'depth_c'},
    "ocean_s_coordinate_g2": {'s', 'C', 'eta', 'depth', 'depth_c'}
 }

def no_missing_terms(formula_name, term_set):
    """
    Returns true if the set is not missing terms corresponding to the
    entries in Appendix D, False otherwise.  The set of terms should be exactly
    equal, and not contain more or less terms than expected.
    """
    reqd_terms = dimless_vertical_coordinates[formula_name]
    def has_all_terms(reqd_termset):
        return len(reqd_termset ^ term_set) == 0

    if isinstance(reqd_terms, set):
        return has_all_terms(reqd_terms)
    # if it's not a set, it's likely some other form of iterable with multiple
    # possible definitions i.e. a/ap are interchangeable in
    else:
        return any(has_all_terms(req) for req in reqd_terms)
