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

# a list of tuples
# Each tuple contains the standard name followed by a regex for the formula
# terms
dimless_vertical_coordinates = [
    ("atmosphere_ln_pressure_coordinate",
     # "p0: var1 lev: var2"
     r'(p0): ([A-Za-z][A-Za-z0-9_]*) (lev): ([A-Za-z][A-Za-z0-9_]*)'),
    ("atmosphere_sigma_coordinate",
     # "sigma: var1 ps: var2 ptop: var3"
     r'(sigma): ([A-Za-z][A-Za-z0-9_]*) (ps): ([A-Za-z][A-Za-z0-9_]*) (ptop): ([A-Za-z][A-Za-z0-9_]*)'),
    ("atmosphere_hybrid_sigma_pressure_coordinate",
     # "a: var1 b: var2 ps: var3 p0: var4"
     r'(a): ([A-Za-z][A-Za-z0-9_]*) (b): ([A-Za-z][A-Za-z0-9_]*) (ps): ([A-Za-z][A-Za-z0-9_]*) (p0): ([A-Za-z][A-Za-z0-9_]*)'),
    ("atmosphere_hybrid_height_coordinate",
     # "a: var1 b: var2 orog: var3"
     r'(a): ([A-Za-z][A-Za-z0-9_]*) (b): ([A-Za-z][A-Za-z0-9_]*) (orog): ([A-Za-z][A-Za-z0-9_]*)'),
    ("atmosphere_sleve_coordinate",
     # "a: var1 b1: var2 b2: var3 ztop: var4 zsurf1: var5 zsurf2: var6"
     r'(a): ([A-Za-z][A-Za-z0-9_]*) (b1): ([A-Za-z][A-Za-z0-9_]*) (b2): ([A-Za-z][A-Za-z0-9_]*) (ztop): ([A-Za-z][A-Za-z0-9_]*) (zsurf1): ([A-Za-z][A-Za-z0-9_]*) (zsurf2): ([A-Za-z][A-Za-z0-9_]*)'),
    ("ocean_sigma_coordinate",
     # "sigma: var1 eta: var2 depth: var3"
     r'(sigma): ([A-Za-z][A-Za-z0-9_]*) (eta): ([A-Za-z][A-Za-z0-9_]*) (depth): ([A-Za-z][A-Za-z0-9_]*)'),
    ("ocean_s_coordinate",
     # "s: var1 eta: var2 depth: var3 a: var4 b: var5 depth_c: var6"
     r'(s): ([A-Za-z][A-Za-z0-9_]*) (eta): ([A-Za-z][A-Za-z0-9_]*) (depth): ([A-Za-z][A-Za-z0-9_]*) (a): ([A-Za-z][A-Za-z0-9_]*) (b): ([A-Za-z][A-Za-z0-9_]*) (depth_c): ([A-Za-z][A-Za-z0-9_]*)'),
    ("ocean_sigma_z_coordinate",
     # "sigma: var1 eta: var2 depth: var3 depth_c: var4 nsigma: var5 zlev: var6"
     r'(sigma): ([A-Za-z][A-Za-z0-9_]*) (eta): ([A-Za-z][A-Za-z0-9_]*) (depth): ([A-Za-z][A-Za-z0-9_]*) (depth_c): ([A-Za-z][A-Za-z0-9_]*) (nsigma): ([A-Za-z][A-Za-z0-9_]*) (zlev): ([A-Za-z][A-Za-z0-9_]*)'),
    ("ocean_double_sigma_coordinate",
     # "sigma: var1 depth: var2 z1: var3 z2: var4 a: var5 href: var6 k_c: var7"
     r'(sigma): ([A-Za-z][A-Za-z0-9_]*) (depth): ([A-Za-z][A-Za-z0-9_]*) (z1): ([A-Za-z][A-Za-z0-9_]*) (z2): ([A-Za-z][A-Za-z0-9_]*) (a): ([A-Za-z][A-Za-z0-9_]*) (href): ([A-Za-z][A-Za-z0-9_]*) (k_c): ([A-Za-z][A-Za-z0-9_]*)'),
    # This comes from CF 1.7 but is used in circulation so we include it
    ("ocean_s_coordinate_g1",
     # "s: var1 C: var2 eta: var3 depth: var4 depth_c: var5"
     r'(s): ([A-Za-z][A-Za-z0-9_]*) (C): ([A-Za-z][A-Za-z0-9_]*) (eta): ([A-Za-z][A-Za-z0-9_]*) (depth): ([A-Za-z][A-Za-z0-9_]*) (depth_c): ([A-Za-z][A-Za-z0-9_]*)'),
    ("ocean_s_coordinate_g2",
     # "s: var1 C: var2 eta: var3 depth: var4 depth_c: var5"
     r'(s): ([A-Za-z][A-Za-z0-9_]*) (C): ([A-Za-z][A-Za-z0-9_]*) (eta): ([A-Za-z][A-Za-z0-9_]*) (depth): ([A-Za-z][A-Za-z0-9_]*) (depth_c): ([A-Za-z][A-Za-z0-9_]*)'),
]
