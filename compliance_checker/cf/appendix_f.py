#!/usr/bin/env python
'''
Appendix F. Grid Mappings
---
Each recognized grid mapping is described in one of the sections below. Each
section contains: the valid name that is used with the grid_mapping_name
attribute; a list of the specific attributes that may be used to assign values
to the mapping's parameters; the standard names used to identify the coordinate
variables that contain the mapping's independent variables; and references to
the mapping's definition or other information that may help in using the
mapping. Since the attributes used to set a mapping's parameters may be shared
among several mappings, their definitions are contained in a table in the final
section. The attributes which describe the ellipsoid and prime meridian may be
included, when applicable, with any grid mapping.

We have used the FGDC "Content Standard for Digital Geospatial Metadata" [FGDC]
as a guide in choosing the values for grid_mapping_name and the attribute names
for the parameters describing map projections.
'''

grid_mapping_names = [
    'albers_conical_equal_area',
    'azimuthal_equidistant',
    'lambert_azimuthal_equal_area',
    'lambert_conformal_conic',
    'lambert_cylindrical_equal_area',
    'latitude_longitude',
    'mercator',
    'orthographic',
    'polar_stereographic',
    'rotated_latitude_longitude',
    'stereographic',
    'transverse_mercator',
    'vertical_perspective' ]

grid_mapping_attrs = [
    'earth_radius',
    'false_easting',
    'false_northing',
    'grid_mapping_name',
    'grid_north_pole_latitude',
    'grid_north_pole_longitude',
    'inverse_flattening',
    'latitude_of_projection_origin',
    'longitude_of_central_meridian',
    'longitude_of_prime_meridian',
    'longitude_of_projection_origin',
    'north_pole_grid_longitude',
    'perspective_point_height',
    'scale_factor_at_central_meridian',
    'scale_factor_at_projection_origin',
    'semi_major_axis',
    'semi_minor_axis',
    'standard_parallel',
    'straight_vertical_longitude_from_pole'
]
