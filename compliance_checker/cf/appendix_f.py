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

grid_mapping_attr_types16 = {
    'earth_radius': 'N',
    'false_easting': 'N',
    'false_northing': 'N',
    'grid_mapping_name': 'S',
    'grid_north_pole_latitude': 'N',
    'grid_north_pole_longitude': 'N',
    'inverse_flattening': 'N',
    'latitude_of_projection_origin': 'N',
    'longitude_of_central_meridian': 'N',
    'longitude_of_prime_meridian': 'N',
    'longitude_of_projection_origin': 'N',
    'north_pole_grid_longitude': 'N',
    'perspective_point_height': 'N',
    'scale_factor_at_central_meridian': 'N',
    'scale_factor_at_projection_origin': 'N',
    'semi_major_axis': 'N',
    'semi_minor_axis': 'N',
    'standard_parallel': 'N',
    'straight_vertical_longitude_from_pole': 'N'
}

grid_mapping_attr_types17 = grid_mapping_attr_types16.copy() # need shallow copy; update() returns None

grid_mapping_attr_types17.update({
    'azimuth_of_central_line': 'N',
    'crs_wkt': 'S',
    'geographic_crs_name': 'S',
    'geoid_name': 'S',
    'geopotential_datum_name': 'S',
    'prime_meridian_name': 'S',
    'projected_crs_name': 'S',
    'reference_ellipsoid_name': 'S',
    'towgs84': 'N'
    })


grid_mapping_dict16 = {
    'albers_conical_equal_area': [
        (
            'longitude_of_central_meridian',
            'latitude_of_projection_origin',
            'false_easting',
            'false_northing'
        ),
        (),
        (
            'projection_x_coordinate',
            'projection_y_coordinate'
        )
    ],
    'azimuthal_equidistant': [
        (
            'longitude_of_projection_origin',
            'latitude_of_projection_origin',
            'false_easting',
            'false_northing'
        ),
        (),
        (
            'projection_x_coordinate',
            'projection_y_coordinate'
        )
    ],
    'lambert_cylindrical_equal_area': [
        (
            'longitude_of_central_meridian',
            'false_easting',
            'false_northing'
        ),
        (),
        (
            'projection_x_coordinate',
            'projection_y_coordinate'
        ),
        (
            'standard_parallel',
            'scale_factor_at_projection_origin'
        )
    ],
    'lambert_azimuthal_equal_area': [
        (
            'longitude_of_projection_origin',
            'latitude_of_projection_origin',
            'false_easting',
            'false_northing'
        ),
        (),
        (
            'projection_x_coordinate',
            'projection_y_coordinate'
        )
    ],
    'lambert_conformal_conic': [
        (
            'standard_parallel',
            'longitude_of_central_meridian',
            'latitude_of_projection_origin',
            'false_easting',
            'false_northing'
        ),
        (),
        (
            'projection_x_coordinate',
            'projection_y_coordinate'
        )
    ],
    'latitude_longitude': [
        (),
        (),
        (
            'longitude',
            'latitude'
        )
    ],
    'mercator': [
        (
            'longitude_of_projection_origin',
            'false_easting',
            'false_northing'
        ),
        (),
        (
            'projection_x_coordinate',
            'projection_y_coordinate'
        ),
        (
            'standard_parallel',
            'scale_factor_at_projection_origin'
        )
    ],
    'orthographic': [
        (
            'longitude_of_projection_origin',
            'latitude_of_projection_origin',
            'false_easting',
            'false_northing'
        ),
        (),
        (
            'projection_x_coordinate',
            'projection_y_coordinate'
        )
    ],
    'polar_stereographic': [
        (
            'straight_vertical_longitude_from_pole',
            'latitude_of_projection_origin',
            'false_easting',
            'false_northing'
        ),
        (),
        (
            'projection_x_coordinate',
            'projection_y_coordinate'
        ),
        (
            'standard_parallel',
            'scale_factor_at_projection_origin'
        )
    ],
    'rotated_latitude_longitude': [
        (
            'grid_north_pole_latitude',
            'grid_north_pole_longitude'
        ),
        (
            'north_pole_grid_longitude'
        ),
        (
            'grid_latitude',
            'grid_longitude'
        )
    ],
    'stereographic': [
        (
            'longitude_of_projection_origin',
            'latitude_of_projection_origin',
            'scale_factor_at_projection_origin',
            'false_easting',
            'false_northing'
        ),
        (),
        (
            'projection_x_coordinate',
            'projection_y_coordinate'
        )
    ],
    'transverse_mercator': [
        (
            'scale_factor_at_central_meridian',
            'longitude_of_central_meridian',
            'latitude_of_projection_origin',
            'false_easting',
            'false_northing'
        ),
        (),
        (
            'projection_x_coordinate',
            'projection_y_coordinate'
        )
    ],
    'vertical_perspective': [
        (
            'longitude_of_projection_origin',
            'latitude_of_projection_origin',
            'perspective_point_height',
            'false_easting',
            'false_northing'
        ),
        (),
        (
            'projection_x_coordinate',
            'projection_y_coordinate'
        )
    ]
}

grid_mapping_dict17 = grid_mapping_dict16.copy() # need shallow copy; update() returns None
grid_mapping_dict17.update({
    'geostationary': [
        (
            'latitude_of_projection_origin',
            'longitude_of_projection_origin',
            'perspective_point_height',
            'false_easting',
            'false_northing'
        )
    ],
    'oblique_mercator': [
        (
            'azimuth',
            'latitude_of_projection_origin',
            'longitude_of_projection_origin',
            'scale_factor_at_projection_origin',
            'false_easting',
            'false_northing'
        )
    ],
    'sinusoidal': [
        (
            'longitude_of_projection_origin',
            'false_easting',
            'false_northing'
        )
    ]
})
