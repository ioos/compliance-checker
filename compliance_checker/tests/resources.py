from pkg_resources import resource_filename
import os
import subprocess


def get_filename(path):
    '''
    Returns the path to a valid dataset
    '''
    filename = resource_filename('compliance_checker', path)
    nc_path = filename.replace('.cdl', '.nc')
    if not os.path.exists(nc_path):
        generate_dataset(filename, nc_path)
    return nc_path


def generate_dataset(cdl_path, nc_path):
    subprocess.call(['ncgen', '-o', nc_path, cdl_path])

STATIC_FILES = {
    'rutgers'                              : get_filename('tests/data/ru07-20130824T170228_rt0.cdl'),
    'conv_multi'                           : get_filename('tests/data/conv_multi.cdl'),
    'conv_bad'                             : get_filename('tests/data/conv_bad.cdl'),
    'example-grid'                         : get_filename('tests/data/example-grid.cdl'),
    'badname'                              : get_filename('tests/data/non-comp/badname.netcdf'),
    'bad'                                  : get_filename('tests/data/non-comp/bad.cdl'),
    'dimensionless'                        : get_filename('tests/data/dimensionless.cdl'),
    '2dim'                                 : get_filename('tests/data/2dim-grid.cdl'),
    'bad2dim'                              : get_filename('tests/data/non-comp/bad2dim.cdl'),
    'rhgrid'                               : get_filename('tests/data/rhgrid.cdl'),
    'bad-rhgrid'                           : get_filename('tests/data/non-comp/bad-rhgrid.cdl'),
    'bad_data_type'                        : get_filename('tests/data/bad_data_type.cdl'),
    'mapping'                              : get_filename('tests/data/mapping.cdl'),
    'bad_region'                           : get_filename('tests/data/bad_region.cdl'),
    'featureType'                          : get_filename('tests/data/example-grid.cdl'),
    'cont_ragged'                          : get_filename('tests/data/cont_ragged.cdl'),
    'index_ragged'                         : get_filename('tests/data/index_ragged.cdl'),
    'bad_missing_data'                     : get_filename('tests/data/bad_missing_data.cdl'),
    'self-referencing-var'                 : get_filename('tests/data/self-referencing-var.cdl'),
    'scalar_coordinate_variable'           : get_filename('tests/data/scalar_coordinate_variable.cdl'),
    'coordinates_and_metadata'             : get_filename('tests/data/coordinates_and_metadata.cdl'),
    'ints64'                               : get_filename('tests/data/ints64.cdl'),
    'units_check'                          : get_filename('tests/data/units_check.cdl'),
    'self_referencing'                     : get_filename('tests/data/non-comp/self_referencing.cdl'),
    'time_units'                           : get_filename('tests/data/non-comp/time_units.cdl'),
    'valid_coordinates'                    : get_filename('tests/data/valid_coordinates.cdl'),
    'point'                                : get_filename('tests/data/point.cdl'),
    'timeseries'                           : get_filename('tests/data/timeseries.cdl'),
    'multi-timeseries-orthogonal'          : get_filename('tests/data/multi-timeseries-orthogonal.cdl'),
    'multi-timeseries-incomplete'          : get_filename('tests/data/multi-timeseries-incomplete.cdl'),
    'trajectory'                           : get_filename('tests/data/trajectory.cdl'),
    'trajectory-single'                    : get_filename('tests/data/trajectory-single.cdl'),
    'profile-orthogonal'                   : get_filename('tests/data/profile-orthogonal.cdl'),
    'profile-incomplete'                   : get_filename('tests/data/profile-incomplete.cdl'),
    'timeseries-profile-single-station'    : get_filename('tests/data/timeseries-profile-single-station.cdl'),
    'timeseries-profile-multi-station'     : get_filename('tests/data/timeseries-profile-multi-station.cdl'),
    'timeseries-profile-single-ortho-time' : get_filename('tests/data/timeseries-profile-single-ortho-time.cdl'),
    'timeseries-profile-multi-ortho-time'  : get_filename('tests/data/timeseries-profile-multi-ortho-time.cdl'),
    'timeseries-profile-ortho-depth'       : get_filename('tests/data/timeseries-profile-ortho-depth.cdl'),
    'timeseries-profile-incomplete'        : get_filename('tests/data/timeseries-profile-incomplete.cdl'),
    'trajectory-profile-orthogonal'        : get_filename('tests/data/trajectory-profile-orthogonal.cdl'),
    'trajectory-profile-incomplete'        : get_filename('tests/data/trajectory-profile-incomplete.cdl'),
    '2d-regular-grid'                      : get_filename('tests/data/2d-regular-grid.cdl'),
    '3d-regular-grid'                      : get_filename('tests/data/3d-regular-grid.cdl'),
}
