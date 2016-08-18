from pkg_resources import resource_filename
import os
import subprocess

def get_filename(path):
    '''
    Returns the path to a valid dataset
    '''
    filename = resource_filename('compliance_checker', path)
    if not os.path.exists(filename):
        cdl_path = filename.replace('.nc', '.cdl')
        generate_dataset(cdl_path, filename)
    return filename

def generate_dataset(cdl_path, nc_path):
    subprocess.call(['ncgen','-o', nc_path, cdl_path])

STATIC_FILES = {
    'rutgers'                              : get_filename('tests/data/ru07-20130824T170228_rt0.nc'),
    'conv_multi'                           : get_filename('tests/data/conv_multi.nc'),
    'conv_bad'                             : get_filename('tests/data/conv_bad.nc'),
    'example-grid'                         : get_filename('tests/data/example-grid.nc'),
    'badname'                              : get_filename('tests/data/non-comp/badname.netcdf'),
    'bad'                                  : get_filename('tests/data/non-comp/bad.nc'),
    'dimensionless'                        : get_filename('tests/data/dimensionless.nc'),
    '2dim'                                 : get_filename('tests/data/2dim-grid.nc'),
    'bad2dim'                              : get_filename('tests/data/non-comp/bad2dim.nc'),
    'rhgrid'                               : get_filename('tests/data/rhgrid.nc'),
    'bad-rhgrid'                           : get_filename('tests/data/non-comp/bad-rhgrid.nc'),
    'bad_data_type'                        : get_filename('tests/data/bad_data_type.nc'),
    'mapping'                              : get_filename('tests/data/mapping.nc'),
    'bad_region'                           : get_filename('tests/data/bad_region.nc'),
    'featureType'                          : get_filename('tests/data/example-grid.nc'),
    'cont_ragged'                          : get_filename('tests/data/cont_ragged.nc'),
    'index_ragged'                         : get_filename('tests/data/index_ragged.nc'),
    'bad_missing_data'                     : get_filename('tests/data/bad_missing_data.nc'),
    'self-referencing-var'                 : get_filename('tests/data/self-referencing-var.nc'),
    'scalar_coordinate_variable'           : get_filename('tests/data/scalar_coordinate_variable.nc'),
    'coordinates_and_metadata'             : get_filename('tests/data/coordinates_and_metadata.nc'),
    'ints64'                               : get_filename('tests/data/ints64.nc'),
    'units_check'                          : get_filename('tests/data/units_check.nc'),
    'self_referencing'                     : get_filename('tests/data/non-comp/self_referencing.nc'),
    'time_units'                           : get_filename('tests/data/non-comp/time_units.nc'),
    'valid_coordinates'                    : get_filename('tests/data/valid_coordinates.nc'),
    'point'                                : get_filename('tests/data/point.nc'),
    'timeseries'                           : get_filename('tests/data/timeseries.nc'),
    'multi-timeseries-orthogonal'          : get_filename('tests/data/multi-timeseries-orthogonal.nc'),
    'multi-timeseries-incomplete'          : get_filename('tests/data/multi-timeseries-incomplete.nc'),
    'trajectory'                           : get_filename('tests/data/trajectory.nc'),
    'trajectory-single'                    : get_filename('tests/data/trajectory-single.nc'),
    'profile-orthogonal'                   : get_filename('tests/data/profile-orthogonal.nc'),
    'profile-incomplete'                   : get_filename('tests/data/profile-incomplete.nc'),
    'timeseries-profile-single-station'    : get_filename('tests/data/timeseries-profile-single-station.nc'),
    'timeseries-profile-multi-station'     : get_filename('tests/data/timeseries-profile-multi-station.nc'),
    'timeseries-profile-single-ortho-time' : get_filename('tests/data/timeseries-profile-single-ortho-time.nc'),
    'timeseries-profile-multi-ortho-time'  : get_filename('tests/data/timeseries-profile-multi-ortho-time.nc'),
    'timeseries-profile-ortho-depth'       : get_filename('tests/data/timeseries-profile-ortho-depth.nc'),
    'timeseries-profile-incomplete'        : get_filename('tests/data/timeseries-profile-incomplete.nc'),
    'trajectory-profile-orthogonal'        : get_filename('tests/data/trajectory-profile-orthogonal.nc'),
    'trajectory-profile-incomplete'        : get_filename('tests/data/trajectory-profile-incomplete.nc'),
    '2d-regular-grid'                      : get_filename('tests/data/2d-regular-grid.nc'),
    '3d-regular-grid'                      : get_filename('tests/data/3d-regular-grid.nc'),
}
