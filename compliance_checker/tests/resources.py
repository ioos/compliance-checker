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
    'rutgers'                    : get_filename('tests/data/ru07-20130824T170228_rt0.nc'),
    'conv_multi'                 : get_filename('tests/data/conv_multi.nc'),
    'conv_bad'                   : get_filename('tests/data/conv_bad.nc'),
    'example-grid'               : get_filename('tests/data/example-grid.nc'),
    'badname'                    : get_filename('tests/data/non-comp/badname.netcdf'),
    'bad'                        : get_filename('tests/data/non-comp/bad.nc'),
    'dimensionless'              : get_filename('tests/data/dimensionless.nc'),
    '2dim'                       : get_filename('tests/data/2dim-grid.nc'),
    'bad2dim'                    : get_filename('tests/data/non-comp/bad2dim.nc'),
    'rhgrid'                     : get_filename('tests/data/rhgrid.nc'),
    'bad-rhgrid'                 : get_filename('tests/data/non-comp/bad-rhgrid.nc'),
    'bad_data_type'              : get_filename('tests/data/bad_data_type.nc'),
    'mapping'                    : get_filename('tests/data/mapping.nc'),
    'bad_region'                 : get_filename('tests/data/bad_region.nc'),
    'featureType'                : get_filename('tests/data/example-grid.nc'),
    'cont_ragged'                : get_filename('tests/data/cont_ragged.nc'),
    'index_ragged'               : get_filename('tests/data/index_ragged.nc'),
    'bad_missing_data'           : get_filename('tests/data/bad_missing_data.nc'),
    'self-referencing-var'       : get_filename('tests/data/self-referencing-var.nc'),
    'scalar_coordinate_variable' : get_filename('tests/data/scalar_coordinate_variable.nc'),
    'coordinates_and_metadata'   : get_filename('tests/data/coordinates_and_metadata.nc'),
    'ints64'                     : get_filename('tests/data/ints64.nc'),
    'units_check'                : get_filename('tests/data/units_check.nc'),
    'self_referencing'           : get_filename('tests/data/non-comp/self_referencing.nc'),
    'time_units'                 : get_filename('tests/data/non-comp/time_units.nc'),
    'valid_coordinates'          : get_filename('tests/data/valid_coordinates.nc')
              }
