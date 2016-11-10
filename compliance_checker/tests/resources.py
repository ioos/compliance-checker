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
    'grid-boundaries'                      : get_filename('tests/data/grid-boundaries.cdl'),
    'ncei_gold_point_1'                    : get_filename('tests/data/ncei_gold_point_1.cdl'),
    'ncei_gold_point_2'                    : get_filename('tests/data/ncei_gold_point_2.cdl'),
    'ghrsst'                               : get_filename('tests/data/20160919092000-ABOM-L3S_GHRSST-SSTfnd-AVHRR_D-1d_dn_truncate.cdl'),
    'climatology'                          : get_filename('tests/data/climatology.cdl'),
    'rotated_pole_grid'                    : get_filename('tests/data/rotated_pole_grid.cdl'),
    'bad_units'                            : get_filename('tests/data/bad_units.cdl'),
    'bad_reference'                        : get_filename('tests/data/bad_reference.cdl'),
    'coordinate_types'                     : get_filename('tests/data/coordinate_types.cdl'),
    'chap2'                                : get_filename('tests/data/chap2.cdl'),
    'vertical_coords'                      : get_filename('tests/data/vertical_coords.cdl'),
    'reduced_horizontal_grid'              : get_filename('tests/data/reduced_horizontal_grid.cdl'),
    'duplicate_axis'                       : get_filename('tests/data/duplicate_axis.cdl'),
    'multi-dim-coordinates'                : get_filename('tests/data/multi-dim-coordinates.cdl'),
    '2d-static-grid'                       : get_filename('tests/data/2d-static-grid.cdl'),
    '3d-static-grid'                       : get_filename('tests/data/3d-static-grid.cdl'),
    'illegal-vertical'                     : get_filename('tests/data/illegal-vertical.cdl'),
    'illegal-aux-coords'                   : get_filename('tests/data/illegal-aux-coords.cdl'),
    'bad-instance'                         : get_filename('tests/data/bad-instance.cdl'),
    'bounds_bad_order'                     : get_filename('tests/data/non-comp/bounds_bad_order.cdl'),
    'bounds_bad_num_coords'                : get_filename('tests/data/non-comp/bounds_bad_num_coords.cdl'),
    '1d_bound_bad'                         : get_filename('tests/data/non-comp/1d_bound_bad.cdl'),
    'cf_example_cell_measures'             : get_filename('tests/data/cf_example_cell_measures.cdl'),
    'h_point'                              : get_filename('tests/data/appendix_h/point.cdl'),
    'h_timeseries-incomplete'              : get_filename('tests/data/appendix_h/timeseries-incomplete.cdl'),
    'h_timeseries-orthogonal'              : get_filename('tests/data/appendix_h/timeseries-orthogonal.cdl'),
    'h_timeseries-single'                  : get_filename('tests/data/appendix_h/timeseries-single.cdl'),
    'sldmb_43093_agg'                      : get_filename('tests/data/examples/sldmb_43093_agg.cdl'),
    'hycom_global'                         : get_filename('tests/data/examples/hycom_global.cdl'),
    'ocos'                                 : get_filename('tests/data/examples/ocos.cdl'),
    'l01-met'                              : get_filename('tests/data/examples/l01-met.cdl'),
    'usgs_dem_saipan'                      : get_filename('tests/data/examples/usgs_dem_saipan.cdl'),
    'sp041'                                : get_filename('tests/data/examples/sp041.cdl'),
    '3mf07'                                : get_filename('tests/data/examples/3mf07.cdl'),
    'ooi_glider'                           : get_filename('tests/data/examples/ooi_glider.cdl'),
    'trajectory-complete'                  : get_filename('tests/data/trajectory-complete.cdl'),
    'trajectory-implied'                   : get_filename('tests/data/trajectory-implied.cdl'),
    'bad-trajectory'                       : get_filename('tests/data/bad-trajectory.cdl'),
    'swan'                                 : get_filename('tests/data/examples/swan.cdl'),
    'kibesillah'                           : get_filename('tests/data/examples/kibesillah.cdl'),
    'cf_example_cell_measures'             : get_filename('tests/data/examples/cf_example_cell_measures.cdl'),
    'pr_inundation'                        : get_filename('tests/data/examples/pr_inundation.cdl'),
    'fvcom'                                : get_filename('tests/data/examples/fvcom.cdl'),
    'ww3'                                  : get_filename('tests/data/examples/ww3.cdl'),
    'glcfs'                                : get_filename('tests/data/examples/glcfs.cdl'),
    'NCEI_profile_template_v2_0'           : get_filename('tests/data/NCEI_profile_template_v2.0_2016-09-22_181835.151325.cdl')
}
