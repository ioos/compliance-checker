import subprocess
from importlib.resources import files


def get_filename(path):
    """
    Returns the path to a valid dataset
    """
    filename = files("compliance_checker") / path
    nc_path = filename.with_suffix(".nc")
    if not nc_path.exists():
        generate_dataset(filename, nc_path)
    return nc_path


def generate_dataset(cdl_path, nc_path):
    subprocess.call(["ncgen", "-4", "-o", nc_path, cdl_path])


STATIC_FILES = {
    # TODO: add defaulltdict implementation for default files, etc
    "bad": get_filename("tests/data/non-comp/bad.cdl"),
    "badname": get_filename("tests/data/non-comp/badname.netcdf"),
    "bad-instance": get_filename("tests/data/bad-instance.cdl"),
    "bad-rhgrid": get_filename("tests/data/non-comp/bad-rhgrid.cdl"),
    "bad-trajectory": get_filename("tests/data/bad-trajectory.cdl"),
    "bad_cell_measure1": get_filename("tests/data/bad_cell_measure1.cdl"),
    "bad_cell_measure2": get_filename("tests/data/bad_cell_measure2.cdl"),
    "bad_cf_role": get_filename("tests/data/bad_cf_role.cdl"),
    "bad_data_type": get_filename("tests/data/bad_data_type.cdl"),
    "bad_missing_data": get_filename("tests/data/bad_missing_data.cdl"),
    "bad_reference": get_filename("tests/data/bad_reference.cdl"),
    "bad_region": get_filename("tests/data/bad_region.cdl"),
    "bad_units": get_filename("tests/data/bad_units.cdl"),
    "bad2dim": get_filename("tests/data/non-comp/bad2dim.cdl"),
    "bounds_bad_order": get_filename("tests/data/non-comp/bounds_bad_order.cdl"),
    "bounds_bad_num_coords": get_filename(
        "tests/data/non-comp/bounds_bad_num_coords.cdl",
    ),
    "cell_measure": get_filename("tests/data/cell_measure.cdl"),
    "cf_example_cell_measures": get_filename(
        "tests/data/examples/cf_example_cell_measures.cdl",
    ),
    "chap2": get_filename("tests/data/chap2.cdl"),
    "climatology": get_filename("tests/data/climatology.cdl"),
    "cont_ragged": get_filename("tests/data/cont_ragged.cdl"),
    "conv_multi": get_filename("tests/data/conv_multi.cdl"),
    "conv_bad": get_filename("tests/data/conv_bad.cdl"),
    "coordinates_and_metadata": get_filename("tests/data/coordinates_and_metadata.cdl"),
    "coordinate_types": get_filename("tests/data/coordinate_types.cdl"),
    "dimensionless": get_filename("tests/data/dimensionless.cdl"),
    "duplicate_axis": get_filename("tests/data/duplicate_axis.cdl"),
    "dimension_order": get_filename("tests/data/dimension_order.cdl"),
    "example-grid": get_filename("tests/data/example-grid.cdl"),
    "featureType": get_filename("tests/data/example-grid.cdl"),
    "forecast_reference": get_filename("tests/data/forecast_reference.cdl"),
    "fvcom": get_filename("tests/data/examples/fvcom.cdl"),
    "ghrsst": get_filename(
        "tests/data/20160919092000-ABOM-L3S_GHRSST-SSTfnd-AVHRR_D-1d_dn_truncate.cdl",
    ),
    "glcfs": get_filename("tests/data/examples/glcfs.cdl"),
    "grid-boundaries": get_filename("tests/data/grid-boundaries.cdl"),
    "grid_mapping_coordinates": get_filename("tests/data/grid_mapping_coordinates.cdl"),
    "hycom_global": get_filename("tests/data/examples/hycom_global.cdl"),
    "h_point": get_filename("tests/data/appendix_h/point.cdl"),
    "h_timeseries-incomplete": get_filename(
        "tests/data/appendix_h/timeseries-incomplete.cdl",
    ),
    "h_timeseries-orthogonal": get_filename(
        "tests/data/appendix_h/timeseries-orthogonal.cdl",
    ),
    "h_timeseries-single": get_filename("tests/data/appendix_h/timeseries-single.cdl"),
    "illegal-vertical": get_filename("tests/data/illegal-vertical.cdl"),
    "illegal-aux-coords": get_filename("tests/data/illegal-aux-coords.cdl"),
    "index_ragged": get_filename("tests/data/index_ragged.cdl"),
    "index_ragged2": get_filename("tests/data/index_ragged2.cdl"),
    "indexed_ragged_domain": get_filename("tests/data/indexed_ragged_domain.cdl"),
    "ints64": get_filename("tests/data/ints64.cdl"),
    "ioos_gold_1_1": get_filename("tests/data/ioos_1_1.cdl"),
    "kibesillah": get_filename("tests/data/examples/kibesillah.cdl"),
    "l01-met": get_filename("tests/data/examples/l01-met.cdl"),
    "lateral_search_example": get_filename(
        "tests/data/examples/lateral_search_example",
    ),
    "line_geometry": get_filename("tests/data/line_geometry.cdl"),
    "mapping": get_filename("tests/data/mapping.cdl"),
    "multi-dim-coordinates": get_filename("tests/data/multi-dim-coordinates.cdl"),
    "multi-timeseries-orthogonal": get_filename(
        "tests/data/multi-timeseries-orthogonal.cdl",
    ),
    "multi-timeseries-incomplete": get_filename(
        "tests/data/multi-timeseries-incomplete.cdl",
    ),
    "ncei_gold_point_1": get_filename("tests/data/ncei_gold_point_1.cdl"),
    "ncei_gold_point_2": get_filename("tests/data/ncei_gold_point_2.cdl"),
    "NCEI_profile_template_v2_0": get_filename(
        "tests/data/NCEI_profile_template_v2.0_2016-09-22_181835.151325.cdl",
    ),
    "ocos": get_filename("tests/data/examples/ocos.cdl"),
    "ooi_glider": get_filename("tests/data/examples/ooi_glider.cdl"),
    "point": get_filename("tests/data/point.cdl"),
    "polygon_geometry": get_filename("tests/data/polygon_geometry.cdl"),
    "profile-orthogonal": get_filename("tests/data/profile-orthogonal.cdl"),
    "profile-incomplete": get_filename("tests/data/profile-incomplete.cdl"),
    "pr_inundation": get_filename("tests/data/examples/pr_inundation.cdl"),
    "reduced_horizontal_grid": get_filename("tests/data/reduced_horizontal_grid.cdl"),
    "rotated_pole_grid": get_filename("tests/data/rotated_pole_grid.cdl"),
    "rhgrid": get_filename("tests/data/rhgrid.cdl"),
    "rutgers": get_filename("tests/data/ru07-20130824T170228_rt0.cdl"),
    "scalar_coordinate_variable": get_filename(
        "tests/data/scalar_coordinate_variable.cdl",
    ),
    "self-referencing-var": get_filename("tests/data/self-referencing-var.cdl"),
    "self_referencing": get_filename("tests/data/non-comp/self_referencing.cdl"),
    "sldmb_43093_agg": get_filename("tests/data/examples/sldmb_43093_agg.cdl"),
    "string": get_filename("tests/data/string_type_variable.cdl"),
    "swan": get_filename("tests/data/examples/swan.cdl"),
    "sp041": get_filename("tests/data/examples/sp041.cdl"),
    "taxonomy_example": get_filename("tests/data/taxonomy_example.cdl"),
    "timeseries": get_filename("tests/data/timeseries.cdl"),
    "timeseries-profile-single-station": get_filename(
        "tests/data/timeseries-profile-single-station.cdl",
    ),
    "timeseries-profile-multi-station": get_filename(
        "tests/data/timeseries-profile-multi-station.cdl",
    ),
    "timeseries-profile-single-ortho-time": get_filename(
        "tests/data/timeseries-profile-single-ortho-time.cdl",
    ),
    "timeseries-profile-multi-ortho-time": get_filename(
        "tests/data/timeseries-profile-multi-ortho-time.cdl",
    ),
    "timeseries-profile-ortho-depth": get_filename(
        "tests/data/timeseries-profile-ortho-depth.cdl",
    ),
    "timeseries-profile-incomplete": get_filename(
        "tests/data/timeseries-profile-incomplete.cdl",
    ),
    "time_units": get_filename("tests/data/non-comp/time_units.cdl"),
    "trajectory-complete": get_filename("tests/data/trajectory-complete.cdl"),
    "trajectory-implied": get_filename("tests/data/trajectory-implied.cdl"),
    "trajectory-profile-orthogonal": get_filename(
        "tests/data/trajectory-profile-orthogonal.cdl",
    ),
    "trajectory-profile-incomplete": get_filename(
        "tests/data/trajectory-profile-incomplete.cdl",
    ),
    "trajectory": get_filename("tests/data/trajectory.cdl"),
    "trajectory-single": get_filename("tests/data/trajectory-single.cdl"),
    "units_check": get_filename("tests/data/units_check.cdl"),
    "usgs_dem_saipan": get_filename("tests/data/examples/usgs_dem_saipan.cdl"),
    "valid_coordinates": get_filename("tests/data/valid_coordinates.cdl"),
    "vertical_coords": get_filename("tests/data/vertical_coords.cdl"),
    "ww3": get_filename("tests/data/examples/ww3.cdl"),
    "1d_bound_bad": get_filename("tests/data/non-comp/1d_bound_bad.cdl"),
    "2dim": get_filename("tests/data/2dim-grid.cdl"),
    "2d-regular-grid": get_filename("tests/data/2d-regular-grid.cdl"),
    "2d-static-grid": get_filename("tests/data/2d-static-grid.cdl"),
    "3d-regular-grid": get_filename("tests/data/3d-regular-grid.cdl"),
    "3d-static-grid": get_filename("tests/data/3d-static-grid.cdl"),
    "3mf07": get_filename("tests/data/examples/3mf07.cdl"),
}
