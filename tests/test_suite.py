from compliance_checker.runner import ComplianceCheckerCheckSuite
from compliance_checker.base import Result
from netCDF4 import Dataset
from pprint import pprint

def test_suite():
    cs = ComplianceCheckerCheckSuite()
    #ds = cs.load_dataset("/Users/asadeveloper/Downloads/hycomglobalnavy_2012120300.nc", ACDDCheck.beliefs())
    #ds = cs.load_dataset("/Users/asadeveloper/Downloads/hycom.ncml", ACDDCheck.beliefs)

    # @TODO obviously need to update to package data'd datasets
    ds = cs.load_dataset("/Users/asadeveloper/Downloads/sresa1b_ncar_ccsm3_0_run1_200001.nc")
    vals = cs.run(ds, 'acdd')

    pprint(vals)
    assert acdd in vals
    assert vals[acdd][0] == (43.5, 78)

def test_scores():
    cs = CheckSuite()

    v2 = [Result(3, False, 'title'),
          Result(3, False, 'summary'),
          Result(3, False, 'keywords'),
          Result(2, False, 'id'),
          Result(2, False, 'naming_authority'),
          Result(2, False, 'keywords_vocabulary'),
          Result(2, False, 'cdm_data_type'),
          Result(2, True, 'history'),
          Result(2, True, 'comment'),
          Result(2, False, 'date_created'),
          Result(2, False, 'creator_name'),
          Result(2, False, 'creator_url'),
          Result(2, False, 'creator_email'),
          Result(2, True, 'institution'),
          Result(2, False, 'project'),
          Result(2, False, 'processing_level'),
          Result(2, False, 'acknowledgement'),
          Result(2, False, 'geospatial_lat_min'),
          Result(2, False, 'geospatial_lat_max'),
          Result(2, False, 'geospatial_lon_min'),
          Result(2, False, 'geospatial_lon_max'),
          Result(2, False, 'geospatial_vertical_min'),
          Result(2, False, 'geospatial_vertical_max'),
          Result(2, False, 'time_coverage_start'),
          Result(2, False, 'time_coverage_end'),
          Result(2, False, 'time_coverage_duration'),
          Result(2, False, 'time_coverage_resolution'),
          Result(2, False, 'standard_name_vocabulary'),
          Result(2, False, 'license'),
          Result(1, False, 'contributor_name'),
          Result(1, False, 'contributor_role'),
          Result(1, False, 'publisher_name'),
          Result(1, False, 'publisher_url'),
          Result(1, False, 'publisher_email'),
          Result(1, False, 'date_modified'),
          Result(1, False, 'date_issued'),
          Result(1, False, 'geospatial_lat_units'),
          Result(1, False, 'geospatial_lat_resolution'),
          Result(1, False, 'geospatial_lon_units'),
          Result(1, False, 'geospatial_lon_resolution'),
          Result(1, False, 'geospatial_vertical_units'),
          Result(1, False, 'geospatial_vertical_resolution'),
          Result(1, False, 'geospatial_vertical_positive'),
          Result(3, False, ('varattr', 'time', 'var_coverage_content_type')),
          Result(3, False, ('varattr', 'tau', 'var_coverage_content_type')),
          Result(3, False, ('varattr', 'depth', 'var_coverage_content_type')),
          Result(3, False, ('varattr', 'lat', 'var_coverage_content_type')),
          Result(3, False, ('varattr', 'lon', 'var_coverage_content_type')),
          Result(3, False, ('varattr', 'water_temp', 'var_coverage_content_type')),
          Result(3, False, ('varattr', 'salinity', 'var_coverage_content_type')),
          Result(3, False, ('varattr', 'surf_el', 'var_coverage_content_type')),
          Result(3, False, ('varattr', 'water_u', 'var_coverage_content_type')),
          Result(3, False, ('varattr', 'water_v', 'var_coverage_content_type')),
          Result(3, True, ('varattr', 'time', 'var_long_name')),
          Result(3, True, ('varattr', 'tau', 'var_long_name')),
          Result(3, True, ('varattr', 'depth', 'var_long_name')),
          Result(3, True, ('varattr', 'lat', 'var_long_name')),
          Result(3, True, ('varattr', 'lon', 'var_long_name')),
          Result(3, True, ('varattr', 'water_temp', 'var_long_name')),
          Result(3, True, ('varattr', 'salinity', 'var_long_name')),
          Result(3, True, ('varattr', 'surf_el', 'var_long_name')),
          Result(3, True, ('varattr', 'water_u', 'var_long_name')),
          Result(3, True, ('varattr', 'water_v', 'var_long_name')),
          Result(3, False, ('varattr', 'time', 'var_std_name')),
          Result(3, False, ('varattr', 'tau', 'var_std_name')),
          Result(3, True, ('varattr', 'depth', 'var_std_name')),
          Result(3, True, ('varattr', 'lat', 'var_std_name')),
          Result(3, True, ('varattr', 'lon', 'var_std_name')),
          Result(3, True, ('varattr', 'water_temp', 'var_std_name')),
          Result(3, True, ('varattr', 'salinity', 'var_std_name')),
          Result(3, False, ('varattr', 'surf_el', 'var_std_name')),
          Result(3, False, ('varattr', 'water_u', 'var_std_name')),
          Result(3, False, ('varattr', 'water_v', 'var_std_name')),
          Result(3, True, ('varattr', 'time', 'var_units')),
          Result(3, True, ('varattr', 'tau', 'var_units')),
          Result(3, True, ('varattr', 'depth', 'var_units')),
          Result(3, True, ('varattr', 'lat', 'var_units')),
          Result(3, True, ('varattr', 'lon', 'var_units')),
          Result(3, True, ('varattr', 'water_temp', 'var_units')),
          Result(3, True, ('varattr', 'salinity', 'var_units')),
          Result(3, True, ('varattr', 'surf_el', 'var_units')),
          Result(3, True, ('varattr', 'water_u', 'var_units')),
          Result(3, True, ('varattr', 'water_v', 'var_units'))]

    val = cs.scores(v2)
    pprint(val)

    assert val[0] == (42.0, 78)
    assert "varattr" in [v[0] for v in val[1]]

    va = [v for v in val[1] if v[0] == "varattr"]
    assert len(va[0][3]) == 10


