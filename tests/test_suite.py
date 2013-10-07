from compliance_checker.suite import CheckSuite
from compliance_checker.acdd import ACDDCheck
from netCDF4 import Dataset
from pprint import pprint

def test_suite():
    cs = CheckSuite()
    #ds = cs.load_dataset("/Users/asadeveloper/Downloads/hycomglobalnavy_2012120300.nc", ACDDCheck.beliefs())
    #ds = cs.load_dataset("/Users/asadeveloper/Downloads/hycom.ncml", ACDDCheck.beliefs)
    acdd = ACDDCheck()

    vals = cs.run("/Users/asadeveloper/Downloads/sresa1b_ncar_ccsm3_0_run1_200001.nc", acdd)

    pprint(vals)
    assert acdd in vals
    assert vals[acdd][0] == (43.5, 78)

def test_scores():
    cs = CheckSuite()

    v2 = [(3, False, 'title'),
          (3, False, 'summary'),
          (3, False, 'keywords'),
          (2, False, 'id'),
          (2, False, 'naming_authority'),
          (2, False, 'keywords_vocabulary'),
          (2, False, 'cdm_data_type'),
          (2, True, 'history'),
          (2, True, 'comment'),
          (2, False, 'date_created'),
          (2, False, 'creator_name'),
          (2, False, 'creator_url'),
          (2, False, 'creator_email'),
          (2, True, 'institution'),
          (2, False, 'project'),
          (2, False, 'processing_level'),
          (2, False, 'acknowledgement'),
          (2, False, 'geospatial_lat_min'),
          (2, False, 'geospatial_lat_max'),
          (2, False, 'geospatial_lon_min'),
          (2, False, 'geospatial_lon_max'),
          (2, False, 'geospatial_vertical_min'),
          (2, False, 'geospatial_vertical_max'),
          (2, False, 'time_coverage_start'),
          (2, False, 'time_coverage_end'),
          (2, False, 'time_coverage_duration'),
          (2, False, 'time_coverage_resolution'),
          (2, False, 'standard_name_vocabulary'),
          (2, False, 'license'),
          (1, False, 'contributor_name'),
          (1, False, 'contributor_role'),
          (1, False, 'publisher_name'),
          (1, False, 'publisher_url'),
          (1, False, 'publisher_email'),
          (1, False, 'date_modified'),
          (1, False, 'date_issued'),
          (1, False, 'geospatial_lat_units'),
          (1, False, 'geospatial_lat_resolution'),
          (1, False, 'geospatial_lon_units'),
          (1, False, 'geospatial_lon_resolution'),
          (1, False, 'geospatial_vertical_units'),
          (1, False, 'geospatial_vertical_resolution'),
          (1, False, 'geospatial_vertical_positive'),
          (3, False, ('varattr', 'time', 'var_coverage_content_type')),
          (3, False, ('varattr', 'tau', 'var_coverage_content_type')),
          (3, False, ('varattr', 'depth', 'var_coverage_content_type')),
          (3, False, ('varattr', 'lat', 'var_coverage_content_type')),
          (3, False, ('varattr', 'lon', 'var_coverage_content_type')),
          (3, False, ('varattr', 'water_temp', 'var_coverage_content_type')),
          (3, False, ('varattr', 'salinity', 'var_coverage_content_type')),
          (3, False, ('varattr', 'surf_el', 'var_coverage_content_type')),
          (3, False, ('varattr', 'water_u', 'var_coverage_content_type')),
          (3, False, ('varattr', 'water_v', 'var_coverage_content_type')),
          (3, True, ('varattr', 'time', 'var_long_name')),
          (3, True, ('varattr', 'tau', 'var_long_name')),
          (3, True, ('varattr', 'depth', 'var_long_name')),
          (3, True, ('varattr', 'lat', 'var_long_name')),
          (3, True, ('varattr', 'lon', 'var_long_name')),
          (3, True, ('varattr', 'water_temp', 'var_long_name')),
          (3, True, ('varattr', 'salinity', 'var_long_name')),
          (3, True, ('varattr', 'surf_el', 'var_long_name')),
          (3, True, ('varattr', 'water_u', 'var_long_name')),
          (3, True, ('varattr', 'water_v', 'var_long_name')),
          (3, False, ('varattr', 'time', 'var_std_name')),
          (3, False, ('varattr', 'tau', 'var_std_name')),
          (3, True, ('varattr', 'depth', 'var_std_name')),
          (3, True, ('varattr', 'lat', 'var_std_name')),
          (3, True, ('varattr', 'lon', 'var_std_name')),
          (3, True, ('varattr', 'water_temp', 'var_std_name')),
          (3, True, ('varattr', 'salinity', 'var_std_name')),
          (3, False, ('varattr', 'surf_el', 'var_std_name')),
          (3, False, ('varattr', 'water_u', 'var_std_name')),
          (3, False, ('varattr', 'water_v', 'var_std_name')),
          (3, True, ('varattr', 'time', 'var_units')),
          (3, True, ('varattr', 'tau', 'var_units')),
          (3, True, ('varattr', 'depth', 'var_units')),
          (3, True, ('varattr', 'lat', 'var_units')),
          (3, True, ('varattr', 'lon', 'var_units')),
          (3, True, ('varattr', 'water_temp', 'var_units')),
          (3, True, ('varattr', 'salinity', 'var_units')),
          (3, True, ('varattr', 'surf_el', 'var_units')),
          (3, True, ('varattr', 'water_u', 'var_units')),
          (3, True, ('varattr', 'water_v', 'var_units'))]

    val = cs.scores(v2)
    pprint(val)

    assert val[0] == (42.0, 78)
    assert "varattr" in [v[0] for v in val[1]]

    va = [v for v in val[1] if v[0] == "varattr"]
    assert len(va[0][3]) == 10


