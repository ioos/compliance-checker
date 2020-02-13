'''
Check for IOOS-approved attributes
'''
from compliance_checker.base import (BaseCheck, BaseNCCheck, BaseSOSGCCheck,
                                     BaseSOSDSCheck, check_has, Result,
                                     attr_check)
from owslib.namespaces import Namespaces
from lxml.etree import XPath
from compliance_checker.acdd import ACDD1_3Check
from compliance_checker.cfutil import (get_geophysical_variables,
                                       get_instrument_variables)
from compliance_checker import base
from compliance_checker.cf.cf import CF1_6Check, CF1_7Check
import validators
import re


class IOOSBaseCheck(BaseCheck):
    _cc_spec = 'ioos'
    _cc_spec_version = '0.1'
    _cc_description = 'IOOS Inventory Metadata'
    _cc_url = 'https://ioos.github.io/ioos-metadata/ioos-metadata-profile-v1-1.html#ioos-netcdf-metadata-profile-attributes'
    _cc_display_headers = {
        3: 'Highly Recommended',
        2: 'Recommended',
        1: 'Suggested'
    }

    @classmethod
    def _has_attr(cls, ds, attr, concept_name, priority=BaseCheck.HIGH):
        """
        Checks for the existance of attr in ds, with the name/message using concept_name.
        """
        val = cls.std_check(ds, attr)
        msgs = []

        if not val:
            msgs.append("Attr '{}' (IOOS concept: '{}') not found in dataset".format(attr, concept_name))

        return Result(priority, val, concept_name, msgs)

    @classmethod
    def _has_var_attr(cls, dataset, vname, attr, concept_name,
                      priority=BaseCheck.HIGH):
        """
        Checks for the existance of an attr on variable vname in dataset, with the name/message using concept_name.
        """
        val = True
        msgs = []
        if vname not in dataset.variables:
            val = False
            msgs.append("Variable '{}' not present while checking for attr '{}' for IOOS concept: '{}'".format(vname, attr, concept_name))
        else:
            v = dataset.variables[vname]
            if attr not in v.ncattrs():
                val = False
                msgs.append("Attr '{}' not present on var '{}' while checking for IOOS concept: '{}'".format(attr, vname, concept_name))

        return Result(priority, val, concept_name, msgs)


class IOOSNCCheck(BaseNCCheck, IOOSBaseCheck):

    def check_time_period(self, ds):
        """
        Check that time period attributes are both set.
        """
        start = self.std_check(ds, 'time_coverage_start')
        end = self.std_check(ds, 'time_coverage_end')

        msgs = []
        count = 2
        if not start:
            count -= 1
            msgs.append("Attr 'time_coverage_start' is missing")
        if not end:
            count -= 1
            msgs.append("Attr 'time_coverage_end' is missing")

        return Result(BaseCheck.HIGH, (count, 2), 'time coverage start/end', msgs)

    def check_station_location_lat(self, ds):
        """
        Checks station lat attributes are set
        """
        gmin = self.std_check(ds, 'geospatial_lat_min')
        gmax = self.std_check(ds, 'geospatial_lat_max')

        msgs = []
        count = 2
        if not gmin:
            count -= 1
            msgs.append("Attr 'geospatial_lat_min' is missing")
        if not gmax:
            count -= 1
            msgs.append("Attr 'geospatial_lat_max' is missing")

        return Result(BaseCheck.HIGH, (count, 2), 'geospatial lat min/max', msgs)

    def check_station_location_lon(self, ds):
        """
        Checks station lon attributes are set
        """
        gmin = self.std_check(ds, 'geospatial_lon_min')
        gmax = self.std_check(ds, 'geospatial_lon_max')

        msgs = []
        count = 2
        if not gmin:
            count -= 1
            msgs.append("Attr 'geospatial_lon_min' is missing")
        if not gmax:
            count -= 1
            msgs.append("Attr 'geospatial_lon_max' is missing")

        return Result(BaseCheck.HIGH, (count, 2), 'geospatial lon min/max', msgs)


class IOOS0_1Check(IOOSNCCheck):
    _cc_spec_version = '0.1'
    _cc_description = 'IOOS Inventory Metadata'
    register_checker = True

    def check_global_attributes(self, ds):
        """
        Check all global NC attributes for existence.

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        return [
            self._has_attr(ds, 'acknowledgement', 'Platform Sponsor'),
            self._has_attr(ds, 'publisher_email', 'Station Publisher Email'),
            self._has_attr(ds, 'publisher_email', 'Service Contact Email', BaseCheck.MEDIUM),
            self._has_attr(ds, 'institution', 'Service Provider Name', BaseCheck.MEDIUM),
            self._has_attr(ds, 'publisher_name', 'Service Contact Name', BaseCheck.MEDIUM),
            self._has_attr(ds, 'Conventions', 'Data Format Template Version', BaseCheck.MEDIUM),
            self._has_attr(ds, 'publisher_name', 'Station Publisher Name', BaseCheck.HIGH),
        ]

    def check_variable_attributes(self, ds):
        """
        Check IOOS concepts that come from NC variable attributes.

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        return [
            self._has_var_attr(ds, 'platform', 'long_name', 'Station Long Name'),
            self._has_var_attr(ds, 'platform', 'short_name', 'Station Short Name'),
            self._has_var_attr(ds, 'platform', 'source', 'Platform Type'),
            self._has_var_attr(ds, 'platform', 'ioos_name', 'Station ID'),
            self._has_var_attr(ds, 'platform', 'wmo_id', 'Station WMO ID'),
            self._has_var_attr(ds, 'platform', 'comment', 'Station Description'),
        ]

    def check_variable_names(self, ds):
        """
        Ensures all variables have a standard_name set.
        """
        msgs = []
        count = 0

        for k, v in ds.variables.items():
            if 'standard_name' in v.ncattrs():
                count += 1
            else:
                msgs.append("Variable '{}' missing standard_name attr".format(k))

        return Result(BaseCheck.MEDIUM, (count, len(ds.variables)), 'Variable Names', msgs)

    def check_altitude_units(self, ds):
        """
        If there's a variable named z, it must have units.

        @TODO: this is duplicated with check_variable_units
        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        if 'z' in ds.variables:
            msgs = []
            val = 'units' in ds.variables['z'].ncattrs()
            if not val:
                msgs.append("Variable 'z' has no units attr")
            return Result(BaseCheck.LOW, val, 'Altitude Units', msgs)

        return Result(BaseCheck.LOW, (0, 0), 'Altitude Units', ["Dataset has no 'z' variable"])

    def check_variable_units(self, ds):
        """
        Ensures all variables have units.
        """
        msgs = []
        count = 0

        for k, v in ds.variables.items():
            if 'units' in v.ncattrs():
                count += 1
            else:
                msgs.append("Variable '{}' missing units attr".format(k))

        return Result(BaseCheck.MEDIUM, (count, len(ds.variables)), 'Variable Units', msgs)


class IOOS1_1Check(IOOSNCCheck):
    '''
    Compliance checker implementation of IOOS Metadata Profile, Version 1.1

    Related links:
    https://ioos.github.io/ioos-metadata/ioos-metadata-profile-v1-1.html#ioos-netcdf-metadata-profile-attributes
    https://github.com/ioos/compliance-checker/issues/69
    https://github.com/ioos/compliance-checker/issues/358
    '''
    _cc_spec_version = '1.1'
    _cc_description = 'IOOS Metadata Profile, Version 1.1'
    _cc_url = 'https://ioos.github.io/ioos-metadata/ioos-metadata-profile-v1-1.html#ioos-netcdf-metadata-profile-attributes'
    register_checker = True

    def __init__(self):
        # Define the global attributes
        self.required_atts = [
            'contributor_name',
            'contributor_role',
            'creator_country',
            'creator_email',
            'creator_sector',
            'featureType',
            'id',
            'institution',
            'naming_authority',
            'platform',
            'platform_vocabulary',
            'publisher_country',
            'publisher_email',
            'publisher_name',
            'standard_name_vocabulary',
            'title'
        ]

        self.rec_atts = [
            'creator_address',
            'creator_city',
            'creator_name',
            'creator_phone',
            'creator_state',
            'creator_url',
            'creator_zipcode',
            'keywords',
            'license',
            'publisher_address',
            'publisher_city',
            'publisher_phone',
            'publisher_state',
            'publisher_url',
            'publisher_zipcode',
            'summary'
        ]

    @check_has(BaseCheck.HIGH)
    def check_high(self, ds):
        '''
        Performs a check on each highly recommended attributes' existence in the dataset

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        return self.required_atts

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        '''
        Performs a check on each recommended attributes' existence in the dataset

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        return self.rec_atts

    def check_platform_variables(self, ds):
        '''
        The value of platform attribute should be set to another variable which
        contains the details of the platform. There can be multiple platforms
        involved depending on if all the instances of the featureType in the
        collection share the same platform or not. If multiple platforms are
        involved, a variable should be defined for each platform and referenced
        from the geophysical variable in a space separated string.

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        platform_names = getattr(ds, 'platform', '').split(' ')
        val = all(platform_name in ds.variables for platform_name in platform_names)
        msgs = []
        if not val:
            msgs = [('The value of "platform" global attribute should be set to another variable '
                     'which contains the details of the platform. If multiple platforms are '
                     'involved, a variable should be defined for each platform and referenced '
                     'from the geophysical variable in a space separated string.')]
        return [Result(BaseCheck.HIGH, val, 'platform variables', msgs)]

    def check_platform_variable_attributes(self, ds):
        '''
        Platform variables must contain the following attributes:
            ioos_code
            long_name
            short_name
            type

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        results = []
        platform_name = getattr(ds, 'platform', '')
        # There can be multiple platforms defined here (space separated)
        for platform in platform_name.split(' '):
            if platform in ds.variables:
                results += [
                    self._has_var_attr(ds, platform, 'long_name', 'Platform Long Name'),
                    self._has_var_attr(ds, platform, 'short_name', 'Platform Short Name'),
                    self._has_var_attr(ds, platform, 'ioos_code', 'Platform IOOS Code'),
                    self._has_var_attr(ds, platform, 'type', 'Platform Type')
                ]
        return results

    def check_geophysical_vars_fill_value(self, ds):
        '''
        Check that geophysical variables contain fill values.

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        results = []
        for geo_var in get_geophysical_variables(ds):
            results.append(
                self._has_var_attr(ds, geo_var, '_FillValue', '_FillValue', BaseCheck.MEDIUM),
            )
        return results

    def check_geophysical_vars_standard_name(self, ds):
        '''
        Check that geophysical variables contain standard names.

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        results = []
        for geo_var in get_geophysical_variables(ds):
            results.append(
                self._has_var_attr(ds, geo_var, 'standard_name', 'geophysical variables standard_name'),
            )
        return results

    def check_units(self, ds):
        '''
        Required for most all variables that represent dimensional quantities.
        The value should come from udunits authoritative vocabulary, which is
        documented in the CF standard name table with it's corresponding
        standard name.

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        cf16 = CF1_6Check()
        return cf16.check_units(ds)

class IOOS1_2_ConventionsValidator(base.RegexValidator):
    validator_regex = r"\bIOOS-1.2\b"
    validator_fail_message = "{} must contain the string \"IOOS 1.2\""

class IOOS1_2Check(IOOSNCCheck):
    """
    Class to implement the IOOS Metadata 1.2 Specification
    """

    _cc_spec_version = '1.2'
    _cc_description = 'IOOS Metadata Profile, Version 1.2'
    _cc_url = 'https://ioos.github.io/ioos-metadata/ioos-metadata-profile-v1-2.html'
    register_checker = True

    def __init__(self):

        # instantiate objects used for delegation
        self.acdd1_6 = ACDD1_3Check()
        self.cf1_7   = CF1_7Check()

        # extend standard_names set to include QARTOD standard_names
        self._qartod_std_names = [
            "aggregate_quality_flag",
            "attenuated_signal_test_quality_flag",
            "climatology_test_quality_flag",
            "flat_line_test_quality_flag",
            "gap_test_quality_flag",
            "gross_range_test_quality_flag",
            "location_test_quality_flag",
            "multi_variate_test_quality_flag",
            "neighbor_test_quality_flag",
            "rate_of_change_test_quality_flag",
            "spike_test_quality_flag",
            "syntax_test_quality_flag"
        ]
        self.cf1_7._std_names._names.extend(self._qartod_std_names)

        # check geophysical variables have the following attrs:
        self.check_var_attrs = (
            ("_FillValue", BaseCheck.MEDIUM),
            ("missing_value", BaseCheck.MEDIUM),
            #( "standard_name", BaseCheck.HIGH # already checked in CF1_7Check.check_standard_name()
            ("standard_name_uri", BaseCheck.MEDIUM),
            #( "units", BaseCheck.HIGH # already checked in CF1_7Check.check_units()
            ("platform", BaseCheck.HIGH),
            #( "wmo_platform_code", BaseCheck.HIGH # only "if applicable", see check_wmo_platform_code()
        )

        self.required_atts = [
            ('Conventions', IOOS1_2_ConventionsValidator),
            'creator_country',
            ('creator_email', base.EmailValidator),
            'creator_institution',
            ('creator_sector', {"gov_state", "nonprofit", "tribal", "other",
                                "unknown", "gov_municipal", "industry",
                                "gov_federal", "academic"}),
            ('creator_url', base.UrlValidator),
            'featureType',
            'id',
            ('infoUrl', base.UrlValidator),
            'license',
            'naming_authority',
            'platform',
            'platform_name',
            'publisher_country',
            ('publisher_email', base.EmailValidator),
            'publisher_institution',
            'publisher_url',
            # TODO: handle standard name table exclusion for v38?
            ('standard_name_vocabulary',
             re.compile(r'^CF Standard Name Table v[1-9]\d*$')),
            'summary',
            'title'
        ]

        self.rec_atts = [
            ('contributor_email', base.EmailValidator),
            'contributor_name',
            'contributor_role',
            'contributor_role_vocabulary',
            ('contributor_url', base.UrlValidator),
            'creator_address',
            'creator_city',
            'creator_name',
            'creator_phone',
            'creator_postalcode',
            'creator_state',
            # checked in check_creator_and_publisher_type
            #'creator_type',
            'institution',
            'instrument',
            'ioos_ingest',
            'keywords',
            'platform_id',
            'publisher_address',
            'publisher_city',
            'publisher_name',
            'publisher_phone',
            'publisher_postalcode',
            'publisher_state',
            # checked in check_creator_and_publisher_type
            #'publisher_type',
            ('references', base.UrlValidator)
        ]

    def setup(self, ds):
        self.platform_vars = self._find_platform_vars(ds)

    def _find_platform_vars(self, ds):
        """
        Finds any variables referenced by 'platform' attribute which exist in
        the dataset.

        Parameters
        ----------
        ds: netCDF4.Dataset
            An open netCDF4 Dataset.

        Returns
        -------
        set of netCDF4.Variable
            Set of variables which are platform variables.
        """
        plat_vars = ds.get_variables_by_attributes(platform=
                                                   lambda p: isinstance(p, str))
        return {ds.variables[var.platform] for var in plat_vars if
                var.platform in ds.variables}

    @check_has(BaseCheck.HIGH)
    def check_high(self, ds):
        '''
        Performs a check on each highly recommended attributes' existence in the dataset

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        return self.required_atts

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        '''
        Performs a check on each recommended attributes' existence in the dataset

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        return self.rec_atts

    def check_standard_name(self, ds):
        """
        Wrapper for checking standard names using the CF module.
        Extends the StandardNameTable to include QARTOD variable
        standard names.
        """

        return self.cf1_7.check_standard_name(ds)

    def check_units(self, ds):
        """
        Wrapper to check units with the CF module.
        """

        return self.cf1_7.check_units(ds)

    def check_vars_have_attrs(self, ds):
        """
        Using the tuples defined in __init__, check that each variable has
        the attributes using the corresponding priority.

        :param netCDF4.Dataset: open netCDF4 dataset
        """

        results = []
        # NOTE should it also find 'geospatial` variables?
        for geo_var in get_geophysical_variables(ds):
            for attr_tuple in self.check_var_attrs:
                results.append(
                    self._has_var_attr(
                        ds,
                        geo_var,
                        attr_tuple[0], # attribute name
                        attr_tuple[0], # attribute name used as 'concept_name'
                        attr_tuple[1]  # priority level
                    )
                )

        return results


    def check_platform_variable_cf_role(self, ds):
        """
        Verify that any platform variables have valid CF roles

        Args:
            ds (netCDF-4 Dataset): open Dataset object

        Returns:
            list of Result objects
        """
        valid_cf_roles = {"timeseries_id", "profile_id", "trajectory_id"}
        cf_role_results = []
        for var in self.platform_vars:
            attr_check(("cf_role", valid_cf_roles), ds, BaseCheck.HIGH,
                        cf_role_results, var_name=var.name)
        return cf_role_results

    def check_creator_and_publisher_type(self, ds):
        """
        Check if global attribute creator_type and publisher_type
        are contained within the values "person", "group", "institution", or
        "position".  If creator_type is not present within the global
        attributes, assume it is set to a value of "person".

        Parameters
        ----------
        ds: netCDF4.Dataset
            An open netCDF4 Dataset

        Returns
        -------
        list of Result
        """
        result_list = []
        for global_att_name in ("creator_type", "publisher_type"):
            messages = []
            try:
                att_value = ds.getncattr(global_att_name)
            except AttributeError:
                # if the attribute isn't found, it's implicitly assigned
                # a value of "person", so it automatically passes.
                pass_stat = True
            else:
                expected_types = {"person", "group", "institution", "position"}
                if att_value in expected_types:
                    pass_stat = True
                else:
                    pass_stat = False
                    messages.append("If specified, {} must be in value list "
                                    "({})".format(global_att_name,
                                                  sorted(expected_types)))

            result_list.append(Result(BaseCheck.MEDIUM, pass_stat,
                                      global_att_name, messages))

        return result_list

    def check_single_platform(self, ds):
        """
        Verify that a dataset only has a single platform attribute. If one exists,
        examine the featureType of the dataset. If the featureType is
        [point, timeSeries, profile, trajectory] and cf_role in [timeseries_id,
        profile_id, trajectory_id] dimensionality of the variable containing
        cf_role must be 1 as "we only want a single glider/auv/ship"; if
        featureType in [timeseries_id, trajectory_id], dimensionality of the
        variable must also be one as  "we only want a single timeSeries aka buoy".
        If cf_role==profile_id, it can have whatever dimension.

        Gridded model datasets are not required to declare a platform
        or platform variables.

        Args:
            ds (netCDF-4 Dataset): open Dataset object

        Returns:
            Result
        """

        results = []
        glb_platform = getattr(ds, "platform", None)

        platform_set = set()
        for v in ds.get_variables_by_attributes(platform=lambda x: x is not None):
            platform_set.add(v.getncattr("platform"))

        num_platforms = len(platform_set)
        if num_platforms > 1 and glb_platform:
            msg = "A dataset may only have one platform; {} found".format(len(platform_set))
            val = False
            results.append(Result(BaseCheck.HIGH, val, "platform", [msg]))


        elif ((not glb_platform) and num_platforms > 0):
            msg = "If platform variables exist, a global attribute \"platform\" must also exist"
            val = False
            results.append(Result(BaseCheck.HIGH, val, "platform", [msg]))

        elif num_platforms == 0 and glb_platform:
            msg = "A dataset with a global \"platform\" attribute must platform have variables"
            val = False
            results.append(Result(BaseCheck.HIGH, val, "platform", [msg]))

        elif num_platforms == 0 and (not glb_platform):
            msg = "Gridded model datasets are not required to declare a platform"
            val = True
            results.append(Result(BaseCheck.HIGH, val, "platform", [msg]))

        else: # num_platforms==1 and glb_platform, test the dimensionality

            num_plat_val = True if num_platforms == 1 else False

            feature_type = getattr(ds, "featureType", "").lower()
            if not feature_type:
                return results

            # filter out cf_role exists
            cf_role_vars = ds.get_variables_by_attributes(cf_role=lambda x: x is not None)
            print(cf_role_vars)
            num_cf_role_vars = len(cf_role_vars)
            msg = "With a single platform provided, the dimension of the cf_role " +\
                  "variable {cf_role_var} (cf_role=={cf_role}) should also " +\
                  "be equal to 1 (it is {dim})"

            for var in cf_role_vars:
                cf_role = getattr(var, "cf_role")
                shp = var.shape[0] if len(var.shape) > 0 else 1
                if (
                       feature_type in ["point", "timeseries", "profile", "trajectory", "timeseriesprofile", "trajectoryprofile"]
                       and
                       cf_role in ["timeseries_id", "profile_id", "trajectory_id"]
                   ):
                    print(var)
                    if (num_cf_role_vars==1) or (num_cf_role_vars>1 and cf_role!="profile_id"):
                        # shape must be 1 (or if no length, that's okay too)
                        _val = shp==1
                    elif (num_cf_role_vars>1) and (cf_role=="profile_id"):
                        # can have any dimension if there are more than one cf_role?
                        _val = True

                    results.append(
                       Result(
                           BaseCheck.HIGH,
                           _val,
                           "platform variables",
                           [msg.format(cf_role_var=var.name, cf_role=cf_role, dim=shp)]
                       )
                    )

        return results

    def check_platform_vocabulary(self, ds):
        """
        The platform_vocabulary attribute is recommended to be a URL to
        http://mmisw.org/ont/ioos/platform or
        http://vocab.nerc.ac.uk/collection/L06/current/. However,
        it is required to at least be a URL.

        Args:
            ds (netCDF4.Dataset): open Dataset

        Returns:
            Result
        """

        m = "platform_vocabulary must be a valid URL"
        pvocab = getattr(ds, "platform_vocabulary", None)
        val = bool(validators.url(getattr(v, "references", "")))
        return Result(BaseCheck.MEDIUM, val, "platform_vocabulary", [m])

    def _check_var_gts_ingest(self, ds, var, var_ingest_val, do_ingest, msg):
        """
        Helper function for check_gts_ingest(). Check that a given variable
          - has a valid CF standard name
          - has a QARTOD aggregates variable
          - has valid units

        Args:
            attr (?): attribute value

        Returns:
            Result
        """

        val = False
        if ingest == "true":

            # should have an ancillary variable with standard_name aggregate_quality_flag
            avar_val = False
            anc_vars = getattr(var, "ancillary_variables", "").split(" ")
            for av in anc_vars:
                if av in ds.variables:
                    if getattr(ds.variables[av], "standard_name", "") == "aggregate_quality_flag":
                        avar_val = True
                        break

            # if variable is flagged for ingest, but no global attr present, error
            val = True if (avar_val and do_ingest) else False

        elif ingest == "false":
            # return positive result
            val = True

        return Result(BaseCheck.HIGH, val, "gts_ingest variable", [msg])

    def check_gts_ingest(self, ds):
        """
        Check if a dataset has a global gts_ingest attribute and if any
        variables also have the gts_ingest attribute.

        According to https://ioos.github.io/ioos-metadata/ioos-metadata-profile-v1-2.html#requirements-for-ioos-dataset-gts-ingest,
        the gts_ingest is "required, if applicable". Because the Compliance Checker
        cannot accurately guess applicability of a certain dataset or variable,
        this check simply verifies that if it exists, the value is a string
        denoting "true" or "false".

        Any variables which a user would like ingested must also contain the
        gts_ingest attribute with a value of true. The variable must:
          - have a valid CF standard_name attribute (already checked)
          - have an ancillary variable reqpresenting QARTOD aggregate flags
          - have a valid udunits units attribute

        Args:
            ds (netCDF4.Dataset): open Dataset

        Returns:
            list of Result objects
        """

        default_pass_result = Result(BaseCheck.HIGH, True, "gts_ingest", ["gts_ingest"])

        results = []

        # check global
        glb_msg = ("If provided, the global attribute \"gts_ingest\" must be a "
                   "string and its value must be one of \"true\" or \"false\"")
        glb_gts_attr = getattr(ds, "gts_ingest", None)
        if glb_gts_attr and (glb_gts_attr=="true" or glb_gts_ingest=="false"):
            do_gts = True if glb_gts_ingest == "true" else False
            results.append(Result(BaseCheck.HIGH, True, "gts_ingest", [glb_msg]))
        else:
            do_gts = False
            results.append(default_pass_result)

        # check variables
        var_msg = ("The attribute \"gts_ingest\" of variable \"{v}\" "
                   " must fulfill the following:\n"
                   "  - must be a string with value \"true\" or \"false\";\n"
                   "  - have a valid CF standard_name attribute (already checked);\n"
                   "  - have an ancillary variable reqpresenting QARTOD aggregate flags;\n"
                   "  - have a valid udunits units attribute\n"
                   "The global attribute \"gts_ingest\" "
                   "must also have a value of \"true\".")

        for v in ds.variables:
            _attr = getattr(ds.variables[v], "gts_ingest", None)
            if _attr:
                results.append(self._check_var_gts_ingest(ds, var, _attr, do_ingest, var_msg.format(v=v)))
            else:
                results.append(default_pass_result)

        return results

    def check_instrument_variables(self, ds):
        """
        If present, the instrument_variable is one that contains additional
        metadata about the instrument the data was collected with.

        Args:
            ds (netCDF4.Dataset): open Dataset

        Returns:
            list of Results
        """

        results = []
        instr_vars = get_instrument_variables(ds)

        # check for component, disciminant
        for instr in instr_vars:
            if instr in ds.variables:
                compnt = getattr(ds.variables[instr], "component", None)
                m = ["component attribute of {} ({}) must be a string".format(instr, compnt)]
                if compnt:
                    results.append(
                        Result(BaseCheck.MEDIUM, isinstance(compnt, str), "instrument_variable", m)
                    )
                else:
                    results.append(Result(BaseCheck.MEDIUM, True, "instrument_variable", m))

                disct = getattr(ds.variables[instr], "discriminant", None)
                m = ["discriminant attribute of {} ({}) must be a string".format(instr, disct)]
                if disct:
                    results.append(
                        Result(BaseCheck.MEDIUM, isinstance(disct, str), "instrument_variable", m)
                    )
                else:
                    results.append(Result(BaseCheck.MEDIUM, True, "instrument_variable", m))

        return results

    def check_qartod_variables_references(self, ds):
        """
        For any variables that are deemed QARTOD variables, check that they
        contain the "references" attribute and that the value of the attribute
        is a valid URL.

        Args:
            ds (netCDF4.Dataset): open Dataset

        Returns:
            list of Results
        """

        results = []
        ctxt = "qartod_variable:references"
        for v in ds.get_variables_by_attributes(standard_name=lambda x: x in self._qartod_std_names):
            msg = "\"references\" attribute for variable \"{}\" must be a valid URL".format(v.name)
            val = bool(validators.url(getattr(v, "references", "")))
            results.append(Result(BaseCheck.MEDIUM, val, ctxt, [msg]))

        return results

    def check_wmo_platform_code(self, ds):
        """
        If a WMO Platform Code is given as a global variable, check that it is
        well-formed. According to the WMO, valid codes are a numeric string
        comprised of 5 or 7 characters.

        Args:
            ds (netCDF4.Dataset): open Dataset

        Returns:
            Result
        """

        valid = True
        ctxt = "wmo_platform_code"
        msg = "The wmo_platform_code must be a numeric string of 5 or 7 characters"

        code = getattr(ds, "wmo_platform_code", None)
        if code:
           if not (
               isinstance(code, str) and
               code.isnumeric() and
               (5 <= len(code) <= 7)
           ):
               valid = False

        return Result(BaseCheck.HIGH, valid, ctxt, [msg])

class IOOSBaseSOSCheck(BaseCheck):
    _cc_spec = 'ioos_sos'
    _cc_spec_version = '0.1'
    _cc_description = ('IOOS Inventory Metadata checks for the Sensor Observation System (SOS). '
                       'Checks SOS functions GetCapabilities and DescribeSensor.')
    register_checker = True
    # requires login
    _cc_url = 'http://sdf.ndbc.noaa.gov/sos/'


class IOOSSOSGCCheck(BaseSOSGCCheck, IOOSBaseSOSCheck):

    # set up namespaces for XPath
    ns = Namespaces().get_namespaces(['sos', 'gml', 'xlink'])
    ns['ows'] = Namespaces().get_namespace('ows110')

    @check_has(BaseCheck.HIGH)
    def check_high(self, ds):
        return [

        ]

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        return [
            ('service_contact_email', XPath("/sos:Capabilities/ows:ServiceProvider/ows:ServiceContact/ows:ContactInfo/ows:Address/ows:ElectronicMailAddress", namespaces=self.ns)),
            ('service_contact_name', XPath("/sos:Capabilities/ows:ServiceProvider/ows:ServiceContact/ows:IndividualName", namespaces=self.ns)),
            ('service_provider_name', XPath("/sos:Capabilities/ows:ServiceProvider/ows:ProviderName", namespaces=self.ns)),

            ('service_title', XPath("/sos:Capabilities/ows:ServiceProvider/ows:ProviderName", namespaces=self.ns)),
            ('service_type_name', XPath("/sos:Capabilities/ows:ServiceIdentification/ows:ServiceType", namespaces=self.ns)),
            ('service_type_version', XPath("/sos:Capabilities/ows:ServiceIdentification/ows:ServiceTypeVersion", namespaces=self.ns)),
            # ds.identification[0].observed_properties has this as well, but
            # don't want to try to shoehorn a function here
            # ('variable_names', len(ds.identification[0].observed_properties) > 0)
            ('variable_names', XPath("/sos:Capabilities/sos:Contents/sos:ObservationOfferingList/sos:ObservationOffering/sos:observedProperty",
             namespaces=self.ns)),
            ('data_format_template_version', XPath("/sos:Capabilities/ows:OperationsMetadata/ows:ExtendedCapabilities/gml:metaDataProperty[@xlink:title='ioosTemplateVersion']/gml:version",
             namespaces=self.ns))
        ]

    @check_has(BaseCheck.LOW)
    def check_suggested(self, ds):
        return [
            'altitude_units'
        ]


class IOOSSOSDSCheck(BaseSOSDSCheck, IOOSBaseSOSCheck):

    # set up namespaces for XPath
    ns = Namespaces().get_namespaces(['sml', 'swe', 'gml', 'xlink'])

    @check_has(BaseCheck.HIGH)
    def check_high(self, ds):
        return [
            ('platform_sponsor', XPath("/sml:SensorML/sml:member/sml:System/sml:classification/sml:ClassifierList/sml:classifier[@name='sponsor']/sml:Term/sml:value", namespaces=self.ns)),
            ('platform_type', XPath("/sml:SensorML/sml:member/sml:System/sml:classification/sml:ClassifierList/sml:classifier[@name='platformType']/sml:Term/sml:value", namespaces=self.ns)),
            ('station_publisher_name', XPath("/sml:SensorML/sml:member/sml:System/sml:contact/sml:ContactList/sml:member[@xlink:role='http://mmisw.org/ont/ioos/definition/publisher']/sml:ResponsibleParty/sml:organizationName", namespaces=self.ns)),
            ('station_publisher_email', XPath("/sml:SensorML/sml:member/sml:System/sml:contact/sml:ContactList/sml:member[@xlink:role='http://mmisw.org/ont/ioos/definition/publisher']/sml:ResponsibleParty/sml:contactInfo/address/sml:electronicMailAddress", namespaces=self.ns)),
            ('station_id', XPath("/sml:SensorML/sml:member/sml:System/sml:identification/sml:IdentifierList/sml:identifier[@name='stationID']/sml:Term/sml:value", namespaces=self.ns)),
            ('station_long_name', XPath("/sml:SensorML/sml:member/sml:System/sml:identification/sml:IdentifierList/sml:identifier[@name='longName']/sml:Term/sml:value", namespaces=self.ns)),
            ('station_short_name', XPath("/sml:SensorML/sml:member/sml:System/sml:identification/sml:IdentifierList/sml:identifier[@name='shortName']/sml:Term/sml:value", namespaces=self.ns)),
            ('station_wmo_id', XPath("/sml:SensorML/sml:member/sml:System/sml:identification/sml:IdentifierList/sml:identifier/sml:Term[@definition=\"http://mmisw.org/ont/ioos/definition/wmoID\"]/sml:value", namespaces=self.ns)),
            ('time_period', XPath("/sml:SensorML/sml:member/sml:System/sml:capabilities[@name='observationTimeRange']/swe:DataRecord/swe:field[@name='observationTimeRange']/swe:TimeRange/swe:value", namespaces=self.ns)),
            ('operator_email', XPath("/sml:SensorML/sml:member/sml:System/sml:contact/sml:ContactList/sml:member[@xlink:role='http://mmisw.org/ont/ioos/definition/operator']/sml:ResponsibleParty/sml:contactInfo/address/sml:electronicMailAddress", namespaces=self.ns)),
            ('operator_name', XPath("/sml:SensorML/sml:member/sml:System/sml:contact/sml:ContactList/sml:member[@xlink:role='http://mmisw.org/ont/ioos/definition/operator']/sml:ResponsibleParty/sml:organizationName", namespaces=self.ns)),
            ('station_description', XPath("/sml:SensorML/sml:member/sml:System/gml:description", namespaces=self.ns)),
            # replaced with lon/lat with point
            ('station_location_point', XPath("/sml:SensorML/sml:member/sml:System/sml:location/gml:Point/gml:pos", namespaces=self.ns))
        ]

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        return [
            ('sensor_descriptions', XPath("/sml:SensorML/sml:member/sml:System/sml:components/sml:ComponentList/sml:component/sml:System/gml:description", namespaces=self.ns)),
            ('sensor_ids', XPath("/sml:SensorML/sml:member/sml:System/sml:components/sml:ComponentList/sml:component/sml:System/@gml:id", namespaces=self.ns)),
            ('sensor_names', XPath("/sml:SensorML/sml:member/sml:System/sml:components/sml:ComponentList/sml:component/@name", namespaces=self.ns)),

            ('data_format_template_version', XPath("/sml:SensorML/sml:capabilities/swe:SimpleDataRecord/swe:field[@name='ioosTemplateVersion']/swe:Text/swe:value", namespaces=self.ns)),

            ('variable_names', XPath("/sml:SensorML/sml:member/sml:System/sml:components/sml:ComponentList/sml:component/sml:System/sml:outputs/sml:OutputList/sml:output/swe:Quantity/@definition", namespaces=self.ns)),
            ('variable_units', XPath("/sml:SensorML/sml:member/sml:System/sml:components/sml:ComponentList/sml:component/sml:System/sml:outputs/sml:OutputList/sml:output/swe:Quantity/swe:uom/@code", namespaces=self.ns)),
            ('network_id', XPath("/sml:SensorML/sml:member/sml:System/sml:capabilities[@name='networkProcedures']/swe:SimpleDataRecord/gml:metaDataProperty/@xlink:href", namespaces=self.ns)),
            ('operator_sector', XPath("/sml:SensorML/sml:member/sml:System/sml:classification/sml:ClassifierList/sml:classifier[@name='operatorSector']/sml:Term/sml:value", namespaces=self.ns)),
        ]

    @check_has(BaseCheck.LOW)
    def check_suggested(self, ds):
        return [
        ]
