"""
Check for IOOS-approved attributes
"""

import re
from numbers import Number

import validators
from cf_units import Unit
from lxml.etree import XPath
from owslib.namespaces import Namespaces

import compliance_checker.cf.util as cfutil
from compliance_checker import base
from compliance_checker.acdd import ACDD1_3Check
from compliance_checker.base import (
    BaseCheck,
    BaseNCCheck,
    BaseSOSDSCheck,
    BaseSOSGCCheck,
    Result,
    TestCtx,
    check_has,
)
from compliance_checker.cf.cf import CF1_6Check, CF1_7Check


class IOOSBaseCheck(BaseCheck):
    _cc_spec = "ioos"
    _cc_spec_version = "0.1"
    _cc_description = "IOOS Inventory Metadata"
    _cc_url = "https://ioos.github.io/ioos-metadata/ioos-metadata-profile-v1-1.html#ioos-netcdf-metadata-profile-attributes"
    _cc_display_headers = {3: "Highly Recommended", 2: "Recommended", 1: "Suggested"}

    @classmethod
    def _has_attr(cls, ds, attr, concept_name, priority=BaseCheck.HIGH):
        """
        Checks for the existence of attr in ds, with the name/message using concept_name.
        """
        val = cls.std_check(ds, attr)
        msgs = []

        if not val:
            msgs.append(
                f"Attr '{attr}' (IOOS concept: '{concept_name}') not found in dataset",
            )

        return Result(priority, val, concept_name, msgs)

    @classmethod
    def _has_var_attr(cls, dataset, vname, attr, concept_name, priority=BaseCheck.HIGH):
        """
        Checks for the existence of an attr on variable vname in dataset, with the name/message using concept_name.
        """
        val = True
        msgs = []
        if vname not in dataset.variables:
            val = False
            msgs.append(
                f"Variable '{vname}' not present while checking for attr '{attr}' for IOOS concept: '{concept_name}'",
            )
        else:
            v = dataset.variables[vname]
            if attr not in v.ncattrs():
                val = False
                msgs.append(
                    f"Attr '{attr}' not present on var '{vname}' while checking for IOOS concept: '{concept_name}'",
                )

        return Result(priority, val, concept_name, msgs)


class IOOSNCCheck(BaseNCCheck, IOOSBaseCheck):
    def check_time_period(self, ds):
        """
        Check that time period attributes are both set.
        """
        start = self.std_check(ds, "time_coverage_start")
        end = self.std_check(ds, "time_coverage_end")

        msgs = []
        count = 2
        if not start:
            count -= 1
            msgs.append("Attr 'time_coverage_start' is missing")
        if not end:
            count -= 1
            msgs.append("Attr 'time_coverage_end' is missing")

        return Result(BaseCheck.HIGH, (count, 2), "time coverage start/end", msgs)

    def check_station_location_lat(self, ds):
        """
        Checks station lat attributes are set
        """
        gmin = self.std_check(ds, "geospatial_lat_min")
        gmax = self.std_check(ds, "geospatial_lat_max")

        msgs = []
        count = 2
        if not gmin:
            count -= 1
            msgs.append("Attr 'geospatial_lat_min' is missing")
        if not gmax:
            count -= 1
            msgs.append("Attr 'geospatial_lat_max' is missing")

        return Result(BaseCheck.HIGH, (count, 2), "geospatial lat min/max", msgs)

    def check_station_location_lon(self, ds):
        """
        Checks station lon attributes are set
        """
        gmin = self.std_check(ds, "geospatial_lon_min")
        gmax = self.std_check(ds, "geospatial_lon_max")

        msgs = []
        count = 2
        if not gmin:
            count -= 1
            msgs.append("Attr 'geospatial_lon_min' is missing")
        if not gmax:
            count -= 1
            msgs.append("Attr 'geospatial_lon_max' is missing")

        return Result(BaseCheck.HIGH, (count, 2), "geospatial lon min/max", msgs)


class IOOS0_1Check(IOOSNCCheck):
    _cc_spec_version = "0.1"
    _cc_description = "IOOS Inventory Metadata"
    register_checker = True

    def check_global_attributes(self, ds):
        """
        Check all global NC attributes for existence.

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        return [
            self._has_attr(ds, "acknowledgement", "Platform Sponsor"),
            self._has_attr(ds, "publisher_email", "Station Publisher Email"),
            self._has_attr(
                ds,
                "publisher_email",
                "Service Contact Email",
                BaseCheck.MEDIUM,
            ),
            self._has_attr(
                ds,
                "institution",
                "Service Provider Name",
                BaseCheck.MEDIUM,
            ),
            self._has_attr(
                ds,
                "publisher_name",
                "Service Contact Name",
                BaseCheck.MEDIUM,
            ),
            self._has_attr(
                ds,
                "Conventions",
                "Data Format Template Version",
                BaseCheck.MEDIUM,
            ),
            self._has_attr(
                ds,
                "publisher_name",
                "Station Publisher Name",
                BaseCheck.HIGH,
            ),
        ]

    def check_variable_attributes(self, ds):
        """
        Check IOOS concepts that come from NC variable attributes.

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        return [
            self._has_var_attr(ds, "platform", "long_name", "Station Long Name"),
            self._has_var_attr(ds, "platform", "short_name", "Station Short Name"),
            self._has_var_attr(ds, "platform", "source", "Platform Type"),
            self._has_var_attr(ds, "platform", "ioos_name", "Station ID"),
            self._has_var_attr(ds, "platform", "wmo_id", "Station WMO ID"),
            self._has_var_attr(ds, "platform", "comment", "Station Description"),
        ]

    def check_variable_names(self, ds):
        """
        Ensures all variables have a standard_name set.
        """
        msgs = []
        count = 0

        for k, v in ds.variables.items():
            if "standard_name" in v.ncattrs():
                count += 1
            else:
                msgs.append(f"Variable '{k}' missing standard_name attr")

        return Result(
            BaseCheck.MEDIUM,
            (count, len(ds.variables)),
            "Variable Names",
            msgs,
        )

    def check_altitude_units(self, ds):
        """
        If there's a variable named z, it must have units.

        @TODO: this is duplicated with check_variable_units
        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        if "z" in ds.variables:
            msgs = []
            val = "units" in ds.variables["z"].ncattrs()
            if not val:
                msgs.append("Variable 'z' has no units attr")
            return Result(BaseCheck.LOW, val, "Altitude Units", msgs)

        return Result(
            BaseCheck.LOW,
            (0, 0),
            "Altitude Units",
            ["Dataset has no 'z' variable"],
        )

    def check_variable_units(self, ds):
        """
        Ensures all variables have units.
        """
        msgs = []
        count = 0

        for k, v in ds.variables.items():
            if "units" in v.ncattrs():
                count += 1
            else:
                msgs.append(f"Variable '{k}' missing units attr")

        return Result(
            BaseCheck.MEDIUM,
            (count, len(ds.variables)),
            "Variable Units",
            msgs,
        )


class IOOS1_1Check(IOOSNCCheck):
    """
    Compliance checker implementation of IOOS Metadata Profile, Version 1.1

    Related links:
    https://ioos.github.io/ioos-metadata/ioos-metadata-profile-v1-1.html#ioos-netcdf-metadata-profile-attributes
    https://github.com/ioos/compliance-checker/issues/69
    https://github.com/ioos/compliance-checker/issues/358
    """

    _cc_spec_version = "1.1"
    _cc_description = "IOOS Metadata Profile, Version 1.1"
    _cc_url = "https://ioos.github.io/ioos-metadata/ioos-metadata-profile-v1-1.html#ioos-netcdf-metadata-profile-attributes"
    register_checker = True

    def __init__(self):
        # Define the global attributes
        self.required_atts = [
            "contributor_name",
            "contributor_role",
            "creator_country",
            "creator_email",
            "creator_sector",
            "featureType",
            "id",
            "institution",
            "naming_authority",
            "platform",
            "platform_vocabulary",
            "publisher_country",
            "publisher_email",
            "publisher_name",
            "standard_name_vocabulary",
            "title",
        ]

        self.rec_atts = [
            "creator_address",
            "creator_city",
            "creator_name",
            "creator_phone",
            "creator_state",
            "creator_url",
            "creator_zipcode",
            "keywords",
            "license",
            "publisher_address",
            "publisher_city",
            "publisher_phone",
            "publisher_state",
            "publisher_url",
            "publisher_zipcode",
            "summary",
        ]

    @check_has(BaseCheck.HIGH)
    def check_high(self, ds):
        """
        Performs a check on each highly recommended attributes' existence in the dataset

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        return self.required_atts

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        """
        Performs a check on each recommended attributes' existence in the dataset

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        return self.rec_atts

    def check_platform_variables(self, ds):
        """
        The value of platform attribute should be set to another variable which
        contains the details of the platform. There can be multiple platforms
        involved depending on if all the instances of the featureType in the
        collection share the same platform or not. If multiple platforms are
        involved, a variable should be defined for each platform and referenced
        from the geophysical variable in a space separated string.

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        platform_names = getattr(ds, "platform", "").split(" ")
        val = all(platform_name in ds.variables for platform_name in platform_names)
        msgs = []
        if not val:
            msgs = [
                (
                    'The value of "platform" global attribute should be set to another variable '
                    "which contains the details of the platform. If multiple platforms are "
                    "involved, a variable should be defined for each platform and referenced "
                    "from the geophysical variable in a space separated string."
                ),
            ]
        return [Result(BaseCheck.HIGH, val, "platform variables", msgs)]

    def check_platform_variable_attributes(self, ds):
        """
        Platform variables must contain the following attributes:
            ioos_code
            long_name
            short_name
            type

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        results = []
        platform_name = getattr(ds, "platform", "")
        # There can be multiple platforms defined here (space separated)
        for platform in platform_name.split(" "):
            if platform in ds.variables:
                results += [
                    self._has_var_attr(ds, platform, "long_name", "Platform Long Name"),
                    self._has_var_attr(
                        ds,
                        platform,
                        "short_name",
                        "Platform Short Name",
                    ),
                    self._has_var_attr(ds, platform, "ioos_code", "Platform IOOS Code"),
                    self._has_var_attr(ds, platform, "type", "Platform Type"),
                ]
        return results

    def check_geophysical_vars_fill_value(self, ds):
        """
        Check that geophysical variables contain fill values.

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        results = []
        for geo_var in cfutil.get_geophysical_variables(ds):
            results.append(
                self._has_var_attr(
                    ds,
                    geo_var,
                    "_FillValue",
                    "_FillValue",
                    BaseCheck.MEDIUM,
                ),
            )
        return results

    def check_geophysical_vars_standard_name(self, ds):
        """
        Check that geophysical variables contain standard names.

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        results = []
        for geo_var in cfutil.get_geophysical_variables(ds):
            results.append(
                self._has_var_attr(
                    ds,
                    geo_var,
                    "standard_name",
                    "geophysical variables standard_name",
                ),
            )
        return results

    def check_units(self, ds):
        """
        Required for most all variables that represent dimensional quantities.
        The value should come from udunits authoritative vocabulary, which is
        documented in the CF standard name table with it's corresponding
        standard name.

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        cf16 = CF1_6Check()
        return cf16.check_units(ds)


class IOOS1_2_ConventionsValidator(base.RegexValidator):
    validator_regex = r"\bIOOS-1.2\b"
    validator_fail_msg = '{} must contain the string "IOOS 1.2"'


class IOOS1_2_PlatformIDValidator(base.RegexValidator):
    validator_regex = r"^[a-zA-Z0-9]+$"
    validator_fail_msg = "{} must be alphanumeric"


class NamingAuthorityValidator(base.UrlValidator):
    """
    Class to check for URL or reversed DNS strings contained within
    naming_authority
    """

    validator_fail_msg = (
        '{} should either be a URL or a reversed DNS name (e.g "edu.ucar.unidata")'
    )

    def validator_func(self, input_value):
        return (
            # also check for reverse DNS strings
            super().validator_func(input_value)
            or validators.domain(".".join(input_value.split(".")[::-1]))
        )


class IOOS1_2Check(IOOSNCCheck):
    """
    Class to implement the IOOS Metadata 1.2 Specification
    """

    _cc_spec_version = "1.2"
    _cc_description = "IOOS Metadata Profile, Version 1.2"
    _cc_url = "https://ioos.github.io/ioos-metadata/ioos-metadata-profile-v1-2.html"
    register_checker = True

    def __init__(self):
        # instantiate objects used for delegation
        self.acdd1_6 = ACDD1_3Check()
        self.cf1_7 = CF1_7Check()

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
            "syntax_test_quality_flag",
        ]
        self.cf1_7._std_names._names.extend(self._qartod_std_names)

        self._default_check_var_attrs = {
            ("_FillValue", BaseCheck.MEDIUM),
            ("missing_value", BaseCheck.MEDIUM),
            # ( "standard_name", BaseCheck.HIGH # already checked in CF1_7Check.check_standard_name()
            # ( "units", BaseCheck.HIGH # already checked in CF1_7Check.check_units()
        }

        # geophysical variables must have the following attrs:
        self.geophys_check_var_attrs = self._default_check_var_attrs.union(
            {
                ("standard_name_url", BaseCheck.MEDIUM),
                # ( "platform", BaseCheck.HIGH) # checked under check_single_platform()
                # ( "wmo_platform_code", BaseCheck.HIGH) # only "if applicable", see check_wmo_platform_code()
                # ( "ancillary_variables", BaseCheck.HIGH) # only "if applicable", see _check_var_gts_ingest()
                # ("accuracy", BaseCheck.MEDIUM), see check_accuracy
                ("precision", BaseCheck.MEDIUM),
                ("resolution", BaseCheck.MEDIUM),
            },
        )

        # valid contributor_role values
        self.valid_contributor_roles = {  # NERC and NOAA
            "author",
            "coAuthor",
            "collaborator",
            "contributor",
            "custodian",
            "distributor",
            "editor",
            "funder",
            "mediator",
            "originator",
            "owner",
            "pointOfContact",
            "principalInvestigator",
            "processor",
            "publisher",
            "resourceProvider",
            "rightsHolder",
            "sponsor",
            "stakeholder",
            "user",
        }

        self.valid_contributor_role_vocabs = {
            "http://vocab.nerc.ac.uk/collection/G04/current/",
            "https://vocab.nerc.ac.uk/collection/G04/current/",
            "http://www.ngdc.noaa.gov/wiki/index.php?title=ISO_19115_and_19115-2_CodeList_Dictionaries#CI_RoleCode",
            "https://www.ngdc.noaa.gov/wiki/index.php?title=ISO_19115_and_19115-2_CodeList_Dictionaries#CI_RoleCode",
        }

        self.required_atts = [
            ("Conventions", IOOS1_2_ConventionsValidator()),
            "creator_country",
            ("creator_email", base.EmailValidator()),
            "creator_institution",
            (
                "creator_sector",
                {
                    "gov_state",
                    "nonprofit",
                    "tribal",
                    "other",
                    "unknown",
                    "gov_municipal",
                    "industry",
                    "gov_federal",
                    "academic",
                },
            ),
            ("creator_url", base.UrlValidator()),
            "featureType",
            "id",
            ("infoUrl", base.UrlValidator()),
            "license",
            ("naming_authority", NamingAuthorityValidator()),
            #'platform', # checked in check_platform_global # noqa
            "platform_name",
            "publisher_country",
            ("publisher_email", base.EmailValidator()),
            "publisher_institution",
            ("publisher_url", base.UrlValidator()),
            # TODO: handle standard name table exclusion for v38?
            (
                "standard_name_vocabulary",
                re.compile(r"^CF Standard Name Table v[1-9]\d*$"),
            ),
            "summary",
            "title",
        ]

        self.rec_atts = [
            ("contributor_email", base.EmailValidator(base.csv_splitter)),
            "contributor_name",
            ("contributor_url", base.UrlValidator(base.csv_splitter)),
            "creator_address",
            "creator_city",
            "creator_name",
            "creator_phone",
            "creator_postalcode",
            "creator_state",
            # checked in check_creator_and_publisher_type
            #'creator_type', # noqa
            "institution",
            "instrument",
            # checked in check_ioos_ingest
            #'ioos_ingest', # noqa
            "keywords",
            ("platform_id", IOOS1_2_PlatformIDValidator()),  # alphanumeric only
            "publisher_address",
            "publisher_city",
            "publisher_name",
            "publisher_phone",
            "publisher_postalcode",
            "publisher_state",
            # checked in check_creator_and_publisher_type
            #'publisher_type', # noqa
            "references",
            "instrument_vocabulary",
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
        plat_vars = ds.get_variables_by_attributes(
            platform=lambda p: isinstance(p, str),
        )
        return {
            ds.variables[var.platform]
            for var in plat_vars
            if var.platform in ds.variables
        }

    @check_has(BaseCheck.HIGH)
    def check_high(self, ds):
        """
        Performs a check on each highly recommended attributes' existence in the dataset

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        return self.required_atts

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        """
        Performs a check on each recommended attributes' existence in the dataset

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        return self.rec_atts

    def check_standard_name(self, ds):
        """
        Wrapper for checking standard names using the CF module.
        Extends the StandardNameTable to include QARTOD variable
        standard names.
        """

        return self.cf1_7.check_standard_name(ds)

    def check_feature_type(self, ds):
        """
        Wrapper for checking featureType global attribute using the CF module.
        """
        return self.cf1_7.check_feature_type(ds)

    def check_units(self, ds):
        """
        Wrapper to check units with the CF module.
        """

        return self.cf1_7.check_units(ds)

    def check_ioos_ingest(self, ds):
        """
        If a dataset contains the global attribute ioos_ingest,
        its value must be "false". All datasets are assumed to be
        ingested except those with this flag. If the dataset should
        be ingested, no flag (or "true") should be present.

        Parameters
        ----------
        ds: netCDF4.Dataset (open)

        Returns
        -------
        Result
        """

        r = True
        m = (
            "To disallow harvest of this dataset to IOOS national products, "
            'global attribute "ioos_ingest" must be a string with value "false"'
        )
        igst = getattr(ds, "ioos_ingest", None)
        if (isinstance(igst, str) and igst.lower() not in ("true", "false")) or (
            not isinstance(igst, str) and igst is not None
        ):
            r = False

        return Result(BaseCheck.MEDIUM, r, "ioos_ingest", None if r else [m])

    def check_contributor_role_and_vocabulary(self, ds):
        """
        Check the dataset has global attributes contributor_role and
        contributor_role_vocabulary. It is recommended to come from
        one of NERC or NOAA-NCEI.

        Parameters
        ----------
        ds: netCDF4.Dataset (open)

        Returns
        -------
        list of Result objects
        """

        role = getattr(ds, "contributor_role", None)
        vocb = getattr(ds, "contributor_role_vocabulary", None)

        role_val = False
        vocb_val = False

        role_msg = "contributor_role '{}' should be from NERC or NOAA-NCEI"
        vocb_msg = "contributor_role_vocabulary '{}' should be one of NERC or NOAA-NCEI"

        role_results = []
        if role:
            # in case it's a CSV, split it and iterate through all
            try:
                _roles = base.csv_splitter(role)
                for _role in _roles:
                    role_val = _role in self.valid_contributor_roles
                    role_results.append(
                        Result(
                            BaseCheck.MEDIUM,
                            role_val,
                            "contributor_role",
                            None if role_val else [role_msg.format(_role)],
                        ),
                    )
            except TypeError:
                role_results.append(
                    Result(
                        BaseCheck.MEDIUM,
                        False,
                        "contributor_role",
                        [f"contributor_role '{role}' must be of type 'string'"],
                    ),
                )
        else:
            role_results.append(
                Result(
                    BaseCheck.MEDIUM,
                    False,
                    "contributor_role",
                    ["contributor_role should be present"],
                ),
            )

        vocb_results = []
        if vocb:
            try:
                _vocbs = base.csv_splitter(vocb)
                for _vocb in _vocbs:
                    vocb_val = _vocb in self.valid_contributor_role_vocabs
                    vocb_results.append(
                        Result(
                            BaseCheck.MEDIUM,
                            vocb_val,
                            "contributor_role_vocabulary",
                            None if vocb_val else [vocb_msg.format(_vocb)],
                        ),
                    )
            except TypeError:
                vocb_results.append(
                    Result(
                        BaseCheck.MEDIUM,
                        False,
                        "contributor_role_vocabulary",
                        [
                            f"contributor_role_vocabulary '{vocb}' must be of type 'string'",
                        ],
                    ),
                )
        else:
            vocb_results.append(
                Result(
                    BaseCheck.MEDIUM,
                    False,
                    "contributor_role_vocabulary",
                    ["contributor_role_vocabulary should be present"],
                ),
            )

        return role_results + vocb_results

    def check_geophysical_vars_have_attrs(self, ds):
        """
        All geophysical variables must have certain attributes.

        Parameters
        ----------
        ds: netCDF4.Dataset

        Returns
        -------
        list: list of Result objects
        """

        # get geophysical variables
        geophys_vars = cfutil.get_geophysical_variables(ds)  # list of str
        results = self._check_vars_have_attrs(  # list
            ds,
            geophys_vars,
            self.geophys_check_var_attrs,
        )

        return results

    def check_accuracy(self, ds):
        """
        Special check for accuracy when in the salinity context.
        https://github.com/ioos/compliance-checker/issues/839

        Parameters
        ----------
        ds: netCDF4.Dataset

        Returns
        -------
        list of Results objects
        """

        results = []
        msg = (
            "Variable '{v}' attribute 'accuracy' should have the " "same units as '{v}'"
        )
        for v in cfutil.get_geophysical_variables(ds):
            _v = ds.variables[v]
            std_name = getattr(_v, "standard_name", None)
            gts_ingest = getattr(_v, "gts_ingest", None)
            if (std_name == "sea_water_practical_salinity") and (gts_ingest == "true"):
                msg = (
                    "Variable '{v}' should have an 'accuracy' attribute "
                    "that is numeric and of the same units as '{v}'"
                )
                r = isinstance(getattr(_v, "accuracy", None), Number)
            else:  # only test if exists
                r = getattr(_v, "accuracy", None) is not None

            results.append(
                Result(
                    BaseCheck.MEDIUM,
                    r,
                    "geophysical_variable:accuracy",
                    [msg.format(v=v)],
                ),
            )

        return results

    def _check_vars_have_attrs(self, ds, vars_to_check, atts_to_check):
        """
        Check that the variables in vars_to_check have the attributes in
        atts_to_check.

        Parameters
        ----------
        ds: netCDF4.Dataset (open)

        Returns
        -------
        list of Result objects
        """

        results = []
        for var in vars_to_check:
            for attr_tuple in atts_to_check:
                results.append(
                    self._has_var_attr(
                        ds,
                        var,
                        attr_tuple[0],  # attribute name
                        attr_tuple[0],  # attribute name used as 'concept_name'
                        attr_tuple[1],  # priority level
                    ),
                )
        return results

    def check_cf_role_variables(self, ds):
        """
        The IOOS-1.2 specification details the following requirements regarding
        the cf_role attribute and its relation to variable dimensionality:

          cf_role may be applied to the "Platform Variable", as indicated by
          geophysical_variable:platform, but it may also be an independent
          variable. To comply with the single platform per dataset rule of
          the IOOS Metadata Profile, the cf_role variable will typically
          have a dimension of 1, unless it is a TimeSeries dataset following
          the 'TimeSeries - multiple station' format.

        To summarize the rules checked in this method:
          - 'timeseries', cf_role var must have dim 1
          - 'timeseriesprofile' must have
            cf_role=timeseries_id variable have dim 1 and dim of cf_role=profile_id
            can be > 1
          - 'trajectory' or 'trajectoryprofile' variable with cf_role=trajectory_id
            must have dim 1, cf_role=profile_id variable can be > 1

        Relevant documentation found in the specification as well as GitHub issues:
        https://github.com/ioos/compliance-checker/issues/748#issuecomment-606659685
        https://github.com/ioos/compliance-checker/issues/828
        """

        feature_type_attr = getattr(ds, "featureType", None)
        # can't do anything, pass
        if not feature_type_attr or not isinstance(feature_type_attr, str):
            return Result(
                BaseCheck.MEDIUM,
                False,
                "CF DSG: Invalid featureType",
                [
                    (
                        f"Invalid featureType '{feature_type_attr}'; please see the "
                        "IOOS 1.2 Profile and CF-1.7 Conformance documents for valid featureType"
                    ),
                ],
            )

        feature_type = feature_type_attr.lower()

        if feature_type == "timeseries":
            return self._check_feattype_timeseries_cf_role(ds)

        elif feature_type == "timeseriesprofile":
            return self._check_feattype_timeseriesprof_cf_role(ds)

        elif feature_type == "trajectory":
            return self._check_feattype_trajectory_cf_role(ds)

        elif feature_type == "trajectoryprofile":
            return self._check_feattype_trajectoryprof_cf_role(ds)

        elif feature_type == "profile":
            return self._check_feattype_profile_cf_role(ds)

        elif feature_type == "point":
            return Result(
                BaseCheck.MEDIUM,
                True,
                "CF DSG: featureType=trajectoryProfile",
            )

        else:
            return Result(
                BaseCheck.MEDIUM,
                False,
                "CF DSG: Unknown featureType",
                [
                    (
                        f"Invalid featureType '{feature_type_attr}'; "
                        "please see the IOOS 1.2 Profile and CF-1.7 "
                        "Conformance documents for valid featureType"
                    ),
                ],
            )

    def _check_feattype_timeseries_cf_role(self, ds):
        ts_msg = (
            "Dimension length of variable with cf_role={cf_role} "
            "(the '{dim_type}' dimension) is {dim_len}. "
            "The IOOS Profile restricts timeSeries "
            "datasets with multiple features to share the same lat/lon position "
            "(i.e. to exist on the same platform). Datasets that include multiple "
            "platforms are not valid and will cause harvesting errors."
        )

        # looking for cf_role=timeseries_id
        cf_role_vars = ds.get_variables_by_attributes(cf_role="timeseries_id")
        if (not cf_role_vars) or (len(cf_role_vars) > 1):
            _val = False
            msgs = [
                (
                    "The IOOS-1.2 Profile specifies a single variable "
                    "must be present with attribute cf_role=timeseries_id"
                ),
            ]

        else:
            _v = cf_role_vars[0]
            _dims = _v.get_dims()
            if not _dims:
                _dimsize = 0
            else:
                _dimsize = _dims[0].size

            # dimension size must be == 1
            _val = _dimsize == 1
            msgs = [
                ts_msg.format(
                    cf_role="timeseries_id",
                    dim_type="station",
                    dim_len=_dimsize,
                ),
            ]

        return Result(
            BaseCheck.HIGH,
            _val,
            "CF DSG: featureType=timeseries",
            msgs,
        )

    def _check_feattype_timeseriesprof_cf_role(self, ds):
        ts_prof_msg = (
            "Dimension length of non-platform variable with cf_role={cf_role} "
            " (the '{dim_type}' dimension) is {dim_len}. "
            "The IOOS profile restricts timeSeriesProfile datasets to a "
            "single platform (ie. station) per dataset "
            "(the profile dimension is permitted to be >= 1."
        )

        # looking for cf_roles timeseries_id and profile_id
        cf_role_vars = []  # extend in specific order for easier checking
        cf_role_vars.extend(ds.get_variables_by_attributes(cf_role="timeseries_id"))
        cf_role_vars.extend(ds.get_variables_by_attributes(cf_role="profile_id"))

        if len(cf_role_vars) != 2:
            _val = False
            msgs = [
                (
                    "Datasets of featureType=timeSeriesProfile must have variables "
                    "containing cf_role=timeseries_id and cf_role=profile_id"
                ),
            ]

        else:
            _ts_id_dims = cf_role_vars[0].get_dims()  # timeseries_id dimensions
            if not _ts_id_dims:
                _ts_id_dimsize = 0
            else:
                _ts_id_dimsize = _ts_id_dims[0].size

            _pf_id_dims = cf_role_vars[1].get_dims()  # profilie_id dimensions
            if not _pf_id_dims:
                _pf_id_dimsize = 0
            else:
                _pf_id_dimsize = _pf_id_dims[0].size

            # timeseries_id must be == 1, profile >= 1
            _val = _ts_id_dimsize == 1 and _pf_id_dimsize >= 1
            msgs = [
                ts_prof_msg.format(
                    cf_role="timeseries_id",
                    dim_type="station",
                    dim_len=_ts_id_dimsize,
                ),
            ]

        return Result(
            BaseCheck.HIGH,
            _val,
            "CF DSG: featureType=timeSeriesProfile",
            msgs,
        )

    def _check_feattype_trajectory_cf_role(self, ds):
        trj_msg = (
            "Dimension length of non-platform variable with cf_role={cf_role} "
            " (the '{dim_type}' dimension) is {dim_len}. "
            "The IOOS profile restricts trjectory "
            "datasets to a single platform (i.e. trajectory) per dataset."
        )

        cf_role_vars = ds.get_variables_by_attributes(cf_role="trajectory_id")

        if len(cf_role_vars) != 1:
            _val = False
            msgs = [
                (
                    "Datasets of featureType=trajectory must have a variable "
                    "containing cf_role=trajectory_id"
                ),
            ]

        else:
            _v = cf_role_vars[0]
            _dims = _v.get_dims()
            if not _dims:
                _dimsize = 0
            else:
                _dimsize = _dims[0].size

            # trajectory dimension must be 1
            _val = _dimsize == 1
            msgs = [
                trj_msg.format(
                    cf_role="trajectory_id",
                    dim_type="station",
                    dim_len=_dimsize,
                ),
            ]

        return Result(BaseCheck.HIGH, _val, "CF DSG: featureType=trajectory", msgs)

    def _check_feattype_trajectoryprof_cf_role(self, ds):
        trj_prof_msg = (
            "Dimension length of non-platform variable with cf_role={cf_role} "
            "(the '{dim_type}' dimension) is {dim_len}. "
            "The IOOS profile restricts trajectory and trajectoryProfile "
            "datasets to a single platform (ie. trajectory) per dataset "
            "(the profile dimension is permitted to be >= 1)."
        )

        # looking for cf_roles trajectory_id and profile_id
        cf_role_vars = []  # extend in specific order for easier checking
        cf_role_vars.extend(ds.get_variables_by_attributes(cf_role="trajectory_id"))
        cf_role_vars.extend(ds.get_variables_by_attributes(cf_role="profile_id"))

        if len(cf_role_vars) != 2:
            _val = False
            msgs = [
                (
                    "Datasets of featureType=trajectoryProfile must have variables "
                    "containing cf_role=trajectory_id and cf_role=profile_id"
                ),
            ]

        else:
            _trj_id_dims = cf_role_vars[0].get_dims()
            if not _trj_id_dims:
                _trj_id_dimsize = 0
            else:
                _trj_id_dimsize = _trj_id_dims[0].size

            _prf_id_dims = cf_role_vars[1].get_dims()

            if not _prf_id_dims:
                _prf_id_dimsize = 0
            else:
                _prf_id_dimsize = _prf_id_dims[0].size

            # trajectory dim must be == 1, profile must be >= 1
            _val = _trj_id_dimsize == 1 and _prf_id_dimsize >= 1
            msgs = [
                trj_prof_msg.format(
                    cf_role="trajectory_id",
                    dim_type="station",
                    dim_len=_trj_id_dimsize,
                ),
            ]

        return Result(
            BaseCheck.HIGH,
            _val,
            "CF DSG: featureType=trajectoryProfile",
            msgs,
        )

    def _check_feattype_profile_cf_role(self, ds):
        prof_msg = (
            "Dimension length of non-platform variable with cf_role={cf_role} "
            " (the '{dim_type}' dimension) is {dim_len}. "
            "The IOOS profile restricts profile datasets to a single "
            "platform (ie. profile) per dataset."
        )

        # looking for cf_role=profile_id
        cf_role_vars = ds.get_variables_by_attributes(cf_role="profile_id")
        if (not cf_role_vars) or (len(cf_role_vars) > 1):
            _val = False
            msgs = [
                "None or multiple variables found with cf_role=profile_id; only one is allowed",
            ]

        else:
            _v = cf_role_vars[0]
            _dims = _v.get_dims()
            if not _dims:
                _dimsize = 0
            else:
                _dimsize = _dims[0].size

            # only one profile valid
            _val = _dimsize == 1
            msgs = [
                prof_msg.format(
                    cf_role="profile_id",
                    dim_type="profile",
                    dim_len=_dimsize,
                ),
            ]

        return Result(BaseCheck.HIGH, _val, "CF DSG: featureType=profile", msgs)

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
                    messages.append(
                        f"If specified, {global_att_name} must be in value list "
                        f"({sorted(expected_types)})",
                    )

            result_list.append(
                Result(BaseCheck.MEDIUM, pass_stat, global_att_name, messages),
            )

        return result_list

    def check_platform_global(self, ds):
        """
        The "platform" attribute must be a single string containing
        no blank characters.

        Parameters
        ----------
        ds: netCDF4.Dataset (open)

        Returns
        -------
        Result
        """

        r = False
        m = (
            'The global attribute "platform" must be a single string '
            + "containing no blank characters; it is {}"
        )
        p = getattr(ds, "platform", None)
        if p:
            if re.match(r"^\S+$", p):
                r = True

        return Result(BaseCheck.HIGH, r, "platform", None if r else [m.format(p)])

    def check_single_platform(self, ds):
        """
        Verify that a dataset only has a single platform attribute, and thus
        a single platform variable. Gridded model datasets are not required
        to declare a platform or platform variables.

        Args:
            ds (netCDF-4 Dataset): open Dataset object

        Returns:
            Result
        """

        glb_platform = getattr(ds, "platform", None)

        platform_set = set()
        for v in ds.get_variables_by_attributes(platform=lambda x: x is not None):
            platform_set.add(v.getncattr("platform"))

        num_platforms = len(platform_set)
        if num_platforms > 1 and glb_platform:
            msg = f"A dataset may only have one platform; {len(platform_set)} found"
            val = False

        elif (not glb_platform) and num_platforms > 0:
            msg = 'If a platform variable exists, a global attribute "platform" must also exist'
            val = False

        elif num_platforms == 0 and glb_platform:
            msg = 'A dataset with a global "platform" attribute must have at least one platform variable'
            val = False

        elif num_platforms == 0 and (not glb_platform):
            msg = "Gridded model datasets are not required to declare a platform"
            val = True

        else:
            val = True

        return Result(BaseCheck.HIGH, val, "platform", None if val else [msg])

    def check_platform_vocabulary(self, ds):
        """
        The platform_vocabulary attribute is recommended to be a URL to
        https://mmisw.org/ont/ioos/platform or
        https://vocab.nerc.ac.uk/collection/L06/current/. However,
        it is required to at least be a URL.

        Args:
            ds (netCDF4.Dataset): open Dataset

        Returns:
            Result
        """

        m = "platform_vocabulary must be a valid URL"
        pvocab = getattr(ds, "platform_vocabulary", "")
        val = bool(validators.url(pvocab))
        return Result(
            BaseCheck.MEDIUM,
            val,
            "platform_vocabulary",
            None if val else [m],
        )

    def _check_gts_ingest_val(self, val):
        """
        Check that `val` is a str and is equal to "true" or "false"

        Parameters
        ----------
        val (?): value to check

        Returns
        -------
        bool
        """

        return isinstance(val, str) and val.lower() in {"true", "false"}

    def check_vertical_coordinates(self, ds):
        """
        Check that vertical units (corresponding to axis "Z") are a unit
        equivalent to one of "meter", "inch", "foot", "yard", "US_survey_foot",
        "mile", or "fathom".  Check that the vertical coordinate variable
        "positive" attribute is either "up" or "down".  Note that unlike the CF
        version of this check, pressure units are not accepted and length units
        are constrained to the aforementioned set instead of accepting any valid
        UDUNITS length unit.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []
        for name in cfutil.get_z_variables(ds):
            variable = ds.variables[name]
            units_str = getattr(variable, "units", None)
            positive = getattr(variable, "positive", None)
            expected_unit_strs = (
                "meter",
                "inch",
                "foot",
                "yard",
                "US_survey_foot",
                "mile",
                "fathom",
            )

            unit_def_set = {
                Unit(unit_str).definition for unit_str in expected_unit_strs
            }

            try:
                units = Unit(units_str)
                pass_stat = units.definition in unit_def_set
            # unknown unit not convertible to UDUNITS
            except ValueError:
                pass_stat = False

            valid_vertical_coord = TestCtx(BaseCheck.HIGH, "Vertical coordinates")
            units_set_msg = (
                f"{name}'s units attribute {units_str} is not equivalent to one "
                f"of {expected_unit_strs}"
            )
            valid_vertical_coord.assert_true(pass_stat, units_set_msg)

            pos_msg = (
                f"{name}: vertical coordinates must include a positive "
                "attribute that is either 'up' or 'down'"
            )
            valid_vertical_coord.assert_true(positive in ("up", "down"), pos_msg)

            ret_val.append(valid_vertical_coord.to_result())

        return ret_val

    def check_gts_ingest_global(self, ds):
        """
        Check if a dataset has the global attribute "gts_ingest" and that
        it matches "true" or "false". This attribute is "required, if applicable".

        Parameters
        ----------
        ds (netCDF4.Dataset): open dataset

        Returns
        -------
        Result
        """

        gts_ingest_value = getattr(ds, "gts_ingest", None)

        is_valid_string = True
        if isinstance(gts_ingest_value, str):
            is_valid_string = self._check_gts_ingest_val(gts_ingest_value)

        fail_message = [
            'Global attribute "gts_ingest" must be a string "true" or "false"',
        ]
        return Result(
            BaseCheck.HIGH,
            is_valid_string,
            "NDBC/GTS Ingest Requirements",
            None if is_valid_string else fail_message,
        )

    def _var_qualifies_for_gts_ingest(self, ds, var):
        """
        Examine a variable to see if it qualifies for GTS Ingest.
        Check that a given variable
          - has a valid CF standard name (checked with check_standard_names())
          - has a QARTOD aggregates variable
          - has valid units (checked with check_units())

        Parameters
        ----------
        ds (netCDF4.Dataset): open Dataset
        var (netCDF4.Variable): variable from dataset

        Returns
        -------
        bool
        """

        # should have an ancillary variable with standard_name aggregate_quality_flag
        avar_val = False
        anc_vars = str(getattr(var, "ancillary_variables", "")).split(" ")
        for av in anc_vars:
            if av in ds.variables:
                if (
                    getattr(ds.variables[av], "standard_name", "")
                    == "aggregate_quality_flag"
                ):
                    avar_val = True
                    break

        # should have compliant standard_name
        # NOTE: standard names are checked extensively in self.check_standard_names()
        # but that method delegates to CF1_7Check.check_standard_name(), which loops through
        # ALL the variables; this takes the absolute core of that check and ASSUMES that the
        # current variable being checked is a coordinate variable, auxiliary coordinate
        # variable, axis variable, flag variable, or geophysical variable
        std_name = getattr(var, "standard_name", False)
        valid_std_name = std_name in self.cf1_7._std_names

        # should have compliant units
        # NOTE: units are checked extensively in self.check_units(), which delegates
        # to CF1_7Check.check_units() --> CF1_6Check.check_units(), which loops through
        # ALL variables; this takes the absolute core and assumes that the variable does
        # not need dimensionless units nor are the units to be compared with any known
        # deprecated ones; it would be nice to reuse machinery, but the similarly convoluted
        # CF1_6Check.check_units() method is too tangled to use directly and would cause a huge
        # time increase
        units = getattr(var, "units", None)
        has_udunits = units is not None and cfutil.units_known(units)

        return avar_val and valid_std_name and has_udunits

    def check_gts_ingest_requirements(self, ds):
        """
        Check which variables qualify for ingest.

        According to https://ioos.github.io/ioos-metadata/ioos-metadata-profile-v1-2.html#requirements-for-ioos-dataset-ndbcgts-ingest,
        the gts_ingest is "required, if applicable". Any variables which a user
        would like ingested must also contain the gts_ingest attribute with a
        value of true. The variable must:
          - have a valid CF standard_name attribute (already checked)
          - have an ancillary variable reqpresenting QARTOD aggregate flags
          - have a valid udunits units attribute (already checked)

        This check will always fail so as to notify the user which variables
        qualified/did not qualify for ingest.
        https://github.com/ioos/compliance-checker/issues/759#issuecomment-629454412

        Parameters
        ----------
        ds (netCDF4.Dataset): open Dataset

        Returns
        -------
        Result
        """

        # is dataset properly flagged for ingest?
        getattr(ds, "gts_ingest", None)

        # check variables
        all_passed_ingest_reqs = True  # default
        var_failed_ingest_msg = (
            "The following variables did not qualify for NDBC/GTS Ingest: {}\n"
        )
        var_passed_ingest_msg = (
            "The following variables qualified for NDBC/GTS Ingest: {}\n"
        )

        var_passed_ingest_reqs = set()
        for v in ds.get_variables_by_attributes(gts_ingest=lambda x: x == "true"):
            var_passed_ingest_reqs.add(
                (v.name, self._var_qualifies_for_gts_ingest(ds, v)),
            )

        # always show which variables have passed
        _var_passed = (y[0] for y in filter(lambda x: x[1], var_passed_ingest_reqs))

        all_passed_ingest_reqs = all(x[1] for x in var_passed_ingest_reqs)
        if not all_passed_ingest_reqs:
            _var_failed = (
                y[0] for y in filter(lambda x: not x[1], var_passed_ingest_reqs)
            )

        return Result(
            BaseCheck.HIGH,
            False,  # always fail
            "NDBC/GTS Ingest Requirements",
            (
                [var_passed_ingest_msg.format(", ".join(_var_passed))]
                if all_passed_ingest_reqs
                else [
                    var_passed_ingest_msg.format(", ".join(_var_passed)),
                    var_failed_ingest_msg.format(", ".join(_var_failed)),
                ]
            ),
        )

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
        instr_vars = cfutil.get_instrument_variables(ds)

        # check for component, disciminant
        for instr in instr_vars:
            if instr in ds.variables:
                compnt = getattr(ds.variables[instr], "component", None)
                m = [
                    f"component attribute of {instr} ({compnt}) must be a string",
                ]
                if compnt:
                    results.append(
                        Result(
                            BaseCheck.MEDIUM,
                            isinstance(compnt, str),
                            "instrument_variable",
                            m,
                        ),
                    )
                else:
                    results.append(
                        Result(BaseCheck.MEDIUM, True, "instrument_variable", m),
                    )

                disct = getattr(ds.variables[instr], "discriminant", None)
                m = [
                    f"discriminant attribute of {instr} ({disct}) must be a string",
                ]
                if disct:
                    results.append(
                        Result(
                            BaseCheck.MEDIUM,
                            isinstance(disct, str),
                            "instrument_variable",
                            m,
                        ),
                    )
                else:
                    results.append(
                        Result(BaseCheck.MEDIUM, True, "instrument_variable", m),
                    )

        return results

    def check_qartod_variables_flags(self, ds):
        """
        https://ioos.github.io/ioos-metadata/ioos-metadata-profile-v1-2.html#quality-controlqartod

        Check that all QARTOD variables have flag_meanings and flag_values attributes.
        Use delegation to methods in the CF module.

        Parameters
        ----------
        ds (netCDF4.Dataset): open dataset

        Returns
        -------
        list of Result objects
        """

        results = []
        # get qartod variables
        for v in ds.get_variables_by_attributes(
            standard_name=lambda x: x in self._qartod_std_names,
        ):
            missing_msg = "flag_{} not present on {}"

            # check if each has flag_values, flag_meanings
            # need isinstance() as can't compare truth value of array
            if getattr(v, "flag_values", None) is None:
                results.append(
                    Result(
                        BaseCheck.MEDIUM,
                        False,
                        "qartod_variables flags",
                        missing_msg.format("values", v.name),
                    ),
                )

            else:  # if exist, test
                results.append(self.cf1_7._check_flag_values(ds, v.name))

            if getattr(v, "flag_meanings", None) is None:
                results.append(
                    Result(
                        BaseCheck.MEDIUM,
                        False,
                        "qartod_variables flags",
                        missing_msg.format("meanings", v.name),
                    ),
                )

            else:  # if exist, test
                results.append(self.cf1_7._check_flag_meanings(ds, v.name))

        # Ensure message name is "qartod_variables flags"
        # NOTE this is a bit of a hack to shove into CF results
        for r in results:
            r.name = "qartod_variables flags"

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
        for v in ds.get_variables_by_attributes(
            standard_name=lambda x: x in self._qartod_std_names,
        ):
            attval = getattr(v, "references", None)
            if attval is None:
                msg = (
                    f'"references" attribute not present for variable {v.name}.'
                    "If present, it should be a valid URL."
                )
                val = False
            else:
                msg = f'"references" attribute for variable "{v.name}" must be a valid URL'
                val = bool(validators.url(attval))

            results.append(
                Result(
                    BaseCheck.MEDIUM,
                    val,
                    "qartod_variable:references",
                    None if val else [msg],
                ),
            )

        return results

    def check_wmo_platform_code(self, ds):
        """
        Per the spec:

        "The WMO identifier for the platform used to measure the data. This
        identifier can be any of the following types:
          - WMO ID for buoys (numeric, 5 digits)
          - WMO ID for gliders (numeric, 7 digits)
          - NWS ID (alphanumeric, 5 digits)"

        This attribute is "required, if applicable" -- a warning message will
        only show up if the attribute is present and does not conform.

        Args:
            ds (netCDF4.Dataset): open Dataset

        Returns:
            Result
        """

        valid = True
        ctxt = "wmo_platform_code"
        msg = (
            "The wmo_platform_code must be an alphanumeric string of 5 "
            "characters or a numeric string of 7 characters"
        )

        code = getattr(ds, "wmo_platform_code", None)
        if code:
            if not (
                isinstance(code, str)
                and (re.search(r"^(?:[a-zA-Z0-9]{5}|[0-9]{7})$", code))
            ):
                valid = False

        return Result(BaseCheck.HIGH, valid, ctxt, None if valid else [msg])

    def check_instrument_make_model_calib_date(self, ds):
        """
        Instrument variables should have attributes make_model and
        calibration_date. Both should be strings, with calibration_date
        following ISO-8601 date format.

        https://github.com/ioos/compliance-checker/issues/839
        """

        results = []
        ivars = cfutil.get_instrument_variables(ds)
        for v in ivars:
            _v = ds.variables[v]

            # make_model
            mm = getattr(_v, "make_model", None)

            valid = isinstance(mm, str)
            results.append(
                Result(
                    BaseCheck.MEDIUM,
                    valid,
                    "instrument_variable:make_model",
                    (
                        None
                        if valid
                        else [f"Attribute {v}:make_model ({mm}) should be a string"]
                    ),
                ),
            )

            # calibration_date
            cd = getattr(_v, "calibration_date", "")
            # thanks folks https://stackoverflow.com/questions/41129921/validate-an-iso-8601-datetime-string-in-python
            valid = bool(
                re.match(
                    r"^(-?(?:[1-9][0-9]*)?[0-9]{4})-(1[0-2]|0[1-9])-(3[01]|0[1-9]|[12][0-9])T(2[0-3]|[01][0-9]):([0-5][0-9]):([0-5][0-9])(\.[0-9]+)?(Z|[+-](?:2[0-3]|[01][0-9]):[0-5][0-9])?$",
                    cd,
                ),
            )
            results.append(
                Result(
                    BaseCheck.MEDIUM,
                    valid,
                    "instrument_variable:calibration_date",
                    (
                        None
                        if valid
                        else [
                            f"Attribute {v}:calibration_date ({cd}) should be an ISO-8601 string",
                        ]
                    ),
                ),
            )

        return results


class IOOSBaseSOSCheck(BaseCheck):
    _cc_spec = "ioos_sos"
    _cc_spec_version = "0.1"
    _cc_description = (
        "IOOS Inventory Metadata checks for the Sensor Observation System (SOS). "
        "Checks SOS functions GetCapabilities and DescribeSensor."
    )
    register_checker = True
    # requires login
    _cc_url = "http://sdf.ndbc.noaa.gov/sos/"


class IOOSSOSGCCheck(BaseSOSGCCheck, IOOSBaseSOSCheck):
    # set up namespaces for XPath
    ns = Namespaces().get_namespaces(["sos", "gml", "xlink"])
    ns["ows"] = Namespaces().get_namespace("ows110")

    @check_has(BaseCheck.HIGH)
    def check_high(self, ds):
        return []

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        return [
            (
                "service_contact_email",
                XPath(
                    "/sos:Capabilities/ows:ServiceProvider/ows:ServiceContact/ows:ContactInfo/ows:Address/ows:ElectronicMailAddress",
                    namespaces=self.ns,
                ),
            ),
            (
                "service_contact_name",
                XPath(
                    "/sos:Capabilities/ows:ServiceProvider/ows:ServiceContact/ows:IndividualName",
                    namespaces=self.ns,
                ),
            ),
            (
                "service_provider_name",
                XPath(
                    "/sos:Capabilities/ows:ServiceProvider/ows:ProviderName",
                    namespaces=self.ns,
                ),
            ),
            (
                "service_title",
                XPath(
                    "/sos:Capabilities/ows:ServiceProvider/ows:ProviderName",
                    namespaces=self.ns,
                ),
            ),
            (
                "service_type_name",
                XPath(
                    "/sos:Capabilities/ows:ServiceIdentification/ows:ServiceType",
                    namespaces=self.ns,
                ),
            ),
            (
                "service_type_version",
                XPath(
                    "/sos:Capabilities/ows:ServiceIdentification/ows:ServiceTypeVersion",
                    namespaces=self.ns,
                ),
            ),
            # ds.identification[0].observed_properties has this as well, but
            # don't want to try to shoehorn a function here
            # ('variable_names', len(ds.identification[0].observed_properties) > 0)
            (
                "variable_names",
                XPath(
                    "/sos:Capabilities/sos:Contents/sos:ObservationOfferingList/sos:ObservationOffering/sos:observedProperty",
                    namespaces=self.ns,
                ),
            ),
            (
                "data_format_template_version",
                XPath(
                    "/sos:Capabilities/ows:OperationsMetadata/ows:ExtendedCapabilities/gml:metaDataProperty[@xlink:title='ioosTemplateVersion']/gml:version",
                    namespaces=self.ns,
                ),
            ),
        ]

    @check_has(BaseCheck.LOW)
    def check_suggested(self, ds):
        return ["altitude_units"]


class IOOSSOSDSCheck(BaseSOSDSCheck, IOOSBaseSOSCheck):
    # set up namespaces for XPath
    ns = Namespaces().get_namespaces(["sml", "swe", "gml", "xlink"])

    @check_has(BaseCheck.HIGH)
    def check_high(self, ds):
        return [
            (
                "platform_sponsor",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:classification/sml:ClassifierList/sml:classifier[@name='sponsor']/sml:Term/sml:value",
                    namespaces=self.ns,
                ),
            ),
            (
                "platform_type",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:classification/sml:ClassifierList/sml:classifier[@name='platformType']/sml:Term/sml:value",
                    namespaces=self.ns,
                ),
            ),
            (
                "station_publisher_name",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:contact/sml:ContactList/sml:member[@xlink:role='http://mmisw.org/ont/ioos/definition/publisher']/sml:ResponsibleParty/sml:organizationName",
                    namespaces=self.ns,
                ),
            ),
            (
                "station_publisher_email",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:contact/sml:ContactList/sml:member[@xlink:role='http://mmisw.org/ont/ioos/definition/publisher']/sml:ResponsibleParty/sml:contactInfo/address/sml:electronicMailAddress",
                    namespaces=self.ns,
                ),
            ),
            (
                "station_id",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:identification/sml:IdentifierList/sml:identifier[@name='stationID']/sml:Term/sml:value",
                    namespaces=self.ns,
                ),
            ),
            (
                "station_long_name",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:identification/sml:IdentifierList/sml:identifier[@name='longName']/sml:Term/sml:value",
                    namespaces=self.ns,
                ),
            ),
            (
                "station_short_name",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:identification/sml:IdentifierList/sml:identifier[@name='shortName']/sml:Term/sml:value",
                    namespaces=self.ns,
                ),
            ),
            (
                "station_wmo_id",
                XPath(
                    '/sml:SensorML/sml:member/sml:System/sml:identification/sml:IdentifierList/sml:identifier/sml:Term[@definition="http://mmisw.org/ont/ioos/definition/wmoID"]/sml:value',
                    namespaces=self.ns,
                ),
            ),
            (
                "time_period",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:capabilities[@name='observationTimeRange']/swe:DataRecord/swe:field[@name='observationTimeRange']/swe:TimeRange/swe:value",
                    namespaces=self.ns,
                ),
            ),
            (
                "operator_email",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:contact/sml:ContactList/sml:member[@xlink:role='http://mmisw.org/ont/ioos/definition/operator']/sml:ResponsibleParty/sml:contactInfo/address/sml:electronicMailAddress",
                    namespaces=self.ns,
                ),
            ),
            (
                "operator_name",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:contact/sml:ContactList/sml:member[@xlink:role='http://mmisw.org/ont/ioos/definition/operator']/sml:ResponsibleParty/sml:organizationName",
                    namespaces=self.ns,
                ),
            ),
            (
                "station_description",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/gml:description",
                    namespaces=self.ns,
                ),
            ),
            # replaced with lon/lat with point
            (
                "station_location_point",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:location/gml:Point/gml:pos",
                    namespaces=self.ns,
                ),
            ),
        ]

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        return [
            (
                "sensor_descriptions",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:components/sml:ComponentList/sml:component/sml:System/gml:description",
                    namespaces=self.ns,
                ),
            ),
            (
                "sensor_ids",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:components/sml:ComponentList/sml:component/sml:System/@gml:id",
                    namespaces=self.ns,
                ),
            ),
            (
                "sensor_names",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:components/sml:ComponentList/sml:component/@name",
                    namespaces=self.ns,
                ),
            ),
            (
                "data_format_template_version",
                XPath(
                    "/sml:SensorML/sml:capabilities/swe:SimpleDataRecord/swe:field[@name='ioosTemplateVersion']/swe:Text/swe:value",
                    namespaces=self.ns,
                ),
            ),
            (
                "variable_names",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:components/sml:ComponentList/sml:component/sml:System/sml:outputs/sml:OutputList/sml:output/swe:Quantity/@definition",
                    namespaces=self.ns,
                ),
            ),
            (
                "variable_units",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:components/sml:ComponentList/sml:component/sml:System/sml:outputs/sml:OutputList/sml:output/swe:Quantity/swe:uom/@code",
                    namespaces=self.ns,
                ),
            ),
            (
                "network_id",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:capabilities[@name='networkProcedures']/swe:SimpleDataRecord/gml:metaDataProperty/@xlink:href",
                    namespaces=self.ns,
                ),
            ),
            (
                "operator_sector",
                XPath(
                    "/sml:SensorML/sml:member/sml:System/sml:classification/sml:ClassifierList/sml:classifier[@name='operatorSector']/sml:Term/sml:value",
                    namespaces=self.ns,
                ),
            ),
        ]

    @check_has(BaseCheck.LOW)
    def check_suggested(self, ds):
        return []
