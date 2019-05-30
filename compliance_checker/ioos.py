'''
Check for IOOS-approved attributes
'''
from __future__ import unicode_literals
from compliance_checker.base import BaseCheck, BaseNCCheck, BaseSOSGCCheck, BaseSOSDSCheck, check_has, Result
from owslib.namespaces import Namespaces
from lxml.etree import XPath
from compliance_checker.cf.cf import CF1_6Check
from compliance_checker.cfutil import get_geophysical_variables, get_instrument_variables
from compliance_checker.cf.cf import CFBaseCheck
import re

try:
    basestring
except NameError:
    basestring = str

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
    def _has_var_attr(cls, dataset, vname, attr, concept_name, priority=BaseCheck.HIGH):
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

    def check_instrument_variables(self, ds):
        '''
        Instrument variables are 'required, if applicable' by the IOOS Profile v1.1.
        If an instrument variable exists, it must have the 'discriminant' attribute.
        Since the Compliance-Checker has no way to discern if it is applicable or not,
        the check is skipped if the variable does not exist.

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        results = []
        instrument_vars = get_instrument_variables(ds)
        if not instrument_vars:
            pass

        for ivar in instrument_vars:
            results.append(
                self._has_var_attr(ds, ivar, 'discriminant', 'instrument_variable:discriminant', BaseCheck.HIGH),
            )
        return results

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


class IOOS1_2Check(IOOS1_1Check):
    """
    Compliance checker implementation of IOOS Metadata Profile, Version 1.2

    Related links:
    https://ioos.github.io/ioos-metadata/ioos-metadata-profile-v1-2.html#ioos-metadata-profile-attributes
    """
    _cc_spec_version = '1.2'
    _cc_description = 'IOOS Metadata Profile, Version 1.2'
    _cc_url = 'https://ioos.github.io/ioos-metadata/ioos-metadata-profile-v1-2.html#ioos-metadata-profile-attributes'
    register_checker = True

    def __init__(self):
        # Define the global attributes
        self.required_atts = [
            'creator_country',
            'creator_email',
            'creator_institution',
            'creator_sector',
            'creator_url',
            'featureType',
            'id',
            'info_url',
            'naming_authority',
            'platform_name',
            'platform_vocabulary',
            'publisher_country',
            'publisher_email',
            'publisher_name',
            'publisher_url',
            'standard_name_vocabulary',
            'title'
        ]

        self.rec_atts = [
            'contributor_email',
            'contributor_name',
            'contributor_role',
            'contributor_role_vocabulary',
            'contributor_url',
            'creator_address',
            'creator_city',
            'creator_name',
            'creator_phone',
            'creator_state',
            'creator_type',
            'creator_postalcode',
            'institution',
            'keywords',
            'license',
            'platform_id',
            'publisher_address',
            'publisher_city',
            'publisher_phone',
            'publisher_state',
            'publisher_type',
            'publisher_postalcode',
            'references',
            'summary'
        ]

    def check_instrument_variables(self, ds):
        '''
        Instrument variables are 'recommended, if applicable' by the IOOS Profile v1.2.
        If an instrument variable exists, it is recommended to have the following attributes:
            component
            discriminant

        Since the Compliance-Checker has no way to discern if it is applicable or not,
        the check is skipped if the variable does not exist.

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        results = []
        instrument_vars = get_instrument_variables(ds)
        if not instrument_vars:
            pass

        for ivar in instrument_vars:
            results.append(
                self._has_var_attr(ds, ivar, 'component', 'instrument_variable:component', BaseCheck.HIGH))
            results.append(
                self._has_var_attr(ds, ivar, 'discriminant', 'instrument_variable:discriminant', BaseCheck.HIGH))
        return results

    def check_platform_variables(self, ds):
        """Since only 1 platform variable is allowed per v1.2, this method
        is overloaded to skip."""
        pass

    def check_platform_variable_attributes(self, ds):
        """Many of the attributes were removed from v1.1 to v1.2, so this
        method is overloaded to skip."""
        pass

    def check_platform_variable(self, ds):
        '''
        Consolidated method to check a singular platform variable exists
        and contains the 'cf_role' attribute.
        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        platform_name = getattr(ds, 'platform', '')
        _multi_platforms = platform_name.split(" ")
        if len(_multi_platforms) > 1:
            return Result(
                BaseCheck.HIGH,
                False,
                "platform variable",
                ["Only one platform variable is allowed per Version 1.2"]
            )
        val = platform_name in ds.variables
        if not val:
            msg = 'The value of "platform" global attribute should be set to another variable '+\
                  'which contains the details of the platform.'
            return Result(BaseCheck.HIGH, val, 'platform variables', [msg])
        return self._has_var_attr(ds, platform_name, 'cf_role', 'platform variable', BaseCheck.HIGH)

    def check_geophysical_vars_standard_name(self, ds):
        """Overloaded implementation to add a check for standard_name_uri"""
        results = []
        for geo_var in get_geophysical_variables(ds):
            results.append(
                self._has_var_attr(ds, geo_var, 'standard_name', 'geophysical variables standard_name'),
            )
            results.append(
                self._has_var_attr(ds, geo_var, 'standard_name_uri', 'geophysical variables standard_name_uri'),
            )
        return results

    def check_wmo_platform_code(self, ds):
        """Check if a dataset is defined as WMO, then it must have the wmo_platform_code
        global attribute. If the attribute does not exist, the check passes without performing
        any further computation; otherwise, it checks if the attribute is a 5 or 7 character alpha-
        numeric string."""

        try:
            wmo_id = ds.getncattr("wmo_platform_code")
            m = "wmo_platform_code must be 5-7 character alphanumeric string"
            if not isinstance(wmo_id, basestring):
                score = False
            else:
                r = re.compile(r"^\w{5}(?:\w{2})?$") # alphanumeric, length 5 or 7 characters
                if r.match(wmo_id):
                    score = True
                else:
                    score = False
        except AttributeError:
            score = True
            m = "wmo_platform_code not found; passing"

        return Result(BaseCheck.HIGH, score, 'Platform', [m])


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
