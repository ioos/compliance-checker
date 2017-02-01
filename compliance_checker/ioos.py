from __future__ import unicode_literals
from compliance_checker.base import BaseCheck, BaseNCCheck, BaseSOSGCCheck, BaseSOSDSCheck, check_has, Result
from owslib.namespaces import Namespaces
from lxml.etree import XPath


class IOOSBaseCheck(BaseCheck):
    register_checker = True
    _cc_spec = 'ioos'
    _cc_spec_version = '0.1'
    _cc_description = 'IOOS Inventory Metadata'
    # requires login
    _cc_url = 'https://docs.google.com/spreadsheets/d/1huUFauh7rPj2oKfiRhLE1ZCsnes8SmAm6fKE95dsybE/'


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


class IOOSNCCheck(BaseNCCheck, IOOSBaseCheck):

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

    def check_global_attributes(self, ds):
        """
        Check all global NC attributes for existence.
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
        """
        return [
            self._has_var_attr(ds, 'platform', 'long_name', 'Station Long Name'),
            self._has_var_attr(ds, 'platform', 'short_name', 'Station Short Name'),
            self._has_var_attr(ds, 'platform', 'source', 'Platform Type'),
            self._has_var_attr(ds, 'platform', 'ioos_name', 'Station ID'),
            self._has_var_attr(ds, 'platform', 'wmo_id', 'Station WMO ID'),
            self._has_var_attr(ds, 'platform', 'comment', 'Station Description'),
        ]

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

        return Result(BaseCheck.HIGH, (count, 2), 'Time Period', msgs)

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

        return Result(BaseCheck.HIGH, (count, 2), 'Station Location Lat', msgs)

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

        return Result(BaseCheck.HIGH, (count, 2), 'Station Location Lon', msgs)

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

    def check_altitude_units(self, ds):
        """
        If there's a variable named z, it must have units.

        @TODO: this is duplicated with check_variable_units
        """
        if 'z' in ds.variables:
            msgs = []
            val = 'units' in ds.variables['z'].ncattrs()
            if not val:
                msgs.append("Variable 'z' has no units attr")
            return Result(BaseCheck.LOW, val, 'Altitude Units', msgs)

        return Result(BaseCheck.LOW, (0, 0), 'Altitude Units', ["Dataset has no 'z' variable"])


class IOOSSOSGCCheck(BaseSOSGCCheck, IOOSBaseCheck):

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
            #('variable_names', len(ds.identification[0].observed_properties) > 0)
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


class IOOSSOSDSCheck(BaseSOSDSCheck, IOOSBaseCheck):

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
