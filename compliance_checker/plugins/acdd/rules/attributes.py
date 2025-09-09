from xrlint.node import DatasetNode
from xrlint.rule import RuleContext, RuleOp

from compliance_checker.plugins.acdd.plugin import plugin

attrs_1_0_high_rec = {
    "title": "A short phrase or sentence describing the dataset. In many discovery systems, the title will be displayed in the results list from a search, and therefore should be human readable and reasonable to display in a list of such names. This attribute is also recommended by the NetCDF Users Guide and the CF conventions.",
    "keywords": "A comma-separated list of key words and/or phrases. Keywords may be common words or phrases, terms from a controlled vocabulary (GCMD is often used), or URIs for terms from a controlled vocabulary (see also 'keywords_vocabulary' attribute).",
    "summary": "A paragraph describing the dataset, analogous to an abstract for a paper.",
}

attrs_1_0_rec = {
    "id": """An identifier for the data set, provided by and unique within its naming authority. The combination of the "naming authority" and the "id" should be globally unique, but the id can be globally unique by itself also. IDs can be URLs, URNs, DOIs, meaningful text strings, a local key, or any other unique string of characters. The id should not include white space characters.""",
    "naming_authority": "The organization that provides the initial id (see above) for the dataset. The naming authority should be uniquely specified by this attribute. We recommend using reverse-DNS naming for the naming authority; URIs are also acceptable. Example: 'edu.ucar.unidata'.",
    "history": "Provides an audit trail for modifications to the original data. This attribute is also in the NetCDF Users Guide: 'This is a character array with a line for each invocation of a program that has modified the dataset. Well-behaved generic netCDF applications should append a line containing: date, time of day, user name, program name and command arguments.' To include a more complete description you can append a reference to an ISO Lineage entity; see NOAA EDM ISO Lineage guidance.",
    "comment": "Miscellaneous information about the data, not captured elsewhere. This attribute is defined in the CF Conventions.",
    "date_created": "The date on which this version of the data was created. (Modification of values implies a new version, hence this would be assigned the date of the most recent values modification.) Metadata changes are not considered when assigning the date_created. The ISO 8601:2004 extended date format is recommended, as described in the Attribute Content Guidance section.",
    "creator_name": "The name of the person (or other creator type specified by the creator_type attribute) principally responsible for creating this data.",
    "creator_url": "The URL of the person (or other creator type specified by the creator_type attribute) principally responsible for creating this data.",
    "creator_email": "The email address of the person (or other creator type specified by the creator_type attribute) principally responsible for creating this data.",
    "institution": "The name of the institution principally responsible for originating this data. This attribute is recommended by the CF convention.",
    "project": "The name of the project(s) principally responsible for originating this data. Multiple projects can be separated by commas, as described under Attribute Content Guidelines. Examples: 'PATMOS-X', 'Extended Continental Shelf Project'.",
    "processing_level": "A textual description of the processing (or quality control) level of the data.",
    "geospatial_bounds": "Describes the data's 2D or 3D geospatial extent in OGC's Well-Known Text (WKT) Geometry format (reference the OGC Simple Feature Access (SFA) specification). The meaning and order of values for each point's coordinates depends on the coordinate reference system (CRS). The ACDD default is 2D geometry in the EPSG:4326 coordinate reference system. The default may be overridden with geospatial_bounds_crs and geospatial_bounds_vertical_crs (see those attributes). EPSG:4326 coordinate values are latitude (decimal degrees_north) and longitude (decimal degrees_east), in that order. Longitude values in the default case are limited to the [-180, 180) range. Example: 'POLYGON ((40.26 -111.29, 41.26 -111.29, 41.26 -110.29, 40.26 -110.29, 40.26 -111.29))'.",
    "geospatial_bounds_vertical_crs": """The vertical coordinate reference system (CRS) for the Z axis of the point coordinates in the geospatial_bounds attribute. This attribute cannot be used if the CRS in geospatial_bounds_crs is 3-dimensional; to use this attribute, geospatial_bounds_crs must exist and specify a 2D CRS. EPSG CRSs are strongly recommended. There is no default for this attribute when not specified. Examples: 'EPSG:5829' (instantaneous height above sea level), "EPSG:5831" (instantaneous depth below sea level), or 'EPSG:5703' (NAVD88 height).""",
    "geospatial_lat_min": "Describes a simple lower latitude limit; may be part of a 2- or 3-dimensional bounding region. Geospatial_lat_min specifies the southernmost latitude covered by the dataset.",
    "geospatial_lat_max": "Describes a simple upper latitude limit; may be part of a 2- or 3-dimensional bounding region. Geospatial_lat_max specifies the northernmost latitude covered by the dataset.",
    "geospatial_lon_min": "Describes a simple longitude limit; may be part of a 2- or 3-dimensional bounding region. geospatial_lon_min specifies the westernmost longitude covered by the dataset. See also geospatial_lon_max.",
    "geospatial_lon_max": "Describes a simple longitude limit; may be part of a 2- or 3-dimensional bounding region. geospatial_lon_max specifies the easternmost longitude covered by the dataset. Cases where geospatial_lon_min is greater than geospatial_lon_max indicate the bounding box extends from geospatial_lon_max, through the longitude range discontinuity meridian (either the antimeridian for -180:180 values, or Prime Meridian for 0:360 values), to geospatial_lon_min; for example, geospatial_lon_min=170 and geospatial_lon_max=-175 incorporates 15 degrees of longitude (ranges 170 to 180 and -180 to -175).",
    "geospatial_vertical_min": "Describes the numerically smaller vertical limit; may be part of a 2- or 3-dimensional bounding region. See geospatial_vertical_positive and geospatial_vertical_units.",
    "geospatial_vertical_max": "Describes the numerically larger vertical limit; may be part of a 2- or 3-dimensional bounding region. See geospatial_vertical_positive and geospatial_vertical_units.",
    "time_coverage_start": "Describes the time of the last data point in the data set. Use ISO 8601:2004 date format, preferably the extended format as recommended in the Attribute Content Guidance section.",
    "time_coverage_end": "Describes the time of the first data point in the data set. Use the ISO 8601:2004 date format, preferably the extended format as recommended in the Attribute Content Guidance section.",
    "time_coverage_duration": "Describes the duration of the data set. Use ISO 8601:2004 duration format, preferably the extended format as recommended in the Attribute Content Guidance section.",
    "time_coverage_resolution": "Describes the targeted time period between each value in the data set. Use ISO 8601:2004 duration format, preferably the extended format as recommended in the Attribute Content Guidance section.",
    "standard_name_vocabulary": "The name and version of the controlled vocabulary from which variable standard names are taken. (Values for any standard_name attribute must come from the CF Standard Names vocabulary for the data file or product to comply with CF.) Example: 'CF Standard Name Table v27'.",
    "license": """Provide the URL to a standard or specific license, enter "Freely Distributed" or "None", or describe any restrictions to data access and distribution in free text.""",
}

attrs_1_0_sug = {
    "contributor_name": "The name of any individuals, projects, or institutions that contributed to the creation of this data. May be presented as free text, or in a structured format compatible with conversion to ncML (e.g., insensitive to changes in whitespace, including end-of-line characters).",
    "contributor_role": "The role of any individuals, projects, or institutions that contributed to the creation of this data. May be presented as free text, or in a structured format compatible with conversion to ncML (e.g., insensitive to changes in whitespace, including end-of-line characters). Multiple roles should be presented in the same order and number as the names in contributor_names.",
    "date_modified": "The date on which the data was last modified. Note that this applies just to the data, not the metadata. The ISO 8601:2004 extended date format is recommended, as described in the Attributes Content Guidance section.",
    "date_issued": "The date on which this data (including all modifications) was formally issued (i.e., made available to a wider audience). Note that these apply just to the data, not the metadata. The ISO 8601:2004 extended date format is recommended, as described in the Attributes Content Guidance section.",
    "geospatial_lat_units": """Units for the latitude axis described in "geospatial_lat_min" and "geospatial_lat_max" attributes. These are presumed to be "degree_north"; other options from udunits may be specified instead.""",
    "geospatial_lat_resolution": "Information about the targeted spacing of points in latitude. Recommend describing resolution as a number value followed by the units. Examples: '100 meters', '0.1 degree'",
    "geospatial_lon_units": """Units for the longitude axis described in "geospatial_lon_min" and "geospatial_lon_max" attributes. These are presumed to be "degree_east"; other options from udunits may be specified instead.""",
    "geospatial_lon_resolution": "Information about the targeted spacing of points in longitude. Recommend describing resolution as a number value followed by units. Examples: '100 meters', '0.1 degree'",
    "geospatial_vertical_units": """Units for the vertical axis described in "geospatial_vertical_min" and "geospatial_vertical_max" attributes. The default is EPSG:4979 (height above the ellipsoid, in meters); other vertical coordinate reference systems may be specified. Note that the common oceanographic practice of using pressure for a vertical coordinate, while not strictly a depth, can be specified using the unit bar. Examples: 'EPSG:5829' (instantaneous height above sea level), 'EPSG:5831' (instantaneous depth below sea level).""",
    "geospatial_vertical_resolution": "Information about the targeted vertical spacing of points. Example: '25 meters'",
}

attrs_1_1_high_rec = {**attrs_1_0_high_rec}

attrs_1_1_rec = {
    **attrs_1_0_rec,
    "keywords_vocabulary": """If you are using a controlled vocabulary for the words/phrases in your "keywords" attribute, this is the unique name or identifier of the vocabulary from which keywords are taken. If more than one keyword vocabulary is used, each may be presented with a prefix and a following comma, so that keywords may optionally be prefixed with the controlled vocabulary key. Example: 'GCMD:GCMD Keywords, CF:NetCDF COARDS Climate and Forecast Standard Names'.""",
}

attrs_1_1_sug = {
    **attrs_1_0_sug,
    "publisher_name": "The name of the person (or other entity specified by the publisher_type attribute) responsible for publishing the data file or product to users, with its current metadata and format.",  # publisher,dataCenter
    "publisher_url": "The URL of the person (or other entity specified by the publisher_type attribute) responsible for publishing the data file or product to users, with its current metadata and format.",  # publisher
    "publisher_email": "The email address of the person (or other entity specified by the publisher_type attribute) responsible for publishing the data file or product to users, with its current metadata and format.",  # publisher
    "geospatial_vertical_positive": "One of 'up' or 'down'. If up, vertical values are interpreted as 'altitude', with negative values corresponding to below the reference datum (e.g., under water). If down, vertical values are interpreted as 'depth', positive values correspond to below the reference datum. Note that if geospatial_vertical_positive is down ('depth' orientation), the geospatial_vertical_min attribute specifies the data's vertical location furthest from the earth's center, and the geospatial_vertical_max attribute specifies the location closest to the earth's center.",
}

attrs_1_3_high_rec = {
    **attrs_1_0_high_rec,
    "conventions": "A comma-separated list of the conventions that are followed by the dataset. For files that follow this version of ACDD, include the string 'ACDD-1.3'. (This attribute is described in the NetCDF Users Guide.)",
}

attrs_1_3_rec = {
    **attrs_1_0_rec,
    "geospatial_vertical_positive": "One of 'up' or 'down'. If up, vertical values are interpreted as 'altitude', with negative values corresponding to below the reference datum (e.g., under water). If down, vertical values are interpreted as 'depth', positive values correspond to below the reference datum. Note that if geospatial_vertical_positive is down ('depth' orientation), the geospatial_vertical_min attribute specifies the data's vertical location furthest from the earth's center, and the geospatial_vertical_max attribute specifies the location closest to the earth's center.",
    "geospatial_bounds_crs": "The coordinate reference system (CRS) of the point coordinates in the geospatial_bounds attribute. This CRS may be 2-dimensional or 3-dimensional, but together with geospatial_bounds_vertical_crs, if that attribute is supplied, must match the dimensionality, order, and meaning of point coordinate values in the geospatial_bounds attribute. If geospatial_bounds_vertical_crs is also present then this attribute must only specify a 2D CRS. EPSG CRSs are strongly recommended. If this attribute is not specified, the CRS is assumed to be EPSG:4326. Examples: 'EPSG:4979' (the 3D WGS84 CRS), 'EPSG:4047'.",
    "geospatial_bounds_vertical_crs": """The vertical coordinate reference system (CRS) for the Z axis of the point coordinates in the geospatial_bounds attribute. This attribute cannot be used if the CRS in geospatial_bounds_crs is 3-dimensional; to use this attribute, geospatial_bounds_crs must exist and specify a 2D CRS. EPSG CRSs are strongly recommended. There is no default for this attribute when not specified. Examples: 'EPSG:5829' (instantaneous height above sea level), "EPSG:5831" (instantaneous depth below sea level), or 'EPSG:5703' (NAVD88 height).""",
    "publisher_name": "The name of the person (or other entity specified by the publisher_type attribute) responsible for publishing the data file or product to users, with its current metadata and format.",  # publisher,dataCenter
    "publisher_url": "The URL of the person (or other entity specified by the publisher_type attribute) responsible for publishing the data file or product to users, with its current metadata and format.",  # publisher
    "publisher_email": "The email address of the person (or other entity specified by the publisher_type attribute) responsible for publishing the data file or product to users, with its current metadata and format.",  # publisher
    "source": "The method of production of the original data. If it was model-generated, source should name the model and its version. If it is observational, source should characterize it. This attribute is defined in the CF Conventions. Examples: 'temperature from CTD #1234'; 'world model v.0.1'.",
}

attrs_1_3_sug = {
    **attrs_1_0_sug,
    "creator_type": "Specifies type of creator with one of the following: 'person', 'group', 'institution', or 'position'. If this attribute is not specified, the creator is assumed to be a person.",
    "creator_institution": "The institution of the creator; should uniquely identify the creator's institution. This attribute's value should be specified even if it matches the value of publisher_institution, or if creator_type is institution.",
    "platform": "Name of the platform(s) that supported the sensor data used to create this data set or product. Platforms can be of any type, including satellite, ship, station, aircraft or other. Indicate controlled vocabulary used in platform_vocabulary.",
    "platform_vocabulary": """Controlled vocabulary for the names used in the "platform" attribute.""",
    "keywords_vocabulary": """If you are using a controlled vocabulary for the words/phrases in your "keywords" attribute, this is the unique name or identifier of the vocabulary from which keywords are taken. If more than one keyword vocabulary is used, each may be presented with a prefix and a following comma, so that keywords may optionally be prefixed with the controlled vocabulary key. Example: 'GCMD:GCMD Keywords, CF:NetCDF COARDS Climate and Forecast Standard Names'.""",
    "instrument": "Name of the contributing instrument(s) or sensor(s) used to create this data set or product. Indicate controlled vocabulary used in instrument_vocabulary.",
    "metadata_link": "A URL that gives the location of more complete metadata. A persistent URL is recommended for this attribute.",
    "product_version": "Version identifier of the data file or product as assigned by the data creator. For example, a new algorithm or methodology could result in a new product_version.",
    "references": "Published or web-based references that describe the data or methods used to produce it. Recommend URIs (such as a URL or DOI) for papers or other references. This attribute is defined in the CF conventions.",
    "publisher_type": "Specifies type of publisher with one of the following: 'person', 'group', 'institution', or 'position'. If this attribute is not specified, the publisher is assumed to be a person.",
    "instrument_vocabulary": """Controlled vocabulary for the names used in the "instrument" attribute.""",
    "date_metadata_modified": "The date on which the metadata was last modified. The ISO 8601:2004 extended date format is recommended, as described in the Attributes Content Guidance section.",
    "program": "The overarching program(s) of which the dataset is a part. A program consists of a set (or portfolio) of related and possibly interdependent projects that meet an overarching objective. Examples: 'GHRSST', 'NOAA CDR', 'NASA EOS', 'JPSS', 'GOES-R'.",
    "publisher_institution": " 	The institution that presented the data file or equivalent product to users; should uniquely identify the institution. If publisher_type is institution, this should have the same value as publisher_name.",
}


class AttributesRule(RuleOp):
    attrs = {}

    def validate_dataset(self, ctx: RuleContext, node: DatasetNode):
        for attr, message in self.attrs.items():
            if attr not in node.dataset.attrs:
                ctx.report(
                    f"Missing attribute '{attr}'",
                    suggestions=[message],
                )


@plugin.define_rule(
    "1.0_attrs_highly_recommended",
    version="1.0",
    docs_url="https://wiki.esipfed.org/Category:Attribute_Conventions_Dataset_Discovery",
)
class Attributes_1_0_Highly_Reccomended(AttributesRule):
    attrs = attrs_1_0_high_rec


@plugin.define_rule(
    "1.0_attrs_recommended",
    version="1.0",
    docs_url="https://wiki.esipfed.org/Category:Attribute_Conventions_Dataset_Discovery",
)
class Attributes_1_0_Reccomended(AttributesRule):
    attrs = attrs_1_0_rec


@plugin.define_rule(
    "1.0_attrs_suggested",
    version="1.0",
    docs_url="https://wiki.esipfed.org/Category:Attribute_Conventions_Dataset_Discovery",
)
class Attributes_1_0_Suggested(AttributesRule):
    attrs = attrs_1_0_sug


@plugin.define_rule(
    "1.3_attrs_highly_recommended",
    version="1.3",
    docs_url="https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3",
)
class Attributes_1_3_Highly_Reccomended(AttributesRule):
    attrs = attrs_1_3_high_rec


@plugin.define_rule(
    "1.3_attrs_recommended",
    version="1.3",
    docs_url="https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3",
)
class Attributes_1_3_Reccomended(AttributesRule):
    attrs = attrs_1_3_rec


@plugin.define_rule(
    "1.3_attrs_suggested",
    version="1.3",
    docs_url="https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3",
)
class Attributes_1_3_Suggested(AttributesRule):
    attrs = attrs_1_3_sug
