# IOOS Compliance Checker

[![Build Status](https://travis-ci.org/ioos/compliance-checker.svg)](https://travis-ci.org/ioos/compliance-checker)

The IOOS Compliance Checker is a Python tool to check local/remote datasets against a variety of compliance standards.
It is primarily a command-line tool (tested on OS X/Linux) and can also be used as a library import.

It currently supports the following sources and standards:

| Standard                                                                                             | Source                                                            | .nc/OPeNDAP | SOS                             |
| ---------------------------------------------------------------------------------------------------- | -----------                                                       | ------      | ------------------------------- |
| [ACDD (1.1)](http://wiki.esipfed.org/index.php/Attribute_Convention_for_Data_Discovery_%28ACDD%29)   | Built-in                                                          | X           | -                               |
| [CF (1.6)](http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html)                       | Built-in                                                          | X           | -                               |
| IOOS Asset Concept                                                                                   | Built-in                                                          | -           | GetCapabilities, DescribeSensor |
| [Glider DAC](https://github.com/ioos/ioosngdac/wiki/NGDAC-NetCDF-File-Format-Version-2)              | [ioos/cc-plugin-glider](https://github.com/ioos/cc-plugin-glider) | X           | -                               |

## Concepts & Terminology

Each compliance standard is executed by a Check Suite,
which functions similar to a Python standard Unit Test.
A Check Suite runs one or more checks against a dataset,
returning a list of Results which are then aggregated into a summary.

Each Result has a (# passed / # total) score, a weight (HIGH/MEDIUM/LOW),
a computer-readable name, an optional list of human-readable messages,
and optionally a list of child Results.

A single score is then calculated by aggregating on the names,
then multiplying the score by the weight and summing them together.

The computer-readable name field controls how Results are aggregated together - in order to prevent the overall score for a Check Suite varying on the number of variables,
it is possible to *group* Results together via the name property.
Grouped results will only add up to a single top-level entry.

See the [Development](//github.com/ioos/compliance-checker/wiki/Development) wiki page for more details on implementation.

## Usage (command line)

The compliance-checker can work against local files (`.nc` files, `.cdl` metadata files, `.xml` files of SOS GetCapabilities/DescribeSensor requests) or against remote URLs (OPeNDAP data URLs, SOS GetCapabilities/DescribeSensor URLs).

> **WARNING** The CF/ACDD checks **will access data**, so if using a remote OPeNDAP URL, please be sure the size is reasonable!

```
usage: compliance-checker [-h] [--test TEST]
                          [--criteria [{lenient,normal,strict}]] [--verbose]
                          [--skip-checks SKIP_CHECKS] [-f {text,html,json}]
                          [-o OUTPUT] [-V] [-l] [-d DOWNLOAD_STANDARD_NAMES]
                          [dataset_location [dataset_location ...]]

positional arguments:
  dataset_location      Defines the location of the dataset to be checked.

optional arguments:
  -h, --help            show this help message and exit
  --test TEST, -t TEST, --test= TEST, -t= TEST
                        Select the Checks you want to perform. Defaults to
                        'acdd' if unspecified
  --criteria [{lenient,normal,strict}], -c [{lenient,normal,strict}]
                        Define the criteria for the checks. Either Strict,
                        Normal, or Lenient. Defaults to Normal.
  --verbose, -v         Increase output. May be specified up to three times.
  --skip-checks SKIP_CHECKS, -s SKIP_CHECKS
                        Specifies tests to skip
  -f {text,html,json}, --format {text,html,json}
                        Output format
  -o OUTPUT, --output OUTPUT
                        Output filename
  -V, --version         Display the IOOS Compliance Checker version
                        information.
  -l, --list-tests      List the available tests
  -d DOWNLOAD_STANDARD_NAMES, --download-standard-names DOWNLOAD_STANDARD_NAMES
                        Specify a version of the cf standard name table to
                        download as packaged version
```

```
$ compliance-checker --test=acdd compliance_checker/tests/data/ru07-20130824T170228_rt0.nc


--------------------------------------------------------------------------------
                     The dataset scored 69 out of 93 points                     
                             during the acdd check                              
--------------------------------------------------------------------------------
                               Scoring Breakdown:                               


                                 High Priority                                  
--------------------------------------------------------------------------------
    Name                            :Priority: Score
Conventions                             :3:     1/2
keywords                                :3:     1/1
summary                                 :3:     1/1
title                                   :3:     1/1
varattr                                 :3:    32/44


                                Medium Priority                                 
--------------------------------------------------------------------------------
    Name                            :Priority: Score
acknowledgment/acknowledgement          :2:     1/1
comment                                 :2:     1/1
creator_email                           :2:     1/1
creator_name                            :2:     1/1
creator_url                             :2:     1/1
date_created                            :2:     1/1
date_created_is_iso                     :2:     0/1
date_issued_is_iso                      :2:     0/1
date_metadata_modified_is_iso           :2:     0/0
date_modified_is_iso                    :2:     0/1
geospatial_bounds                       :2:     0/1
geospatial_bounds_crs                   :2:     0/1
geospatial_bounds_vertical_crs          :2:     0/1
geospatial_lat_extents_match            :2:     1/2
geospatial_lat_max                      :2:     1/1
geospatial_lat_min                      :2:     1/1
geospatial_lon_extents_match            :2:     1/2
geospatial_lon_max                      :2:     1/1
geospatial_lon_min                      :2:     1/1
geospatial_vertical_extents_match       :2:     0/2
geospatial_vertical_max                 :2:     1/1
geospatial_vertical_min                 :2:     1/1
geospatial_vertical_positive            :2:     1/1
history                                 :2:     1/1
id                                      :2:     1/1
institution                             :2:     1/1
license                                 :2:     1/1
naming_authority                        :2:     1/1
no_blanks_in_id                         :2:     1/1
processing_level                        :2:     1/1
project                                 :2:     1/1
publisher_email                         :2:     1/1
publisher_name                          :2:     1/1
publisher_url                           :2:     1/1
source                                  :2:     1/1
standard_name_vocabulary                :2:     1/1
time_coverage_duration                  :2:     0/1
time_coverage_end                       :2:     1/1
time_coverage_extents_match             :2:     2/2
time_coverage_resolution                :2:     1/1
time_coverage_start                     :2:     1/1


--------------------------------------------------------------------------------
                  Reasoning for the failed tests given below:                   


Name                             Priority:     Score:Reasoning
--------------------------------------------------------------------------------
Conventions                            :3:     1/ 2 : Attr Conventions does not
                                                      contain 'ACDD-1.3'
varattr                                :3:    32/44 :  
    conductivity                       :3:     3/ 4 :  
        coverage_content_type          :3:     0/ 1 : Var conductivity missing
                                                      attr coverage_content_type
    density                            :3:     3/ 4 :  
        coverage_content_type          :3:     0/ 1 : Var density missing attr
                                                      coverage_content_type
    profile_id                         :3:     1/ 4 :  
        coverage_content_type          :3:     0/ 1 : Var profile_id missing
                                                      attr coverage_content_type
        var_std_name                   :3:     1/ 2 : Var profile_id missing
                                                      attr standard_name
        var_units                      :3:     0/ 1 : Var profile_id missing
                                                      attr units
    salinity                           :3:     3/ 4 :  
        coverage_content_type          :3:     0/ 1 : Var salinity missing attr
                                                      coverage_content_type
    segment_id                         :3:     1/ 4 :  
        coverage_content_type          :3:     0/ 1 : Var segment_id missing
                                                      attr coverage_content_type
        var_std_name                   :3:     1/ 2 : Var segment_id missing
                                                      attr standard_name
        var_units                      :3:     0/ 1 : Var segment_id missing
                                                      attr units
    temperature                        :3:     3/ 4 :  
        coverage_content_type          :3:     0/ 1 : Var temperature missing
                                                      attr coverage_content_type
    u                                  :3:     3/ 4 :  
        coverage_content_type          :3:     0/ 1 : Var u missing attr
                                                      coverage_content_type
    v                                  :3:     3/ 4 :  
        coverage_content_type          :3:     0/ 1 : Var v missing attr
                                                      coverage_content_type
date_created_is_iso                    :2:     0/ 1 : Datetime provided is not
                                                      in a valid ISO 8601 format
date_issued_is_iso                     :2:     0/ 1 : Datetime provided is not
                                                      in a valid ISO 8601 format
date_modified_is_iso                   :2:     0/ 1 : Datetime provided is not
                                                      in a valid ISO 8601 format
geospatial_bounds                      :2:     0/ 1 : Attr geospatial_bounds not
                                                      present
geospatial_bounds_crs                  :2:     0/ 1 : Attr geospatial_bounds_crs
                                                      not present
geospatial_bounds_vertical_crs         :2:     0/ 1 : Attr
                                                      geospatial_bounds_vertical
                                                      _crs not present
geospatial_lat_extents_match           :2:     1/ 2 : Data for possible latitude
                                                      variables ({u'lat':
                                                      9.9692099683868702e+36,
                                                      u'lat_uv':
                                                      9.9692099683868702e+36})
                                                      did not match
                                                      geospatial_lat_max value
                                                      (34.85172)
geospatial_lon_extents_match           :2:     1/ 2 : Data for possible
                                                      longitude variables
                                                      ({u'lon_uv':
                                                      9.9692099683868702e+36,
                                                      u'lon':
                                                      9.9692099683868702e+36})
                                                      did not match
                                                      geospatial_lon_max value
                                                      (-120.78092)
geospatial_vertical_extents_match      :2:     0/ 2 : geospatial_vertical_min !=
                                                      min(depth) values, 1.1 !=
                                                      0.11,
                                                      geospatial_vertical_max !=
                                                      max(depth) values, 1.1 !=
                                                      58.9
time_coverage_duration                 :2:     0/ 1 : Attr
                                                      time_coverage_duration not
                                                      present
creator_institution                    :1:     0/ 1 : Attr creator_institution
                                                      not present
creator_type                           :1:     0/ 2 : Attr creator_type not
                                                      present
date_metadata_modified                 :1:     0/ 1 : Attr
                                                      date_metadata_modified not
                                                      present
instrument                             :1:     0/ 1 : Attr instrument not
                                                      present
instrument_vocabulary                  :1:     0/ 1 : Attr instrument_vocabulary
                                                      not present
metadata_link                          :1:     0/ 1 : Attr metadata_link is
                                                      empty or completely
                                                      whitespace
metadata_link_valid                    :1:     0/ 1 : Metadata URL should
                                                      include http:// or
                                                      https://
platform                               :1:     0/ 1 : Attr platform not present
platform_vocabulary                    :1:     0/ 1 : Attr platform_vocabulary
                                                      not present
product_version                        :1:     0/ 1 : Attr product_version not
                                                      present
program                                :1:     0/ 1 : Attr program not present
publisher_institution                  :1:     0/ 1 : Attr publisher_institution
                                                      not present
publisher_type                         :1:     0/ 2 : Attr publisher_type not
                                                      present
references                             :1:     0/ 1 : Attr references is empty
                                                      or completely whitespace

```

```
$ compliance-checker --test=cf compliance_checker/tests/data/examples/hycom_global.nc


--------------------------------------------------------------------------------
                    The dataset scored 113 out of 122 points                    
                              during the cf check                               
--------------------------------------------------------------------------------
                               Scoring Breakdown:                               


                                 High Priority                                  
--------------------------------------------------------------------------------
    Name                            :Priority: Score
§2.2 Valid netCDF data types            :3:     6/6
§2.4 Unique dimensions                  :3:     6/6
§3.1 Variable depth contains valid CF u :3:     3/3
§3.1 Variable lat contains valid CF uni :3:     3/3
§3.1 Variable lon contains valid CF uni :3:     3/3
§3.1 Variable time contains valid CF un :3:     3/3
§3.1 Variable water_u contains valid CF :3:     3/3
§3.1 Variable water_v contains valid CF :3:     3/3
§3.3 Variable time has valid standard_n :3:     0/1
§4 depth contains a valid axis          :3:     2/2
§4 lat contains a valid axis            :3:     2/2
§4 lon contains a valid axis            :3:     2/2
§4.1 Latitude variable lat has required :3:     1/1
§4.1 Longitude variable lon has require :3:     1/1
§4.3.1 depth is a valid vertical coordi :3:     1/2
§4.4 Time coordinate variable and attri :3:     2/2
§5.0 Variable water_u does not contain  :3:     4/4
§5.0 Variable water_v does not contain  :3:     4/4
§5.6 Grid Feature water_u is associated :3:     2/2
§5.6 Grid Feature water_v is associated :3:     2/2
§9.1 Dataset contains a valid featureTy :3:     1/1
§9.1 Feature Types are all the same     :3:     1/1


                                Medium Priority                                 
--------------------------------------------------------------------------------
    Name                            :Priority: Score
cell_methods                            :2:     0/0
§2.3 Naming Conventions for attributes  :2:    27/27
§2.3 Naming Conventions for dimensions  :2:     4/4
§2.3 Naming Conventions for variables   :2:     6/6
§2.3 Unique variable names              :2:     6/6
§2.4 Dimension Order                    :2:     6/6
§2.5.1 Fill Values should be outside th :2:     0/0
§2.6.1 Global Attribute Conventions inc :2:     0/1
§2.6.2 Recommended Attributes           :2:     0/3
§2.6.2 Recommended Global Attributes    :2:     1/2
§4.1 Latitude variable lat defines eith :2:     1/1
§4.1 Latitude variable lat uses recomme :2:     1/1
§4.1 Longitude variable lon defines eit :2:     1/1
§4.1 Longitude variable lon uses recomm :2:     1/1
§5.0 multidimensional coordinate lat sh :2:     1/1
§5.0 multidimensional coordinate lon sh :2:     1/1
§8.1 Packed Data defined by water_u con :2:     1/1
§8.1 Packed Data defined by water_u con :2:     0/1
§8.1 Packed Data defined by water_v con :2:     1/1
§8.1 Packed Data defined by water_v con :2:     0/1


--------------------------------------------------------------------------------
                  Reasoning for the failed tests given below:                   


Name                             Priority:     Score:Reasoning
--------------------------------------------------------------------------------
§3.3 Variable time has valid standard_n:3:     0/ 1 : variable time's attribute
                                                      standard_name must be a
                                                      non-empty string or it
                                                      should define a long_name
                                                      attribute.
§4.3.1 depth is a valid vertical coordi:3:     1/ 2 : vertical coordinates not
                                                      defining pressure must
                                                      include a positive
                                                      attribute that is either
                                                      'up' or 'down'
§2.6.1 Global Attribute Conventions inc:2:     0/ 1 : Conventions global
                                                      attribute does not contain
                                                      "CF-1.6"
§2.6.2 Recommended Attributes          :2:     0/ 3 : institution should be
                                                      defined, source should be
                                                      defined, references should
                                                      be defined
§2.6.2 Recommended Global Attributes   :2:     1/ 2 : global attribute history
                                                      should exist and be a non-
                                                      empty string
§8.1 Packed Data defined by water_u con:2:     0/ 1 : Attributes add_offset and
                                                      scale_factor are not of
                                                      type float or double.
§8.1 Packed Data defined by water_v con:2:     0/ 1 : Attributes add_offset and
                                                      scale_factor are not of
                                                      type float or double.



$ compliance-checker -d 35

Downloading cf-standard-names table version 35 from: http://cfconventions.org/Data/cf-standard-names/35/src/cf-standard-name-table.xml
```

## Installation

### Conda users

`compliance-checker` depends on many external C libraries,
so the easiest way to install it on MS-Windows/OS X/Linux is with `conda`.

```shell
$ conda install -c conda-forge compliance-checker
```

For more information on `conda` and installing the IOOS software stack see:

https://github.com/ioos/notebooks_demos/wiki/Installing-Conda-Python-with-the-IOOS-environment

### Pip users

When using pip the user must install the non-Python dependencies first.
Known C library dependencies include:
  - UDUNITS (2.x series)
  - HDF5
  - NetCDF4
  - libxml2/libxslt

Installation for these libraries will vary depending on your choice of operating system and method of installation
(i.e. binary packages versus compiling from source).
For more information on installing these libraries,
reference the documentation from the individual libraries.
(Check our [Ubuntu Install Guide](doc/ubuntu-install-guide.md) out.)

To install locally, set up a virtual environment (recommend using
[virtualenv-burrito](https://github.com/brainsik/virtualenv-burrito);

```
$ mkvirtualenv --no-site-packages compliance-checker
$ workon compliance-checker
```

The Python dependencies require several underlying system packages that most package managers should have.
See the [Installation](//github.com/ioos/compliance-checker/wiki/Installation) wiki page for more information.

```shell
$ pip install compliance-checker
```

## Usage (from Python code)

```python
from compliance_checker.runner import ComplianceChecker, CheckSuite

# Load all available checker classes
check_suite = CheckSuite()
check_suite.load_all_available_checkers()

# Run cf and adcc checks with normal strictness, verbose text format to stdout
return_value, errors = ComplianceChecker.run_checker(
                                '/path/or/url/to/your/dataset',
                                ['cf', 'acdd'], 0, 'normal', '-', 'text')
```


## Usage (Command Line)

```
compliance-checker <data-source> -t <test>[ <test>...]
```

The compliance checker command line tool will print out (to STDOUT) the test
results. The command line tool will also return 0 for a successful run and
non-0 for a failure.

## Available Test Suites

- [CF 1.6](http://cfconventions.org/)
- [ACDD 1.1, 1.3](http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/metadata/DataDiscoveryAttConvention.html)
- [IOOS](http://ioos.noaa.gov/data/contribute-data/)

## Available Test Suites as Plugins

- [GliderDAC](https://github.com/ioos/ioosngdac/wiki/NGDAC-NetCDF-File-Format-Version-2) - [link](https://github.com/ioos/cc-plugin-glider)

## Development

The compliance-checker is designed to be simple and hackable to edit existing compliance suites or introduce new ones. See the [Development](https://github.com/ioos/compliance-checker/wiki/Development) wiki page for more information.

#### Testing

Please run any new code through `flake8` and `pep8`. You can easily do this by using `pylama` (config is already setup in `pytest.ini`).

```bash
$ pip install pylama
$ py.test --pylama
```

Take a look at the failed tests and fix accordingly. Travis does not run with the `--pylama` flag so you MUSt do this yourself!


## Resource materials

[Compliance checker Webinar](
https://mmancusa.webex.com/mmancusa/ldr.php?RCID=e5e6fc5b6d218307f9eec863111e6034)


## Roadmap

- Improved text output (\#12)
- UGRID compliance (\#33)

## Contributors

- [Dave Foster](https://github.com/daf) &lt;dave@axiomdatascience.com&gt;
- [Dan Maher](https://github.com/danieljmaher) &lt;daniel.maher@rpsgroup.com&gt;
- [Luke Campbell](https://github.com/lukecampbell) &lt;luke.campbell@rpsgroup.com&gt;
- [Kyle Wilcox](https://github.com/kwilcox) &lt;kyle@axiomdatascience.com&gt;
- [Ben Adams](https://github.com/benjwadams) &lt;ben.adams@rpsgroup.com&gt;

And many more testers!

Portions of the CF checker are based on Michael Decker's work, http://repositories.iek.fz-juelich.de/hg/CFchecker/
