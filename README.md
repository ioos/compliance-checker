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
usage: compliance-checker [-h] [--test <test>[ <test>...]]
                   [--criteria [{lenient,normal,strict}]] [--verbose]
                   [-f {stdout,html}] [-o OUTPUT] [-d]
                   dataset_location [dataset_location ...]

positional arguments:
  dataset_location      Defines the location of the dataset to be checked.

optional arguments:
  -h, --help            show this help message and exit
  --test {gliderdac|acdd|cf|ioos[:version]}[ gliderdac|acdd|cf|ioos[:version]...] -t {gliderdac|acdd|cf|ioos[:version]}[ gliderdac|acdd|cf|ioos[:version]...] --test= {gliderdac|acdd|cf|ioos[:version]}[ gliderdac|acdd|cf|ioos[:version]...], -t= {gliderdac|acdd|cf|ioos[:version]}[ gliderdac|acdd|cf|ioos[:version]...]
                        Select the Checks you want to perform. Versions may be specified by using a colon followed by a version number.  Using a specification without a version number or with `:latest` will select the latest version of the standard.  Multiple tests may be specified by separating with a comma
  --criteria [{lenient,normal,strict}], -c [{lenient,normal,strict}]
                        Define the criteria for the checks. Either Strict,
                        Normal, or Lenient. Defaults to Normal.
  --verbose, -v         Increase output. May be specified up to three times.
  -f {stdout,html}, --format {stdout,html}
                        Output format
  -o OUTPUT, --output OUTPUT
                        Output filename
  -s check1 [-s check2] ..., --skip-checks check1 [--skip-checks check2] ...
                        Skips the any check functions contained within the
                        check suites selected which have the same function
                        name(s) as the string(s) specified.  May be specified
                        multiple times in order to skip multiple checks.
  -d DOWNLOAD_STANDARD_NAMES, --download-standard-names DOWNLOAD_STANDARD_NAMES
                        Specify a version of the cf standard name table to
                        download as packaged version
```

```
$ compliance-checker --test=acdd compliance_checker/tests/data/ru07-20130824T170228_rt0.nc
Running Compliance Checker on the dataset from: compliance_checker/tests/data/ru07-20130824T170228_rt0.nc


--------------------------------------------------------------------------------
                    The dataset scored 100 out of 188 points
                             during the acdd check
--------------------------------------------------------------------------------
                               Scoring Breakdown:


                                 High Priority
--------------------------------------------------------------------------------
    Name                            :Priority: Score
keywords                                :3:     1/1
summary                                 :3:     1/1
title                                   :3:     1/1
varattr                                 :3:    69/150


                                Medium Priority
--------------------------------------------------------------------------------
    Name                            :Priority: Score
acknowledgement                         :2:     0/1
cdm_data_type                           :2:     1/2
comment                                 :2:     1/1
creator_email                           :2:     1/1
creator_name                            :2:     1/1
creator_url                             :2:     1/1
date_created                            :2:     1/1
geospatial_lat_extents_match            :2:     1/2
geospatial_lat_max                      :2:     1/1
geospatial_lat_min                      :2:     1/1
geospatial_lon_extents_match            :2:     1/2
geospatial_lon_max                      :2:     1/1
geospatial_lon_min                      :2:     1/1
geospatial_vertical_extents_match       :2:     0/2
geospatial_vertical_max                 :2:     1/1
geospatial_vertical_min                 :2:     1/1
history                                 :2:     1/1
id                                      :2:     1/1
institution                             :2:     1/1
keywords_vocabulary                     :2:     1/1
license                                 :2:     1/1
naming_authority                        :2:     1/1
processing_level                        :2:     1/1
project                                 :2:     1/1
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
varattr                                :3:    69/150 :
    conductivity                       :3:     3/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var conductivity missing
                                                      attr coverage_content_type
    conductivity_qc                    :3:     2/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var conductivity_qc
                                                      missing attr
                                                      coverage_content_type
        var_units                      :3:     0/ 1 : Var conductivity_qc
                                                      missing attr units
    density                            :3:     3/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var density missing attr
                                                      coverage_content_type
    density_qc                         :3:     2/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var density_qc missing
                                                      attr coverage_content_type
        var_units                      :3:     0/ 1 : Var density_qc missing
                                                      attr units
    depth                              :3:     3/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var depth missing attr
                                                      coverage_content_type
    depth_qc                           :3:     2/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var depth_qc missing attr
                                                      coverage_content_type
        var_units                      :3:     0/ 1 : Var depth_qc missing attr
                                                      units
    instrument_ctd                     :3:     1/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var instrument_ctd missing
                                                      attr coverage_content_type
        var_std_name                   :3:     0/ 1 : Var instrument_ctd missing
                                                      attr standard_name
        var_units                      :3:     0/ 1 : Var instrument_ctd missing
                                                      attr units
    lat                                :3:     3/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var lat missing attr
                                                      coverage_content_type
    lat_qc                             :3:     2/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var lat_qc missing attr
                                                      coverage_content_type
        var_units                      :3:     0/ 1 : Var lat_qc missing attr
                                                      units
    lat_uv                             :3:     3/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var lat_uv missing attr
                                                      coverage_content_type
    lon                                :3:     3/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var lon missing attr
                                                      coverage_content_type
    lon_qc                             :3:     2/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var lon_qc missing attr
                                                      coverage_content_type
        var_units                      :3:     0/ 1 : Var lon_qc missing attr
                                                      units
    lon_uv                             :3:     3/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var lon_uv missing attr
                                                      coverage_content_type
    platform                           :3:     1/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var platform missing attr
                                                      coverage_content_type
        var_std_name                   :3:     0/ 1 : Var platform missing attr
                                                      standard_name
        var_units                      :3:     0/ 1 : Var platform missing attr
                                                      units
    pressure                           :3:     3/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var pressure missing attr
                                                      coverage_content_type
    pressure_qc                        :3:     2/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var pressure_qc missing
                                                      attr coverage_content_type
        var_units                      :3:     0/ 1 : Var pressure_qc missing
                                                      attr units
    profile_id                         :3:     1/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var profile_id missing
                                                      attr coverage_content_type
        var_std_name                   :3:     0/ 1 : Var profile_id missing
                                                      attr standard_name
        var_units                      :3:     0/ 1 : Var profile_id missing
                                                      attr units
    salinity                           :3:     3/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var salinity missing attr
                                                      coverage_content_type
    salinity_qc                        :3:     2/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var salinity_qc missing
                                                      attr coverage_content_type
        var_units                      :3:     0/ 1 : Var salinity_qc missing
                                                      attr units
    segment_id                         :3:     1/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var segment_id missing
                                                      attr coverage_content_type
        var_std_name                   :3:     0/ 1 : Var segment_id missing
                                                      attr standard_name
        var_units                      :3:     0/ 1 : Var segment_id missing
                                                      attr units
    temperature                        :3:     3/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var temperature missing
                                                      attr coverage_content_type
    temperature_qc                     :3:     2/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var temperature_qc missing
                                                      attr coverage_content_type
        var_units                      :3:     0/ 1 : Var temperature_qc missing
                                                      attr units
    time                               :3:     3/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var time missing attr
                                                      coverage_content_type
    time_qc                            :3:     2/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var time_qc missing attr
                                                      coverage_content_type
        var_units                      :3:     0/ 1 : Var time_qc missing attr
                                                      units
    time_uv                            :3:     3/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var time_uv missing attr
                                                      coverage_content_type
    trajectory                         :3:     1/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var trajectory missing
                                                      attr coverage_content_type
        var_std_name                   :3:     0/ 1 : Var trajectory missing
                                                      attr standard_name
        var_units                      :3:     0/ 1 : Var trajectory missing
                                                      attr units
    u                                  :3:     3/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var u missing attr
                                                      coverage_content_type
    u_qc                               :3:     2/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var u_qc missing attr
                                                      coverage_content_type
        var_units                      :3:     0/ 1 : Var u_qc missing attr
                                                      units
    v                                  :3:     3/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var v missing attr
                                                      coverage_content_type
    v_qc                               :3:     2/ 5 :
        var_coverage_content_type      :3:     0/ 2 : Var v_qc missing attr
                                                      coverage_content_type
        var_units                      :3:     0/ 1 : Var v_qc missing attr
                                                      units
acknowledgement                        :2:     0/ 1 : Attr acknowledgement not
                                                      present
cdm_data_type                          :2:     1/ 2 : Attr cdm_data_type present
                                                      but not in expected value
                                                      list (['vector', 'grid',
                                                      'textTable', 'tin',
                                                      'stereoModel', 'video'])
geospatial_lat_extents_match           :2:     1/ 2 : Data for possible latitude
                                                      variables ({u'lat':
                                                      9.969209968386869e+36,
                                                      u'lat_uv':
                                                      9.969209968386869e+36})
                                                      did not match
                                                      geospatial_lat_max value
                                                      (34.85172)
geospatial_lon_extents_match           :2:     1/ 2 : Data for possible
                                                      longitude variables
                                                      ({u'lon_uv':
                                                      9.969209968386869e+36,
                                                      u'lon':
                                                      9.969209968386869e+36})
                                                      did not match
                                                      geospatial_lon_max value
                                                      (-120.78092)
geospatial_vertical_extents_match      :2:     0/ 2 : Data for possible vertical
                                                      variables ({u'pressure':
                                                      0.10999999999999999,
                                                      u'depth':
                                                      0.10999999999999999}) did
                                                      not match
                                                      geospatial_vertical_min
                                                      value (1.1) Data for
                                                      possible vertical
                                                      variables ({u'pressure':
                                                      9.969209968386869e+36,
                                                      u'depth':
                                                      9.969209968386869e+36})
                                                      did not match
                                                      geospatial_vertical_max
                                                      value (589.0)
time_coverage_duration                 :2:     0/ 1 : Attr
                                                      time_coverage_duration not
                                                      present
contributor_role                       :1:     1/ 2 : Attr contributor_role
                                                      present but not in
                                                      expected value list
                                                      (['principalInvestigator',
                                                      'author'])



$ compliance-checker --test=cf compliance_checker/tests/data/sss20140107.v2.0cap.nc

Running Compliance Checker on the dataset from: compliance_checker/tests/data/sss20140107.v2.0cap.nc


--------------------------------------------------------------------------------
                     The dataset scored 12 out of 14 points
                              during the cf check
--------------------------------------------------------------------------------
                               Scoring Breakdown:


                                 High Priority
--------------------------------------------------------------------------------
    Name                            :Priority: Score
Variable names                          :3:     3/3
conventions                             :3:     0/1
data_types                              :3:     3/3
dimension_names                         :3:     3/3
units                                   :3:     0/1


                                Medium Priority
--------------------------------------------------------------------------------
    Name                            :Priority: Score
all_features_are_same_type              :2:     0/0
contiguous_ragged_array                 :2:     0/0
coordinate_type                         :2:     2/2
coordinates_and_metadata                :2:     0/0
feature_type                            :2:     0/0
incomplete_multidim_array               :2:     0/0
indexed_ragged_array                    :2:     0/0
missing_data                            :2:     0/0
orthogonal_multidim_array               :2:     0/0
var                                     :2:     1/1


--------------------------------------------------------------------------------
                  Reasoning for the failed tests given below:


Name                             Priority:     Score:Reasoning
--------------------------------------------------------------------------------
conventions                            :3:     0/ 1 : Conventions field is not
                                                      present
units                                  :3:     0/ 1 :
    sss_cap                            :3:     0/ 1 :
        known                          :3:     0/ 1 : unknown units type (PSU)



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
                                ['cf', 'acdd'], True, 'normal', '-', 'text')
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
- [IOOS](http://www.ioos.noaa.gov/data/contribute_data.html)

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

- Dave Foster <dave@axiomdatascience.com>
- Dan Maher <dmaher@asascience.com>
- Luke Campbell <lcampbell@asascience.com>
- [Kyle Wilcox](https://github.com/kwilcox) <kyle@axiomdatascience.com>
- Ben Adams <ben.adams@rpsgroup.com>

And many more testers!

Portions of the CF checker are based on Michael Decker's work, http://repositories.iek.fz-juelich.de/hg/CFchecker/
