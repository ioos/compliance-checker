# IOOS Compliance Checker

[![Build Status](https://travis-ci.org/ioos/compliance-checker.svg)](https://travis-ci.org/ioos/compliance-checker)
[![Build status](https://ci.appveyor.com/api/projects/status/lcc9co38pi6o45ho/branch/master?svg=true)](https://ci.appveyor.com/project/ocefpaf/compliance-checker/branch/master)

The IOOS Compliance Checker is a python based tool for data providers to check
for completeness and community standard compliance of local or remote
[netCDF](https://en.wikipedia.org/wiki/NetCDF) files against
[CF](https://en.wikipedia.org/wiki/NetCDF) and
[ACDD](http://wiki.esipfed.org/index.php/Attribute_Convention_for_Data_Discovery_1-3)
file standards. The python module can be used as a command-line tool or as a
library that can be integrated into other software.

A [web-based version](https://data.ioos.us/compliance/index.html) of the Compliance
Checker was developed to enable a broader audience and improve accessibility for the
checker. With the web version, providers can simply provide a link or upload their
datasets and get the full suite of capabilities that Compliance Checker offers.


It currently supports the following sources and standards:

| Standard                                                                                                                            | Source                                                            | .nc/OPeNDAP/.cdl | SOS                             |
| ----------------------------------------------------------------------------------------------------                                | -----------                                                       | ------           | ------------------------------- |
| [ACDD (1.1, 1.3)](http://wiki.esipfed.org/index.php/Attribute_Convention_for_Data_Discovery_1-3)                                    | Built-in                                                          | X                | -                               |
| [CF (1.6)](http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html)                                                      | Built-in                                                          | X                | -                               |
| IOOS SOS                                                                                                                            | Built-in                                                          | -                | GetCapabilities, DescribeSensor |
| [IOOS (1.1)](https://ioos.github.io/ioos-netcdf/ioos-netcdf-metadata-description-v1-1.html#ioos-netcdf-metadata-profile-attributes) | Built-in                                                          | X                | -                               |
| [Glider DAC](https://github.com/ioos/ioosngdac/wiki/NGDAC-NetCDF-File-Format-Version-2)                                             | [ioos/cc-plugin-glider](https://github.com/ioos/cc-plugin-glider) | X                | -                               |
| [NCEI (1.1, 2.0)](https://www.nodc.noaa.gov/data/formats/netcdf/v2.0/)                                                              | [ioos/cc-plugin-ncei](https://github.com/ioos/cc-plugin-ncei)     | X                | -                               |


## Advice to data providers

While the command-line version of this tool can be run in a loop, it is not necessary to check
every file if they are all created the same way. In short, this tool is not meant for
identifying bugs in your data processing stream. It is, however, intended to help you identify
your process procedure compliance to the standards.  If you change your processing procedure
for any reason it would be worth your while to run one file through the Compliance Checker to
insure you procedure change does not impact your file’s compliance.

If you feel you will need to run a batch of files through the Compliance Checker, please contact
the IOOS Program Office Operations Division for assistance.


# [The Compliance Checker Web Tool](https://data.ioos.us/compliance/)

The IOOS Compliance Checker front end companion.

[https://data.ioos.us/compliance/](https://data.ioos.us/compliance/)

Source Code is available on GitHub:

[https://github.com/ioos/compliance-checker-web](https://github.com/ioos/compliance-checker-web)

## Usage
Select the test you want to run from the dropdown menu. Then, either upload your dataset or provide a url to a
remote dataset (OPeNDAP) and click 'Submit'.

The output of the Compliance Checker will give you an overall score and a scoring breakdown.
You may download the Compliance Checker report as a text file by clicking the 'Download Report' button

![Compliance-Checker-Web](https://user-images.githubusercontent.com/5702672/30527267-b4bb136c-9bf4-11e7-8345-dd9b8e2e859f.png)

## API

In addition to a web-based front-end for the IOOS Compliance Checker project, an API is provided for
users interested in batch processing files hosted via OPeNDAP. Details on how to use the API are
available on the Compliance Checker Web [wiki page](https://github.com/ioos/compliance-checker-web/wiki/API).

Here are a couple examples:

**HTML Output**

https://data.ioos.us/compliance/api/run?report_format=html&test=acdd&url=http://sos.maracoos.org/stable/dodsC/hrecos/stationHRMARPH-agg.ncml

**JSON Output**

https://data.ioos.us/compliance/api/run?report_format=json&test=acdd&url=http://sos.maracoos.org/stable/dodsC/hrecos/stationHRMARPH-agg.ncml


# The Compliance Checker Command Line Tool


## Concepts & Terminology

Each compliance standard is executed by a Check Suite,
which functions similar to a Python standard Unit Test.
A Check Suite runs checks against a dataset based on a metadata standard,
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

## Installation

Check out the [Installation wiki](https://github.com/ioos/compliance-checker/wiki/Installation) for instructions on how to install.

## Command Line Usage

The compliance-checker can work against local files (`.nc` files, `.cdl`
metadata files, .xml files of SOS GetCapabilities/DescribeSensor requests)
or against remote URLs (OPeNDAP data URLs, SOS GetCapabilities/DescribeSensor URLs).

If you are aiming to check a netCDF-dump, also known as a CDL file, the file
must be named to end with a `.cdl` for the check-suite to be able to correctly
parse it's contents.

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
                        Output filename(s)
  -V, --version         Display the IOOS Compliance Checker version
                        information.
  -l, --list-tests      List the available tests
  -d DOWNLOAD_STANDARD_NAMES, --download-standard-names DOWNLOAD_STANDARD_NAMES
                        Specify a version of the cf standard name table to
                        download as packaged version
```

## Examples

### Check a local file against CF 1.6
```
$ compliance-checker --test=cf:1.6 compliance_checker/tests/data/examples/hycom_global.nc


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

```

### Check a remote file against ACDD 1.3

The remote dataset url is taken from the Data URL section of an OPeNDAP endpoint.

```shell
$ compliance-checker --test=acdd:1.3 "http://sos.maracoos.org/stable/dodsC/hrecos/stationHRMARPH-agg.ncml"
```

### Write results to text file

```shell
$ compliance-checker --test=acdd:1.3 --format=text --output=/tmp/report.txt compliance_checker/tests/data/examples/hycom_global.nc
```

### Write results to JSON file

```shell
$ compliance-checker --test=acdd:1.3 --format=json --output=/tmp/report.json compliance_checker/tests/data/examples/hycom_global.nc
```

### Write results to HTML file

```shell
$ compliance-checker --test=acdd:1.3 --format=html --output=/tmp/report.html compliance_checker/tests/data/examples/hycom_global.nc
```

### Download a particular CF standard names table for use in the test

**Note**
During the CF test, if a file has a particular version of the cf standard name table specified in the global attributes
(i.e. ```:standard_name_vocabulary = "CF Standard Name Table v30" ;```) that doesn't match the packaged version, it will
try to download the specified version. If it fails, it will fall back to packaged version.

```
$ compliance-checker -d 35

Downloading cf-standard-names table version 35 from: http://cfconventions.org/Data/cf-standard-names/35/src/cf-standard-name-table.xml
```



## Python Usage

If you are interested in incorporating the IOOS Compliance Checker into your own python projects, check out the following python code example:
```python
from compliance_checker.runner import ComplianceChecker, CheckSuite

# Load all available checker classes
check_suite = CheckSuite()
check_suite.load_all_available_checkers()

# Run cf and adcc checks
path = '/path/or/url/to/your/dataset'
checker_names = ['cf', 'acdd']
verbose = 0
criteria = 'normal'
output_filename = '/output/report.json'
output_format = 'json'
"""
Inputs to ComplianceChecker.run_checker

path            Dataset location (url or file)
checker_names   List of string names to run, should match keys of checkers dict (empty list means run all)
verbose         Verbosity of the output (0, 1, 2)
criteria        Determines failure (lenient, normal, strict)
output_filename Path to the file for output
output_format   Format of the output

@returns                If the tests failed (based on the criteria)
"""
return_value, errors = ComplianceChecker.run_checker(path,
                                                     checker_names,
                                                     verbose,
                                                     criteria,
                                                     output_filename=output_filename,
                                                     output_format=output_format)

# Open the JSON output and get the compliance scores
with open(output_filename, 'r') as fp:
    cc_data = json.load(fp)
    scored = cc_data[cc_test[0]]['scored_points']
    possible = cc_data[cc_test[0]]['possible_points']
    log.debug('CC Scored {} out of {} possible points'.format(scored, possible))
```

## Compliance Checker Plug-Ins

Separate Plug-ins have been developed to complement the master Compliance Checker tool with
specifications for preparing data to be submitted to different data assembly centers.
The version numbering of these plug-ins are not necessarily link to the version of the
master Compliance Checker, but they are all designed to run with the master Compliance Checker tool.

### Current Plug-in Releases:

- [GliderDAC](https://github.com/ioos/cc-plugin-glider/releases)

This is a checker for [GliderDAC](https://github.com/ioos/ioosngdac/wiki/NGDAC-NetCDF-File-Format-Version-2) files

- [NCEI](https://github.com/ioos/cc-plugin-ncei/releases) - [link](https://github.com/ioos/cc-plugin-ncei)

This is a checker for NCEI netCDF Templates [v1.1](https://www.nodc.noaa.gov/data/formats/netcdf/v1.1/) and [v2.0](https://www.nodc.noaa.gov/data/formats/netcdf/v2.0/) files.

These plug-ins must be installed separately but work on top of the base compliance checker software.

```
pip install cc-plugin-ncei
```

Check to see if it installed correctly, list the tests:

```
compliance-checker -l
```

You should see

```
 IOOS compliance checker available checker suites (code version):
 - ncei-grid (2.1.0)
 - ncei-grid:1.1 (2.1.0)
 - ncei-grid:2.0 (2.3.0)
 - ncei-grid:latest (2.1.0)
 - ncei-point (2.3.0)
 - ncei-point:1.1 (2.1.0)
 - ncei-point:2.0 (2.3.0)
 etc ....
```

Once installing the plug-in the usage is similar to the built in checkers.

### Examples of how to use the Plug-Ins

1. Run the NCEI Point check on a THREDDS endpoint

```python
compliance-checker -t ncei-point -v "https://data.nodc.noaa.gov/thredds/dodsC/testdata/mbiddle/GOLD_STANDARD_NETCDF/1.1/NODC_point_template_v1.1_2016-06-15_133710.844375.nc"
```

2. Run NCEI Trajectory Profile Orthogonal Check on local dataset

```python
compliance-checker -t ncei-trajectory-profile-orthogonal -v ~/data/sample-trajectory-profile.nc

```

3. Outputting JSON from a gridded file check
```
compliance-checker -t ncei-grid -f json -o ~/Documents/sample_grid_report.json ~/Documents/sample_grid_report.nc
```



## Contributors

- [Dave Foster](https://github.com/daf) &lt;dave@axiomdatascience.com&gt;
- [Dan Maher](https://github.com/danieljmaher) &lt;daniel.maher@gdit.com&gt;
- [Luke Campbell](https://github.com/lukecampbell) &lt;luke.campbell@gdit.com&gt;
- [Kyle Wilcox](https://github.com/kwilcox) &lt;kyle@axiomdatascience.com&gt;
- [Ben Adams](https://github.com/benjwadams) &lt;ben.adams@rpsgroup.com&gt;
- [Bob Fratantonio](https://github.com/bobfrat) &lt;robert.fratantonio@rpsgroup.com&gt;

And many more testers!

Portions of the CF checker are based on Michael Decker's work, http://repositories.iek.fz-juelich.de/hg/CFchecker/
