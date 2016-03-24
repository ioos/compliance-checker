# IOOS Compliance Checker

[![Build Status](https://travis-ci.org/ioos/compliance-checker.svg)](https://travis-ci.org/ioos/compliance-checker)

The IOOS Compliance Checker is a Python tool to check local/remote datasets against a variety of compliance standards. It is primarily a command-line tool (tested on OSX/Linux) and can also be used as a library import.

It currently supports the following sources and standards:

| Standard                                                                                             | Source                                                            | .nc/OPeNDAP | SOS                             |
| ---------------------------------------------------------------------------------------------------- | -----------                                                       | ------      | ------------------------------- |
| [ACDD (1.1)](http://wiki.esipfed.org/index.php/Attribute_Convention_for_Data_Discovery_%28ACDD%29)   | Built-in                                                          | X           | -                               |
| [CF (1.6)](http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html)                       | Built-in                                                          | X           | -                               |
| IOOS Asset Concept                                                                                   | Built-in                                                          | -           | GetCapabilities, DescribeSensor |
| [Glider DAC](https://github.com/ioos/ioosngdac/wiki/NGDAC-NetCDF-File-Format-Version-2)              | [ioos/cc-plugin-glider](https://github.com/ioos/cc-plugin-glider) | X           | -                               |

## Concepts & Terminology

Each compliance standard is executed by a Check Suite, which functions similar to a Python standard Unit Test. A Check Suite runs one or more checks against a dataset, returning a list of Results which are then aggregated into a summary.

Each Result has a (# passed / # total) score, a weight (HIGH/MEDIUM/LOW), a computer-readable name, an optional list of human-readable messages, and optionally a list of child Results.

A single score is then calculated by aggregating on the names, then multiplying the score by the weight and summing them together.

The computer-readable name field controls how Results are aggregated together - in order to prevent the overall score for a Check Suite varying on the number of variables, it is possible to *group* Results together via the name property. Grouped results will only add up to a single top-level entry.

See the [Development](//github.com/ioos/compliance-checker/wiki/Development) wiki page for more details on implementation.

## Usage (command line)

The compliance-checker can work against local files (.nc files, .xml files of SOS GetCapabilities/DescribeSensor requests) or against remote URLs (OPeNDAP data URLs, SOS GetCapabilities/DescribeSensor URLs).

> **WARNING** The CF/ACDD checks **will access data**, so if using a remote OPenDAP URL, please be sure the size is reasonable!

```
usage: compliance-checker [-h] [--test {gliderdac,acdd,cf,ioos}]
                   [--criteria [{lenient,normal,strict}]] [--verbose]
                   [-f {stdout,html}] [-o OUTPUT]
                   dataset_location [dataset_location ...]

positional arguments:
  dataset_location      Defines the location of the dataset to be checked.

optional arguments:
  -h, --help            show this help message and exit
  --test {gliderdac,acdd,cf,ioos}, -t {gliderdac,acdd,cf,ioos}, --test= {gliderdac,acdd,cf,ioos}, -t= {gliderdac,acdd,cf,ioos}
                        Select the Checks you want to perform.
  --criteria [{lenient,normal,strict}], -c [{lenient,normal,strict}]
                        Define the criteria for the checks. Either Strict,
                        Normal, or Lenient. Defaults to Normal.
  --verbose, -v         Increase output. May be specified up to three times.
  -f {stdout,html}, --format {stdout,html}
                        Output format
  -o OUTPUT, --output OUTPUT
                        Output filename
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
```

## Installation

### Non-Python dependencies

A number of the Python modules used by the compliance checker such as NumPy
require linking against external C or Fortran libraries.  If using an
 installation method for Python dependencies such as pip below, these external
 libraries will need to be installed beforehand prior to the installation of
 the Python dependencies.
Known C library dependencies include:
  - NetCDF4
  - HDF5
  - UDUNITS (2.x series)
  - libxml2/libxslt

Installation for these libraries will vary depending on your choice of
operating system and method of installation (i.e. binary packages versus
compiling from source).  For more information on installing these libraries,
reference the documentation from the individual libraries.

### Python dependencies with `pip`

To install locally, set up a virtual environment (recommend using
[virtualenv-burrito](https://github.com/brainsik/virtualenv-burrito);
however, do not use `virtualenv-burrito` if you are using Anaconda python.
[Instructions for Anaconda](#anaconda) appear below):

```
$ mkvirtualenv --no-site-packages compliance-checker
$ workon compliance-checker
```

The Python dependencies require several underlying system packages that most package managers should have. See the [Installation](//github.com/ioos/compliance-checker/wiki/Installation) wiki page for more information.

Install dependencies, numpy must be installed on its own:

```
$ pip install numpy
$ pip install compliance-checker
```

###<a name="anaconda"></a> Anaconda on linux

Assuming you have a standard Anaconda installation, run the following:

1. From the compliance-checker top-level directory create the virtual
environment (stay as far away from virtualenv's burrito as you possibly
can):    
	```     
conda create --file requirements.txt -n compliance_checker    
source activate compliance_checker    
conda install compliance-checker    
	```    
     
2. Test the installation using the following commands:     
	```    
	compliance-checker --help    
	compliance-checker --test=acdd compliance_checker/tests/data/ru07-20130824T170228_rt0.nc    
	compliance-checker --test=cf compliance_checker/tests/data/sss20140107.v2.0cap.nc     
	```    

### MS-Windows

Before [creating the virtual environment and installing compliance-checker](https://github.com/ioos/compliance-checker#-anaconda-on-linux), the channel locations should be updated in the Anaconda configuration file `.condarc`; without that the `conda` command fails to find the required packages. If that file does not exist it should be created in the Anaconda root directory with the following text (the comments with "#" prefix are optional):

```
# channel locations. These override conda defaults, i.e., conda will
# search *only* the channels listed here, in the order given. Use "defaults" to
# automatically include all default channels. Non-url channels will be
# interpreted as binstar usernames (this can be changed by modifying the
# channel_alias key; see below).
channels:
  - defaults
  - IOOS
  
```

If after that conda still cannot find the required packages, the `binstar` package should be installed:

```
conda install binstar
```

### Ubuntu

Installing on Ubuntu? Check out the [Ubuntu Install Guide](doc/ubuntu-install-guide.md)

## Usage (from Python code)

```python
from compliance_checker.runner import ComplianceCheckerCheckSuite

cs = ComplianceCheckerCheckSuite()
dataset = cs.load_dataset("/path/or/url/to/your/dataset")
groups = cs.run(dataset, 'acdd')
scores = groups['acdd']
```


## Usage (Command Line)

```
compliance-checker <data-source> -t [cf|acdd|ioos|gliderdac]
```

The compliance checker command line tool will print out (to STDOUT) the test
results. The command line tool will also return 0 for a successful run and
non-0 for a failure.

## Available Test Suites

- [CF 1.6](http://cfconventions.org/)
- [ACDD](http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/metadata/DataDiscoveryAttConvention.html)
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


To run tests for both Python 2 and 3, install `tox` globally via
```
pip install tox
```

then run

`tox -c tox.ini`

to run the tests

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

And many more testers!    

Portions of the CF checker are based on Michael Decker's work, http://repositories.iek.fz-juelich.de/hg/CFchecker/    
 
<!--
   15  vim --version
   17  python --version
   23  git --version
   27  wget http://09c8d0b2229f813c1b93-c95ac804525aac4b6dba79b00b39d1d3.r79.cf1.rackcdn.com/Anaconda-2.1.0-Linux-x86_64.sh
   33  bash Anaconda-2.1.0-Linux-x86_64.sh
   67  wget https://github.com/spf13/hugo/releases/download/v0.12/hugo_0.12_linux_amd64.tar.gz
   70  git clone ssh://git@github.com/duncombe/dotfiles.git
   82  git clone ssh://git@github.com/duncombe/centos-instance.git
  106  git config --global user.name "Christopher Duncombe Rae"
  108  git config --global user.email "christopher.duncombe.rae@noaa.gov"
  278  git remote add origin ssh://git@github.com/duncombe/PROGS.git 
  344  git clone ssh://git@github.com/duncombe/system-test.git
  348  git remote set-url upstream ssh://git@github.com/ioos/system-test.git
  449  git submodule add ssh://git@github.com/duncombe/bash-git-prompt.git
  749  git clone ssh://git@github.com/duncombe/unidata-python-workshop.git
  750  git remote add upstream  ssh://git@github.com/Unidata/unidata-python-workshop.git
  783  git clone https://github.com/Unidata/unidata-python-workshop
  784  conda config --add channels https://conda.binstar.org/rsignell
  785  conda config --add channels https://conda.binstar.org/Unidata
  791  conda create -n workshop python=2 numpy matplotlib cartopy ipython ipython-notebook     netcdf4 owslib pyudl networkx basemap
  836  compliance-checker --test=acdd compliance_checker/tests/data/test-data/ru07-20130824T170228_rt0.nc
  837  compliance-checker --test=acdd compliance_checker/tests/data/ru07-20130824T170228_rt0.nc 
  838  compliance-checker --test=acdd compliance_checker/tests/data/ru07-20130824T170228_rt0.nc | less -S
  839  compliance-checker --test=acdd test-data/ru07-20130824T170228_rt0.nc
  
-->


