# IOOS Compliance Checker

The IOOS Compliance Checker is a Python tool to check local/remote datasets against a variety of compliance standards. It is primarily a command-line tool (tested on OSX/Linux) and can also be used as a library import.

It currently supports the following sources and standards:


| Standard                                                                                            | .nc/OPeNDAP             | SOS                             |
| --------------------------------------------------------------------------------------------------- | ----------------------- | ------------------------------- |
| [ACDD (1.1)](http://wiki.esipfed.org/index.php/Attribute_Convention_for_Data_Discovery_%28ACDD%29)  | Complete                | -                               |
| IOOS Asset Concept                                                                                  | -                       | GetCapabilities, DescribeSensor |
| [CF (1.6)](http://cfconventions.org/Data/cf-convetions/cf-conventions-1.6/build/cf-conventions.html)                                                 | Complete                | -                               |

### Concepts & Terminology

Each compliance standard is executed by a Check Suite, which functions similar to a Python standard Unit Test. A Check Suite runs one or more checks against a dataset, returning a list of Results which are then aggregated into a summary.

Each Result has a (# passed / # total) score, a weight (HIGH/MEDIUM/LOW), a computer-readable name, an optional list of human-readable messages, and optionally a list of child Results.

A single score is then calculated by aggregating on the names, then multiplying the score by the weight and summing them together.

The computer-readable name field controls how Results are aggregated together - in order to prevent the overall score for a Check Suite varying on the number of variables, it is possible to *group* Results together via the name property. Grouped results will only add up to a single top-level entry.

See the [Development](//github.com/ioos/compliance-checker/wiki/Development) wiki page for more details on implementation.

### Usage (command line)

The compliance-checker can work against local files (.nc files, .xml files of SOS GetCapabilities/DescribeSensor requests) or against remote URLs (OPeNDAP data URLs, SOS GetCapabilities/DescribeSensor URLs).

> **WARNING** The CF/ACDD checks **will access data**, so if using a remote OPenDAP URL, please be sure the size is reasonable!

```
$ compliance-checker --help
usage: compliance-checker [-h] [--test {acdd,cf,ioos} [{acdd,cf,ioos} ...]]
                          [--criteria [{lenient,normal,strict}]] [--verbose]
                          dataset_location

positional arguments:
  dataset_location      Defines the location of the dataset to be checked.

optional arguments:
  -h, --help            show this help message and exit
  --test {acdd,cf,ioos} [{acdd,cf,ioos} ...], -t {acdd,cf,ioos} [{acdd,cf,ioos} ...], --test= {acdd,cf,ioos} [{acdd,cf,ioos} ...], -t= {acdd,cf,ioos} [{acdd,cf,ioos} ...]
                        Select the Checks you want to perform. Either all
                        (default), cf, ioos, or acdd.
  --criteria [{lenient,normal,strict}], -c [{lenient,normal,strict}]
                        Define the criteria for the checks. Either Strict,
                        Normal, or Lenient. Defaults to Normal.
  --verbose, -v         Increase output. May be specified up to three times.
```

```
$ compliance-checker --test=acdd test-data/ru07-20130824T170228_rt0.nc
Running Compliance Checker on the dataset from: test-data/ru07-20130824T170228_rt0.nc


-------------------------------------------------------
   The dataset scored 95 out of 149 required points
            during the acdd check
      This test has passed under normal critera
-------------------------------------------------------

$ compliance-checker --test=cf sss20140107.v2.0cap.nc
Running Compliance Checker on the dataset from: sss20140107.v2.0cap.nc


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

### Installation

To install locally, set up a virtual environment (recommend using [virtualenv-burrito](https://github.com/brainsik/virtualenv-burrito)):

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

#### MS-Windows

IOOS Compliance Checker cannot be installed on MS-Windows systems.  The required `udunitspy` and `swig` packages are either broken or unavailable for that OS (2014-08-18).

#### Ubuntu

Installing on Ubuntu? Check out the [Ubuntu Install Guide](doc/ubuntu-install-guide.md)

### Usage (from Python code)

```python
from compliance_checker.runner import ComplianceCheckerCheckSuite

cs = ComplianceCheckerCheckSuite()
dataset = cs.load_dataset("/path/or/url/to/your/dataset")
groups = cs.run(dataset, 'acdd')
scores = groups['acdd']
```

### Development

The compliance-checker is designed to be simple and hackable to edit existing compliance suites or introduce new ones. See the [Development](https://github.com/ioos/compliance-checker/wiki/Development) wiki page for more information.

### Roadmap

- Improved text output (#12)
- UGRID compliance (#33)

### Contributors

- Dave Foster <dfoster@asascience.com>
- Dan Maher <dmaher@asascience.com>
- Luke Campbell <lcampbell@asascience.com>

And many more testers!

Portions of the CF checker are based on Michael Decker's work, http://repositories.iek.fz-juelich.de/hg/CFchecker/
