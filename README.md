# IOOS Compliance Checker

The IOOS Compliance Checker is a Python tool to check local/remote datasets against a variety of compliance standards. It is primarily a command-line tool (tested on OSX/Linux) and can also be used as a library import.

It currently supports the following sources and standards:


| Standard                                                                                            | .nc/OPeNDAP             | SOS                             |
| --------------------------------------------------------------------------------------------------- | ----------------------- | ------------------------------- |
| [ACDD (1.1)](http://wiki.esipfed.org/index.php/Attribute_Convention_for_Data_Discovery_%28ACDD%29)  | Complete                | -                               |
| IOOS Asset Concept                                                                                  | -                       | GetCapabilities, DescribeSensor |
| [CF (1.6)](http://cf-convention.github.io/1.6.html)                                                 | Partial (chs 2-5)       | -                               |

### Concepts & Terminology

Each compliance standard is executed by a Check Suite, which functions similar to a Python standard Unit Test. A Check Suite runs one or more checks against a dataset, returning a list of Results which are then aggregated into a summary.

Each Result has a (# passed / # total) score, a weight (HIGH/MEDIUM/LOW), a computer-readable name, an optional list of human-readable messages, and optionally a list of child Results.

A single score is then calculated by aggregating on the names, then multiplying the score by the weight and summing them together.

The computer-readable name field controls how Results are aggregated together - in order to prevent the overall score for a Check Suite varying on the number of variables, it is possible to *group* Results together via the name property. Grouped results will only add up to a single top-level entry.

For example, ...

See the Development section for more details on implementation.

### Usage

```
$ ./cchecker.py --help
usage: cchecker.py [-h] [--test {acdd,cf,ioos} [{acdd,cf,ioos} ...]]
                   [--criteria [{lenient,normal,strict}]] [--verbose]
                   dataset_location

positional arguments:
  dataset_location      Defines the location of the dataset to be checked.

optional arguments:
  -h, --help            show this help message and exit
  --test {acdd,cf,ioos} [{acdd,cf,ioos} ...], -t {acdd,cf,ioos} [{acdd,cf,ioos} ...], --test= {acdd,cf,ioos} [{acdd,cf,ioos} ...], -t= {acdd,cf,ioos} [{acdd,cf,ioos} ...]
                        Select the Checks you want to perform. Either All, CF,
                        or ACDD. Defaults to All.
  --criteria [{lenient,normal,strict}], -c [{lenient,normal,strict}]
                        Define the criteria for the checks. Either Strict,
                        Normal, or Lenient. Defaults to Normal.
  --verbose, -v         Increase Output Verbosity
```

```
$ ./cchecker.py --test=acdd test-data/ru07-20130824T170228_rt0.nc
Running Compliance Checker on the dataset from: test-data/ru07-20130824T170228_rt0.nc


-------------------------------------------------------
   The dataset scored 95 out of 149 required points
            during the acdd check
      This test has passed under normal critera
-------------------------------------------------------

$ ./cchecker.py -v --test=acdd test-data/ru07-20130824T170228_rt0.nc
Running Compliance Checker on the dataset from: test-data/ru07-20130824T170228_rt0.nc

-------------------------------------------------------
The following tests failed:
----High priority tests failed-----
    Name                            :Priority: Score
varattr                                 :3:    69/120
----Medium priority tests failed-----
    Name                            :Priority: Score
acknowledgement                         :2:     0/1
cdm_data_type                           :2:     0/1
time_coverage_duration                  :2:     0/1
```

### Installation

To install locally, set up a virtual environment (recommend using virtualenv-burrito):

```
$ mkvirtualenv --no-site-packages compliance-checker
$ workon compliance-checker
```

Get compliance checker source:

```
$ git clone git@github.com:ioos/compliance-checker.git
$ cd compliance-checker
```

Install dependencies (you may need C dependencies for netCDF-python):

```
$ pip install numpy
$ pip install -r requirements.txt
```

### Development

The compliance-checker is designed to be simple and hackable to edit existing compliance suites or introduce new ones. See the [Development](https://github.com/ioos/compliance-checker/wiki/Development) wiki page for more information.

### Roadmap

- Complete CF 1.6 checks
- Improve text output
- Improve integration usage with other Python applications
 

