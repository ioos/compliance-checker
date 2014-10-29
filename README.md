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

### Installation

To install locally, set up a virtual environment (recommend using
[virtualenv-burrito](https://github.com/brainsik/virtualenv-burrito);
however, do not use virtualenv-burrito if you are using Anaconda python.
Instructions for Anaconda appear below):

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

#### Anaconda on linux

Assuming you have a standard Anaconda installation, run the following:

1. From the compliance-checker top-level directory create the virtual
environment (stay as far away from virtualenv's burrito as you possibly
can):
	```
	 conda create --file requirements.txt -n compliance-checker
	 conda install compliance-checker
	source activate compliance_checker
	```
2. Test the installation using the following commands:
	```
  compliance-checker --help
	 compliance-checker --test=acdd compliance_checker/tests/data/ru07-20130824T170228_rt0.nc
  compliance-checker --test=cf compliance_checker/tests/data/sss20140107.v2.0cap.nc 
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

<!--
    1  ls -ltr
    2  ls /data
    3  mount 
    4  man man
    5  ntpq -pn
    6  ls -ltr /etc
    7  ls -ltr /etc/ntp
    8  cat /etc/ntp.conf 
    9  w
   10  ls -ltr
   11  which wget
   12  git status
   13  which git
   14  which vi
   15  vim --version
   16  which python
   17  python --version
   18  ps ax | grep ntp
   19  less /etc/init.d/ntpd
   20  service ntpd status
   21  date
   22  vim --version | less
   23  git --version
   24  which curl
   25  wget http://cf1.rackcdn.com/Anaconda-2.1.0-Linux-x86_64.sh
   26  wget http://rackcdn.com/Anaconda-2.1.0-Linux-x86_64.sh
   27  wget http://09c8d0b2229f813c1b93-c95ac804525aac4b6dba79b00b39d1d3.r79.cf1.rackcdn.com/Anaconda-2.1.0-Linux-x86_64.sh
   28  ls -altr
   29  mkdir software
   30  cd software/
   31  mv -i ../Anaconda-2.1.0-Linux-x86_64.sh  .
   32  wget http://09c8d0b2229f813c1b93-c95ac804525aac4b6dba79b00b39d1d3.r79.cf1.rackcdn.com/Anaconda-2.1.0-Linux-x86_64.sh
   33  bash Anaconda-2.1.0-Linux-x86_64.sh
   34  ls -ltr
   35  rm Anaconda-2.1.0-Linux-x86_64.sh
   36  mv Anaconda-2.1.0-Linux-x86_64.sh.1 Anaconda-2.1.0-Linux-x86_64.sh
   37  bash Anaconda-2.1.0-Linux-x86_64.sh 
   38  ls 
   39  cd ..
   40  df -h
   41  ls -altr /data
   42  ls -ltr
   43  mount 
   44  mkdir  IOOS
   45  ls -latr
   46  vi .bashrc
   47  vim .bashrc
   48  ls -ltr
   49  ls
   50  git clone ssh://git@github.com/duncombe/dotfiles.git
   51  ssh git@github.com
   52  ssh -v git@github.com
   53  ls -ltr .ssh
   54  cd .ssh
   55  ssh-keygen 
   56  ls -ltr
   57  cat id_rsa.pub 
   58  ssh-rsa AAAAB3NzaC1yc2EAAAABIwAAAQEAlqDGWPzUuJ6W+7tqp2HHChRX9rNFN1u94hbdKYGw5XPr1JT01dSf54e4G8IeIAYDwf5+qqQJeMQjDFHvX1kY0g4l+T7SCUfy2qZGCGMAf//wsp4aJONaQ3scboHGIqIj88CqFTaMsnYjQGNby8+t82024o3fQ/mhU70bsq/EPI0Zqj1VPKSOWoL2YKgQsKtmE1P1YocQFZvrb0uITQlatgarwXcgBQQorFZ3fVyyolm53VcO93sJM8xJYliOzlpa1g9cAJop8QXfGk6LtyIxJzmgf7ZYaYj4SDwV2fLmqyf1NVb8eoK3n1qQ73fTWLdt+wLrll+mP3aVUVv9OAVCmw== cduncomberae@ip-10-191-198-218
   59  ssh -v git@github.com
   60  ssh git@github.com
   61  w
   62  last
   63  ls -ltr
   64  ls
   65  ls IOOS/
   66  cd software/
   67  wget https://github.com/spf13/hugo/releases/download/v0.12/hugo_0.12_linux_amd64.tar.gz
   68  ls
   69  cd ..
   70  git clone ssh://git@github.com/duncombe/dotfiles.git
   71  cd dotfiles/
   72  ls -ltr
   73  cd duncombe/
   74  ls
   75  ls -altr
   76  history
   77  cd ..
   78  ls -ltr
   79  cd ..
   80  ls
   81  git clone ssh://git@github.com/duncombe/centos-install.git
   82  git clone ssh://git@github.com/duncombe/centos-instance.git
   83  cd centos-instance/
   84  ls -ltr
   85  git config
   86  git config --get-all
   87  git config -l
   88  cd ../dotfiles
   89  ls
   90  git status
   91  cd duncombe/
   92  ls -altr
   93  ls -altr ~
   94  vim bash_profile 
   95  ls -ltr
   96  history >> commands-file
   97  fg
   98  ls ~/.gitconfig
   99  fg
  100  git config --get user.name
  101  git config --set user.name="Christopher Duncombe Rae"
  102  git config --global user.name="Christopher Duncombe Rae"
  103  git config --get user.name
  104  git config -l
  105  git config --global add user.name="Christopher Duncombe Rae"
  106  git config --global user.name "Christopher Duncombe Rae"
  107  git config -l
  108  git config --global user.email "christopher.duncombe.rae@noaa.gov"
  109  git config -l
  110  fg
  111  ls -ltr
  112  git status
  113  cd ../../centos-instance/
  114  ls
  115  ls -ltr
  116  git status
  117  git add setup.sh 
  118  git commit -m "Begin a set script."
  119  git config --global -l
  120  git remote -v
  121  git push
  122  git commit -m "Begin a set script."
  123  git fetch
  124  git status
  125  git log --oneline --graph
  126  git merge 
  127  git merge origin/master
  128  git push 
  129  cd ..
  130  ls
  131  ls IOO
  132  ls IOOS/
  133   
  134  ls -ltr
  135  w
  136  ls -ltr
  137  cd anaconda/
  138  ls
  139  cd ..
  140  df -h
  141  df .
  142  df -h .
  143  du -s .
  144  du -sh .
  145  cd dotfiles/
  146  ls -ltr
  147  cd duncombe
  148  ls -ltr
  149  man ln
  150  pwd
  151  cd ..
  152  vi setup-duncombe.sh
  153  w
  154  ls -ltr
  155  git status
  156  vi .gitignore
  157  mv -i duncombe/commands-file ~
  158  git add setup-duncombe.sh 
  159  w
  160  uptime
  161  ls -ltr
  162  history
  163  cd dotfiles/
  164  ls
  165  git status
  166  git commit -m "Creating a setup script." 
  167  git push
  168  cd ..
  169  ls
  170  vi test-unpushed
  171  bash test-unpushed 
  172  man find
  173  cd dotfiles/
  174  ls -ltr
  175  cat setup-duncombe.sh 
  176  ln -sb ~/dotfiles/duncombe/bashrc ~/.bashrc
  177  ls -altr ~
  178  ln -sb ~/dotfiles/duncombe/bash_logout .bash_logout
  179  ls -altr ~
  180  xterm
  181  cat setup-duncombe.sh 
  182  ln -sb ~/dotfiles/duncombe/bashrc_custom .bashrc_custom
  183  ln -sb ~/dotfiles/duncombe/vimrc .vimrc
  184  ls -ltr
  185  vi setup-duncombe.sh 
  186  bash setup-duncombe.sh 
  187  vi setup-duncombe.sh 
  188  man test
  189  vi setup-duncombe.sh 
  190  bash setup-duncombe.sh 
  191  vi setup-duncombe.sh 
  192  bash setup-duncombe.sh 
  193  vi setup-duncombe.sh 
  194  bash setup-duncombe.sh 
  195  ls -altr ~
  196  git status
  197  which git
  198  git status
  199  rehahs
  200  rehash
  201  hash -r
  202  git status
  203  ls -altr
  204  rm .bash_logout .bashrc_custom .vimrc 
  205  ls -ltr
  206  git status
  207  git add setup-duncombe.sh 
  208  git commit -m "Make the script do something."
  209  git push
  210  git config --global push.default current
  211  ls -ltr
  212  cd ..
  213  ls
  214  cd dotfiles/
  215  git status
  216  git status -uno
  217  git status -u no
  218  git fetch --dry-run
  219  man find
  220  git help remote
  221  bash
  222  source ~/.bashrc
  223  git help remote
  224  man ssh
  225  ls -ltr
  226  vi ~/.bashrc
  227  vi ~/.bashrc_custom 
  228  grep PATH ~/.*
  229  grep -il PATH ~/.*
  230  vi ~/.bash_profile
  231  echo $PATH
  232  export PATH=~/PROGS:$PATH
  233  ls ~/PROGS/
  234  chmod +x ~/PROGS/test-unpushed 
  235  (cd ~; test-unpushed )
  236  vi ~/PROGS/test-unpushed 
  237  (cd ~; test-unpushed )
  238  pushd ~/PROGS/
  239  ls -altr
  240  git status
  241  git status -u no
  242  git status
  243  git help status
  244  git status -s
  245  cd ..
  246  test-unpushed 
  247  git help update
  248  git help remote 
  249  gitk 
  250  git status
  251  test-unpushed 
  252  cd dotfiles/
  253  git status
  254  git commit -a -m "duncombe/bash_profile updated." 
  255  git push
  256  ls -ltr
  257  rm mooky~ mooky hooky 
  258  ls -altr
  259  ls -tlr
  260  vi test-unpushed 
  261  git --version
  262  find . -name ".git" -type d -exec dirname \{\} \; 
  263  find . -iname "*.swp" -print
  264  vi test-unpushed 
  265  bash test-unpushed 
  266  vi test-unpushed 
  267  bash test-unpushed 
  268  vi test-unpushed 
  269  bash test-unpushed 
  270  vi test-unpushed 
  271  bash test-unpushed 
  272  vi test-unpushed 
  273  git init PROGS
  274  cd PROGS/
  275  mv -i ../test-unpushed .
  276  git add test-unpushed 
  277  git commit -m "Created test-unpushed."
  278  git remote add origin ssh://git@github.com/duncombe/PROGS.git 
  279  git push origin master
  280  vi test-unpushed 
  281  git add test-unpushed 
  282  git commit -m "test-unpushed updated." 
  283  git remote -v
  284  git remote update origin
  285  git status
  286  git log --oneline --graph
  287  git status
  288  git remote -v
  289  git branch -ra
  290  git fetch origin
  291  git branch -ra
  292  git log --oneline --graph
  293  git status -u no
  294  git show-branch master
  295  git show-branch 
  296  git diff
  297  git diff origin/master
  298  git log --oneline --graph
  299  git rev-parse
  300  git rev-parse --help
  301  git rev-parse @
  302  git rev-parse @{u}
  303  git merge-base --help
  304  git branch --help
  305  git merge-base @ @{u}
  306  git branch --set-upstream master origin/master
  307  git help branch
  308  git status
  309  git status -s
  310  git status -u no
  311  vi test-unpushed 
  312  ./test-unpushed 
  313  git diff 
  314  git diff origin
  315  git diff origin/master
  316  git diff origin/master master
  317  git commit -a -m "test-unpushed updated." 
  318  git push
  319  ./test-unpushed 
  320  cd ..
  321  test-unpushed 
  322  ls -altr
  323  ls 
  324  ls -altr dotfiles/
  325  ls -ltr
  326  ls -altr dotfiles/
  327  touch mooky
  328  ln -sb hooky mooky
  329  ls -ltr
  330  touch hooky
  331  ln -sb hooky mooky
  332  ls -ltr
  333  ls -altr
  334  bash
  335  bash
  336  ls
  337  cd IOOS
  338  ls
  339  ls ../PROGS/
  340  less ../PROGS/test-unpushed 
  341  ls ../dotfiles/
  342  ls -altr
  343  ls ../.ssh/
  344  git clone ssh://git@github.com/duncombe/system-test.git
  345  cd system-test/
  346  git remote add upstream ssh://git@github.com/ioos/systemtest.git
  347  git fetch --all
  348  git remote set-url upstream ssh://git@github.com/ioos/system-test.git
  349  git fetch --all
  350  git status
  351  git pull --rebase upstream master
  352  vi README.md 
  353  git checkout upstream/master README.md
  354  ls -ltr
  355  git rebase --continue
  356  git status
  357  git rebase --continue
  358  git add  README.md 
  359  git status
  360  git rebase --continue
  361  git rebase --skip
  362  vi README.md 
  363  git status
  364  git add README.md 
  365  git rebase --continue
  366  git commit -m "Fixed conflicts."
  367  git rebase --continue
  368  git rebase --skip
  369  git rebase --abort
  370  git status
  371  git fetch --all
  372  git reset --hard upstream/master
  373  git log --oneline --graph
  374  git status
  375  git pull origin master
  376  vi README.md 
  377  git branch -ra
  378  git reset --hard origin/master
  379  git status
  380  git checkout -b master-upstream upstream/master
  381  git branch rename master messedup-master
  382  git branch -m master messedup-master
  383  git status
  384  git uptime
  385  uptime
  386  w
  387  pwd
  388  pushd ~/dotfiles/
  389  ls
  390  ls -altr
  391  git status
  392  git fetch --all
  393  git merge 
  394  ls -altr
  395  ls -altr duncombe/
  396  man history
  397  history -a
  398  man history
  399  ls -altr
  400  git status
  401  git fetch --all
  402  git merge 
  403  ls -altr
  404  ls -altr duncombe/
  405  man history
  406  history -a
  407  git clone ssh://git@github.com/duncombe/compliance-checker.git
  408  cd compliance-checker/
  409  git remote add upstream ssh://git@github.com/ioos/compliance-checker.git
  410  git fetch --all
  411  git status
  412  git pull --rebase upstream master
  413  git status
  414  git push origin master
  415  vi .bashrc
  416  ls -altr 
  417  man ln
  418  ls -ltr
  419  .  /home/cduncomberae/.bash-git-prompt/gitprompt.sh
  420  pwd
  421  whatami
  422  uname -a
  423  PROGS/whatami
  424  ls ~/PROGS
  425  pwd
  426  cd dotfiles
  427  ls
  428  cd ..
  429  ls
  430  cd PROGS
  431  git clone ssh://git@github.com/duncombe/WHATAMI.git
  432  cd WHATAMI/
  433  ls
  434  ./whatami 
  435  gmt
  436  GMT
  437  ntpq -pn
  438  ls /home/
  439  date
  440  git status
  441  pwd
  442  cd dotfiles/
  443  ls
  444  cd duncombe/
  445  ls
  446  ls -alt
  447  ls -altr
  448  cd ..
  449  git submodule add ssh://git@github.com/duncombe/bash-git-prompt.git
  450  ls -ltr
  451  git status
  452  cat .gitmodules 
  453  git commit -a -m "Added bash-git-prompt as a submodule."
  454  vi setup-duncombe.sh 
  455  ls -latr
  456  . setup-duncombe.sh 
  457  vi setup-duncombe.sh 
  458  . setup-duncombe.sh 
  459  vi setup-duncombe.sh 
  460  . setup-duncombe.sh 
  461  vi setup-duncombe.sh 
  462  ln -sb /home/cduncomberae/dotfiles/bash-git-prompt /home/cduncomberae/.bash-git-prompt
  463  ls -ltr
  464  vu ~/.bashrc
  465  grep gitprom ~/.* 
  466  .  /home/cduncomberae/.bash-git-prompt/gitprompt.sh
  467  git status
  468  git log --oneline --graph
  469  pwd
  470  git add setup-duncombe.sh 
  471  git commit -m "Updated setup-duncombe.sh."
  472  git push origin master
  473  w
  474  ls -ltr
  475  ls centos-instance/
  476  less centos-instance/setup.sh 
  477  PROGS/test-unpushed 
  478  ls
  479  PROGS/test-unpushed | less
  480  PROGS/test-unpushed 2>&1 | less
  481  pwd
  482  cd IOOS/system-test/
  483  ls -altr
  484  git status
  485  git remote -v
  486  git checkout master
  487  git status
  488  git branch -m master-upstream master
  489  git log --oneline --graph
  490  git log -p --oneline --graph
  491  git log -p --oneline --graph | less -r
  492  git log -p --oneline --graph | less -R
  493  git log--oneline --graph | less -R
  494  git log --oneline --graph | less -R
  495  git log --oneline --graph 
  496  git config --global -l
  497  git config -l
  498  git log --oneline --graph 
  499  git config -l
  500  git log --oneline --graph 
  501  git branch -ra
  502  gitk
  503  git push origin master
  504  git push origin master -f
  505  git status
  506  git help track
  507  git help branch
  508  git branch -u origin/master 
  509  git branch -l 
  510  git branch -l upstream
  511  git branch -ra 
  512  git checkout upstream
  513  git remote -v
  514  git status
  515  git push origin messedup-master
  516  git log --oneline --graph
  517  git status
  518  git diff master
  519  git diff master upstream
  520  git branch -vv
  521  git branch -vv -ra
  522  git branch -vv -ra | less -R
  523  git branch -vv -ra
  524  git branch -u upstream/master upstream
  525  git branch -vv -ra
  526  git branch --unset-upstream messedup-master
  527  git branch -vv -ra | cut --80
  528  cut --help
  529  git branch -vv -ra | cut -c -80
  530  ls -ltr
  531  git checkout master
  532  ls -ltr
  533  echo $PATH
  534  cd 
  535  PROGS/test-unpushed 
  536  PROGS/test-unpushed  | cut -c -80
  537  PROGS/WHATAMI/whatami 
  538  which lsb-release
  539  man which
  540  man which > /dev/null
  541  which lsb-release > /dev/null
  542  which lsb-release @> /dev/null
  543  which lsb-release 2> /dev/null
  544  which cat 2> /dev/null
  545  [ -x `which lsb-release` ] && echo hello
  546  [ -x "`which lsb-release`" ] && echo hello
  547  [ -x "`which lsb-release 2> /dev/null`" ] && echo hello
  548  [ -x "`which cat 2> /dev/null`" ] && echo hello
  549  PROGS/WHATAMI/whatami 
  550  PROGS/test-unpushed  | cut -c -80
  551  (cd PROGS; PROGS/test-unpushed  | cut -c -80 )
  552  (cd PROGS; ~/PROGS/test-unpushed  | cut -c -80 )
  553  git log --stat 
  554  pd
  555  cd PROGS/
  556  git log --stat 
  557  git log --stat  WHATAMI/
  558  cd WHATAMI/
  559  git log --stat  
  560  git log origin/master..HEAD
  561  git diff origin/master..HEAD 
  562  git log origin/master..HEAD --oneline
  563  git log origin/master..HEAD --oneline | wc -l
  564  git log origin/master..HEAD --oneline 
  565  git log origin/master..HEAD --oneline | cat
  566  git --nopager log origin/master..HEAD --oneline | cat
  567  git --no-pager log origin/master..HEAD --oneline | cat
  568  git help alias
  569  git alias
  570  git --help
  571  git config --help
  572  git config alias.ahead="git config --help --nopager log origin/master..HEAD --oneline"
  573  git config alias.ahead "git config --help --nopager log origin/master..HEAD --oneline"
  574  git ahead
  575  git config alias.ahead "git config alias.ahead "git config --help --nopager log origin/master..HEAD --oneline" config --help --nopager log origin/master..HEAD --oneline"
  576  git config alias.ahead "\!git config --help --nopager log origin/master..HEAD --oneline"
  577  git ahead
  578  git config alias.ahead "\!git config --nopager log origin/master..HEAD --oneline"
  579  git ahead
  580  git config alias.ahead '!git config --nopager log origin/master..HEAD --oneline'
  581  git ahead
  582  git config alias.ahead '!git --nopager log origin/master..HEAD --oneline'
  583  git ahead
  584  history | less
  585  git config alias.ahead '!git --no-pager log origin/master..HEAD --oneline'
  586  git ahead
  587  cd ..
  588  git ahead
  589  git status
  590  cd WHATAMI/
  591  git status
  592  git push
  593  git config
  594  git config -l
  595  git config --global -l
  596  git config --global user.name "christopher.duncombe.rae@noaa.gov"
  597  git config --global user.email "christopher.duncombe.rae@noaa.gov"
  598  git config --global user.name "Christopher Duncombe Rae"
  599  git ahead
  600  git help status
  601  ls -ltr
  602  git status
  603  git status -s
  604  git ahead
  605  cd ..
  606  ls 
  607  git sttus
  608  git status
  609  git add test-unpushed 
  610  git commit -m "Updated test-unpushed."
  611  git push
  612  git submodule --help
  613  git submodule status
  614  git submodule --help
  615  history | less
  616  git status
  617  git push
  618  w
  619  cd dotfiles/
  620  vi setup-duncombe.sh 
  621  pwd
  622  ls -ltr
  623  cd 
  624  git status
  625  ls -altr
  626  vi .gitconfig 
  627  mv -i .gitconfig gitconfig.old
  628  cp -i ~/dotfiles/duncombe/gitconfig .gitconfig
  629  ls -ltra 
  630  vi .gitconfig
  631  ls -altr
  632  vi .gitconfig
  633  vi .gitconfig_extras
  634  vi .gitconfig
  635  vi .gitconfig_extras
  636  PROGS/test-unpushed 
  637  vi PROGS/test-unpushed 
  638  PROGS/test-unpushed 
  639  vi PROGS/test-unpushed 
  640  PROGS/test-unpushed 
  641  cd dotfiles/
  642  ls -ltr
  643  cd duncombe/
  644  grep git *
  645  less gitconfig 
  646  git config --global core.pager less -R
  647  git  config -l
  648  git config --global core.pager "less -R"
  649  git  config -l
  650  pd
  651  pushd ~
  652  PROGS/test-unpushed 
  653  PROGS/test-unpushed  | less
  654  vi PROGS/test-unpushed 
  655  PROGS/test-unpushed  | less
  656  git status
  657  pwd
  658  cd PROGS/
  659  git status
  660  git remote -v
  661  git add test-unpushed 
  662  git commit -m "Update test-unpushed to output stderr to stdout."
  663  git push
  664  ls -ltr
  665  ls WHATAMI/
  666  ls
  667  git help branch 
  668  ls -altr
  669  cd WHATAMI/
  670  ls -ltr
  671  vi whatami 
  672  git add whatami 
  673  git commit -m "Stop stderr output from lsb-release."
  674  git fetch 
  675  git status
  676  ls -altr ~
  677  pd
  678  pushd 
  679  ls 
  680  ls -altr
  681  diff gitconfig ~/.gitconfig 
  682  git log origin/master..HEAD --oneline 
  683  git help log
  684  git log --help origin/master..HEAD --oneline 
  685  git log --help
  686  git log --help origin/master..HEAD --oneline | cat
  687  git log origin/master..HEAD --oneline | cat
  688  git log origin/master..HEAD --oneline 
  689  vi gitconfig 
  690  ls -la ~
  691  vi gitconfig 
  692  cd ..
  693  ls
  694  cd duncombe/
  695  ls -ltr
  696  git mv gitconfig gitconfig_extras
  697  ls -ltr
  698  cp -i ~/.gitconfig .
  699  mv -i .gitconfig gitconfig
  700  vi gitconfig
  701  git status
  702  git add gitconfig_extras ../setup-duncombe.sh gitconfig
  703  git commit -m "Including gitconfig files in a careful way."
  704  git status
  705  git push
  706  cd ..
  707  ln -sb duncombe/gitconfig_extras ~/.gitconfig_extras
  708  ln -sb `/bin/pwd`/duncombe/gitconfig_extras ~/.gitconfig_extras
  709  vi gitconfig
  710  vi duncombe/gitconfig
  711  git status
  712  git ahead
  713  vi setup-duncombe.sh 
  714  git ahead
  715  git status
  716  git add setup-duncombe.sh 
  717  git status
  718  git ahead
  719  git commit -m "Update setup-duncombe.sh."
  720  git ahead
  721  git push
  722  vu duncombe/gitconfig_extras 
  723  vi duncombe/gitconfig_extras 
  724  vi duncombe/bashrc
  725  vi duncombe/bashrc_custom 
  726  vi ~/.bashrc_custom 
  727  git ahead
  728  git status -s
  729  ls -ltr
  730  git status
  731  git add duncombe/bashrc_custom
  732  git commit -m "Updated duncombe/bashrc_custom."
  733  git push
  734  ifconfig
  735  source .bashrc_local 
  736  echo $PATH
  737  cd IOOS/
  738  sl
  739  ls
  740  cd compliance-checker/
  741  git status
  742  git branch -ra
  743  git branch -b my-test
  744  git checkout -b my-test
  745  conda 
  746  conda -V
  747  pwd
  748  pushd ..
  749  git clone ssh://git@github.com/duncombe/unidata-python-workshop.git
  750  git remote add upstream  ssh://git@github.com/Unidata/unidata-python-workshop.git
  751  cd unidata-python-workshop/
  752  git remote add upstream  ssh://git@github.com/Unidata/unidata-python-workshop.git
  753  git fetch --all
  754  git pull --rebase upstream master
  755  git push origin master
  756  git checkout my-workshop
  757  git merge master
  758  vi CompositeRadar.ipynb
  759  ipython trust CompositeRadar.ipynb 
  760  git status
  761  git add CompositeRadar.ipynb 
  762  git status
  763  git commit
  764  vi ~/.gitconfig
  765  vi ~/.gitconfig_extras
  766  git commit
  767  git push 
  768  locate activate
  769  find . -iname "*activate*" -print
  770  find . -iname "*source*" -print
  771  which source
  772  man source
  773  source activate workshop
  774  alias
  775  man source
  776  echo $sourcepath
  777  man source
  778  set 
  779  set sourcepath
  780  source activate workshop
  781  locate activate
  782  ls -l /home/cduncomberae/anaconda/bin/activate
  783  git clone https://github.com/Unidata/unidata-python-workshop
  784  conda config --add channels https://conda.binstar.org/rsignell
  785  conda config --add channels https://conda.binstar.org/Unidata
  786  ls -ltr
  787  rm -rf unidata-python-workshop/
  788  ls -ltr
  789  conda config --add channels https://conda.binstar.org/rsignell
  790  conda config --add channels https://conda.binstar.org/Unidata
  791  conda create -n workshop python=2 numpy matplotlib cartopy ipython ipython-notebook     netcdf4 owslib pyudl networkx basemap
  792  source activate workshop
  793  ipython notebook
  794  (cd ~; ~/PROGS/test-unpushed )
  795  git status
  796  (cd ~; ~/PROGS/test-unpushed )
  797  git ahead
  798  view ~/.gitconfig
  799  view ~/.gitconfig_extras 
  800  git status
  801  git fetch --all
  802  git branch -ra
  803  git ahead
  804  git branch my-workshop
  805  view .git/refs/heads/my-workshop 
  806  git lg
  807  ls -ltr
  808  pwd
  809  pd
  810  pushd
  811  ls -ltr
  812  git lg
  813  it lg
  814  git lg
  815  view requirements.txt 
  816  conda --help
  817  conda --help install
  818  conda install --help
  819  conda install --file requirements.txt 
  820  ls
  821  curl -sL https://raw.githubusercontent.com/brainsik/virtualenv-burrito/master/virtualenv-burrito.sh | $SHELL
  822  mkvirtualenv --no-site-packages compliance-checker
  823  source /home/cduncomberae/.venvburrito/startup.sh
  824  mkvirtualenv --no-site-packages compliance-checker
  825  pip install numpy
  826  conda install pip
  827  pip install numpy
  828  pip install compliance-checker
  829  compliance-checker --help
  830  compliance-checker --test=acdd test-data/ru07-20130824T170228_rt0.nc
  831  ls test-data
  832  ls
  833  ls compliance_checker/
  834  ls compliance_checker/tests/
  835  ls compliance_checker/tests/data/
  836  compliance-checker --test=acdd compliance_checker/tests/data/test-data/ru07-20130824T170228_rt0.nc
  837  compliance-checker --test=acdd compliance_checker/tests/data/ru07-20130824T170228_rt0.nc 
  838  compliance-checker --test=acdd compliance_checker/tests/data/ru07-20130824T170228_rt0.nc | less -S
  839  compliance-checker --test=acdd test-data/ru07-20130824T170228_rt0.nc
  840  mkvirtualenv --no-site-packages compliance-checker
  841  conda install --file requirements.txt 
  842  locate sss20140107
  843  compliance-checker --test=acdd compliance_checker/tests/data/sss20140107.v2.0cap.nc 
  844  compliance-checker --test=acdd compliance_checker/tests/data/sss20140107.v2.0cap.nc | less
  845  mkvirtualenv --no-site-packages compliance-checker
  846  ls -ltr
  847  cd ..
  848  ls -ltr
  849  mkdir trash
  850  mv -i compliance-checker/ trash/
  851  git clone ssh://git@github.com/duncombe/compliance-checker.git
  852  cd compliance-checker/
  853  source deactivate workshop
  854  source deactivate 
  855  deactivate workshop
  856  deactivate 
  857  pushd ../trash/compliance-checker/
  858  deactivate 
  859  source deactivate 
  860  source --help
  861  vi ~/.bashrc
  862  echo $PATH
  863  grep burri ~/.*
  864  grep -il burri ~/.*
  865  vi /home/cduncomberae/.bash_profile
  866  vi $HOME/.venvburrito/startup.sh
  867  orkon compliance_checker
  868  workon compliance_checker
  869  workon compliance-checker
  870  view /home/cduncomberae/.virtualenvs/compliance-checker
  871  history 
  872  history  | less
  873  cd IOOS/compliance-checker/
  874  ls -ltr
  875  cd 
  876  ls -ltr
  877  ls IOOS/trash/compliance-checker/
  878  ls -l IOOS/trash/compliance-checker/
  879  ls -altr 
  880  ls .virtualenvs/
  881  ls .virtualenvs/compliance-checker/
  882  ls .virtualenvs/compliance-checker/bin/
  883  vi .venvburrito/
  884  curl -sL https://raw.githubusercontent.com/brainsik/virtualenv-burrito/master/virtualenv-burrito.sh | less
  885  ls .virtualenvs/
  886  ls .virtualenvs/compliance-checker/
  887  ls .virtualenvs/compliance-checker/bin/
  888  ls .virtualenvs/compliance-checker/bin/python
  889  view .virtualenvs/compliance-checker/bin/python
  890  ls -l .virtualenvs/compliance-checker/bin/python2.7 
  891  which python
  892  ls -ltr
  893  ls -altr
  894  curl -sL https://raw.githubusercontent.com/brainsik/virtualenv-burrito/master/virtualenv-burrito.sh > virtualenv-burrito.sh
  895  view virtualenv-burrito.sh 
  896  w
  897  ls -ltr
  898  cd IOOS/
  899  ls
  900  ls -ltr
  901  cd compliance-checker/
  902  conda create --help
  903  conda create -f requirements.txt -n compliance-checker
  904  conda create --file requirements.txt -n compliance-checker
  905  source activate compliance_checker
  906  source activate compliance-checker
  907  compliance-checker --help
  908  echo $PATH
  909  ls /home/cduncomberae/anaconda/envs/compliance-checker/bin
  910  conda install numpy
  911  conda install compliance-checker
  912  compliance-checker --help
  913  compliance-checker --test=acdd compliance_checker/tests/data/ru07-20130824T170228_rt0.nc 
  914  compliance-checker --test=acdd compliance_checker/tests/data/ru07-20130824T170228_rt0.nc | less
  915  compliance-checker --test=cf compliance_checker/tests/data/sss20140107.v2.0cap.nc 
  916  git branch -ra
  917  git fetch --all
  918  git checkout -b vm-run 
  919  vi README.md 
  920  history >> README.md 
-->
