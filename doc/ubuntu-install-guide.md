# Installing Compliance-Checker on Ubuntu


## Installing SWIG on Ubuntu

```
sudo apt-get install swig
```

## Installing UDUnits on Ubuntu

```
wget 'ftp://ftp.unidata.ucar.edu/pub/udunits/udunits-2.1.24.tar.gz'
tar -zxvf udunits-2.1.24.tar.gz
cd udunits-2.1.24.tar.gz
./configure
make
sudo make install
cd ..
```

Linux infrequently updates the shared library cache (ldconfig). To force the cache to update:
```
sudo ldconfig -v
```

To ensure that UDUnits is properly installed and recognized by the operating system as a registered shared library:

```
ldconfig -p | grep udunits
```

You should see:

```
	libudunits2.so.0 (libc6,x86-64) => /usr/local/lib/libudunits2.so.0
	libudunits2.so (libc6,x86-64) => /usr/local/lib/libudunits2.so
```

## Installing Udunitspy

In your Python virtual environment that you set up for Compliance-Checker

```
pip install udunitspy
```

## Testing your UDUnitspy installation

In your Python virtual environment that you set up for Compliance-Checker

Run Python with the `python` command

```
Python 2.7.3 (default, Apr 10 2013, 06:20:15) 
[GCC 4.6.3] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> from udunitspy import Unit
>>> 
```

If you don't see an error then udunitspy was succesfully installed, and you can now finish the installation of the compliance checker.

## Installing lxml on Ubuntu

### Get the libxml2 and libxslt packages

```
sudo apt-get install libxml2-dev
sudo apt-get install libxslt1-dev
pip install lxml
```

## Installing Compliance Checker

```
pip install compliance-checker
```

## Testing your compliance checker installation:

Run python on your virtual environment and try:

```
from compliance_checker.runner import ComplianceCheckerCheckSuite
```

If it succeeds, then the Compliance Checker should be working correctly.


