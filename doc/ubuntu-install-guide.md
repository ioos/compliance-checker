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


