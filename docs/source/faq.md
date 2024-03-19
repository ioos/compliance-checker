# Frequently Asked Questions

## What is the Compliance Checker?

The IOOS Compliance Checker is a Python-based tool for data providers to check for completeness and community standard compliance of local or remote netCDF files against CF and ACDD file standards.
The Python module can be used as a command-line tool or as a library that can be integrated into other software.

You are currently viewing the web-based version of the Compliance Checker.
It was developed to enable a broader audience and improve accessibility for the checker.
With the web version,
providers can simply provide a link or upload their datasets and get the full suite of capabilities that Compliance Checker offers.

## What does the Compliance Checker check?

It currently supports the following sources and standards:
- [ACDD (1.1, 1.3)](https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3)
- [CF (1.6, 1.7)](http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.html)
- [IOOS (1.1, 1.2)](https://ioos.github.io/ioos-metadata/)
- [Glider DAC](https://github.com/ioos/ioosngdac/wiki/NetCDF-Specification)
- [NCEI (1.1, 2.0)](https://www.ncei.noaa.gov/data/oceans/ncei/formats/netcdf/v2.0/index.html)

## Can I test an ERDDAP dataset with the Compliance Checker?

Yes.
When testing an ERDDAP dataset,
please supply the URL to the dataset without a file extension.
For example,
to test this [Morro Bay dataset](https://standards.sensors.ioos.us/erddap/tabledap/morro-bay-bs1-met.html),
you should supply the URL like so:
"https://standards.sensors.ioos.us/erddap/tabledap/morro-bay-bs1-met".


## What version of the Compliance Checker is run on [compliance.ioos.us](https://compliance.ioos.us/index.html)?

This web site is using [version 5.0.0](https://pypi.org/project/compliance-checker/) of the Compliance Checker.

## Is there an API?

There sure is.
Check out the details on how to use the web [API here](https://github.com/ioos/compliance-checker-web/wiki/API).

## Where can I find more information about the Compliance Checker?

The Compliance Checker is completely open-source and available on [GitHub](https://github.com/ioos/compliance-checker).

## Disclaimer

The objective of the IOOS Compliance Checker is to check your file against
our interpretation of select dataset metadata standards to use as a
guideline in generating compliant files.  The compliance checker should
not be considered the authoritative source on whether your file is 100% "compliant".
Instead, we recommend that users use the results as a guide to work towards compliance.
