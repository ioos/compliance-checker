from __future__ import with_statement
from setuptools import setup, find_packages
import versioneer

def readme():
    with open('README.md') as f:
        return f.read()


reqs = [line.strip() for line in open('requirements.txt')]

setup(
    name                 = "compliance-checker",
    version              = versioneer.get_version(),
    description          = "Checks Datasets and SOS endpoints for standards compliance",
    long_description     = readme(),
    license              = 'Apache License 2.0',
    author               = "Dave Foster",
    author_email         = "dave@axiomdatascience.com",
    url                  = "https://github.com/ioos/compliance-checker",
    packages             = find_packages(),
    install_requires     = reqs,
    tests_require        = ['pytest'],
    classifiers          = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
    ],
    include_package_data = True,
    scripts=['cchecker.py'],

    # Note: Do not use colons in the entry-point keys. Python 3 reserves
    # portions of the key after a colon for special use.

    # Note: The entry point names are not used at all. All methods in the
    # compliance checker use class attributes to determine the checker's name
    # and version. But, an entry point must be defined for each plugin to be
    # loaded.

    entry_points         = {
        'console_scripts': [
            'compliance-checker = cchecker:main'
        ],
        'compliance_checker.suites': [
            'cf = compliance_checker.cf.cf:CFBaseCheck',
            'acdd-1.1 = compliance_checker.acdd:ACDD1_1Check',
            'acdd-1.3 = compliance_checker.acdd:ACDD1_3Check',
            'ioos_sos = compliance_checker.ioos:IOOSBaseSOSCheck',
            'ioos-0.1 = compliance_checker.ioos:IOOS0_1Check',
            'ioos-1.1 = compliance_checker.ioos:IOOS1_1Check',
        ]
    },
    package_data         = {
        'compliance_checker': ['data/*.xml', 'tests/data/*.nc', 'tests/data/*.cdl', 'tests/data/non-comp/*.cdl', 'data/templates/*.j2'],
    },
    cmdclass=versioneer.get_cmdclass(),
)
