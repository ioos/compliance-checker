from __future__ import with_statement
from setuptools import setup, find_packages
from compliance_checker import __version__


def readme():
    with open('README.md') as f:
        return f.read()

reqs = [line.strip() for line in open('requirements.txt')]

setup(
    name                 = "compliance-checker",
    version              = __version__,
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
    entry_points         = {
        'console_scripts': [
            'compliance-checker = cchecker:main'
        ],
        'compliance_checker.suites': [
            'cf = compliance_checker.cf.cf:CFBaseCheck',
            'acdd = compliance_checker.acdd:ACDD1_1Check',
            'acdd:1.1 = compliance_checker.acdd:ACDD1_1Check',
            'acdd:1.3 = compliance_checker.acdd:ACDD1_3Check',
            'ioos = compliance_checker.ioos:IOOSBaseCheck',
        ]
    },
    package_data         = {
        'compliance_checker': ['data/*.xml', 'tests/data/*.nc', 'tests/data/non-comp/*.nc', 'data/templates/*.j2'],
    }
)
