from __future__ import with_statement
import sys

from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

from compliance_checker import __version__

def readme():
    with open('README.md') as f:
        return f.read()

reqs = [line.strip() for line in open('requirements.txt')]

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True
    def run_tests(self):
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)

setup(name                 = "compliance-checker",
    version              = __version__,
    description          = "Checks Datasets and SOS endpoints for standards compliance",
    long_description     = readme(),
    license              = 'Apache License 2.0',
    author               = "Dave Foster",
    author_email         = "dfoster@asascience.com",
    url                  = "https://github.com/asascience-open/compliance-checker",
    packages             = find_packages(),
    install_requires     = reqs,
    tests_require        = ['pytest'],
    cmdclass             = {'test': PyTest},
    classifiers          = [
            'Development Status :: 3 - Alpha',
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
        ]
    },
    package_data         = {
        'compliance_checker':['data/*.json', 'data/*.xml', 'tests/data/*.nc', 'tests/data/non-comp/*.nc'],
    }
)

