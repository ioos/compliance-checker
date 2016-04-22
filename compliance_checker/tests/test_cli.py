#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Tests for command line output and parsing

'''

from unittest import TestCase
from compliance_checker.runner import ComplianceChecker, CheckSuite
import tempfile
import os


class TestCLI(TestCase):
    '''
    Tests various functions and aspects of the command line tool and runner
    '''

    def shortDescription(self):
        return None

    # override __str__ and __repr__ behavior to show a copy-pastable nosetest name for ion tests
    #  ion.module:TestClassName.test_function_name
    def __repr__(self):
        name = self.id()
        name = name.split('.')
        if name[0] not in ["ion", "pyon"]:
            return "%s (%s)" % (name[-1], '.'.join(name[:-1]))
        else:
            return "%s ( %s )" % (name[-1], '.'.join(name[:-2]) + ":" + '.'.join(name[-2:]))
    __str__ = __repr__

    def setUp(self):
        check_suite = CheckSuite()
        check_suite.load_all_available_checkers()

    def test_unicode_html(self):
        '''
        Tests that the checker is capable of producing HTML with unicode characters
        '''

        fd, path = tempfile.mkstemp()
        os.close(fd)
        self.addCleanup(os.remove, path)
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc='compliance_checker/tests/data/2dim-grid.nc',
            verbose=0,
            criteria='strict',
            checker_names=['acdd'],
            output_filename=path,
            output_format='html'
        )

        assert os.stat(path).st_size > 0

