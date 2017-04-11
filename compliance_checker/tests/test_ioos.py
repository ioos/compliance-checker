import unittest
from compliance_checker.suite import CheckSuite
from compliance_checker.runner import ComplianceChecker
import os
import sys

if sys.version_info < (3,0):
    import httpretty
else:
    httpretty = False

# TODO: Use inheritance to eliminate redundant code in test setup, etc
if httpretty:
    class TestIOOSSOSGetCapabilities(unittest.TestCase):

        def setUp(self):
            with open(os.path.join(os.path.dirname(__file__),
                        'data/http_mocks/ncsos_getcapabilities.xml')) as f:
                self.resp = f.read()
            # need to monkey patch checkers prior to running tests, or no checker
            # classes will show up
            CheckSuite().load_all_available_checkers()


        @httpretty.activate
        def test_retrieve_getcaps(self):
            """Method that simulates retrieving SOS GetCapabilities"""
            url = "http://data.oceansmap.com/thredds/sos/caricoos_ag/VIA/VIA.ncml"
            httpretty.register_uri(httpretty.GET, url, body=self.resp)
            # need to mock out the HEAD response so that compliance checker
            # recognizes this as some sort of XML doc instead of an OPeNDAP
            # source
            httpretty.register_uri(httpretty.HEAD, url, status=200,
                                    content_type='text/xml')
            ComplianceChecker.run_checker(url, ['ioos'], 1, 'normal')

    class TestIOOSSOSDescribeSensor(unittest.TestCase):

        def setUp(self):
            with open(os.path.join(os.path.dirname(__file__),
                        'data/http_mocks/ncsos_describesensor.xml')) as f:
                self.resp = f.read()
            # need to monkey patch checkers prior to running tests, or no checker
            # classes will show up
            CheckSuite().load_all_available_checkers()


        @httpretty.activate
        def test_retrieve_describesensor(self):
            """Method that simulates retrieving SOS DescribeSensor"""
            url = "http://data.oceansmap.com/thredds/sos/caricoos_ag/VIA/VIA.ncml?request=describesensor&service=sos&procedure=urn:ioos:station:ncsos:VIA&outputFormat=text/xml%3Bsubtype%3D%22sensorML/1.0.1/profiles/ioos_sos/1.0%22&version=1.0.0"
            httpretty.register_uri(httpretty.GET, url, body=self.resp)
            # need to mock out the HEAD response so that compliance checker
            # recognizes this as some sort of XML doc instead of an OPeNDAP
            # source
            httpretty.register_uri(httpretty.HEAD, url, status=200,
                                    content_type='text/xml')
            ComplianceChecker.run_checker(url, ['ioos'], 1, 'normal')
