import os

from mocket import mocketize
from mocket.mockhttp import Entry

from compliance_checker.runner import ComplianceChecker
from compliance_checker.suite import CheckSuite


# TODO: Use inheritance to eliminate redundant code in test setup, etc
class TestIOOSSOSGetCapabilities:
    def setup_method(self):
        with open(
            os.path.join(
                os.path.dirname(__file__),
                "data/http_mocks/ncsos_getcapabilities.xml",
            ),
        ) as f:
            self.resp = f.read()
        # need to monkey patch checkers prior to running tests, or no checker
        # classes will show up
        CheckSuite().load_all_available_checkers()

    @mocketize
    def test_retrieve_getcaps(self):
        """Method that simulates retrieving SOS GetCapabilities"""
        url = "http://data.oceansmap.com/thredds/sos/caricoos_ag/VIA/VIA.ncml"
        Entry.single_register(
            method=Entry.GET,
            uri=url,
            body=self.resp,
            headers={"content-type": "text/xml"},
        )
        Entry.single_register(
            Entry.HEAD,
            url,
            headers={"content-type": "text/xml"},
            body="HTTP/1.1 200",
        )
        ComplianceChecker.run_checker(url, ["ioos_sos"], 1, "normal")


class TestIOOSSOSDescribeSensor:
    def setup_method(self):
        with open(
            os.path.join(
                os.path.dirname(__file__),
                "data/http_mocks/ncsos_describesensor.xml",
            ),
        ) as f:
            self.resp = f.read()
        # need to monkey patch checkers prior to running tests, or no checker
        # classes will show up
        CheckSuite().load_all_available_checkers()

    @mocketize
    def test_retrieve_describesensor(self):
        """Method that simulates retrieving SOS DescribeSensor"""
        url = (
            "http://data.oceansmap.com/thredds/sos/caricoos_ag/VIA/VIA.ncml?"
            "request=describesensor"
            "&service=sos"
            "&procedure=urn:ioos:station:ncsos:VIA"
            "&outputFormat=text/xml%3Bsubtype%3D%22sensorML/1.0.1/profiles/ioos_sos/1.0%22"
            "&version=1.0.0"
        )
        Entry.single_register(
            method=Entry.GET,
            uri=url,
            body=self.resp,
            headers={"content-type": "text/xml"},
        )
        Entry.single_register(
            method=Entry.HEAD,
            uri=url,
            body="HTTP/1.1 200",
            headers={"content-type": "text/xml"},
        )
        # need to mock out the HEAD response so that compliance checker
        # recognizes this as some sort of XML doc instead of an OPeNDAP
        # source
        ComplianceChecker.run_checker(url, ["ioos_sos"], 1, "normal")
