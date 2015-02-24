#!/usr/bin/env python
'''
compliance_checker.glider_dac

Compliance Test Suite for the IOOS National Glider Data Assembly Center
https://github.com/ioos/ioosngdac/wiki
'''

from compliance_checker.base import BaseCheck, BaseNCCheck, Result

class GliderCheck(BaseNCCheck):
    @classmethod
    def beliefs(cls): 
        '''
        Not applicable for gliders
        '''
        return {}

    @classmethod
    def make_result(cls, level, score, out_of, name, messages):
        return Result(level, (score, out_of), name, messages)

    def setup(self, ds):
        pass

    def check_nothing(self, ds):
        level = BaseCheck.HIGH
        score = 16
        out_of = 20
        name = 'Check Units'
        messages = []
        return self.make_result(level, score, out_of, name, messages)

