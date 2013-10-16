#!/usr/bin/env python

import argparse

from compliance_checker.suite import CheckSuite
from compliance_checker.acdd import ACDDCheck
from compliance_checker.cf import CFCheck

def main(ds_loc):
    cs = CheckSuite()
    out = cs.run(ds_loc, ACDDCheck(), CFCheck())

    import pprint
    pprint.pprint({str(type(k)): (v[0], [vv for vv in v[1] if vv.value != (0, 0)]) for k,v in out.iteritems()})

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_location', nargs='?')

    args = parser.parse_args()

    main(args.dataset_location)

