#!/usr/bin/env python

import argparse

from compliance_checker.suite import CheckSuite
from compliance_checker.acdd import ACDDCheck

def main(ds_loc):
    cs = CheckSuite()
    out = cs.run(ds_loc, ACDDCheck())

    import pprint
    pprint.pprint({str(type(k)): v[0] for k,v in out.iteritems()})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_location', nargs='?')

    args = parser.parse_args()

    main(args.dataset_location)

