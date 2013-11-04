"""
Compliance Checker suite runner
"""

import inspect
import itertools
from netCDF4 import Dataset
from lxml import etree as ET
from wicken.netcdf_dogma import NetCDFDogma
from compliance_checker.base import BaseCheck, fix_return_value, Result

class DSPair(object):
    '''
    Structure to hold a dataset and dogma pairing
    '''
    dataset = None 
    dogma   = None
    def __init__(self, ds, dogma):
        self.dataset = ds
        self.dogma = dogma

class CheckSuite(object):

    def _get_checks(self, checkclass):
        """
        Helper method to retreive check methods from a Checker class.

        The name of the methods in the Checker class should start with "check_" for this
        method to find them.
        """
        meths = inspect.getmembers(checkclass, inspect.ismethod)
        return [x[1] for x in meths if x[0].startswith("check_")]

    def _run_check(self, check_method, ds):
        val = check_method(ds)

        if isinstance(val, list):
            return [fix_return_value(v, check_method.im_func.func_name) for v in val]

        return [fix_return_value(val, check_method.im_func.func_name)]

    def run(self, dataset_location, *args):
        """
        Runs this CheckSuite on the dataset with all the passed Checker instances.

        Returns a dictionary mapping Checkers to their grouped scores.
        """

        ret_val = {}

        for a in args:

            ds = self.load_dataset(dataset_location, a.beliefs())

            a.setup(ds)
            checks = self._get_checks(a)

            vals = list(itertools.chain.from_iterable(map(lambda c: self._run_check(c, ds), checks)))

            # remove Nones? (aka skips)

            # transform to scores
            scores = self.scores(vals)

            ret_val[a] = scores

        return ret_val


    def load_dataset(self, ds_str, belief_map):
        """
        Helper method to load a dataset.
        """
        ds = Dataset(ds_str)
        data_object = NetCDFDogma('ds', belief_map, ds)

        return self.DSPair(ds, data_object)

    def scores(self, raw_scores):
        """
        Transforms raw scores from a single checker into a fully tallied and grouped scoreline.
        """
        grouped = self._group_raw(raw_scores)

        # score each top level item in groups
        weights = [x.weight for x in grouped]
        scores = [x.value for x in grouped]

        def score(acc, cur):
            # cur is ((x, x), weight)
            cs, weight = cur

            # acc is (score, total possible)
            cas, catotal = acc

            pts = 0
            if cs[1] == 0:
                return acc  # effectively a skip @TODO keep track of that?

            if cs[0] == cs[1]:
                return (cas + (1 * weight), catotal + weight)
            elif cs[0] == 0:
                return (cas, catotal + weight)

            # @TODO: could do cs[0] / cs[1] for exact pct points
            return (cas + (0.5 * weight), catotal + weight)

        cur_score = reduce(score, zip(scores, weights), (0, 0))

        return (cur_score, grouped)

    def _group_raw(self, raw_scores, cur=None, level=1):
        """
        Internal recursive method to group raw scores into a cascading score summary.

        Only top level items are tallied for scores.
        """
        #print "_group_raw", level, cur, len(raw_scores), raw_scores[0]

        def build_group(label=None, weight=None, value=None, sub=None):
            label  = label
            weight = weight
            value  = self._translate_value(value)
            sub    = sub or []

            return Result(weight=weight,
                          value=value,
                          name=label,
                          children=sub)

        def trim_groups(r):
            if isinstance(r.name, tuple) or isinstance(r.name, list):
                new_name = r.name[1:]
            else:
                new_name = []

            return Result(r.weight, r.value, new_name, r.msgs)

        # CHECK FOR TERMINAL CONDITION: all raw_scores.name are single length
        # @TODO could have a problem here with scalar name, but probably still works
        terminal = map(lambda x: len(x.name), raw_scores)
        if terminal == [0] * len(raw_scores):
            return []

        def group_func(r):
            """
            Slices off first element (if list/tuple) of classification or just returns it if scalar.
            """
            if isinstance(r.name, tuple) or isinstance(r.name, list):
                return r.name[0:1][0]
            return r.name

        grouped = itertools.groupby(sorted(raw_scores, key=group_func), key=group_func)
        ret_val = []

        for k, v in grouped:
            v = list(v)
            #print "group", k, level

            cv = self._group_raw(map(trim_groups, v), k, level+1)
            if len(cv):
                # if this node has children, max weight of children + sum of all the scores
                max_weight = max(map(lambda x: x.weight, cv))
                sum_scores = tuple(map(sum, zip(*(map(lambda x: x.value, cv)))))
                msgs = []
            else:
                max_weight = max(map(lambda x: x.weight, v))
                sum_scores = tuple(map(sum, zip(*(map(lambda x: self._translate_value(x.value), v)))))
                msgs = sum(map(lambda x: x.msgs, v), [])

            ret_val.append(Result(name=k, weight=max_weight, value=sum_scores, children=cv, msgs=msgs))

        return ret_val

    def _translate_value(self, val):
        """
        Turns shorthand True/False/None checks into full scores (1, 1)/(0, 1)/(0, 0).
        Leaves full scores alone.
        """
        if val == True:
            return (1, 1)
        elif val == False:
            return (0, 1)
        elif val is None:
            return (0, 0)

        return val

