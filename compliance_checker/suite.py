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

    def run(self, dataset_location, criteria, tests_to_run, verbose, *args):
        """
        Runs this CheckSuite on the dataset with all the passed Checker instances.

        Returns a dictionary mapping Checkers to their grouped scores.
        """

        ret_val = {}
        check_number = 0
        for a in args:

            ds = self.load_dataset(dataset_location, a.beliefs())
            a.setup(ds)
            checks = self._get_checks(a)

            vals = list(itertools.chain.from_iterable(map(lambda c: self._run_check(c, ds), checks)))
            groups = self.scores(vals)
            #Calls output routine to display results in terminal, including scoring.  Goes to verbose function if called by user.
            score_list, fail_flag, check_number, limit = self.standard_output(criteria, check_number, groups, tests_to_run)
            

            if verbose:
                self.verbose_output_generation(groups, verbose, score_list, limit)
        return fail_flag

    def standard_output(self, criteria, check_number, groups, tests_to_run):
        '''
        Generates the Terminal Output for Standard cases

        Returns the dataset needed for the verbose output, as well as the failure flags.
        '''
        if criteria == 'normal':
            limit = 2
        elif criteria == 'strict':
            limit = 1
        elif criteria == 'lenient':
            limit = 3

        score_list = []
        score_only_list= []

        for v in range(len(groups)):
            score_list.append([groups[v].name, groups[v].weight, groups[v].value, groups[v].children])
            if groups[v].weight >= limit:
                score_only_list.append(groups[v].value)

        points = [x[0] for x in score_only_list]
        out_of = [x[1] for x in score_only_list]

        points = sum(points)
        out_of = sum(out_of)

        score_list.sort(key=lambda x: x[1], reverse=True)

        fail_flag = 0

        if points < out_of:
            fail_flag = limit
            print '\n'
            print "-"*55
            print "   The dataset scored %r out of %r required points" % (points, out_of)
            print "            during the %s check" % tests_to_run[check_number]
            print "      This test has passed under %s critera" % criteria
            print "-"*55
        check_number = check_number +1      
        return [score_list, fail_flag, check_number, limit]


    def verbose_output_generation(self, groups, verbose, score_list, limit):
        '''
        Generates the Terminal Output for Verbose cases
        '''
        
        sub_tests = []
        if verbose == 1:
            print "\n"+"-"*55
            print "The following tests failed:" 
            priority_flag = 3

            for x in range(len(score_list)):
                if score_list[x][1] == 3 and limit <= 3 and priority_flag == 3:
                    print '----High priority tests failed-----'
                    print '%-36s:%8s:%6s' % ('    Name', 'Priority', 'Score')
                    priority_flag -= 1
                elif score_list[x][1] == 2 and limit <= 2 and priority_flag == 2:
                    print '----Medium priority tests failed-----'
                    print '%-36s:%8s:%6s' % ('    Name', 'Priority', 'Score')
                    priority_flag -=1
                elif score_list[x][1] == 1 and limit <= 1 and priority_flag == 1:
                    print '----Low priority tests failed-----'
                    print '%-36s:%8s:%6s' % ('    Name', 'Priority', 'Score')
                    priority_flag -= 1
                if score_list[x][2][0] < score_list[x][2][1] and score_list[x][1] >= limit:
                    print '%-40s:%s:%6s/%1s'  % (score_list[x][0], score_list[x][1], score_list[x][2][0], score_list[x][2][1])
                    print '-'*55

        if verbose >= 2:
            print "-"*55
            print "Summary of all the checks performed:" 
            
            priority_flag = 3

            self.print_routine(groups, 0, verbose, priority_flag)

        pass


    def print_routine(self, list_of_results, indent, verbose, priority_flag):
        """
        print routine performed
        """
        def weight_func(r):
            """
            Function that returns the weight, used for sorting by priority
            """
            return r.weight

        #Sorting method used to properly sort the output by priority.
        grouped_sorted = []
        grouped_sorted = sorted(list_of_results, key=weight_func, reverse=True)        


        #Loop over inoput
        for res in grouped_sorted:
            
            #If statements to print the proper Headings
            if res.weight == 3 and indent == 0 and priority_flag == 3:
                print "\nHigh Priority"
                print "-------------"
                print '%-36s:%8s:%6s' % ('    Name', 'Priority', 'Score')

                priority_flag -= 1
            if res.weight == 2 and indent == 0 and priority_flag == 2:
                print "\nMedium Priority"
                print "---------------"
                print '%-36s:%8s:%6s' % ('    Name', 'Priority', 'Score')

                priority_flag -= 1
            if res.weight ==1 and indent ==0 and priority_flag == 1:
                print "\nLow Priority"
                print "------------"
                print '%-36s:%8s:%6s' % ('    Name', 'Priority', 'Score')
                priority_flag -= 1


            print '%-40s:%s:%6s/%1s' % (indent*'    '+res.name, res.weight, res.value[0], res.value[1])
            if res.children and verbose >1:
                self.print_routine(res.children, indent+1, verbose-1, priority_flag)
            if indent == 0:
                print '-'*55


    def load_dataset(self, ds_str, belief_map):
        """
        Helper method to load a dataset.
        """
        ds = Dataset(ds_str)
        data_object = NetCDFDogma('ds', belief_map, ds)

        return DSPair(ds, data_object)

    def scores(self, raw_scores):
        """
        Transforms raw scores from a single checker into a fully tallied and grouped scoreline.
        """
        grouped = self._group_raw(raw_scores)

        return (grouped)


    def _group_raw(self, raw_scores, cur=None, level=1):
        """
        Internal recursive method to group raw scores into a cascading score summary.

        Only top level items are tallied for scores.
        """

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

