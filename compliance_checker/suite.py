"""
Compliance Checker suite runner
"""

import sys
import inspect
import itertools
from netCDF4 import Dataset
from lxml import etree as ET
from compliance_checker.base import BaseCheck, BaseNCCheck, fix_return_value, Result
from owslib.sos import SensorObservationService
from owslib.swe.sensor.sml import SensorML
from urlparse import urlparse
import requests
import textwrap
import pkg_resources

class CheckSuite(object):

    checkers = {}       # Base dict of checker names to BaseCheck derived types, override this in your CheckSuite implementation

    @classmethod
    def load_all_available_checkers(cls):
        """
        Helper method to retrieve all sub checker classes derived from various
        base classes.
        """
        from pkg_resources import working_set
        for x in working_set.iter_entry_points('compliance_checker.suites'):
            try:
                xl = x.load()

                cls.checkers[xl.name] = xl
            except Exception as e:
                print >>sys.stderr, "Could not load", x, ":", e

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
            return [fix_return_value(v, check_method.im_func.func_name, check_method, check_method.im_self) for v in val]

        return [fix_return_value(val, check_method.im_func.func_name, check_method, check_method.im_self)]

    def _get_valid_checkers(self, ds, checker_names):
        """
        Returns a filtered list of 2-tuples: (name, valid checker) based on the ds object's type and
        the user selected names.
        """
        if len(checker_names) == 0:
            checker_names = self.checkers.keys()

        args = [(name, self.checkers[name]) for name in checker_names if name in self.checkers]
        valid = []

        all_checked = set([a[1] for a in args])  # only class types
        checker_queue = set(args)

        while len(checker_queue):
            name, a = checker_queue.pop()
            if type(ds) in a().supported_ds:
                valid.append((name, a))

            # add all to queue
            for subc in a.__subclasses__():
                if subc not in all_checked:
                    all_checked.add(subc)
                    checker_queue.add((name, subc))

        return valid

    def run(self, ds, *checker_names):
        """
        Runs this CheckSuite on the dataset with all the passed Checker instances.

        Returns a dictionary mapping checker names to a 2-tuple of their grouped scores and errors/exceptions while running checks.
        """

        ret_val      = {}
        fail_flag    = False

        checkers     = self._get_valid_checkers(ds, checker_names)

        if len(checkers) == 0:
            print "No valid checkers found for tests '%s'" % ",".join(checker_names)

        for checker_name, checker_class in checkers:

            checker            = checker_class()   # @TODO: combine with load_datapair/setup
            dsp                = checker.load_datapair(ds)
            checker.setup(dsp)

            checks             = self._get_checks(checker)
            vals               = []
            errs               = {}   # check method name -> (exc, traceback)

            for c in checks:
                try:
                    vals.extend(self._run_check(c, dsp))
                except Exception as e:
                    errs[c.im_func.func_name] = (e, sys.exc_info()[2])

            # score the results we got back
            groups = self.scores(vals)

            ret_val[checker_name] = groups, errs

        return ret_val

    @classmethod
    def passtree(cls, groups, limit):
        for r in groups:
            if r.children:
                x = cls.passtree(r.children, limit)
                if r.weight >= limit and x == False:
                    return False

            if r.weight >= limit and r.value[0] != r.value[1]:
                return False

        return True

    def html_output(self, check_name, groups, file_object, source_name):
        '''
        Renders an HTML file using Jinja2 and saves the output to the file specified.

        @param check_name      The test which was run
        @param groups          List of results from compliance checker
        @param output_filename Path to file to save output
        @param source_name     Source of the dataset, used for title
        '''
        from jinja2 import Environment, PackageLoader
        self.j2 = Environment(loader=PackageLoader('compliance_checker', 'data/templates'))
        template = self.j2.get_template('ccheck.html.j2')

        template_vars = {}

        template_vars['scored_points'] = 0
        template_vars['possible_points'] = 0
        high_priorities   = []
        medium_priorities = []
        low_priorities    = []
        all_priorities    = []

        template_vars['high_count']   = 0
        template_vars['medium_count'] = 0
        template_vars['low_count']    = 0

        def named_function(result):
            for child in result.children:
                template_vars['scored_points'] += child.value[0]
                template_vars['possible_points'] += child.value[1]
                all_priorities.append(child)
                named_function(child)



        # For each result, bin them into the appropriate category, put them all
        # into the all_priorities category and add up the point values
        for res in groups:
            template_vars['scored_points'] += res.value[0]
            template_vars['possible_points'] += res.value[1]
            if res.weight == 3:
                high_priorities.append(res)
                if res.value[0] < res.value[1]:
                    template_vars['high_count'] += 1
            elif res.weight == 2:
                medium_priorities.append(res)
                if res.value[0] < res.value[1]:
                    template_vars['medium_count'] += 1
            else:
                low_priorities.append(res)
                if res.value[0] < res.value[1]:
                    template_vars['low_count'] += 1
            all_priorities.append(res)
            # Some results have children
            # We don't render children inline with the top three tables, but we
            # do total the points and display the messages
            named_function(res)

        template_vars['high_priorities']   = high_priorities
        template_vars['medium_priorities'] = medium_priorities
        template_vars['low_priorities']    = low_priorities
        template_vars['all_priorities']    = all_priorities
        template_vars['testname']          = check_name
        template_vars['source_name']       = source_name

        buf = template.render(**template_vars)
        file_object.write(buf)

    def get_points(self, groups, limit):
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

        return score_list, points, out_of


    def standard_output(self, limit, check_name, groups):
        """
        Generates the Terminal Output for Standard cases

        Returns the dataset needed for the verbose output, as well as the failure flags.
        """
        score_list, points, out_of = self.get_points(groups, limit)

        fail_flag = 0

        fail_flag = limit
        print '\n'
        print "-"*80
        print '{:^80}'.format("The dataset scored %r out of %r points" % (points, out_of))
        print '{:^80}'.format("during the %s check" % check_name)
        print "-"*80

        return  [score_list, points, out_of]

    def non_verbose_output_generation(self, score_list, groups, limit, points, out_of):

        if points < out_of:
            print '{:^80}'.format("Scoring Breakdown:")
            print '\n'
            priority_flag = 3
            for x in range(len(score_list)):
                if score_list[x][1] == 3 and limit <= 3 :
                    if priority_flag ==3:
                        print '{:^80}'.format("High Priority")
                        print "-"*80
                        print '%-36s:%8s:%6s' % ('    Name', 'Priority', 'Score')
                        priority_flag -= 1
                    print '%-40s:%s:%6s/%1s'  % (score_list[x][0][0:39], score_list[x][1], score_list[x][2][0], score_list[x][2][1])
                
                elif score_list[x][1] == 2 and limit <= 2 :
                    if priority_flag ==2:
                        print '\n'
                        print '{:^80}'.format("Medium Priority")
                        print "-"*80
                        print '%-36s:%8s:%6s' % ('    Name', 'Priority', 'Score')
                        priority_flag -= 1
                    print '%-40s:%s:%6s/%1s'  % (score_list[x][0][0:39], score_list[x][1], score_list[x][2][0], score_list[x][2][1])
                    
                elif score_list[x][1] == 1 and limit == 1 :
                    if priority_flag ==1:
                        print '\n'
                        print '{:^80}'.format("Low Priority")
                        print "-"*80
                        print '%-36s:%8s:%6s' % ('    Name', 'Priority', 'Score')
                        priority_flag -= 1
                    print '%-40s:%s:%6s/%1s'  % (score_list[x][0][0:39], score_list[x][1], score_list[x][2][0], score_list[x][2][1])
                    
                elif score_list[x][1] == 1 and limit == 1 and priority_flag == 2:
                    print '{:^80}'.format('No medium priority tests present')
                    print '-'*80
                    priority_flag -= 1
            #Catch All for pretty presentation
            if priority_flag == 2 and limit == 2:
                print '{:^80}'.format('No Medium priority tests present')
                print '-'*80

            if priority_flag == 2 and limit == 1:
                print '{:^80}'.format('No Medium priority tests present')
                print '-'*80
                print ''
                print '{:^80}'.format('No Low priority tests present')
                print '-'*80

            if priority_flag == 1 and limit == 1:
                print '{:^80}'.format('No Low priority tests present')
                print '-'*80

            print "\n"+"\n"+'-'*80
            print '{:^80}'.format('Reasoning for the failed tests given below:')
            print '\n'
            print '%s%37s:%10s:%8s' % ('Name','Priority', '  Score', 'Reasoning')
            print "-"*80
            self.reasoning_routine(groups, 0)

        else: 
            print "All tests passed!"



    def verbose_output_generation(self, groups, limit, points, out_of):
        '''
        Generates the Terminal Output for Verbose cases
        '''
        priority_flag = 3
        print '{:^80}'.format("Verbose Scoring Breakdown:"),
        self.print_routine(groups, 0, priority_flag)
        if points < out_of:
            print "\n"+"\n"+'-'*80
            print '{:^80}'.format('Reasoning for the failed tests given below:')
            print '\n'
            print '%s%37s:%10s:%8s' % ('Name','Priority', '  Score', 'Reasoning')
            print "-"*80
            self.reasoning_routine(groups, 0)

        pass


    def print_routine(self, list_of_results, indent, priority_flag):
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

        #Loop over input
        for res in grouped_sorted:
            #If statements to print the proper Headings
            if res.weight == 3 and indent == 0 and priority_flag == 3:
                print '\n'
                print '{:^80}'.format("High Priority")
                print "-"*80
                print '%-36s:%8s:%6s' % ('    Name', 'Priority', 'Score')

                priority_flag -= 1
            if res.weight == 2 and indent == 0 and priority_flag == 2:
                print '\n'
                print '{:^80}'.format("Medium Priority")
                print "-"*80
                print '%-36s:%8s:%6s' % ('    Name', 'Priority', 'Score')

                priority_flag -= 1
            if res.weight ==1 and indent ==0 and priority_flag == 1:
                print '\n'
                print '{:^80}'.format("Low Priority")
                print "-"*80
                print '%-36s:%8s:%6s' % ('    Name', 'Priority', 'Score')
                priority_flag -= 1


            print '%-40s:%s:%s%6s/%1s' % ((indent*'    '+res.name)[0:39], res.weight, indent*'  ', res.value[0], res.value[1])
            if res.children:
                self.print_routine(res.children, indent+1, priority_flag)


    def reasoning_routine(self, list_of_results, indent, line = True):
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

        wrapper = textwrap.TextWrapper(initial_indent = '', width = 80, subsequent_indent = ' '*54)
        for res in grouped_sorted:
            if (res.value[0] != res.value[1]) and not res.msgs:
                print '%-39s:%1s:%6s/%2s : %s' %(str(indent*'    '+res.name)[0:39], res.weight, str(res.value[0]), str(res.value[1]), ' ')
            
            if (res.value[0] != res.value[1]) and res.msgs:
                print wrapper.fill('%-39s:%1s:%6s/%2s : %s' %(str(indent*'    '+res.name)[0:39], res.weight, str(res.value[0]), str(res.value[1]), str(", ".join(res.msgs))))

            if res.children:
                self.reasoning_routine(res.children, indent+1, False)

    def load_dataset(self, ds_str):
        """
        Helper method to load a dataset or SOS GC/DS url.
        """
        ds = None

        # try to figure out if this is a local NetCDF Dataset, a remote one, or an SOS GC/DS url
        doc = None
        pr = urlparse(ds_str)
        if pr.netloc:       # looks like a remote url
            rhead = requests.head(ds_str)

            # if we get a 400 here, it's likely a Dataset openable OpenDAP url
            if rhead.status_code == 400:
                pass
            elif rhead.status_code == 200 and rhead.headers['content-type'] == 'text/xml':
                # probably interesting, grab it
                r = requests.get(ds_str)
                r.raise_for_status()

                doc = r.text
            else:
                raise StandardError("Could not understand response code %s and content-type %s" % (rhead.status_code, rhead.headers.get('content-type', 'none')))
        else:
            # do a cheap imitation of libmagic
            # http://stackoverflow.com/a/7392391/84732
            textchars = ''.join(map(chr, [7,8,9,10,12,13,27] + range(0x20, 0x100)))
            is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))

            with open(ds_str) as f:
                first_chunk = f.read(1024)
                if is_binary_string(first_chunk):
                    # likely netcdf file
                    pass
                else:
                    f.seek(0)
                    doc = "".join(f.readlines())

        if doc is not None:
            xml_doc = ET.fromstring(str(doc))
            if xml_doc.tag == "{http://www.opengis.net/sos/1.0}Capabilities":
                ds = SensorObservationService(ds_str, xml=str(doc))

            elif xml_doc.tag == "{http://www.opengis.net/sensorML/1.0.1}SensorML":
                ds = SensorML(xml_doc)
            else:
                raise StandardError("Unrecognized XML root element: %s" % xml_doc.tag)
        else:
            # no doc? try the dataset constructor
            ds = Dataset(ds_str)

        return ds

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

