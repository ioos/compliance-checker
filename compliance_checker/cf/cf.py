import re
from functools import wraps
from collections import defaultdict, OrderedDict
import numpy as np

from compliance_checker.base import BaseCheck, BaseNCCheck, check_has, score_group, Result
from compliance_checker.cf.appendix_d import dimless_vertical_coordinates
from compliance_checker.cf.util import NCGraph, StandardNameTable, units_known, units_convertible, units_temporal

from netCDF4 import Dimension, Variable


# copied from paegan
# paegan may depend on these later
_possiblet = ["time", "TIME", "Time",
           "t", "T",
           "ocean_time", "OCEAN_TIME",
           "jd", "JD",
           "dn", "DN",
           "times", "TIMES", "Times",
           "mt", "MT",
           "dt", "DT",
          ]
_possiblez = ["depth", "DEPTH",
           "depths", "DEPTHS",
           "height", "HEIGHT",
           "altitude", "ALTITUDE",
           "alt", "ALT", 
           "Alt", "Altitude",
           "h", "H",
           "s_rho", "S_RHO",
           "s_w", "S_W",
           "z", "Z",
           "siglay", "SIGLAY",
           "siglev", "SIGLEV",
           "sigma", "SIGMA",
           "vertical", "VERTICAL", "lev", "LEV", "level", "LEVEL"
          ]
_possiblex = ["x", "X",
           "lon", "LON",
           "xlon", "XLON",
           "lonx", "lonx",
           "lon_u", "LON_U",
           "lon_v", "LON_V",
           "lonc", "LONC",
           "Lon", "Longitude",
           "longitude", "LONGITUDE",
           "lon_rho", "LON_RHO",
           "lon_psi", "LON_PSI",

          ]
_possibley = ["y", "Y",
           "lat", "LAT",
           "ylat", "YLAT",
           "laty", "laty",
           "lat_u", "LAT_U",
           "lat_v", "LAT_V",
           "latc", "LATC",
           "Lat", "Latitude",
           "latitude", "LATITUDE",
           "lat_rho", "LAT_RHO",
           "lat_psi", "LAT_PSI",

          ]

_possibleaxis = _possiblet + _possiblez + _possiblex + _possibley


_possiblexunits = ['degrees_east',
                    'degree_east',
                    'degrees_E',
                    'degree_E',
                    'degreesE',
                    'degreeE'
                    ]

_possibleyunits = ['degrees_north',
                'degree_north',
                'degrees_N',
                'degree_N',
                'degreesN',
                'degreeN'
                    ]

_possibletunits = ['day', 
                'days', 
                'd', 
                'hour', 
                'hours', 
                'hr', 
                'hrs', 
                'h', 
                'year', 
                'years', 
                'minute', 
                'minutes', 
                'm', 
                'min', 
                'mins', 
                'second', 
                'seconds', 
                's', 
                'sec', 
                'secs'
                ]

_possibleaxisunits =  _possiblexunits + _possibleyunits +_possibletunits

def guess_dim_type(dimension):
    """
    Guesses the type of dimension of a variable X/Y/Z/T

    If can't figure it out, None is returned.
    """

    dimclasses = {'T':_possiblet,
                  'Z':_possiblez,
                  'Y':_possibley,
                  'X':_possiblex}

    for dcname, dcvals in dimclasses.iteritems():
        if dimension in dcvals:
            return dcname

    return None

def guess_coord_type(units, positive=None):
    """
    Guesses coordinate variable type (X/Y/Z/T) from units and positive attrs
    """
    coord_types = {'X':['degrees_east', 'degree_east', 'degrees_E', 'degree_E', 'degreesE', 'degreeE'],
                   'Y':['degrees_north', 'degree_north', 'degrees_N', 'degree_N', 'degreesN', 'degreeN']}

    deprecated = ['level', 'layer', 'sigma_level']

    if not units or not isinstance(units, basestring) or not (units_known(units) or units in deprecated):
        return None

    if isinstance(positive, basestring) and positive.lower() in ['up', 'down']:
        # only Z if units without positive deos not result in something else
        if guess_coord_type(units, None) in [None, 'Z']:
            return 'Z'
        else:
            # they differ, then we can't conclude
            return None

    if units in deprecated or units_convertible(units, 'hPa', reftimeistime=True):
        return 'Z'

    if positive:
        return None

    for ctype, unitsposs in coord_types.iteritems():
        if units in unitsposs:
            return ctype

    if units_convertible(units, 'days since 1970-01-01', reftimeistime=False):
        return 'T'

    return None

def is_variable(name, var):
    dims = var.dimensions
    if (name,) == dims:
        # Coordinate Type
        return False
    # Probably a variable
    return True


def map_axes(dim_vars, reverse_map=False):
    """
    axis name       -> [dimension names]
    dimension name  -> [axis_name], length 0 if reverse_map
    """
    ret_val = defaultdict(list)
    axes = ['X', 'Y', 'Z', 'T']

    for k, v in dim_vars.iteritems():
        axis = getattr(v, 'axis', '')
        if not axis:
            continue

        axis = axis.upper()
        if axis in axes:
            if reverse_map:
                ret_val[k].append(axis)
            else:
                ret_val[axis].append(k)

    return dict(ret_val)

# helper to see if we should do DSG tests
def is_likely_dsg(func):
    @wraps(func)
    def _dec(s, ds):
        if hasattr(ds.dataset, 'featureType'):
            return func(s, ds)

        # @TODO: skips if we have formalized skips
        return None

    return _dec

class CFBaseCheck(BaseCheck):
    @classmethod
    def beliefs(cls): # @TODO
        return {}
    """
    CF Convention Checker (1.6)

    These checks are translated documents: 
        http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/cf-conventions.html
        http://cf-pcmdi.llnl.gov/conformance/requirements-and-recommendations/1.6/
    """

    def __init__(self):
        self._coord_vars     = defaultdict(list)
        self._ancillary_vars = defaultdict(list)
        self._clim_vars      = defaultdict(list)
        self._boundary_vars  = defaultdict(dict)

        self._std_names      = StandardNameTable('cf-standard-name-table.xml')

    ################################################################################
    #
    # Helper Methods - var classifications, etc
    #
    ################################################################################

    def setup(self, ds):
        self._find_coord_vars(ds)
        self._find_ancillary_vars(ds)
        self._find_clim_vars(ds)
        self._find_boundary_vars(ds)

    def _find_coord_vars(self, ds, refresh=False):
        """
        Finds all coordinate variables in a dataset.

        A variable with the same name as a dimension is called a coordinate variable.

        The result is cached by the passed in dataset object inside of this checker. Pass refresh=True
        to redo the cached value.
        """
        if ds in self._coord_vars and not refresh:
            return self._coord_vars[ds]

        for d in ds.dataset.dimensions:
            if d in ds.dataset.variables and ds.dataset.variables[d].dimensions == (d,):
                self._coord_vars[ds].append(ds.dataset.variables[d])

        return self._coord_vars[ds]

    def _find_ancillary_vars(self, ds, refresh=False):
        """
        Finds all ancillary variables in a dataset.

        TODO: fully define

        An ancillary variable generally is a metadata container and referenced from
        other variables via a string reference in an attribute.

        - via ancillary_variables (3.4)
        - "grid mapping var" (5.6)
        - TODO: more?

        The result is cached by the passed in dataset object inside of this checker. Pass refresh=True
        to redo the cached value.
        """
        if ds in self._ancillary_vars and not refresh:
            return self._ancillary_vars[ds]

        for name, var in ds.dataset.variables.iteritems():
            if hasattr(var, 'ancillary_variables'):
                for anc_name in var.ancillary_variables.split(" "):
                    if anc_name in ds.dataset.variables:
                        self._ancillary_vars[ds].append(ds.dataset.variables[anc_name])

            if hasattr(var, 'grid_mapping'):
                gm_name = var.grid_mapping
                if gm_name in ds.dataset.variables:
                    self._ancillary_vars[ds].append(ds.dataset.variables[gm_name])

        return self._ancillary_vars

    def _find_data_vars(self, ds):
        """
        Finds all variables that could be considered Data variables.

        Returns a dictionary mapping name -> variable.

        Excludes variables that are:
            - coordinate variables
            - ancillary variables
            - no dimensions

        Results are NOT CACHED.
        """
        return {k:v for k, v in ds.dataset.variables.iteritems() if v not in self._find_coord_vars(ds) \
                                                                  and v not in self._find_ancillary_vars(ds) \
                                                                  and v.dimensions}

    def _find_clim_vars(self, ds, refresh=False):
        """
        Finds all climatology variables in a dataset.

        Climatology variables (7.4)

        Cached.
        """

        if ds in self._clim_vars and not refresh:
            return self._clim_vars[ds]

        c_time = set()  # set of climatological time axes
        c_vars = set()  # set of climatological variables

        coord_vars = self._find_coord_vars(ds)

        # find all time dimension variables
        time_vars = [v for v in coord_vars if guess_dim_type(v.dimensions[0]) == 'T']

        for k, v in ds.dataset.variables.iteritems():
            is_cvar = False

            if v in coord_vars:
                if hasattr(v, 'climatology') and not hasattr(v, 'bounds'):
                    is_cvar = True

                if k in time_vars and hasattr(v, 'units'):
                    units_split = v.units.split()
                    if len(units_split) == 3 and units_split[1:] == ['since', '0-1-1']:
                        is_cvar = True

                if is_cvar:
                    c_time.add(v)

            else:
                # check cell_methods
                if hasattr(v, 'cell_methods'):
                    try:
                        cell_methods = parse_cell_methods(v.cell_methods)
                    except:
                        pass

        dvars = set(ds.dataset.variables.itervalues()) - set(coord_vars)
        for dim in c_time:
            for v in dvars:
                if dim in v.dimensions:
                    c_vars.add(v)

        self._clim_vars[ds] = c_vars

        return c_vars

    def _find_boundary_vars(self, ds, refresh=False):
        """
        Returns dict of coordinates vars -> associated boundary vars

        Cached.
        """
        if ds in self._boundary_vars and not refresh:
            return self._boundary_vars[ds]

        b_vars = {}

        for k, v in ds.dataset.variables.iteritems():
            bounds = getattr(v, 'bounds', None)
            if bounds is not None and isinstance(bounds, basestring) and bounds in ds.dataset.variables:
                b_vars[v] = ds.dataset.variables[bounds]

        self._boundary_vars[ds] = b_vars
        return b_vars

    ###############################################################################
    #
    # CHAPTER 2: NetCDF Files and Components
    #
    ###############################################################################



    def check_data_types(self, ds):
        """
        2.2 The netCDF data types char, byte, short, int, float or real, and double are all acceptable
        """
        fails = []
        total = len(ds.dataset.variables)

        for k, v in ds.dataset.variables.iteritems():
            if v.datatype not in [np.character,
                                  np.dtype('c'),
                                  np.dtype('b'),
                                  np.dtype('i4'),
                                  np.int32,
                                  np.float32,
                                  np.double,
                                  'int16',
                                  'float32'
                                  ]:

                fails.append(('The variable %s failed because the datatype is %s' %(k, v.datatype)))
        return Result(BaseCheck.HIGH, (total - len(fails), total), msgs=fails)

    def check_naming_conventions(self, ds):
        """
        2.3 Variable, dimension and attribute names should begin with a letter and be composed of letters, digits, and underscores. 
        """
        fails = []
        total = len(ds.dataset.variables)

        rname = re.compile("[A-Za-z][A-Za-z0-9_]*")

        for k, v in ds.dataset.variables.iteritems():
            if not rname.match(k):
                fails.append('Variable %s failed because it does not start with a letter, digit, or an underscore' %k)

        return Result(BaseCheck.HIGH, (total - len(fails), total), 'Variable names', fails)

    def check_names_unique(self, ds):
        """
        2.3 names should not be distinguished purely by case, i.e., if case is disregarded, no two names should be the same.
        """
        fails = []
        total = len(ds.dataset.variables)
        names = defaultdict(int)

        for k in ds.dataset.variables:
            names[k.lower()] += 1

        fails = ['Variables are not case sensitive.  Duplicate variables named: %s' %k for k,v in names.iteritems() if v > 1]
        return Result(BaseCheck.LOW, (total - len(fails), total), msgs=fails)

    def check_dimension_names(self, ds):
        """
        2.4 A variable may have any number of dimensions, including zero, and the dimensions must all have different names.
        """
        fails = []
        total = len(ds.dataset.variables)

        for k, v in ds.dataset.variables.iteritems():
            dims = defaultdict(int)
            for d in v.dimensions:
                dims[d] += 1

            cur_fails = [(k, kk) for kk, vv in dims.iteritems() if vv > 1]
            fails.extend(cur_fails)

        return Result(BaseCheck.HIGH, (total - len(fails), total), msgs=fails)

    def check_dimension_order(self, ds):
        """
        2.4 If any or all of the dimensions of a variable have the interpretations of "date or time" (T), "height or depth" (Z),
        "latitude" (Y), or "longitude" (X) then we recommend, those dimensions to appear in the relative order T, then Z, then Y,
        then X in the CDL definition corresponding to the file. All other dimensions should, whenever possible, be placed to the
        left of the spatiotemporal dimensions.
        """
        fails = []
        total = len(ds.dataset.variables)

        expected = ['T', 'Z', 'Y', 'X']

        for k, v in ds.dataset.variables.iteritems():

            dclass = map(guess_dim_type, v.dimensions)
            # any nones should be before classified ones
            nones    = [i for i, x in enumerate(dclass) if x is None]
            nonnones = [i for i, x in enumerate(dclass) if x is not None]

            if len(nones) and len(nonnones) and max(nones) > min(nonnones):
                fails.append("Variable %s has a non-space-time dimension after space-time-dimensions"%k)

            # classified ones should be in correct order
            nonnones = [expected.index(x) for x in dclass if x is not None]
            nonnones_sorted = sorted(nonnones)

            if nonnones != nonnones_sorted:
                fails.append("The dimensions for %s are not in T Z Y X order"%k)

        # there are two checks here per variable so totals must be doubled
        return Result(BaseCheck.LOW, (total*2 - len(fails), total*2), msgs=fails)

    #def check_dimension_single_value_applicable(self, ds):
        """
        2.4 When a single value of some coordinate applies to all the values in a variable, the recommended means of attaching this
        information to the variable is by use of a dimension of size unity with a one-element coordinate variable. It is also
        acceptable to use a scalar coordinate variable which eliminates the need for an associated size one dimension in the data
        variable.
        """


        #TODO: We need to identify a non-compliant example of this that can be verified, but I believe
        #      that if the file is netCDF then this requirement may be met.  When we do we can reinsert this check
        #pass

    def check_fill_value_outside_valid_range(self, ds):
        """
        2.5.1 The _FillValue should be outside the range specified by valid_range (if used) for a variable.
        """
        fails = []
        checked = 0

        for k, v in ds.dataset.variables.iteritems():
            if hasattr(v, '_FillValue'):
                attrs = v.ncattrs()

                if 'valid_range' in attrs:
                    rmin, rmax = v.valid_range
                elif 'valid_min' in attrs and 'valid_max' in attrs:
                    rmin = v.valid_min
                    rmax = v.valid_max
                else:
                    continue

                checked += 1

                if v._FillValue >= rmin and v._FillValue <= rmax:
                    fails.append((k, "%s is between %s and %s" % (v._FillValue, rmin, rmax)))
        if checked >= 1:
            return Result(BaseCheck.HIGH, (checked - len(fails), checked), msgs=fails)
        else:
            return []

    def check_conventions_are_cf_16(self, ds):
        """
        2.6.1 the NUG defined global attribute Conventions to the string value "CF-1.6"
        """
        if hasattr(ds.dataset, 'Conventions'):
            if getattr(ds.dataset, 'Conventions', '') == 'CF-1.6':
                valid = True
                reasoning = ['Conventions field is "CF-1.6"']
            else:
                valid = False
                reasoning = ['Conventions field is not "CF-1.6"']
        else:
            valid = False
            reasoning = ['Conventions field is not present']
        return Result(BaseCheck.HIGH, valid, 'conventions', msgs=reasoning)

    @score_group('convention_attrs')
    def check_convention_globals(self, ds):
        """
        2.6.2 title/history global attributes, must be strings. Do not need to exist.
        """
        attrs = ['title', 'history']
        ret = []

        for a in attrs:
            if hasattr(ds.dataset, a):
                ret.append(Result(BaseCheck.HIGH, isinstance(getattr(ds.dataset, a), basestring), ('global', a)))

        return ret

    @score_group('convention_attrs')
    def check_convention_possibly_var_attrs(self, ds):
        """
        2.6.2 institution, source, references, and comment, either global or assigned to individual variables.
        When an attribute appears both globally and as a variable attribute, the variable's version has precedence.
        Must be strings.
        """
        attrs = ['institution', 'source', 'references', 'comment']
        ret = []

        # check attrs on global ds

        # can't predetermine total - we only report attrs we find
        for k, v in ds.dataset.variables.iteritems():
            vattrs = v.ncattrs()
            for a in attrs:
                if a in vattrs:
                    ret.append(Result(BaseCheck.HIGH, isinstance(getattr(v, a), basestring), (k, a)))

        return ret

    ###############################################################################
    #
    # CHAPTER 3: Description of the Data
    #
    ###############################################################################

    def check_units(self, ds):
        """
        3.1 The units attribute is required for all variables that represent dimensional quantities
        (except for boundary variables defined in Section 7.1, "Cell Boundaries" and climatology variables
        defined in Section 7.4, "Climatological Statistics").

        Units are not required for dimensionless quantities. A variable with no units attribute is assumed
        to be dimensionless. However, a units attribute specifying a dimensionless unit may optionally be
        included.

        - units required
        - type must be recognized by udunits
        - if std name specified, must be consistent with standard name table, must also be consistent with a
          specified cell_methods attribute if present
        """
        ret_val = []

        deprecated = ['level', 'layer', 'sigma_level']

        for k, v in ds.dataset.variables.iteritems():

            # skip climatological vars, boundary vars
            if v in self._find_clim_vars(ds) or \
               v in self._find_boundary_vars(ds).itervalues() or \
               v.shape == ():
               continue

            # skip string type vars
            if v.dtype.char == 'S':
                continue

            units = getattr(v, 'units', None)

            # 1) "units" attribute must be present
            presence = Result(BaseCheck.HIGH, units is not None, ('units', k, 'present'))
            if not presence.value:
                presence.msgs = ['units attribute required']
                ret_val.append(presence)
                continue

            # 2) units attribute must be a string
            astring = Result(BaseCheck.HIGH, isinstance(units, basestring), ('units', k, 'string'))
            if not astring.value:
                astring.msgs = ["units not a string (%s)" % type(units)]
                ret_val.append(astring)
                continue

            # now, units are present and string
            # 3) units are not deprecated
            resdeprecated = Result(BaseCheck.LOW, not units in deprecated, ('units', k, 'deprecated'))
            if not resdeprecated.value:
                resdeprecated.msgs = ['units (%s) is deprecated' % units]
                ret_val.append(resdeprecated)
                continue

            # 4) units are known

            knownu = Result(BaseCheck.HIGH, units_known(units), ('units', k, 'known'))
            if not knownu.value:
                knownu.msgs = ['unknown units type (%s)' % units]
                ret_val.append(knownu)
                #continue
            # units look ok so far, check against standard name / cell methods
            std_name = getattr(v, 'standard_name', None)
            std_name_modifier = None

            if isinstance(std_name, basestring):
                if ' ' in std_name:
                    std_name, std_name_modifier = std_name.split(' ', 1)


            # if no standard name or cell_methods, nothing left to do
            if std_name is None and not hasattr(v, 'cell_methods'):
                #ret_val.append(Result(BaseCheck.HIGH, True, ('units', k, 'ok')))
                continue

            # 5) if a known std_name, use the units provided
            if std_name is not None and std_name in self._std_names:

                std_units = self._std_names[std_name].canonical_units

                #@TODO modifiers changes units
                msgs = []
                valid = True
                if units is not None:
                    if units.split(" ")[0] in ['day', 'days', 'd', 'hour', 'hours', 'hr', 'hrs', 'h', 'year', 'years', 'minute', 'minutes', 'm', 'min', 'mins', 'second', 'seconds', 's', 'sec', 'secs']:
                        if len(units.split(" "))>1:
                            if units.split(" ")[1] == 'since':
                                std_units = units
                        else:
                            std_units = units
    
                    if std_units == 'm' and units in ['meter', 'meters']:
                        std_units = units
    
                    if units != std_units and units not in ['degrees_north', 'degree_N', 'degreeN', 'degreesN', 'degrees_east', 'degree_E', 'degreeE', 'degreesE'] :
                        msgs = ['units are %s, standard_name units should be %s' % (units, std_units)]
                        valid = False
                else:
                    valid = False
                    msgs = ['The unit for variable %s in of type None.'%name]

                ret_val.append(Result(BaseCheck.HIGH, valid, ('units', k, 'standard_name'), msgs))

            # 6) cell methods @TODO -> Isnt this in the check_cell_methods section?
            #if hasattr(v, 'cell_methods'):
            #    cell_methods = v.cell_methods
#
            #    # placemarker for future check
            #    ret_val.append(Result(BaseCheck.HIGH, False, ('units', k, 'cell_methods'), ['TODO: implement cell_methods check']))

        return ret_val

    def check_standard_name(self, ds):
        """
        3.3 A standard name is associated with a variable via the attribute standard_name which takes a
        string value comprised of a standard name optionally followed by one or more blanks and a
        standard name modifier
        """
        ret_val = []

        for k, v in ds.dataset.variables.iteritems():
            std_name = getattr(v, 'standard_name', None)

            std_name_modifier = None

            # no standard name? is ok by the letter of the law
            if std_name is None:
                continue

            if isinstance(std_name, basestring):
                if ' ' in std_name:
                    std_name, std_name_modifier = std_name.split(' ', 1)

            # 1) standard name is a string and in standard name table
            msgs = []
            is_str = isinstance(std_name, basestring)
            in_table = std_name in self._std_names

            if not is_str:
                msgs.append("The standard name '%s' is not of type string.  It is type %s" % (std_name, type(std_name)))
            if not in_table:
                msgs.append("The standard name '%s' is not in standard name table" % std_name)

            ret_val.append(Result(BaseCheck.HIGH, is_str and in_table, ('std_name', k, 'legal'), msgs))

            # 2) optional - if modifiers, should be in table
            if std_name_modifier:
                allowed = ['detection_minimum',
                           'number_of_observations',
                           'standard_error',
                           'status_flag']

                msgs = []
                if not std_name_modifier in allowed:
                    msgs.append("modifier (%s) not allowed" % std_name_modifier)

                ret_val.append(Result(BaseCheck.HIGH, std_name_modifier in allowed, ('std_name', k, 'modifier'), msgs))

        return ret_val

    def check_ancillary_data(self, ds):
        """
        3.4 It is a string attribute whose value is a blank separated list of variable names.
        The nature of the relationship between variables associated via ancillary_variables must
        be determined by other attributes. The variables listed by the ancillary_variables attribute
        will often have the standard name of the variable which points to them including a modifier
        (Appendix C, Standard Name Modifiers) to indicate the relationship.
        """
        ret_val = []

        for k, v in ds.dataset.variables.iteritems():
            anc = getattr(v, 'ancillary_variables', None)
            if anc is None:
                continue

            # should be a string, splittable, and each should exist
            anc_result = Result(BaseCheck.HIGH, name=('ancillary', k))
            msgs = []

            if not isinstance(anc, basestring):
                anc_result.value = False
                anc_result.msgs = ["ancillary_variables is not a string"]
                ret_val.append(anc_result)
                continue

            ancs = anc.split()
            existing = 0

            for a in ancs:
                if a in ds.dataset.variables:
                    existing += 1
                else:
                    msgs.append("ancillary var %s does not exist" % a)

            anc_result.value = (existing, len(ancs))
            anc_result.msgs = msgs

            ret_val.append(anc_result)

        return ret_val

    def check_flags(self, ds):
        """
        3.5 The attributes flag_values, flag_masks and flag_meanings are intended to make variables
        that contain flag values self describing. Status codes and Boolean (binary) condition flags may be
        expressed with different combinations of flag_values and flag_masks attribute definitions.

        The flag_values and flag_meanings attributes describe a status flag consisting of mutually exclusive coded values.

        The flag_meanings attribute is a string whose value is a blank separated list of descriptive words
        or phrases, one for each flag value. Each word or phrase should consist of characters from
        the alphanumeric set and the following five: '_', '-', '.', '+', '@'. 

        The flag_masks and flag_meanings attributes describe a number of independent Boolean conditions
        using bit field notation by setting unique bits in each flag_masks value. 

        The flag_masks, flag_values and flag_meanings attributes, used together, describe a blend of
        independent Boolean conditions and enumerated status codes. A flagged condition is identified
        by a bitwise AND of the variable value and each flag_masks value; a result that matches the
        flag_values value indicates a true condition. 
        """
        ret_val = []

        for k, v in ds.dataset.variables.iteritems():

            flag_values   = getattr(v, "flag_values", None)
            flag_masks    = getattr(v, "flag_masks", None)
            flag_meanings = getattr(v, "flag_meanings", None)

            if not (flag_values is not None or flag_masks is not None):
                continue

            # 1) flags_values attribute must have same type as variable to which it is attached
            if flag_values is not None:
                fvr = Result(BaseCheck.HIGH, flag_values.dtype == v.dtype, name=('flags', k, 'flag_values_type'))
                if not fvr.value:
                    fvr.msgs = ['flag_values attr does not have same type as var (fv: %s, v: %s)' % (flag_values.dtype, v.dtype)]

                ret_val.append(fvr)

                # 2) if flag_values, must have flag_meanings
                fmr = Result(BaseCheck.HIGH, flag_meanings is not None, name=('flags', k, 'flag_meanings_present'))
                if not fmr.value:
                    fmr.msgs = ['flag_meanings must be present']

                ret_val.append(fmr)

                # 8) flag_values attribute values must be mutually exclusive
                fvset = set(flag_values)
                fvsr = Result(BaseCheck.HIGH, len(fvset) == len(flag_values), ('flags', k, 'flag_values_mutually_exclusive'))
                if not fvsr.value:
                    fvsr.msgs = ['repeated items in flag_values']

                ret_val.append(fvsr)

            # 3) type of flag_meanings is a string, blank separated list of words
            if flag_meanings is not None:
                fmt = Result(BaseCheck.HIGH, isinstance(flag_meanings, basestring), name=('flags', k, 'flag_meanings_type'))
                if not fmt.value:
                    fmt.msgs = ['flag_meanings must be a string']

                ret_val.append(fmt)

                # split and check each word
                rflags = re.compile("^[0-9A-Za-z_\-.+@]+$")
                meanings = flag_meanings.split()
                msgs = []
                ok_count = 0

                for fm in meanings:
                    if rflags.match(fm) is not None:
                        ok_count += 1
                    else:
                        msgs.append("flag_meaning %s of var %s is incorrectly named" % (fm, k))

                ret_val.append(Result(BaseCheck.HIGH, (ok_count, len(meanings)), name=('flags', k, 'flag_meanings_names'), msgs=msgs))

                # now that we've split meanings up, check length vs values/masks

                # 4) number of flag_values must equal number of meanings
                if flag_values is not None:
                    fvfmr = Result(BaseCheck.HIGH, len(flag_values) == len(meanings), ('flags', k, 'flag_values_equal_meanings'))
                    if not fvfmr.value:
                        fvfmr.msgs = ['flag_values length (%d) not equal to flag_meanings length (%d)' % (len(flag_values), len(meanings))]

                    ret_val.append(fvfmr)

                # 5) number of flag_masks must equal number of meanings
                if flag_masks is not None:
                    fmfmr = Result(BaseCheck.HIGH, len(flag_masks) == len(meanings), ('flags', k, 'flag_masks_equal_meanings'))
                    if not fmfmr.value:
                        fmfmr.msgs = ['flag_masks length (%d) not equal to flag_meanings length (%d)' % (len(flag_masks), len(meanings))]

                    ret_val.append(fmfmr)

            # 6) flag_masks must have same type as var and those vars must be compatible with bit field expressions
            if flag_masks is not None:
                msgs = []
                ok_count = 0

                same_type = flag_masks.dtype == v.dtype
                type_ok = v.dtype in [np.character,
                                      np.dtype('b'),
                                      np.dtype('i4'),
                                      np.int32]

                if same_type:
                    ok_count += 1
                else:
                    msgs.append("flag_masks is not same type as v (fm: %s, v: %s)" % (flag_masks.dtype, v.dtype))

                if type_ok:
                    ok_count += 1
                else:
                    msgs.append("variable not of appropriate type to have flag_masks (%s)" % (v.dtype))

                ret_val.append(Result(BaseCheck.HIGH, (ok_count, 2), ('flags', k, 'flag_masks_type'), msgs=msgs))

                # 7) the flag_masks attribute values must be non-zero
                zeros = [x for x in flag_masks if x == 0]
                msgs = []
                if len(zeros):
                    msgs = ['flag_masks attribute values contains a zero']

                ret_val.append(Result(BaseCheck.HIGH, len(zeros) != 0, ('flags', k, 'flag_masks_zeros'), msgs=msgs))

            # 9) when both defined, boolean AND of each entry in flag_values with corresponding entry in flag_masks
            #    should equal the flags_value entry
            if flag_values is not None and flag_masks is not None:
                allv = map(lambda a, b: a & b == a, zip(flag_values, flag_masks))

                allvr = Result(BaseCheck.MEDIUM, all(allv), ('flags', k, 'flag_masks_with_values'))
                if not allvr.value:
                    allvr.msgs = ["flag masks and flag values combined don't equal flag value"]

                ret_val.append(allvr)

        return ret_val

    ###############################################################################
    #
    # CHAPTER 4: Coordinate Types
    #
    ###############################################################################

    def check_coordinate_axis_attr(self, ds):
        """
        4 The attribute axis may be attached to a coordinate variable and given one of the values X, Y, Z or T
        which stand for a longitude, latitude, vertical, or time axis respectively. Alternatively the standard_name
        attribute may be used for direct identification.
        """
        ret_val = []
        dim_to_axis = map_axes({k:v for k,v in ds.dataset.variables.iteritems() if v in self._find_coord_vars(ds)}, reverse_map=True)
        data_vars = {k:v for k,v in ds.dataset.variables.iteritems() if v not in self._find_coord_vars(ds)}

        for k, v in ds.dataset.variables.iteritems():
            axis = getattr(v, 'axis', None)

            if axis is None:
                continue

            # 1) axis must be X, Y, Z, or T
            axis_valid = axis in ['X', 'Y', 'Z', 'T']

            avr = Result(BaseCheck.HIGH, axis_valid, ('axis', k, 'valid_value'))
            if not axis_valid:
                avr.msgs = ['axis value (%s) is not valid' % axis]

            ret_val.append(avr)

            # 2) only coordinate vars are allowed to have axis set
            acvr = Result(BaseCheck.HIGH, v in self._find_coord_vars(ds), ('axis', k, 'is_coordinate_var'))
            if not acvr.value:
                acvr.msgs = ['%s is not allowed to have an axis attr as it is not a coordinate var' % k]

            ret_val.append(acvr)

            # 3) must be consistent with coordinate type deduced from units and positive
            axis_type = guess_coord_type(getattr(v, 'units', None), getattr(v, 'positive', None))
            if axis_type is not None:
                atr = Result(BaseCheck.HIGH, axis_type == axis, ('axis', k, 'consistent_with_coord_type'))
                if not atr.value:
                    atr.msgs = ['%s guessed type (%s) is not consistent with coord type (%s)' % (k, axis_type, axis)]

                ret_val.append(atr)

            # 4) a data variable must not have more than one coordinate variable with a particular value of the axis attribute
            if k in data_vars:
                dep_axes = [(dim_to_axis[d], d) for d in v.dimensions if d in dim_to_axis]
                dups = defaultdict(int)
                for d in dep_axes:
                    dups[d[0][0]] += 1

                dups = {kk:vv for kk,vv in dups.iteritems() if vv > 1}

                coores = Result(BaseCheck.HIGH, len(dups) == 0, ('axis', k, 'does_not_depend_on_mult_coord_vars'))
                if not coores.value:
                    coores.msgs = []
                    for kk, vv in dups.iteritems():
                        same_axis = [item[1] for item in dep_axes if item[0] == kk]
                        coores.msgs.append('%s depends on multiple coord vars with axis attribute (%s): %s' % (k, kk, ','.join(same_axis)))

                ret_val.append(coores)

        return ret_val

    def check_coordinate_vars_for_all_coordinate_types(self, ds):
        """
        4 We strongly recommend that coordinate variables be used for all coordinate types whenever they are applicable.
        """
        ret_val = []
        # 1. Verify that for any known or common coordinate name as a dmension
        #    there is a coordinate variable for that dimension.
        known_coordinate_names = ('longitude', 'lon'   , 'x',
                                  'latitude' , 'lat'   , 'y',
                                  'vertical' , 'height', 'z',
                                  'time'               , 't')
        for k,v in ds.dataset.dimensions.iteritems():
            if k.lower() in known_coordinate_names:
                valid = k in ds.dataset.variables
                result = Result(BaseCheck.MEDIUM, valid, ('coordinate_type', k, 'var_for_coordinate_type'))
                if not valid:
                    result.msgs = ['No coordinate variable for coordinate type %s' % k]

                ret_val.append(result)

        #@TODO: Additional verifiable requirements

        return ret_val


    def _coord_has_units(self, name,coordinate, var, recommended, acceptable):
        ret_val = []
        has_units = hasattr(var, 'units')
        result = Result(BaseCheck.HIGH, has_units, (coordinate, name, 'has_units'))
        ret_val.append(result)


        # 0 - does not have units
        # 1 - incorrect units
        # 2 - also acceptable units
        # 3 - recommend units
        if not has_units:
            result = Result(BaseCheck.MEDIUM, (0, 3), (coordinate, name, 'correct_units'), ['%s does not have units'%name])
            ret_val.append(result)
        elif has_units and var.units == recommended:
            result = Result(BaseCheck.MEDIUM, (3, 3), (coordinate, name, 'correct_units'))
            ret_val.append(result)
        elif has_units and var.units in acceptable:
            result = Result(BaseCheck.MEDIUM, (2, 3), (coordinate, name, 'correct_units'), ['%s units are acceptable, but not recommended'%name])
            ret_val.append(result)
        else:
            result = Result(BaseCheck.MEDIUM, (1, 3), (coordinate, name, 'correct_units'), ['%s does not have units'%name])
            ret_val.append(result)
        return ret_val

    def check_latitude(self, ds):
        """
        4.1 Variables representing latitude must always explicitly include the units attribute; there is no default value.
        The recommended unit of latitude is degrees_north. Also acceptable are degree_north, degree_N, degrees_N, degreeN, and degreesN.

        Optionally, the latitude type may be indicated additionally by providing the standard_name attribute with the
        value latitude, and/or the axis attribute with the value Y.
        """
        ret_val = []

        recommended = 'degrees_north'
        acceptable = ['degree_north', 'degree_N', 'degrees_N', 'degreeN', 'degreesN']
    
        for k,v in ds.dataset.variables.iteritems():
            if k == 'latitude' or getattr(v, 'standard_name', None) == 'latitude':
                results = self._coord_has_units(k, 'latitude', v, recommended, acceptable)
                ret_val.extend(results)


        return ret_val

    def check_longitude(self, ds):
        """
        4.2 Variables representing longitude must always explicitly include the units attribute; there is no default value.
        The recommended unit of longitude is degrees_east. Also acceptable are degree_east, degree_E, degrees_E, degreeE, and degreesE.

        Optionally, the longitude type may be indicated additionally by providing the standard_name attribute with the
        value longitude, and/or the axis attribute with the value X.
        """
        ret_val = []

        recommended = 'degrees_east'
        acceptable = ['degree_east', 'degree_E', 'degrees_E', 'degreeE', 'degreesE']
    
        for k,v in ds.dataset.variables.iteritems():
            if k == 'longitude' or getattr(v, 'standard_name', None) == 'longitude':
                results = self._coord_has_units(k, 'longitude', v, recommended, acceptable)
                ret_val.extend(results)


        return ret_val

    def _is_vertical_coordinate(self, var_name, var):
        '''
        Determines if a variable is a vertical coordinate variable
        
        4.3
        A vertical coordinate will be identifiable by: units of pressure; or the presence of the positive attribute with a
        value of up or down (case insensitive).  Optionally, the vertical type may be indicated additionally by providing
        the standard_name attribute with an appropriate value, and/or the axis attribute with the value Z.
        '''
        # Known name
        satisfied = var_name.lower() in _possiblez 
        satisfied |= getattr(var, 'standard_name', '') in _possiblez
        # Is the axis set to Z?
        satisfied |= getattr(var, 'axis', '').lower() == 'z'
        is_pressure = units_convertible(getattr(var, 'units', '1'), 'dbar')
        # Pressure defined or positive defined
        satisfied |= is_pressure
        if not is_pressure:
            satisfied |= getattr(var,'positive', '').lower() in ('up', 'down')
        return satisfied




    def check_vertical_coordinate(self, ds):
        """
        4.3 Variables representing dimensional height or depth axes must always
        explicitly include the units attribute; there is no default value.

        The attribute positive is required if the vertical axis units are not a
        valid unit of pressure. The positive attribute may have the value up or
        down (case insensitive). This attribute may be applied to either
        coordinate variables or auxillary coordinate variables that contain
        vertical coordinate data.
        
        """
        
        ret_val = [] 
        for k,v in ds.dataset.variables.iteritems(): 
            if self._is_vertical_coordinate(k,v):
                # Vertical variables MUST have units
                has_units = hasattr(v, 'units') 
                result = Result(BaseCheck.HIGH, \
                                has_units,      \
                                ('vertical', k, 'has_units'))
                ret_val.append(result)

                # If it's not pressure then it must have positive defined
                if not has_units:
                    result = Result(BaseCheck.HIGH, \
                                    False,          \
                                    ('vertical', k, 'correct_units'),['%s does not have units'%k])
                    ret_val.append(result)
                    continue

                # Do we have pressure?
                is_pressure = units_convertible('dbar', v.units)
                if is_pressure: 
                    result = Result(BaseCheck.HIGH, \
                                    True,           \
                                    ('vertical', k, 'correct_units'))
                # What about positive?
                elif getattr(v,'positive', '').lower() in ('up', 'down'):
                    result = Result(BaseCheck.HIGH, \
                                    True,           \
                                    ('vertical', k, 'correct_units'))
                # Not-compliant
                else:
                    result = Result(BaseCheck.HIGH,                   \
                                    False,                            \
                                    ('vertical', k, 'correct_units'), \
                                    ['vertical variable needs to define positive attribute'])
                ret_val.append(result)
        return ret_val

        

    def check_dimensional_vertical_coordinate(self, ds):
        """
        4.3.1 The units attribute for dimensional coordinates will be a string
        formatted as per the udunits.dat file.

        The acceptable units for vertical (depth or height) coordinate variables
        are:
        - units of pressure as listed in the file udunits.dat. For vertical axes
          the most commonly used of these include include bar, millibar,
          decibar, atmosphere (atm), pascal (Pa), and hPa.
        - units of length as listed in the file udunits.dat. For vertical axes
          the most commonly used of these include meter (metre, m), and
          kilometer (km).
        - other units listed in the file udunits.dat that may under certain
          circumstances reference vertical position such as units of density or
          temperature.

        Plural forms are also acceptable.
        """
        ret_val = []
        for k,v in ds.dataset.variables.iteritems():
            # If this is not a vertical coordinate
            if not self._is_vertical_coordinate(k,v):
                continue

            # If this is not height or depth
            vertical_coordinates = ('height', 'depth')
            if k not in vertical_coordinates and \
                    getattr(v, 'standard_name', '') not in vertical_coordinates:
                continue

            # Satisfies 4.3.1
            # Pressure or length is okay
            is_pressure = units_convertible(getattr(v, 'units', '1'), 'dbar')
            is_length   = units_convertible(getattr(v, 'units', '1'), 'm')
            is_temp     = units_convertible(getattr(v, 'units', '1'), 'degrees_C')
            is_density  = units_convertible(getattr(v, 'units', '1'), 'kg m-3')

            if is_pressure or is_length:
                result = Result(BaseCheck.HIGH, True,                     \
                            ('dimensional_vertical', k, 'correct_units'), \
                            ['dimensional vertical coordinate is pressure or length'])

            # Temperature or Density are okay as well
            elif is_temp or is_density:
                result = Result(BaseCheck.HIGH, True,                     \
                            ('dimensional_vertical', k, 'correct_units'), \
                            ['dimensional vertical coordinate is temp or density'])
            else:
                result = Result(BaseCheck.HIGH, False,                    \
                            ('dimensional_vertical', k, 'correct_units'), \
                            ['incorrect vertical units'])
            ret_val.append(result)

        return ret_val
                    

    def check_dimensionless_vertical_coordinate(self, ds):
        """
        4.3.2 The units attribute is not required for dimensionless coordinates.

        The standard_name attribute associates a coordinate with its definition
        from Appendix D, Dimensionless Vertical Coordinates. The definition
        provides a mapping between the dimensionless coordinate values and
        dimensional values that can positively and uniquely indicate the
        location of the data.

        A new attribute, formula_terms, is used to associate terms in the
        definitions with variables in a netCDF file.  To maintain backwards
        compatibility with COARDS the use of these attributes is not required,
        but is strongly recommended.  
        """
        ret_val = []

        dimless = dict(dimless_vertical_coordinates)
        for k,v in ds.dataset.variables.iteritems():
            std_name = getattr(v, 'standard_name', '')
            if std_name not in dimless:
                continue
            # Determine if the regex matches for formula_terms
            valid_formula = re.match(dimless[std_name],  \
                                     getattr(v, 'formula_terms', ''))

            if valid_formula is not None:
                result = Result(BaseCheck.MEDIUM, \
                                True,             \
                                ('dimensionless_vertical', k, 'formula_terms'))
            else:
                result = Result(BaseCheck.MEDIUM,                              \
                                False,                                         \
                                ('dimensionless_vertical', k, 'formula_terms'),\
                                ['formula_terms missing'])
            ret_val.append(result)

            # Determine that each of the terms actually exists
            # If formula_terms wasn't defined then this fails
            if not valid_formula:
                result = Result(BaseCheck.MEDIUM,                             \
                                False,                                        \
                                ('dimensionless_vertical', k, 'terms_exist'), \
                                ['formula_terms not defined'])
                ret_val.append(result)
                continue

            # Check the terms
            missing_terms = []
            groups = valid_formula.groups()
            for i in xrange(1, len(groups), 2):
                varname = groups[i]
                if varname not in ds.dataset.variables:
                    missing_terms.append(varname)
            # Report the missing terms
            result = Result(BaseCheck.MEDIUM,                             \
                            not missing_terms,                            \
                            ('dimensionless_vertical', k, 'terms_exist'), \
                            ['%s missing' % i for i in missing_terms])

            ret_val.append(result)
                

        return ret_val

    def _is_time_variable(self, varname, var):
        '''
        Identifies if a variable is represents time
        '''
        satisfied = varname.lower() == 'time'
        satisfied |= getattr(var, 'standard_name', '') == 'time'
        satisfied |= getattr(var, 'axis', '') == 'T'
        satisfied |= units_convertible('seconds since 1900-01-01', getattr(var, 'units', ''))
        return satisfied

    def check_time_coordinate(self, ds):
        """
        4.4 Variables representing time must always explicitly include the units
        attribute; there is no default value.

        The units attribute takes a string value formatted as per the
        recommendations in the Udunits package.

        The acceptable units for time are listed in the udunits.dat file. The
        most commonly used of these strings (and their abbreviations) includes
        day (d), hour (hr, h), minute (min) and second (sec, s). Plural forms
        are also acceptable. The reference time string (appearing after the
        identifier since) may include date alone; date and time; or date, time,
        and time zone. The reference time is required. A reference time in year
        0 has a special meaning (see Section 7.4, "Climatological Statistics").

        Recommend that the unit year be used with caution. It is not a calendar
        year.  For similar reasons the unit month should also be used with
        caution.

        A time coordinate is identifiable from its units string alone.
        Optionally, the time coordinate may be indicated additionally by
        providing the standard_name attribute with an appropriate value, and/or
        the axis attribute with the value T.  
        """ 

        ret_val = []
        for k,v in ds.dataset.variables.iteritems():
            if not self._is_time_variable(k,v):
                continue
            # Has units
            has_units = hasattr(v, 'units')
            if not has_units:
                result = Result(BaseCheck.HIGH, \
                                False,          \
                                ('time', k, 'has_units'),['%s does not have the correct units'%k])
                ret_val.append(result)
                result = Result(BaseCheck.HIGH, \
                                False,          \
                                ('time', k, 'correct_units'),['%s does not have the correct units'%k])
                ret_val.append(result)
                continue
            # Correct and identifiable units
            result = Result(BaseCheck.HIGH, \
                            True,           \
                            ('time', k, 'has_units'))
            ret_val.append(result)
            correct_units = units_temporal(v.units)
            result = Result(BaseCheck.HIGH, \
                            correct_units,  \
                            ('time', k, 'correct_units'))
            ret_val.append(result)

        return ret_val


    def check_calendar(self, ds):
        """
        4.4.1 In order to calculate a new date and time given a base date, base
        time and a time increment one must know what calendar to use.

        The values currently defined for calendar are:
        - gregorian or standard
        - proleptic_gregorian
        - noleap or 365_day
        - all_leap or 366_day
        - 360_day
        - julian
        - none

        The calendar attribute may be set to none in climate experiments that
        simulate a fixed time of year.
        The time of year is indicated by the date in the reference time of the
        units attribute.

        If none of the calendars defined above applies, a non-standard calendar
        can be defined. The lengths of each month are explicitly defined with
        the month_lengths attribute of the time axis.

        If leap years are included, then two other attributes of the time axis
        should also be defined:

        leap_year, leap_month

        The calendar attribute is not required when a non-standard calendar is
        being used. It is sufficient to define the calendar using the
        month_lengths attribute, along with leap_year, and leap_month as
        appropriate. However, the calendar attribute is allowed to take
        non-standard values and in that case defining the non-standard calendar
        using the appropriate attributes is required.
        """
        valid_calendars = [
            'gregorian',
            'standard',
            'proleptic_gregorian',
            'noleap',
            '365_day',
            'all_leap',
            '366_day',
            '360_day',
            'julian',
            'none'
        ]

        ret_val = []

        for k,v in ds.dataset.variables.iteritems():
            if not self._is_time_variable(k,v):
                continue
            has_calendar = hasattr(v, 'calendar')
            result = Result(BaseCheck.LOW,  \
                            has_calendar,   \
                            ('time', k, 'has_calendar'))
            ret_val.append(result)
            valid_calendar = has_calendar and v.calendar in valid_calendars
            result = Result(BaseCheck.LOW,  \
                            valid_calendar, \
                            ('time', k, 'valid_calendar'))
            ret_val.append(result)

        return ret_val
        


    ###############################################################################
    #
    # CHAPTER 5: Coordinate Systems
    #
    ###############################################################################

    def _get_coord_vars(self, ds):
        coord_vars = []
        for name,var in ds.dataset.variables.iteritems():
            if (name,) == var.dimensions:
                coord_vars.append(name)
        return coord_vars

    #def check_coordinate_systems(self, ds):
        """
        5 All of a variable's spatiotemporal dimensions that are not latitude,
        longitude, vertical, or time dimensions are required to be associated
        with the relevant latitude, longitude, vertical, or time coordinates via
        the new coordinates attribute of the variable. The value of the
        coordinates attribute is a blank separated list of the names of
        auxiliary coordinate variables.

        The dimensions of an auxiliary coordinate variable must be a subset of
        the dimensions of the variable with which the coordinate is associated,
        with two exceptions:
        - String-valued coordinates (Section 6.1, "Labels") have a dimension for
          maximum string length
        - In the ragged array representations of data (Chapter 9, Discrete
          Sampling Geometries), special methods are needed to connect the data
          and coordinates

        Recommend that the name of a multidimensional coordinate variable should
        not match the name of any of its dimensions because that precludes
        supplying a coordinate variable for the dimension.

        Auxiliary coordinate variables may not be used as the only way to
        identify latitude and longitude coordinates that could be identified
        using coordinate variables.

        An application that is trying to find the latitude coordinate of a
        variable should always look first to see if any of the variable's
        dimensions correspond to a latitude coordinate variable. If the latitude
        coordinate is not found this way, then the auxiliary coordinate
        variables listed by the coordinates attribute should be checked. Note
        that it is permissible, but optional, to list coordinate variables as
        well as auxiliary coordinate variables in the coordinates attribute.

        It is not permissible for a data variable to have both a coordinate
        variable and an auxiliary coordinate variable, or more than one of
        either type of variable, having an axis attribute with any given value
        e.g. there must be no more than one axis attribute for X for any data
        variable.

        """

        #pass

    def check_independent_axis_dimensions(self, ds):
        """
        5.1 When each of a variable's spatiotemporal dimensions is a latitude,
        longitude, vertical, or time dimension, then each axis is identified by
        a coordinate variable.

        """
        ret_val = []

        space_time_coord_var = []
        #Check to find all space-time coordinate variables (Lat/Lon/Time/Height)
        for each in  self._find_coord_vars(ds):
            if str(each._name) in _possibleaxis \
               or (hasattr(each, 'units') and (each.units in _possibleaxisunits or each.units.split(" ")[0]  in _possibleaxisunits)) \
               or hasattr(each,'positive'):
                space_time_coord_var.append(each._name)

        #Find all all space-time variables that are not coordinate variables
        space_time_non_coord_var=[]
        space_time_non_coord_var_dim = []
        for name,var in ds.dataset.variables.iteritems():
            if hasattr(var,'units'):
                if (var  in _possibleaxis or var.units in _possibleaxisunits or var.units.split(" ")[0]  in _possibleaxisunits or hasattr(var,'positive')) and name not in space_time_coord_var:
                    space_time_non_coord_var.append(name)
                    for every in var.dimensions:
                        space_time_non_coord_var_dim.append(every)

        #Looks to ensure that every dimension of each variable that is a space-time dimension has associated coordinate variables
        for name,var in ds.dataset.variables.iteritems():
            valid = ''
            for each in var.dimensions:
                if each in space_time_non_coord_var_dim:
                    valid = False
                    dim_name = each
                    break
                elif each in space_time_coord_var:
                    valid = True

            if valid == False :
                ret_val.append(Result(BaseCheck.MEDIUM, \
                               valid, \
                               ('var', name, 'check_independent_axis_dimensions'),['The %s dimension for the variable %s does not have an associated coordinate variable, but is a Lat/Lon/Time/Height dimension.'%(dim_name,name)]))

            if valid == True and name not in space_time_coord_var:
                ret_val.append(Result(BaseCheck.MEDIUM, \
                               valid, \
                               ('var', name, 'check_independent_axis_dimensions')))
        return ret_val


    def check_two_dimensional(self, ds):
        """
        5.2 The latitude and longitude coordinates of a horizontal grid that was
        not defined as a Cartesian product of latitude and longitude axes, can
        sometimes be represented using two-dimensional coordinate variables.
        These variables are identified as coordinates by use of the coordinates
        attribute.

        For each variable, if the variable has a coordinates attribute:
          for each coordinate defined, verify that the coordinate:
            is either a coordinate variable OR comprises coordinate variables
          
        """

        ret_val = []
        for name,var in ds.dataset.variables.iteritems():
            g = NCGraph(ds.dataset, name, var)
            #Determine if 2-D coordinate variables (Lat and Lon are of shape (i,j)
            for each in g.coords:
                valid = g.coords[each].ndim == 2

            if len(g.coords) == 2 and valid:
                #------------------------------------------------------------
                # Check all the dims are coordinate variables
                #------------------------------------------------------------
                valid_dims = True
                reasoning = []
                for dim in g.dims.iterkeys():
                    if dim not in ds.dataset.variables:
                        valid_2d = False
                        reasoning.append("Variable %s's dimension %s is not a coordinate variable" % (name, dim))
    
                result = Result(BaseCheck.HIGH,                             \
                                valid_dims,                                 \
                                ('var', name, '2d_hgrid_valid_dimensions'), \
                                reasoning)
                ret_val.append(result)
                
                #------------------------------------------------------------
                # Check that the coordinates are correct
                #------------------------------------------------------------
                valid_2d = True
                reasoning = []
                for cname, coord in g.coords.iteritems():
                    if coord is None:
                        valid_2d = False
                        reasoning.append("Variable %s's coordinate, %s, is not a coordinate or auxiliary variable" %(name, cname))
                        continue
                    for dim in coord.dims.iterkeys():
                        if dim not in g.dims:
                            valid_2d = False
                            reasoning.append("Variable %s's coordinate, %s, does not share dimension %s with the variable" % (name, cname, dim))
                result = Result(BaseCheck.MEDIUM,                   \
                                valid_2d,                           \
                                ('var', name, 'valid_coordinates'), \
                                reasoning)
                ret_val.append(result)
    
                
                #------------------------------------------------------------
                # Can make lat/lon?
                #------------------------------------------------------------
    
                lat_check = False
                lon_check = False

                for cname, coord in g.coords.iteritems():

                    if coord.units in ['degrees_north', 'degree_north', 'degrees_N', 'degree_N', 'degreesN', 'degreeN']:
                        lat_check = True
                    elif coord.units in ['degrees_east', 'degree_east', 'degrees_E', 'degree_E', 'degreesE', 'degreeE']:
                        lon_check = True
    
                result = Result(BaseCheck.HIGH,          \
                                lat_check and lon_check, \
                                ('var', name, 'lat_lon_correct'))
                ret_val.append(result)
            else:
                continue # Not a 2d horizontal grid


        return ret_val


    def check_reduced_horizontal_grid(self, ds):
        """
        5.3 A "reduced" longitude-latitude grid is one in which the points are
        arranged along constant latitude lines with the number of points on a
        latitude line decreasing toward the poles.

        Recommend that this type of gridded data be stored using the compression
        scheme described in Section 8.2, "Compression by Gathering". The
        compressed latitude and longitude auxiliary coordinate variables are
        identified by the coordinates attribute.

        """
        ret_val = []
        coord_vars = self._get_coord_vars(ds)

        for name, var in ds.dataset.variables.iteritems():
            if name in coord_vars:
                continue
            if not hasattr(var, 'coordinates'):
                continue

            valid = True
            reasoning = []

            coords = var.coordinates.split(' ')
            for coord in coords:
                if coord not in ds.dataset.variables:
                    valid = False
                    reasoning.append("Coordinate %s is not a proper variable" % coord)
                    continue

                for dim_name in ds.dataset.variables[coord].dimensions:

                    if dim_name not in var.dimensions:
                        valid = False
                        reasoning.append("Coordinate %s's dimension, %s, is not a dimension of %s" %(coord, dim_name, name))
                        continue

                    if dim_name not in coord_vars:
                        valid = False
                        reasoning.append("Coordinate %s's dimension, %s, is not a coordinate variable" % (coord, dim_name))
                        continue

                    dim = ds.dataset.variables[dim_name]
                    if not hasattr(dim, 'compress'):
                        valid = False
                        reasoning.append("Coordinate %s's dimension, %s, does not define compress" % (coord, dim_name))
                        continue

                    compress_dims = dim.compress.split(' ')
                    for cdim in compress_dims:
                        if cdim not in ds.dataset.dimensions:
                            valid = False
                            reasoning.append("Dimension %s compresses non-existent dimension, %s" % (dim_name, cdim))
                            continue
            result = Result(BaseCheck.MEDIUM,                            \
                            valid,                                       \
                            ('var', name, 'is_reduced_horizontal_grid'), \
                            reasoning)
            if valid == True:
                ret_val.append(result)
        
        return ret_val

    # grid mapping dictionary, appendix F
    grid_mapping_dict = {
        'albers_conical_equal_area': [('longitude_of_central_meridian', 'latitude_of_projection_origin', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate')],
        'azimuthal_equidistant': [('longitude_of_projection_origin', 'latitude_of_projection_origin', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate')],
        'lambert_cylindrical_equal_area': [('longitude_of_central_meridian', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate'), ('standard_parallel','scale_factor_at_projection_origin')],
        'lambert_azimuthal_equal_area': [('longitude_of_projection_origin', 'latitude_of_projection_origin', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate')],
        'lambert_conformal_conic': [('standard_parallel', 'longitude_of_central_meridian', 'latitude_of_projection_origin', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate')],
        'latitude_longitude': [(),(),('longitude', 'latitude')],
        'mercator': [('longitude_of_projection_origin', 'false_easting', 'false_northing'),(),('projection_x_coordinate','projection_y_coordinate'), ('standard_parallel', 'scale_factor_at_projection_origin')],
        'orthographic': [('longitude_of_projection_origin', 'latitude_of_projection_origin', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate')],
        'polar_stereographic': [('straight_vertical_longitude_from_pole', 'latitude_of_projection_origin', 'false_easting', 'false_northing'),(),('projection_x_coordinate', 'projection_y_coordinate'), ('standard_parallel', 'scale_factor_at_projection_origin')],
        'rotated_latitude_longitude': [('grid_north_pole_latitude', 'grid_north_pole_longitude'),('north_pole_grid_longitude'),('grid_latitude', 'grid_longitude')],
        'stereographic':[('longitude_of_projection_origin', 'latitude_of_projection_origin', 'scale_factor_at_projection_origin', 'false_easting', 'false_northing'),(),('projection_x_coordinate', 'projection_y_coordinate')],
        'transverse_mercator': [('scale_factor_at_central_meridian', 'longitude_of_central_meridian', 'latitude_of_projection_origin', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate')],
        'vertical_perspective': [('longitude_of_projection_origin', 'latitude_of_projection_origin', 'perspective_point_height', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate')]
    }


    def check_horz_crs_grid_mappings_projections(self, ds):
        """
        5.6 When the coordinate variables for a horizontal grid are not
        longitude and latitude, it is required that the true latitude and
        longitude coordinates be supplied via the coordinates attribute. If in
        addition it is desired to describe the mapping between the given
        coordinate variables and the true latitude and longitude coordinates,
        the attribute grid_mapping may be used to supply this description.

        This attribute is attached to data variables so that variables with
        different mappings may be present in a single file. The attribute takes
        a string value which is the name of another variable in the file that
        provides the description of the mapping via a collection of attached
        attributes. This variable is called a grid mapping variable and is of
        arbitrary type since it contains no data. Its purpose is to act as a
        container for the attributes that define the mapping.

        The one attribute that all grid mapping variables must have is
        grid_mapping_name which takes a string value that contains the mapping's
        name. The other attributes that define a specific mapping depend on the
        value of grid_mapping_name. The valid values of grid_mapping_name along
        with the attributes that provide specific map parameter values are
        described in Appendix F, Grid Mappings.

        When the coordinate variables for a horizontal grid are longitude and
        latitude, a grid mapping variable with grid_mapping_name of
        latitude_longitude may be used to specify the ellipsoid and prime
        meridian.


        In order to make use of a grid mapping to directly calculate latitude
        and longitude values it is necessary to associate the coordinate
        variables with the independent variables of the mapping. This is done by
        assigning a standard_name to the coordinate variable. The appropriate
        values of the standard_name depend on the grid mapping and are given in
        Appendix F, Grid Mappings.  
        """ 
        
        
        ret_val = []
        reasoning = []

        
        for name, var in ds.dataset.variables.iteritems():
            valid_mapping_count = 0
            total_mapping_count = 0
            if hasattr(var, 'grid_mapping_name'):
                total_mapping_count = 1
                
                mapping = getattr(var, 'grid_mapping_name', '')
                if mapping in self.grid_mapping_dict.iterkeys():
                    valid_mapping_count = valid_mapping_count +1
                else:
                    reasoning.append('The grid_mapping_name attribute is not an accepted value.  See Appendix F.')
                    
                    
                for each in self.grid_mapping_dict[mapping][0]:
                    total_mapping_count = total_mapping_count + 1
                    if each in dir(var):
                        valid_mapping_count = valid_mapping_count +1
                    else:
                        reasoning.append('The map parameters are not accepted values.  See Appendix F.')
                
                if len(self.grid_mapping_dict[mapping]) >=4:
                    for each in self.grid_mapping_dict[mapping][3:]:
                        every_flag = 0
                        total_mapping_count = total_mapping_count + 1
                        for every in each:
                            if every in dir(var):
                                valid_mapping_count = valid_mapping_count + 1
                                every_flag = every_flag +1
                        
                        if every_flag == 0:
                            reasoning.append('Neither of the "either/or" parameters are present')
                        if every_flag == 2:
                            valid_mapping_count = valid_mapping_count - 2
                
                total_mapping_count = total_mapping_count + len(self.grid_mapping_dict[mapping][2])
                for name_again, var_again in ds.dataset.variables.iteritems():
                    if hasattr(var_again,'standard_name'):
                        if var_again.standard_name in self.grid_mapping_dict[mapping][2]:
                            valid_mapping_count = valid_mapping_count + 1
                
                result = Result(BaseCheck.MEDIUM,                            \
                        (valid_mapping_count, total_mapping_count),                                       \
                        ('var', name, 'horz_crs_grid_mappings_projections'), \
                        reasoning)

                ret_val.append(result)

        return ret_val
                


    def check_scalar_coordinate_system(self, ds):

        """
        5.7 When a variable has an associated coordinate which is single-valued, that coordinate may be represented as a
        scalar variable. Since there is no associated dimension these scalar coordinate variables should be attached to a
        data variable via the coordinates attribute.

        Once a name is used for a scalar coordinate variable it can not be used for a 1D coordinate variable. For this
        reason we strongly recommend against using a name for a scalar coordinate variable that matches the name of any
        dimension in the file.
        """
        ret_val = []

        for name, var in ds.dataset.variables.iteritems():
            valid_scalar_coordinate_var = 0
            total_scalar_coordinate_var = 0
            reasoning = []

            if not hasattr(var, 'coordinates'):
                continue

            for coordinate in getattr(var, 'coordinates', '').split(" "):
                if coordinate in ds.dataset.variables:
                    if ds.dataset.variables[coordinate].shape == (1,):
                        total_scalar_coordinate_var += 1
                        if coordinate not in ds.dataset.dimensions.keys():
                            valid_scalar_coordinate_var += 1
                        else:
                            reasoning.append('Scalar coordinate var (%s) of var (%s) is correct size but is present in the dimensions list, which is not allowed.'% (coordinate, name))

            if total_scalar_coordinate_var > 0:
                result = Result(BaseCheck.MEDIUM,
                                (valid_scalar_coordinate_var, total_scalar_coordinate_var),
                                ('var', name, 'scalar_coordinates'),
                                reasoning)
                ret_val.append(result)

        return ret_val

    ###############################################################################
    #
    # CHAPTER 6: Labels and Alternative Coordinates
    #
    ###############################################################################

    def check_geographic_region(self, ds):
        """
        6.1.1 When data is representative of geographic regions which can be identified by names but which have complex
        boundaries that cannot practically be specified using longitude and latitude boundary coordinates, a labeled
        axis should be used to identify the regions. 

        Recommend that the names be chosen from the list of standardized region names whenever possible. To indicate
        that the label values are standardized the variable that contains the labels must be given the standard_name
        attribute with the value region.
        """
        ret_val = []
        reasoning = []
        region_list = [
                        'africa',
                        'antarctica',
                        'arabian_sea',
                        'aral_sea',
                        'arctic_ocean',
                        'asia',
                        'atlantic_ocean',
                        'australia',
                        'baltic_sea',
                        'barents_opening',
                        'barents_sea',
                        'beaufort_sea',
                        'bellingshausen_sea',
                        'bering_sea',
                        'bering_strait',
                        'black_sea',
                        'canadian_archipelago',
                        'caribbean_sea',
                        'caspian_sea',
                        'central_america',
                        'chukchi_sea',
                        'contiguous_united_states',
                        'denmark_strait',
                        'drake_passage',
                        'east_china_sea',
                        'english_channel',
                        'eurasia',
                        'europe',
                        'faroe_scotland_channel',
                        'florida_bahamas_strait',
                        'fram_strait',
                        'global',
                        'global_land',
                        'global_ocean',
                        'great_lakes',
                        'greenland',
                        'gulf_of_alaska',
                        'gulf_of_mexico',
                        'hudson_bay',
                        'iceland_faroe_channel',
                        'indian_ocean',
                        'indonesian_throughflow',
                        'indo_pacific_ocean',
                        'irish_sea',
                        'lake_baykal',
                        'lake_chad',
                        'lake_malawi',
                        'lake_tanganyika',
                        'lake_victoria',
                        'mediterranean_sea',
                        'mozambique_channel',
                        'north_america',
                        'north_sea',
                        'norwegian_sea',
                        'pacific_equatorial_undercurrent',
                        'pacific_ocean',
                        'persian_gulf',
                        'red_sea',
                        'ross_sea',
                        'sea_of_japan',
                        'sea_of_okhotsk',
                        'south_america',
                        'south_china_sea',
                        'southern_ocean',
                        'taiwan_luzon_straits',
                        'weddell_sea',
                        'windward_passage',
                        'yellow_sea'
                        ]

        
        for name, var in ds.dataset.variables.iteritems():
            if getattr(var, 'standard_name', '') == 'region':
                if ''.join(var[:]).lower() in region_list:
                    reasoning.append('The Region Value is from the allowable list.')
                    result = Result(BaseCheck.LOW,                            \
                            True,                                       \
                            ('var', name, 'geographic_region'), \
                            reasoning)
                else:
                    reasoning.append('The Region Value is not from the allowable list.')
                    result = Result(BaseCheck.LOW,                            \
                            False,                                       \
                            ('var', name, 'geographic_region'), \
                            reasoning)
                ret_val.append(result)
        return ret_val

    def check_alternative_coordinates(self, ds):
        """
        6.2 In some situations a dimension may have alternative sets of coordinates values. Since there can only be
        one coordinate variable for the dimension (the variable with the same name as the dimension), any alternative
        sets of values have to be stored in auxiliary coordinate variables. For such alternative coordinate variables,
        there are no mandatory attributes, but they may have any of the attributes allowed for coordinate variables.
        """
        ret_val = []
        reasoning = []
        valid_alt_coordinate_var = 0
        total_alt_coordinate_var = 0


        for name, var in ds.dataset.variables.iteritems():
            valid_alt_coordinate_var = 0
            total_alt_coordinate_var = 0
            if getattr(var, 'coordinates', ''):
                for coordinate in getattr(var, 'coordinates', '').split(' '):
                    if coordinate in ds.dataset.variables and coordinate not in ds.dataset.dimensions:
                        reasoning.append('The Alternative Coordinate system for variable %s coordinate %s is not derived from the Coordinate Variables and is not Dimensionalized by them.' %(name, coordinate))
                        total_alt_coordinate_var = total_alt_coordinate_var + 1
                        valid_alt_coordinate_var = valid_alt_coordinate_var + 1
                        continue
                    elif coordinate not in ds.dataset.variables:
                        reasoning.append('The Alternative Coordinate system for variable %s coordinate %s is derived from the Coordinate Variables and is Dimensionalized by them.' %(name, coordinate))
                        total_alt_coordinate_var = total_alt_coordinate_var + 1


                result = Result(BaseCheck.MEDIUM,                            \
                            (valid_alt_coordinate_var,total_alt_coordinate_var),                                       \
                            ('var', name, 'alternative_coordinates'), \
                            reasoning)
                ret_val.append(result)

        return ret_val

    ###############################################################################
    #
    # CHAPTER 7: Data Representative of Cells
    #
    ###############################################################################

    def check_cell_boundaries(self, ds):
        """
        7.1 To represent cells we add the attribute bounds to the appropriate coordinate variable(s). The value of bounds
        is the name of the variable that contains the vertices of the cell boundaries. We refer to this type of variable as
        a "boundary variable." A boundary variable will have one more dimension than its associated coordinate or auxiliary
        coordinate variable. The additional dimension should be the most rapidly varying one, and its size is the maximum
        number of cell vertices. 

        Applications that process cell boundary data often times need to determine whether or not adjacent cells share an
        edge. In order to facilitate this type of processing the following restrictions are placed on the data in boundary
        variables:

        Bounds for 1-D coordinate variables

            For a coordinate variable such as lat(lat) with associated boundary variable latbnd(x,2), the interval endpoints
            must be ordered consistently with the associated coordinate, e.g., for an increasing coordinate, lat(1) > lat(0)
            implies latbnd(i,1) >= latbnd(i,0) for all i

            If adjacent intervals are contiguous, the shared endpoint must be represented indentically in each instance where
            it occurs in the boundary variable. For example, if the intervals that contain grid points lat(i) and lat(i+1) are
            contiguous, then latbnd(i+1,0) = latbnd(i,1).
            
        Bounds for 2-D coordinate variables with 4-sided cells
            
            In the case where the horizontal grid is described by two-dimensional auxiliary coordinate variables in latitude
            lat(n,m) and longitude lon(n,m), and the associated cells are four-sided, then the boundary variables are given
            in the form latbnd(n,m,4) and lonbnd(n,m,4), where the trailing index runs over the four vertices of the cells. 

        Bounds for multi-dimensional coordinate variables with p-sided cells

            In all other cases, the bounds should be dimensioned (...,n,p), where (...,n) are the dimensions of the auxiliary
            coordinate variables, and p the number of vertices of the cells. The vertices must be traversed anticlockwise in the
            lon-lat plane as viewed from above. The starting vertex is not specified.
        """
        ret_val = []
        reasoning = []
        valid = ' '

        for cvar, bvar in self._find_boundary_vars(ds).iteritems():
                valid = True
                if bvar.ndim !=cvar.ndim + 1:
                    valid = False
                    reasoning.append('The number of dimensions of the Coordinate Variable is %s, but the number of dimensions of the Boundary Variable is %s.'%(cvar.ndim, bvar.ndim))

                result = Result(BaseCheck.MEDIUM,                          
                            valid,                                       
                            ('var', cvar._name, 'cell_boundaries'), 
                            reasoning)
                ret_val.append(result)
                reasoning = []

        return ret_val

    def check_cell_measures(self, ds):
        """
        7.2 To indicate extra information about the spatial properties of a variable's grid cells, a cell_measures attribute may
        be defined for a variable. This is a string attribute comprising a list of blank-separated pairs of words of the form
        "measure: name". "area" and "volume" are the only defined measures.

        The "name" is the name of the variable containing the measure values, which we refer to as a "measure variable". The
        dimensions of the measure variable should be the same as or a subset of the dimensions of the variable to which they are
        related, but their order is not restricted.

        The variable must have a units attribute and may have other attributes such as a standard_name.
        """
        ret_val = []
        reasoning = []
        paragraph = []
        for name, var in ds.dataset.variables.iteritems():
            for dim in var.dimensions:
                if getattr(var, 'cell_measures', ''):
                    measures = getattr(var,'coordinates','')
                    measures = measures.split(': ')
                    if measures[0] not in ['area', 'volume']:
                        reasoning.append("The 'measures' field is not equal to 'area' or 'volume'.")
                        return Result(BaseCheck.MEDIUM,                            \
                                    False,                                       \
                                    ('var', name, 'cell_measures'), \
                                    reasoning)
                    for every, attri in ds.dataset.variables.iteritems():
                        if every == measures[1]:
                            for dimi in attri.dimensions:
                                if dimi in var.dimensions:
                                    reasoning.append('The measure variable dimensions are a set or subset of the cell_measure variable.')
                                    valid = True
                                else:
                                    reasoning.append('The measure variable dimensions are not a set or subset of the cell_measure variable.')
                                    valid == False

    
                    result = Result(BaseCheck.MEDIUM,                            \
                                    valid,                                       \
                                    ('var', name, 'cell_measures'), \
                                    reasoning)
                    ret_val.append(result)


        return ret_val


    def check_cell_methods(self, ds):
        """
        7.3 To describe the characteristic of a field that is represented by cell values, we define the cell_methods attribute
        of the variable. This is a string attribute comprising a list of blank-separated words of the form "name: method". Each
        "name: method" pair indicates that for an axis identified by name, the cell values representing the field have been
        determined or derived by the specified method.

        name can be a dimension of the variable, a scalar coordinate variable, a valid standard name, or the word "area"

        values of method should be selected from the list in Appendix E, Cell Methods, which includes point, sum, mean, maximum,
        minimum, mid_range, standard_deviation, variance, mode, and median. Case is not significant in the method name. Some
        methods (e.g., variance) imply a change of units of the variable, as is indicated in Appendix E, Cell Methods.

        Because the default interpretation for an intensive quantity differs from that of an extensive quantity and because this
        distinction may not be understood by some users of the data, it is recommended that every data variable include for each
        of its dimensions and each of its scalar coordinate variables the cell_methods information of interest (unless this
        information would not be meaningful). It is especially recommended that cell_methods be explicitly specified for each
        spatio-temporal dimension and each spatio-temporal scalar coordinate variable.
        """
        
        
        _areatype_names     = ["bare_ground",
                                "all_area_types", 
                                "burnt_vegetation",
                                "c3_plant_functional_types",
                                "c4_plant_functional_types",
                                "clear_sky", 
                                "cloud", 
                                "crops",
                                "floating_ice", 
                                "ice_free_land",
                                "ice_free_sea", 
                                "lake_ice_or_sea_ice", 
                                "land", 
                                "land_ice",
                                "natural_grasses",
                                "pastures",
                                "primary_deciduous_trees",
                                "primary_evergreen_trees",
                                "sea",
                                "sea_ice", 
                                "secondary_deciduous_trees",
                                "secondary_evergreen_trees",
                                "shrubs"
                                "snow", 
                                "trees"
                                "vegetation"]
        methods = [ 'point',
                    'sum',
                    'mean',
                    'maximum',
                    'minimum',
                    'mid_range',
                    'standard_deviation',
                    'variance',
                    'mode',
                    'median']

        ret_val = []
        reasoning = []
        paragraph = ''
        named = ''
        
        names = list(ds.dataset.variables.iterkeys())
        for name, var in ds.dataset.variables.iteritems():
            named_dict = OrderedDict()
            if getattr(var, 'cell_methods', '') :
                method = getattr(var, 'cell_methods', '')
                #Parse the cell_methods format.
                paragraph = re.split(r"(?<=:)\s|(?<=\s)\(",method)
                
                i = 0
                while i < (len(paragraph)):
                    if paragraph[i][-10:] == " interval:":
                        paragraph[i] = paragraph[i][0:-10]
                        paragraph.insert(i+1,"interval:")
                        i = i+1
                    if paragraph[i][-6:] == " area:":
                        paragraph[i] = paragraph[i][0:-6]
                        paragraph.insert(i+1,"area:")
                        i = i+1
                    if paragraph[i][-9:] == " comment:":
                        paragraph[i] = paragraph[i][0:-9]
                        paragraph.insert(i+1,"comment:")
                        i = i+1
                    i = i+1
        
                
                dict_count = 0
                #Create a dict that holds each statement in the cell_method field
                for i in range(len(paragraph)-1):
        
                    if paragraph[i][-1] == ":":
                        named = named +paragraph[i]
                    if paragraph[i+1][-1] != ":":
                        named_dict[str(dict_count)+named] = []
                    if paragraph[i][-1] != ":":
                        named_dict[str(dict_count)+named].append(paragraph[i])
        
                        dict_count = dict_count+1
                        named = ''
                named_dict[str(dict_count)+named].append(paragraph[i+1])
                
                
    
                valid_name_count = 0
                total_name_count = 0
                
                #Checks if the name value of the 'name: method' pair is either a dimension of the variable, a standard_name, or a scalar_variable (which is         trumped by the standard_name requirement_
                
                for each in named_dict.iterkeys():
                    title = each[1:].split(':')
                    for i in range(len(title)):
                        if title[i].lower() in _areatype_names:
                            valid_name_count = valid_name_count +1
                            reasoning.append('The name field appears in the allowable types.')
                        else:
                            reasoning.append('The name field does not appear in the allowable types.')
                                
                for var_dim in ds.dataset.variables[name].dimensions:
                    for each in named_dict.iterkeys():
                        title = each[1:].split(':')
                        for i in range(len(title)):
                            if title[i].lower() == var_dim:
                                valid_name_count = valid_name_count +1
                                reasoning.append('The name field matches the dimension.') 
                            else:
                                reasoning.append('The name field does not match the dimension.') 

                for each in named_dict.iterkeys():
                    title = each[1:].split(':')
                
                    for i in range(len(title)):
                        if title[i] != '':
                            total_name_count = total_name_count +1
                        if title[i].lower() in ["interval", "area", "comment"] :
                            reasoning.append('The name field matches the reserved words "interval", "area", or "comment".')
                            valid_name_count = valid_name_count +1
                        else:
                            reasoning.append('The name field does not match the reserved words "interval", "area", or "comment".')

                result = Result(BaseCheck.MEDIUM,                            \
                        (valid_name_count, total_name_count),                                       \
                        ('var', name, 'cell_methods_name'), \
                        reasoning)
                ret_val.append(result)

                reasoning = []
                            
                #Checks if the method value of the 'name: method' pair is acceptable
                methods = ['point', 'sum', 'mean', 'maximum', 'minimum', 'mid_range', 'standard_deviation', 'variance', 'mode', 'median']            
                valid_method_count = 0
                total_method_count = 0
                
                for each in named_dict.iterkeys():
                    title = each[1:].split(':')
                
                    for i in range(len(title)):
                        if title[i] not in ['interval', 'comment', 'area', '']:
                            if title[i] != '':
                                total_method_count = total_method_count +1
                            if named_dict[each][0].strip() in methods:
                                valid_method_count = valid_method_count+1
                                reasoning.append('The method field matches a valid method value.')
                            else:
                                reasoning.append('The method field does not match a valid method value.')
                result = Result(BaseCheck.MEDIUM,                            \
                        (valid_method_count, total_method_count),                                       \
                        ('var', name, 'cell_methods_method'), \
                        reasoning)
                ret_val.append(result)            
                
                #Checks the format of the interval field
                reasoning = []
                valid_interval_count = 0
                total_interval_count = 0
                
                for each in named_dict.iterkeys():
                    title = each[1:].split(':')
                
                    for i in range(len(title)):
                        if title[i] == 'interval':
                            total_interval_count = total_interval_count +1
                            if len(named_dict[each][0].split(" ")) == 2:
                                valid_interval_count = valid_interval_count+1
                                reasoning.append('The "interval: value units" format is the correct length.')
                            else:
                                reasoning.append('The "interval: value units" format is not the correct length.')

                result = Result(BaseCheck.MEDIUM,                            \
                        (valid_interval_count, total_interval_count),                                       \
                        ('var', name, 'cell_methods_interval'), \
                        reasoning)
                ret_val.append(result)    

                #Checks the 'method where' formats
                reasoning = []
                valid_area_count = 0
                total_area_count = 0
                for each in named_dict.iterkeys():
                    title = each[1:].split(':')
                
                    for i in range(len(title)):
                        if title[i] == 'area':
                            total_area_count = total_area_count +1
                            area_data = named_dict[each][0].split(" ")
                            if len(area_data) == 4:
                                if area_data[0] in methods and area_data[1] == "where" and area_data[2] in _areatype_names:
                                    valid_area_count = valid_area_count+1
                                    reasoning.append('The "name: method where _areatype_names" format is correct.')
                                else:
                                    reasoning.append('The "name: method where _areatype_names" format is not correct.')
                            elif len(area_data) == 6:
                                if area_data[0] in methods and area_data[1] == "where" and area_data[2] in _areatype_names and area_data[3] == "over" and       area_data[4] in _areatype_names :
                                    valid_area_count = valid_area_count+1
                                    reasoning.append('The "name: method where type over _areatype_names" format is correct.')
                                else:
                                    reasoning.append('The "name: method where type over _areatype_names" format is not correct.')

                result = Result(BaseCheck.MEDIUM,                            \
                        (valid_area_count, total_area_count),                                       \
                        ('var', name, 'cell_methods_area'), \
                        reasoning)
                ret_val.append(result)   
        
                #Checks the no coordinate case
                reasoning = []
                valid_no_coord_count = 0
                total_no_coord_count = 0
                for each in named_dict.iterkeys():
                    title = each[1:].split(':')
                
                    for i in range(len(title)):
                        if title[i].lower() in self._std_names and title[i].lower() not in ds.dataset.variables[name].dimensions and title[i].lower() not in getattr(var,"coordinates",""):
                            if title[i] != '':
                                total_no_coord_count = total_no_coord_count +1
                            if named_dict[each][0].strip() in methods:
                                valid_no_coord_count = valid_no_coord_count+1
                                reasoning.append('The method is in the a value provided in the allowable method list.')
                            else:
                                reasoning.append('The method is not in the a value provided in the allowable method list.')
                                
                result = Result(BaseCheck.MEDIUM,                            \
                        (valid_no_coord_count, total_no_coord_count),                                       \
                        ('var', name, 'cell_methods_no_coord'), \
                        reasoning)
                ret_val.append(result)   
        #Checks the Climatology Variables - 7.4                               
        reasoning = []
        paragraph = []
        total_climate_count = 0
        valid_climate_count = 0
        for name, var in ds.dataset.variables.iteritems():
            if getattr(var, 'climatology', ''):
                climate_dim = ds.dataset.variables[name].dimensions
                clim_method = getattr(var, 'climatology', '')
                
        
                
                for each in climate.split(" "):
                    paragraph.append(each)

                total_climate_count = total_climate_count+ 1
                for name_again, var_again in ds.dataset.variables.iteritems():
                    if getattr(var_again,"cell_methods",""):
                        climate = getattr(var, 'cell_methods', '')
                        name_dim = ds.dataset.variables[name_again].dimensions
                        if len(climate_dim)>0:
                            if climate_dim[0] in name_dim:
                                case1 = re.search(r"time: \w* within years time: \w* over years",climate)
                                case2 = re.search(r"time: \w* within days time: \w* over days$",climate)
                                case3 = re.search(r"time: \w* within days time: \w* over days time: \w* over years",climate)
                        
                        if (case1 or case2 or case3) and len(ds.dataset.variables[clim_method].shape) == 2 and ds.dataset.variables[clim_method].shape[1] == 2 and ds.dataset.variables[clim_method].shape[0] == ds.dataset.variables[name_again].shape[0] :
                            
                            valid_climate_count = 1
                            reasoning.append('The "time: method within years/days over years/days" format is correct.')
                        if not (case1 or case2 or case3):
                            reasoning.append('The "time: method within years/days over years/days" format is not correct.')

                        if not (len(ds.dataset.variables[clim_method].shape) == 2 and ds.dataset.variables[clim_method].shape[1] == 2 and ds.dataset.variables[clim_method].shape[0] == ds.dataset.variables[name_again].shape[0]):
                            reasoning.append('The dimensions of the climatology varaible is incorrect.')


                result = Result(BaseCheck.MEDIUM,                            \
                        (valid_climate_count, total_climate_count),                                       \
                        ('var', name, 'cell_methods_climatology'), \
                        reasoning)
                ret_val.append(result)
        
        return ret_val



    #def check_cell_methods_for_multi_axes(self, ds):
        """
        7.3.1 If a data value is representative of variation over a combination of axes, a single method should be prefixed by the
        names of all the dimensions involved (listed in any order, since in this case the order must be immaterial). 
        
        There is no way to check this.  A warning should be posted explaining this method to the user!"

        """



    #def check_spacing_and_extra_info(self, ds):
        """
        7.3.2 To indicate more precisely how the cell method was applied, extra information may be included in parentheses ()
        after the identification of the method. This information includes standardized and non-standardized parts.

        The only standardized information is to provide the typical interval between the original data values to which the method
        was applied, in the situation where the present data values are statistically representative of original data values which
        had a finer spacing.

        The syntax is (interval: value unit), where value is a numerical value and unit is a string that can be recognized by
        UNIDATA's Udunits package.

        If the cell method applies to a combination of axes, they may have a common original interval. Alternatively, they may have
        separate intervals, which are matched to the names of axes by position.

        If there is both standardized and non-standardized information, the non-standardized follows the standardized information
        and the keyword comment:. If there is no standardized information, the keyword comment: should be omitted.

        A dimension of size one may be the result of "collapsing" an axis by some statistical operation, for instance by
        calculating a variance from time series data. We strongly recommend that dimensions of size one be retained (or scalar
        coordinate variables be defined) to enable documentation of the method (through the cell_methods attribute) and its
        domain (through the cell_bounds attribute).
        """



    #def check_stats_applying_to_portions_of_cells(self, ds):
        """
        7.3.3 By default, the statistical method indicated by cell_methods is assumed to have been evaluated over the entire
        horizontal area of the cell. Sometimes, however, it is useful to limit consideration to only a portion of a cell.

        One of two conventions may be used.

        The first convention is a method that can be used for the common case of a single area-type. In this case, the
        cell_methods attribute may include a string of the form "name: method where type".

        The second convention is the more general. In this case, the cell_methods entry is of the form "name: method where
        _areatype_names". Here _areatype_names is a string-valued auxiliary coordinate variable or string-valued scalar coordinate variable
        with a standard_name of area_type. The variable _areatype_names contains the name(s) of the selected portion(s) of the grid
        cell to which the method is applied. 

        If the method is mean, various ways of calculating the mean can be distinguished in the cell_methods attribute with
        a string of the form "mean where type1 [over type2]". Here, type1 can be any of the possibilities allowed for _areatype_names
        or type (as specified in the two paragraphs preceding above Example). The same options apply to type2, except it is
        not allowed to be the name of an auxiliary coordinate variable with a dimension greater than one (ignoring the
        dimension accommodating the maximum string length)
        """

    #def check_cell_methods_with_no_coords(self, ds):
        """
        7.3.4 To provide an indication that a particular cell method is relevant to the data without having to provide a
        precise description of the corresponding cell, the "name" that appears in a "name: method" pair may be an
        appropriate standard_name (which identifies the dimension) or the string, "area" (rather than the name of a scalar
        coordinate variable or a dimension with a coordinate variable). This convention cannot be used, however, if the name
        of a dimension or scalar coordinate variable is identical to name. 

        Recommend that whenever possible, cell bounds should be supplied by giving the variable a dimension of size one
        and attaching bounds to the associated coordinate variable.
        """
        

    #def check_climatological_statistics(self, ds):
        """
        7.4 A climatological time coordinate variable does not have a bounds attribute. Instead, it has a climatology
        attribute, which names a variable with dimensions (n,2), n being the dimension of the climatological time axis.
        Using the units and calendar of the time coordinate variable, element (i,0) of the climatology variable specifies
        the beginning of the first subinterval and element (i,1) the end of the last subinterval used to evaluate the
        climatological statistics with index i in the time dimension. The time coordinates should be values that are
        representative of the climatological time intervals, such that an application which does not recognise climatological
        time will nonetheless be able to make a reasonable interpretation.

        Valid values of the cell_methods attribute must be in one of the forms from the following list.

        - time: method1 within years   time: method2 over years
        - time: method1 within days   time: method2 over days
        - time: method1 within days   time: method2 over days   time: method3 over years

        The methods which can be specified are those listed in Appendix E, Cell Methods and each entry in the cell_methods
        attribute may also, contain non-standardised information in parentheses after the method. 
        """
        

    ###############################################################################
    #
    # CHAPTER 8: Reduction of Dataset Size
    #
    ###############################################################################

    def check_packed_data(self, ds):
        """
        8.1 Simple packing may be achieved through the use of the optional NUG defined attributes scale_factor and
        add_offset. After the data values of a variable have been read, they are to be multiplied by the scale_factor,
        and have add_offset added to them. 

        The units of a variable should be representative of the unpacked data.

        If the scale_factor and add_offset attributes are of the same data type as the associated variable, the unpacked
        data is assumed to be of the same data type as the packed data. However, if the scale_factor and add_offset
        attributes are of a different data type from the variable (containing the packed data) then the unpacked data
        should match the type of these attributes, which must both be of type float or both be of type double. An additional
        restriction in this case is that the variable containing the packed data must be of type byte, short or int. It is
        not advised to unpack an int into a float as there is a potential precision loss.

        When data to be packed contains missing values the attributes that indicate missing values (_FillValue, valid_min,
        valid_max, valid_range) must be of the same data type as the packed data.
        """
        
        ret_val = []
        reasoning = []
        names = list(ds.dataset.variables.iterkeys())
        for name, var in ds.dataset.variables.iteritems():
            if hasattr(var, 'add_offset') or hasattr(var, 'scale_factor'):
                valid = False
                if hasattr(var, 'add_offset') and hasattr(var, 'scale_factor'):
                    offset = getattr(var, 'add_offset', '')
                    scale = getattr(var, 'scale_factor', '')

                    data_size = np.s_[ds.dataset.variables[name].size]
                    data_type_check = np.reshape(ds.dataset.variables[name],data_size)[0]
                    #, type(ds.dataset.variables[name][0,0,0]) ,type(scale), type(offset)
                    for x in xrange(data_size):
                        if type(np.reshape(ds.dataset.variables[name],data_size)[x]) != 'numpy.ma.core.MaskedConstant':
                            data_type_check = np.reshape(ds.dataset.variables[name],data_size)[x]
                            break
                    if type(data_type_check) == type(scale) == type(offset):
                        valid = True
                        #reasoning.append("'add_offset' and 'scale_factor' both exist and match data type to the variable.")
                        result = Result(BaseCheck.MEDIUM,                            \
                                valid,                                       \
                                ('var', name, 'packed_data'), \
                                reasoning)
                        reasoning = []
                    elif type(scale) == type(offset) != type(np.reshape(ds.dataset.variables[name],data_size)[x]):
                        if type(scale) in [float(), type(np.float32(1.)), int(), type(np.int16(1))]  and (type(np.reshape(ds.dataset.variables[name],data_size)[x])) in [type(np.int8(1)), type(np.int16(1)), int()]:
                           valid = True
                           #reasoning.append("'add_offset' and 'scale_factor' are both of type float or int and the data variable is of type byte, short, or int.")
                           result = Result(BaseCheck.MEDIUM,                            \
                                valid,                                       \
                                ('var', name, 'packed_data'), \
                                reasoning)
                           ret_val.append(result)
                           reasoning = []
                        else: 
                           valid = False

                           reasoning.append("'add_offset' and 'scale_factor' are not of type float or int, or the data variable is not of type byte, short, or int.")
                           result = Result(BaseCheck.MEDIUM,                            \
                                valid,                                       \
                                ('var', name, 'packed_data'), \
                                reasoning)
                           ret_val.append(result) 
                           reasoning = []    
                    else:
                        valid = False
                        reasoning.append("'add_offset' and 'scale_factor' are not both of type float or int or the data variable is not of type byte, short, or int.")
                        result = Result(BaseCheck.MEDIUM,                            \
                                valid,                                       \
                                ('var', name, 'packed_data'), \
                                reasoning)
                        ret_val.append(result) 
                        reasoning = []



        return ret_val

        

    def check_compression(self, ds):
        """
        8.2 To save space in the netCDF file, it may be desirable to eliminate points from data arrays that are invariably
        missing. Such a compression can operate over one or more adjacent axes, and is accomplished with reference to a list
        of the points to be stored. 

        The list is stored as the coordinate variable for the compressed axis of the data array. Thus, the list variable and
        its dimension have the same name. The list variable has a string attribute compress, containing a blank-separated
        list of the dimensions which were affected by the compression in the order of the CDL declaration of the uncompressed
        array. 
        """
        ret_val = []
        reasoning = []

        for name, var in ds.dataset.variables.iteritems():
            if getattr(var, 'compress', ''):
                valid_form = False
                valid_dim = False

                if len(getattr(var, 'compress', '').split(" ")) >= 1:
                    valid_form = True
                    reasoning.append("The 'compress' attribute is in the form of a coordinate.")
                else: 
                    reasoning.append("The 'compress' attribute is not in the form of a coordinate.")

                for name_again, var_again in ds.dataset.variables.iteritems():
                
                    if name in ds.dataset.variables[name].dimensions:
                        valid_dim = True
                        reasoning.append("The 'compress' attribute is a referenced dimension.")
                    else:
                        reasoning.append("The 'compress' attribute is not a referenced dimension.")

                result = Result(BaseCheck.MEDIUM,                            \
                                valid_form and valid_dim,                                       \
                                ('var', name, 'compressed_data'), \
                                reasoning)
                ret_val.append(result)

        return ret_val
    ###############################################################################
    #
    # CHAPTER 9: Discrete Sampling Geometries
    #
    ###############################################################################

    @is_likely_dsg
    def check_all_features_are_same_type(self, ds):
        """
        9.1 The features contained within a collection must always be of the same type; and all the collections in a CF file
        must be of the same feature type. 

        point, timeSeries, trajectory, profile, timeSeriesProfile, trajectoryProfile.

        The space-time coordinates that are indicated for each feature are mandatory.  However a featureType may also include
        other space-time coordinates which are not mandatory (notably the z coordinate).
        """
        flag = 0
        x = ''
        y = ''
        z = ''
        t = ''

        flag = 0    
        for name,var in ds.dataset.variables.iteritems():
            if getattr(var,"grid_mapping_name", ""):
                #DO GRIDMAPPING CHECKS FOR X,Y,Z,T
                flag = 1
                for name_again, var_again in ds.dataset.variables.iteritems():
                    if getattr(var_again,"standard_name","") == self.grid_mapping_dict[getattr(var,"grid_mapping_name", "")][2][0]:
                        x = name_again
                    if getattr(var_again,"standard_name","") == self.grid_mapping_dict[getattr(var,"grid_mapping_name", "")][2][1]:
                        y = name_again
        
                        
                    
                
 
        for name,var in ds.dataset.variables.iteritems():
            #DO STANDARD SEARCH
            if getattr(var,'units','').lower() in ['pa', 'kpa', 'mbar', 'bar', 'atm', 'hpa', 'dbar'] or getattr(var,'positive','') or getattr(var,'standard_name','') == 'z' or getattr(var,'axis','') ==  'z':
                z = name
            if name.lower() in ['lon', 'longitude'] and flag == 0:
                x = name
            elif name.lower()in ['lat', 'latitude'] and flag == 0:
                y = name
            elif name.lower() == 'time':
                t = name
                    
            if getattr(var, '_CoordinateAxisType', ''):
                axis_type = getattr(var, '_CoordinateAxisType', '')
                if axis_type.lower() in ['lon', 'longitude'] and flag == 0:
                    x = name
                elif axis_type.lower()in ['lat', 'latitude'] and flag == 0:
                    y = name
                elif axis_type.lower() == 'time':
                    t = name
                                      
        valid = False                
        feature_types = []


        #create shape size tuple
        if x =='' or y == '' or t == '':
            return
        elif z == '':
            feature_tuple = (len(ds.dataset.variables[x].shape), len(ds.dataset.variables[y].shape), len(ds.dataset.variables[t].shape))
        else:
            feature_tuple = (len(ds.dataset.variables[x].shape), len(ds.dataset.variables[y].shape), len(ds.dataset.variables[t].shape), len(ds.dataset.variables[z].shape))



        #point
        if feature_tuple == (0,0,0):
            for name,var in ds.dataset.variables.iteritems():
                if name not in [x,y,t] and ds.dataset.variables[name].shape == ds.dataset.variables[t].shape:
                    feature_types.append('point')
                    
        #timeSeries
        if feature_tuple == (0,0,1) or feature_tuple == (1,1,2) :
            for name,var in ds.dataset.variables.iteritems():
                if name not in [x,y,t] and ds.dataset.variables[name].shape == ds.dataset.variables[t].shape:
                    feature_types.append('timeSeries')

        #trajectory or profile
        if feature_tuple == (1,1,1) :
            for name,var in ds.dataset.variables.iteritems():
                if name not in [x,y,t] and ds.dataset.variables[name].shape == ds.dataset.variables[t].shape:
                    feature_types.append('point-or-trajectory')

        #trajectory
        if feature_tuple == (2,2,2) :
            for name,var in ds.dataset.variables.iteritems():
                if name not in [x,y,t] and ds.dataset.variables[name].shape == ds.dataset.variables[t].shape:
                    feature_types.append('trajectory')
                    
        #profile
        if feature_tuple == (0,0,1,0) :
            for name,var in ds.dataset.variables.iteritems():
                if name not in [x,y,t,z] and ds.dataset.variables[name].shape == ds.dataset.variables[z].shape:
                    feature_types.append('profile')
        
        #timeseriesProfile
        if feature_tuple == (0,0,2,1) or feature_tuple == (1,1,3,2) :
            for name,var in ds.dataset.variables.iteritems():
                if name not in [x,y,t,z] and ds.dataset.variables[name].shape == ds.dataset.variables[z].shape:
                   feature_types.append('timeSeriesProfile')  
              
        #profile or trajectoryProfile
        if feature_tuple == (1,1,2,1) :
            for name,var in ds.dataset.variables.iteritems():
                if name not in [x,y,t,z] and ds.dataset.variables[name].shape == ds.dataset.variables[z].shape:
                   feature_types.append('profile-or-trajectoryProfile')

        #trajectoryProfile
        if  feature_tuple == (2,2,3,2) :
            for name,var in ds.dataset.variables.iteritems():
                if name not in [x,y,t,z] and ds.dataset.variables[name].shape == ds.dataset.variables[z].shape:
                   feature_types.append('trajectoryProfile')  

        valid = all(x == feature_types[0] for x in feature_types)

        return Result(BaseCheck.HIGH, valid)

    @is_likely_dsg
    def check_orthogonal_multidim_array(self, ds):
        """
        9.3.1 The orthogonal multidimensional array representation, the simplest representation, can be used if each feature
        instance in the collection has identical coordinates along the element axis of the features. 

        Data variables have both an instance dimension and an element dimension.  The dimensions may be given in any order. 
        If there is a need for either the instance or an element dimension to be the netCDF unlimited dimension (so that more
        features or more elements can be appended), then that dimension must be the outer dimension of the data variable
        i.e. the leading dimension in CDL.
        """
        ret_val = []
        reasoning = []

        for name,var in ds.dataset.variables.iteritems():
            if not getattr(var,'count_variable','') and not getattr(var,'index_variable',''):
                if getattr(var, '_FillValue', ''):
                    if getattr(var, '_FillValue', '') in ds.dataset.variables[name][:]:
                        pass
                    else:
                        reasoning.append('No _FillValue was used and there is no count_variable or index_variable field, so therefore it must be an orthogonal multidimensional array.')
                        result = Result(BaseCheck.MEDIUM,                            \
                                True,                                       \
                                ('var', name, 'orthogonal_multidimensional'), \
                                reasoning)
                        ret_val.append(result)
        return ret_val
        
    @is_likely_dsg
    def check_incomplete_multidim_array(self, ds):
        """
        9.3.2 The incomplete multidimensional array representation can used if the features within a collection do not all have
        the same number of elements, but sufficient storage space is available to allocate the number of elements required by
        the longest feature to all features.  That is, features that are shorter than the longest feature must be padded with
        missing values to bring all instances to the same storage size.

        Data variables have both an instance dimension and an element dimension.  The dimensions may be given in any order. 
        If there is a need for either the instance or an element dimension to be the netCDF unlimited dimension (so that more
        features or more elements can be appended), thlen that dimension must be the outer dimension of the data variable
        i.e. the leading dimension in CDL.
        """

        ret_val = []
        reasoning = []
        for name,var in ds.dataset.variables.iteritems():
            if not getattr(var,'count_variable','') and not getattr(var,'index_variable',''):
                if getattr(var, '_FillValue', ''):
                    if getattr(var, '_FillValue', '') in ds.dataset.variables[name][:]:
                        reasoning.append('A _FillValue was given and used and there is no count_variable or index_variable field, so therefore it must be an orthogonal multidimensional array.')
                        result = Result(BaseCheck.MEDIUM,                            \
                                True,                                       \
                                ('var', name, 'ragged_multidimensional'), \
                                reasoning)
                        ret_val.append(result)
                    else:
                        return []

        return ret_val

    @is_likely_dsg
    def check_contiguous_ragged_array(self, ds):
        """
        9.3.3 The contiguous ragged array representation can be used only if the size of each feature is known at the time
        that it is created. 

        In this representation, the file contains a count variable, which must be of type integer and must have the instance
        dimension as its sole dimension.  The count variable contains the number of elements that each feature has. This
        representation and its count variable are identifiable by the presence of an attribute, sample_dimension, found on
        the count variable, which names the sample dimension being counted. For indices that correspond to features, whose
        data have not yet been written, the count variable should  have a value of zero or a missing value. 

        In the ragged array representations, the instance dimension (i), which sequences the individual features within the
        collection, and the element dimension, which sequences the data elements of each feature (o and p), both occupy the
        same dimension (the sample dimension).   If the sample dimension is the netCDF unlimited dimension, new data can be
        appended to the file.  
        """
        ret_val = []
        reasoning = []
        for name,var in ds.dataset.variables.iteritems():
            if getattr(var,'count_variable',''):
                reasoning.append('The count_variable is present, therefore the data is a contiguous ragged array.')
                result = Result(BaseCheck.MEDIUM,                            \
                        True,                                       \
                        ('var', name, 'continuous_ragged'), \
                        reasoning)
                ret_val.append(result)
            else:
                return []

        return ret_val

    @is_likely_dsg
    def check_indexed_ragged_array(self, ds):
        """
        9.3.4 The indexed ragged array representation stores the features interleaved along the sample dimension in the data
        variable.

        In this representation, the file contains an index variable, which must be of type integer, and must have the sample
        dimension as its single dimension. The index variable contains the zero-based index of the feature to which each
        element belongs. This representation is identifiable by the presence of an attribute, instance_dimension, on the index
        variable, which names the dimension of the instance variables. For those indices of the sample dimension, into which
        data have not yet been written, the index variable should be pre-filled with missing values.

        In the ragged array representations, the instance dimension (i), which sequences the individual features within the
        collection, and the element dimension, which sequences the data elements of each feature (o and p), both occupy the
        same dimension (the sample dimension).   If the sample dimension is the netCDF unlimited dimension, new data can be
        appended to the file.  
        """
        ret_val = []
        reasoning = []
        for name,var in ds.dataset.variables.iteritems():
            if getattr(var,'index_variable',''):
                reasoning.append('The index_variable is present, therefore the data is an indexed ragged array.')
                result = Result(BaseCheck.MEDIUM,                            \
                        True,                                       \
                        ('var', name, 'continuous_ragged'), \
                        reasoning)
                ret_val.append(result)
            else:
                return []

        return ret_val

    @is_likely_dsg
    def check_feature_type(self, ds):
        """
        9.4 A global attribute, featureType, is required for all Discrete Geometry representations except the orthogonal
        multidimensional array representation, for which it is highly recommended.

        The value assigned to the featureType attribute is case-insensitive.
        """
        ret_val = []
        reasoning=[]
        feature_list = ['point', 'timeseries','trajectory','profile', 'timeseriesprofile','trajectoryprofile']

        if getattr(ds.dataset, 'featureType', '').lower() in feature_list:
            reasoning.append('The featureType is provided and is from the featureType list.')
            result = Result(BaseCheck.MEDIUM,                            \
                        True, 'feature_type',  
                        reasoning)
            ret_val.append(result)
        elif getattr(ds.dataset, 'featureType', ''):
            reasoning.append('The featureType is provided and is not from the featureType list.')
            result = Result(BaseCheck.MEDIUM, 'feature_type',                            \
                        False,
                        reasoning)
            ret_val.append(result)
        return ret_val

    @is_likely_dsg
    def check_coordinates_and_metadata(self, ds):
        """
        9.5 Every feature within a Discrete Geometry CF file must be unambiguously associated with an extensible collection
        of instance variables that identify the feature and provide other metadata as needed to describe it.  Every element of
        every feature must be unambiguously associated with its space and time coordinates and with the feature that contains
        it.

        The coordinates attribute must be attached to every data variable to indicate the spatiotemporal coordinate variables
        that are needed to geo-locate the data.


        Where feasible a variable with the attribute cf_role should be included.  The only acceptable values of cf_role for
        Discrete Geometry CF data sets are timeseries_id, profile_id, and trajectory_id.   The variable carrying the cf_role
        attribute may have any data type.  When a variable is assigned this attribute, it must provide a unique identifier
        for each feature instance.

        CF files that contain timeSeries, profile or trajectory featureTypes, should include only a single occurrence of a
        cf_role attribute;  CF files that contain timeSeriesProfile or trajectoryProfile may contain two occurrences,
        corresponding to the two levels of structure in these feature types.

        CF Discrete Geometries provides a mechanism to encode both the nominal and the precise positions, while retaining the
        semantics of the idealized feature type. Only the set of coordinates which are regarded as the nominal (default or
        preferred) positions should be indicated by the attribute axis, which should be assigned string values to indicate
        the orientations of the axes (X, Y, Z, or T).

        Auxiliary coordinate variables containing the nominal and the precise positions should be listed in the relevant
        coordinates attributes of data variables. In orthogonal representations the nominal positions could be  coordinate
        variables, which do not need to be listed in the coordinates attribute, rather than auxiliary coordinate variables.

        Coordinate bounds may optionally be associated with coordinate variables and auxiliary coordinate variables using
        the bounds attribute.

        If there is a vertical coordinate variable or auxiliary coordinate variable, it must be identified by the means
        specified in section 4.3.   The use of the attribute axis=Z is recommended for clarity.  A standard_name attribute
        that identifies the vertical coordinate is recommended.
        """
        ret_val = []
        reasoning = []

        name_list = []
        non_data_list = []
        data_list = []

        for name,var in ds.dataset.variables.iteritems():


            if var.dimensions and not hasattr(var, 'cf_role'):
                if var.dimensions != (name,):
                    name_list.append(name)

        for name,var in ds.dataset.variables.iteritems():
            if hasattr(var, 'coordinates'):
                for each in getattr(var, 'coordinates', '').split(' '):
                    if each in name_list:
                        non_data_list.append(each)

        data_list = [each for each in name_list if each not in non_data_list]


        for each in data_list:
            if getattr(ds.dataset.variables[each], 'coordinates', ''):
                reasoning.append('The data variable %s has associated coordinates' %each)
                result = Result(BaseCheck.MEDIUM,                            \
                            True,                                       \
                            ('var', each, 'check_coordinates'), \
                            reasoning)
                ret_val.append(result)
                reasoning = []
            else:
                reasoning.append('The variable %s does not have associated coordinates' %each)
                result = Result(BaseCheck.MEDIUM,                            \
                            False,                                       \
                            ('var', each, 'check_coordinates'), \
                            reasoning)
                ret_val.append(result)
                reasoning = []

        role_list = [getattr(var, 'cf_role', '').split(' ') for name,var in ds.dataset.variables.iteritems() if hasattr(var, 'cf_role')]
        single_role = ['timeseries', 'profile', 'trajectory']
        dual_role = ['timeseries', 'profile', 'trajectory','timeSeriesProfile', 'trajectoryProfile']
        if getattr(ds.dataset, 'featureType', '').lower() in single_role and len(np.ravel(role_list)) == 1:
            reasoning = []
            valid = True
        elif getattr(ds.dataset, 'featureType', '').lower() in dual_role and len(np.ravel(role_list)) in [1,2]:
            reasoning = []
            valid = True
        else:
            valid = False
            reasoning = []
            reasoning.append('The cf_role featureType is not properly defined.')

        result = Result(BaseCheck.MEDIUM,                            \
                valid,                                       \
                'check_cf_role', \
                reasoning)
        ret_val.append(result)

        return ret_val

    @is_likely_dsg
    def check_missing_data(self, ds):
        """
        9.6 Auxiliary coordinate variables (spatial and time) must contain missing values to indicate a void in data storage
        in the file but must not have missing data for any other reason.

        It is not permitted for auxiliary coordinate variables to have missing values for elements where there is non-missing
        data. Where any auxiliary coordinate variable contains a missing value, all other coordinate, auxiliary coordinate
        and data values corresponding to that element should also contain missing values. Data variables should (as usual)
        also contain missing values to indicate when there is no valid data available for the element, although the
        coordinates are valid.

        Similarly, for indices where the instance variable identified by cf_role contains a missing value indicator, all other
        instance variable should also contain missing values.
        """


        ret_val = []
        
        
        name_list = ds.dataset.variables.keys()
        dim_list = ds.dataset.dimensions.keys()

        aux_index_dict = {}
        dim_index_dict = {}


        for name, var in ds.dataset.variables.iteritems():
            reasoning = []
            valid = False
            aux_valid = False
            
            if hasattr(var, 'coordinates'):
                for coordinate in getattr(var, 'coordinates', '').split(" "):
                    indices = []
                    if coordinate in name_list and coordinate not in dim_list:
                        try:
                            indices = np.where(ds.dataset.variables[coordinate] == var._FillValue).tolist()
                        except:
                            indices = np.where(ds.dataset.variables[coordinate] == var._FillValue)[0].tolist()

                        dim_index_dict[name+'-'+coordinate] = indices
                        aux_index_dict[name+'-'+coordinate] = indices
                            
                    elif coordinate in name_list and coordinate in dim_list:
                        try:
                            indices = np.where(ds.dataset.variables[coordinate] == var._FillValue).tolist()
                        except:
                            indices = np.where(ds.dataset.variables[coordinate] == var._FillValue)[0].tolist()
                        dim_index_dict[name+'-'+coordinate] = indices
                    else:
                        dim_index_dict[name+'-'+coordinate] = []
    

            #Check to see that all coordinate variable mising data locations are the same
            aux_index_list = []
            for each in aux_index_dict:
                aux_index_list.append(aux_index_dict[each])
            if aux_index_list != []:  
                aux_valid = all(x == aux_index_list[0] for x in aux_index_list)
            else: 
                aux_valid = True
            
            #Check to see that all auxilliary coordinate variable missing data appears in the coordinate variables
            dim_index_list = []
            for each in dim_index_dict:
                dim_index_list.append(dim_index_dict[each])  
            if dim_index_list != []:
                valid = all(x == dim_index_list[0] for x in dim_index_list)
            else:
                valid = True

                
            if aux_valid == False:
                reasoning.append('The auxillary coordinates do not have the same missing data locations')
            if valid == False:
                reasoning.append('The coordinate variables do not have the same missing data locations as the auxillary coordinates')
            
            #Check to see that all coordinate variable mising data is reflceted in the dataset
            valid_missing = True
            count = 0

            if hasattr(var, '_FillValue'):
                try:
                    x_indices = np.where(var==var._FillValue).tolist()
                except:
                    x_indices = np.where(var==var._FillValue)[0].tolist()

                if hasattr(var, 'coordinates'):
                    for coordinate in var.coordinates.split(" "):
                        coordinate_ind_list = dim_index_dict[name+'-'+coordinate]
                        valid_missing = all(each in x_indices[count] for each in coordinate_ind_list)
                        count = count+1

                if valid_missing == False:
                    reasoning.append('The data does not have the same missing data locations as the coordinates')
            

            result = Result(BaseCheck.MEDIUM,                            \
                valid and aux_valid and valid_missing,                                       \
                ('var', name, 'missing_data'), \
                reasoning)
            ret_val.append(result)
        return ret_val

class CFNCCheck(BaseNCCheck, CFBaseCheck):
    @classmethod
    def beliefs(cls): # @TODO
        return {}

