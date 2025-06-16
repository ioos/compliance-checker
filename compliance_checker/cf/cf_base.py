#!/usr/bin/env python
import logging
import os
import sys
from collections import OrderedDict, defaultdict
from warnings import warn

import numpy as np
import regex

import compliance_checker.cf.util as cfutil
from compliance_checker.base import BaseCheck, BaseNCCheck, Result, TestCtx
from compliance_checker.cf import util
from compliance_checker.cf.appendix_d import no_missing_terms

logger = logging.getLogger(__name__)


class CFBaseCheck(BaseCheck):
    """
    CF Convention Checker Base
    """

    def __init__(self, options=None):
        # The compliance checker can be run on multiple datasets in a single
        # instantiation, so caching values has be done by the unique identifier
        # for each dataset loaded.

        # Each default dict is a key, value mapping from the dataset object to
        # a list of variables
        super().__init__(options)
        self._coord_vars = defaultdict(list)
        self._ancillary_vars = defaultdict(list)
        self._clim_vars = defaultdict(list)
        self._metadata_vars = defaultdict(list)
        self._boundary_vars = defaultdict(list)
        self._geophysical_vars = defaultdict(list)
        self._aux_coords = defaultdict(list)

        self._std_names = util.StandardNameTable()

        self.section_titles = {  # dict of section headers shared by grouped checks
            "1.2": "§1.2 Terminology",
            "2.1": "§2.1 Filename",
            "2.2": "§2.2 Data Types",
            "2.3": "§2.3 Naming Conventions",
            "2.4": "§2.4 Dimensions",
            "2.5": "§2.5 Variables",
            "2.5.1": "§2.5.1. Missing data, valid and actual range of data",
            "2.6": "§2.6 Attributes",
            "2.6.3": "§2.6.3 External Variables",
            "2.7.1": "§2.7.1. Scope",
            "3.1": "§3.1 Units",
            "3.2": "§3.2 Long Name",
            "3.3": "§3.3 Standard Name",
            "3.4": "§3.4 Ancillary Data",
            "3.5": "§3.5 Flags",
            "4": "§4 Coordinate Types",
            "4.1": "§4.1 Latitude Coordinate",
            "4.2": "§4.2 Longitude Coordinate",
            "4.3": "§4.3 Vertical Coordinate",
            "4.4": "§4.4 Time Coordinate",
            "4.4.1": "§4.4.1 Calendar",
            "4.5": "§4.5 Discrete Axis",
            "5": "§5 Coordinate Systems",
            "5.1": "§5.1 Independent Latitude, Longitude, Vertical, and Time Axes",
            "5.2": "§5.2 2-D Latitude, Longitude, Coordinate Variables",
            "5.3": "§5.3 Reduced Horizontal Grid",
            "5.4": "§5.4 Timeseries of Station Data",
            "5.5": "§5.5 Trajectories",
            "5.6": "§5.6 Horizontal Coordinate Reference Systems, Grid Mappings, Projections",
            "5.7": "§5.7 Scalar Coordinate Variables",
            "6.1": "§6.1 Labels",
            "6.2": "§6.2 Alternative Coordinates",
            "7.1": "§7.1 Cell Boundaries",
            "7.2": "§7.2 Cell Measures",
            "7.3": "§7.3 Cell Methods",
            "7.4": "§7.4 Climatological Statistics",
            "8.1": "§8.1 Packed Data",
            "8.2": "§8.2 Compression by Gathering",
            "9.1": "§9.1 Features and feature types",
            "9.2": "§9.2 Collections, instances, and elements",
            "9.3": "§9.3 Representations of Collections of features in data variables",
            "9.4": "§9.4 The featureType attribute",
            "9.5": "§9.5 Coordinates and metadata",
            "9.6": "§9.6 Missing Data",
        }

    ################################################################################
    # Helper Methods - var classifications, etc
    ################################################################################

    def setup(self, ds):
        """
        Initialize various special variable types within the class.
        Mutates a number of instance variables.

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        self.coord_vars = self._find_coord_vars(ds)
        self._find_aux_coord_vars(ds)
        self._find_ancillary_vars(ds)
        self._find_clim_vars(ds)
        self._find_boundary_vars(ds)
        self._find_metadata_vars(ds)
        self._find_cf_standard_name_table(ds)
        self._find_geophysical_vars(ds)
        coord_containing_vars = ds.get_variables_by_attributes(
            coordinates=lambda val: isinstance(val, str),
        )

        # coordinate data variables

        # Excerpt from "§1.3 Overview" on coordinate data
        # There are two methods used to identify variables that contain
        # coordinate data. The first is to use the NUG-defined "coordinate
        # variables." The use of coordinate variables is required for all
        # dimensions that correspond to one dimensional space or time
        # coordinates . In cases where coordinate variables are not applicable,
        # the variables containing coordinate data are identified by the
        # coordinates attribute.

        # first read in variables referred to in coordinates which exist
        # in the dataset
        self.coord_data_vars = set()
        for var in coord_containing_vars:
            for coord_var_name in var.coordinates.strip().split(" "):
                if coord_var_name in ds.variables:
                    self.coord_data_vars.add(coord_var_name)
        # then add in the NUG coordinate variables -- single dimension with
        # dimension name the same as coordinates
        self.coord_data_vars.update(self.coord_vars)

    def check_grid_mapping(self, ds):
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

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """

        ret_val = OrderedDict()
        grid_mapping_variables = cfutil.get_grid_mapping_variables(ds)

        # Check the grid_mapping attribute to be a non-empty string and that its reference exists
        for variable in ds.get_variables_by_attributes(
            grid_mapping=lambda x: x is not None,
        ):
            grid_mapping = getattr(variable, "grid_mapping", None)
            defines_grid_mapping = self.get_test_ctx(
                BaseCheck.HIGH,
                self.section_titles["5.6"],
                variable.name,
            )
            defines_grid_mapping.assert_true(
                (isinstance(grid_mapping, str) and grid_mapping),
                f"{variable.name}'s grid_mapping attribute must be a "
                "space-separated non-empty string",
            )
            if isinstance(grid_mapping, str):
                # TODO (badams): refactor functionality to split functionality
                #                into requisite classes
                if ":" in grid_mapping and self._cc_spec_version >= "1.7":
                    colon_count = grid_mapping.count(":")
                    re_all = regex.findall(
                        r"(\w+):\s*((?:\w+\s+)*(?:\w+)(?![\w:]))",
                        grid_mapping,
                    )
                    if colon_count != len(re_all):
                        defines_grid_mapping.out_of += 1
                        defines_grid_mapping.messages.append(
                            "Could not consume entire grid_mapping expression, please check for well-formedness",
                        )
                    else:
                        for grid_var_name, coord_var_str in re_all:
                            defines_grid_mapping.assert_true(
                                grid_var_name in ds.variables,
                                f"grid mapping variable {grid_var_name} must exist in this dataset",
                            )
                            for ref_var in coord_var_str.split():
                                defines_grid_mapping.assert_true(
                                    ref_var in ds.variables,
                                    f"Coordinate-related variable {ref_var} referenced by grid_mapping variable {grid_var_name} must exist in this dataset",
                                )

                else:
                    for grid_var_name in grid_mapping.split():
                        defines_grid_mapping.assert_true(
                            grid_var_name in ds.variables,
                            f"grid mapping variable {grid_var_name} must exist in this dataset",
                        )
            ret_val[variable.name] = defines_grid_mapping.to_result()

        # Check the grid mapping variables themselves
        for grid_var_name in grid_mapping_variables:
            valid_grid_mapping = self.get_test_ctx(
                BaseCheck.HIGH,
                self.section_titles["5.6"],
                grid_var_name,
            )
            grid_var = ds.variables[grid_var_name]

            grid_mapping_name = getattr(grid_var, "grid_mapping_name", None)

            # Grid mapping name must be in appendix F
            valid_grid_mapping.assert_true(
                grid_mapping_name in self.grid_mapping_dict,
                f"{grid_mapping_name} is not a valid grid_mapping_name."
                + " See Appendix F for valid grid mappings",
            )

            # The self.grid_mapping_dict has a values of:
            # - required attributes
            # - optional attributes (can't check)
            # - required standard_names defined
            # - at least one of these attributes must be defined

            # We can't do any of the other grid mapping checks if it's not a valid grid mapping name
            if grid_mapping_name not in self.grid_mapping_dict:
                ret_val[grid_mapping_name] = valid_grid_mapping.to_result()
                continue

            grid_mapping = self.grid_mapping_dict[grid_mapping_name]
            required_attrs = grid_mapping[0]
            # Make sure all the required attributes are defined
            for req in required_attrs:
                valid_grid_mapping.assert_true(
                    hasattr(grid_var, req),
                    f"{req} is a required attribute for grid mapping {grid_mapping_name}",
                )

            # Make sure that exactly one of the exclusive attributes exist
            if len(grid_mapping) == 4:
                at_least_attr = grid_mapping[3]
                number_found = 0
                for attr in at_least_attr:
                    if hasattr(grid_var, attr):
                        number_found += 1
                valid_grid_mapping.assert_true(
                    number_found == 1,
                    f"grid mapping {grid_mapping_name}"
                    + "must define exactly one of these attributes: "
                    + "{}".format(" or ".join(at_least_attr)),
                )

            # Make sure that exactly one variable is defined for each of the required standard_names
            expected_std_names = grid_mapping[2]
            for expected_std_name in expected_std_names:
                found_vars = ds.get_variables_by_attributes(
                    standard_name=expected_std_name,
                )
                valid_grid_mapping.assert_true(
                    len(found_vars) == 1,
                    f"grid mapping {grid_mapping_name} requires exactly "
                    + "one variable with standard_name "
                    + f"{expected_std_name} to be defined",
                )

            ret_val[grid_var_name] = valid_grid_mapping.to_result()

        return ret_val

    def check_conventions_version(self, ds):
        """
        CF §2.6.1 the NUG defined global attribute Conventions to the string
        value "CF-<version_number>"; check the Conventions attribute contains
        the appropriate string.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        """

        valid = False
        reasoning = []
        correct_version_string = f"{self._cc_spec}-{self._cc_spec_version}".upper()
        if hasattr(ds, "Conventions"):
            conventions = regex.split(r",|\s+", getattr(ds, "Conventions", ""))
            for convention in conventions:
                if convention == correct_version_string:
                    valid = True
                    break
            else:
                reasoning = [
                    "§2.6.1 Conventions global attribute does not contain "
                    f'"{correct_version_string}"',
                ]
        else:
            valid = False
            reasoning = ["§2.6.1 Conventions field is not present"]
        return Result(
            BaseCheck.MEDIUM,
            valid,
            self.section_titles["2.6"],
            msgs=reasoning,
        )

    def _check_dimensionless_vertical_coordinates(
        self,
        ds,
        deprecated_units,
        version_specific_check,
        version_specific_dimless_vertical_coord_dict,
    ):
        """
        Check the validity of dimensionless coordinates under CF

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param list deprecated_units: list of string names of deprecated units
        :param function version_specific_check: version-specific implementation to check dimensionless vertical coord
        :param dict version_specific_dimless_coord_dict: version-specific dict of dimensionless vertical coords and computed standard names
        :return: List of results
        """
        ret_val = []

        z_variables = cfutil.get_z_variables(ds)

        # call version-specific implementation
        for name in z_variables:
            version_specific_check(
                ds,
                name,
                deprecated_units,
                ret_val,
                version_specific_dimless_vertical_coord_dict,
            )

        return ret_val

    def _check_formula_terms(self, ds, coord, dimless_coords_dict):
        """
        Checks a dimensionless vertical coordinate contains valid formula_terms

        - formula_terms is a non-empty string
        - formula_terms matches regdimless_coords_dictx
        - every variable defined in formula_terms exists

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        """
        variable = ds.variables[coord]
        standard_name = getattr(variable, "standard_name", None)
        formula_terms = getattr(variable, "formula_terms", None)
        valid_formula_terms = TestCtx(BaseCheck.HIGH, self.section_titles["4.3"])

        valid_formula_terms.assert_true(
            isinstance(formula_terms, str) and formula_terms,
            f"§4.3.2: {coord}'s formula_terms is a required attribute and must be a non-empty string"
            "",
        )
        # We can't check any more
        if not formula_terms:
            return valid_formula_terms.to_result()

        # check that the formula_terms are well formed and are present
        # The pattern for formula terms is always component: variable_name
        # the regex grouping always has component names in even positions and
        # the corresponding variable name in odd positions.
        poorly_formed_formula_terms = ("Attribute formula_terms is not well-formed",)
        matches = list(
            regex.finditer(
                r"(\w+):\s+(\w+)(?:\s+(?!$)|$)",
                variable.formula_terms,
            ),
        )
        if not matches:
            valid_formula_terms.add_failure(poorly_formed_formula_terms)
            return valid_formula_terms.to_result()

        terms = {m.group(1) for m in matches}
        # get the variables named in the formula terms and check if any
        # are not present in the dataset
        missing_vars = sorted({m.group(2) for m in matches} - set(ds.variables))
        missing_fmt = "The following variable(s) referenced in {}:formula_terms are not present in the dataset: {}"
        valid_formula_terms.assert_true(
            len(missing_vars) == 0,
            missing_fmt.format(coord, ", ".join(missing_vars)),
        )
        # try to reconstruct formula_terms by adding space in between the regex
        # matches.  If it doesn't exactly match the original, the formatting
        # of the attribute is incorrect
        reconstructed_formula = "".join(m.group(0) for m in matches)
        valid_formula_terms.assert_true(
            reconstructed_formula == formula_terms,
            "Attribute formula_terms is not well-formed",
        )

        valid_formula_terms.assert_true(
            standard_name in dimless_coords_dict,
            f"unknown standard_name '{standard_name}' for dimensionless vertical coordinate {coord}"
            "",
        )
        if standard_name not in dimless_coords_dict:
            return valid_formula_terms.to_result()

        valid_formula_terms.assert_true(
            no_missing_terms(standard_name, terms, dimless_coords_dict),
            f"{coord}'s formula_terms are invalid for {standard_name}, please see appendix D of CF 1.6"
            "",
        )

        return valid_formula_terms.to_result()

    def _check_grid_mapping_attr_condition(self, attr, attr_name, ret_val):
        """
        Evaluate a condition (or series of conditions) for a particular
        attribute. Designed to be overloaded in subclass implementations.

        :param attr: attribute to teset condition for
        :param str attr_name: name of the attribute
        :param list ret_val: list of results to append to
        :rtype None
        :return None
        """
        raise NotImplementedError

    def _dims_in_order(self, dimension_order):
        """
        :param list dimension_order: A list of axes
        :rtype: bool
        :return: Returns True if the dimensions are in order U*, T, Z, Y, X,
                 False otherwise
        """
        regx = regex.compile(r"^[^TZYX]*T?Z?Y?X?$")
        dimension_string = "".join(dimension_order)
        return regx.match(dimension_string) is not None

    def _parent_var_attr_type_check(self, attr_name, var, ctx):
        """
        Checks that an attribute has an equivalent value to a parent variable.
        Takes an attribute name, variable, and test context on which to operate.
        :param str attr_name: The name of the attribute to be checked
        :param netCDF4.Variable var: The variable against which to be checked
        :param compliance_checker.base.TestCtx ctx: The associated test context to modify
        :rtype None
        :return None
        """
        attr_val = var.getncattr(attr_name)

        if isinstance(attr_val, (str, bytes)):
            type_match = (var.dtype is str) or (var.dtype.kind == "S")
            val_type = type(attr_val)
        else:
            val_type = attr_val.dtype.type
            type_match = val_type == var.dtype.type

        ctx.assert_true(
            type_match,
            f"Attribute '{attr_name}' (type: {val_type}) and parent variable '{var.name}' (type: {var.dtype.type}) "
            "must have equivalent datatypes",
        )

    def _find_aux_coord_vars(self, ds, refresh=False):
        """
        Returns a list of auxiliary coordinate variables

        An auxiliary coordinate variable is any netCDF variable that contains
        coordinate data, but is not a coordinate variable (in the sense of the term
        defined by CF).

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param bool refresh: if refresh is set to True, the cache is
                             invalidated.
        :rtype: list
        :return: List of variable names (str) that are defined to be auxiliary
                 coordinate variables.
        """
        ds_str = ds.__str__()
        if self._aux_coords.get(ds, None) and refresh is False:
            return self._aux_coords[ds_str]

        self._aux_coords[ds_str] = cfutil.get_auxiliary_coordinate_variables(ds)
        return self._aux_coords[ds_str]

    def _find_boundary_vars(self, ds, refresh=False):
        """
        Returns dictionary of boundary variables mapping the variable instance
        to the name of the variable acting as a boundary variable.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param bool refresh: if refresh is set to True, the cache is
                             invalidated.
        :rtype: list
        :return: A list containing strings with boundary variable names.
        """
        ds_str = ds.__str__()
        if self._boundary_vars.get(ds, None) and refresh is False:
            return self._boundary_vars[ds_str]

        self._boundary_vars[ds_str] = cfutil.get_cell_boundary_variables(ds)

        return self._boundary_vars[ds_str]

    def _find_ancillary_vars(self, ds, refresh=False):
        """
        Returns a list of variable names that are defined as ancillary
        variables in the dataset ds.

        An ancillary variable generally is a metadata container and referenced
        from other variables via a string reference in an attribute.

        - via ancillary_variables (3.4)
        - "grid mapping var" (5.6)
        - TODO: more?

        The result is cached by the passed in dataset object inside of this
        checker. Pass refresh=True to redo the cached value.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param bool refresh: if refresh is set to True, the cache is
                             invalidated.
        :rtype: list
        :return: List of variable names (str) that are defined as ancillary
                 variables in the dataset ds.
        """
        ds_str = ds.__str__()
        # Used the cached version if it exists and is not empty
        if self._ancillary_vars.get(ds, None) and refresh is False:
            return self._ancillary_vars[ds_str]

        # Invalidate the cache at all costs
        self._ancillary_vars[ds_str] = []

        for var in ds.variables.values():
            if hasattr(var, "ancillary_variables"):
                for anc_name in var.ancillary_variables.split(" "):
                    if anc_name in ds.variables:
                        self._ancillary_vars[ds_str].append(anc_name)

            if hasattr(var, "grid_mapping"):
                gm_name = var.grid_mapping
                if gm_name in ds.variables:
                    self._ancillary_vars[ds_str].append(gm_name)

        return self._ancillary_vars[ds_str]

    def _find_clim_vars(self, ds, refresh=False):
        """
        Returns a list of variables that are likely to be climatology variables based on CF §7.4

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param bool refresh: if refresh is set to True, the cache is
                             invalidated.
        :rtype: list
        :return: A list containing strings with geophysical variable
                 names.
        """
        ds_str = ds.__str__()
        if self._clim_vars.get(ds, None) and refresh is False:
            return self._clim_vars[ds_str]

        climatology_variable = cfutil.get_climatology_variable(ds)
        if climatology_variable:
            self._clim_vars[ds_str].append(climatology_variable)

        return self._clim_vars[ds_str]

    def _find_cf_standard_name_table(self, ds):
        """
        Parse out the `standard_name_vocabulary` attribute and download that
        version of the cf standard name table.  If the standard name table has
        already been downloaded, use the cached version.  Modifies `_std_names`
        attribute to store standard names.  Returns True if the file exists and
        False if it fails to download.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: bool
        """
        # Get the standard name vocab
        standard_name_vocabulary = getattr(ds, "standard_name_vocabulary", "")

        # Try to parse this attribute to get version
        version = None
        try:
            if "cf standard name table" in standard_name_vocabulary.lower():
                version = [
                    s.strip("(").strip(")").strip("v").strip(",")
                    for s in standard_name_vocabulary.split()
                ]
                # This assumes that table version number won't start with 0.
                version = [
                    s
                    for s in version
                    if s.isdigit() and len(s) <= 2 and not s.startswith("0")
                ]
                if len(version) > 1:
                    return False
                else:
                    try:
                        version = version[0]
                    except IndexError:
                        warn(
                            "Cannot extract CF standard name version number "
                            "from standard_name_vocabulary string",
                            stacklevel=2,
                        )
                        return False
            else:
                # Can't parse the attribute, use the packaged version
                return False
        # usually raised from .lower() with an incompatible (non-string)
        # data type
        except AttributeError:
            warn(
                "Cannot convert standard name table to lowercase.  This can "
                "occur if a non-string standard_name_vocabulary global "
                "attribute is supplied",
                stacklevel=2,
            )
            return False

        if version.startswith("v"):  # i.e 'v34' -> '34' drop the v
            version = version[1:]

        # If the packaged version is what we're after, then we're good
        if version == self._std_names._version:
            print(
                f"Using packaged standard name table v{version}",
                file=sys.stderr,
            )
            return False

        # Try to download the version specified
        try:
            data_directory = util.create_cached_data_dir()
            location = os.path.join(
                data_directory,
                f"cf-standard-name-table-test-{version}.xml",
            )
            # Did we already download this before?
            if not os.path.isfile(location):
                util.download_cf_standard_name_table(version, location)
                print(
                    f"Using downloaded standard name table v{version}",
                    file=sys.stderr,
                )
            else:
                print(
                    f"Using cached standard name table v{version} from {location}",
                    file=sys.stderr,
                )

            self._std_names = util.StandardNameTable(location)
            return True
        except Exception as e:
            # There was an error downloading the CF table. That's ok, we'll just use the packaged version
            warn(
                f"Problem fetching standard name table:\n{e}\n"
                f"Using packaged v{self._std_names._version}",
                stacklevel=2,
            )
            return False

    def _find_coord_vars(self, ds, refresh=False):
        """
        Returns a list of variable names that identify as coordinate variables.

        The result is cached by the passed in dataset object inside of this
        checker. Pass refresh=True to redo the cached value.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param bool refresh: if refresh is set to True, the cache is
                             invalidated.
        :rtype: list
        :return: A list of variables names (str) that are defined as coordinate
                 variables in the dataset ds.
        """
        ds_str = ds.__str__()
        if ds_str in self._coord_vars and refresh is False:
            return self._coord_vars[ds_str]

        self._coord_vars[ds_str] = cfutil.get_coordinate_variables(ds)

        return self._coord_vars[ds_str]

    def _find_geophysical_vars(self, ds, refresh=False):
        """
        Returns a list of geophysical variables.  Modifies
        `self._geophysical_vars`

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param bool refresh: if refresh is set to True, the cache is
                             invalidated.
        :rtype: list
        :return: A list containing strings with geophysical variable
                 names.
        """
        ds_str = ds.__str__()
        if self._geophysical_vars.get(ds, None) and refresh is False:
            return self._geophysical_vars[ds_str]

        self._geophysical_vars[ds_str] = cfutil.get_geophysical_variables(ds)

        return self._geophysical_vars[ds_str]

    def _find_metadata_vars(self, ds, refresh=False):
        """
        Returns a list of netCDF variable instances for those that are likely metadata variables

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param bool refresh: if refresh is set to True, the cache is
                             invalidated.
        :rtype: list
        :return:   List of variable names (str) that are likely metadata
                   variable candidates.

        """
        ds_str = ds.__str__()
        if self._metadata_vars.get(ds, None) and refresh is False:
            return self._metadata_vars[ds_str]

        self._metadata_vars[ds_str] = []
        for name, var in ds.variables.items():
            if name in self._find_ancillary_vars(ds) or name in self._find_coord_vars(
                ds,
            ):
                continue

            if name in (
                "platform_name",
                "station_name",
                "instrument_name",
                "station_id",
                "platform_id",
                "surface_altitude",
            ):
                self._metadata_vars[ds_str].append(name)

            elif getattr(var, "cf_role", "") != "":
                self._metadata_vars[ds_str].append(name)

            elif (
                getattr(var, "standard_name", None) is None and len(var.dimensions) == 0
            ):
                self._metadata_vars[ds_str].append(name)

        return self._metadata_vars[ds_str]

    def _get_coord_axis_map(self, ds):
        """
        Returns a dictionary mapping each coordinate to a letter identifier
        describing the _kind_ of coordinate.

        :param netCDF4.Dataset ds: An open netCDF dataset

        :rtype: dict
        :return: A dictionary with variable names mapped to axis abbreviations,
                 i.e. {'longitude': 'X', ... 'pressure': 'Z'}
        """
        expected = ["T", "Z", "Y", "X"]
        coord_vars = self._find_coord_vars(ds)
        coord_axis_map = {}

        # L - Unlimited Coordinates
        # T - Time coordinates
        # Z - Depth/Altitude Coordinate
        # Y - Y-Coordinate (latitude)
        # X - X-Coordinate (longitude)
        # A - Auxiliary Coordinate
        # I - Instance Coordinate

        time_variables = cfutil.get_time_variables(ds)
        lat_variables = cfutil.get_latitude_variables(ds)
        lon_variables = cfutil.get_longitude_variables(ds)
        z_variables = cfutil.get_z_variables(ds)

        for coord_name in coord_vars:
            coord_var = ds.variables[coord_name]
            axis = getattr(coord_var, "axis", None)
            standard_name = getattr(coord_var, "standard_name", None)

            # Unlimited dimensions must come first
            if ds.dimensions[coord_name].isunlimited():
                coord_axis_map[coord_name] = "L"
            # axis takes precedence over standard_name
            elif axis in expected:
                coord_axis_map[coord_name] = axis
            elif standard_name == "time":
                coord_axis_map[coord_name] = "T"
            elif standard_name == "longitude":
                coord_axis_map[coord_name] = "X"
            elif standard_name == "latitude":
                coord_axis_map[coord_name] = "Y"
            elif standard_name in ["height", "depth", "altitude"]:
                coord_axis_map[coord_name] = "Z"
            elif cfutil.is_compression_coordinate(ds, coord_name):
                coord_axis_map[coord_name] = "C"
            elif coord_name in time_variables:
                coord_axis_map[coord_name] = "T"
            elif coord_name in z_variables:
                coord_axis_map[coord_name] = "Z"
            elif coord_name in lat_variables:
                coord_axis_map[coord_name] = "Y"
            elif coord_name in lon_variables:
                coord_axis_map[coord_name] = "X"
            else:
                # mark the coordinate variable as unknown
                coord_axis_map[coord_name] = "U"

        for dimension in self._get_instance_dimensions(ds):
            if dimension not in coord_axis_map:
                coord_axis_map[dimension] = "I"

        # Dimensions of auxiliary coordinate variables will be marked with A.
        # This is useful to help determine if the dimensions are used like a
        # mapping from grid coordinates to physical lat/lon
        for coord_name in self._find_aux_coord_vars(ds):
            coord_var = ds.variables[coord_name]
            # Skip label auxiliary coordinates
            if hasattr(coord_var.dtype, "char") and coord_var.dtype.char == "S":
                continue
            elif coord_var.dtype == str:
                continue
            for dimension in coord_var.dimensions:
                if dimension not in coord_axis_map:
                    coord_axis_map[dimension] = "A"

        # If a dimension does not have a coordinate variable mark it as unknown
        # 'U'
        for dimension in ds.dimensions:
            if dimension not in coord_axis_map:
                coord_axis_map[dimension] = "U"

        return coord_axis_map

    def _get_coord_vars(self, ds):
        coord_vars = []
        for name, var in ds.variables.items():
            if (name,) == var.dimensions:
                coord_vars.append(name)
        return coord_vars

    def _get_dimension_order(self, ds, name, coord_axis_map):
        """
        Returns a list of strings corresponding to the named axis of the dimensions for a variable.

        Example::
            self._get_dimension_order(ds, 'temperature', coord_axis_map)
            --> ['T', 'Y', 'X']

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: Name of the variable
        :param dict coord_axis_map: A dictionary mapping each coordinate variable and dimension to a named axis

        :rtype: list
        :return: A list of strings corresponding to the named axis of the dimensions for a variable
        """

        retval = []
        variable = ds.variables[name]
        for dim in variable.dimensions:
            retval.append(coord_axis_map[dim])
        return retval

    def _get_instance_dimensions(self, ds):
        """
        Returns a list of dimensions marked as instance dimensions

        :param netCDF4.Dataset ds: An open netCDF dataset

        :rtype: list
        :returns: A list of variable dimensions
        """
        ret_val = []
        for variable in ds.get_variables_by_attributes(
            cf_role=lambda x: isinstance(x, str),
        ):
            if variable.ndim > 0:
                ret_val.append(variable.dimensions[0])
        return ret_val

    def _get_pretty_dimension_order(self, ds, name):
        """
        Returns a comma separated string of the dimensions for a specified
        variable

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: A string with a valid NetCDF variable name for the
                         dataset
        :rtype: str
        :return: A comma separated string of the variable's dimensions
        """
        dim_names = []
        for dim in ds.variables[name].dimensions:
            dim_name = dim
            if ds.dimensions[dim].isunlimited():
                dim_name += " (Unlimited)"
            dim_names.append(dim_name)
        return ", ".join(dim_names)

    def _get_pretty_dimension_order_with_type(self, ds, name, dim_types):
        """
        Returns a comma separated string of the dimensions for a specified
        variable of format "DIMENSIONS_NAME (DIMENSION_TYPE[, unlimited])"

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: A string with a valid NetCDF variable name for the
                         dataset
        :param list dim_types: A list of strings returned by
                               _get_dimension_order for the same "name"
        :rtype: str
        :return: A comma separated string of the variable's dimensions
        """
        dim_names = []
        for dim, dim_type in zip(ds.variables[name].dimensions, dim_types):
            dim_name = f"{dim} ({dim_type}"
            if ds.dimensions[dim].isunlimited():
                dim_name += ", unlimited)"
            else:
                dim_name += ")"
            dim_names.append(dim_name)
        return ", ".join(dim_names)

    def _is_station_var(self, var):
        """
        Returns True if the NetCDF variable is associated with a station, False
        otherwise.

        :param netCDF4.Variable var: a variable in an existing NetCDF dataset
        :rtype: bool
        :return: Status of whether variable appears to be associated with a
                 station
        """

        if getattr(var, "standard_name", None) in (
            "platform_name",
            "station_name",
            "instrument_name",
        ):
            return True
        return False

    def _split_standard_name(self, standard_name):
        """
        Returns a tuple of the standard_name and standard_name modifier

        Nones are used to represent the absence of a modifier or standard_name

        :rtype: tuple
        :return: 2-tuple of standard_name and modifier as strings
        """

        if isinstance(standard_name, str) and " " in standard_name:
            return standard_name.split(" ", 1)
        # if this isn't a string, then it doesn't make sense to split
        # -- treat value as standard name with no modifier
        else:
            return standard_name, None

    def check_appendix_a(self, ds):
        """
        Validates a CF dataset against the contents of its Appendix A table for
        attribute types and locations. Returns a list of results with the
        outcomes of the Appendix A validation results against the existing
        attributes in the docstring.

        :param netCDF4.Variable var: a variable in an existing NetCDF dataset
        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: A list of results corresponding to the results returned
        """
        # if 'enable_appendix_a_checks' isn't specified in the checks,
        # don't do anything on this check
        results = []
        if "enable_appendix_a_checks" not in self.options:
            return results
        possible_global_atts = set(ds.ncattrs()).intersection(self.appendix_a.keys())
        attr_location_ident = {
            "G": "global attributes",
            "C": "coordinate data",
            "D": "non-coordinate data",
        }

        def att_loc_print_helper(att_letter):
            """
            Returns a string corresponding to attr_location ident in
            human-readable form.  E.g. an input of 'G' will return
            "global attributes (G)"

            :param str att_letter: An attribute letter corresponding to the
                                   "Use" column in CF Appendix A
            :rtype: str
            :return: A string with a human-readable name followed by the input
                     letter specified
            """

            return "{} ({})".format(
                attr_location_ident.get(att_letter, "other"),
                att_letter,
            )

        def _att_loc_msg(att_loc):
            """
            Helper method for formatting an error message when an attribute
            appears in the improper location corresponding to the "Use" column
            in CF Appendix A.

            :param set att_loc: A set with the possible valid locations of the
                                attribute corresponding to the "Use" column
                                in CF Appendix A
            :rtype: str
            :return: A human-readable string with the possible valid locations
                     of the attribute
            """
            att_loc_len = len(att_loc)
            # this is a fallback in case an empty att_loc is passed
            # it generally should not occur
            valid_loc = "no locations in the dataset"
            loc_sort = sorted(att_loc)
            if att_loc_len == 1:
                valid_loc = att_loc_print_helper(loc_sort[0])
            elif att_loc_len == 2:
                valid_loc = f"{att_loc_print_helper(loc_sort[0])} and {att_loc_print_helper(loc_sort[1])}"
            # shouldn't be reached under normal circumstances, as any attribute
            # should be either G, C, or D but if another
            # category is added, this will be useful.
            else:
                valid_loc = (
                    ", ".join(loc_sort[:-1])
                    + f", and {att_loc_print_helper(loc_sort[-1])}"
                )
            return f"This attribute may only appear in {valid_loc}."

        for global_att_name in possible_global_atts:
            global_att = ds.getncattr(global_att_name)
            att_dict = self.appendix_a[global_att_name]
            att_loc = att_dict["attr_loc"]
            valid_loc_warn = _att_loc_msg(att_loc)
            if att_dict["cf_section"] is not None:
                subsection_test = ".".join(att_dict["cf_section"].split(".")[:2])

                section_loc = self.section_titles.get(
                    subsection_test,
                    att_dict["cf_section"],
                )
            else:
                section_loc = None
            test_ctx = TestCtx(BaseCheck.HIGH, section_loc)

            test_ctx.out_of += 1
            if "G" not in att_loc:
                test_ctx.messages.append(
                    f'[Appendix A] Attribute "{global_att_name}" should not be present in global (G) '
                    f"attributes. {valid_loc_warn}",
                )
            else:
                result = self._handle_dtype_check(global_att, global_att_name, att_dict)
                if not result[0]:
                    test_ctx.messages.append(result[1])
                else:
                    test_ctx.score += 1
            results.append(test_ctx.to_result())

        noncoord_vars = set(ds.variables) - set(self.coord_data_vars)
        for var_set, coord_letter in (
            (self.coord_data_vars, "C"),
            (noncoord_vars, "D"),
        ):
            for var_name in var_set:
                var = ds.variables[var_name]
                possible_attrs = set(var.ncattrs()).intersection(self.appendix_a.keys())
                for att_name in possible_attrs:
                    att_dict = self.appendix_a[att_name]
                    if att_dict["cf_section"] is not None:
                        subsection_test = ".".join(
                            att_dict["cf_section"].split(".")[:2],
                        )

                        section_loc = self.section_titles.get(
                            subsection_test,
                            att_dict["cf_section"],
                        )
                    else:
                        section_loc = None
                    test_ctx = TestCtx(BaseCheck.HIGH, section_loc, variable=var_name)
                    att_loc = att_dict["attr_loc"]
                    valid_loc_warn = _att_loc_msg(att_loc)
                    att = var.getncattr(att_name)
                    test_ctx.out_of += 1
                    if coord_letter not in att_loc:
                        test_ctx.messages.append(
                            f'[Appendix A] Attribute "{att_name}" should not be present in {att_loc_print_helper(coord_letter)} '
                            f'variable "{var_name}". {valid_loc_warn}',
                        )
                    else:
                        result = self._handle_dtype_check(att, att_name, att_dict, var)
                        if not result[0]:
                            test_ctx.messages.append(result[1])
                        else:
                            test_ctx.score += 1
                    results.append(test_ctx.to_result())

        return results

    def _check_attr_type(self, attr_name, attr_type, attribute, variable=None):
        """
        Check if an attribute `attr` is of the type `attr_type`. Upon getting
        a data type of 'D', the attr must have the same data type as the
        variable it is assigned to.

        Attributes designated type 'S' must be of type `str`. 'N' require
        numeric types, and 'D' requires the attribute type match the type
        of the variable it is assigned to.

        :param str attr_name: name of attr being checked (to format message)
        :param str attr_type: the correct type of the attribute
        :param attribute: attribute to check
        :param variable: if given, type should match attr
        :rtype tuple
        :return A two-tuple that contains pass/fail status as a boolean and
                a message string (or None if unset) as the second element.
        """

        if attr_type == "S":
            if not isinstance(attribute, str):
                return [False, f"{attr_name} must be a string"]
        else:
            # if it's not a string, it should have a numpy dtype
            underlying_dtype = getattr(attribute, "dtype", None)

            # TODO check for np.nan separately
            if underlying_dtype is None:
                return [False, f"{attr_name} must be a numeric type"]

            # both D and N should be some kind of numeric value
            is_numeric = np.issubdtype(underlying_dtype, np.number)
            if attr_type == "N":
                if not is_numeric:
                    return [False, f"{attr_name} must be a numeric type"]
            elif attr_type == "D":
                # TODO: handle edge case where variable is unset here
                temp_ctx = TestCtx()
                self._parent_var_attr_type_check(attr_name, variable, temp_ctx)
                var_dtype = getattr(variable, "dtype", None)
                if temp_ctx.messages:
                    return (
                        False,
                        f"{attr_name} must be numeric and must be equivalent to {var_dtype} dtype",
                    )
            else:
                # If we reached here, we fell off with an unrecognized type
                return (
                    False,
                    f"{attr_name} has unrecognized type '{attr_type}'",
                )
        # pass if all other possible failure conditions have been evaluated
        return (True, None)

    def _handle_dtype_check(self, attribute, attr_name, attr_dict, variable=None):
        """
        Helper function for Appendix A checks.

        :param attribute: The value of the attribute being checked
        :param str attr_name: The name of the attribute being processed
        :param dict attr_dict: The dict entry with type and attribute location
                               information corresponding to this attribute
        :param variable: if given, the variable whose type to check against
        :rtype: tuple
        :return: A two-tuple that contains pass/fail status as a boolean and
                 a message string (or None if unset) as the second element.
        """
        attr_type = attr_dict["Type"]
        if variable is None and "G" not in attr_dict["attr_loc"]:
            raise ValueError(
                "Non-global attributes must be associated with a " " variable",
            )
        attr_str = (
            f"Global attribute {attr_name}"
            if "G" in attr_dict["attr_loc"] and variable is None
            else f"Attribute {attr_name} in variable {variable.name}"
        )

        # check the type
        return_value = self._check_attr_type(attr_name, attr_type, attribute, variable)

        # if the second element is a string, format it
        if isinstance(return_value[1], str):
            return_value[1] = return_value[1].format(attr_str)

        # convert to tuple for immutability and return
        return tuple(return_value)


class CFNCCheck(BaseNCCheck, CFBaseCheck):
    """Inherits from both BaseNCCheck and CFBaseCheck to support
    checking netCDF datasets. Must inherit in this order, or certain
    attributes from BaseNCCheck (like supported_ds) will not be passed to
    CFNCCheck."""


appendix_a_base = {
    "Conventions": {"Type": "S", "attr_loc": {"G"}, "cf_section": None},
    "_FillValue": {"Type": "D", "attr_loc": {"D", "C"}, "cf_section": None},
    "add_offset": {"Type": "N", "attr_loc": {"D"}, "cf_section": "8.1"},
    "ancillary_variables": {"Type": "S", "attr_loc": {"D"}, "cf_section": "3.4"},
    "axis": {"Type": "S", "attr_loc": {"C"}, "cf_section": "4"},
    "bounds": {"Type": "S", "attr_loc": {"C"}, "cf_section": "7.1"},
    "calendar": {"Type": "S", "attr_loc": {"C"}, "cf_section": "4.4.1"},
    "cell_measures": {"Type": "S", "attr_loc": {"D"}, "cf_section": "7.2"},
    "cell_methods": {"Type": "S", "attr_loc": {"D"}, "cf_section": "7.3"},
    # cf_role type is "C" in document, which does not correspond
    # to types used, replaced with "S"
    "cf_role": {"Type": "S", "attr_loc": {"C"}, "cf_section": "9.5"},
    "climatology": {"Type": "S", "attr_loc": {"C"}, "cf_section": "7.4"},
    # comment was removed in this implementation
    "compress": {"Type": "S", "attr_loc": {"C"}, "cf_section": "8.2"},
    "coordinates": {"Type": "S", "attr_loc": {"D"}, "cf_section": "5"},
    # featureType type is "C" in document, which does not
    # correspond to types used, replaced with "S"
    "featureType": {"Type": "S", "attr_loc": {"G"}, "cf_section": "9.4"},
    "flag_masks": {"Type": "D", "attr_loc": {"D"}, "cf_section": "3.5"},
    "flag_meanings": {"Type": "S", "attr_loc": {"D"}, "cf_section": "3.5"},
    "flag_values": {"Type": "D", "attr_loc": {"D"}, "cf_section": "3.5"},
    "formula_terms": {"Type": "S", "attr_loc": {"C"}, "cf_section": "4.3.2"},
    "grid_mapping": {"Type": "S", "attr_loc": {"D"}, "cf_section": "5.6"},
    "history": {"Type": "S", "attr_loc": {"G"}, "cf_section": None},
    "institution": {"Type": "S", "attr_loc": {"G", "D"}, "cf_section": "2.6.2"},
    "leap_month": {"Type": "N", "attr_loc": {"C"}, "cf_section": "4.4.1"},
    "leap_year": {"Type": "N", "attr_loc": {"C"}, "cf_section": "4.4.1"},
    "long_name": {"Type": "S", "attr_loc": {"D", "C"}, "cf_section": "3.2"},
    "missing_value": {"Type": "D", "attr_loc": {"D", "C"}, "cf_section": "2.5.1"},
    "month_lengths": {"Type": "N", "attr_loc": {"C"}, "cf_section": "4.4.1"},
    "positive": {"Type": "S", "attr_loc": {"C"}, "cf_section": None},
    "references": {"Type": "S", "attr_loc": {"G", "D"}, "cf_section": "2.6.2"},
    "scale_factor": {"Type": "N", "attr_loc": {"D"}, "cf_section": "8.1"},
    "source": {"Type": "S", "attr_loc": {"G", "D"}, "cf_section": "2.6.2"},
    "standard_error_multiplier": {"Type": "N", "attr_loc": {"D"}, "cf_section": None},
    "standard_name": {"Type": "S", "attr_loc": {"D", "C"}, "cf_section": "3.3"},
    "title": {"Type": "S", "attr_loc": {"G"}, "cf_section": None},
    "units": {"Type": "S", "attr_loc": {"D", "C"}, "cf_section": "3.1"},
    "valid_max": {"Type": "N", "attr_loc": {"D", "C"}, "cf_section": None},
    "valid_min": {"Type": "N", "attr_loc": {"D", "C"}, "cf_section": None},
    "valid_range": {"Type": "N", "attr_loc": {"D", "C"}, "cf_section": None},
}
