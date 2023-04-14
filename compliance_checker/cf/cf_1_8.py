"""
What's new in CF-1.8
--------------------
2.7. Groups
  2.7.1. Scope
  2.7.2. Application of attributes

6.1.2. Taxon Names and Identifiers

7.5. Geometries
"""

import itertools
import re
import warnings

import numpy as np
import requests
from lxml import etree
from netCDF4 import Dataset
from shapely.geometry import Polygon

from compliance_checker import MemoizedDataset
from compliance_checker.base import BaseCheck, TestCtx
from compliance_checker.cf.cf_1_7 import CF1_7Check
from compliance_checker.cf.util import reference_attr_variables, string_from_var_type


class CF1_8Check(CF1_7Check):
    """Implementation for CF v1.8. Inherits from CF1_7Check."""

    # things that are specific to 1.8
    _cc_spec_version = "1.8"
    _cc_url = "http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html"

    ROOT_GROUP_ONLY_ATTRS = ["Conventions", "external_variables"]
    NON_ROOT_GROUP_OPT = ["title", "history"]

    def __init__(self, options=None):
        super(CF1_8Check, self).__init__(options)
        self.section_titles.update(
            {
                "2.7": "§2.7 Groups",
                "6.1.2": "§6.1.2 Taxon Names and Identifiers",
                "7.5": "§7.5 Geometries",
            }
        )

    def check_groups(self, ds: MemoizedDataset):
        """
        2.7.2. Application of attributes

        The following attributes are optional for non-root groups. They are allowed in order to
        provide additional provenance and description of the subsidiary data. They do not override
        attributes from parent groups.

        - title
        - history

        If these attributes are present, they may be applied additively to the parent attributes of
        the same name. If a file containing groups is modified, the user or application need only
        update these attributes in the root group, rather than traversing all groups and updating
        all attributes that are found with the same name. In the case of conflicts, the root group
        attribute takes precedence over per-group instances of these attributes.

        The following attributes MAY ONLY be used in the root group and SHALL NOT be duplicated or
        overridden in child groups:

        - Conventions
        - external_variables

        Furthermore, per-variable attributes MUST be attached to the variables to which they refer.
        They MAY NOT be attached to a group, even if all variables within that group use the same
        attribute and value.

        If attributes are present within groups without being attached to a variable, these
        attributes apply to the group where they are defined, and to that group’s descendants, but
        not to ancestor or sibling groups. If a group attribute is defined in a parent group, and
        one of the child group redefines the same attribute, the definition within the child group
        applies for the child and all of its descendants.
        """

        results = []

        ctx_hi = TestCtx(BaseCheck.HIGH, self.section_titles["2.7"])
        ctx_lo = TestCtx(BaseCheck.LOW, self.section_titles["2.7"])

        # Make sure `Conventions` & `external_variables` attributes are only present in the
        # root group.
        for gname in ds.groups:
            ginstance = ds.createGroup(
                gname
            )  # returns existing Group; doesn't create a new one

            for attr in ginstance.ncattrs():
                if attr in CF1_8Check.ROOT_GROUP_ONLY_ATTRS:

                    ctx_hi.messages.append(
                        f'§2.7.2 Attribute "{ attr }" MAY ONLY be used in the root group '
                        "and SHALL NOT be duplicated or overridden in child groups."
                    )

                    results.append(ctx_hi.to_result())

                elif attr in CF1_8Check.NON_ROOT_GROUP_OPT:

                    ctx_lo.messages.append(
                        f"§2.7.2 Note: attribute '{ attr }' found on non-root group '{ gname }'. "
                        "This is optional for non-root groups. It is allowed in order to provide additional "
                        "provenance and description of the subsidiary data. It does not override "
                        "attributes from parent groups."
                    )
                    results.append(ctx_lo.to_result())

        return results

    def check_geometry(self, ds: Dataset):
        """Runs any necessary checks for geometry well-formedness
        :param netCDF4.Dataset ds: An open netCDF dataset
        :returns list: List of error messages
        """
        vars_with_geometry = ds.get_variables_by_attributes(
            geometry=lambda g: g is not None
        )
        results = []
        unique_geometry_var_names = {var.geometry for var in vars_with_geometry}
        if unique_geometry_var_names:
            geom_valid = TestCtx(BaseCheck.MEDIUM, self.section_titles["7.5"])
            geom_valid.out_of += 1
        for geometry_var_name in unique_geometry_var_names:
            if geometry_var_name not in ds.variables:
                geom_valid.messages.append(
                    "Cannot find geometry variable " f"named {geometry_var_name}"
                )
                results.append(geom_valid.to_result())
                continue
            else:
                geometry_var = ds.variables[geometry_var_name]

            geometry_type = getattr(geometry_var, "geometry_type")
            try:
                node_coord_var_names = geometry_var.node_coordinates
            except AttributeError:
                geom_valid.messages.append(
                    "Could not find required attribute "
                    '"node_coordinates" in geometry '
                    f'variable "{geometry_var_name}"'
                )
                results.append(geom_valid.to_result())
            if not isinstance(node_coord_var_names, str):
                geom_valid.messages.append(
                    'Attribute "node_coordinates" in geometry '
                    f'variable "{geometry_var_name}" must be '
                    "a string"
                )
                results.append(geom_valid.to_result())
                continue
            split_coord_names = node_coord_var_names.strip().split(" ")
            node_coord_vars, not_found_node_vars = [], []
            for coord_var_name in split_coord_names:
                try:
                    node_coord_vars.append(ds.variables[coord_var_name])
                except KeyError:
                    not_found_node_vars.append(coord_var_name)
            # If any variables weren't found, we can't continue
            if not_found_node_vars:
                geom_valid.messages.append(
                    "The following referenced node coordinate"
                    "variables for geometry variable"
                    f'"{geometry_var_name}" were not found: '
                    f"{not_found_node_vars}"
                )
                results.append(geom_valid.to_result())
                continue

            node_count = reference_attr_variables(
                ds, getattr(geometry_var, "node_count", None)
            )
            # multipart lines and polygons only
            part_node_count = reference_attr_variables(
                ds, getattr(geometry_var, "part_node_count", None)
            )
            # polygons with interior geometry only
            interior_ring = reference_attr_variables(
                ds, getattr(geometry_var, "interior_ring", None)
            )

            if geometry_type == "point":
                geometry = PointGeometry(node_coord_vars, node_count)
            elif geometry_type == "line":
                geometry = LineGeometry(node_coord_vars, node_count, part_node_count)
            elif geometry_type == "polygon":
                geometry = PolygonGeometry(
                    node_coord_vars, node_count, part_node_count, interior_ring
                )
            else:
                geom_valid.messages.append(
                    f'For geometry variable "{geometry_var_name}'
                    'the attribute "geometry_type" must exist'
                    "and have one of the following values:"
                    '"point", "line", "polygon"'
                )
                results.append(geom_valid.to_result())
                continue
            # check geometry
            messages = geometry.check_geometry()
            if not messages:
                geom_valid.score += 1
            else:
                geom_valid.messages.extend(messages)
            results.append(geom_valid.to_result())
        return results

    def check_taxa(self, ds: Dataset):
        """
        6.1.2. Taxon Names and Identifiers

        A taxon is a named level within a biological classification, such as a class, genus and
        species. QUANTITIES DEPENDENT ON TAXA HAVE GENERIC STANDARD NAMES CONTAINING THE PHRASE
        "organisms_in_taxon", AND THE TAXA ARE IDENTIFIED BY AUXILIARY COORDINATE VARIABLES.

        The taxon auxiliary coordinate variables are string-valued. The plain-language name of the
        taxon MUST be contained in a variable with standard_name of 'biological_taxon_name'. A Life
        Science Identifier (LSID) may be contained in a variable with standard_name of
        biological_taxon_lsid. This is a URN with the syntax
        "urn:lsid:<Authority>:<Namespace>:<ObjectID>[:<Version>]". This includes the reference
        classification in the <Authority> element and these are restricted by the LSID governance.
        It is strongly recommended in CF that the authority chosen is World Register of Marine
        Species (WoRMS) for oceanographic data and Integrated Taxonomic Information System (ITIS)
        for freshwater and terrestrial data. WoRMS LSIDs are built from the WoRMS AphiaID taxon
        identifier such as "urn:lsid:marinespecies.org:taxname:104464" for AphiaID 104464. This may
        be converted to a URL by adding prefixes such as http://www.lsid.info/. ITIS LSIDs are
        built from the ITIS Taxonomic Serial Number (TSN), such as
        "urn:lsid:itis.gov:itis_tsn:180543".

        The biological_taxon_name auxiliary coordinate variable included for human readability is
        MANDATORY. The biological_taxon_lsid auxliary coordinate variable included for software
        agent readability is optional, but strongly recommended. If both are present then each
        biological_taxon_name coordinate must exactly match the name resolved from the
        biological_taxon_lsid coordinate. If LSIDs are available for some taxa in a dataset then
        the biological_taxon_lsid auxiliary coordinate variable should be included and missing data
        given for those taxa that do not have an identifier.
        """
        ret_val = []

        def match_taxa_standard_names(standard_name_string):
            """
            Match variables which are standard_names related to taxa, but
            are not the taxon identifiers or LSIDs themselves.
            """
            return (
                standard_name_string is not None
                and "taxon" in standard_name_string
                and standard_name_string  # exclude the identifiers we just looked at
                not in {"biological_taxon_lsid", "biological_taxon_name"}
                and standard_name_string in self._std_names
            )

        taxa_quantifier_variables = ds.get_variables_by_attributes(
            standard_name=match_taxa_standard_names
        )
        # If there are no matches, there either are no taxa variables
        # or the standard names are not appropriate, which will be picked up
        # by the standard_name check
        if not taxa_quantifier_variables:
            return

        for taxon_quantifier_variable in taxa_quantifier_variables:
            valid_taxa = TestCtx(BaseCheck.HIGH, self.section_titles["6.1.2"])
            if not isinstance(
                getattr(taxon_quantifier_variable, "coordinates", None), str
            ):
                valid_taxa.add_failure(
                    f'{taxon_quantifier_variable.name} must have a string valued "coordinates" attribute'
                )
                continue

            coordinate_var_names = taxon_quantifier_variable.coordinates.split(" ")
            invalid_coord_vars = set(coordinate_var_names) - ds.variables.keys()
            if invalid_coord_vars:
                valid_taxa.add_failure(
                    'The following values for "coordinates" attributes were not found in the dataset\'s variables '
                    f"{invalid_coord_vars}"
                )

            if len(coordinate_var_names) > 2:
                valid_taxa.add_failure(
                    "coordinates attribute for taxon data must either reference one or two variable names"
                )
                continue

            coordinate_vars = [
                ds.variables[var_name] for var_name in coordinate_var_names
            ]

            coord_var_standard_names = {
                var: getattr(var, "standard_name", None) for var in coordinate_vars
            }
            # if we have no authority, we can't check validity of the name -- assume it's OK
            standard_name_set = set(coord_var_standard_names.values())
            if set(coord_var_standard_names.keys()) == {"biological_taxon_name"}:
                # TODO: Check for at least binomial nomenclature?
                continue
            # check against WoRMS or ITIS if applicable
            elif standard_name_set == {
                "biological_taxon_name",
                "biological_taxon_lsid",
            }:
                inverted_dict = {v: k for k, v in coord_var_standard_names.items()}
                taxon_lsid_var = inverted_dict["biological_taxon_lsid"]
                taxon_name_var = inverted_dict["biological_taxon_name"]
                lsid_messages = self.handle_lsid(taxon_lsid_var, taxon_name_var)
                valid_taxa.out_of += 1
                if lsid_messages:
                    valid_taxa.messages.extend(lsid_messages)
                else:
                    valid_taxa.score += 1
            else:
                valid_taxa.add_failure(
                    f"coordinates attribute for variable {taxon_quantifier_variable} must consist of "
                    'variables containing standard names of either just "biological_taxon_name", or "biological_taxon_name" and "biological_taxon_identifier"'
                )
            ret_val.append(valid_taxa.to_result())

        return ret_val

    def handle_lsid(self, taxon_lsid_variable, taxon_name_variable):
        """
        Checks if LSID is well formed and present in the LSID database,
        and then attempts to delegate to WoRMS or ITIS, the LSID is applicable.
        If the LSID does not check the above authorities, it is not
        currently checked for correctness.
        """
        messages = []
        match_str = (
            r"(?:http://(?:www\.)?lsid.info/)?urn:lsid:"
            r"(?P<authority>[^:]+):(?P<namespace>[^:]+):"
            r"(?P<object_id>\w+)(?::(?P<version>\w+))?"
        )
        for taxon_lsid, taxon_name in zip(
            taxon_lsid_variable[:], taxon_name_variable[:]
        ):
            # TODO: handle case where LSID is not present.  This can happen
            #       if the species is not present in the database desired.
            taxon_name_str = string_from_var_type(taxon_name)
            lsid_str = string_from_var_type(taxon_lsid)
            # if nodata/empty string for LSID, skip validity check
            if lsid_str == "":
                continue
            taxon_match = re.fullmatch(match_str, lsid_str)
            if not taxon_match:
                messages.append(
                    "Taxon id must match one of the following forms:\n"
                    "- urn:lsid:<authority>:<namespace>:<object_id>\n"
                    "- urn:lsid:<authority>:<namespace>:<object_id>:<version>\n"
                    "- www.lsid.info/urn:lsid.info:<authority>:<namespace>/<object_id>\n"
                    "- www.lsid.info/urn:lsid.info:<authority>:<namespace>/<object_id>:<version>\n"
                    "- lsid.info/urn:lsid.info:<authority>:<namespace>/<object_id>\n"
                    "- lsid.info/urn:lsid.info:<authority>:<namespace>/<object_id>:<version>\n"
                    "- http://lsid.info/urn:lsid.info:<authority>:<namespace>/<object_id>\n"
                    "- http://lsid.info/urn:lsid.info:<authority>:<namespace>/<object_id>:<version>\n"
                    "- http://www.lsid.info/urn:lsid.info:<authority>:<namespace>/<object_id>\n"
                    "- http://www.lsid.info/urn:lsid.info:<authority>:<namespace>/<object_id>:<version>"
                )
                continue
            if lsid_str.startswith("urn"):
                lsid_url = f"http://www.lsid.info/{lsid_str}"
            else:
                lsid_url = lsid_str

            try:
                response = requests.get(lsid_url, timeout=10)
                response.raise_for_status()
            except requests.exceptions.RequestException as e:
                # 400 error code indicates something is malformed on client's
                # end
                if response.status_code == 400:
                    tree = etree.HTML(response.text)
                    problem_text = tree.find("./body/p").text
                    messages.append(
                        "http://lsid.info returned an error message "
                        f"for submitted LSID string '{lsid_str}': "
                        f"{problem_text}"
                    )
                else:
                    messages.append(
                        "Error occurred attempting to check LSID "
                        f"'{lsid_str}': {str(e)}"
                    )
                continue

            # WoRMS -- marine bio data
            if (
                taxon_match["authority"] == "marinespecies.org"
                and taxon_match["namespace"] == "taxname"
            ):
                try:
                    response = requests.get(
                        f"http://www.marinespecies.org/rest/AphiaRecordByAphiaID/{taxon_match['object_id']}",
                        timeout=15,
                    )
                    response.raise_for_status()
                except requests.exceptions.RequestException as e:  # noqa: F841
                    messages.append(
                        "Aphia ID {taxon_match['object_id'] returned "
                        "other error: {str(e)}"
                    )
                # record not found in database
                if response.status_code == 204:
                    messages.append(
                        "Aphia ID {taxon_match['object_id'] "
                        "not found in WoRMS database"
                    )
                # good case, parse JSON
                elif response.status_code == 200:
                    valid_name = response.json()["valid_name"]
                    if valid_name != taxon_name_str:
                        messages.append(
                            "Supplied taxon name and WoRMS valid name do not match. "
                            f"Supplied taxon name is '{taxon_name_str}', WoRMS valid name "
                            f"is '{valid_name}.'"
                        )
                # Misc non-error code.  Should not reach here.
                else:
                    messages.append(
                        f"Aphia ID {taxon_match['object_id']}"
                        "returned an unhandled HTTP status "
                        f"code {response.status_code}"
                    )
                    continue

            # ITIS -- freshwater bio data
            elif (
                taxon_match["authority"] == "itis.gov"
                and taxon_match["namespace"] == "itis_tsn"
            ):
                itis_url = f"https://www.itis.gov/ITISWebService/jsonservice/getFullRecordFromTSN?tsn={taxon_match['object_id']}"
                try:
                    itis_response = requests.get(itis_url, timeout=15)
                    itis_response.raise_for_status()
                except requests.exceptions.RequestException as e:
                    if itis_response.status_code == 404:
                        messages.append(
                            "itis.gov TSN " f"{taxon_match['object_id']} not found."
                        )
                        continue
                    else:
                        messages.append(
                            "itis.gov identifier returned other " f"error: {str(e)}"
                        )
                        continue
                json_contents = itis_response.json()
                combined_name = json_contents["scientificName"]["combinedName"]

                if taxon_name_str != combined_name:
                    messages.append(
                        "Supplied taxon name and ITIS scientific name do not match. "
                        f"Supplied taxon name is '{taxon_name_str}', ITIS scientific name "
                        f"for TSN {taxon_match['object_id']} is '{combined_name}.'"
                    )

            else:
                warnings.warn(
                    "Compliance checker only supports checking valid "
                    "LSID URNs of the form "
                    "'urn:lsid:marinespecies.org:taxname:<AphiaID>' or "
                    "'urn:lsid:itis.gov:itis_tsn:<TSN>'.  Assuming "
                    "pass condition"
                )

        return messages


class GeometryStorage(object):
    """Abstract base class for geometries"""

    def __init__(self, coord_vars, node_count):
        self.coord_vars = coord_vars
        self.node_count = node_count
        self.errors = []
        # geometry is later parsed after sanity checks are run
        self.geometry = None

    def check_geometry(self):
        invalid_vars = []
        for coord_var in self.coord_vars:
            if not np.issubdtype(coord_var, np.float):
                invalid_vars.append(coord_var.name)
        # can't continue if the geometry variables are not the correct size
        if invalid_vars:
            self.errors.append(
                "The following geometry variables "
                f"have non-numeric contents: {invalid_vars}"
            )

    def _split_mulitpart_geometry(self):
        arr_extents_filt = self.part_node_count[self.part_node_count > 0]
        splits = np.split(np.vstack(self.coord_vars).T, arr_extents_filt.cumsum()[:-1])
        return splits


class PointGeometry(GeometryStorage):
    """Class for validating Point/MultiPoint geometries"""

    def check_geometry(self):
        super().check_geometry()

        if all(len(cv.dimensions) != 0 for cv in self.coord_vars):
            same_dim_group = itertools.groupby(self.coord_vars, lambda x: x.dimensions)
            same_dim = next(same_dim_group, True) and not next(same_dim_group, False)
            if not same_dim:
                self.errors.append(
                    "For a point geometry, coordinate "
                    "variables must be the same length as "
                    "node_count defined, or must be "
                    "length 1 if node_count is not set"
                )
        return self.errors


class LineGeometry(GeometryStorage):
    """Class for validating Line/MultiLine geometries"""

    def __init__(self, coord_vars, node_count, part_node_count):
        super().__init__(coord_vars, node_count)
        self.part_node_count = part_node_count
        if not np.issubdtype(self.node_count.dtype, np.integer):
            raise TypeError("For line geometries, node_count must be an integer")

    def check_geometry(self):
        geom_errors = []
        same_dim_group = itertools.groupby(self.coord_vars, lambda x: x.dimensions)
        same_dim = next(same_dim_group, True) and not next(same_dim_group, False)
        if not same_dim:
            raise IndexError(
                "Coordinate variables must be the same length. "
                "If node_count is specified, this value must "
                "also sum to the length of the coordinate "
                "variables."
            )
        # if a multipart
        if self.node_count is not None:
            same_length = len(self.coord_vars[0]) == self.node_count[:].sum()
            if not same_length:
                geom_errors.append(
                    "Coordinate variables must be the same "
                    "length. If node_count is specified, this "
                    "value must also sum to the length of the "
                    "coordinate variables."
                )
        if self.part_node_count is not None:
            if not np.issubdtype(self.part_node_count.dtype, np.integer):
                geom_errors.append(
                    "when part_node_count is specified, it must "
                    "be an array of integers"
                )
            same_node_count = len(self.coord_vars[0]) == self.node_count[:].sum()
            if not same_node_count:
                geom_errors.append(
                    "The sum of part_node_count must be equal "
                    "to the value of node_count"
                )
        return geom_errors


class PolygonGeometry(LineGeometry):
    """Class for validating Line/MultiLine geometries"""

    # TODO/clarify: Should polygons be simple, i.e. non-self intersecting?
    # Presumably
    def __init__(self, coord_vars, node_count, part_node_count, interior_ring):
        super().__init__(coord_vars, node_count, part_node_count)
        self.part_node_count = part_node_count
        self.interior_ring = interior_ring

    def check_polygon_orientation(self, transposed_coords, interior=False):
        """
        Checks that the polygon orientation is counter-clockwise if an
        exterior ring, otherwise clockwise if an interior ring.  Orientation
        is indicated by the `interior` boolean variable with False for an
        exterior ring and True for an interior ring (hole), defaulting to False.
        This function operates piecewise on individual interior/exterior
        polygons as well as multipart polygons
        :param np.array transposed_coords: A 2-by-n array of x and y coordinates
        :param bool interior: A boolean defaulting to False which has False
        indicating a counter-clockwise or exterior polygon, and True
        indicating a clockwise or interior polygon.
        :rtype bool:
        :returns: True if the polygon follows the proper orientation,
                  False if it fails the orientation test.
        """

        try:
            polygon = Polygon(transposed_coords.tolist())
        except ValueError:
            raise ValueError(
                "Polygon contains too few points to perform orientation test"
            )

        ccw = polygon.exterior.is_ccw
        return not ccw if interior else ccw

    def check_geometry(self):
        messages = super().check_geometry()
        # If any errors occurred within the preliminary checks, they preclude
        # running checks against the geometry here.
        if messages:
            return messages
        if self.part_node_count is not None:
            extents = np.concatenate([np.array([0]), self.part_node_count[:].cumsum()])
            if self.interior_ring is not None:
                ring_orientation = self.interior_ring[:].astype(bool)
            else:
                ring_orientation = np.zeros(len(self.part_count), dtype=bool)
            node_indexer_len = len(self.part_node_count)
        else:
            extents = np.concatenate([np.array([0]), self.node_count[:].cumsum()])
            node_indexer_len = len(self.node_count)
            ring_orientation = np.zeros(node_indexer_len, dtype=bool)
        # TODO: is it necessary to check whether part_node_count "consumes"
        #       node_count in the polygon, i.e. first (3, 3, 3) will consume
        #       a node part of 9, follow by next 3 will consume a node part of
        #       3 after consuming
        for i in range(node_indexer_len):
            extent_slice = slice(extents[i], extents[i + 1])
            poly_sliced = np.vstack([cv[extent_slice] for cv in self.coord_vars]).T
            pass_orientation = self.check_polygon_orientation(
                poly_sliced, ring_orientation[i]
            )
            if not pass_orientation:
                orient_fix = (
                    ("exterior", "counterclockwise")
                    if not ring_orientation[i]
                    else ("interior", "clockwise")
                )
                message = (
                    f"An {orient_fix[0]} polygon referred to by "
                    f"coordinates ({poly_sliced}) must have coordinates "
                    f"in {orient_fix[1]} order"
                )
                messages.append(message)
        return messages
