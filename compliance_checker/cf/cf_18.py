
from compliance_checker import MemoizedDataset
from compliance_checker.cf.cf import CF1_7Check
from netCDF4 import Dataset
from compliance_checker.base import BaseCheck, BaseNCCheck, Result, TestCtx
from shapely.geometry import (MultiPoint, LineString, MultiLineString, Polygon,
                              MultiPolygon)
import numpy as np
from compliance_checker.cf.util import reference_attr_variables
import itertools

'''
What's new in CF-1.8
--------------------
2.7. Groups
  2.7.1. Scope
  2.7.2. Application of attributes

6.1.2. Taxon Names and Identifiers

7.5. Geometries
'''

class CF1_8Check(CF1_7Check):
    """ Implementation for CF v1.8. Inherits from CF1_7Check. """

    # things that are specific to 1.8
    _cc_spec_version = "1.8"
    _cc_url = "http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html"

    def __init__(self, options=None):
        super(CF1_8Check, self).__init__(options)
        self.section_titles.update({"7.5": "§7.5 Geometries"})

    def check_groups(self, ds: MemoizedDataset):
        '''
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
        '''

        # TODO make sure `Conventions` & `external_variables` attributes are only present in the
        # root group.

        print("GROUPS")
        for g in ds.groups:
            print(g)
        print()

        print("DIMENSIONS")
        for d in ds.dimensions:
            print(d)
        print()

        print("VARIABLES")
        for v in ds.variables:
            print(v)
        print()

    def check_taxa(self, ds: MemoizedDataset):
        '''
        6.1.2. Taxon Names and Identifiers

        A taxon is a named level within a biological classification, such as a class, genus and
        species. Quantities dependent on taxa have generic standard names containing the phrase
        "organisms_in_taxon", and the taxa are identified by auxiliary coordinate variables.

        The taxon auxiliary coordinate variables are string-valued. The plain-language name of the
        taxon must be contained in a variable with standard_name of biological_taxon_name. A Life
        Science Identifier (LSID) may be contained in a variable with standard_name of
        biological_taxon_lsid. This is a URN with the syntax
        "urn:lsid:<Authority>:<Namespace>:<ObjectID>[:<Version>]". This includes the reference
        classification in the <Authority> element and these are restricted by the LSID governance.
        It is strongly recommended in CF that the authority chosen is World Register of Marine
        Species (WoRMS) for oceanographic data and Integrated Taxonomic Information System (ITIS)
        for freshwater and terrestrial data. WoRMS LSIDs are built from the WoRMS AphiaID taxon
        identifier such as "urn:lsid:marinespecies.org:taxname:104464" for AphiaID 104464. This may
        be converted to a URL by adding prefixes such as ​http://www.lsid.info/. ITIS LSIDs are
        built from the ITIS Taxonomic Serial Number (TSN), such as
        "urn:lsid:itis.gov:itis_tsn:180543".

        The biological_taxon_name auxiliary coordinate variable included for human readability is
        mandatory. The biological_taxon_lsid auxliary coordinate variable included for software
        agent readability is optional, but strongly recommended. If both are present then each
        biological_taxon_name coordinate must exactly match the name resolved from the
        biological_taxon_lsid coordinate. If LSIDs are available for some taxa in a dataset then
        the biological_taxon_lsid auxiliary coordinate variable should be included and missing data
        given for those taxa that do not have an identifier.
        '''
        pass

    def check_geometry(self, ds: Dataset):
        """Runs any necessary checks for geometry well-formedness
        :param netCDF4.Dataset ds: An open netCDF dataset
        :returns list: List of error messages

        """
        vars_with_geometry = ds.get_variables_by_attributes(
                                geometry=lambda g: g is not None)
        results = []
        unique_geometry_var_names = {var.geometry for var in vars_with_geometry}
        if unique_geometry_var_names:
            geom_valid = TestCtx(BaseCheck.MEDIUM, self.section_titles["7.5"])
            geom_valid.out_of += 1
        for geometry_var_name in unique_geometry_var_names:
            if geometry_var_name not in ds.variables:
                geom_valid.messages.append("Cannot find geometry variable "
                                           f"named {geometry_var_name}")
                results.append(geom_valid.to_result())
                continue
            else:
                geometry_var = ds.variables[geometry_var_name]

            geometry_type = getattr(geometry_var, "geometry_type")
            valid_geometry_types = {"point", "line", "polygon"}
            try:
                node_coord_var_names = geometry_var.node_coordinates
            except AttributeError as e:
                geom_valid.messsages.append('Could not find required attribute '
                                            '"node_coordinates" in geometry '
                                            f'variable "{geometry_var_name}"')
                results.append(geom_valid.to_result())
            if not isinstance(node_coord_var_names, str):
                geom_valid.messages.append(
                                  'Attribute "node_coordinates" in geometry '
                                  f'variable "{geometry_var_name}" must be '
                                  'a string')
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
                                  f'{not_found_node_vars}')
                results.append(geom_valid.to_result())
                continue
                return error_msgs

            node_count = reference_attr_variables(ds,
                               getattr(geometry_var, "node_count", None))
            # multipart lines and polygons only
            part_node_count = reference_attr_variables(ds,
                                getattr(geometry_var, "part_node_count", None))
            # polygons with interior geometry only
            interior_ring = reference_attr_variables(ds,
                               getattr(geometry_var, "interior_ring", None))

            if geometry_type == "point":
                geometry = PointGeometry(node_coord_vars, node_count)
            elif geometry_type == "line":
                geometry = LineGeometry(node_coord_vars, node_count,
                                        part_node_count)
            elif geometry_type == "polygon":
                geometry = PolygonGeometry(node_coord_vars, node_count,
                                                    part_node_count,
                                                    interior_ring)
            else:
                geom_valid.messages.append(
                                  f'For geometry variable "{geometry_var_name}'
                                  'the attribute "geometry_type" must exist'
                                  'and have one of the following values:'
                                  '"point", "line", "polygon"')
                results.append(geom_valid.to_result())
                continue
            # check geometry
            geometry.check_geometry()
            if not geometry.errors: #geom_valid.messages:
                geom_valid.score += 1
            results.append(geom_valid.to_result())
        return results

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
            self.errors.append("The following geometry variables "
                               f"have non-numeric contents: {invalid_vars}")

    def _split_mulitpart_geometry(self):
        arr_extents_filt = self.part_node_count[self.part_node_count > 0]
        splits = np.split(np.vstack(self.coord_vars).T,
                          arr_extents_filt.cumsum()[:-1])
        return splits

class PointGeometry(GeometryStorage):
    """Class for validating Point/MultiPoint geometries"""

    def check_geometry(self):
        super().check_geometry()
        # non-multipoint should have exactly one feature
        if self.node_count is None:
            expected_node_count = 1
        else:
            expected_node_count = self.node_count

        if all(len(cv.dimensions) != 0 for cv in self.coord_vars):
            same_dim_group = itertools.groupby(self.coord_vars,
                                            lambda x: x.dimensions)
            same_dim = (next(same_dim_group, True) and
                             not next(same_dim_group, False))
            if not same_dim:
                self.errors.append("For a point geometry, coordinate "
                                   "variables must be the same length as "
                                   "node_count defined, or must be "
                                   "length 1 if node_count is not set")
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
        same_dim_group = itertools.groupby(self.coord_vars,
                                           lambda x: x.dimensions)
        same_dim = (next(same_dim_group, True) and
                         not next(same_dim_group, False))
        if not same_dim:
            raise IndexError("Coordinate variables must be the same length. "
                            "If node_count is specified, this value must "
                            "also sum to the length of the coordinate "
                            "variables.")
        # if a multipart
        if self.node_count is not None:
            same_length = len(self.coord_vars[0]) == self.node_count[:].sum()
            if not same_length:
                geom_errors.append("Coordinate variables must be the same "
                                   "length. If node_count is specified, this "
                                   "value must also sum to the length of the "
                                    "coordinate variables.")
        if self.part_node_count is not None:
            if not np.issubdtype(self.part_node_count.dtype, np.integer):
                geom_errors.append("when part_node_count is specified, it must "
                                   "be an array of integers")
            same_node_count = len(self.coord_vars[0]) == self.node_count[:].sum()
            if not same_node_count:
                geom_errors.append("The sum of part_node_count must be equal "
                                   "to the value of node_count")
        return geom_errors

class PolygonGeometry(LineGeometry):
    """Class for validating Line/MultiLine geometries"""
    # TODO/clarify: Should polygons be simple, i.e. non-self intersecting?
    # Presumably
    def __init__(self, coord_vars, node_count, part_node_count,
                    interior_ring):
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
            raise ValueError("Polygon contains too few points to perform orientation test")

        ccw = polygon.exterior.is_ccw
        return not ccw if interior else ccw

    def check_geometry(self):
        messages = super().check_geometry()
        # If any errors occurred within the preliminary checks, they preclude
        # running checks against the geometry here.
        if messages:
            return messages
        if self.part_node_count is not None:
            extents = np.concatenate([np.array([0]),
                                     self.part_node_count[:].cumsum()])
            if self.interior_ring is not None:
                ring_orientation = self.interior_ring[:].astype(bool)
            else:
                ring_orientation = np.zeros(len(self.part_count), dtype=bool)
            current_node_count = self.node_count[:].copy()
            node_indexer_len = len(self.part_node_count)
        else:
            extents = np.concatenate([np.array([0]),
                                     self.node_count[:].cumsum()])
            node_indexer_len = len(self.node_count)
            ring_orientation = np.zeros(node_indexer, dtype=bool)
        # TODO: is it necessary to check whether part_node_count "consumes"
        #       node_count in the polygon, i.e. first (3, 3, 3) will consume
        #       a node part of 9, follow by next 3 will consume a node part of
        #       3 after consuming
        for i in range(node_indexer_len):
            extent_slice = slice(extents[i], extents[i+1])
            poly_sliced = np.vstack([cv[extent_slice] for cv in
                                     self.coord_vars]).T
            pass_orientation = (self.check_polygon_orientation(
                                            poly_sliced,
                                            ring_orientation[i]))
            if not pass_orientation:
                orient_fix = (("exterior", "counterclockwise")
                                 if not ring_orientation[i] else
                               ("interior", "clockwise"))
                message = (f"An {orient_fix[0]} polygon referred to by "
                           f"coordinates ({poly_sliced}) must have coordinates "
                           f"in {orient_fix[1]} order")
                messages.append(message)
        return messages


if __name__ == "__main__":

    import os
    checker = CF1_8Check()

    _dir = '/Users/ben/Downloads/netcdf/examples/'
    ncfs = [ f for f in os.listdir(_dir) if f.endswith('.nc') ]

    for ncf in ncfs:
        ds = Dataset(_dir + ncf)

        print('-' * 210)
        print(ncf)
        print()

        checker.check_groups(ds)
