
from compliance_checker import MemoizedDataset
from compliance_checker.cf.cf import CF1_7Check
from netCDF4 import Dataset

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
