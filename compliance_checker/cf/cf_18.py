
from compliance_checker.base import BaseCheck, TestCtx, Result
from compliance_checker import MemoizedDataset
from compliance_checker.cf.cf import CF1_7Check
from netCDF4 import Dataset, Variable

import sys

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

    ROOT_GROUP_ONLY_ATTRS = ['Conventions', 'external_variables']
    NON_ROOT_GROUP_OPT = ['title', 'history']

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

        results = []

        ctx_hi = TestCtx(BaseCheck.HIGH, self.section_titles["2.7"])
        ctx_lo = TestCtx(BaseCheck.LOW, self.section_titles["2.7"])

        # Make sure `Conventions` & `external_variables` attributes are only present in the
        # root group.
        for gname in ds.groups:
            ginstance = ds.createGroup(gname) # returns existing Group; doesn't create a new one

            for attr in ginstance.ncattrs():
                if attr in CF1_8Check.ROOT_GROUP_ONLY_ATTRS:

                    ctx_hi.messages.append(
                        f'§2.7.2 Attribute "{ attr }" MAY ONLY be used in the root group '
                        'and SHALL NOT be duplicated or overridden in child groups.'
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


    def check_taxa(self, ds: MemoizedDataset):
        '''
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
        be converted to a URL by adding prefixes such as ​http://www.lsid.info/. ITIS LSIDs are
        built from the ITIS Taxonomic Serial Number (TSN), such as
        "urn:lsid:itis.gov:itis_tsn:180543".

        The biological_taxon_name auxiliary coordinate variable included for human readability is
        MANDATORY. The biological_taxon_lsid auxliary coordinate variable included for software
        agent readability is optional, but strongly recommended. If both are present then each
        biological_taxon_name coordinate must exactly match the name resolved from the
        biological_taxon_lsid coordinate. If LSIDs are available for some taxa in a dataset then
        the biological_taxon_lsid auxiliary coordinate variable should be included and missing data
        given for those taxa that do not have an identifier.
        '''

        ## TODO
        # find standard names contains 'organisms_in_taxon'
        # find aux. coord. vars with standard_name='biological_taxon_name'
        # LSID,'biological_taxon_lsid' = "urn:lsid:<Authority>:<Namespace>:<ObjectID>[:<Version>]"

        results = []

        ctx_hi = TestCtx(BaseCheck.HIGH, self.section_titles["6.1"])
        ctx_med = TestCtx(BaseCheck.MEDIUM, self.section_titles["6.1"])
        ctx_lo = TestCtx(BaseCheck.LOW, self.section_titles["6.1"])

        for var in [ds[name] for name in ds.variables]:
            if 'organisms_in_taxon' in var.standard_name:

                # taxa identified by aux. coordinate vars
                for aux_var in [ds[name] for name in var.coordinates.split()]:
                    print(ds[aux_var.name])



class LSID():
    PRE = 'urn:lsid'
    OCEAN = ':'.join([PRE, 'marinespecies.org', 'taxname'])
    TERRA = ':'.join([PRE, 'itis.gov', 'itis_tsn'])

