#!/usr/bin/env python
# -*- coding: utf-8 -*-

import compliance_checker.cf.util as cfutil  # noqa: F401
from compliance_checker.base import (  # noqa: F401
    BaseCheck,
    BaseNCCheck,
    Result,
    TestCtx,
)
from compliance_checker.cf import util  # noqa: F401
from compliance_checker.cf.appendix_d import (  # noqa: F401
    dimless_vertical_coordinates_1_6,
    dimless_vertical_coordinates_1_7,
    no_missing_terms,
)
from compliance_checker.cf.appendix_e import cell_methods16  # noqa: F401
from compliance_checker.cf.appendix_e import cell_methods17  # noqa: F401
from compliance_checker.cf.appendix_f import ellipsoid_names17  # noqa: F401
from compliance_checker.cf.appendix_f import (  # noqa: F401
    grid_mapping_attr_types16,
    grid_mapping_attr_types17,
    grid_mapping_dict16,
    grid_mapping_dict17,
    horizontal_datum_names17,
    prime_meridian_names17,
)

# Version specific checkers organized in other modules
from compliance_checker.cf.cf_1_6 import CF1_6Check  # noqa: F401
from compliance_checker.cf.cf_1_7 import CF1_7Check  # noqa: F401
from compliance_checker.cf.cf_1_8 import CF1_8Check  # noqa: F401
from compliance_checker.cf.cf_1_9 import CF1_9Check  # noqa: F401
from compliance_checker.cf.cf_1_10 import CF1_10Check  # noqa: F401
from compliance_checker.cf.cf_1_11 import CF1_11Check  # noqa: F401
