from compliance_checker.cf.cf_1_9 import CF1_9Check


# no significant features for code implementation between CF 1.9 and CF 1.10
class CF1_10Check(CF1_9Check):
    _cc_spec_version = "1.10"
    _cc_url = "http://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html"
