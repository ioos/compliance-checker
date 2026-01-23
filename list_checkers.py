from compliance_checker.suite import CheckSuite
cs = CheckSuite()
cs.load_all_available_checkers()
print("Available checkers:", list(cs.checkers.keys()))
