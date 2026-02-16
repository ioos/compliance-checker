import pytest
from compliance_checker.suite import CheckSuite
from compliance_checker import run_checker


class TestSimpleAPI:
    
    def teardown_method(self):
        # Reset state after each test
        CheckSuite._checkers_loaded = False
        CheckSuite._instance = None

    def test_lazy_loading(self):
        """Test that checkers are not loaded initially, but load on demand."""
        # Reset to ensure clean state
        CheckSuite._checkers_loaded = False
        CheckSuite._instance = None
        
        # Verify not loaded yet
        assert CheckSuite._checkers_loaded is False
        
        # Trigger loading via factory
        CheckSuite._get_instance()
        
        # Verify loaded
        assert CheckSuite._checkers_loaded is True
        assert CheckSuite._instance is not None

    def test_run_checker_auto_loads(self):
        """Test that run_checker automatically loads checkers."""
        """Test that run_checker automatically loads checkers."""
        CheckSuite._checkers_loaded = False
        CheckSuite._instance = None
        
        # We try to run the checker. It will fail because the file doesn't exist,
        # but we just want to verify that it *tried* (which means it loaded checkers).
        try:
             run_checker("dummy.nc", "acdd")
        except Exception:
             # We expect a crash (FileNotFoundError), so we catch it and ignore it.
             pass
             
        # The key test: Did checkers get loaded?
        assert CheckSuite._checkers_loaded is True

    def test_backward_compatibility(self):
        """Test that the old explicit loading method still works."""
        CheckSuite._checkers_loaded = False
        CheckSuite._instance = None
        
        cs = CheckSuite()
        cs.load_all_available_checkers()
        
        assert CheckSuite._checkers_loaded is True