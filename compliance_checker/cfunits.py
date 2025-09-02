try:
    import cf_units
except ImportError:
    cf_units = False

UT_SYSTEM = None

if cf_units:
    PyUdunits2 = cf_units.Unit
else:
    import pyudunits2

    class PyUdunits2:
        """Workaround for the differences in pyudunits2 and cf-units.

        NB: Some of these may change and/or get implemented upstream. Pyudunits2 is new and in-flux.

        1/4 Raise the same ValueError to match cf-unit errors.
        2/4 Creates an empty unit from None to mimic cf-unit's Unit('unknown')
        3/4 Add a definition object that is ust units.expanded()
        """

        def __init__(self, units: str | None):
            """Keep unit system so we can convert from string later."""
            global UT_SYSTEM
            if UT_SYSTEM is None:
                UT_SYSTEM = pyudunits2.UnitSystem.from_udunits2_xml()

            self.ut_system = UT_SYSTEM

            if units is None:
                units = ""

            try:
                self.units = self.ut_system.unit(units)
            except (SyntaxError, pyudunits2.UnresolvableUnitException) as err:
                raise ValueError from err
            self.definition = self.units.expanded()

        def __eq__(self, other):
            return self.units == other

        def is_dimensionless(self):
            return self.units.is_dimensionless()

        def is_time_reference(self):
            return isinstance(self.units, pyudunits2.DateUnit)

        def is_convertible(self, other):
            if isinstance(other, str):
                other = self.ut_system.unit(other)
            elif isinstance(other, (PyUdunits2)):
                other = other.units
            else:
                msg = f"Expected valid unit string or pyudunits2 unit object. Got {other}."
                raise ValueError(msg)

            # FIXME: cf-units Workaround 1/4 -> cf_units.Unit(None) -> Unit('unknown').
            if "" in (self.units.expanded(), other.expanded()):
                return False

            convertible = self.units.is_convertible_to(other)
            # FIXME: cf-units Workaround 2/4 -> time is not convertible to time reference.

            # Both are time reference confirm.
            if self.is_time_reference() and isinstance(other, pyudunits2.DateUnit):
                convertible = True

            return convertible


class Unit(PyUdunits2):
    def __init__(self, units):
        super().__init__(units)
