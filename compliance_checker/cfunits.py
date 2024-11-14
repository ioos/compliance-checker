from pyudunits2 import UnitSystem, UnresolvableUnitException

try:
    import cf_units
except ImportError:
    cf_units = False


class PyUdunits2:
    """Workaround for the differences in pyudunits2 and cf-units.

    NB: Some of these may change and/or get implemented upstream. Pyudunits2 is new and in-flux.

    1/4 Raise the same ValueError to match cf-unit errors.
    2/4 Creates an empty unit from None to mimic cf-unit's Unit('unknown')
    3/4 Add a definition object that is ust units.expanded()
    """

    def __init__(self, units):
        """Keep unit system so we can convert from string later."""
        self.ut_system = UnitSystem.from_udunits2_xml()

        if units is None:
            units = ""

        try:
            self.units = self.ut_system.unit(units)
        except (SyntaxError, UnresolvableUnitException) as err:
            raise ValueError from err
        self.definition = self.units.expanded()

    def __eq__(self, other):
        return self.units == other

    def is_dimensionless(self):
        return self.units.is_dimensionless()

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
        if _is_time_reference(self.units) and _is_time_reference(other):
            convertible = True
        # One is time, the other is not, change it to False.
        if sum((_is_time_reference(self.units), _is_time_reference(other))) == 1:
            convertible = False

        return convertible

    def is_time_reference(self):
        return _is_time_reference(self.units)


def _is_time_reference(self):
    # FIXME: cf-units Workaround 4/4 -> cf_units can differentiante between time reference and time units.
    is_time_reference = False
    try:
        if hasattr(self._definition, "shift_from"):
            is_time_reference = True
    except KeyError:
        # FIXME: hasattr should return None in that case.
        # pyudunits2/_expr_graph.py:27, in Node.__getattr__(self, name)
        #      25 def __getattr__(self, name):
        #      26     # Allow the dictionary to raise KeyError if the key doesn't exist.
        # ---> 27     return self._attrs[name]
        # KeyError: 'shift_from'
        pass
    return is_time_reference


if cf_units:
    PyUdunits2 = cf_units.Unit


class Unit(PyUdunits2):
    def __init__(self, units):
        super().__init__(units)
