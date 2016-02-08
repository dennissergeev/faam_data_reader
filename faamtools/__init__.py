from .version import __version__

__all__ = ['ObsData',
           'DsFld',
           'FaamFld'
          ]


class ObsData:
    """Generic class for storing several fields of observational data.

    Contains methods `__init__` and `__call__` for initialising and adding
    fields to this class.
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
    def __call__(self,**kwds):
        self.__dict__.update(kwds)


class DsFld:
    """A class for storing dropsonde data.

    Contains several attributes:
        raw: array-like, 'raw' data
        fil: array-like, 'filtered' data
        units: string, units
        long_name: string, name
    """
    def __init__(self, raw=None, fil=None, units='', long_name=''):
        self.raw = raw
        self.fil = fil
        self.units = units
        self.long_name = long_name


class FaamFld:
    """A class for storing core FAAM aircraft data.

    Contains several attributes:
        val: array-like, data values
        units: str, units
        long_name: str, name
    """
    def __init__(self, val=None, units='', long_name=''):
        self.val = val
        self.units = units
        self.long_name = long_name
