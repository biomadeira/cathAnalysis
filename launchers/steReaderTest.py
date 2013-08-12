import unittest
from steReader import STE
class steReaderTestCase(unittest.TestCase):
    def testLoader(self):
        fileName = 'state.ste'
        ste = STE(fileName)
        assert  ste.runMmult is False