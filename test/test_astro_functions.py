import unittest
import doctest
import os

import astrolyze.functions.astro_functions

class Test(unittest.TestCase):
    """Unit tests for astro_functions."""

    def test_doctests(self):
        """Run astro_functions doctests"""
        doctest.testmod(astrolyze.functions.astro_functions)
        os.system('/usr/bin/convert black_body.eps testfigures/black_body.jpg')
        os.system('cp testfigures/black_body.jpg ../doc/figures/')
        os.system('rm black_body.eps')
if __name__ == "__main__":
    unittest.main()
