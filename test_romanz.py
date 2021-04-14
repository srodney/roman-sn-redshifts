import unittest
from . import romanz

class TestCatalogBasedRedshiftSim(unittest.TestCase):
    """Test class for projecting redshift completeness from an input
    galaxy catalog.
    """

    def test_read_catalog(self):
        """Read in a catalog of galaxy properties"""
        # TODO : convert to a unittest setup step?
        self.assertEqual(True, False)

    def test_read_specz_completeness_map(self):
        """Read in a 'map' for spectroscopic redshift completeness, which
        maps from one or more galaxy properties (mag, SFR, z...) onto a
        probability of getting a spec-z.
        """
        # TODO : convert to a unittest setup step?
        self.assertEqual(True, False)

    def test_apply_specz_completeness_map(self):
        """Apply a specz completeness map (already read in) to a catalog
        of host galaxy properties (already read in) to define a random
        sampling of the galaxies that are SN hosts, and define exactly which
        of those host galaxies gets a redshift.
        """
        self.assertEqual(True, False)

    def test_make_photoz_accuracy_map(self):
        """For every galaxy in the catalog of galaxy properties, apply a
        photo-z function that defines the 'measured' photo-z value and
        uncertainty (photoz pdf).  Includes catastrophic outliers.
        """
        self.assertEqual(True, False)

    def test_report_redshift_completeness(self):
        """Produce a report of the overall redshift completeness, accuracy
        and precision, based on multiple spectroscopic 'filters' and the
        random assignment of photo-z values.
        """
        self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
