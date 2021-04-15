import unittest
from romanz import romanz

_TEST_HOSTLIB_ = "romanz/data/Akari_Hosts_subset_SNR_v7.HOSTLIB"
_TEST_GALCAT_ = "romanz/data/Akari_input_catalog_20210415.fits"
_TEST_SPECZCOMPLETENESSMAP_ = "romanz/data/specz_completeness_map_preromandata_cosmos.ecsv"

class TestCatalogBasedRedshiftSim(unittest.TestCase):
    """Test class for projecting redshift completeness from an input
    galaxy catalog.
    """

    #def setUpClass(cls) -> None:
    #    """ Called once, before all tests"""

    def setUp(self) -> None:
        """Called before each and every test"""
        self.romanz_sim = romanz.CatalogBasedRedshiftSim()
        self.romanz_sim.read_galaxy_catalog(_TEST_GALCAT_)
        return

    def test_read_catalog(self):
        """Read in a catalog of galaxy properties"""
        # TODO : convert to a unittest setup step?
        self.assertTrue(self.romanz_sim.galaxies is not None)
        self.assertTrue(len(self.romanz_sim.galaxies) > 0)

    def test_pick_host_galaxies(self):
        """Use a SN rate function to define the relative probability of each
        galaxy hosting a SN in any given year.
        """
        # Assign probabilities based on spec. SFR
        self.romanz_sim.assign_snhost_prob(
            snr_model='AH18PW', logmasscolname='logmass',
            logsfrcolname='logsfr')
        # Random sampling to pick host galaxies
        self.romanz_sim.pick_host_galaxies(
            100, snrcolname='snr_AH18_piecewise', replace=False)

    def test_apply_specz_completeness_map(self):
        """ Read in a 'map' for spectroscopic redshift completeness, which
        maps from one or more galaxy properties (mag, SFR, z...) onto a
        probability of getting a spec-z.

        Then apply this specz completeness map to a catalog
        of SN host galaxy properties (already read in) to define exactly which
        of those host galaxies gets a redshift.
        """
        interpolatordict = self.romanz_sim.apply_specz_completeness_map(
            _TEST_SPECZCOMPLETENESSMAP_,
            defining_columns_galcat=['g_LSST', 'r_LSST', 'i_LSST', 'z_LSST'],
            defining_columns_speczmap=['mag', 'mag', 'mag', 'mag'],
            efficiency_columns_speczmap=['efficiency_g', 'efficiency_r',
                                        'efficiency_i', 'efficiency_z'],
        )


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
