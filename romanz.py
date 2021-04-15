import numpy as np
from astropy import units as u
from astropy.table import Table
from scipy import interpolate as scinterp

class SNANAHostLib():
    """Class for parsing a SNANA HOSTLIB file.
    The file may contain a weight map, and must contain host galaxy data,
    following the standard SNANA HOSTLIB format.
    """
    def __init__(self, filename):
        """Read in a SNANA HOSTLIB file"""
        # find the 'VARNAMES' line, use it to define the start of the hostlib
        # section (as opposed to the wgtmap section)
        ngaldatastart = -1
        nwgtmapstart = -1
        iline = 0
        with open(filename, 'r') as read_obj:
            for line in read_obj:
                if len(line.strip().lstrip('#'))==0:
                    continue
                if line.strip().startswith('NVAR_WGTMAP:'):
                    wgtmaphdrline = line.split()
                    varnames_wgtmap = wgtmaphdrline[3:]
                if line.strip().startswith('WGT:') and nwgtmapstart<0:
                    nwgtmapstart = iline
                if line.strip().startswith('GAL:') and ngaldatastart<0:
                    ngaldatastart = iline
                if line.strip().startswith('VARNAMES:'):
                    ngaldataheader = iline
                iline += 1
        if nwgtmapstart >= 0:
            self.wgtmaptable = Table.read(
                filename, format='ascii.basic',
                names=['label']+varnames_wgtmap+['wgt','snmagshift'],
                data_start=nwgtmapstart-1,
                data_end=ngaldataheader-2,
                comment='#'
                )
        else:
            self.wgtmaptable = None

        galdatatable = Table.read(filename, format='ascii.basic',
                                  header_start=ngaldataheader-1,
                                  data_start=ngaldatastart-1, comment='#'
                                  )
        galdatatable.remove_columns(['VARNAMES:'])
        self.galdatatable = galdatatable
        return


class CatalogBasedRedshiftSim():
    """Class for projecting redshift completeness from an input
    galaxy catalog.
    """

    def __init__(self):
        self.postsurvey = False
        self.galaxies = None

    def read_galaxy_catalog(self, filename):
        """Read in a catalog of galaxy properties

        Parameters
        ----------
        filename : str
          full path to the file containing galaxy properties (e.g. Mass, SFR,
          magnitudes, etc.).  May be a SNANA HOSTLIB file, or any formtat that
          can be auto-parsed by astropy.table.Table.read()
        """
        # check if it is a hostlib file


        self.galaxies = Table.read(filename)
        return

    def pick_host_galaxies(self):
        """Use a SN rate function to define a random
        sampling of the galaxies that are SN hosts.

        Alternatively, read in a SNANA output file (.dump file maybe?) that
        has already run a survey simulation and picked host galaxies.
        """
        pass

    def apply_specz_completeness_map(self, filename,
                                     defining_columns_galcat,
                                     defining_columns_speczmap,
                                     efficiency_column_speczmap,
                                     fill_value = np.nan
                                     ):
        """Read in a 'map' for spectroscopic redshift completeness, which
        maps from one or more galaxy properties (mag, SFR, z...) onto a
        probability of getting a spec-z.

        Preferred format of the input file is a .ecsv file, but anything
        that astropy.table can read is OK in principle.

        Then apply the specz completeness map to the catalog
        of host galaxy properties (already read in) to define exactly which
        of the galaxies gets a redshift.

        If the method 'pick_host_galaxies' has already been run
        (so the flag postsurvey == True), then only galaxies defined as SN
        hosts are assigned a redshift.

        Parameters
        ----------
        filename : str
           path to astropy-readable file

        defining_columns_galcat : listlike
           list of strings specifying the column names in the galaxy catalog
           (self.galaxies) for parameters that are used to define the specz
           efficiency (e.g. if this is a SFR-based specz map then this may
           be ['logSFR'])

        defining_columns_speczmap : listlike
           list of strings specifying the corresponding column names in the
           specz map file (given by 'filename').  Must be the same length as
           defining_columns_galcat, giving corresponding column names in the
           same order.

        efficiency_column_speczmap : str
           Name of the column that contains the specz efficiency (or
           completeness fraction) for each row in the specz map file.
        """
        if len(defining_columns_galcat)!=len(defining_columns_speczmap):
            raise RuntimeError(
                'You must specify the same number of columns from the '
                'galaxy catalog and the specz efficiency catalog.'
            )

        speczmap = Table.read(filename)

        # build an interpolating function from the specz efficiency map
        points = np.array(
            [speczmap[colname] for colname in defining_columns_speczmap])
        points = points.reshape(
            (len(defining_columns_speczmap), len(speczmap)))

        # TODO: reduce dimensionality of indata when possible
        values = speczmap[efficiency_column_speczmap]

        # TODO? : choose interpolation function based on dimensionality?
        interpolator = scinterp.LinearNDInterpolator(
            points, values, fill_value=fill_value)

        return(interpolator)

    def make_photoz_accuracy_map(self):
        """For every galaxy in the catalog of galaxy properties, apply a
        photo-z function that defines the 'measured' photo-z value and
        uncertainty (photoz pdf).  Includes catastrophic outliers.
        """
        pass

    def report_redshift_completeness(self):
        """Produce a report of the overall redshift completeness, accuracy
        and precision, based on multiple spectroscopic 'filters' and the
        random assignment of photo-z values.
        """
        pass

