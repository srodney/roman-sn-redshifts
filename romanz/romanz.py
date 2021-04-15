import numpy as np
from astropy import units as u
from astropy import table
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
        nwgtmapstart = -1
        ngaldataheader = -1
        ngaldatastart = -1
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
        if ngaldataheader < 0:
            raise RuntimeError(r"{filename} is not an SNANA HOSTLIB file")

        if nwgtmapstart >= 0:
            self.wgtmaptable = table.Table.read(
                filename, format='ascii.basic',
                names=['label']+varnames_wgtmap+['wgt','snmagshift'],
                data_start=nwgtmapstart-1,
                data_end=ngaldataheader-2,
                comment='#'
                )
        else:
            self.wgtmaptable = None

        galdatatable = table.Table.read(filename, format='ascii.basic',
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
        # TODO: check if it is a hostlib without try/except
        try :
            self.galaxies = table.Table.read(filename)
        except:
            try:
                hostlib = SNANAHostLib(filename)
                self.galaxies = hostlib.galdatatable
            except:
                raise RuntimeError(
                    f"Can't read in {filename}. "
                    "It may not be a valid hostlib or astropy-readable table.")
        return


    def assign_snhost_prob(self, snr_model='AH18S',
                           logmasscolname='logmass',
                           logsfrcolname='logsfr',
                           verbose=True):
        """Add a column to the 'galaxies' catalog that gives the relative
        probability for each galaxy hosting a SN in any given observer-frame
        year.  This is computed based on the predicted SN rate (number of SN
        explosions per observer-frame year) of each galaxy, adopting the
        specified SN rate model.

        Parameters
        ----------
        snr_model : str
           'A+B' : SNR = A*M + B*SFR   (Scannapieco & Bildsten 2005)
           'AH18S' : the smooth logarithmic sSFR model (Andersen & Hjorth 2018)
           'AH18PW' : the piecewise sSFR model (Andersen & Hjorth 2018)

        logmasscolname : str
           name of column in the galaxies Table containing the log10(Mass)

        logsfrcolname : str
           name of column in the galaxies Table containing the
           log10(StarFormationRate)

        verbose : bool
            Set to True to print messages.
        """
        if self.galaxies is None:
            print("No 'galaxies' catalog loaded. Use 'read_galaxy_catalog()'")

        if snr_model.lower()=='a+b':
            # Note: adopting the A and B values from Andersen & Hjorth 2018
            # but dividing by 1e-4 (so the SNR below actually counts the number
            # of SN explodiing per 10000 yrs)
            A = 4.66 * 1e-10
            B = 4.88
            snr = A * 10 ** self.galaxies[logmasscolname] + B * 10 ** self.galaxies[logsfrcolname]
            # divide by the total snr to get relative probabilities
            snr /= np.nanmax(snr)
            snrcolname = 'snr_A+B'
            snrcol = table.Column(data=snr, name='snr_A+B')
        elif snr_model.lower() == 'ah18s':
            logssfr = self.galaxies[logsfrcolname] - self.galaxies[logmasscolname]
            ssnr = ssnr_ah18_smooth(logssfr)
            snr = ssnr * 10 ** self.galaxies[logmasscolname]
            snr /= np.nanmax(snr)
            snrcolname = 'snr_AH18_smooth'
            snrcol = table.Column(data=snr, name=snrcolname)
        elif snr_model.lower() == 'ah18pw':
            logssfr = self.galaxies[logsfrcolname] - self.galaxies[logmasscolname]
            ssnr = ssnr_ah18_piecewise(logssfr)
            snr = ssnr * 10 ** self.galaxies[logmasscolname]
            snr /= np.nanmax(snr)
            snrcolname = 'snr_AH18_piecewise'
        else:
            raise RuntimeError(r"{snr_model} is not a know SN rate model.")

        snrcol = table.Column(data=snr, name=snrcolname)
        if snrcolname in self.galaxies.colnames:
            self.galaxies[snrcolname] = snr
        else:
            self.galaxies.add_column(snrcol)
        if verbose:
            print(f"Added/updated relative SN rate column using {snr_model} model")
        return


    def pick_host_galaxies(self, nsn, snrcolname='snr_AH18_piecewise',
                           replace=False, verbose=True):
        """Do a random draw to assign 'nsn' supernovae to galaxies in the
        galaxies catalog, based on the (pre-defined) relative SN rates.

        TODO: (Alternatively, read in a SNANA output file (.dump file maybe?)
        that has already run a survey simulation and picked host galaxies.)

        Parameters
        ----------
        replace
        nsn : int
          number of SN to assign to host galaxies

        snrcolname : str
           name of the column in the galaxies catalog that gives the relative
           SN rate (or 'weight') for each galaxy.  This may be created by the
           assign_snhost_prob() method.

        replace : bool
           Whether to sample with replacement.  If True, a galaxy may host
           more than one SN. If False, then assign no more than one SN to
           each galaxy (requires nsn<len(galaxies))
        """
        if ~replace and nsn > len(self.galaxies):
            raise RuntimeError(
                r'Picking hosts without replacement, but Nsn > len(galaxies)')

        # Pick SN host galaxies
        galindices = np.arange(len(self.galaxies))
        psnhost = self.galaxies[snrcolname]/np.sum(self.galaxies[snrcolname])
        snindices = np.random.choice(
            galindices, nsn, replace=replace, p=psnhost)

        # Add a boolean 'host' column to the galaxies catalog
        ishost = np.zeros(len(self.galaxies), dtype=bool)
        ishost[snindices] = True
        hostcol = table.Column(name='host', data=ishost)
        if 'host' in self.galaxies.colnames:
            self.galaxies['host'] = hostcol
        else:
            self.galaxies.add_column(hostcol, index=1)

        # TODO: Alternate approach:  read in a SNANA output file (.dump file
        #  maybe?) that has already run a survey simulation and picked hosts.

        if verbose:
            print(f"Assigned {nsn} SNe to hosts using {snrcolname} probabilities.")
        return



    def apply_specz_completeness_map(self, filename,
                                     defining_columns_galcat,
                                     defining_columns_speczmap,
                                     efficiency_columns_speczmap,
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

        defining_columns_speczmap : listlike, same length as above
           list of strings specifying the corresponding column names in the
           specz map file (given by 'filename').  Must be the same length as
           defining_columns_galcat, giving corresponding column names in the
           same order.

        efficiency_columns_speczmap : listlike, same length as above
           list of column names giving the specz
           efficiency (or completeness fraction) for each row in the specz
           map file.
        """
        if (len(defining_columns_galcat)!=len(defining_columns_speczmap) or
            len(defining_columns_galcat)!=len(efficiency_columns_speczmap)):
            raise RuntimeError(
                'You must specify the same number of columns from the '
                'galaxy catalog and the specz efficiency catalog.')

        # TODO : make a masked array to remove NaNs ? ?
        speczmap = table.Table.read(filename)

        # TODO : build a separate interpolating function for each of
        #  the possible input parameters ?
        interpolatordict = {}
        for i in range(len(defining_columns_galcat)):
            colname_galcat = defining_columns_galcat[i]
            xobs = self.galaxies[colname_galcat]
            colname_param = defining_columns_speczmap[i]
            x = speczmap[colname_param]
            colname_efficiency = efficiency_columns_speczmap[i]
            y = speczmap[colname_efficiency]
            interpolator = scinterp.interp1d(
                x, y, bounds_error=False, fill_value=fill_value)
            interpolatordict[colname_galcat] = interpolator

        return(interpolatordict)

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



def ssnr_ah18_smooth(logssfr):
    """ Returns the Type Ia specific SN rate per Tyr
    (number of SN Ia exploding per 10^12 yr per solar mass)
    for a galaxy, using the model of Andersen & Hjorth 2018, which is based
    on the specific star formation rate, given as log10(SSFR).
    """
    a = (1.5)*1e-13 # (1.12)*1e-13
    b = 0.5 # 0.73
    k = 0.4 # 0.49
    ssfr0 = 1.7e-10# 1.665e-10
    # logssfr0 = -9.778585762157661    # log10(ssfr0)
    ssfr = np.power(10.,logssfr)
    ssnr = (a + (a/k) * np.log10(ssfr/ssfr0 + b)) * 1e12
    #ssnr = np.max(ssnr, 0.7)
    return(ssnr)


def ssnr_ah18_piecewise(logssfr):
    """ Returns the Type Ia specific SN rate per Tyr
    (number of SN Ia exploding per 10^12 yr per solar mass)
    for a galaxy, using the piecwise linear model
    of Andersen & Hjorth 2018, which is based
    on the specific star formation rate, given as log10(SSFR).
    """
    # Note that the alpha scaling parameter
    # has been multiplied by 1e12 to get units of Tyr-1
    alpha = (1.12)* 1e5
    beta = 0.586
    ssfr2 = 1.01e-11
    ssfr1 = 1.04e-9

    S1 = np.power(ssfr1, beta)
    S2 = np.power(ssfr2, beta)

    if not np.iterable(logssfr):
        logssfr = np.array([logssfr])

    ssfr = np.power(10.,logssfr)

    ilow = np.where(ssfr<=ssfr2)[0]
    imid = np.where((ssfr>ssfr2) & (ssfr<ssfr1))[0]
    ihi = np.where(ssfr>=ssfr1)[0]

    ssnrmid = alpha * np.power(ssfr[imid], beta)

    ssnr = alpha * np.where(ssfr<=ssfr2, S2,
                            np.where(ssfr>=ssfr1, S1,
                                     np.power(ssfr, beta)))
    if len(ssnr)==1:
        ssnr = ssnr[0]
    return(ssnr)

