#!/usr/bin/env python
# coding: utf-8
# (1) Read in the master catalog (a SNANA HOSTLIB file) generated in Part I
# (2) For each simulated SN host galaxy, use the EAZY code to make a simulated host galaxy spectrum (from the best-fitting photoz template)
# (3) DOABLE, BUT NOT YET DONE: Store each simulated spectrum as an ascii .dat file with wavelength in nm and AB mag (suitable for input to the Subaru ETC).
# (4) STILL TBD :  Store the revised master catalog (now updated with SED .dat file names) as a modified SNANA HOSTLIB file (ascii text)


import os
import numpy as np
from matplotlib import pyplot as plt

from astropy.table import Table, Column
# from astropy.io import fits, ascii
#from astropy.coordinates import SkyCoord
#from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

import sncosmo

#from scipy.interpolate import interp1d
#from scipy.integrate import trapz


# TODO: Adjust the `datadir` variable to use a local path from the install dir
datadir = os.path.abspath("./data")

#Data_overview_table_name = os.path.join(datadir, "gal_lib_short1.dat") #name of table which contains galaxy ID information and location for DEIMOS and vUDS spectra
hostlib_filename = os.path.join(datadir,'cosmos_example_hostlib.txt') #name of SNANA HOSTLIB file (includes observed and rest-frame-synthetic photometry)
eazy_templates_filename = os.path.join(datadir,"Akari_Hosts_subset_SNR_v7.HOSTLIB")

vUDS_spec_location = os.path.join(datadir,"vUDS_spec/") 
DEIMOS_spec_location = os.path.join(datadir,'deimos_spec/') 
sim_spec_location = os.path.join(datadir,'sim_spec/') 

HST_table_name = os.path.join(datadir,'cosmos_match.dat') #name of table with HST ID information to locate and match spec files
HST_cosmos_folder_location = os.path.join(datadir,'COSMOS_3DHST_SPECTRA/') #location of all cosmos tile folders (i.e. directory which contains cosmos-02, cosmos-03, etc.)


medsmooth = lambda f,N : np.array(
    [np.median( f[max(0,i-N):min(len(f),max(0,i-N)+2*N)])
     for i in range(len(f))])

flcdm = FlatLambdaCDM(H0=73, Om0=0.27)


def load_eazypy_templates(eazytemplatefilename,
                          format='ascii.commented_header',
                          verbose=True,
                          **kwargs):
    """Read in the galaxy SED templates (basis functions for the
    eazypy SED fitting / simulation) and store as the 'eazytemplatedata'
    property.

    We read in an astropy Table object with N rows and M+1 columns, where
    N is the number of wavelength steps and M is the
    number of templates (we expect 13).
    The first column is the  wavelength array, common to all templates.

    We translate the Nx(M+1) Table data into a np structured array,
    then reshape as a (M+1)xN numpy ndarray, with the first row giving
    the wavelength array and each subsequent row giving a single
    template flux array.
    See the function simulate_eazy_sed_from_coeffs() to construct
    a simulated galaxy SED with a linear combination from this matrix.
    """
    eazytemplates = Table.read(eazytemplatefilename,
                               format=format, **kwargs)
    tempdata = eazytemplates.as_array()
    eazytemplatedata = tempdata.view(np.float64).reshape(
        tempdata.shape + (-1,)).T
    if verbose:
        print("Loaded Eazypy template SEDs from {0}".format(
            eazytemplatefilename))
    return eazytemplatedata


def scale_to_match_imag(wave, flam, imag, medsmooth_window=20):
    """KLUDGE!!  Using sncosmo to make this galaxy SED into a Source so 
    we can integrate into mags using the sncosmo bandmag, and rescale 
    to match a pre-defined mag
    
    wave: wavelength in angstroms
    flam: flambda in erg/s/cm2/A
    imag: sdss i band magnitude to scale to
    """    
    # check that we cover the i band
    if wave[0]>6600:
        wave = np.append([6580], wave)
        flam = np.append([1e-20], flam)
    if wave[-1]<8380:
        wave = np.append(wave, [8400])
        flam = np.append(flam, [1e-20])
    
    if medsmooth_window>1:
        # If a smoothing window size is given, use only the smoothed flux
        flam = medsmooth(flam, medsmooth_window)

    # Make a dummy sncosmo Source and scale it to the given sdss i band mag
    phase = np.array([-1, 0, 1, 2]) # need at least 4 phase positions for a source
    flux = np.array([flam, flam, flam, flam]) 
    galsource = sncosmo.TimeSeriesSource(phase, wave, flux)    
    galsource.set_peakmag(imag, 'sdssi', 'ab')

    fout = galsource.flux(0,wave)
    
    return(wave, fout)


def simulate_eazy_sed_from_coeffs(eazycoeffs, eazytemplatedata, z,
            returnfluxunit='flambda', returnwaveunit='A',
            limitwaverange=True, savetofile='', **outfile_kwargs):
    """
    Generate a simulated SED from a given set of input eazy-py coefficients
    and eazypy templates.

    NB: Requires the eazy-py package to apply the IGM absorption!
    (https://github.com/gbrammer/eazy-py)

    Optional Args:
    returnfluxunit: ['AB', 'flambda', 'fnu'] TODO: add Jy
        'AB'= return log(flux) as monochromatic AB magnitudes
        'AB25' = return AB mags, rescaled to a zeropoint of 25:  m=-2.5*log10(fnu)+25
        'flambda' = return flux density in erg/s/cm2/A
        'fnu' = return flux density in erg/s/cm2/Hz
    returnwaveunit: ['A' or 'nm'] limitwaverange: limit the output
    wavelengths to the range covered by PFS savetofile: filename for saving
    the output spectrum as a two-column ascii data file (suitable for use
    with the SubaruPFS ETC from C. Hirata.

    Returns
    -------
        wave : observed-frame wavelength, Angstroms or  nm
        flux : flux density of best-fit template, erg/s/cm2/A or AB mag
    """
    # the input data units are Angstroms for wavelength
    # and cgs for flux (flambda): erg s-1 cm-2 Ang-1
    wave_em = eazytemplatedata[0]  # rest-frame (emitted) wavelength
    wave_obs = wave_em * (1 + z)  # observer-frame wavelength
    obsfluxmatrix = eazytemplatedata[1:]
    flam = np.dot(eazycoeffs, obsfluxmatrix) # flux in erg/s/cm2/A
    
    if limitwaverange:
        # to simplify things, we only work with data over the Subaru PFS
        # + WFIRST prism wavelength range, from 200 to 2500 nm
        # (2000 to 25000 Angstroms)
        iuvoir = np.where((wave_obs>2000) & (wave_obs<25000))[0]
        wave_obs = wave_obs[iuvoir]
        wave_em = wave_em[iuvoir]
        flam = flam[iuvoir]
    
    # convert flux units to fnu using :  fnu=(lam^2/c)*flam  ;  c = 3.e18 A/s
    fnu = (wave_em * wave_em / 3.e18) * flam  # flux in erg/s/cm2/Hz
    
    # Confusing previous setup from GB, used to convert to AB mags w/ zpt=25
    #fnu_factor = 10 ** (-0.4 * (25 + 48.6))
    # flam_spec = 1. / (1 + z) ** 2
    # obsflux = sedsimflux * fnu_factor * flam_spec   
    
    try:
        import eazy.igm
        igmz = eazy.igm.Inoue14().full_IGM(z, wave_obs)
        fnu *= igmz
    except:
        pass
    
    if returnfluxunit=='AB':
        # convert from flux density fnu into monochromatic AB mag:
        returnflux = -2.5 * np.log10(fnu) - 48.6
    elif returnfluxunit=='AB25':
        # convert from flux density fnu into AB mags for zpt=25:
        returnflux = -2.5 * np.log10(fnu) + 25
    elif returnfluxunit=='fnu':
        returnflux = fnu
    elif returnfluxunit.startswith('flam'):
        returnflux = flam
    else:
        print("I don't recognize flux unit {}".format(returnfluxunit))
        return None,None
        
    if returnwaveunit=='nm':
        returnwave = wave_obs / 10.
    elif returnwaveunit.startswith('A'):
        returnwave = wave_obs
    else:
        print("I don't recognize wave unit {}".format(returnwaveunit))
        return None,None

    if savetofile:
        out_table = Table()
        outcol1 = Column(data=wave_obs, name='Angstroms')
        outcol2 = Column(data=flam, name='flambda')
        out_table.add_columns([outcol1, outcol2])
        out_table.write(savetofile, **outfile_kwargs)

    return returnwave, returnflux


def mAB_from_flambda(flambda, wave):
    """ Convert from flux density f_lambda in erg/s/cm2/A 
    into AB mag
    
    flambda: flux density f_lambda (erg/s/cm2/A)
    wave : wavelength in angstroms
    
    (see https://en.wikipedia.org/wiki/AB_magnitude)
    """
    return(-2.5 * np.log10(3.34e4 * wave * wave * (flambda / 3631)))



def plot_spec_comparison(galid, showphot=True, showvuds=True, showdeimos=True,
                         showhst=True, showeazy=True,
                         medsmooth_deimos=20, medsmooth_vuds=20,
                         medsmooth_hst=20,
                         rescaledeimos=True, rescalevuds=False, ax=None):
    """Plot flux vs wavelength for the given galaxy ID, showing the observed 
    the Eazy-simulated spectrum.

    TBD : also plot the simulated photometry from the Akari catalog.
    """
    if ax is None:
        fig = plt.figure(figsize=[12,4])
        ax = fig.add_subplot(1,1,1)

    # read in the eazy spectral templates data
    # NOTE: could do this without loading the whole hostlib as a SnanaSimData object, would just need to grab
    # the code from snhostspec 
    #sim1 = snhostspec.SnanaSimData()
    #sim1.load_hostlib_catalog("DATA/cosmos_example_hostlib.txt")
    #sim1.
    eazytemplatedata = load_eazypy_templates(eazy_templates_filename)

    # ---------------------------------
    # Simulated and Observed photometry :
    # --------------------------------


    # plot the EAZY simulated spectrum
    eazycoeffs = np.array([mastercat[col][ithisgal_mastercat]
                           for col in mastercat.colnames
                           if col.startswith('coeff_specbasis')])
    outfilename = "DATA/cosmos_example_spectra/cosmos_example_host_simspec_" +\
                  "{:6d}.fits".format(galid)
    wobs, mobs = simulate_eazy_sed_from_coeffs(
        eazycoeffs, eazytemplatedata, z,
        returnwaveunit='A', returnfluxunit='AB25',
        savetofile=outfilename, overwrite=True)
    if showeazy:
        ax.plot(wobs, mobs, label='EAZY SED fit', color='0.5', zorder=10)
    
    ax.set_xlim(3000,19000)
    #ax.set_ylim(-0.25*1e-16,0.3*1e-16)
    #ax.set_ylim(27, 20)
    ax.text(0.95,0.95, galid, ha='right', va='top', transform=ax.transAxes)
    ax.text(0.95,0.88, "z={0}".format(z), ha='right', va='top', transform=ax.transAxes)

    ax = plt.gca()
    ax.set_xlim(3000, 19000)
    ax.set_ylim(magmin-2,magmax+1)

    ax.legend(loc='upper left')
    ax.invert_yaxis()
    ax.grid()
    ax.set_xlabel('Observed Wavelength (Angstroms)')
    ax.set_ylabel("AB mag")
    plt.tight_layout()
    #plt.savefig("cosmos_example_spec_eazysims.pdf")

    return





