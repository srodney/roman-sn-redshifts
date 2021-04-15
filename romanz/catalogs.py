# /usr/bin/env python
# 2017.03.10  S. Rodney
# Reading in host galaxy data from WFIRST simulations
# and computing an estimated exposure time for Subaru+PFS
# to get a redshift for that host.

from numbers import Number
from astropy.io import fits
from astropy.table import Table, Column, Row, MaskedColumn
from astropy import table
from astropy.io import ascii
import os
import subprocess
#import exceptions
import numpy as np
from glob import glob
from matplotlib import pyplot as plt
from io import StringIO
import time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck15 as cosmo

__3DHST_DATADIR__ =  '/Users/rodney/src/wfirst/3DHST/'
__CANDELS_DATADIR__ = '/Users/rodney/src/wfirst/CANDELS/'

__CANDELS_FIELDS_LIST__ = ['AEGIS', 'UDS', 'COSMOS', 'GOODS-N', 'GOODS-S']

__CANDELS_COLUMN_KEEPLIST__ = [
    'NUMBER', 'X_IMAGE', 'Y_IMAGE', 'CXX_IMAGE', 'CYY_IMAGE', 'CXY_IMAGE',
    'THETA_IMAGE', 'THETA_WORLD', 'ELLIPTICITY', 'KRON_RADIUS',
    'A_IMAGE', 'B_IMAGE', 'ALPHA_J2000','DELTA_J2000',
    'BACKGROUND','FLUX_BEST', 'FLUXERR_BEST', 'MAG_BEST', 'MAGERR_BEST',
    'FLUX_RADIUS', 'ISOAREA_IMAGE', 'FWHM_IMAGE', 'FLAGS', 'CLASS_STAR',
    'NITER_GALFIT', 'NEIGH_GALFIT', 'CHISQ_GALFIT', 'NFREE_GALFIT',
    'NFIX_GALFIT', 'NDOF_GALFIT', 'CHI2NU_GALFIT', 'FLAG_GALFIT',
    'RE_GALFIT_BAND','REERR_GALFIT_BAND',
    'N_GALFIT_BAND','NERR_GALFIT_BAND',
    'Q_GALFIT_BAND','QERR_GALFIT_BAND',
    'PA_GALFIT_BAND','PAERR_GALFIT_BAND'
     ]



class Supernova(object):
    """ A single SN object.  Observed in an HST survey
    """
    def __init__(self):
        self.id = None
        self.ra = 0 * u.degree
        self.dec = 0 * u.degree
        self.z = -9

    def load_from_catalog(self, catalogfile, snid):
        catalog = ascii.read(catalogfile)
        idx = np.where(catalog['id'] == snid)[0][0]
        self.load_from_catalog_row(catalog[idx])


    def load_from_catalog_row(self, rowdata):
        """load the data for a single SN from a single row of a catalog (or
        any dict object).

        :param rowdata: dict containing required keys id, RA, DEC, field, z.
           Any other keys are optional.
        :return:
        """
        # TODO : make this more general ?
        required_keys = ['id', 'z', 'field', 'ra', 'dec']
        if isinstance(rowdata, table.Row):
            colnames = rowdata.colnames
        elif isinstance(rowdata, table.Table):
            rowdata = rowdata[0]
            colnames = rowdata.colnames
        elif isinstance(rowdata, dict):
            colnames = list(rowdata.keys())
        for key in required_keys:
            assert key in colnames

        self.ra =  rowdata['ra'] * u.degree
        self.dec =  rowdata['dec'] * u.degree
        self.z =  rowdata['z']
        self.field = rowdata['field']
        self.id = rowdata['id']
        self.location_3D = SkyCoord(
            ra = self.ra, dec = self.dec,
            distance = cosmo.luminosity_distance(self.z))
        self.location_sky = SkyCoord(ra = self.ra, dec = self.dec)

        for key in colnames:
            self.__dict__[key] = rowdata[key]

    @property
    def summary(self):
        summary_string = ""
        for key in ['id', 'nickname', 'field']:
            id_row = '{:s} = {:s}\n'.format(
                key.upper(), self.__dict__[key])
            summary_string += id_row

        if 'zerrplus' in list(self.__dict__.keys()):
            zrow = 'z = {:.3f} +{:.3f} {:.3f}\n'.format(
                self.z, self.zerrplus, self.zerrminus)
        else:
            zrow = 'z = {:.3f}\n'.format(self.z)
        summary_string += zrow

        if self.z > 0:
            location_row = self.location_3D.__str__()
        else:
            location_row = self.location_sky.__str__()
        summary_string += location_row

        return summary_string


    def print_summary(self):
        print((self.summary))


    def find_nearest_galaxy(self, galaxy_catalog):
        """ Find galaxies in the given catalog that are potential host galaxies
        for this SN.  Returns a Table object in which the first row gives all
        available information for the galaxy with the smallest angular
        separation from the SN.   A second row will report the galaxy with the
        smallest physical distance, whenever it is possible to calculate (i.e.,
        when the catalog has redshift information).

        :param galaxy_catalog: a catalog with locations of galaxies (3D
        locations with
          distances defined, or 2D locations, as projected on the sky)
        :return: Table object with one row for each nearest galaxy.
        """
        idx_skymatch, d2d_skymatch, d3d_skymatch = \
            self.location_sky.match_to_catalog_sky(
                galaxy_catalog.locations)

        match_row = galaxy_catalog.catalog_data[idx_skymatch]
        output_table = Table(rows=match_row)
        d2d_data = [d2d_skymatch]
        d3d_data = [d3d_skymatch]

        if 'z' in list(galaxy_catalog.__dict__.keys()):
            idx_3dmatch, d2d_3dmatch, d3d_3dmatch = \
                self.location_3D.match_to_catalog_3d(
                    galaxy_catalog.locations)

            if idx_3dmatch != idx_skymatch:
                match_row_3d = galaxy_catalog.catalog_data[idx_3dmatch]
                output_table.add_row(match_row_3d)
                d2d_data.append(d2d_3dmatch)
                d3d_data.append(d3d_3dmatch)
        output_table.add_column(Column(data=d2d_data, name='sep_to_sn_arcsec',
                                       unit='degrees'))
        output_table.add_column(Column(data=d3d_data, name='dist_to_sn_Mpc',
                                       unit='Mpc'))
        return(output_table)


    def find_nearby_galaxies(self, galaxy_catalog, search_radius_arcsec=0,
                             search_radius_Mpc=0, zmin=0, zmax=5):
        """ Find galaxies in the given catalog that are near to the location
        of this SN.  The user may provide a search radius in arcseconds for
        finding all galaxies with a projected separation on the sky within
        that radius.  Alternatively, or in addtion, the user may provide a
        search radius in Mpc to find galaxies within a physical separation
        distance.  This latter option only works if the galaxy_catalog
        includes distance information (i.e., from redshifts).

        Returns a Table object in which each row gives all available
        information for one of the galaxies that satisfies either of the
        user-specified separation criteria.

        :param galaxy_catalog: a catalog with locations of galaxies (3D
        locations with distances defined, or 2D locations, as projected on
        the sky)
        :param search_radius_arcsec: angular search radius in arcsec
        :param search_radius_Mpc: physical search radius in arcsec
        :param zmin, zmax: limit the catalog search to only galaxies within
         this redshift range
        :return: Table object with one row for each nearby galaxy.
        """
        assert search_radius_arcsec + search_radius_Mpc > 0

        if zmax>0 and 'z' in galaxy_catalog.catalog_data.colnames:
            izrange = np.where((zmin < galaxy_catalog.catalog_data['z']) &
                               (galaxy_catalog.catalog_data['z']< zmax))[0]
            galcatdata = galaxy_catalog.catalog_data[izrange]
            galcatloc = galaxy_catalog.locations[izrange]
        else:
            galcatdata = galaxy_catalog.catalog_data
            galcatloc = galaxy_catalog.locations

        output_table_angular = None
        output_table_physical = None

        if search_radius_arcsec > 0:
            d2d = self.location_sky.separation(galcatloc)
            catalogmsk = d2d < search_radius_arcsec * u.arcsec
            output_table_angular = galcatdata[catalogmsk]
            output_table_angular.add_column(
                Column(data=d2d[catalogmsk].to(u.arcsec),
                       name='angular_distance_to_sn',
                       unit=u.arcsec))
            if search_radius_Mpc <= 0:
                return(output_table_angular)

        if search_radius_Mpc > 0:
            d3d = self.location_3D.separation_3d(galcatloc)
            catalog3dmsk = d3d < search_radius_Mpc * u.Mpc
            output_table_physical = galcatdata[catalog3dmsk]
            output_table_physical.add_column(
                Column(data=d3d[catalog3dmsk].to(u.Mpc),
                       name='physical_distance_to_sn',
                       unit=u.Mpc))
            if search_radius_arcsec <= 0:
                return(output_table_physical)

        return(output_table_angular, output_table_physical)




class GalaxyCatalog(object):
    """ A generic catalog of galaxies, constructed from any one of the
    many CANDELS/CLASH catalogs.

    Necessary properties:
    ra, dec, field, catalogfile, idx_catalog

    Optional:
    z, zerr, photometry, mass, luminosity, restframe_photometry

    """
    def find_nearest_galaxy(self, location, tolerance_arcsec=0):
        """ Find the galaxy in this catalog that is closest to the given
        position.  Returns a Table object with a single row, giving all
        available information for the galaxy with the smallest angular
        separation from the specified location.

        :param location: a SkyCoord coordinate location to search near
        :param tolerance_arcsec: return None if there are no matches within
           this radius, given in arcseconds
        :return: Table object with one row.
        """
        idx_match, angsep_match, distance_skymatch = \
            location.match_to_catalog_sky(self.locations)
        if isinstance(idx_match, np.ndarray):
            idx_match = int(idx_match)

        if ((tolerance_arcsec > 0) and
            (angsep_match.to(u.arcsec) > tolerance_arcsec * u.arcsec)):
            return None

        match_row = self.catalog_data[idx_match]
        output_table = Table(rows=match_row)
        # angsep_column = Column(data=[angsep_match],
        #                        name='angular_distance',
        #                        unit='degrees')
        # output_table.add_column(angsep_column)
        return(output_table)



class SNCatalog(object):
    """A catalog of SN observed with HST"""

    def __init__(self, txtfilename):
        self.catalog_3dhst = {}
        self.catalog_candels = {}
        self.catalog_clash = {}

        self.catalog_sn = ascii.read(txtfilename)
        self.supernova_list = []
        for rowdata in self.catalog_sn:
            sn = Supernova()
            sn.load_from_catalog_row(rowdata)
            self.supernova_list.append(sn)

    @property
    def fieldlist(self):
        return(np.array([sn.field for sn in self.supernova_list]))

    @property
    def idlist(self):
        return(np.array([sn.id for sn in self.supernova_list]))

    @property
    def ralist(self):
        return(np.array([sn.ra for sn in self.supernova_list]))

    @property
    def declist(self):
        return(np.array([sn.dec for sn in self.supernova_list]))

    @property
    def zlist(self):
        return(np.array([sn.z for sn in self.supernova_list]))


    def load_galaxy_catalogs(self, field=None):
        if field is None:
            fieldlist = np.unique(self.fieldlist)
        elif isinstance(field, str):
            fieldlist = [field]
        else:
            fieldlist = field

        for field in fieldlist:
            if field in ['EGS', 'AEGIS', 'COSMOS', 'UDS',
                             'GOODS-S', 'GOODS-N']:
                if field in list(self.catalog_3dhst.keys()):
                    continue

                self.catalog_3dhst[field] = Catalog3DHST(field)
                self.catalog_candels[field] = CatalogCANDELS(field)
            else:
                if field in list(self.catalog_clash.keys()):
                    continue
                self.catalog_clash[field] = CatalogCLASH(field)


    def find_nearby_galaxies(self, search_radius_arcsec=10.0,
                             galcatmatch_tolerance_arcsec=0.3):
        """ find nearby galaxies with similar redshift.  make a table with the
        SN ID and a unique "host candidate" identifier
        :return:
        """
        if len(list(self.catalog_3dhst.keys()))<5:
            self.load_galaxy_catalogs()

        output_table_list = []
        for sn in self.supernova_list:
            if sn.field not in __CANDELS_FIELDS_LIST__:
                # TODO: get CLASH catalogs including parallel fields
                continue

            # make a table of nearby galaxies from the latest CANDELS master
            # photometry catalogs (via Boris Hauessler)
            catphot = self.catalog_candels[sn.field]

            # For now, no z restrictions applied
            zmin = 0 # sn.z - max(3*np.abs(sn.zerrminus), 0.1)
            zmax = 0 # sn.z + max(3*np.abs(sn.zerrplus), 0.1)
            nearby_galaxies_table = sn.find_nearby_galaxies(
                catphot, search_radius_arcsec=search_radius_arcsec,
                search_radius_Mpc=0, zmin=zmin, zmax=zmax)

            # Add a unique ID for every nearby galaxy in the table
            nearby_galaxy_id_list = [
                '{:s}-gal{:02d}'.format(sn.id, ihostcand)
                for ihostcand in range(len(nearby_galaxies_table))]
            nearby_galaxies_table.add_column(
                Column(data=nearby_galaxy_id_list, name='id_gal',
                       dtype=str))

            # Extend the table of nearby galaxies to include masked
            # columns for all the data from the 3DHST catalog.  The values are
            # initialized to None in every cell, and masked.
            cat3d = self.catalog_3dhst[sn.field]
            n_nearby = len(nearby_galaxies_table)
            for cat3dcolname in cat3d.catalog_data.colnames:
                cat3dcolumn = cat3d.catalog_data[cat3dcolname]
                if isinstance(cat3dcolumn[0], Number):
                    fillval = 0
                else:
                    fillval = ''
                if cat3dcolname in nearby_galaxies_table.colnames:
                    cat3dcolname += "_3dhst"
                newphotcol = MaskedColumn(
                    data = [fillval for i in range(n_nearby)],
                    name = cat3dcolname, dtype=cat3dcolumn.dtype,
                    mask = [True for i in range(n_nearby)],
                    fill_value = fillval
                )
                nearby_galaxies_table.add_column(newphotcol)

            # For each nearby galaxy found in the CANDELS photometry catalog,
            # seek it out in the 3DHST catalog and add any matching 3DHST data
            # to our host galaxy candidate info table
            for irow in range(len(nearby_galaxies_table)):
                row = nearby_galaxies_table[irow]
                ra_gal = row['ALPHA_J2000']
                dec_gal = row['DELTA_J2000']
                loc_gal = SkyCoord(ra=ra_gal * u.deg, dec=dec_gal * u.deg)
                nearest_galaxy_3dhst_table = cat3d.find_nearest_galaxy(
                    loc_gal, tolerance_arcsec=galcatmatch_tolerance_arcsec)
                if nearest_galaxy_3dhst_table is None:
                    continue
                for colname in nearest_galaxy_3dhst_table.colnames:
                    nearby_galaxies_table[colname][irow] = \
                        nearest_galaxy_3dhst_table.columns[colname][0]
                    nearby_galaxies_table[colname].mask[irow] = False

            sn_data_table = Table()
            isn = np.where(self.catalog_sn['id'] == sn.id)[0][0]
            for i in range(len(nearby_galaxies_table)):
                if i==0:
                    sn_data_table = self.catalog_sn[[isn]]
                else:
                    sn_data_table.add_row(self.catalog_sn[isn])
            for colname in sn_data_table.colnames:
                sn_data_table[colname].name = sn_data_table[colname].name + '_sn'

            output_table = table.hstack(
                [sn_data_table, nearby_galaxies_table])
            output_table_list.append(output_table)
        output_table = table.vstack(output_table_list)
        return output_table


class CatalogCLASH(GalaxyCatalog):
    """ Galaxy photometry from the CLASH team"""

    def __init__(self, fieldname, verbose=False):
        self.fieldname = fieldname


class Catalog3DHST(GalaxyCatalog):
    """ Galaxy photometry, redshifts and derived parameters 
    from 3DHST for a single field
    """

    def __init__(self, fieldname, verbose=False):
        if fieldname == 'EGS':
            fieldname = 'AEGIS'

        self.fieldname = fieldname.upper()
        self.fieldnameshort = fieldname.replace('-','').lower()
        self.verbose = verbose

        zoutglob = os.path.join(__3DHST_DATADIR__,
                                '{:s}_3dhst.*.zout.fits'.format(
                                    self.fieldnameshort))
        filelist = glob(zoutglob)
        if not len(filelist):
            raise RuntimeError(
                'Catalog of 3DHST data for' +
                ' {:s} does not exists'.format(self.fieldname))
        self.catalog_data = Table.read(filelist[0])
        zerrminus = Column( data=self.catalog_data['z160'].data -
                                 self.catalog_data['z_phot'].data,
                            name='zerrminus')
        zerrplus = Column( data=self.catalog_data['z840'].data -
                                self.catalog_data['z_phot'].data,
                            name='zerrplus')
        z = Column(data=np.where(self.catalog_data['z_spec'] > 0,
                                 self.catalog_data['z_spec'],
                                 self.catalog_data['z_phot']),
                   name='z')

        # read in the 3DHST master photometry catalog to get RA and Dec
        masterphotfile = os.path.join(__3DHST_DATADIR__,
                                      '3dhst_master.phot.v4.1.cat.FITS')
        masterphottable = Table.read(masterphotfile)
        ithisfield = np.where(masterphottable['field'] == self.fieldname)[0]
        photkeylist = [k for k in list(masterphottable.keys())
                       if k.startswith('f_') or k.startswith('e_')]
        self.photometry = {}
        for photkey in photkeylist:
            self.photometry[photkey] = masterphottable[photkey][ithisfield]

        idmaster = masterphottable['id'][ithisfield]
        isort = np.searchsorted(self.catalog_data['id'], idmaster)
        ra = masterphottable['ra'][ithisfield][isort]
        dec = masterphottable['dec'][ithisfield][isort]

        self.locations = SkyCoord(
            ra = ra * u.degree, dec = dec * u.degree,
            distance = cosmo.luminosity_distance(z))
        self.catalog_data.add_columns([z, zerrplus, zerrminus, ra, dec])

        for colname in self.catalog_data.colnames:
            if not colname.endswith('_3dhst'):
                self.catalog_data[colname].name = colname + '_3dhst'


    def load_eazy_fitparam(self):
        """ Load the metadata from the EAZY fits """
        eazydatglob = os.path.join(__3DHST_DATADIR__,
                                   '{:s}_3dhst.*.eazypy.data.fits'.format(
                                       self.fieldnameshort))
        filelist = glob(eazydatglob)
        if not len(filelist):
            raise exceptions.RuntimeError(
                'Catalog of EAZY fits to 3DHST' +
                ' for {:s} does not exists'.format(self.fieldname))
        eazydatfile = filelist[0]
        hdulist = fits.open(eazydatfile)
        self.eazyfitdict = dict([(hdu.name, hdu.data) for hdu in hdulist])
        hdulist.close()




class CatalogCANDELS(GalaxyCatalog):

    def __init__(self, fieldname, verbose=False):
        if fieldname == 'EGS':
            fieldname = 'AEGIS'

        self.fieldname = fieldname.upper()
        self.fieldnameshort = fieldname.replace('-','').lower()
        self.verbose = verbose
        if fieldname == 'AEGIS':
            self.fieldnameshort = 'EGS'

        candelsdatglob = os.path.join(__CANDELS_DATADIR__,
                                      'CANDELS_{:s}*.fits'.format(
                                          self.fieldname.upper()))
        filelist = glob(candelsdatglob)
        if not len(filelist):
            raise RuntimeError(
                'Catalog of CANDELS data for' +
                ' {:s} does not exist'.format(self.fieldname))
        candelsdatfile = filelist[0]
        self.catalog_data = Table.read(candelsdatfile)
        for colname in self.catalog_data.colnames:
            if (colname not in __CANDELS_COLUMN_KEEPLIST__):
                self.catalog_data.remove_column(colname)
                continue
            if (len(self.catalog_data.columns[colname].shape) > 1):
                column_single = self.catalog_data.columns[colname][:,-1]
                self.catalog_data.replace_column(colname, column_single)
                if 'BAND' in colname:
                    modname = colname.replace('BAND','F160W')
                self.catalog_data.columns[colname].name = modname
        fieldcolumn = Column(
            data=[fieldname for i in range(len(self.catalog_data))],
            name='FIELD')
        self.catalog_data.add_column(fieldcolumn)
        self.ra = self.catalog_data['ALPHA_J2000']
        self.dec = self.catalog_data['DELTA_J2000']
        self.locations = SkyCoord(ra = self.ra * u.degree,
                                  dec = self.dec * u.degree)

class CatalogGeneric(GalaxyCatalog):

    def __init__(self, filename, fieldname=None, verbose=False):

        self.fieldname = fieldname
        self.fieldnameshort = fieldname.replace('-','').lower()
        self.verbose = verbose
        self.catalog_data = Table.read(filename)
        for colname in self.catalog_data.colnames:
            if (colname not in __CANDELS_COLUMN_KEEPLIST__):
                self.catalog_data.remove_column(colname)
                continue
            if (len(self.catalog_data.columns[colname].shape) > 1):
                column_single = self.catalog_data.columns[colname][:,-1]
                self.catalog_data.replace_column(colname, column_single)
                if 'BAND' in colname:
                    modname = colname.replace('BAND','F160W')
                self.catalog_data.columns[colname].name = modname
        fieldcolumn = Column(
            data=[fieldname for i in range(len(self.catalog_data))],
            name='FIELD')
        self.catalog_data.add_column(fieldcolumn)
        self.ra = self.catalog_data['ALPHA_J2000']
        self.dec = self.catalog_data['DELTA_J2000']
        self.locations = SkyCoord(ra = self.ra * u.degree,
                                  dec = self.dec * u.degree)





class CatalogMerged(object):
    """ a catalog that has merged together SNe and galaxies"""
    def __init__(self, filename, **kwargs):
        self.input_filename = filename
        self.catalog = Table.read(filename, **kwargs)
        if 'gal_relation_to_sn' not in self.catalog.colnames:
            galrelationcol = Column(data=['unknown' for
                                          i in range(len(self.catalog))],
                                    name='gal_relation_to_sn')
            self.catalog.add_column(galrelationcol, 0)

    def print_sn_summary(self, nickname=None, id=None):
        if nickname is not None:
            ithissn = np.where(self.catalog['nickname_sn']==nickname)[0]
        elif id is not None:
            ithissn = np.where(self.catalog['id_sn']==id)[0]
        else:
            return
        return(self.catalog[ithissn])


    def check_nearby_galaxies(self, name, revise=True, pixscale=0.03):
        """ make a figure showing a SN position and the relations of nearby
        galaxies.
        :return:
        """
        from matplotlib.patches import Ellipse
        if name in self.catalog['nickname_sn']:
            id_sn = self.catalog['id_sn'][np.where(
                self.catalog['nickname_sn'] == name)[0]][0]
        else:
            id_sn = name

        fig = plt.gcf()
        fig.clf()
        ax = fig.add_subplot(1, 1, 1)
        i_this_sn = np.where(self.catalog['id_sn'] == id_sn)[0]
        rasn = self.catalog['ra_sn'][i_this_sn][0]
        decsn = self.catalog['dec_sn'][i_this_sn][0]
        zsn = self.catalog['z_sn'][i_this_sn][0]
        thisnick = self.catalog['nickname_sn'][i_this_sn][0]
        ax.set_title("%s z=%.2f" % (thisnick, zsn))
        ax.plot(0, 0, marker='x', color='k', ms=10, mew=2)

        galrelationguess = Column(data=self.catalog['gal_relation_to_sn'].data,
                                  name='gal_relation_guess')
        coord_sn = SkyCoord(rasn * u.degree, decsn * u.degree)
        for igal in i_this_sn:
            ragal = self.catalog['ALPHA_J2000'][igal]
            decgal = self.catalog['DELTA_J2000'][igal]
            coord_gal = SkyCoord(ragal * u.degree, decgal * u.degree)
            # sep = coord_sn.separation(coord_gal)
            dra, ddec = coord_sn.spherical_offsets_to(coord_gal)
            # r_arcsec = nearby_galaxies_table['KRON_RADIUS'][igal] * pixscale
            a_arcsec = self.catalog['A_IMAGE'][igal] * pixscale
            b_arcsec = self.catalog['B_IMAGE'][igal] * pixscale
            theta_deg = -self.catalog['THETA_WORLD'][igal]
            zspec_3dhst = self.catalog['z_spec_3dhst'][igal]
            zphot_3dhst = self.catalog['z_phot_3dhst'][igal]
            zgal_3dhst = self.catalog['z_3dhst'][igal]
            Mgal_3dhst = self.catalog['mass_3dhst'][igal]
            mag_gal = self.catalog['MAG_BEST'][igal]
            galrelationold = self.catalog['gal_relation_to_sn'][igal]

            hostsummary = """%i
    zsp,ph=%.2f, %.2f
    log(M)=%.2f""" % (
                igal, zspec_3dhst, zphot_3dhst, np.log10(Mgal_3dhst))
            if not zgal_3dhst > 0:
                hostsummary = "%i" % igal

            if np.abs(zgal_3dhst - zsn) < 0.1:
                ec = 'g'
                galrelationguess[igal] = 'host'
            elif 0 < zgal_3dhst < zsn:
                ec = 'b'
                galrelationguess[igal] = 'fg'
            elif zgal_3dhst > zsn:
                ec = 'r'
                galrelationguess[igal] = 'bg'
            else:
                ec = 'k'
                galrelationguess[igal] = 'unknown'

            if galrelationold == 'unknown':
                fc = 'w'
            elif galrelationold == 'host':
                fc = 'g'
                galrelationguess[igal] = 'host'
            elif galrelationold == 'bg':
                fc = 'r'
                galrelationguess[igal] = 'bg'
            elif galrelationold == 'fg':
                fc = 'b'
                galrelationguess[igal] = 'fg'
            elif galrelationold == 'sn':
                fc = 'c'
                galrelationguess[igal] = 'sn'

            alpha_m = -0.1 * mag_gal + 2.8
            alpha = max(min(alpha_m, 0.8), 0.1)

            xgal = dra.to(u.arcsec).value
            ygal = ddec.to(u.arcsec).value

            ellipse = Ellipse(xy=(xgal, ygal),
                              width=4 * a_arcsec, height=4 * b_arcsec,
                              angle=theta_deg, ec=ec, fc=fc, alpha=alpha,
                              ls='-', lw=2)
            ax.add_artist(ellipse)
            ax.set_xlim(-3., 3.)
            ax.set_ylim(-3., 3.)
            ax.text(xgal, ygal, hostsummary)
        ax.invert_xaxis()
        plt.draw()

        if revise:
            userdone = 'no'
            print(("  Revising gal labels for %s" % thisnick))
            while len(userdone) > 0:
                for igal in i_this_sn:
                    userin = input(
                        "    Galaxy %i is b:bg, f:fg, h:host, s:sn, u:unknown."%igal +
                        " [Return = %s]" % galrelationguess[igal])
                    if userin.lower().startswith('f'):
                        galrelation = 'fg'
                    elif userin.lower().startswith('b'):
                        galrelation = 'bg'
                    elif userin.lower().startswith('h'):
                        galrelation = 'host'
                    elif userin.lower().startswith('s'):
                        galrelation = 'sn'
                    elif userin.lower().startswith('u'):
                        galrelation = 'unknown'
                    else:
                        galrelation = galrelationguess[igal]
                    self.catalog['gal_relation_to_sn'][igal] = galrelation
                self.check_nearby_galaxies(id_sn, revise=False)
                userdone = input(
                    "  Revised %s. Return to continue." % thisnick +
                    "Any other key to repeat. ")
        return


    def checkall_nearby_galaxies(self, istart=0, clobber=False):
        """ review the galaxy relations for all SNe in the catalog
        :return:
        """
        id_sn_list = np.unique(self.catalog['id_sn'])
        nsn = len(id_sn_list)
        isnlist = list(range(nsn))[istart:]
        for isn in isnlist:
            id_sn = id_sn_list[isn]
            print(("SN %i of %i" % (isn, nsn)))
            self.check_nearby_galaxies(id_sn)
            userin = input("Done with SN %i of %i" % (isn, nsn) +
                               " Q to quit, Return to continue.  ")
            if userin.lower().startswith('q'):
                break
        userin = input("Save catalog to %s ? [y]/n  "% self.input_filename)
        if len(userin) == 0 or userin.lower().startswith('y'):
            self.write(clobber=clobber)


    def write(self, filename=None, clobber=False, **kwargs):
        if filename is None:
            filename = self.input_filename

        if os.path.exists(filename) and not clobber:
            userin = input('output table %s exists. ' % filename +
                               'Clobber? [y]/n  ')
            if len(userin)==0 or userin.lower().startswith('y'):
                clobber = True
        if os.path.exists(filename) and clobber:
            os.remove(filename)
        if not os.path.isfile(filename):
            self.catalog.write(filename, **kwargs)
        else:
            print(("file %s exists. Not clobbering." % filename))
        return


    def host_galaxies_table(self):
        """ return a view of the table with only host galaxies"""
        ihost = np.where(self.catalog['gal_relation_to_sn']=='host')[0]
        return self.catalog[ihost]

    def makeall_eazy_seds(self, output_directory='HOSTGAL.SEDS', istart=0,
                          showplot=False):
        id_sn_list = np.unique(self.catalog['id_sn'])
        nsn = len(id_sn_list)
        isnlist = list(range(nsn))[istart:]
        for isn in isnlist:
            id_sn = id_sn_list[isn]

            self.make_eazy_sed(id_sn, output_directory=output_directory,
                               showplot=showplot)
            if showplot:
                plt.draw()
                userin =  input(
                    "Showing host SED for SN %i of %i." % (isn, nsn) +
                    "  Return to continue.")



    def make_eazy_sed(self, name, output_directory='HOSTGAL.SEDS',
                      showplot=False):
        from . import snhostspec
        if name in self.catalog['nickname_sn']:
            ithisobj = np.where(self.catalog['nickname_sn'] == name)[0][0]
        else:
            ithisobj = np.where(self.catalog['id_sn'] == name)[0][0]

        fieldid_3dhst = "%s.%i" % (
            self.catalog['FIELD'][ithisobj],
            self.catalog['id_3dhst'][ithisobj])
        if not os.path.isdir(output_directory):
            os.makedirs(output_directory)
        outfilename = ('%s_host_bestfit_sed.txt' %
                       self.catalog['id_sn'][ithisobj])
        savetofile = os.path.join(output_directory, outfilename)
        output = snhostspec.simulate_eazy_sed(
            fieldid_3dhst, savetofile=savetofile)

        if showplot:
            eazysed = snhostspec.EazySpecSim(savetofile)
            plt.clf()
            eazysed.plot()
            return eazysed
        return output







