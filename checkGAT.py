import os
import sys
import pyvo as vo
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from astroquery.simbad import Simbad
import argparse
from astropy.time import Time
import urllib.request

from ztfquery import lightcurve

from astropy.io import ascii



from gaiaxpy import convert,plot_spectra,calibrate
from astroquery.gaia import Gaia

# To avoid warnings
import warnings
warnings.simplefilter("ignore")

'''

Author: A.Ederoclite

Description:
This script queries VSX (via VizieR), Simbad, Gaia and WISE to find counterparts
of transients, given the coordinates of an object and a search radius.

Dependencies:
- numpy
- matplotlib
- astropy
- astroquery
- pyvo

Todo:
- add Gaia light curve
- add Gaia BP/RP spectra

- add GALEX

- add eROSITA

- add light curve plot: AAVSO, ASAS-SN, ZTF

- add J-surveys

- make more elegant to avoid plotting in case no objects are found


'''

def GATname(mySkycoords):
    myRAHH = str(int(mySkycoords.ra.hms[0])).zfill(2)
    myRAMM = str(int(mySkycoords.ra.hms[1])).zfill(2)
    myRASS = str(int(mySkycoords.ra.hms[2])).zfill(2)
    myDEDD = str(int(abs(mySkycoords.dec.dms[0]))).zfill(2)
    myDEMM = str(int(abs(mySkycoords.dec.dms[1]))).zfill(2)
    myDESS = str(int(abs((mySkycoords.dec.dms[2])))).zfill(2)

    mysign = ''
    if mySkycoords.dec.degree >= 0:
        mysign = '+'
    else:
        mysign = '-'

    myGATname = "GAT" + str(myRAHH) + str(myRAMM)  + str(myRASS) + str(mysign) + str(myDEDD) + str(myDEMM) + str(myDESS)
    return myGATname

def checkSimbad(mySkycoords,myRadius):
    print('== Explore Simbad ==')
    ## https://astroquery.readthedocs.io/en/latest/simbad/simbad.html
    #result_table = Simbad.query_region(mySkycoords, radius=myRadius * u.arcsec)
    ##result_table.to_table()
    #print(result_table.info())
    simbad_service = vo.dal.TAPService("http://simbad.cds.unistra.fr/simbad/sim-tap")
    ra = mySkycoords.ra.degree  # Right Ascension in degrees
    dec = mySkycoords.dec.degree   # Declination in degrees
    radius = myRadius/3600.  # Search radius in degrees
    query = (
    f"SELECT * FROM basic "
    f"WHERE 1=CONTAINS(POINT('ICRS',basic.ra,basic.dec), "
    f"CIRCLE('ICRS',{ra},{dec},{radius}))"
    )
    result_table = simbad_service.search(query)
    result_table = result_table.to_table()
    #print(result_table.info())    
    return result_table

def checkVSX(mySkycoords,myRadius):
    print('== Explore VSX (via VizieR) ==')
    ## https://astroquery.readthedocs.io/en/latest/simbad/simbad.html
    #result_table = Simbad.query_region(mySkycoords, radius=myRadius * u.arcsec)
    ##result_table.to_table()
    #print(result_table.info())
    vsx_service = vo.dal.TAPService("http://TAPVizieR.u-strasbg.fr/TAPVizieR/tap/")
    ra = mySkycoords.ra.degree  # Right Ascension in degrees
    dec = mySkycoords.dec.degree   # Declination in degrees
    radius = myRadius/3600.  # Search radius in degrees
    query = (
    f"SELECT * FROM \"B/vsx/vsx\" "
    f"WHERE 1=CONTAINS(POINT('ICRS',\"B/vsx/vsx\".RAJ2000,\"B/vsx/vsx\".DEJ2000), "
    f"CIRCLE('ICRS',{ra},{dec},{radius}))"
    )
    result_table = vsx_service.search(query)
    result_table = result_table.to_table()
    #print(result_table.info())    
    return result_table

def checkGaia(myRA,myDec,myRadius):
    print('== Explore Gaia ==')

    gaia_service = vo.dal.TAPService("https://gea.esac.esa.int/tap-server/tap")
    # Define the position and search radius
    ra = myRA  # Right Ascension in degrees
    dec = myDec   # Declination in degrees
    radius = myRadius/3600.  # Search radius in degrees

    query = (
    f"SELECT * FROM gaiadr3.gaia_source "
    f"WHERE 1=CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec), "
    f"CIRCLE('ICRS',{ra},{dec},{radius}))"
    )
    print(query)
    result = gaia_service.search(query)
    #print(result)
    print(len(result))
    result = result.to_table()

    bck_result = ''
    if len(result) > 0:
        print('== Explore Gaia Background ==')
        bck_result = ''
        while len(bck_result) < 1000 :
            # Formulate the ADQL query
            query = (
            f"SELECT * FROM gaiadr3.gaia_source "
            f"WHERE 1=CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec), "
            f"CIRCLE('ICRS',{ra},{dec},{radius}))"
            )
            bck_result = gaia_service.search(query)
            #  print(result)
            print(len(bck_result))
            radius = radius * 2.
        bck_result = bck_result.to_table()
        #print(bck_result.info())
    return result, bck_result

def plotGaia(myname,myRA,myDec,myGaiaObjects,myGaia_bck_Objects,mySimbadResults,myVSXResults):
    print('== Plot Gaia ==')

    if len(myGaiaObjects) == 0:
        print('No Gaia objects')
    else:
        print('Gaia objects')

    # Create the figure
    fig = plt.figure(figsize=(8.27, 11.69))

    # Create a GridSpec with 3 rows and 3 columns
    gs = gridspec.GridSpec(3, 3)

    # First row: 3 plots
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    # Second row: 2 plots (spanning over 3 columns)
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[1, 2])

    # Third row: 2 plots (first two spanning 2 columns)
    ax7 = fig.add_subplot(gs[2, :2])
    ax8 = fig.add_subplot(gs[2, 2])

    # Add some data or customize the axes here
    #ax1.plot([1, 2, 3], [1, 4, 9])
    ax1.text(0.2, 0.8, 'Name:' + myname, horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)
    ax1.text(0.2, 0.6, 'RA:' + str(myRA), horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)
    ax1.text(0.2, 0.4, 'Dec:' + str(myDec), horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)
    ax1.set_axis_off()
    #ax2.plot([1, 2, 3], [1, 4, 9])
    if len(mySimbadResults) > 0:
        ax2.text(0.2, 0.8, 'Simbad Name: ' + mySimbadResults['main_id'][0], horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes)
        ax2.text(0.2, 0.6, 'Simbad Type: ' + mySimbadResults['otype_txt'][0], horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes)
        ax2.text(0.2, 0.4, 'Simbad Sp.Type: ' + mySimbadResults['sp_type'][0], horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes)
    else:
        ax2.text(0.2, 0.6, 'No Simbad detection', horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes)
    ax2.set_axis_off()

    if len(myVSXResults) > 0:
        ax3.text(0.2, 0.8, 'VSX Name: ' + myVSXResults['Name'][0], horizontalalignment='left', verticalalignment='center', transform=ax3.transAxes)
        ax3.text(0.2, 0.6, 'VSX Type: ' + myVSXResults['Type'][0], horizontalalignment='left', verticalalignment='center', transform=ax3.transAxes)
        ax3.text(0.2, 0.4, 'VSX Period [d]: ' + str(myVSXResults['Period'][0]), horizontalalignment='left', verticalalignment='center', transform=ax3.transAxes)
    else:
        ax3.text(0.2, 0.6, 'No VSX detection', horizontalalignment='left', verticalalignment='center', transform=ax3.transAxes)
    ax3.set_axis_off()

    #ax3.plot([1, 2, 3], [1, 4, 9])
    #tmp = myGaiaObjects[ myGaiaObjects['ruwe'] < 1.4 ]
    #tmp = tmp[ tmp['parallax_over_error'] > 5. ]
    ax4.scatter(myGaia_bck_Objects['bp_rp'],myGaia_bck_Objects['phot_g_mean_mag'])
    ax4.scatter(myGaiaObjects['bp_rp'],myGaiaObjects['phot_g_mean_mag'])
    ax4.invert_yaxis()
    ax4.set_xlabel('BP - RP')
    ax4.set_ylabel('Gmag')

    tmp = myGaiaObjects[ myGaiaObjects['ruwe'] < 1.4 ]
    tmp = tmp[ tmp['parallax_over_error'] > 5. ]

    tmp_bck = myGaia_bck_Objects[ myGaia_bck_Objects['ruwe'] < 1.4 ]
    tmp_bck = tmp_bck[ tmp_bck['parallax_over_error'] > 5. ]

    ax5.scatter(tmp_bck['bp_rp'],tmp_bck['phot_g_mean_mag'] + 5 - 5 * np.log10(1000./tmp_bck['parallax']))
    ax5.scatter(tmp['bp_rp'],tmp['phot_g_mean_mag'] + 5 - 5 * np.log10(1000./tmp['parallax']))
    ax5.invert_yaxis()
    ax5.set_xlabel('BP - RP')
    ax5.set_ylabel('Gmag [abs]')

    ax6.scatter(myGaia_bck_Objects['pmra'],myGaia_bck_Objects['pmdec'])
    ax6.scatter(myGaiaObjects['pmra'],myGaiaObjects['pmdec'])
    ax6.set_xlabel('proper motion RA')
    ax6.set_ylabel('proper motion DEC')


    if len(myGaiaObjects) == 1:
        sources = ['679528804789642240']
        print(sources)
        sources = [myGaiaObjects['source_id'][0]]
        print(sources)
        calibrated_Gaia_spectrum, sampling = calibrate(sources,save_file=False)
        #plot_spectra(converted_data, sampling=sampling, multi=False, show_plot=True)
        ax7.plot(sampling,calibrated_Gaia_spectrum['flux'][0])


        # https://www.cosmos.esa.int/web/gaia-users/archive/datalink-products#datalink_jntb_get_all

        # https://gea.esac.esa.int/data-server/datalink/links?ID=Gaia+DR3+30343944744320
        # https://gea.esac.esa.int/data-server/datalink/links?ID=Gaia+DR3+679528804789642240

        retrieval_type = 'EPOCH_PHOTOMETRY'  # Options are: 'EPOCH_PHOTOMETRY', 'MCMC_GSPPHOT', 'MCMC_MSC', 'XP_SAMPLED', 'XP_CONTINUOUS', 'RVS', 'ALL'
        data_structure = 'INDIVIDUAL'   # Options are: 'INDIVIDUAL', 'COMBINED', 'RAW'
        data_release   = 'Gaia DR3'     # Options are: 'Gaia DR3' (default), 'Gaia DR2'
        datalink = Gaia.load_data(ids=sources, data_release = data_release, retrieval_type=retrieval_type, data_structure = data_structure, verbose = False, output_file = None)
        dl_keys  = [inp for inp in datalink.keys()]
        dl_keys.sort()
        print()
        print(f'The following Datalink products have been downloaded:')
        for dl_key in dl_keys:
            print(f' * {dl_key}')
            inp_table  = datalink[dl_key][0].to_table()
            ax8.set_xlabel   = f'JD date [{inp_table["time"].unit}]'
            ax8.set_ylabel   = f'magnitude [{inp_table["mag"].unit}]'
            gbands   = ['G', 'RP', 'BP']
            colours  = iter(['green', 'red', 'blue'])
            for band in gbands:
                phot_set = inp_table[inp_table['band'] == band]
                ax8.plot(phot_set['time'], phot_set['mag'], 'o', label = band, color = next(colours))
            ax8.legend()

    #ax6.scatter(myGaiaObjects['bp_rp'],myGaiaObjects['phot_g_mean_mag'])


    #ax4.plot([1, 2, 3], [1, 4, 9])
    #ax5.plot([1, 2, 3], [1, 4, 9])
    #ax6.plot([1, 2, 3], [1, 4, 9])
    #ax7.plot([1, 2, 3], [1, 4, 9])

    ## Adjust layout to prevent overlap
    #plt.tight_layout()

    # Save the figure
    if not os.path.exists('./'+myname):
        os.mkdir('./'+myname)
    plt.savefig('./'+ myname + '/' + myname + '_Gaia.png')


    # Show the figure
    #plt.show()

def checkWISE(myRA,myDec,myRadius):#, max_rows=10):
    print('== Explore WISE ==')
    """
    Query the AllWISE catalog using a cone search.

    Parameters:
    - myRA: Right Ascension in degrees
    - myDec: Declination in degrees
    - myRadius: Search radius in degrees (default is 1.0 degree)

    Returns:
    - result: astropy table result containing the query results
    """

    # Define the position and search radius
    ra = myRA  # Right Ascension in degrees
    dec = myDec   # Declination in degrees
    radius = myRadius/3600.  # Search radius in degrees

    ## AllWISE
    print('== Explore AllWISE ==')

    # Formulate the ADQL query
    query = (
        f"SELECT * FROM allwise_p3as_psd WHERE "
        f"CONTAINS(POINT('ICRS', allwise_p3as_psd.ra, allwise_p3as_psd.dec), "
        f"CIRCLE('ICRS', {ra}, {dec}, {radius}))"
    )

    # Set up the AllWISE TAP service
    allwise_service = vo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")

    # Execute the query
    allwise_result = allwise_service.search(query)#, maxrec=max_rows)

    print(len(allwise_result))
    allwise_bck_result = allwise_result

    if len(allwise_result) >  0:

        current_radius = radius
        while len(allwise_bck_result) < 100:
                query = (
                f"SELECT * FROM allwise_p3as_psd WHERE "
                f"CONTAINS(POINT('ICRS', allwise_p3as_psd.ra, allwise_p3as_psd.dec), "
                f"CIRCLE('ICRS', {ra}, {dec}, {current_radius}))"
                )
                allwise_bck_result = allwise_service.search(query)#, maxrec=max_rows)
                current_radius = current_radius * 2.
                print(round(current_radius*3600,2),len(allwise_bck_result))
        allwise_bck_result = allwise_bck_result.to_table()
    #print(bck_result.info())

    ## CatWISE
    print('== Explore CatWISE ==')

    # Formulate the ADQL query
    query = (
        f"SELECT * FROM catwise_2020 WHERE "
        f"CONTAINS(POINT('ICRS', catwise_2020.ra, catwise_2020.dec), "
        f"CIRCLE('ICRS', {ra}, {dec}, {radius}))"
    )

    # Set up the AllWISE TAP service
    catwise_service = vo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")

    # Execute the query
    catwise_result = catwise_service.search(query)#, maxrec=max_rows)

    print(len(catwise_result))
    catwise_bck_result = catwise_result

    if len(catwise_result) >  0:

        current_radius = radius
        while len(catwise_bck_result) < 100:
            query = (
            f"SELECT * FROM catwise_2020 WHERE "
            f"CONTAINS(POINT('ICRS', catwise_2020.ra, catwise_2020.dec), "
            f"CIRCLE('ICRS', {ra}, {dec}, {current_radius}))"
            )
            catwise_bck_result = catwise_service.search(query)#, maxrec=max_rows)
            current_radius = current_radius * 2.
            print(round(current_radius*3600,2),len(catwise_bck_result))
        catwise_bck_result = catwise_bck_result.to_table()

    return allwise_result, allwise_bck_result, catwise_result, catwise_bck_result

def plotWISE(myname,myRA,myDec,mySimbad_result,myVSX_result,myWISE_result, myWISE_bck_result,myCatWISE_result, myCatWISE_bck_result):
    print('== Plot WISE ==')

    # Create the figure
    fig = plt.figure(figsize=(8.27, 11.69))

    fig.suptitle(myname + ' WISE')

    # Create a GridSpec with 3 rows and 3 columns
    gs = gridspec.GridSpec(3, 3)

    # First row: 3 plots
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    # Second row: 3 plots
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[1, 2])

    # Third row: 2 plots (first two spanning 2 columns)
    ax6 = fig.add_subplot(gs[2, 0])
    #ax7 = fig.add_subplot(gs[2, 2])

    ax1.text(0.2, 0.8, 'Name:' + myname, horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)
    ax1.text(0.2, 0.6, 'RA:' + str(myRA), horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)
    ax1.text(0.2, 0.4, 'Dec:' + str(myDec), horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)
    ax1.set_axis_off()

    if len(mySimbad_result) > 0:
        ax2.text(0.2, 0.8, 'Simbad Name: ' + mySimbad_result['main_id'][0], horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes)
        ax2.text(0.2, 0.6, 'Simbad Type: ' + mySimbad_result['otype_txt'][0], horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes)
        ax2.text(0.2, 0.4, 'Simbad Sp.Type: ' + mySimbad_result['sp_type'][0], horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes)
    else:
        ax2.text(0.2, 0.6, 'No Simbad detection', horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes)
    ax2.set_axis_off()

    if len(myVSX_result) > 0:
        ax3.text(0.2, 0.8, 'VSX Name: ' + myVSX_result['Name'][0], horizontalalignment='left', verticalalignment='center', transform=ax3.transAxes)
        ax3.text(0.2, 0.6, 'VSX Type: ' + myVSX_result['Type'][0], horizontalalignment='left', verticalalignment='center', transform=ax3.transAxes)
        ax3.text(0.2, 0.4, 'VSX Period [d]: ' + str(myVSX_result['Period'][0]), horizontalalignment='left', verticalalignment='center', transform=ax3.transAxes)
    else:
        ax3.text(0.2, 0.6, 'No VSX detection', horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes)
    ax3.set_axis_off()

    if len(myWISE_result) > 0:
        ax4.scatter(myWISE_result['w1mag']-myWISE_result['w2mag'],myWISE_result['w2mag']-myWISE_result['w3mag'])
    ax4.scatter(myWISE_bck_result['w1mag']-myWISE_bck_result['w2mag'],myWISE_bck_result['w2mag']-myWISE_bck_result['w3mag'])
    ax4.set_xlabel('W1 - W2 (AllWISE)')
    ax4.set_ylabel('W2 - W3 (AllWISE)')

    ax6.hist(myCatWISE_bck_result['w1mag']-myCatWISE_bck_result['w2mag'])
    if len(myCatWISE_result) > 0:
        tmp = myCatWISE_result['w1mag']-myCatWISE_result['w2mag']
        for eachcolour in tmp:
            ax6.axvline(eachcolour,color='r')
    ax6.set_xlabel('W1 - W2 (CatWISE)')

    ## Adjust layout to prevent overlap
    #plt.tight_layout()

    # Save the figure
    if not os.path.exists('./'+myname):
        os.mkdir('./'+myname)
    plt.savefig('./'+ myname + '/' + myname + '_WISE.png')

    # Show the figure
    #plt.show()

def checkAAVSOLC(myVSXResults):
    
    # AAVSO

    # https://www.aavso.org/vsx/index.php?view=api.object&ident=000-BCQ-471&data=50000&fromjd=2459752.5&tojd=2460483.5&csv&band=Vis.,V&mtype=std

    # https://www.aavso.org/vsx/index.php?view=api.object&ident=T-Pyx&data=50000&fromjd=2459752.5&tojd=2460483.5&csv&band=Vis.,V&mtype=std

    # https://www.aavso.org/programmatic-access-aavso-light-curves

    jds = []
    mags = []
    if len(myVSXResults) > 0:
        t = Time.now()
        t_jd = round(t.jd,2)
        print(t_jd)
        t0 = t_jd - 1000.
        myurl = "https://www.aavso.org/vsx/index.php?view=api.object&ident=" + myVSXResults['Name'][0] + "&data=50000&fromjd="
        myurl = myurl + str(t0) + "&tojd=" + str(t_jd) + "&csv&band=Vis.,V&mtype=std"
        #print(myurl)
        myurl= myurl.replace(' ','-')
        # https://www.aavso.org/vsx/index.php?view=api.object&ident=AT Cnc&data=50000&fromjd=2459558.07&tojd=2460558.07&csv&band=Vis.,V&mtype=std

        try: 
            urllib.request.urlopen(myurl)
            urllib.request.urlretrieve(myurl, 'tmp.votable')
        except urllib.error.URLError as e:
            print(e.reason)

        myfile = open('tmp.votable','r')
        for index,eachline in enumerate(myfile):
            if index > 1 :
                print(index, eachline)
                tmp = eachline.strip().split(',')
                try :
                    jds.append(float(tmp[0]))
                    mags.append(float(tmp[1]))
                except:
                    pass
        myfile.close()
        os.system('rm tmp.votable')
    return [jds,mags]

# ASAS-SN

# ZTF

def checkZTF(mySkycoords,myRadius):
    print('== Explore ZTF (via ZTFQuery) ==')
    print('you may need to be patiend with this')

    ra = mySkycoords.ra.degree  # Right Ascension in degrees
    dec = mySkycoords.dec.degree   # Declination in degrees

    # https://github.com/MickaelRigault/ztfquery/blob/master/doc/lightcurve.md
    lcq = lightcurve.LCQuery.from_position(ra, dec, myRadius)
    #print(lcq.data)

    return lcq.data

def plotLCs(myname,myRA,myDec,mySimbad_result,myVSX_result,myAAVSOLC,myZTFLC):

    print('== Plot LCs ==')

    # Create the figure
    fig = plt.figure(figsize=(8.27, 11.69))

    fig.suptitle(myname + ' Light Curves')

    # Create a GridSpec with 3 rows and 3 columns
    gs = gridspec.GridSpec(4, 3)

    # First row: 3 plots
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    # Second row: 3 plots
    ax4 = fig.add_subplot(gs[1, :2])
    #ax5 = fig.add_subplot(gs[1, 1])
    #ax6 = fig.add_subplot(gs[1, 2])

    # Third row: 2 plots (first two spanning 2 columns)
    ax6 = fig.add_subplot(gs[2, 0])
    #ax7 = fig.add_subplot(gs[2, 2])

    ax1.text(0.2, 0.8, 'Name:' + myname, horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)
    ax1.text(0.2, 0.6, 'RA:' + str(myRA), horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)
    ax1.text(0.2, 0.4, 'Dec:' + str(myDec), horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)
    ax1.set_axis_off()

    if len(mySimbad_result) > 0:
        ax2.text(0.2, 0.8, 'Simbad Name: ' + mySimbad_result['main_id'][0], horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes)
        ax2.text(0.2, 0.6, 'Simbad Type: ' + mySimbad_result['otype_txt'][0], horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes)
        ax2.text(0.2, 0.4, 'Simbad Sp.Type: ' + mySimbad_result['sp_type'][0], horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes)
    else:
        ax2.text(0.2, 0.6, 'No Simbad detection', horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes)
    ax2.set_axis_off()

    if len(myVSX_result) > 0:
        ax3.text(0.2, 0.8, 'VSX Name: ' + myVSX_result['Name'][0], horizontalalignment='left', verticalalignment='center', transform=ax3.transAxes)
        ax3.text(0.2, 0.6, 'VSX Type: ' + myVSX_result['Type'][0], horizontalalignment='left', verticalalignment='center', transform=ax3.transAxes)
        ax3.text(0.2, 0.4, 'VSX Period [d]: ' + str(myVSX_result['Period'][0]), horizontalalignment='left', verticalalignment='center', transform=ax3.transAxes)
    else:
        ax3.text(0.2, 0.6, 'No VSX detection', horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes)
    ax3.set_axis_off()

    if len(myAAVSOLC[0]) > 0:
        #print(jds,mags)
        #plt.scatter(jds,mags)
        #plt.show()
        ax4.scatter(myAAVSOLC[0],myAAVSOLC[1])
        ax4.invert_yaxis()

    if len(myZTFLC) > 0:
        #print(jds,mags)
        #plt.scatter(jds,mags)
        #plt.show()
        ax6.scatter(myZTFLC['mjd'],myZTFLC['mag'])
        ax6.invert_yaxis()

    # Save the figure
    if not os.path.exists('./'+myname):
        os.mkdir('./'+myname)
    plt.savefig('./'+ myname + '/' + myname + '_LCs.png')


# eROSITA
# https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/A+A/682/A34/erass1-m    

# J-PLUS

# Construct the argument parser
ap = argparse.ArgumentParser()
ap.add_argument("ra", help="right ascension")
ap.add_argument("dec", help="declination")
ap.add_argument("-f", required=True, default='x',help="format of coordinates x for hexadecimal, d for decimal degrees")
ap.add_argument("-r", required=True, help="the search radius (in arcseconds)")
args = vars(ap.parse_args())

# ============================================================================= #

inputra = args['ra']
inputdec = args['dec']
inputRadius = float(args['r'])

inputUnitType = args['f']
if inputUnitType == 'x':
    myUnitType = (u.hourangle, u.deg)
elif inputUnitType == 'd':
    myUnitType = (u.deg, u.deg)
else:
    sys.exit('invalid format of coordinates, it must be either x (hexadecimal) or d (decimal degrees)')

# ============================================================================= #

'''
python checkGAT.py -f d -r 5 303.8972306 -14.2790126

python checkGAT.py -f x -r 5 07:25:53.03 +29:04:10.2
'''

# AT Cnc
# python checkGAT.py -f x -r 5 08:28:36.92 +25:20:03.04

# JVAR binary
# python checkGAT.py -f x -r 5 14:05:03.34 +33:16:13.68

c = SkyCoord(inputra, inputdec, frame='icrs',unit=myUnitType)
print(c.ra.degree,c.dec.degree)
name = GATname(c)
Simbad_result = checkSimbad(c,5)
VSX_result = checkVSX(c,5)

GaiaResults,Gaia_bck_results = checkGaia(c.ra.degree,c.dec.degree,inputRadius)
if len(GaiaResults) > 0:
    plotGaia(name,c.ra,c.dec,GaiaResults,Gaia_bck_results,Simbad_result,VSX_result)
else:
    print('no Gaia results')
allwise_result, allwise_bck_result, catwise_result, catwise_bck_result = checkWISE(c.ra.degree,c.dec.degree,5.0)
if len(allwise_result) == 0 and len(catwise_result) == 0:
    print('no WISE results')
else:
    plotWISE(name,c.ra,c.dec,Simbad_result,VSX_result,allwise_result, allwise_bck_result, catwise_result, catwise_bck_result)

AAVSOLC = checkAAVSOLC(VSX_result)
ZTFLC = checkZTF(c,inputRadius)
plotLCs(name,c.ra,c.dec,Simbad_result,VSX_result,AAVSOLC,ZTFLC)
