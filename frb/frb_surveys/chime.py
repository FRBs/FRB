""" CHIME/FRB calculations """
import os
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

from importlib import resources
#%matplotlib inline
import json
import pandas

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18 as cosmo

from frb.galaxies import hosts

try:
    from zdm.chime import grids
except ImportError:
    print('WARNING:  zdm.chime.grids not available')
    grids = None

#import dustmaps.sfd
#dustmaps.sfd.fetch()

from IPython import embed

def check_frb_mr(outfile:str='mr_pdf.png',
                 pdf_file:str=None):
    """ Generate a plot of the host galaxy M_r PDF

    Args:
        outfile (str, optional): _description_. Defaults to 'mr_pdf.png'.
        pdf_file (str, optional): _description_. Defaults to 'pdf_Mr.npy'.
    """

    #host galaxy M_r
    xvals, prob1 = hosts.load_Mr_pdf(pdf_file)

    plt.clf()
    ax = plt.gca()

    ax.plot(xvals, prob1, color='k', lw=2, 
            label='Host Galaxy M_r PDF')

    ax.legend()

    ax.set_xlabel(r'$M_r$')
    ax.set_ylabel(r'PDF')

    #plt.tight_layout(pad=0.5, h_pad=0.5, w_pad=0.5)
    plt.show()
    plt.savefig(outfile, dpi=300)
    plt.close()
    print('Wrote {:s}'.format(outfile))


    

def calc_mr_dist(catalog_file:str=None, 
                 pdf_file:str=None,
                 figfile:str=None,
                 tblfile:str='CHIME_mr_5Jyms.csv',
                 dm_mw_host:float=200.):
    """ Generate a distribution of m_r values for bright CHIME/FRBs

    Args:
        catalog_file (str, optional): 
            Filename for the CHIME/FRB catalog. Defaults to None.
        pdf_file (str, optional): 
            Filename for the host galaxy M_r PDF. Defaults to None.
        figfile (str, optional): 
            Filename for the output figure. Defaults to None.
        tblfile (str, optional): 
            Filename for the output table. Defaults to 'CHIME_mr_5Jyms.csv'
        dm_mw_host (float, optional): 
            Aveage DM contribution from the Milky Way and Host. Defaults to 200.
    """
    # Hiding this here to avoid a dependency
    from dustmaps.sfd import SFDQuery

    # Generate p(z|DM)
    #  Requires zdm
    dmvals, zvals, all_rates, all_singles, all_reps = grids.load()

    # Load up
    if catalog_file is None:
        catalog_file = os.path.join(resources.files('frb'), 'data', 
                                    'FRBs', 'CHIME_catalog-2021-1-27.json')

    #host galaxy M_r
    xvals, prob1 = hosts.load_Mr_pdf(pdf_file)

    dm_excess = []
    fluence = []
    ra = []
    dec = []
    with open(catalog_file) as json_file:
        data = json.load(json_file)
        for i in data:
            rep = i['repeater_of']   
            if len(rep) == 0:
                fluence.append(i['fluence'])
                ra.append(i['ra'])
                dec.append(i['dec'])
                dm_excess.append((i['dm_excess_ne2001']+i['dm_excess_ymw16'])/2)
    dm_excess = np.array(dm_excess)
    fluence = np.array(fluence)
    ra = np.array(ra)
    dec = np.array(dec)
    index_5_Jyms = np.where(fluence >=5)[0]

    n_samples = len(index_5_Jyms) #total number of bright FRBs in CHIME/FRB catalog 1

    # Deal with extinction
    coords = SkyCoord(ra*u.degree, dec*u.degree, frame='icrs')[index_5_Jyms]
    sfd = SFDQuery()
    Ar = sfd(coords)*2.285 # SDSS r-band
    index_Ar = np.arange(len(Ar))
    Index_Ar = np.random.shuffle(index_Ar)
    mw_extinction = np.squeeze(Ar[Index_Ar]) #don't know why but np.random.shuffle(mw_extinction) is not working

    # This needs to be replaced with a more sophisticated calculation from zdm
    z = (dm_excess[index_5_Jyms] - dm_mw_host) / 935
    z[z < 0.01] = (dm_excess[index_5_Jyms][z < 0.01] - 30) / 935
    redshift = np.array(z)

    dist_mod = cosmo.distmod(redshift).value

    r_mag_1 = 20
    r_mag_2 = 22
    r_mag_3 = 24
    r_mag_4 = 26

    samples = 100000
    frac = []
    m_rs= []
    Dist_s= []
    for i in range(samples):
        M_r = np.random.choice(xvals, n_samples, p=prob1)
        Index_ = np.random.shuffle(np.arange(len(Ar)))
        mw_extinction = np.squeeze(Ar[Index_])
        Dist_mod = dist_mod[Index_]
        host_m_r = Dist_mod + M_r + mw_extinction
        frac1 = len(np.where(host_m_r <= r_mag_1)[0])/n_samples
        frac2 = len(np.where((host_m_r > r_mag_1) & (host_m_r <= r_mag_2))[0])/n_samples
        frac3 = len(np.where((host_m_r > r_mag_2) & (host_m_r <= r_mag_3))[0])/n_samples
        frac4 = len(np.where((host_m_r > r_mag_3) & (host_m_r <= r_mag_4))[0])/n_samples
        frac5 = len(np.where(host_m_r > r_mag_4)[0])/n_samples
        val = np.squeeze(np.array([frac1,frac2,frac3,frac4,frac5]))
        # Save
        frac.append(val)
        m_rs.append(host_m_r.flatten())
        Dist_s.append(Dist_mod.flatten())

    frac_ = np.round(np.mean(frac,axis=0),2)

    if figfile:
        plt.style.use('seaborn-poster')
        fig, ax = plt.subplots(figsize=(8,6))

        # Save the chart so we can loop through the bars below.
        bars = ax.bar(
            x= [0,1,2,3,4],
            height= frac_ ,
            tick_label=['$<$ 20','20 $-$ 22','22 $-$ 24','24 $-$ 26','$>$ 26']
        )

        # Axis formatting.
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_color('#DDDDDD')
        ax.tick_params(bottom=False, left=False)
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color='#EEEEEE')
        ax.xaxis.grid(False)

        # Add text annotations to the top of the bars.
        bar_color = bars[0].get_facecolor()
        for bar in bars:
            ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.01,
            round(bar.get_height(), 2),
            horizontalalignment='center',
            color=bar_color,
            weight='bold'
        )

        ax.set_xlabel('R-band Apparent magnitude [AB] ', labelpad=15, color='#333333')
        ax.set_ylabel('Fraction of CHIME FRBs (1 YR)', labelpad=15, color='#333333')
        ax.set_title('R-band Magnitude of Bright CHIME FRBs (Fluence > 5 Jy ms)', pad=15, color='#333333',
                    weight='bold')

        fig.tight_layout()

        #plt.tight_layout(pad=0.5, h_pad=0.5, w_pad=0.5)
        plt.savefig(figfile, dpi=300)
        plt.close()
        print('Wrote {:s}'.format(figfile))

    # Save the m_r distribution
    df = pandas.DataFrame(
        dict(m_r=np.concatenate(m_rs),
             Dist=np.concatenate(Dist_s)))

    df.to_parquet(tblfile)

        
# Command line execution
#if __name__ == '__main__':
#    check_frb_mr()
#
#    # Runs
#    #calc_mr_dist()
#    #calc_mr_dist(figfile='mr_dist_150.png',
#    #             tblfile='CHIME_mr_5Jyms_150.parquet',
#    #             dm_mw_host=150.)
