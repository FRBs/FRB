""" Module for Tables for the Alopeke paper """
# Imports
from signal import raise_signal
import numpy as np
import os, sys
import copy

import pandas

from astropy import units
from astropy.table import Table 


# Local
sys.path.append(os.path.abspath("../../Analysis/py"))
import alopeke_defs
import alopeke_utils2
import alopeke_analy

from IPython import embed

# Summary table of results
def mktab_photom(camera, outroot='tab_photom_', sub=False):

    # Open
    if sub:
        outfile = outroot+camera+'.tex'
    else:
        outfile = outroot+camera+'_sub.tex'
    tbfil = open(outfile, 'w')
    lbl = 'Red' if camera == 'red' else 'Blue'
    filt = '$i$' if camera == 'red' else '$r$'
    gain = alopeke_defs.gain_red if camera == 'red' else alopeke_defs.gain_blue

    # Load
    #data = Table.read('../Data/master_table_redcam.fits')
    data_dict = alopeke_utils2.load_camera(camera)

    # Header
    #tbfil.write('\\clearpage\n')
    tbfil.write('\\begin{deluxetable}{cccccccccccccccc}\n')
    #tbfil.write('\\rotate\n')
    tbfil.write('\\tablewidth{0pc}\n')
    tbfil.write('\\tablecaption{Photometry for '+lbl+' Camera \\label{tab:photom_'+camera+'}}\n')
    tbfil.write('\\tabletypesize{\\footnotesize}\n')
    tbfil.write('\\tablehead{\\colhead{MJD} & \\colhead{Filter} \n')
    #tbfil.write('& \\colhead{$x_1$} & \\colhead{$y_1$} & \\colhead{FWHM$_1$}\n')
    tbfil.write('& \\colhead{\cbkg} & \\colhead{\cstar} & \\colhead{\ctfrb} \n')
    #tbfil.write("& (deg) & (deg) & ($''$) & & (deg) & (deg) & ($''$) & ($''$) & ($''$) & ($''$) & (mag)\n")
    #tbfil.write("\\\\ (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) & (10) & (11) & (12) & (13) & (14) & (15)")
    tbfil.write('} \n')

    tbfil.write('\\startdata \n')

    # Loop 
    for kk in range(len(data_dict['C_star1_full'])):
        if sub and kk > 10:
            continue
        #if data_dict['frame'][kk] == 1:
        #    raise ValueError("Bad frame value!!")

        sline = ''

        # Name and filter
        sline += '{}'.format(data_dict['MJD_star1'][kk]) + '& '+filt

        # x1, y1
        #sline += '&' + '{:0.2f}'.format(row['x_star1'])
        #sline += '&' + '{:0.2f}'.format(row['y_star1'])

        # FWHM
        #sline += '&' + 'TBD'

        # C_bkg
        sline += '&' + '{:0.2f}'.format(data_dict['C_bkg_full'][kk])
        # C_star
        sline += '&' + '{:0.2f}'.format(data_dict['C_star1_full'][kk])
        # C FRB
        sline += '&' + '{:0.2f}'.format(data_dict['C_gal_full'][kk])

        tbfil.write(sline + '\\\\ \n')

    # End end
    tbfil.write('\\hline \n')


    tbfil.write('\\enddata \n')

    #tbfil.write('\\tablenotetext{a}{Spectroscopic redshifts are reported to 4 significant digits.  Photometric to 2.} \n')
    #tbfil.write('the gas below the line $\\rm DEC_{\\rm off} = \\aslope RA_{\\rm off} \\ayoff$}\n')
    #tbfil.write('\\tablecomments{Column 1: FRB source. Columns 2 and 3: R.A. and Decl. of the FRB (J2000). Column 4: FRB error ellipse. Column 5: FRB classication. Repeating = yes(y)/no(n). Column 6 and 7: R.A. and Dec. of the associated host galaxy (J2000). Column 8: projected angular offset of the FRB to the host galaxy center. Column 9: association radius $\delta x$ \citep{Tunnicliffe14}. Column 10: angular effective radius of the host measured from a sersic model using GALFIT \citep{galfit} on the $i$-band images (or equivalent). Column 11: effective search radius \citep{Bloom02}. Column 12: measured apparent magnitude of the host. Column 13: filter used for the magnitude measurement. Column 14: probability of chance coincidence using the \citet{Bloom02} formalism. Column 15: sample designations following the criteria outlined in $\S$~\\ref{ssec:associate}.}\n')
    # End
    tbfil.write('\\end{deluxetable} \n')

    tbfil.close()
    print('Wrote {:s}'.format(outfile))



# Summary table of results
def mktab_summary(outfile='tab_summary.tex'):

    # Load up
    summary = alopeke_analy.summary_dict(calc_upper=True)

    # Items to show
    items = ['A', 'C_1', 'C_gal', 'sC_gal', #'dCdt',
             'mxFRB', 'uppFRB', 'uppFlu']

    # Header
    #tbfil.write('\\clearpage\n')
    tbfil = open(outfile, 'w')
    tbfil.write('\\begin{deluxetable*}{cccccccccccccccc}\n')
    #tbfil.write('\\rotate\n')
    tbfil.write('\\tablewidth{0pc}\n')
    tbfil.write('\\tablecaption{Summary of Results\\label{tab:summary}}\n')
    tbfil.write('\\tabletypesize{\\footnotesize}\n')
    tbfil.write('\\tablehead{\\colhead{Item} & \\colhead{Camera} \n')
    tbfil.write('& \\colhead{Value} & \\colhead{Error} & \\colhead{Unit}\n')
    tbfil.write('& \\colhead{Desc.} \n')
    #tbfil.write("& (deg) & (deg) & ($''$) & & (deg) & (deg) & ($''$) & ($''$) & ($''$) & ($''$) & (mag)\n")
    #tbfil.write("\\\\ (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) & (10) & (11) & (12) & (13) & (14) & (15)")
    tbfil.write('} \n')

    tbfil.write('\\startdata \n')

    # Loop 
    for item in items:
        for camera in ['blue', 'red']:
            sline = ''

            # Name and filter
            sline += summary[camera][item]['variable'] + '&'+camera
 
            # Value 
            sline += '&' + summary[camera][item]['vformat'].format(
                summary[camera][item]['value'])
            # Error
            sline += '&' + summary[camera][item]['eformat'].format(
                summary[camera][item]['error'])
            # Unit
            sline += '&' + summary[camera][item]['units']
            # Description
            sline += '&' + summary[camera][item]['desc']

            tbfil.write(sline + '\\\\ \n')

    # End end
    tbfil.write('\\hline \n')


    tbfil.write('\\enddata \n')

    #tbfil.write('\\tablenotetext{a}{Spectroscopic redshifts are reported to 4 significant digits.  Photometric to 2.} \n')
    #tbfil.write('the gas below the line $\\rm DEC_{\\rm off} = \\aslope RA_{\\rm off} \\ayoff$}\n')
    #tbfil.write('\\tablecomments{Column 1: FRB source. Columns 2 and 3: R.A. and Decl. of the FRB (J2000). Column 4: FRB error ellipse. Column 5: FRB classication. Repeating = yes(y)/no(n). Column 6 and 7: R.A. and Dec. of the associated host galaxy (J2000). Column 8: projected angular offset of the FRB to the host galaxy center. Column 9: association radius $\delta x$ \citep{Tunnicliffe14}. Column 10: angular effective radius of the host measured from a sersic model using GALFIT \citep{galfit} on the $i$-band images (or equivalent). Column 11: effective search radius \citep{Bloom02}. Column 12: measured apparent magnitude of the host. Column 13: filter used for the magnitude measurement. Column 14: probability of chance coincidence using the \citet{Bloom02} formalism. Column 15: sample designations following the criteria outlined in $\S$~\\ref{ssec:associate}.}\n')
    # End
    tbfil.write('\\end{deluxetable*} \n')

    tbfil.close()
    print('Wrote {:s}'.format(outfile))



#### ########################## #########################
#### ########################## #########################
#### ########################## #########################

# Command line execution
if __name__ == '__main__':

    #mktab_photom('blue', sub=True)
    #mktab_photom('red', sub=True)
    mktab_summary()
