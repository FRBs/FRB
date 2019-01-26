from frb.halos import build_grid

# Command line execution
if __name__ == '__main__':
    #build_grid(outfile='z1_mNFW_10000', ntrial=10000)
    #build_grid(outfile='z1_mNFW_10000_21dec2018', ntrial=10000)
    #build_grid(outfile='test', ntrial=500, r_max=1.)
    #build_grid(outfile='test', ntrial=100, r_max=1.)
    #build_grid(outfile='test', ntrial=10)

    # Fiducial model
    build_grid(outfile='z1_mNFW_10000_rmax1', ntrial=10000, r_max=1., f_hot=0.75)
