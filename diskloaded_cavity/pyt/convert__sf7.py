import os, sys
import numpy as np

# ========================================================= #
# ===  convert__sf7.py                                  === #
# ========================================================= #

def convert__sf7():

    # ------------------------------------------------- #
    # --- [1] load config & sf7 file                --- #
    # ------------------------------------------------- #

    import nkUtilities.load__constants as lcn
    cnfFile = "dat/parameter.conf"
    const   = lcn.load__constants( inpFile=cnfFile )
    
    with open( const["sf7File"], "r" ) as f:
        lines = f.readlines()

    # ------------------------------------------------- #
    # --- [2] search for the data start line        --- #
    # ------------------------------------------------- #

    nLine      = int( const["sf7_xMinMaxNum"][2] * const["sf7_yMinMaxNum"][2] )
    searchline = "Electromagnetic fields for a rectangular area with corners at:"
    offset     = 7
    DataStart  = None
    for iL,line in enumerate(lines):
        if ( line.strip() == searchline ):
            DataStart = iL + offset
            break
    if ( DataStart is None ):
        sys.exit( "[convert__sf7.py] cannot find searchline in {0}".format( const["sf7File"] ) )
        
    # ------------------------------------------------- #
    # --- [3] fetch Data from outsf7.txt            --- #
    # ------------------------------------------------- #
    
    # -- [3-1] load file contents                   --  #
    with open( const["sf7File"], "r" ) as f:
        Data = np.loadtxt( f, skiprows=DataStart, max_rows=nLine )

    # -- [3-2] unit conversion                      --  #
    Data[:,0] = Data[:,0] * 1.e-3 #  Z  :: (mm)   -> (m)
    Data[:,1] = Data[:,1] * 1.e-3 #  R  :: (mm)   -> (m)
    Data[:,2] = Data[:,2] * 1.e+6 #  Ez :: (MV/m) -> (V/m)
    Data[:,3] = Data[:,3] * 1.e+6 #  Er :: (MV/m) -> (V/m)
    Data[:,4] = Data[:,4] * 1.e+6 # |E| :: (MV/m) -> (V/m)
    Data[:,5] = Data[:,5]         #  H  :: (A/m)
    
    # ------------------------------------------------- #
    # --- [4] save as a pointData                   --- #
    # ------------------------------------------------- #
    
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=const["efdFile"], Data=Data )


# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    convert__sf7()
