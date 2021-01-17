import numpy as np


# ========================================================= #
# ===  generate sample wave for travelling wave         === #
# ========================================================= #

def generate__samplewave():

    # ------------------------------------------------- #
    # --- [1] Load Config                           --- #
    # ------------------------------------------------- #

    import nkUtilities.load__constants as lcn
    cnsFile  = "dat/parameter.conf"
    const = lcn.load__constants( inpFile=cnsFile )

    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = const["x1MinMaxNum"]
    x2MinMaxNum = const["x2MinMaxNum"]
    x3MinMaxNum = const["x3MinMaxNum"]
    grid        = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "point" )
    kx, ky      = const["sample_kx"], const["sample_ky"]
    wave1       = np.cos( grid[:,0] * 2.0*np.pi*kx ) + np.sin( grid[:,1] *2.0*np.pi*ky )
    wave2       = np.sin( grid[:,0] * 2.0*np.pi*kx ) + np.cos( grid[:,1] *2.0*np.pi*ky )
    wave1_      = np.concatenate( (grid,wave1[:,None]), axis=1 )
    wave2_      = np.concatenate( (grid,wave2[:,None]), axis=1 )
    size        = (int(x3MinMaxNum[2]),int(x2MinMaxNum[2]),int(x1MinMaxNum[2]),4)
    wave1       = np.reshape( wave1_, size )
    wave2       = np.reshape( wave2_, size )
    
    # ------------------------------------------------- #
    # --- [2] save in File                          --- #
    # ------------------------------------------------- #

    import nkUtilities.save__pointFile as spf
    wavFile1   = "dat/sample_wave1.dat"
    wavFile2   = "dat/sample_wave2.dat"
    spf.save__pointFile( outFile=wavFile1, Data=wave1 )
    spf.save__pointFile( outFile=wavFile2, Data=wave2 )
    
    # ------------------------------------------------- #
    # --- [3] display eigen mode                    --- #
    # ------------------------------------------------- #

    import nkUtilities.cMapTri as cmt
    pngFile1    = "png/sample_wave1.png"
    pngFile2    = "png/sample_wave2.png"
    cmt.cMapTri( xAxis=wave1_[...,0], yAxis=wave1_[...,1], cMap=wave1_[...,3], pngFile=pngFile1 )
    cmt.cMapTri( xAxis=wave2_[...,0], yAxis=wave2_[...,1], cMap=wave2_[...,3], pngFile=pngFile2 )

    return()


# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    generate__samplewave()
