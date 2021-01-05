import sys
import numpy                      as np
import nkUtilities.load__config   as lcf
import nkUtilities.cMapTri        as cmt
import nkUtilities.configSettings as cfs


# ========================================================= #
# ===  display__sf7                                     === #
# ========================================================= #
def display__sf7():
    
    # ------------------------------------------------- #
    # --- [1] Arguments                             --- #
    # ------------------------------------------------- #
    import nkUtilities.load__constants as lcn
    cnfFile = "dat/parameter.conf"
    const   = lcn.load__constants( inpFile=cnfFile )

    config  = lcf.load__config()
    datFile = const["efdFile"]
    pngFile = const["efdFile"].replace( "dat", "png" )
    pngFile = pngFile.replace( ".png", "_{0}.png" )

    # ------------------------------------------------- #
    # --- [2] Fetch Data                            --- #
    # ------------------------------------------------- #
    import nkUtilities.load__pointFile as lpf
    Data  = lpf.load__pointFile( inpFile=datFile, returnType="point" )
    xAxis = Data[:,0]
    yAxis = Data[:,1]
    Ez    = Data[:,2]
    Er    = Data[:,3]
    Ea    = Data[:,4]
    Hp    = Data[:,5]
    
    # ------------------------------------------------- #
    # --- [3] config Settings                       --- #
    # ------------------------------------------------- #
    cfs.configSettings( configType="cMap_def", config=config )
    config["FigSize"]        = (8,3)
    config["cmp_position"]   = [0.16,0.12,0.97,0.88]
    config["xTitle"]         = "Z (m)"
    config["yTitle"]         = "R (m)"
    config["xMajor_Nticks"]  = 8
    config["yMajor_Nticks"]  = 3
    config["cmp_xAutoRange"] = True
    config["cmp_yAutoRange"] = True
    config["cmp_xRange"]     = [-5.0,+5.0]
    config["cmp_yRange"]     = [-5.0,+5.0]

    # ------------------------------------------------- #
    # --- [4] plot Figure                           --- #
    # ------------------------------------------------- #
    cmt.cMapTri( xAxis=xAxis, yAxis=yAxis, cMap=Ez, pngFile=pngFile.format( "Ez" ), config=config )
    cmt.cMapTri( xAxis=xAxis, yAxis=yAxis, cMap=Er, pngFile=pngFile.format( "Er" ), config=config )
    cmt.cMapTri( xAxis=xAxis, yAxis=yAxis, cMap=Ea, pngFile=pngFile.format( "Ea" ), config=config )
    cmt.cMapTri( xAxis=xAxis, yAxis=yAxis, cMap=Hp, pngFile=pngFile.format( "Hp" ), config=config )
    
    cfs.configSettings( configType="vector_def", config=config )
    cmt.cMapTri( xAxis=xAxis, yAxis=yAxis, cMap=Ea, \
                 uvec =Ez   , vvec =Er   ,
                 pngFile=pngFile.format( "Ea" ), config=config )
    
    

# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    display__sf7()
