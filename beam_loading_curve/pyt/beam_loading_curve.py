import numpy as np

# ========================================================= #
# ===  beam_loading_curve                               === #
# ========================================================= #

def beam_loading_curve():

    # ------------------------------------------------- #
    # --- [1] make beam loading curve               --- #
    # ------------------------------------------------- #
    import nkUtilities.load__constants as lcn
    cnsFile    = "dat/parameter.conf"
    const      = lcn.load__constants( inpFile=cnsFile )
    Ibeam      = np.linspace( const["Ibeam_Min"], const["Ibeam_Max"], const["Ibeam_num"] )
    P0_RFs     = np.linspace( const["P0_RF_Min"], const["P0_RF_Max"], const["P0_RF_num"] )
    BeamCurve  = np.zeros( (const["Ibeam_num"],const["P0_RF_num"]) )
    tau_attenuation = const["alpha_attenuation"] * const["length_cavity"]
    
    for ik,P0_RF in enumerate( P0_RFs ):
        term1  = np.sqrt( const["shunt_impedance"]*const["length_cavity"]*P0_RF*( 1.0 - np.exp( -2.0*tau_attenuation ) ) )
        term2  = Ibeam*const["shunt_impedance"]*const["length_cavity"] / 2.0 * ( 1.0 - 2.0*tau_attenuation*np.exp( -2.0*tau_attenuation ) / ( 1.0 - np.exp( -2.0*tau_attenuation ) ) )
        BeamCurve[:,ik] = ( term1 - term2 )

    outFile    = "dat/beam_loading_curve.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=BeamCurve )


    # ------------------------------------------------- #
    # --- [2] config settings                       --- #
    # ------------------------------------------------- #
    import nkUtilities.load__config as lcf
    config                   = lcf.load__config()
    config["xTitle"]         = "Beam Current (mA)"
    config["yTitle"]         = "Beam Energy (MeV)"
    config["plt_xAutoRange"] = False
    config["plt_yAutoRange"] = False
    config["plt_xRange"]     = [ 0.0, 300  ]
    config["plt_yRange"]     = [ 0.0, 30.0 ]
    pngFile                  = "png/beam_loading_curve.png"

    
    # ------------------------------------------------- #
    # --- [3] plot                                  --- #
    # ------------------------------------------------- #
    Ibeam     = Ibeam     / const["current_unit"]
    BeamCurve = BeamCurve / const["energy_unit"]
    import nkUtilities.plot1D as pl1
    fig = pl1.plot1D( config=config, pngFile=pngFile )
    for ik in range( const["P0_RF_num"] ):
        fig.add__plot( xAxis=Ibeam, yAxis=BeamCurve[:,ik] )
    fig.set__axis()
    fig.save__figure()


    

# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    beam_loading_curve()
