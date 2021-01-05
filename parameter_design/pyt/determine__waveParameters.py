import os, sys
import numpy                       as np
import scipy.special               as spc
import nkUtilities.load__constants as lcn


# ========================================================= #
# ===  determine wave parameters                        === #
# ========================================================= #

def determine__waveParameters():

    # ------------------------------------------------- #
    # --- [1] constants                             --- #
    # ------------------------------------------------- #

    prmFile   = "dat/parameter.conf"
    outFile   = "dat/wave_parameters.dat"
    params    = lcn.load__constants( inpFile=prmFile )

    mu0       = 4.0 * np.pi * 1.0e-7
    
    # ------------------------------------------------- #
    # --- [2] calculation of parameters             --- #
    # ------------------------------------------------- #
    if ( str( params["cavity_mode"] ) == "2/3" ):
        params["n_cell_per_wavelength"] = 3.0
        params["cavity_mode_float"]     = 2.0 / 3.0 * np.pi
    else:
        print( "[determine__waveParameters] cavity_mode == ??? " )
        sys.exit()

    # ------------------------------------------------- #
    # --- [3] calculation of the parameters         --- #
    # ------------------------------------------------- #

    # -- [3-1] velocity / wavelength                --  #
    eth                             = np.abs( params["qe"] ) * params["pEnergy"] / ( params["me"] * params["cv"]**2 )
    params["beta_max"]              = np.sqrt( 1.0 - 1.0 / ( 1.0 + eth )**2 )
    params["vph"]                   = params["cv"]  * params["beta_p"]
    params["wavelength_incavity"]   = params["vph"] / params["fres"]
    params["cavity_length"]         = params["wavelength_incavity"] / 3.0
    params["propagation_const"]     = 2.0 * np.pi / params["wavelength_incavity"]

    # -- [3-2] geometry of cavity                   --  #
    params["2b_waveguide_approx"]   = 2.0 * params["bessel_root"] * params["cv"] / ( 2.0 * np.pi * params["fres"] )

    # -- [3-3] cavity parameters                    --  #
    params["kappa_coupling_coef"]   = 4.0 * params["r_disk"]**3 / ( 3.0 * np.pi * ( spc.jv( 1.0, params["bessel_root"] ) )**2 * params["r_cavity"]**2 * ( params["z_cavity"] - params["z_disk"] ) )
    params["xi_attenate_coef"]      = np.sqrt( ( params["bessel_root"] / params["r_disk"] )**2 - ( 2.0 * np.pi * params["fres"] / params["cv"] )**2 )
    params["fres_corrected"]        = params["bessel_root"] * params["cv"] / ( 2.0 * np.pi * params["r_cavity"] ) * np.sqrt( 1.0 + params["kappa_coupling_coef"] * ( 1.0 - np.exp( (-1.0) * params["xi_attenate_coef"] * params["z_disk"] ) * np.cos( params["cavity_mode_float"] ) ) )
    params["skin_depth"]            = np.sqrt( ( 2.0 ) / ( 2.0 * np.pi * params["fres"] * mu0 * params["conductivity"] ) )
    params["Qvalue"]                = ( params["wavelength_incavity"] / params["skin_depth"] ) * ( params["beta_p"] * ( 1.0 - params["z_disk"] / params["z_cavity"] ) ) / ( params["n_cell_per_wavelength"] + 2.61 * params["beta_p"] * ( 1.0 - params["z_disk"] / params["z_cavity"] ) )
    params["beta_g"]                = ( 2.0 * params["bessel_root"] ) / ( 3.0 * np.pi * ( spc.jv( 1.0, params["bessel_root"] ) )**2 ) * ( params["r_disk"] / params["r_cavity"] )**3 * np.sin( params["cavity_mode_float"] ) * np.exp( (-1.0)*params["xi_attenate_coef"] * params["z_disk"] )
    params["attenation_coef"]       = ( 2.0*np.pi * params["fres"] ) / ( 2.0 * params["beta_g"] * params["cv"] * params["Qvalue"] )
    params["surface_resistivity"]   = np.sqrt( np.pi * params["Z0_vacuum"] / ( params["conductivity"] * params["wavelength_incavity"] ) )
    params["transit_time_factor"]   = ( np.sin( params["cavity_mode_float"]*0.5 * ( ( params["z_cavity"] - params["z_disk"] ) / params["z_cavity"] ) ) ) / ( params["cavity_mode_float"]*0.5 )
    params["rsh_shunt_impedance"]   = params["Z0_vacuum"]**2 * params["z_cavity"] * params["transit_time_factor"]**2 * ( spc.jv( 0.0, 2.0*np.pi*params["r_disk"] / params["wavelength_incavity"] ) )**2 / ( np.pi * params["surface_resistivity"] * params["r_cavity"] * ( params["r_cavity"] + params["z_cavity"] - params["z_disk"] ) * ( spc.jv( 1.0, params["bessel_root"] ) )**2 )

    
    # ------------------------------------------------- #
    # --- [4] save as text output                   --- #
    # ------------------------------------------------- #

    with open( outFile, "w" ) as f:
        f.write( "\n" )
        f.write( "[wave_parametrs.py] output text file " + "\n" + "\n" )
        for item,value in params.items():
            f.write( " {0:30}, {1:30}\n".format( item, value ) )
        f.write( "\n" )
    print( "[determine__waveParameters] outFile :: {0} ".format( outFile ) )

    
# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    determine__waveParameters()
