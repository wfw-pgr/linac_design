import numpy as np

# ========================================================= #
# ===  calculation of einstein Beta Routine             === #
# ========================================================= #

def einsteinBeta():

    # ------------------------------------------------- #
    # --- [1] constants                             --- #
    # ------------------------------------------------- #
    
    mp     =  1.6726219e-27
    me     =  9.1093837e-31
    qp     =  1.6021766e-19
    qe     = -1.6021766e-19
    cv     =  3.0e8
    Ek     =  40.0e6

    # ------------------------------------------------- #
    # --- [2] velocity calculation                  --- #
    # ------------------------------------------------- #

    Th_p   = np.abs( qp ) * Ek / ( mp * cv**2 )
    Th_e   = np.abs( qe ) * Ek / ( me * cv**2 )
    
    beta_p = np.sqrt( 1.0 - 1.0 / ( 1.0 + Th_p )**2 )
    velo_p = beta_p * cv
    
    beta_e = np.sqrt( 1.0 - 1.0 / ( 1.0 + Th_e )**2 )
    velo_e = beta_e * cv

    
    # ------------------------------------------------- #
    # --- [3] print out                             --- #
    # ------------------------------------------------- #

    print()
    print( "[einsteinBeta.py] calculation of beta & velocity of {0:20} (eV) proton.".format( Ek ) )
    print( "[einsteinBeta.py]          Thermal Energy :: {0:20}".format( Th_p   ) )
    print( "[einsteinBeta.py]          beta           :: {0:20}".format( beta_p ) )
    print( "[einsteinBeta.py]          velocity (m/s) :: {0:20}".format( velo_p ) )
    print()
    print( "[einsteinBeta.py] calculation of beta & velocity of {0:20} (eV) electron.".format( Ek ) )
    print( "[einsteinBeta.py]          Thermal Energy :: {0:20}".format( Th_e   ) )
    print( "[einsteinBeta.py]          beta           :: {0:20}".format( beta_e ) )
    print( "[einsteinBeta.py]          velocity (m/s) :: {0:20}".format( velo_e ) )
    print()


# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):

    einsteinBeta()
