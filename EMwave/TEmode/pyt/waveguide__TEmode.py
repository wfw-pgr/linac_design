import os, sys
import numpy                     as np
import nkUtilities.equiSpaceGrid as esg


# ========================================================= #
# ===  waveguide__TEmode.py                             === #
# ========================================================= #

def waveguide__TEmode( time=0.0, kstep=0 ):

    # ------------------------------------------------- #
    # --- [1] parameters                            --- #
    # ------------------------------------------------- #

    import nkUtilities.load__constants as lcn
    inpFile = "dat/parameter.conf"
    const   = lcn.load__constants( inpFile=inpFile )
    mu      = 4.0*np.pi*1e-7
    cv      = 3.0e8
    omega   = 2.0*np.pi*const["freq"]
    Hmn     = const["Hmn"]

    # ------------------------------------------------- #
    # --- [2] grid making                           --- #
    # ------------------------------------------------- #

    x_,y_,z_    = 0, 1, 2
    x1MinMaxNum = [ 0.0, const["wg_a"], const["LI"] ]
    x2MinMaxNum = [ 0.0, const["wg_b"], const["LJ"] ]
    x3MinMaxNum = [ 0.0, const["wg_c"], const["LK"] ]
    grid        = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "structured" )
    xg          = grid[...,0]
    yg          = grid[...,1]
    zg          = grid[...,2]
    
    k0          = omega / cv
    kx          = const["mmode"] * np.pi / const["wg_a"]
    ky          = const["nmode"] * np.pi / const["wg_b"]
    kc          = np.sqrt( kx**2 + ky**2 )
    beta        = np.sqrt( k0**2 - kc**2 )

    wt          = omega * time
    kxx         = kx * xg
    kyy         = ky * yg
    exp_wt_bz_r = np.cos( wt - beta*zg )
    exp_wt_bz_i = np.sin( wt - beta*zg )

    lambda_g    = 2.0 * np.pi / beta
    
    print( " omega    :: {0}".format( omega    ) )
    print( " k0       :: {0}".format( k0       ) )
    print( " kx       :: {0}".format( kx       ) )
    print( " ky       :: {0}".format( ky       ) )
    print( " kc       :: {0}".format( kc       ) )
    print( " beta     :: {0}".format( beta     ) )
    print( " lambda_g :: {0}".format( lambda_g ) )
    
    # ------------------------------------------------- #
    # --- [3] TE mode wave                          --- #
    # ------------------------------------------------- #
    
    Ez_Re  =  0.0  * exp_wt_bz_r 
    Hz_Re  =                           Hmn * np.cos( kxx ) * np.cos( kyy ) * exp_wt_bz_r
    Ex_Re  = omega * ky * mu / kc**2 * Hmn * np.cos( kxx ) * np.sin( kyy ) * exp_wt_bz_i * (-1)
    Ey_Re  = omega * kx * mu / kc**2 * Hmn * np.sin( kxx ) * np.cos( kyy ) * exp_wt_bz_i
    Hx_Re  =  beta * kx      / kc**2 * Hmn * np.sin( kxx ) * np.cos( kyy ) * exp_wt_bz_i * (-1)
    Hy_Re  =  beta * ky      / kc**2 * Hmn * np.cos( kxx ) * np.sin( kyy ) * exp_wt_bz_i * (-1)
    

    # ------------------------------------------------- #
    # --- [4] colormap of the Amplitude             --- #
    # ------------------------------------------------- #

    import nkBasicAlgs.pileupArray as pil
    Data        = pil.pileupArray( (xg,yg,zg,Ex_Re,Ey_Re,Ez_Re,Hx_Re,Hy_Re,Hz_Re) )
    DataLabel   = ["xg","yz","zg","Ex","Ey","Ez","Hx","Hy","Hz"]
    
    import nkVTKRoutines.vtkDataConverter as vdc
    outFile     = "dat/TEwave_{0:04}.vts".format( kstep )
    cvt1        = vdc.vtkDataConverter( vtkFile=outFile, Data=Data, \
                                        tag="data", DataType="structured", DataLabel=DataLabel )


# ========================================================= #
# ===  time sequence of the TEmode                      === #
# ========================================================= #

def sequencial__TEmode():

    # ------------------------------------------------- #
    # --- [1] load constants                        --- #
    # ------------------------------------------------- #
    import nkUtilities.load__constants as lcn
    inpFile = "dat/parameter.conf"
    const   = lcn.load__constants( inpFile=inpFile )

    # ------------------------------------------------- #
    # --- [2] calculate TE mode at each time        --- #
    # ------------------------------------------------- #
    
    timearr = np.linspace( const["tMin"], const["tMax"], const["LT"] )
    
    for ik in range( const["LT"] ):
        waveguide__TEmode( time=timearr[ik], kstep=ik )
    

        
# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    
    # waveguide__TEmode()

    sequencial__TEmode()
