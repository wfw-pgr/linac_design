import numpy                     as np
import nkUtilities.equiSpaceGrid as esg

# ------------------------------------------------- #
# --- [1] parameters                            --- #
# ------------------------------------------------- #

mmode       = 0
nmode       = 1
wg_a        = 0.10
tMax        = 0.0
LI          = 51
LJ          = 51
LK          =  1
LT          =  1
omega       = 2.856e8 * 2.0 * np.pi
mu          = 4.0*np.pi*1e-7
epsilon     = 8.85e-12
cv          = 3.0e8
Emn         = 1.0
rho_mn      = 2.405

# ------------------------------------------------- #
# --- [2] grid making                           --- #
# ------------------------------------------------- #

x1MinMaxNum = [ -wg_a, wg_a, LI ]
x2MinMaxNum = [ -wg_a, wg_a, LJ ]
x3MinMaxNum = [ 0.0  , 0.0 , LK ]
grid        = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                 x3MinMaxNum=x3MinMaxNum, returnType = "point" )
radii       = np.sqrt( grid[:,0]**2 + grid[:,1]**2 )
phi         = np.arctan2( grid[:,1], grid[:,0] )
index       = np.where( radii <= wg_a )
radii       = radii[index]
phi         =   phi[index]
xg          = ( grid[:,0] )[index]
yg          = ( grid[:,1] )[index]
zg          = ( grid[:,2] )[index]

kc          = rho_mn / wg_a
kcr         = kc     * radii
mphi        = mmode * phi
gamma       = np.sqrt( kc**2 - omega**2 / cv**2 )

# ------------------------------------------------- #
# --- [3] TE mode wave                          --- #
# ------------------------------------------------- #

import scipy.special         as special
import nkBasicAlgs.robustInv as inv
rinv        = inv.robustInv( radii )
Ez_Amp      =                                     Emn * special.jv( mmode  , kcr ) * np.cos( mphi )
Hz_Amp      = 0.0 * kcr
Er_Amp      = -           gamma         / kc    * Emn * special.jv( mmode+1, kcr ) * np.cos( mphi )
Ep_Amp      =             gamma * mmode / kc**2 * Emn * special.jv( mmode  , kcr ) * np.sin( mphi ) * rinv
Hr_Amp      =   epsilon * omega * mmode / kc**2 * Emn * special.jv( mmode  , kcr ) * np.sin( mphi ) * rinv
Hp_Amp      = - epsilon * omega         / kc    * Emn * special.jv( mmode+1, kcr ) * np.cos( mphi )


# ------------------------------------------------- #
# --- [4] colormap of the Amplitude             --- #
# ------------------------------------------------- #

import nkUtilities.cMapTri as cmt
cmt.cMapTri( xAxis=xg, yAxis=yg, cMap=Er_Amp, pngFile="Er_Amp.png" )
cmt.cMapTri( xAxis=xg, yAxis=yg, cMap=Ep_Amp, pngFile="Ep_Amp.png" )
cmt.cMapTri( xAxis=xg, yAxis=yg, cMap=Ez_Amp, pngFile="Ez_Amp.png" )
cmt.cMapTri( xAxis=xg, yAxis=yg, cMap=Hr_Amp, pngFile="Hr_Amp.png" )
cmt.cMapTri( xAxis=xg, yAxis=yg, cMap=Hp_Amp, pngFile="Hp_Amp.png" )
cmt.cMapTri( xAxis=xg, yAxis=yg, cMap=Hz_Amp, pngFile="Hz_Amp.png" )
