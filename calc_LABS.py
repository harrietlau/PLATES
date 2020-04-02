#! /usr/bin/python

## we calculate all parameters required for plotting

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import ticker,cm

plt.rc('font', family='sans serif')
mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13

dir = '../../../HL_out/HOT/xfit_mxw/'
fout = 'zlab_hot_TK1.dat'


J1 = np.loadtxt(dir+'J1_zf.dat')
J2 = np.loadtxt(dir+'J2_zf.dat')
z = np.loadtxt(dir+'Z_km.dat')
f = np.loadtxt(dir+'freqs.dat')
Qi = J1/J2
eta_ss_vbr = np.loadtxt(dir+'eta_diff_HK.dat')
M0 = np.loadtxt(dir+'Gu_z.dat')
J0 = 1./M0
om = 2.*np.pi*f
nz = len(z)
nf = len(f)
J = J1-1.j*J2
M = 1./J
maxf = 1./(eta_ss_vbr/M0)
tauM = 1./maxf

eta = np.zeros((nz,nf),dtype=complex)
eta_abs = np.zeros((nz,nf))
for i in range(nz):
    eta[i,:] = -1.j*M[i,:]/om
    eta_abs[i,:] = np.abs(eta[i,:])
eta_mx = np.zeros((nz,nf),dtype=complex)
for i in range(nz):
    eta_mx[i,:] = eta_ss_vbr[i]/(1.+1.j*om*eta_ss_vbr[i]/M0[i])
eta_norm = eta_abs/np.abs(eta_mx)

##################################################
################ FIND LITHOSPHERE ################
##################################################

# first time running? run the following and find
# where zmax is at high frequency
#### <---
#zlab = np.zeros((nf))
#for i in range(nf):
#    maxf_trim = maxf
#    ival = np.argmin(np.abs(maxf_trim-f[i]))
#    zlab[i] = z[ival]
#print 'pick largest zlab for zmax at high frequency:'
#print zlab
#### --->

# know zmax? run the followng and comment out above
#### <---
zmax = 110.0
zlab = np.zeros((nf))
for i in range(nf):
    ivals = np.where(z<=zmax)[0]
    maxf_trim = maxf[ivals]
    ival = np.argmin(np.abs(maxf_trim-f[i]))
    zlab[i] = z[ival]
newdata = np.zeros((nf,2))
newdata[:,0] = f
newdata[:,1] = zlab
np.savetxt(fout,newdata,fmt='%e')
#### --->


















