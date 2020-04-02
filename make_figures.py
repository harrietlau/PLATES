#! /usr/bin/python

## we calculate all parameters required for plotting

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import ticker,cm

plt.rc('font', family='sans serif')
mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13

dir = '../../../HL_out/COLD/xfit_mxw/'

J1 = np.loadtxt(dir+'J1_zf.dat')
J2 = np.loadtxt(dir+'J2_zf.dat')
z = np.loadtxt(dir+'Z_km.dat')
f = np.loadtxt(dir+'freqs.dat')
Q = J1/J2
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

##################################################
### CACULATE NEW J2 WITH DISLOCATION VISCOSITY ###
##################################################

# for Takei's model:
# J2 = Ju * (B + p/(2pi*tauM))      
# J2_tot = Ju * (B + p/(2pi*tauM_tot) )
# where p = 1/f, tauM = Maxwell time from VBR, 
#       tauM_tot = Maxwell time with dislocation 
#       creep

# for Extended Burgers model:
# J2 = Ju * (w*A + 1/(w*tauM))
# J2_tot = Ju * (w*A + 1/(w*tauM_tot) )
# where w = angular frequency, tauM = Maxwell time from diff
#       tuaM_tot = Maxwell time with disl
#       A = A(w)
#       A = (J2/Ju - 1/(w*tauM))/w

## do this for each MPa.
eta_ss_tot = np.loadtxt(dir+'eta_ss_sigvar/eta_tot_HK03.dat')
sig = np.loadtxt(dir+'eta_ss_sigvar/sig_vec_MPa.dat')
nsig = len(sig)
J_tot_sig = np.zeros((nz,nf,nsig),dtype=complex)
M_tot_sig = np.zeros((nz,nf,nsig),dtype=complex)
maxf_tot = np.zeros((nz,nsig),dtype=complex)

for isig in range(nsig):

    ##### FOR TAKEI #####
    ## first calculate B
    p = 1./f
    B = np.zeros((nz,nf))
    for i in range(nz):
        B[i,:] = J2[i,:]/J0[i] - p/(2.*np.pi*tauM[i])
    tauM_tot = eta_ss_tot[:,isig]/M0 
    maxf_tot[:,isig] = 1./tauM_tot
    J2_tot = np.zeros((nz,nf))
    for i in range(nz):
        J2_tot[i,:] = J0[i]*(B[i,:]+p/(2.*np.pi*tauM_tot[i]))
    #####################

    ###### FOR J&F ######
    #A = np.zeros((nz,nf))
    #for i in range(nz):
    #    A[i,:] = (J2[i,:]/J0[i] - 1./(om*tauM[i]))/om 
    #tauM_tot = eta_ss_tot[:,isig]/M0
    #maxf_tot[:,isig] = 1./tauM_tot
    #J2_tot = np.zeros((nz,nf))
    #for i in range(nz):
    #    J2_tot[i,:] = J0[i]*(om*A[i,:] + 1./(om*tauM_tot[i]))
    #####################

    J_tot_sig[:,:,isig] = J1-1.j*J2_tot
    M_tot_sig[:,:,isig] = 1./J_tot_sig[:,:,isig]

ipick = 6 # where sig = 1MPa
grays = np.linspace(0.2,0.8,nsig)

##################################################
######### CALCULATE COMPLEX VISCOSITIES ##########
##################################################

eta = np.zeros((nz,nf),dtype=complex)
eta_abs = np.zeros((nz,nf))
for i in range(nz):
    eta[i,:] = -1.j*M[i,:]/om
    eta_abs[i,:] = np.abs(eta[i,:])
eta_mx = np.zeros((nz,nf),dtype=complex)
for i in range(nz):
    eta_mx[i,:] = eta_ss_vbr[i]/(1.+1.j*om*eta_ss_vbr[i]/M0[i])
eta_norm = eta_abs/np.abs(eta_mx)

# With stress dependent dislocation creep:
eta_tot = np.zeros((nz,nf,nsig),dtype=complex)
eta_tot_abs = np.zeros((nz,nf,nsig))
eta_tot_norm = np.zeros((nz,nf,nsig))
for isig in range(nsig):
    for i in range(nz):
        eta_tot[i,:,isig] = -1.j*M_tot_sig[i,:,isig]/om
        eta_tot_abs[i,:,isig] = np.abs(eta_tot[i,:,isig])
    eta_tot_norm[:,:,isig] = eta_tot_abs[:,:,isig]/np.abs(eta_mx)

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
zmax = 240.0
zlab = np.zeros((nf))
for i in range(nf):
    ivals = np.where(z<=zmax)[0]
    maxf_trim = maxf[ivals]
    ival = np.argmin(np.abs(maxf_trim-f[i]))
    zlab[i] = z[ival]
newdata = np.zeros((nf,2))
newdata[:,0] = f
newdata[:,1] = zlab
np.savetxt('zlab.dat',newdata,fmt='%e')
#### --->

##################################################
################### PLOT FIGURES #################
##################################################

##### FIGURE 1 ######
# Temperature Profile
[P,T,rho1] = np.loadtxt('cold_tkmxw.dat',unpack=True)
i1200 = np.argmin(np.abs(T-(1100)))
i1350 = np.argmin(np.abs(T-(1425)))


plt.figure(1,figsize=(2.4,4.5))
z1 = np.linspace(0,350,len(T))
ax1=plt.gca()
ax1.plot(T,z1,'k-')
ax1.set_yticks([0,50,100,150,200,250,300])
ax1.plot([0,1600],[z[i1200],z[i1200]],':',color='tab:red')
ax1.plot([0,1600],[z[i1350],z[i1350]],':',color='tab:red')
ax1.set_ylim([0,350])
ax1.set_ylabel(r'$z$ (km)')
ax1.set_xlim([0,1500])
ax1.set_xticks([0,800,1600])
ax1.set_xlabel(r'$T$ ($^\circ$C)')
ax1.invert_yaxis()

ax2 = ax1.twiny()
ax2.set_xticks([-2,-6,-8,-10,-13])
ax2.set_xticklabels(['hr','week','y',r'10$^3$ y',r'10$^6$ y'],fontsize=11)
ax2.set_xlabel('time scale')
ax2.set_xlim([-15,0])

plt.tight_layout()
plt.savefig('fig1_a.pdf')


##### FIGURE 2 ######
# CONTOUR PLOT WITH 1-MINUS-eta-star-bar

plt.figure(2,figsize=(5,4.5))
levels = np.arange(0,0.55,0.05)
val=1.-eta_norm
[r,c]=np.where(val>np.max(levels))
val[r,c]=np.max(levels)
[r,c]=np.where(val<np.min(levels))
val[r,c]=np.min(levels)
ff,zz=np.meshgrid(np.log10(f),z)
ax1=plt.gca()
ax1.contourf(ff,zz,val,levels=levels,cmap='summer')
ax1.set_ylim([0,350])
ax1.set_yticks([0,50,100,150,200,250,300])
ax1.set_ylabel(r'$z$ (km)')
ax1.set_xticks([-15,-10,-5,0])
ax1.set_xticklabels(['$10^{-15}$','$10^{-10}$','$10^{-5}$','$10^{0}$'])
ax1.set_xlabel(r'$f$ (Hz)')
ax1.set_xlim([-15,0])
ax1.invert_yaxis()

ax2 = ax1.twiny()
ax2.plot(np.log10(f),zlab,color='orangered',lw=6)
ax2.plot(np.log10(maxf),z,'o',mew=1,mec='k',ms=5,mfc='orangered',zorder=10)
ax2.set_xticks([-2,-6,-8,-10,-13])
ax2.set_xticklabels(['hr','week','y',r'10$^3$ y',r'10$^6$ y'],fontsize=11)
ax2.set_xlabel('time scale')
ax2.set_xlim([-15,0])
plt.tight_layout()

plt.savefig('fig1_b.pdf')

##### FIGURE ** #####
# DUMMY FIGURE FOR COLORBAR 

plt.figure(15,figsize=(5,5))
plt.contourf(ff,zz,val,levels=levels,cmap='summer')
plt.yticks([0,50,100,150,200,250,300])
plt.ylabel(r'$z$ (km)')
cbar=plt.colorbar(orientation='horizontal')
cbar.set_label(r'$1-\bar{\eta}^*$')
plt.tight_layout()

plt.savefig('fig1_a_cb.pdf')

##### FIGURE 3 #####
### DETERMINE VISCOSITY-DEFINED PLATE THICKNESS
imin = np.argmin(eta_ss_vbr)
visc_adiabat = np.mean(eta_ss_vbr[imin:-1])
visc_lab = 10.*visc_adiabat
iilab = np.argmin( np.abs(eta_ss_vbr-visc_lab) )


# VISCOSITY WITH DEPTH
ichoose = [0,9,17,25,33] #1e-15,1e-13,1e-11,1e-9,1e-7
nchoose = len(ichoose)
colorsz2 = plt.cm.Oranges_r(np.linspace(0.2,0.8,nchoose))
plt.figure(3,figsize=(2.8,4.5))
for isig in range(nsig):
    if (isig==ipick):
        print 'do nothing'
    else:
        plt.semilogx(eta_ss_tot[:,isig],z,'-',color=str(grays[isig]),lw=1)
        
for i in range(nchoose):
    iz = np.argmin(np.abs(z-zlab[ichoose[i]]))
    plt.semilogx(eta_abs[iz:nz,ichoose[i]],z[iz:nz],color=colorsz2[i],lw=2)
    plt.semilogx([eta_abs[iz,ichoose[i]],eta_abs[iz,ichoose[i]]],[z[iz],z[iz]],'o',mew=0,mfc=colorsz2[i],ms=8)

plt.semilogx(eta_ss_vbr,z,color='k',lw=2)
plt.semilogx(eta_ss_tot[:,ipick],z,'--',color='k',lw=2)
plt.semilogx([1e17,1e26],[z[i1200],z[i1200]],':',color='tab:red')
plt.semilogx([1e17,1e26],[z[i1350],z[i1350]],':',color='tab:red')
plt.semilogx([1e17,1e26],[z[iilab],z[iilab]],'-',color='orangered',lw=3)
plt.yticks([0,50,100,150,200,250,300])
plt.ylim([0,350])
plt.ylabel(r'$z$ (km)')
plt.xlim([1e17,1e26])
plt.xticks([1e17,1e20,1e23,1e26])
plt.xlabel(r'$\|\eta^*\|$ (Pa s)')
plt.gca().invert_yaxis()
plt.tight_layout()

plt.savefig('fig1_c.pdf')

##### FIGURE 4 ######
# Q with depth

ichoose = [-10,-9,-8,-5,-1] #200s,100s,50s,10s,1s
nchoose = len(ichoose)

# LAB with Q criterion
imin = np.argmin(Q[:,ichoose[3]])
Q_adiabat = np.mean(Q[imin:-1,ichoose[3]])
Q_lab = 10.*Q_adiabat
iilab = np.argmin( np.abs(Q[:,ichoose[3]]-Q_lab) )

colorsz2 = plt.cm.Greens(np.linspace(0.2,0.8,nchoose))
plt.figure(4,figsize=(2.4,4.5))
for i in range(nchoose):
    plt.plot(Q[:,ichoose[i]],z,color=colorsz2[i],lw=2)
plt.plot([0,900],[z[i1200],z[i1200]],':',color='tab:red')
plt.plot([0,900],[z[i1350],z[i1350]],':',color='tab:red')
plt.plot([0,900],[z[iilab],z[iilab]],'-',color='orangered',lw=3)
plt.yticks([0,50,100,150,200,250,300])
plt.ylim([0,350])
plt.ylabel(r'$z$ (km)')
plt.xlim([0,900])
plt.xticks([0,300,600,900])
plt.xlabel(r'$Q$')
plt.gca().invert_yaxis()
plt.tight_layout()

plt.savefig('fig1_d.pdf')

##### FIGURE 5 ######
# vs with depth
rho = np.concatenate(([rho1[0]],rho1))
vs = np.zeros(np.shape(M))
for i in range(nf):
    vs[:,i] = np.sqrt(np.real(M[:,i])/rho)

# LAB with vs criterion (maximum negative gradient - Fischer et al., 2010)
vgrad = np.zeros(len(z))
for i in range(len(z)-1):
    vgrad[i] = (vs[i+1,ichoose[3]]-vs[i,ichoose[3]])/(z[i+1]-z[i])
vgrad[-1] = vgrad[-2]
# erase top bit:
vgrad[0:15] = 100.
iilab = np.argmin(vgrad)

plt.figure(5,figsize=(2.4,4.5))
for i in range(nchoose):
    plt.plot(vs[:,ichoose[i]]/1000.,z,color=colorsz2[i],lw=2)
plt.plot([4.6,4.9],[z[i1200],z[i1200]],':',color='tab:red')
plt.plot([4.6,4.9],[z[i1350],z[i1350]],':',color='tab:red')
plt.plot([4.6,4.9],[z[iilab],z[iilab]],'-',color='orangered',lw=3)
plt.yticks([0,50,100,150,200,250,300])
plt.ylim([0,350])
plt.ylabel(r'$z$ (km)')
plt.xlim([4.6,4.9])
plt.xticks([4.6,4.7,4.8,4.9])
plt.xlabel(r'$v_{\rm s}$ (km/s)')
plt.gca().invert_yaxis()
plt.tight_layout()

plt.savefig('fig1_e.pdf')

##### FIGURE ** #####
# DUMMY FIGURE FOR COLORBAR FOR ETA
plt.figure(16)
val = np.zeros((nz,nf))
val[:,0] = 0.2
val[:,1:8] = 0.8
levels=np.arange(0,1.2,0.2)
plt.contourf(ff,zz,val,levels=levels,cmap=cm.Oranges_r)
cbar=plt.colorbar(orientation='horizontal',ticks=[-0.9,-0.7,-0.5,-.3,-.1,.1,0.3,0.5,0.7,0.9,1.1])
cbar.ax.set_xticklabels([r'$10^{15}$',r'$10^{13}$',r'$10^{11}$',r'$10^{9}$',r'$10^{7}$'],\
    fontsize=20)
cbar.set_label(r'$\tau$ (s)',fontsize=20)
plt.tight_layout()
plt.savefig('fig1_c_cb.pdf')

##### FIGURE ** #####
# DUMMY FIGURE FOR COLORBAR FOR Q and Vs
plt.figure(17)
val = np.zeros((nz,nf))
val[:,0] = 0.2
val[:,1:8] = 0.8
levels=np.arange(0,1.2,0.2)
plt.contourf(ff,zz,val,levels=levels,cmap=cm.Greens)
cbar=plt.colorbar(orientation='horizontal',ticks=[-0.9,-0.7,-0.5,-.3,-.1,.1,0.3,0.5,0.7,0.9,1.1],)
cbar.ax.set_xticklabels([r'200',r'100',r'50',r'10',r'1'],fontsize=20)
cbar.set_label(r'$\tau$ (s)',fontsize=20)
plt.tight_layout()
plt.savefig('fig1_de_cb.pdf')


##### FIGURE 6 ######
# ETA STAR BAR WITH DEPTH
colorsz1 = plt.cm.winter(np.linspace(0,1,nz))
plt.figure(6,figsize=(3,8.5))
r = z[-1]-z
fac = 15.
for i in range(0,nz,6):

    plt.semilogx(f,np.ones(nf)*fac+r[i]-fac,'k:',zorder=0)

    plt.semilogx(f,eta_norm[i,:]*fac+r[i]-fac,color=colorsz1[i],\
                 lw=2,zorder=1)
    ival = np.argmin(np.abs(maxf[i]-f))
    dr = eta_norm[i,ival]*fac-fac
    plt.plot([maxf[i],maxf[i]],[r[i]+dr,r[i]+dr],'o',color='orangered',\
             mew=0,ms=8,zorder=10)

    plt.semilogx(f,eta_tot_norm[i,:,ipick]*fac+r[i]-fac,'--',color=colorsz1[i],\
                 lw=2,zorder=1)
    ival = np.argmin(np.abs(maxf_tot[i,ipick]-f))
    dr = eta_tot_norm[i,ival,ipick]*fac-fac
    plt.plot([maxf_tot[i,ipick],maxf_tot[i,ipick]],[r[i]+dr,r[i]+dr],'^',color='orangered',\
             mew=0,ms=8,zorder=10)

plt.xlabel(r'$f$ (Hz)')
plt.ylabel(r'$z$ (km)')
plt.xlim([1.e-15,1.e-0])
plt.xticks([1e-15,1e-10,1e-5,1e-0])
plt.yticks([50,100,150,200,250,300,350])
plt.gca().set_yticklabels(['300','250','200','150','100','50','0'])
plt.ylim([0,350])
plt.tight_layout()
colorsz1 = plt.cm.winter(np.linspace(0,1,nz))
plt.figure(6,figsize=(3,8.5))
r = z[-1]-z
fac = 15.
for i in range(0,nz,6):

    plt.semilogx(f,np.ones(nf)*fac+r[i]-fac,'k:',zorder=0)

    plt.semilogx(f,eta_norm[i,:]*fac+r[i]-fac,color=colorsz1[i],\
                 lw=2,zorder=1)
    ival = np.argmin(np.abs(maxf[i]-f))
    dr = eta_norm[i,ival]*fac-fac
    plt.plot([maxf[i],maxf[i]],[r[i]+dr,r[i]+dr],'o',color='orangered',\
             mew=0,ms=8,zorder=10)

    plt.semilogx(f,eta_tot_norm[i,:,ipick]*fac+r[i]-fac,'--',color=colorsz1[i],\
                 lw=2,zorder=1)
    ival = np.argmin(np.abs(maxf_tot[i,ipick]-f))
    dr = eta_tot_norm[i,ival,ipick]*fac-fac
    plt.plot([maxf_tot[i,ipick],maxf_tot[i,ipick]],[r[i]+dr,r[i]+dr],'^',color='orangered',\
             mew=0,ms=8,zorder=10)

plt.xlabel(r'$f$ (Hz)')
plt.ylabel(r'$z$ (km)')
plt.xlim([1.e-15,1.e-0])
plt.xticks([1e-15,1e-10,1e-5,1e-0])
plt.yticks([50,100,150,200,250,300,350])
plt.gca().set_yticklabels(['300','250','200','150','100','50','0'])
plt.ylim([0,350])
plt.tight_layout()
plt.savefig('fig1_f.pdf')


##### FIGURE 7 #####
# INSET OF FIG 6
plt.figure(7,figsize=(3,2))
idx = 81 # 250 km
for isig in range(nsig):
    if (isig==ipick):
        print 'do nothing'
    else: 
        plt.semilogx(f,eta_tot_norm[idx,:,isig],'-',color=str(grays[isig]))
        ival = np.argmin(np.abs(maxf_tot[idx,isig]-f))
        plt.semilogx([maxf_tot[idx,isig],maxf_tot[idx,isig]],\
            [eta_tot_norm[idx,ival,isig],eta_tot_norm[idx,ival,isig]],\
            '^',color=str(grays[isig]),mew=0,ms=8,zorder=10)
plt.semilogx(f,eta_norm[idx,:],color=colorsz1[idx],lw=2)
ival = np.argmin(np.abs(maxf[idx]-f))
plt.semilogx([maxf[idx],maxf[idx]],[eta_norm[idx,ival],eta_norm[idx,ival]],\
	'o',color='orangered',mew=0,ms=8,zorder=10)
plt.semilogx(f,eta_tot_norm[idx,:,ipick],'--',color=colorsz1[idx],lw=2)
ival = np.argmin(np.abs(maxf_tot[idx,ipick]-f))
plt.semilogx([maxf_tot[idx,ipick],maxf_tot[idx,ipick]],\
    [eta_tot_norm[idx,ival,ipick],eta_tot_norm[idx,ival,ipick]],\
    '^',color='orangered',mew=0,ms=8,zorder=10)

plt.xlim([1.e-15,1e0])
plt.xticks([1.e-15,1e-10,1e-5,1e0])
plt.ylim([0,1.1])
plt.yticks([0.5,1.0])
plt.gca().yaxis.set_label_position("right")
plt.gca().yaxis.tick_right()
plt.xlabel(r'$f$ (Hz)')
plt.ylabel(r'$\bar{\eta}^*$')
plt.tight_layout()
plt.savefig('fig1_g.pdf')


##### FIGURE 8 #####
# INSET OF FIG 6
plt.figure(8,figsize=(3,2))
idx = 99 # 300 km
for isig in range(nsig):
    if (isig==ipick):
        print 'do nothing'
    else:
        plt.semilogx(f,eta_tot_norm[idx,:,isig],'-',color=str(grays[isig]))
        ival = np.argmin(np.abs(maxf_tot[idx,isig]-f))
        plt.semilogx([maxf_tot[idx,isig],maxf_tot[idx,isig]],\
            [eta_tot_norm[idx,ival,isig],eta_tot_norm[idx,ival,isig]],\
            '^',color=str(grays[isig]),mew=0,ms=8,zorder=10)
plt.semilogx(f,eta_norm[idx,:],color=colorsz1[idx],lw=2)
ival = np.argmin(np.abs(maxf[idx]-f))
plt.semilogx([maxf[idx],maxf[idx]],[eta_norm[idx,ival],eta_norm[idx,ival]],\
        'o',color='orangered',mew=0,ms=8,zorder=10)
plt.semilogx(f,eta_tot_norm[idx,:,ipick],'--',color=colorsz1[idx],lw=2)
ival = np.argmin(np.abs(maxf_tot[idx,ipick]-f))
plt.semilogx([maxf_tot[idx,ipick],maxf_tot[idx,ipick]],\
    [eta_tot_norm[idx,ival,ipick],eta_tot_norm[idx,ival,ipick]],\
    '^',color='orangered',mew=0,ms=8,zorder=10)

plt.xlim([1.e-15,1e0])
plt.xticks([1.e-15,1e-10,1e-5,1e0])
plt.ylim([0,1.1])
plt.yticks([0.5,1.0])
plt.gca().yaxis.set_label_position("right")
plt.gca().yaxis.tick_right()
plt.xlabel(r'$f$ (Hz)')
plt.ylabel(r'$\bar{\eta}^*$')
plt.tight_layout()

plt.savefig('fig1_h.pdf')

##### FIGURE 9 #####
# INSET OF FIG 6
plt.figure(9,figsize=(3,2))
idx = nz-1 # 350 km
for isig in range(nsig):
    if (isig==ipick):
        print 'do nothing'
    else:
        plt.semilogx(f,eta_tot_norm[idx,:,isig],'-',color=str(grays[isig]))
        ival = np.argmin(np.abs(maxf_tot[idx,isig]-f))
        plt.semilogx([maxf_tot[idx,isig],maxf_tot[idx,isig]],\
            [eta_tot_norm[idx,ival,isig],eta_tot_norm[idx,ival,isig]],\
            '^',color=str(grays[isig]),mew=0,ms=8,zorder=10)
plt.semilogx(f,eta_norm[idx,:],color=colorsz1[idx],lw=2)
ival = np.argmin(np.abs(maxf[idx]-f))
plt.semilogx([maxf[idx],maxf[idx]],[eta_norm[idx,ival],eta_norm[idx,ival]],\
        'o',color='orangered',mew=0,ms=8,zorder=10)
plt.semilogx(f,eta_tot_norm[idx,:,ipick],'--',color=colorsz1[idx],lw=2)
ival = np.argmin(np.abs(maxf_tot[idx,ipick]-f))
plt.semilogx([maxf_tot[idx,ipick],maxf_tot[idx,ipick]],\
    [eta_tot_norm[idx,ival,ipick],eta_tot_norm[idx,ival,ipick]],\
    '^',color='orangered',mew=0,ms=8,zorder=10)

plt.xlim([1.e-15,1e0])
plt.xticks([1.e-15,1e-10,1e-5,1e0])
plt.ylim([0,1.1])
plt.yticks([0.5,1.0])
plt.gca().yaxis.set_label_position("right")
plt.gca().yaxis.tick_right()
plt.xlabel(r'$f$ (Hz)')
plt.ylabel(r'$\bar{\eta}^*$')
plt.tight_layout()

plt.savefig('fig1_i.pdf')

##### FIGURE 10 #####
# LEGEND FOR INSETS
plt.figure(10,figsize=(3,7))
plt.semilogx(f,eta_norm[idx,:],'k-',lw=2)
for isig in range(nsig):
    if (isig==ipick):
        plt.semilogx(f,eta_norm[idx,:],'k--',lw=2)
    else:
        plt.semilogx(f,eta_tot_norm[idx,:,isig],'-',color=str(grays[isig]))

plt.xlim([1.e-15,1e0])
plt.xticks([1.e-15,1e-10,1e-5,1e0])
plt.ylim([0,6.1])
plt.legend(['0.0','.01','.02','.05','0.1','0.2','0.5','1.0','2.0','5.0','10'],\
    loc='upper right')


plt.tight_layout()
plt.savefig('fig1_i_legend.pdf')


plt.show()




















