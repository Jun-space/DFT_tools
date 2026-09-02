#!/usr/bin/env python
# This script was created by zqj
# <Tue May 26 09:05:25 CST 2015>

##################################### NOTES ##################################### 
# 1. Set NBLOCK = 1 in the INCAR, so that all the configuration is wrtten to
# XDATCAR.
################################################################################# 

import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftfreq, fftshift
import time
import scipy
import copy
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import sys

import matplotlib as mpl                                                                                                                                                         
mpl.use('Agg')


# Boltzmann Constant in [eV/K]
kB = 8.617332478E-5
# electron volt in [Joule]
ev = 1.60217733E-19
# Avogadro's Constant
Navogadro = 6.0221412927E23

################################################################################      
class xdatcar:
    """ Python Class for VASP XDATCAR """

    def __init__(self, File=None):
        if File is None:
            self.xdatcar = 'XDATCAR'
        else:
            self.xdatcar = File

        # time step of MD
        self.potim = None
        # mass per type
        self.mtype = None
        self.readoutcar()

        self.TypeName = None
        self.ChemSymb = None
        self.Ntype = None
        self.Nions = None
        self.Nelem = None
        self.Niter = None

        # position in Direct Coordinate
        self.position = None
        # position in Cartesian Coordinate
        self.positionC = None
        # Velocity in Angstrom per Femtosecond
        self.velocity = None
        self.readxdat()

        self.mass_and_name_per_ion()
        # Temperature
        self.Temp = np.zeros(self.Niter-1)
        # Kinetic Energy
        self.Ken = np.zeros(self.Niter-1)
        # Time in femtosecond
        self.Time = np.arange(self.Niter-1) * self.potim
        self.getTemp()

        # Velocity Autocorrelation Function 
        self.VAF = None
        self.VAF2= None
        # Pair Correlation Function
        # self.PCF = None

    def mass_and_name_per_ion(self):
        # mass per ion
        self.mions = []
        self.ChemSymb = []

        if self.TypeName is None:
            self.TypeName = [chr(i) for i in range(65,91)][:self.Ntype]

        for i in range(self.Ntype):
            self.mions += [np.tile(self.mtype[i], self.Nelem[i])]
            self.ChemSymb += [np.tile(self.TypeName[i], self.Nelem[i])]

        self.mions = np.concatenate(self.mions)
        self.ChemSymb = np.concatenate(self.ChemSymb)

    def readxdat(self):
        """ Read VASP XDATCAR """
    
        inp = [line for line in open(self.xdatcar) if line.strip()]
    
    
        ta = inp[5].split()
        tb = inp[6].split()
        # print ta, tb
        if ta[0].isalpha():
            self.TypeName = ta
            self.Ntype = len(ta)
            self.Nelem = np.array(tb, dtype=int)
            self.Nions = self.Nelem.sum()
        else:
            print( "VASP 4.X Format encountered...")
            self.Nelem = np.array(tb, type=int)
            self.Nions = self.Nelem.sum()
            self.Ntype = len(tb)
            self.TypeName = None
            
        if self.isif==3: 
            niter=int(len(inp)/(8+self.Nions)) 
            scale = np.array([float(inp[1+(8+self.Nions)*m]) for m in range(niter)])
            self.cell = np.array([[line.split() for line in inp[2+(8+self.Nions)*m:5+(8+self.Nions)*m]] for m in range(niter)],
                              dtype=float)
            self.cell = np.array([self.cell[m] *scale[m] for m in range(niter)])
        
            pos = np.array([[line.split() for line in inp[8+(8+self.Nions)*m:(8+self.Nions)*(m+1)]]
                        for m in range(niter)],
                        dtype=float)
        else:
            scale = float(inp[1])
            self.cell = np.array([line.split() for line in inp[2:5]],
                              dtype=float)
            self.cell *= scale
            pos = np.array([line.split() for line in inp[7:]
                        if not line.split()[0].isalpha()],
                        dtype=float)            
        self.position = pos.ravel().reshape((-1,self.Nions,3))
        self.Niter = self.position.shape[0]
    
        dpos = np.diff(self.position, axis=0)
        self.positionC = np.zeros_like(self.position)
        # apply periodic boundary condition
        dpos[dpos > 0.5] -= 1.0
        dpos[dpos <-0.5] += 1.0
        # Velocity in Angstrom per femtosecond
        if self.isif==3:
            for i in range(self.Niter-1):
                self.positionC[i,:,:] = np.dot(self.position[i,:,:], self.cell[i])
                dpos[i,:,:] = np.dot(dpos[i,:,:], self.cell[i]) / self.potim
    
            self.positionC[-1,:,:] = np.dot(self.position[-1,:,:], self.cell[-1])
        else:
            for i in range(self.Niter-1):
                self.positionC[i,:,:] = np.dot(self.position[i,:,:], self.cell)
                dpos[i,:,:] = np.dot(dpos[i,:,:], self.cell) / self.potim
    
            self.positionC[-1,:,:] = np.dot(self.position[-1,:,:], self.cell)
    
        self.velocity = dpos


    def readoutcar(self):
        """ read POTIM and POMASS from OUTCAR """

        if os.path.isfile("OUTCAR"):
            # print "OUTCAR found!"
            # print "Reading POTIM & POMASS from OUTCAR..."
            
            outcar = [line.strip() for line in open('OUTCAR')]
            lm = 0;
            for ll, line in enumerate(outcar):
                if 'POTIM  =' in line:
                    # lp = ll
                    self.potim = float(line.split()[2])

                # For ISIF >= 3, VASP output CELLs for each step
                if 'ISIF   =' in line:
                    self.isif = int(line.split()[2])

                if 'Mass of Ions in am' in line:
                    lm = ll + 1

                if lm:
                    break

            # In case Masses not written in OUTCAR
            if lm == 0:
                raise ValueError("Masses for atoms NOT found! Check OUTCAR to see if 'POMASS' for atoms are present!")

            # For heavy atoms, digits for atomic masses may stick together,
            # resulting in cases like: "POMASS =  95.94 32.07183.85"

            pomass_line = outcar[lm]
            pomass_tmp  = pomass_line.split()[2:]
            # Count the number of decimal points, which should equal to the
            # number of types of elements
            if len(pomass_tmp) != pomass_line.count('.'):
                # Fortunately, VASP use fixed-format for printing the atomic
                # masses, i.e. the number of decimal digits for all the floats
                # are the same. Check the last float for this number.

                nd   = pomass_line[::-1].index('.')
                # Find the positions for the decimal points, and add "nd"
                dpos = [ii+nd for ii,xx in enumerate(pomass_line) if xx == '.']
                # Add extra space to the end the number and rejoin the string
                pomass_new_line = ''.join(
                        [xx + ' ' if ii in dpos
                         else xx
                         for ii, xx in
                         enumerate(pomass_line)]
                    )
                self.mtype = np.array(pomass_new_line.split()[2:], dtype=float)
            else:
                self.mtype = np.array(pomass_tmp, dtype=float)

    def getTemp(self, Nfree=None):
        """ Temp vs Time """

        for i in range(self.Niter-1):
            ke = np.sum(np.sum(self.velocity[i,:,:]**2, axis=1) * self.mions / 2.)
            self.Ken[i] = ke * 1E7 / Navogadro / ev
            if Nfree is None:
                Nfree = 3 * (self.Nions - 1)
            self.Temp[i] = 2 * self.Ken[i] / (kB * Nfree)
    def get_center_of_mass(self):
        """
        Compute the center-of-mass position and velocity for each frame.
    
        Returns
        -------
        com_pos : ndarray, shape (Nt, 3)
            Center-of-mass positions in Å
        com_vel : ndarray, shape (Nt, 3)
            Center-of-mass velocities in Å/fs
        """
        masses = self.mions[:, None]  # shape (Nions, 1)
        com_pos = np.sum(self.positionC * masses[None, :, :], axis=1) / np.sum(masses)
        com_vel = np.sum(self.velocity * masses[None, :, :], axis=1) / np.sum(masses)
        return com_pos, com_vel
    
    
    def remove_com_drift(self):
        """
        Remove the center-of-mass velocity from atomic velocities.
        After this, velocity[:,i,:] -= v_com for all atoms i at each frame.
        """
        _, com_vel = self.get_center_of_mass()
        self.velocity -= com_vel[:, None, :]


    import numpy as np


    def find_correlation_cutoff(self, corr,  threshold_ratio=0.02,  min_index=5,  max_fit_points=10,  min_prominence_ratio=0.05 ):
        """
        Determine truncation index for oscillatory correlation functions (JAF / VAF)
        using exponential decay of the peak envelope.
    
        Parameters
        ----------
        corr : array_like (1D)
            Correlation function (JAF or VAF)
        threshold_ratio : float
            Cutoff when envelope decays to this fraction of initial value
        min_index : int
            Minimum index allowed for cutoff
        max_fit_points : int
            Maximum number of early peaks used for exponential fitting
        min_prominence_ratio : float
            Minimum peak prominence relative to corr[0]
    
        Returns
        -------
        cutoff : int
            Index at which correlation function should be truncated
        """
    
        corr = np.asarray(corr, dtype=float)
    
        if corr.ndim != 1:
            raise ValueError("corr must be a 1D array")
    
        n = len(corr)
        if n < min_index + 5:
            raise ValueError("corr is too short to determine cutoff")
    
        c0 = np.abs(corr[0])
        if c0 == 0:
            raise ValueError("corr[0] is zero; cannot normalize envelope")
    
        # -------- step 1: find positive peaks (envelope) --------
        peaks, props = find_peaks(
            np.abs(corr),
            prominence=min_prominence_ratio * c0
        )
    
        if len(peaks) < 3:
            # fallback: simple threshold crossing
            for i in range(min_index, n):
                if np.abs(corr[i]) < threshold_ratio * c0:
                    return i
            return n
    
        # use only early peaks
        peaks = peaks[:max_fit_points]
        peak_vals = np.abs(corr[peaks])
    
        # -------- step 2: exponential envelope fit --------
        def exp_decay(t, A, tau):
            return A * np.exp(-t / tau)
    
        t = peaks.astype(float)
    
        try:
            popt, _ = curve_fit(
                exp_decay,
                t,
                peak_vals,
                p0=(peak_vals[0], t[-1]),
                maxfev=5000
            )
        except RuntimeError:
            # fitting failed → conservative cutoff
            return peaks[-1]
    
        A, tau = popt
        if tau <= 0:
            return peaks[-1]
    
        # -------- step 3: compute cutoff index --------
        t_cut = -tau * np.log(threshold_ratio)
    
        cutoff = int(np.ceil(t_cut))
        cutoff = max(cutoff, min_index)
        cutoff = min(cutoff, n)
    
        return np.min([cutoff,200])






    def getVAF(self):
        """ Velocity Autocorrelation Function """

        # VAF definitions
        # VAF(t) = Natoms^-1 * \sum_i <V_i(0) V_i(t)>
        ############################################################
        # Fast Fourier Transform Method to calculate VAF
        ############################################################
        # The cross-correlation theorem for the two-sided correlation:
        # corr(a,b) = ifft(fft(a)*fft(b).conj()

        # If a == b, then this reduces to the special case of the
        # Wiener-Khinchin theorem (autocorrelation of a):

        # corr(a,a) = ifft(abs(fft(a))**2)
        # where the power spectrum of a is simply:
        # fft(corr(a,a)) == abs(fft(a))**2
        ############################################################
        # in this function, numpy.correlate is used to calculate the VAF

        """
        Velocity Autocorrelation Function (VAF)
    
        Returns
        -------
        time : ndarray
            time array (fs)
        VAF : ndarray
            velocity autocorrelation function
            units: (Å/fs)^2
        """
    
        Nt = self.Niter - 1
        VAF = np.zeros(Nt)
    
        for i in range(self.Nions):
            for a in range(3):
                v = self.velocity[:, i, a]
                corr = scipy.signal.correlate(v, v, mode='full', method='fft')
                corr = corr[Nt-1:]            # t >= 0
                corr /= np.arange(Nt, 0, -1) # (N - t)
                VAF += corr
        
        
        VAF /= self.Nions * 3.0
    
        time = np.arange(Nt) * self.potim  # fs
    
        self.VAF = VAF
        self.Time = time
    
        
        VAF_type=[]
        for Ntype in range(self.Ntype):
            idx = np.where(self.ChemSymb == self.TypeName[Ntype])[0]
            VAF = np.zeros(Nt)
            for i in idx:
                for a in range(3):
                    v = self.velocity[:, i, a]
                    corr = scipy.signal.correlate(v, v, mode='full', method='fft')
                    corr = corr[Nt-1:]            # t >= 0
                    corr /= np.arange(Nt, 0, -1) # (N - t)
                    VAF += corr
            VAF/=self.Nelem[Ntype]*3
            VAF_type.append(VAF)
        self.VAF_type=VAF_type
        
        
    def get_diffusion_coefficient(self,cut=500):
        Ntype=self.Ntype
        d_types=[]
        for i in range(Ntype):
            cut=self.find_correlation_cutoff(self.VAF_type[i])
            tmp=np.trapz(self.VAF_type[i][:cut],self.Time[:cut])*0.1  # A^2/fs to cm^2/s
            d_types.append(tmp)
        self.d_types=d_types        
        
        
    
    def get_J(self,q, remove_drift=True):
        """
        velocity : (Nt, Natoms, 3) in Å/fs
        q        : (Natoms,) charges in e
        return   : J(t) (Nt, 3) in e·Å/fs
        """
        # 核心：原子求和
        velocity=self.velocity
        q = np.asarray(q, dtype=float)

        J = np.sum(velocity * q[None, :, None], axis=1)
        
        if remove_drift:
            J -= J.mean(axis=0)


        self.J = J
        return J
    
    
    # def get_J_trans(self, q, remove_com_drift=True):
    #     """
    #     Transport current excluding polarization contribution (numerically consistent).
    
    #     Parameters
    #     ----------
    #     q : (Natoms,) array
    #         Charges in units of e
    #     remove_com_drift : bool
    #         Remove center-of-mass velocity drift
    
    #     Returns
    #     -------
    #     J_trans : (Nt, 3) ndarray
    #         Transport current in e·Å/fs
    #     """
    
    #     q = np.asarray(q, dtype=float)
    #     Nt, Nat, _ = self.velocity.shape
    #     dt = self.potim  # fs
    
    #     # ---------- velocities (Å/fs) ----------
    #     vel = self.velocity.copy()
    
    #     # ---------- remove COM drift ----------
    #     if remove_com_drift:
    #         v_com = vel.mean(axis=1, keepdims=True)  # (Nt,1,3)
    #         vel -= v_com
    
    #     # ---------- total current ----------
    #     # J_total(t) = Σ_i q_i v_i(t)
    #     J_total = np.sum(vel * q[None, :, None], axis=1)  # (Nt,3)
    
    #     # ---------- reconstruct positions from corrected velocity ----------
    #     pos = np.zeros((Nt, Nat, 3), dtype=float)
    #     pos[0] = self.positionC[0]
    
    #     for t in range(1, Nt):
    #         pos[t] = pos[t-1] + vel[t-1] * dt
    
    #     # ---------- polarization ----------
    #     # P(t) = Σ_i q_i r_i(t)
    #     P = np.sum(pos * q[None, :, None], axis=1)  # (Nt,3)
    
    #     # ---------- polarization current ----------
    #     J_pol = np.zeros_like(P)
    #     J_pol[1:-1] = (P[2:] - P[:-2]) / (2 * dt)
    #     J_pol[0]    = (P[1]  - P[0])   / dt
    #     J_pol[-1]   = (P[-1] - P[-2])  / dt
    
    #     # ---------- transport current ----------
    #     J_trans = J_total - J_pol
    
    #     # ---------- remove residual DC ----------
    #     J_trans -= J_trans.mean(axis=0)
    
    #     self.J_trans = J_trans
    #     return J_trans

    
    
    
    def get_J_trans(self, q, remove_com_drift=True, smooth=False, sigma=5):
        """
        Transport current excluding polarization contribution.
    
        Parameters
        ----------
        q : (Natoms,) array
            charges in units of e
        remove_com_drift : bool
            remove center-of-mass velocity drift before computing current
        smooth : bool
            whether to smooth J_trans(t)
        sigma : float
            Gaussian smoothing width (fs index)
    
        Returns
        -------
        J_trans : (Nt, 3) ndarray
            transport current in e·Å/fs
        """
    
        q = np.asarray(q, dtype=float)
        Nt, Nat, _ = self.velocity.shape
        dt = self.potim  # fs
    
        # ---------- velocities (Å/fs) ----------
        vel = self.velocity.copy()
    
        # ---------- remove COM drift ----------
        # if remove_com_drift:
        #     v_com = vel.mean(axis=1, keepdims=True)   # (Nt,1,3)
        #     vel -= v_com
    
        # ---------- total current ----------
        # J_total(t) = Σ_i q_i v_i(t)
        J_total = np.sum(vel * q[None, :, None], axis=1)  # (Nt,3)
    
        # ---------- polarization ----------
        # P(t) = Σ_i q_i r_i(t)   (Å·e)
        pos = self.positionC[:Nt]   # Cartesian positions (Å)
        P = np.sum(pos * q[None, :, None], axis=1)  # (Nt,3)
    
        # ---------- polarization current ----------
        # central finite difference for stability
        J_pol = np.zeros_like(P)
        J_pol[1:-1] = (P[2:] - P[:-2]) / (2 * dt)
        J_pol[0]    = (P[1]  - P[0])   / dt
        J_pol[-1]   = (P[-1] - P[-2])  / dt
    
        # ---------- transport current ----------
        J_trans = J_total - J_pol
    
        # ---------- optional smoothing ----------
        # if smooth:
        #     from scipy.ndimage import gaussian_filter1d
        #     J_trans = gaussian_filter1d(J_trans, sigma=sigma, axis=0, mode='nearest')
    
        self.J_trans = J_trans
        return J_trans

    def get_JAF(self):
        """
        J : (Nt, 3) electric current
        return : JAF(t) (Nt,)  in (e·Å/fs)^2
        """
        
        J=self.J
        Nt = J.shape[0]
        JAF = np.zeros(Nt)
    
        for i in range(3):
            c = scipy.signal.correlate(J[:, i], J[:, i], mode='full', method='fft')
            JAF += c[Nt-1:]  # 取正时间
    
        # 归一化（与 getVAF 一致的统计意义）
        JAF /= np.arange(Nt, 0, -1)

        self.JAF = JAF
        
        JAF_smooth = copy.deepcopy(JAF)        
        JAF_smooth [200:] = scipy.ndimage.gaussian_filter1d(JAF[200:],sigma=20,mode='nearest')
        self.JAF_smooth=JAF_smooth
        
        return JAF
    
    def get_sigma_gk(self, T, volume_A3, acf_fraction=0.9):
        """
        JAF       : (Nt,) in (e·Å/fs)^2
        dt_fs     : timestep in fs
        T         : K
        V_A3      : volume in Å^3
        acf_fraction : fallback cutoff if ACF never crosses zero
        """
        dt_fs = self.potim
        JAF = self.JAF
        # JAF = self.JAF_smooth
        
        # ---------- constants ----------
        kB = 1.380649e-23      # J/K
        e  = 1.602176634e-19   # C
        A  = 1e-10             # m
        fs = 1e-15             # s
    
        # ---------- unit conversion ----------
        JAF_SI = JAF * (e*A/fs)**2
        dt_s   = dt_fs * fs
        V_m3   = volume_A3 * A**3
        
        # ---------- automatic truncation ----------
        # JAF_norm = JAF / JAF[0]
    
        # idx = np.where(JAF_norm < 0)[0]
        # if len(idx) > 0:
        #     cut = idx[0]
        # else:
        #     cut = int(len(JAF) * acf_fraction)
        # cut = int(len(JAF) * acf_fraction)
        # ---------- cumulative integral ----------
        cut=self.find_correlation_cutoff(self.JAF)
        sigma_t = np.zeros(cut)
        integ = 0.0
        for i in range(cut):
            integ += JAF_SI[i] * dt_s
            sigma_t[i] = integ / (3 * kB * T * V_m3)
        
        self.sigma_t=sigma_t
        return sigma_t

    
    
    def phononDos(self, unit='THz', sigma=5):
        """ Phonon DOS from VAF """

        N = self.Niter - 1
        # Frequency in THz
        omega = fftfreq(2*N-1, self.potim) * 1E3
        # Frequency in cm^-1
        if unit.lower() == 'cm-1':
            omega *= 33.35640951981521
        if unit.lower() == 'mev':
            omega *= 4.13567
        # from scipy.ndimage.filters import  gaussian_filter1d as gaussian
        # smVAF = gaussian(self.VAF2, sigma=sigma)
        # pdos = np.abs(fft(smVAF))**2
        if self.VAF2 is None:
            self.getVAF()
        pdos = np.abs(fft(self.VAF2 - np.average(self.VAF2)))**2

        return omega[:N], pdos[:N]

    def PCF(self, bins=50, Niter=10, A='', B=''):
        """ Pair Correlation Function """

        if not A:
            A = self.TypeName[0]
        if not B:
            B = A

        whichA = self.ChemSymb == A
        whichB = self.ChemSymb == B
        indexA = np.arange(self.Nions)[whichA]
        indexB = np.arange(self.Nions)[whichB]
        posA = self.position[:,whichA,:]
        posB = self.position[:,whichB,:]

        steps = range(0, self.Niter, Niter)
        rABs = np.array([posA[i,k,:]-posB[i,j,:]
                         for k in range(indexA.size)
                         for j in range(indexB.size)
                         for i in steps
                         if indexA[k] != indexB[j]])
        # periodic boundary condition
        rABs[rABs > 0.5] -= 1.0
        rABs[rABs <-0.5] += 1.0
        # from direct to cartesian coordinate
        rABs = np.linalg.norm(np.dot(self.cell, rABs.T), axis=0)
        # histogram of pair distances
        val, b = np.histogram(rABs, bins=bins)
        # density of the system
        rho = self.Nions / np.linalg.det(self.cell)
        # Number of A type atom
        Na = self.Nelem[self.TypeName.index(A)]
        # Number of B type atom
        Nb = self.Nelem[self.TypeName.index(B)]
        dr = b[1] - b[0]
        val = val * self.Nions / (4*np.pi*b[1:]**2 * dr) / (Na * Nb * rho) / len(steps)

        return val, b[1:]

################################################################################      

# test code of the above class
if __name__ == '__main__':
    
    t_start=time.time()
    
    Temperature = float(sys.argv[1])
    
    inp = xdatcar()
    inp.remove_com_drift()

    q = np.zeros(inp.Nions)
    q[inp.ChemSymb == 'H'] = 1
    q[inp.ChemSymb == 'C'] = 4
    q[inp.ChemSymb == 'O'] = -2


    inp.get_J(q)
    # JAF_smooth = inp.JAF_smooth
    # J_trans = inp.get_J_trans(q, remove_com_drift=True)
    # inp.J = J_trans
    JAF = inp.get_JAF()
    
   
    t = inp.Time
    plt.figure()    
    plt.plot(t, JAF)
    # plt.plot(t, JAF_smooth)
    jaf_cut=inp.find_correlation_cutoff(inp.JAF)
    plt.xlabel('Time (fs)')
    plt.ylabel(r'JAF (e$^2$·$\rm \AA^2$/fs$^2$)')
    plt.vlines(jaf_cut,ymin=np.min(inp.JAF),ymax=np.max(inp.JAF),color='r')
    # plt.show()
    plt.xlim(-100,3000)
    plt.savefig('JAF.png',dpi=600,bbox_inches='tight',transparent=True)


    header="{:^20}{:^20}".format('  time(fs)',  'JAF(e^2·A^2/fs^2)')
    data=np.column_stack([t,JAF])
    np.savetxt('JAF.dat',data,fmt='%20.8e',header=header)



    V = np.linalg.det(inp.cell)  #Angstrom^3
    # conductivity
    sigma_t = inp.get_sigma_gk(T = Temperature, volume_A3=V, acf_fraction=0.6)
    t = inp.Time[:len(sigma_t)]
    plt.figure() 
    plt.plot(t, sigma_t * 10)  # S/m → mS/cm
    plt.xlabel('Time (fs)')
    plt.ylabel(r'$\sigma$ (mS/cm)')
    # plt.show()
    # plt.xlim(-100,3000)
    plt.savefig('sigma.png',dpi=600,bbox_inches='tight',transparent=True)
    
    header="{:^20}{:^20}".format('  time(fs)',  'sigma(mS/cm)')
    data=np.column_stack([t,sigma_t*10])
    np.savetxt('sigma_t.dat',data,fmt='%20.8e',header=header)
    
    inp.getVAF()
    inp.get_diffusion_coefficient()
    plt.figure()
    vaf_cut1=inp.find_correlation_cutoff(inp.VAF_type[0])
    vaf_cut2=inp.find_correlation_cutoff(inp.VAF_type[1])
    vaf_cut3=inp.find_correlation_cutoff(inp.VAF_type[2])
    
    plt.vlines(vaf_cut1,ymin=np.min(inp.VAF_type[0]),ymax=np.max(inp.VAF_type[0]),color='r')
    plt.vlines(vaf_cut2,ymin=np.min(inp.VAF_type[1]),ymax=np.max(inp.VAF_type[1]),color='g')
    plt.vlines(vaf_cut3,ymin=np.min(inp.VAF_type[2]),ymax=np.max(inp.VAF_type[2]),color='b')
    
    plt.plot(inp.Time,inp.VAF_type[0],label='H')
    plt.plot(inp.Time,inp.VAF_type[1],label='C')
    plt.plot(inp.Time,inp.VAF_type[2],label='O')
    plt.legend()
    
    plt.xlim(-100,3000)
    plt.savefig('vacf.png',dpi=600,bbox_inches='tight',transparent=True)
    
    
    header="{:^20}{:^20}{:^20}{:^20}{:^20}".format('  time(fs)',  'vacf[H]','vacf[C]','vacf[O]','unit: A^2/fs^2')
    data=np.column_stack([inp.Time,inp.VAF_type[0],inp.VAF_type[1],inp.VAF_type[2]])
    np.savetxt('vacf.dat',data,fmt='%20.8e',header=header)
    
    
    print('d(H):', inp.d_types[0],'cm2/s')
    print('d(C):', inp.d_types[1],'cm2/s')
    print('d(O):', inp.d_types[2],'cm2/s')
    
    header="{:^20}{:^20}{:^20}{:^20}".format('sigma (mS/cm)',  'd_H (cm^2/s)','d_C (cm^2/s)','d_O (cm^2/s)')
    data=np.array([sigma_t[-1]*10,inp.d_types[0],inp.d_types[1],inp.d_types[2]]).reshape(1,4)
    np.savetxt('diffusion_conductivity.dat',data,fmt='%20.8e',header=header)

    t_end=time.time()
    print('total time: ',t_end-t_start,'s')

    # velocity=inp.velocity
    
    # vel_H=[]
    # vel_C=[]
    # vel_O=[]
    # for i in range(len(velocity[:,0,0])):
    #     tmp_H=[]
    #     tmp_C=[]
    #     tmp_O=[]
    #     for j in range(0,inp.Nelem[0]):
    #         v=np.linalg.norm(velocity[i,j,:])
    #         tmp_H.append(v)
    #     for j in range(inp.Nelem[0],inp.Nelem[0]+inp.Nelem[1]):            
    #         v=np.linalg.norm(velocity[i,j,:])
    #         tmp_C.append(v)
    #     for j in range(inp.Nelem[0]+inp.Nelem[1],inp.Nelem[0]+inp.Nelem[1]+inp.Nelem[2]):            
    #         v=np.linalg.norm(velocity[i,j,:])
    #         tmp_O.append(v)
    #     v_H_avg=np.average(np.array(tmp_H))
    #     v_C_avg=np.average(np.array(tmp_C))
    #     v_O_avg=np.average(np.array(tmp_O))
    #     vel_H.append(v_H_avg)
    #     vel_C.append(v_C_avg)
    #     vel_O.append(v_O_avg)
    # plt.figure()
    # plt.hist(vel_H)
    # plt.hist(vel_C)
    # plt.hist(vel_O)   
    
    
    # plt.figure()
    
    # vel_H=[]
    # vel_C=[]
    # vel_O=[]
    # for i in range(len(velocity[:,0,0])):
    #     tmp_Hx=np.sum(velocity[i,0:inp.Nelem[0],0])
    #     tmp_Hy=np.sum(velocity[i,0:inp.Nelem[0],1])
    #     tmp_Hz=np.sum(velocity[i,0:inp.Nelem[0],2])
    #     tmp_Cx=np.sum(velocity[i,inp.Nelem[0]:inp.Nelem[0]+inp.Nelem[1],0])
    #     tmp_Cy=np.sum(velocity[i,inp.Nelem[0]:inp.Nelem[0]+inp.Nelem[1],1])
    #     tmp_Cz=np.sum(velocity[i,inp.Nelem[0]:inp.Nelem[0]+inp.Nelem[1],2])
    #     tmp_Ox=np.sum(velocity[i,inp.Nelem[0]+inp.Nelem[1]:,0])
    #     tmp_Oy=np.sum(velocity[i,inp.Nelem[0]+inp.Nelem[1]:,1])
    #     tmp_Oz=np.sum(velocity[i,inp.Nelem[0]+inp.Nelem[1]:,2])

    #     vel_H.append([tmp_Hx,tmp_Hy,tmp_Hz])
    #     vel_C.append([tmp_Cx,tmp_Cy,tmp_Cz])
    #     vel_O.append([tmp_Ox,tmp_Oy,tmp_Oz])
    
    # vel_H=np.array(vel_H)
    # vel_C=np.array(vel_C)
    # vel_O=np.array(vel_O)
    # plt.hist(vel_H[:,1])
    # plt.hist(vel_C[:,1])
    # plt.hist(vel_O[:,1])  
    

