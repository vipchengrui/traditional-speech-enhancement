import numpy as np
import scipy.interpolate as itp
import math

# Weighted spectral average
class Init_Weight(object):

	def __init__(self,ns_ps,fs):
		self.ns_ps = ns_ps
		self.fs = fs
	
	def info(self):
		parameters = {'ass':0.85 , 'beta':1.5 , 'noise_ps':self.ns_ps , 'P':self.ns_ps}
		return parameters

# Continuous minimal tracking
class Init_ConMinTrack(object):
	
	def __init__(self,ns_ps,fs):
		self.ns_ps = ns_ps
		self.fs = fs
	
	def info(self):
		len_val = len(self.ns_ps)
		parameters = {'n':2,'leng':len_val,'alpha':0.7,'beta':0.96,'gamma':0.998,\
		'noise_ps':self.ns_ps,'pxk_old':self.ns_ps,'pxk':self.ns_ps,'pnk_old':self.ns_ps,'pnk':self.ns_ps}
		return parameters

# MCRA algorithm
class Init_MCRA(object):

	def __init__(self,ns_ps,fs):
		self.ns_ps = ns_ps
		self.fs = fs

	def info(self):
		len_val = len(self.ns_ps)
		parameters = {'n':2,'ad':0.95,'ass':0.8,'L':1000*2//20,'delta':5,'ap':0.2,'leng':len_val,\
		'P':self.ns_ps,'Pmin':self.ns_ps,'Ptmp':self.ns_ps,'pk':np.zeros(len_val),'noise_ps':self.ns_ps}
		return parameters

# MCRA2 algorithm
class Init_MCRA2(object):

	def __init__(self,ns_ps,fs):
		self.ns_ps = ns_ps
		self.fs = fs

	def info(self):
		len_val = len(self.ns_ps)
		freq_res = self.fs / len_val
		k_1khz = int(1000 // freq_res)
		k_3khz = int(3000 // freq_res)
		
		# [9.60] delta
		low_1 = 2*np.ones(k_1khz,dtype=np.int),
		low_2 = 2*np.ones(k_3khz-k_1khz,dtype=np.int),
		high = 5*np.ones(len_val//2-k_3khz,dtype=np.int),
		delta_val = np.append(np.append(np.append(np.append(np.append(low_1,low_2),high),high),low_2),low_1)

		parameters = {'n':2,'leng':len_val,'ad':0.95,'ass':0.8,'ap':0.2,'beta':0.8,'beta1':0.98,'gamma':0.998,'alpha':0.7,\
		'delta':delta_val,'pk':np.zeros(len_val),'noise_ps':self.ns_ps,'pxk_old':self.ns_ps,'pxk':self.ns_ps,'pnk_old':self.ns_ps,'pnk':self.ns_ps}

		return parameters

'''
# Minimum Statistics Algorithm
class Init_MS(object):

	def __init__(self,ns_ps):
		self.ns_ps = ns_ps

	def info(self):
		len_val = len(self.ns_ps)
		L_val = len_val
		R_val = len_val / 2
		D_val = 150
		V_val = 15
		Um_val = 10
		Av_val = 2.12
		alpha_max_val = 0.96
		alpha_min_val = 0.3
		beta_max_val = 0.8
		x_val = [1 , 2 , 5 , 8 , 10 , 15 , 20 , 30 , 40 , 60 , 80 , 120 , 140 , 160]
		Y_M_val = [0 , 0.26 , 0.48 , 0.58 , 0.61 , 0.668 , 0.705 , 0.762 , 0.8 , 0.841 , 0.865 , 0.89 , 0.9 , 0.91]
		Y_H_val = [0 , 0.15 , 0.48 , 0.78 , 0.98 , 1.55 , 2.0 , 2.3 , 2.52 , 2.9 , 3.25 , 4.0 , 4.1 , 4.1]
		xi_val = D_val
		M_D_val = itp.spline(x_val , Y_M_val , xi_val)		#interpolate
		H_D_val = itp.spline(x_val , Y_H_val , xi_val)
		xi_val = V_val
		M_V_val = itp.spline(x_val , Y_M_val , xi_val)
		H_V_val = itp.spline(x_val , Y_H_val , xi_val)
		
		#minact_val[1 : L_val , 1 : Um_val]=np.maximum(self.ns_ps)
		
		parameters = {'n' : 2 , 'len' : len_val , 'alpha_corr' : 0.96 , 'alpha' : 0.96*np.ones(len_val) , 'P' : self.ns_ps , \
		'noise_ps' : self.ns_ps , 'Pbar' : self.ns_ps , 'Psqbar' : self.ns_ps , 'actmin' : self.ns_ps , 'actmin_sub' : self.ns_ps ,\
		'Pmin_u' : self.ns_ps , 'subwc' : 2 , 'u' : 1 , 'lmin_flag' : np.zeros(len_val) , 'L' : L_val ,\
		'R' : R_val , 'D' : D_val ,'V' : V_val , 'Um' : Um_val , 'Av' : Av_val , 'alpha_max' : alpha_max_val ,\
		'alpha_min' : alpha_min_val , 'beta_max' : beta_max_val , 'Y_M' : Y_M_val , 'Y_H' : Y_H_val , 'M_D' : M_D_val ,\
		'H_D' : H_D_val , 'M_V' : M_V_val , 'H_V' : H_V_val}

		return parameters
'''