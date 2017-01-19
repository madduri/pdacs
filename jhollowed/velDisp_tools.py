'''
Joe Hollowed
Last edited 1/16/2017

Module prodiving functions relevant to velocity dispersion analysis of simulated 
dark matter halos and observational galaxy clusters, intended for use with the 
PDACS platform (pdacs.mcs.anl.gov , portal-auth.nersc.gov/galaxy-pdacs-dev)
'''

import dtk
import glob
import numpy as np
import matplotlib.pyplot as plt
import pdb
import csv
from astropy.constants import M_sun
from astropy.cosmology import WMAP7 as cosmo
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as plticker
import datetime
import os

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

def velDisp_mass_rel_z0(sod_file):
	'''
	Function to gather mass and dark matter velocity dispersion data columns from
	sodproperties halo files, scale masses by redshift, and return a .csv of the data
	
	:param sod_file: path to a sodproperties file
	:return: path to .csv file with mass and DM velocity dispersion columns
	'''

	masses = dtk.gio_read(sod_file, 'sod_halo_mass')
	vDisp = dtk.gio_read(sod_file, 'sod_halo_vel_disp')
	h_1 = 1/(cosmo.h)
	mask = (masses*h_1) >= 1e14
	masses_big = masses[mask] * h_1
	vDisp_big = vDisp[mask]

	today = str(datetime.date.today())
	with open('velDisp_vs_SodMass_{}.csv'.format(today), 'wb') as output:
		writer = csv.writer(output, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		for row in zip(masses_big, vDisp_big): 
			writer.writerow(list(row))
			output.flush()
		return os.path.abspath(output.name)
	

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

def velDisp_mass_rel_z0_plot(sig_m_file):
	'''
	Function to plot the mass vs dark matter velocity dispersion of halos as
	a scatter plot, including the Evrard et al. (2007) relation, and a best
	fit to the input data.

	:param sig_m_file: .csv file handle with columns of redshift-scaled mass and 
			   velocity dispersion (each row one halo)
	:return: path to matplotlib figure .svg file
	'''
	masses = []
	vDisp = []
	with open(sig_m_file, 'rb') as sig_m_data:	
		reader = csv.reader(sig_m_data, delimiter=',', quotechar='|')
		for row in reader: 
			masses.append(float(row[0]))
			vDisp.append(float(row[1]))
			
	# Overplotting Evrard relation
	t_x = np.linspace(7e13, 2e15, 300)
	h = cosmo.h
	sig_dm15 = 1082.9
	alpha = 0.3361
	t_y = sig_dm15 * (( (h*t_x) / (1e15) )**alpha)

	# preforming Least Squares on AlphaQuadrant Data
	# (fitting to log linear form, as in Evrard et al.)
	# X = feature matrix (masses)
	# P = parameter matrix (sig_dm15(log intercept) and alpha(log slope))
	X = np.array([np.log(mi / 1e15) for mi in masses])
	X = np.vstack([X, np.ones(len(X))]).T
	P = np.linalg.lstsq(X, np.log(vDisp))[0]
	alpha_fit = P[0]
	sig_dm15_fit = np.e**P[1]

	fit_x = np.linspace(7e13, 2e15, 300)
	fit_y = sig_dm15_fit * (( (fit_x) / (1e15) )**alpha_fit)


	# plotting
	fig = plt.figure()
	plt.rc('text', usetex=True)
	plt.rc('font', **{'family' : "sans-serif"})
	params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
	plt.rcParams.update(params)

	ax = fig.add_subplot(1,1,1)
	p = ax.loglog(masses, vDisp, 'rx', markersize=9, label=r'$\text{halos with } M_{200}h(z) > 1e14$')
	ax.hold(True)

	t = ax.loglog(t_x, t_y, '--b', linewidth = 1.2, label=r'$\text{Evrard et al.(2007) relation}$')
	fit = ax.loglog(fit_x, fit_y, 'b', linewidth = 1.2, label=r'$\text{least squares fit}$')
	
	ax.set_ylim([300, 1400])
	ax.set_xlim([7e13, 2e15])
	ax.yaxis.set_major_formatter(ScalarFormatter())
	loc = plticker.MultipleLocator(base=100)
	ax.yaxis.set_major_locator(loc)
	ax.legend(loc=4)
	ax.set_ylabel(r'$\sigma_{DM}\>\>(\text{km s}^{-1})$', fontsize=16)
	ax.set_xlabel(r'$M_{200}\>\>(10^{15}h^{-1}M_{\odot})$', fontsize=16)
	fig_name = '{}_plot.svg'.format(sig_m_file.split('.')[-2])
	plt.savefig(fig_name, format='svg')
	
	fig_svg = open(fig_name, 'rb')
	fig_path = os.path.abspath(fig_svg.name)
	fig_svg.close()
	return(fig_path)

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
