#!/usr/bin/env python

# Dominique Marion, Diego F. Gauto, Isabel Ayala, Audrey Hessel, Paul Schanda 
# Univ. Grenoble Alpes, CEA CNRS, Institute for Structural Biology (IBS) 71 avenue des martyrs,
# 38044 Grenoble (France); E-mail: dominique.marion@ibs.fr, paul.schanda@ibs.fr
 


#_______________________________________________________________________
#PACKAGES
#
import os
import numpy as np
import scipy.linalg as sci

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,show,clf,close,figure,hold,savefig, errorbar, title,semilogx
from math import sqrt,exp,pi,e
import scipy.io as sio
from pylab import show, fft, rfft, fftshift
from scipy.optimize import fmin
import math
import sys
from random import gauss
from time import gmtime, strftime
from numpy.lib.scimath import logn

#_______________________________________________________________________
#FUNCTIONS PART
##_______________________________________________________________________


# monoexp_fitting fuction
def fit(fix_ori, xaxis):
	fix=np.zeros(len(fix_ori))
	l1 = len(fix_ori)
	l2 = 3*l1/4
	l3 = l1/4
	for i in range(len(fix_ori)):
		fix[i]=fix_ori[i]/fix_ori[0]

#	v0= [1.4,fix[0]]
#	print 'Initial value: used      ' , '[%7.3f %3.1f]' % (v0[0], v0[1])
#	v0= [(logn(e,fix[5])-logn(e,fix[l3]))/(xaxis[l3]-xaxis[5]),fix[0]]
#	print 'Initial value: computed   5-', l3 , '[%7.3f %3.1f]' % (v0[0], v0[1])
	v0= [(logn(e,fix[30])-logn(e,fix[l2]))/(xaxis[l2]-xaxis[30]),fix[0]]
	print 'Initial value: computed  30-',l2 , '[%7.3f %3.1f]' % (v0[0], v0[1])
	
	vv = fmin(expdecay, v0, args=(fix, xaxis),maxiter=1000, maxfun=1000,full_output=True,disp=True)
	v=vv[0]
	chi_mono=vv[1]
	fitparam = [v[0],v[1],vv[1]]
	print '[%7.3f %3.1f]' % (v[0], v[1]) , ' chi_mono %6.2e' % chi_mono
	return fitparam
#


# double exp_fitting fuction
def fit2(fix_ori, xaxis):
	fix=np.zeros(len(fix_ori))
	l1 = len(fix_ori)
	l2 = 3*l1/4
	l3 = l1/4
	for i in range(len(fix_ori)):
		fix[i]=fix_ori[i]/fix_ori[0]
	v0= [1.4,fix[0], 1.4,fix[0]]
#	print 'Initial value: used      ' , '[%7.3f %3.1f]' % (v0[0], v0[1])
#	v0= [(logn(e,fix[5])-logn(e,fix[l3]))/(xaxis[l3]-xaxis[5]),fix[0]]
#	print 'Initial value: computed   5-', l3 , '[%7.3f %3.1f]' % (v0[0], v0[1])
	decay = (logn(e,fix[30])-logn(e,fix[l2]))/(xaxis[l2]-xaxis[30])
	v0= [decay*2.0,fix[0]/2.0,decay/2.0,fix[0]/2.0]
	print 'Initial value: computed  30-',l2 , '[%7.3f %3.1f ][ %7.3f %3.1f]' % (v0[0], v0[1], v0[2], v0[3])
	
	vv = fmin(double_expdecay, v0, args=(fix, xaxis),maxiter=1000, maxfun=1000,full_output=True,disp=True)
	v=vv[0]
	chi_mono=vv[1]
	fitparam = [v[0],v[1],v[2],v[3],vv[1]]
	print '[%7.3f %3.1f ][%7.3f %3.1f]' % (v[0], v[1], v[2], v[3]), ' chi_bi %6.2e' % chi_mono
	return fitparam
#_______________________________________________________________________

#_______________________________________________________________________
# fmin function expdecay

def expdecay(v, fix, xaxis):
	allchi2=[]
	
	for i in range(len(xaxis)):
		chi2=(v[1]*np.exp(-1.*v[0]*xaxis[i])-fix[i])**2
		allchi2.append(chi2)
	allchi2sum=sum(allchi2)/len(xaxis)
	return allchi2sum
#_________________________________________________________________
# fmin function double_expdecay

def double_expdecay(v, fix, xaxis):
	allchi2=[]
	
	for i in range(len(xaxis)):
		chi2=(v[1]*np.exp(-1.*v[0]*xaxis[i])  + v[3]*np.exp(-1.*v[2]*xaxis[i])-fix[i])**2
		allchi2.append(chi2)
	allchi2sum=sum(allchi2)/len(xaxis)
	return allchi2sum
#_________________________________________________________________


#___________________________________________________________________
#MAIN PART
#
#_______________________________________________________________________
#
# simulation of the calculation of Rr1rho with spinlock sequences
#_______________________________________________________________________
# initialisation

if len(sys.argv) < 2:
    print '-->  Missing file name'
    print 'Syntax ', sys.argv[0], ' filename'
    sys.exit(0)
filename = sys.argv[1]
comment = sys.argv[2]
#print filename

start=strftime('%d.%m - %H:%M:%S', gmtime())
print '\n Script start at:', start


######################### often changed parameters ################################
output_folder=filename + '.fit'





w1_gamma=[12000]

poplist=['2spin_50perc.sys']
#pop=['3perc','10perc']
pop=['50perc']


############## less changed parameters ############################
fieldlist=[60]
dppmlist=[0.0]
D=20000 

simname=[]
simfile=[]
SxIa_gamma_real=[]
SxIb_gamma_real=[]
temps_gamma=[]
temps_gamma1=[]


########################create folders##################################
if not os.path.exists(output_folder): os.makedirs(output_folder)
if not os.path.exists(output_folder+'/matrices'): os.makedirs(output_folder+'/matrices')
if not os.path.exists(output_folder+'/graphs_fit_2spins'): os.makedirs(output_folder+'/graphs_fit_2spins')


simname.append(filename)

for i in range(len(simname)):
  simfile.append('./'+str(simname[i])+'.mat')

#Recovery of the datas from the matlab files

for i in range(len(simfile)):                                          
  SxIa_gamma_real.append(np.real(sio.loadmat(simfile[i])[str(simname[i])][0]))
  temps_gamma1.append(np.real(sio.loadmat(simfile[i])[str(simname[i])][1]))

temps_gamma=temps_gamma1[0]

#skip_points=int(5.*len(Mx_gamma_real[0][0])/200.)

#temps_gamma=np.arange(0.,0.0001*len(SxIa_gamma_real[0]),0.0001)
print 'Nb points  ', len(temps_gamma)
print '#0 %7.2f #1 %7.2f .. #%3d .. %7.2f (us)' %(temps_gamma[0]*1e6,temps_gamma[1]*1e6,len(temps_gamma)-1, temps_gamma[-1]*1e6)

#print 'Data length',len(SxIa_gamma_real[0])
number_sims=len(SxIa_gamma_real)

###########################################################################################################################
# Calculation of r1rho
###########################################################################################################################

#AllRr1rho_calc_gamma=np.zeros((number_sims,len(Mx_gamma_real[0])-1))   #initialisation
rate_SxIa=np.zeros(number_sims)   #initialisation
monoexp_chi_SxIa=np.zeros(number_sims)

rate_SxIa_1=np.zeros(number_sims)  
rate_SxIa_2=np.zeros(number_sims)  
ampl_SxIa_1=np.zeros(number_sims)  
ampl_SxIa_2=np.zeros(number_sims)  
biexp_chi_SxIa=np.zeros(number_sims)


biexp_chi_SxIa=np.zeros(number_sims)

for i in range(number_sims):                       #loop for the number of simulations #for ii in range(len(SxIa_gamma_real[0])-1):                  #loop for the aimantation of each f1

	fitparam2_SxIa=fit2(SxIa_gamma_real[i][0:],temps_gamma[0:])
	rate_SxIa_1[i]=fitparam2_SxIa[0]
	rate_SxIa_2[i]=fitparam2_SxIa[2]
	ampl_SxIa_1[i]=fitparam2_SxIa[1]/ (fitparam2_SxIa[1] + fitparam2_SxIa[3])
	ampl_SxIa_2[i]=fitparam2_SxIa[3]/ (fitparam2_SxIa[1] + fitparam2_SxIa[3])
	biexp_chi_SxIa[i]=fitparam2_SxIa[4]
    
	fitparam_SxIa=fit(SxIa_gamma_real[i][0:],temps_gamma[0:])
	rate_SxIa[i]=fitparam_SxIa[0]
	monoexp_chi_SxIa[i]=fitparam_SxIa[2]


###########################################################################################################################
# Create figures
################################## monoexp ##############################################################################

	Evo=plt.figure()
	plot(temps_gamma,SxIa_gamma_real[i]/SxIa_gamma_real[i][0],color='r')
	plt.hold(True)

	plot(temps_gamma,fitparam_SxIa[1]*np.exp(-fitparam_SxIa[0]*np.array(temps_gamma)),color='b')
	plot(temps_gamma,fitparam2_SxIa[1]*np.exp(-fitparam2_SxIa[0]*np.array(temps_gamma))+		\
	         fitparam2_SxIa[3]*np.exp(-fitparam2_SxIa[2]*np.array(temps_gamma)),color='g')
 

	rate1= '(1) %7.2f  [1.0]\n        (2) %7.2f %7.2f  [%3.2f :%3.2f]' % (rate_SxIa[i] , rate_SxIa_1[i], rate_SxIa_2[i], \
			ampl_SxIa_1[i], ampl_SxIa_2[i]) 
	
	str_rate1 = "{:.6e}    {:.6e} {:.6e} {:.2f} {:.2f}  ".format(rate_SxIa[i], rate_SxIa_1[i], rate_SxIa_2[i],  \
			ampl_SxIa_1[i], ampl_SxIa_2[i])
	plt.figtext(0.25, 0.80, 'Rate '+rate1 + '    ' + comment, size=10)
	 
#	plt.ylim([-0.1, 1.1])
	plt.title('M vs time')
#	plt.axis([0,0.16,0,1.0])
	plt.savefig(str(output_folder)+'/graphs_fit_2spins/'+str(simname[i])+'monoexponential.pdf',format='pdf', dpi=300)



###########################################################################################################################
# Save data in MATLAB format
################################### monoexp ##################################################################################

sio.savemat(str(output_folder)+"/matrices/monoexp_rate_SxIa",{'monoexp_rate_SxIa':rate_SxIa})
sio.savemat(str(output_folder)+"/matrices/monoexp_chi_SxIa",{'monoexp_chi_SxIa':monoexp_chi_SxIa})



fname = "Summary"
with open(fname, "a+") as f:

  print >> f, str_rate1, filename
f.close()


close()

print '\n Script end at:  ', strftime('%d.%m - %H:%M:%S', gmtime()), '\n'
