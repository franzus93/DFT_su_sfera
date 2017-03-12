#!/usr/bin/python2.7                           
                                               
import matplotlib.pyplot as plt			#
import numpy as np				#
import sys					#
						#
						#                                               
import pyshtools as sht                        	#
from pyshtools.expand import SHExpandDH;       	#
from pyshtools.expand import MakeGridDH;       	#
from scipy.special import sph_harm             	#
					       	#
import dft_su_sfera as dfts		       	#
					       	#
					       	#
#######################################################################
####			FUNCTIONS				#######
#######################################################################
		#
		#
########################################
####		INPUT		########
########################################

def get_pot(filename,lmax):
	"reads file formatted as: #legendre coefficient# ... is the legendre trasform of the potential"
	INPUT = open(filename,'r')
	
	W = np.zeros(lmax+1,order='F')
	
	l = 0
	
	for line in INPUT:
		W[l] = float(line)
		l += 1
		if l>lmax:
			break
		
	INPUT.close()
	
	return W
	
def get_distro(filename, lmax):
	"reads file formatted as: # cos/sin ## l ## m ## coeff # ... coefficients for the starting distribution"
	INPUT = open(filename,'r')
	
	Clm = np.zeros((2,lmax+1,lmax+1), order='F')
	
	for line in INPUT:
		fields = line.split()
		u = int(fields[0])
		l = int(fields[1])
		m = int(fields[2])
		C = float(fields[3])
		if l < lmax+1:
			Clm[u,l,m] = C
	
	INPUT.close()
	
	return Clm

def get_parm(filename):
	"reads parameter file and changes the values"
	INPUT = open(filename,'r')
	
	#define global parameters
	
	global BETA
	global L_MAX
	global TRIAL_RHO
	global FILE_POT
	global FILE_DISTRO
	global FILE_OUTPUT_DISTRO
	global EPSILON_CHECK
	global EPSILON
	
	
	#define dictionary with default values, 
	#these are the values I get if I don't 
	#set a given parameter in the parm file
	
	parm_dict ={
			'BETA' : 		    1.,
			'L_MAX': 		   128,
			'TRIAL_RHO':		    1.,
			'FILE_POT':		   " ",
			'FILE_DISTRO':		   " ",
			'FILE_OUTPUT_DISTRO':	   " ",
			'EPSILON_CHECK': 	 10E-8,
			'EPSILON':		 10E-3,
					}
	
	#fill dictionary with parameters
	
	for line in INPUT:
		fields = line.split()
		parm_dict[fields[0]] = fields[1]
		
	#assign values to parameters from dictionary
	
	BETA = float(parm_dict['BETA'])
	L_MAX = int(parm_dict['L_MAX'])
	TRIAL_RHO = float(parm_dict['TRIAL_RHO'])
	FILE_POT = parm_dict['FILE_POT']
	FILE_DISTRO = parm_dict['FILE_DISTRO']
	FILE_OUTPUT_DISTRO = parm_dict['FILE_OUTPUT_DISTRO']
	EPSILON_CHECK = float(parm_dict['EPSILON_CHECK'])
	EPSILON = float(parm_dict['EPSILON'])
	
	INPUT.close()
		
########################################
####		OUTPUT		########
########################################	

def write_distro(filename,distro_):
	"writes spherical harmonics coefficient of distro to filename"
	OUTPUT = open(filename,'w')
	
	lmax = distro_.shape[1]
	
	for u in range(2):
		for l in range(0,lmax):
			for m in range(0,l+1):
				OUTPUT.write("%d\t%d\t%d\t%f\n" %(u,l,m,distro_[u,l,m]))
					
	
	OUTPUT.close()	

########################################
####	      WRAPPERS 		########
########################################

def compute_wz(pot_):
	"pythonic wrapper of the fortran routine."
	
	wz = np.array([0.], order='F')
	dfts.compute_wz(pot_,wz)
	return wz[0]
	
def compute_zeta(precon_new, precon_old = None):
	"pythonic wrapper of the fortran routine."
	
	zeta = np.array([0.], order='F')
	dfts.compute_zeta(precon_new, precon_old, zeta)
	return zeta[0]
	
def convo_product(pot_, distroA_, distroB_, BETA):
	"pythonic wrapper of the fortran routine."
	
	convo = np.array([0.], order='F')
	dfts.convolution_product(pot_, distroA_, distroB_, convo, BETA)
	return convo[0]
	
def integrate(field):
	"pythonic wrapper of the fortran routine."
	
	integral = np.array([0.], order='F')
	dfts.integrate(field,integral)
	
	return integral[0]
	
def convolve_product(pot_, fieldA_, fieldB_, BETA):
	"pythonic wrapper of the fortran routine."
	
	convo = np.array([0.], order='F')
	dfts.convolve_product(pot_,fieldA_,fieldB_,convo,BETA)
	
	return convo[0]

########################################
####	 MINIMIZE - UTILITY	########
########################################

def compute_muex(pot_, TRIAL_RHO):
	"computes the value of mu_ex: the calculation has been carried out numerically, see notes."
	return TRIAL_RHO*np.sqrt(4.*np.pi)*pot_[0]

def compute_grad(pot_, rho, TRIAL_RHO, BETA, MUEX):
	"computes the gradient of the grand canonical potential."
	
	rho_ = SHExpandDH(rho, sampling=2, norm=4)			# compute distro expansion in spherical harmonics
									#		     v
	grad_id = np.zeros(rho.shape, order='F')			#       define grad_id and grad_ex_
	grad_ex_= np.zeros(rho_.shape, order='F')			#		     v
	
	dfts.compute_grad_id(rho, grad_id, TRIAL_RHO, BETA, MUEX)	# compute grad_id = ln(rho) - ( ln(rho_trial) + beta*muex )
	dfts.convolve(pot_, rho_, grad_ex_, BETA)			# compute grad_ex_ by convolution theorem
	grad_ex = MakeGridDH(grad_ex_, sampling=2, norm=4)		# compute grad_ex from its spherical harmonics coefficients
									#		     v
	return grad_id + grad_ex					# return the total gradient as the sum of the two terms
	
def compute_psi(rho, grad, wz, BETA, precon_old = None, psi_old = None):
	"computes the conjugate gradient psi using the old preconditioners and the old psi. If these are not given, returns the preconditioners."
	
	if precon_old is None:
		precon_old = np.array([], order='F')
		psi_old = np.array([], order='F')
	
	precon_new = np.zeros(rho.shape, order='F')			# create vector for the new preconditioners
	dfts.compute_precon(rho, grad, wz, precon_new, BETA)		# compute precon in fortran routine
									#		v
	if precon_old.shape == (0,):					# if it's the first step (on the Qth step) psi_new = precon_new, 
		return precon_new, precon_new				# then remember precon_new as precon_old for the next step.
									#		x
	else:								#
		zeta = compute_zeta(precon_new, precon_old)		# compute the zeta offset if step is not first or Qth.
									#		v
	psi_new = precon_new + zeta*psi_old				# compute the new conjugate gradient psi
									#		v
	return psi_new, precon_new					# returns psi_new and precon_new to be stored in precon_old

def compute_step(pot_, rho, psi, grad, BETA):
	"computes the best step to use."
	
	fieldA = psi*grad						# compute vector of the products of psi and grad
	numerator = integrate(fieldA)					# compute the first element of the step
									#		v
	fieldB = psi*psi/rho						# compute the vector of the products psi^2/rho
	denominator = integrate(fieldB)					# compute the first term of the denominator
									#		v
									#		v	
	psi_ = SHExpandDH(psi, sampling=2, norm=4)			# coefficients of the SH expansion of psi
	denominator += convolve_product(pot_,psi_,psi_,BETA)		# compute the last term of the denominator
									#		v
	return numerator/denominator					# return the best step a num/den
	
def compute_check(grad, EPSILON_CHECK):
	"checks how close to grad=0 and returns true if we are far from EPSILON_CHECK."
	
	sq_grad = np.sum(grad*grad)
	return sq_grad > EPSILON_CHECK

def compute_grand_pot(rho,pot_,BETA,MUEX,TRIAL_RHO):
	"computes the grand potential value from the density profile."
	
	rho_ = SHExpandDH(rho, sampling=2, norm=4)
	
	Fid = integrate(rho*( np.log(rho/TRIAL_RHO) - BETA*MUEX ))

	Fex = convolve_product(pot_, rho_, rho_, BETA)
	
	return Fid + Fex

########################################
####	      MINIMIZE 		########
########################################

def minimize(pot_, rho, TRIAL_RHO, BETA, MUEX, EPSILON, EPSILON_CHECK ):
	"find optimal distribution"

	Q = 10										# define the number of step of conjugate gradient
	counter = 0									# define the counter
											#		v
	wz = compute_wz(pot_)								# compute w(0) once and for all
											#		v
											#		v
	grad = compute_grad(pot_,rho,TRIAL_RHO, BETA, MUEX)				# compute gradient at first step
	check= compute_check(grad, EPSILON_CHECK)					# check is true if norm squared of grad is > EPSILON_CHECK
											#		v
	while check:									# continue loop until check is false, that is grad is smaller than EPS_CHECK
		if counter%Q == 0:							# see if we should do a normal step instead of conj grad
			psi, precon = compute_psi(rho, grad, wz, BETA)			# compute psi and precon without conj grad
		else:									#		v
			psi, precon = compute_psi(rho, grad, wz, BETA, precon, psi)	# compute conj grad
											#		v
		counter += 1								# update counter
											#		v
		step = compute_step(pot_, rho, psi, grad, BETA)				# compute the best step
											#		v
		dfts.update_distro(rho, psi, step, EPSILON)				# update the distro
											#		v
		grad = compute_grad(pot_,rho,TRIAL_RHO, BETA, MUEX)			#compute gradient
		check= compute_check(grad, EPSILON_CHECK)				#check is true if norm squared of grad is > EPSILON_CHECK
											#		v
	return rho									#return the optimized rho
	
########################################
####		PAINT		########
########################################

def show(F):

	N = 2*L_MAX

	theta = np.linspace(0.,np.pi,N)
	phi   = np.linspace(0.,2*np.pi,2*N)

	PHI,THETA = np.meshgrid(phi, theta)
	
	fig, ax = plt.subplots(1,1,figsize=(10,5))
	ax.imshow(F, extent=(0,2*np.pi,0,np.pi), cmap='magma')
	ax.set(xlabel='longitude', ylabel='latitude');
	
	plt.show()
		
#######################################################################
####			MAIN					#######
#######################################################################
#	notice that _ indicates the spherical harmonics transform of 
#	the variable. So if a is a function of theta and phi, a_ will
#	be its spherical harmonic transform.
#######################################################################

def main():
	#####################INPUT#####################################
	try:
		FILE_PARM = sys.argv[1]
	except IndexError:
		print "A parameter file must be given from the command line."
		sys.exit(1)

	get_parm(FILE_PARM)	#initialize parameters

	#####################INIT######################################

	W_   = get_pot(FILE_POT,L_MAX)			

	rho_ = get_distro(FILE_DISTRO,L_MAX)
	
	MUEX = compute_muex(W_,TRIAL_RHO)
	
	N = L_MAX*2
	
	#####################MINIMIZE##################################

	rho  = MakeGridDH(rho_, sampling=2, norm=4)
	
	show(rho)
	
	rho  = minimize(W_, rho, TRIAL_RHO, BETA, MUEX, EPSILON, EPSILON_CHECK )

	show(rho)
	
	rho_ = SHExpandDH(rho, sampling=2, norm=4)
	
	C = sht.SHCoeffs.from_array(rho_) 

	fig, ax = C.plot_spectrum2d()
	
	write_distro(FILE_OUTPUT_DISTRO, rho_)


###############################RUN######################################

main()









