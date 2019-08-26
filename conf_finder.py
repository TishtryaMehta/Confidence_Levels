
import numpy as np
from matplotlib import pyplot as plt
import colorednoise as cn
from scipy.fftpack import fft
from scipy.stats import norm
from astroML.fourier import PSD_continuous
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=8, usetex=True)
from scipy import signal
from scipy.stats import chi2
from matplotlib.gridspec import GridSpec
from colorednoise import *
from scipy.optimize import curve_fit
import math


#------------------------------------------------------------
def broken_conf(x,step_len, *args, **kwargs):
	
	"""
	Generate confidence intervals for a signal containing white and 
	colored noise.
    Based on the algorithm in:
    S. Vaughan: 
    "A simple test for periodic signals in red noise" (2005)
    10.1051/0004-6361:20041453 

    Parameters:
    -----------
    x : float.
        The input signal. 
        Needs to be evenly spaced in time. 
        1-dimensional array. 
        Labelling assumes the signal to be on the scale of days, but can 
        have any units. 
    step_len : integer.
        Cadence of the input signal. 
        Need to implement a non-integer
        cadence in future. For now, only input integers. 
    P : float, optional.
        Probability wanted from confidence level.
        Should range from 0-1; i.e P = 0.85 is an 85% significance. 
        Default level is 0.95 (95%).
    show: bool,
		If 'True', returns a figure with the input signal, 
		power distribution of the periodogram, power against frequency,
		power against period, Nyquist freq, and associated 
		confidence level.
		The fit is given by a solid red line, and the confidence level 
		is given as a dashed red line. 
		If 'False', produces no illustrative output. 
    
    Returns
    -------
    out : array, array, array, float
         f				 The frequency array from the periodogram
         power_sig  	 The power array from the periodogram
         fit_model_ci 	 The fitted confidence interval array
         output_per 	 The largest period which exceeds the confidence
						 level
         
    Issues to be resolved:
    ---------
	If multiple peaks exceed the confidence level, output_per 
	automatically selects the largest periodicity of the set, rather the
	one which exceeds the confidence level by the greatest amount. 
         
    Example:
    ---------
    # Find the (most) significant peak in an input signal at a 95% 
    # confidence level.
    # Input signal is an array x, with cadence of 5 units. 
    
    >>> from conf_finder import broken_conf
    >>> broken_conf(x,5,P=0.95, show=True )
    
    Created 25th August 2019 by Tishtrya Mehta
    """
    
    
	sig_input=x
	step=step_len
	P=kwargs.get('P', None)
	show= kwargs.get('show', None)
	freq_max=[]	
	t=np.linspace(0, len(sig_input)*step, num=len(sig_input))
	t_hist=np.linspace(0,9,9000)										# NEED A SEPARATE TIME ARRAY FOR CADENCE ON HISTOGRAM
	if P == None:														# SET A CONFIDENCE LEVEL (DEFAULT 95)
		P=0.95
	if P>1.:
		print('The input confidence level should be between 0 and 1.0. For example, for a 99% level input: "P=0.99"')
		quit()							
	#-------------------------------------------------------------
	# COMPUTATION

	# CREATE SCIPY PERIODOGRAM
	f, power_sig = signal.periodogram(sig_input, fs=1./step, scaling='spectrum')
	power_sig=power_sig[0:int(len(sig_input)/2)]						# DISCARD THE NEGATIVE FREQUENCIES
	f=f[0:int(len(sig_input)/2)]
	#power_sig= power_sig/np.mean(power_sig)							# UNCOMMENT TO NORMALISE THE MEAN POWER TO EQUAL 1.

	# CREATE FITTING						
	f_log,power_log = zip(*[(math.log10(n), math.log10(g)) for n,g in zip(f, power_sig) if n>0.]) # DISCARD THE FIRST FREQUENCY (NYQUIST) AND ANY NEGATIVE VALUES

	# FIT THE POWER SPECTRA IN LINSPACE (SEE NOTES FOR LOGIC)
	# (SURPRESSES BIAS FOR THE LEAST SQUARE FITTING FOR THE LARGEST VALUES)
	def linbreak(x,m,C1,C2):
		return np.piecewise(x,[x<(C1-C2)/m,x>=(C1-C2)/m],[lambda x: C1 - m*x,C2])
		
	popt,pcov= curve_fit( linbreak, f_log, power_log, p0=[3.5,1.5,-0.5])# INITIAL GUESSES ARE BASED ON SIMULATED DATA
	fit_model_log=linbreak(f_log, *popt)								# I DONT KNOW WHAT POPT ACTUALLY MEANS, BUT I LIKE THE WORD. POPT! POP! POPT!
	fit_model=[ 10.**(n) for n in fit_model_log]						# MOVE MODEL BACK TO 'LOG' SCALES
	f=list([10.**(n) for n in f_log])							
	power_sig=[10.**(n) for n in power_log]
	
	nf= (len(sig_input))/2.	
	# FIND CONFIDENCE LEVELS
	ci= -2.*np.log(1.-(P)**(1/nf)) 										# USING P=1. -(1-np.exp(-x/2))**N. (SEE PAGE 6 OF VAUGHAN 2005 OR NOTES)
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------
	# IMPORTANT! HERE WE USE A CONF. LEVEL FOR A CHI-SQR 2 DOF DISTRIBUTION. THIS IS NOT NECESSAIRILY THE DISTRIBUTION OF AN ARTIFICIAL INPUT SIGNAL.
	# THIS IS THE ASSUMPTION WE TEST AGAINST TO ESTIMATE THE SIGNIFICANCE OF PEAKS AGAINST THE NULL HYPOTHESIS
	#---------------------------------------------------------------------------------------------------------------------------------------------------
	
	fit_model_ci=[i*ci*0.5 for i in fit_model]							# THE CI IS CALCULATED FOR A CHI-SQR 2 DOF DISTRIBUTION (WHICH HAS A MEAN OF 2). TO APPLY
																		# THE CI TO OUR DISTRIBUTION, WE MUST RENORMALISE SUCH THAT MEAN =1 (WHICH IS THE 0.5 FACTOR) 
																		# THEN MULTIPLY BY YFIT. SEE YOUR NOTES ( OR REREAD VAUGHAN OR PUGH 2017 (SEE BELOW) ) 
										
	for n in range(len(fit_model_ci)):									# COLLECT THE ARRAY OF POINTS WHICH HAVE POWER > CI. 
		if power_sig[n]> fit_model_ci[n]:
			freq_max.append(1./f[power_sig.index(power_sig[n])])

	#-------------------------------------------------------------
	#PLOTTING
	props = dict(boxstyle='square', facecolor='white', alpha=0.9)
	fig = plt.figure(figsize=(10,7))
	fig.subplots_adjust(hspace=0.45)
	gs=GridSpec(3,2)
	fs=11
	ax1 = fig.add_subplot(gs[0,:])
	ax2 = fig.add_subplot(gs[1,0])
	ax3 = fig.add_subplot(gs[1,1])
	ax4 = fig.add_subplot(gs[2,:])
	#-------------------------------------------------------------
	# PLOT INPUT SIGNAL
	ax1.set_title(r'Input signal ', fontsize=14)
	ax1.plot( t, sig_input,  color='k')
	#ax1.scatter( t, sig_input,  color='k', s=2)
	ax1.set_xlabel('$time \ [days]$' ,fontsize=fs)
	ax1.set_xlim(0,t[-1])
	ax1.set_ylabel(r'Amplitude' ,fontsize=fs)
	ax1.text(0.95,.8,r'$Nyquist \ Freq: {:04.2f} \ {:s}$'.format(1./(2.*step), 'days^{-1}'),fontsize=16, horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes,  bbox=props)

	#-------------------------------------------------------------
	# PLOT PERIODOGRAM (WITH CI)
	ax2.set_title(r'SciPy: Periodogram Routine', fontsize=14)
	ax2.loglog(f, power_sig,'-k')
	ax2.loglog(f, fit_model, 'r')
	ax2.loglog(f, fit_model_ci, 'r:', label='CI: {:04.1f}$\%$'.format(100.*P))
	ax2.set_ylim(min(power_sig[1:]), 10.*max(power_sig))
	ax2.set_xlabel(r'$Frequency \ [days^{-1}]$' ,fontsize=fs)
	ax2.set_ylabel(r'$Power \ [V^{2}]$' ,fontsize=fs)
	ax2.legend()

	#-------------------------------------------------------------
	# PLOT HISTOGRAM OF POWERS
	ax3.hist(power_sig, density=True,bins=int(len(x)/5.))
	ax3.set_xlabel(r'$Power \ [V^{2}]$' ,fontsize=fs)
	ax3.set_xlim(0,5)
	ax3.plot(t_hist,0.5*np.exp(-0.5*t_hist), 'r--',label=r'$\chi^{2}_{2} = \frac{1}{2} e^{\frac{-x}{2}}$' ) 
	# AGAIN, HERE, NOTE WN FOLLOWS A 0.5EXP(-0.5X) DISTRIBTION - NOT NECESSAIRILY THE CASE FOR ALL WN
	# INPUT AS SOMETIMES IT FOLLOWS THE EXP(-X') DISTRIBUTION- DEPENDS ON HOW THEYVE BEEN NORMALISED.
	ax3.legend()

	#-------------------------------------------------------------
	# PLOT PERIOD SPECTRA (WITH CI) 
	per=[float(1./n) for n in f]
	
	ax4.loglog(per, power_sig,'-k')
	#ax4.set_xlim(2,len(sig_input)*step/2)
	ax4.set_xlabel(r'$Period \ [days]$' ,fontsize=fs)
	ax4.set_ylabel(r'$Power \ [V^{2}]$' ,fontsize=fs)
	ax4.loglog(per, fit_model_ci, ':r')
	ax4.set_xlim(min(per),max(per))
	ax4.set_ylim(0.1*min(power_sig), 10*max(power_sig))
	ax4.plot(per, fit_model, 'r')
	if (len(freq_max) <1):
		print('No peaks above the {:04.1f}th percentile are detected'.format(P*100))
		output_per= 0.0
	else:	
		output_per= max(freq_max)
		print('The periods above the {:04.1f}th percentile are: {}'.format(P*100,freq_max))
		#textstrout = '\n'.join((
		#	'Input period: {:04.2f} days'.format(period),				# UNCOMMENT IF YOU KNOW THE INPUT PERIOD
		#	'Output period: {:04.2f} days'.format(max(freq_max))		
		#	))
		textstrout = 'Output period: {:04.2f} days'.format(output_per)# ONLY SELECTS THE GREATEST PERIODICITY
		ax4.text(0.70,.9,textstrout,fontsize=12, horizontalalignment='left', verticalalignment='top', transform=ax4.transAxes,  bbox=props)
		
	#-------------------------------------------------------------
	# SHOW FIGURE		
	if show==None or show==True:
		plt.show()
	return f, power_sig, fit_model_ci,output_per 
