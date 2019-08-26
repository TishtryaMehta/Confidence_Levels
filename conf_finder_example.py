""" Example: Generate confidence intervals for a signal containing white and colored noise """
import numpy as np
from colorednoise import *
import math
from conf_finder import broken_conf

#------------------------------------------------------------

# GENERATE SIGNAL WITH WHITE AND COLORED NOISE
# CONTAINS A PERIODICITY OF 20 UNITS, CADENCE OF 3 UNITS. 

np.random.seed(0)
period=20.
end=500
step=3
t=np.linspace(0, end, num=int(end/step))
white_noise = np.random.normal( 0,1,size=t.shape) 			
sig_raw = np.sin((2*math.pi)*t/period) 
col_noise= powerlaw_psd_gaussian(2, len(t))
sig_input=.5*sig_raw + .5*col_noise + 0.5*white_noise 


y= broken_conf(sig_input,step,P=0.95, show=True )
