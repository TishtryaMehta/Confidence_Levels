# Confidence_Levels
My tailored Python code for testing the significance of peaks in a Fourier spectrum by fitting of a broken power law.

This code builds on the techniques outlined by S. Vaughan ("A simple test for periodic signals in red noise" (2005), 10.1051/0004-6361:20041453 ), and C. Pugh ("Significance testing for quasi-periodic pulsations in solar and stellar flares" (2017), 10.1051/0004-6361/201730595), and mirrors the FFT_alpha.pro IDL code, written by S. Afinogentov (github.com/Sergey-Anfinogentov -See the EMD_conf repos.) 

The code tests the significance of peaks against the null hypothesis; that the peak is the result of chi-squared 2 d.o.f white noise. This uses the relation x = -2* log(1 - P^(1/N)), where P is the level of confidence sought. 

As white noise dominates over high frequencies, and colored noise over lower frequencies, the spectra is fitted with a broken power law, which the confidence level is scaled by. 

The example code makes use of the "colorednoise.py" code found here (https://github.com/felixpatzelt). 
