Steps to produce Figure 0.2 in Appendix D (important - must complete steps in the order listed)

1. Run the scripts "K_function_osc_1.m" and "K_function_osc_2.m" to generate the K functions for oscillators 1 and 2.
   - output: (a) plots of the K functions
	     (b) .mat files containing the K function data (up to order 15)

2. Run the scripts "ResponseCurves_osc_1.m" and "ResponseCurves_osc_2.m" to generate the PRC and IRC for oscillators 1 and 2.
   - output: (a) plots of the PRC and IRC functions
	     (b) .mat files containing the PRC and IRC data (up to order 8)

3. Run the script "Frequency_matching_curve_synch.m" to produce the plots in Figure 0.2 (left panel) 
   - output: (a) left panel of Figure 0.2 in Appendix D

4. To generate the right panel of Figure 0.2, run the script "Period_curve_synch.m"
	- NOTE: this script takes several minutes to run
	- NOTE: the main figure and inset were generate separately
		- the default code generates the main figure
		- to generate the inset, change L25 to "delta = linspace(-0.02,0.02,10);"
	- output: (a) right panel of Figure 0.2 in Appendix D
