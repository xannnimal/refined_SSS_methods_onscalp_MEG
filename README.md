# refined_SSS_methods_onscalp_MEG
Analysis scripts for "Refined SSS Methods for on-scalp MEG systems"

"main_methods.mat" calculates metric comparisons and current dipole simulation results for each helmet system. 
"Kernel_system.mat" constructs each SSS method variant and compares the reconstruction of evoked data collectected at University of Washington and processed using MNE-Python. Data is available online at OSF: https://osf.io/teygz/
"UCL_system.mat" constructs each SSS method variant and compares the reconstruction of raw data collectected at UCL, data and corresponding tutorial is publicaly available online on MNE-Python "Preprocessing optically pumped magnetometer (OPM) MEG data"
"sensor_deviation.mat" uses the MEGIN/Elekta Neuromag system to simulated a deviated on-scalp OPM system of magnetometers with all SSS variants and simulated data
"SSS_function" folder contains necessary function for file i/o as well as for SSS, mSSS calculation
