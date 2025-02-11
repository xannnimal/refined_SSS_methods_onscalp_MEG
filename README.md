# refined_SSS_methods_onscalp_MEG
Analysis scripts for "Refined SSS Methods for on-scalp MEG systems"

Description of MATLAB files
* "main_methods.mat" calculates metric comparisons and current dipole simulation results for each helmet system. 
* "Kernel_system.mat" constructs each SSS method variant and compares the reconstruction of evoked data collectected at University of Washington and processed using MNE-Python
* "UCL_system.mat" constructs each SSS method variant and compares the reconstruction of raw data collectected at UCL
* "sensor_deviation.mat" uses the MEGIN/Elekta Neuromag system to simulated a deviated on-scalp OPM system of magnetometers with all SSS variants and simulated data
* "SSS_function" folder contains necessary function for file i/o as well as for SSS,
* "multi_sss.mat" is implementation of novel mSSS

Necessary packages and other files
* MNE-Matlab
* Kernel Flux OPM raw and evoked data from UW is available online at OSF: https://osf.io/teygz/
* Spheroidal harmonic functions "spm_ipharm.mat" and "spm_epharm.mat" available from Tim Tierney et al. through SPM on GitHub at https://github.com/spm
* For the MEGIN/Elekta Neuromag system, the raw file "sample_audvis_raw.fif" for SQUID sensor geometry is obtained on MNE-Python
* Data and sensor geometry for QuSpin helmet at Nottingham can be obtained via Niall Holmes et al. with corresponding reference https://ieeexplore.ieee.org/document/10685146
*UCL data and corresponding tutorial is publicaly available online on MNE-Python at https://mne.tools/stable/auto_tutorials/preprocessing/80_opm_processing.html
