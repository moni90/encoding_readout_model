# encoding_readout_model
Supplementary code for the paper "Correlations enhance the behavioral readout of neural population activity in association cortex" by M. Valente, G. Pica, G. Bondanelli, M. Moroni, C. A. Runyan, A. S. Morcos, C. D. Harvey, S. Panzeri, (2021).

The repository contains the code to reproduce Figure 2 d, e, f right, g, h, I right, j, k, Figure 3 and Extended Data Figure 4.
It is structured as follows:
-	Encoding_readout_model: scripts in this folder are the main scripts for the generation of simulated data. In particular, the main scripts to run is main_model_pipeline.m;
  - Utilities: auxiliary scripts necessary to run the analyses and generate figures, but not specific to the encoding-readout model are saved in this folder and its subfolders
  - Simulation_output_files: simulated data are automatically saved in this subfolder as .fig, .png or .svg files;
  - Simulation_figures: users can save the final figures in this subfolder;
     - Data_for_figures: data necessary to generated figures can be saved in this subfolder in .mat or .csv format;
  

**main_model_pipeline.m**

Running this script users can generate simulated data. Parameters to be set are:
-	save_figs, save_figs_data: binary flags to save the generated figures and their data
-	n_simu: number of simulations
-	n_samples: number of samples for each simulation
-	decoder_type: decoder type to compute the encoded stimulus ('lda')
-	stdev: standard deviation of each neuron (to set the diagonal values of the covariance matrix)
-	rho_within: correlation strength of neurons in the same feature
-	gamma_range: range of signal-to-noise angle to explore (in \pi units)
-	Ni_range: range of correlation strength to explore
-	Equalization_type: type of equalization ('','joint_posterior', 'corr_con_decoding', 'correct&error_trials')
-	Dist_stim: distance between the average response to the two stimuli
-	N_workers: number of workers for the parallel computing
-	N_range: range of number of neurons in each feature to explore.
-	Consistency_idx: consistency index to use to generate figures

For Figure 3 and Extended Data Figure 4
-	n_simu: 100
-	n_samples: 300,000
-	decoder_type: 'lda'
-	stdev: 0.2
-	rho_within: 0.8
-	gamma_range: [0:.01:.5]'
-	Ni_range: [0.5:0.1:0.8]
-	Equalization_type: ''
-	Dist_stim: 0.15
-	N_workers: 30
-	N_range: 20 (Fig. 3) or 10 (ED Fig. 4)
-	Consistency_idx: 3

**main_generate_Fig2.m**

Running this script users can reproduce Figure 2 d, e, f right, g, h, I right, j, k.
