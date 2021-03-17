utilities_path = './utilities';
addpath(genpath(utilities_path));

%folder to save simulations
data_folder = './simulation_output_files';

%folder and flags to save figures and figures data
save_figs_data = 1;
save_figs = 1;
save_path = './simulation_figures/';
%% code for models

n_simu = 10;%100;
n_samples = 50;%300000;
decoder_type = 'lda';
stdev = 0.2; %std of single neurons
rho_within = 0.8;
gamma_range = [0:.01:.5]'; %range of signal to noise angle
Ni_range = [0.5:0.1:0.8];%[0.1:0.05:0.7]; %populationwise correlations strength
equalization_type='';
dist_stim = 0.15; %half distance between stimuli
n_workers = 1; %worker number for parallel computing
%%
N_range = 2;%[5 10 15 20]; %number of neurons per sample

for idx_N = 1:length(N_range)
    N=N_range(idx_N);
    rho_range = ((2*N)*Ni_range - 1 )/(2*N-1) ; %model v1
    rho_range = rho_range(rho_range>0);
    id_string = ['_v3_N' num2str(N) '_rho_within_' num2str(rho_within) '_d_' num2str(dist_stim)];
    run_simu_suppl_model_figures_N_v3(N,n_simu,n_samples,decoder_type,dist_stim,stdev,gamma_range, rho_range, equalization_type, id_string, n_workers);
    
    file_name = ['suppl_model_figure_simulation_output_lda_decoder_' num2str(n_simu) 'simu_' num2str(n_samples) 'samples' id_string '.mat'];
    
    consistency_idx = 3; %for which consistency figures are generated. 1 for cons-enhanced readout with highest consistency
    for which_panel = 1:9
        generate_figure_simulations(which_panel, data_folder, file_name, N, consistency_idx, save_path, save_figs_data, save_figs);
    end
end

