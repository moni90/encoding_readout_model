utilities_path = './utilities';
addpath(genpath(utilities_path));

%folder to save simulations
data_folder = './simulation_output_files';

%folder and flags to save figures and figures data
save_figs_data = 1;
save_figs = 1;
save_path = './simulation_figures/';
%% code for models
n_simu = 10;
n_samples = 50;
decoder_type = 'lda';
stdev = 0.2; %std of single neurons
gamma_range = [.02:.02:.1];%[.02:.02:.5]; %range of signal to noise angle
rho_range = [.1:.1:.9];%[.1:.1:.9]; %correlations strength
equalization_type='';
dist_stim = sqrt(0.02)/2; %half distance between stimuli
n_workers = 1; %worker number for parallel computing
%%
N_range = 1; %number of neurons per sample

for idx_N = 1:length(N_range)
    N=N_range(idx_N);
    id_string = ['_v3_N' num2str(N) '_sigma_' num2str(stdev) '_d_' num2str(dist_stim*2)];
    tic;
    run_simu_suppl_model_figures_N_v3(N,n_simu,n_samples,decoder_type,dist_stim,stdev,gamma_range, rho_range, equalization_type, id_string, n_workers);
    toc;
    
    file_name = ['suppl_model_figure_simulation_output_lda_decoder_' num2str(n_simu) 'simu_' num2str(n_samples) 'samples' id_string '.mat'];
    
    %for which correlations strenght, signal to noise angle and consistency figures are generated
    rho_idx = 8;
    gamma_idx = 4;
    consistency_idx = 2;
    for which_panel = 4:11
        generate_figure_2(which_panel, data_folder, file_name, rho_idx, gamma_idx, consistency_idx, save_path, save_figs_data, save_figs);
    end
end
