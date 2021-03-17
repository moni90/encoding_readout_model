utilities_path = './utilities';
addpath(genpath(utilities_path));

%% code for models
N=3;
n_simu = 10;
n_samples = 50000;
decoder_type = 'lda';
stdev = 0.2;  %std of single neurons
rho_within = 0.8;
gamma_range = [0:.01:.25]'; %range of signal to noise angle
rho_range = [0.8]; %populationwise correlations strength
equalization_type='';
dist_stim = sqrt(0.02); %half distance between stimuli
n_workers = 1; %worker number for parallel computing
%%
cons_range = [.25 .5 .75 1];
for idx_cons = 1:length(cons_range)
    consistency_thr = cons_range(idx_cons);
    rho_range = rho_range(rho_range>0);
    id_string = ['_v7_N' num2str(N) '_rho_within_' num2str(rho_within) '_d_' num2str(dist_stim) '_cons_thr_' num2str(consistency_thr) '_temp'];
    tic;
    run_simu_suppl_model_figures_N_v7(N,n_simu,n_samples,decoder_type,dist_stim,stdev,gamma_range, rho_range, consistency_thr, equalization_type, id_string, n_workers);
    toc;
end
