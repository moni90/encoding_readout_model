function out=simu_consistency_readout_N_v3(N,n_samples,mu_s0,mu_s1,stdev1_s0,stdev2_s0,stdev1_s1,stdev2_s1,rho_s0,rho_s1,rho_within,decoder_type,fitting_fraction,...
    alpha_readout,readout_strength,readout_strength_params,equalization_type,display_figures)

% n_samples: number of samples per stimulus; mu_s0, mu_s1, stdev_r1,
% stdev2, rho_s0, rho_s1: parameters defining the bivariate gaussian
% distributions to be simulated; decoder_type: 'analytical','lda','svm';
% fitting_fraction: fraction of trials used for decoder fitting;
% alpha_readout: angle the type of readout (pi/2: squared readout,
% pi:optimal readout); readout_strength: 'ideal','strong','custom' defines
% the strength of stimulus to choice conversion in each readout region;
% readout_strength_params (only for 'custom' readout strength): user
% defined readout strength; equalization_type: 'joint_posterior',
% 'corr_con_decoding', 'correct&error_trials' different types of
% equalization of stimulus information in consistent and inconsistent
% trials

%% Simulate Gaussian neural responses
% Assumption of the model: stdev1 (standard deviation along feature 1) and
% stdev2 (standard deviation along feature 2) are the same for both stimuli
% s0 and s1 (not necessarily equal between each other)
sigma_s0 = [rho_within*stdev1_s0*stdev1_s0*ones(N,N) rho_s0*stdev1_s0*stdev2_s0*ones(N,N);...
    rho_s0*stdev1_s0*stdev2_s0*ones(N,N) rho_within*stdev2_s0*stdev2_s0*ones(N,N)].*(1-eye(2*N));
sigma_s0 = diag([(stdev1_s0^2)*ones(N,1); (stdev2_s0^2)*ones(N,1)])+sigma_s0;%/(N);
sigma_s1 = [rho_within*stdev1_s1*stdev1_s1*ones(N,N) rho_s1*stdev1_s1*stdev2_s1*ones(N,N);...
    rho_s1*stdev1_s1*stdev2_s1*ones(N,N) rho_within*stdev2_s1*stdev2_s1*ones(N,N)].*(1-eye(2*N));
sigma_s1 = diag([(stdev1_s1^2)*ones(N,1); (stdev2_s1^2)*ones(N,1)])+sigma_s1;%/(N);

sigma_sh_s0 = rho_within*[stdev1_s0*stdev1_s0*ones(N,N) zeros(N,N);...
    zeros(N,N) stdev2_s0*stdev2_s0*ones(N,N)].*(1-eye(2*N));
sigma_sh_s0 = diag([(stdev1_s0^2)*ones(N,1); (stdev2_s0^2)*ones(N,1)])+sigma_sh_s0;%/(N);
sigma_sh_s1 = rho_within*[stdev1_s1*stdev1_s1*ones(N,N) zeros(N,N);...
    zeros(N,N) stdev2_s1*stdev2_s1*ones(N,N)].*(1-eye(2*N));
sigma_sh_s1 = diag([(stdev2_s1^2)*ones(N,1); (stdev2_s1^2)*ones(N,1)])+sigma_sh_s1;%/(N);

% sigma_s0 = build_sigma_N(stdev1_s0,stdev2_s0,rho_s0);
% sigma_s1 = build_sigma_N(stdev1_s1,stdev2_s1,rho_s1);
% sigma_sh_s0 = build_sigma_N(stdev1_s0,stdev2_s0,0);
% sigma_sh_s1 = build_sigma_N(stdev1_s1,stdev2_s1,0);

r_s0 = mvnrnd(mu_s0,sigma_s0,n_samples); % response to stimulus 1
r_s1 = mvnrnd(mu_s1,sigma_s1,n_samples); % response to stimulus 2

r_all=[r_s1;r_s0];
s_all=[ones(n_samples,1);zeros(n_samples,1)];

train_idxs = get_training_trials_indices(s_all,fitting_fraction);
test_idxs = setdiff([1:2*n_samples]',train_idxs);

train_trials = false(2*n_samples,1); train_trials(train_idxs)=true;
test_trials = false(2*n_samples,1); test_trials(test_idxs)=true;

%project on stimulus axis
%population 8
r_pool1 = [r_s1(:,1:N); r_s0(:,1:N)];
stim_axis_1 = nanmean(r_s1(:,1:N),1)-nanmean(r_s0(:,1:N),1);
stim_axis_1 = stim_axis_1'/norm(stim_axis_1);
r_pool1 = r_pool1*stim_axis_1;
r_pool2 = [r_s1(:,N+1:end); r_s0(:,N+1:end)];
stim_axis_2 = nanmean(r_s1(:,N+1:end),1)-nanmean(r_s0(:,N+1:end),1);
stim_axis_2 = stim_axis_2'/norm(stim_axis_2);
r_pool2 = r_pool2*stim_axis_2;
r_all_2d=[r_pool1 r_pool2];

%% Shuffle data
r_all_sh=nan*ones(size(r_all));

for this_s = 0:1
    
    this_class_trials = (s_all==this_s) & train_trials;
    
    r_all_sh(this_class_trials,1:N) = datasample(r_all(this_class_trials,1:N),nnz(this_class_trials),'Replace',false);
    r_all_sh(this_class_trials,N+1:end) = datasample(r_all(this_class_trials,N+1:end),nnz(this_class_trials),'Replace',false);
    
    this_class_trials = (s_all==this_s) & test_trials;
    
    r_all_sh(this_class_trials,1:N) = datasample(r_all(this_class_trials,1:N),nnz(this_class_trials),'Replace',false);
    r_all_sh(this_class_trials,N+1:end) = datasample(r_all(this_class_trials,N+1:end),nnz(this_class_trials),'Replace',false);
    
end

%% Run 2D and 1D decoder

switch decoder_type
    
    case 'analytical'
        
        % the decoding boundary is defined analytically from the
        % distributions parameters (not from the available samples)
        
        % when selecting this type of decoding, consistency defined on
        % posterior probs is the same as consistency defined on responses
        % directly
        
        if not(isequal(rho_s0,rho_s1) & isequal(stdev1_s0,stdev1_s1) & isequal(stdev2_s0,stdev2_s1))
            error('Error: analytical decoder requires equal covariance matrix for the two stimuli');
        end
        
        class1_idx = 2;
        
        [posterior_probs_j,posterior_probs_1,posterior_probs_2]=...
            compute_stimulus_posterior_analytical(r_all(test_idxs,:),mu_s0,mu_s1,sigma_s0);
        
        [K,L]=compute_lda_theoretical_coeffs(mu_s0',mu_s1',sigma_s0);
        
        [posterior_probs_j_sh,posterior_probs_1_sh,posterior_probs_2_sh]=...
            compute_stimulus_posterior_analytical(r_all_sh(test_idxs,:),mu_s0,mu_s1,sigma_sh_s0);
        
        [K_sh,L_sh]=compute_lda_theoretical_coeffs(mu_s0',mu_s1',sigma_sh_s0);
        
        decoding_dir=L;
        decoding_dir_sh=L_sh;
        
        f_db = @(x,y) L(1)*x+L(2)*y+K;
        f_db_sh = @(x,y) L_sh(1)*x+L_sh(2)*y+K_sh;
        
    case 'lda'
        
        class1_idx = 2;
        
        [label, posterior_probs_j, ~,  K, L]=lda_classify(r_all, s_all, [], train_idxs);
        [~, posterior_probs_1]=lda_classify(r_all(:,1:N), s_all, [], train_idxs);
        [~, posterior_probs_2]=lda_classify(r_all(:,N+1:end), s_all, [], train_idxs);
        
        [~, posterior_probs_j_sh, ~,  K_sh, L_sh]=lda_classify(r_all_sh, s_all, [], train_idxs);
        [~, posterior_probs_1_sh]=lda_classify(r_all_sh(:,1:N), s_all, [], train_idxs);
        [~, posterior_probs_2_sh]=lda_classify(r_all_sh(:,N+1:end), s_all, [], train_idxs);
        
        decoding_dir=L;
        decoding_dir_sh=L_sh;
        
        f_db = @(x,y) L(1)*x+L(2)*y+K;
        f_db_sh = @(x,y) L_sh(1)*x+L_sh(2)*y+K_sh;
        
    case 'svm'
        
        class1_idx = 1;
        
        [~, ~, ~, posterior_probs_j,model_opt]=svm_classification_linear(r_all, s_all, test_idxs);
        [~, ~, ~, posterior_probs_1]=svm_classification_linear(r_all(:,1), s_all, test_idxs);
        [~, ~, ~, posterior_probs_2]=svm_classification_linear(r_all(:,2), s_all, test_idxs);
        
        [~, ~, ~, posterior_probs_j_sh,model_opt_sh]=svm_classification_linear(r_all_sh, s_all, test_idxs);
        [~, ~, ~, posterior_probs_1_sh]=svm_classification_linear(r_all_sh(:,1), s_all, test_idxs);
        [~, ~, ~, posterior_probs_2_sh]=svm_classification_linear(r_all_sh(:,2), s_all, test_idxs);
        
        w = compute_decoding_boundary_normal(model_opt);
        w_sh = compute_decoding_boundary_normal(model_opt_sh);
        b = -model_opt.rho;
        
        decoding_dir = w;
        decoding_dir_sh = w_sh;
        
        f_db = @(x,y) w(1)*x+w(2)*y+b;
        f_db_sh = @(x,y) w(1)*x+w(2)*y+b;
        
end

% select posterior probabilities of s=1
ps_joint = posterior_probs_j(:,class1_idx);
ps_1 = posterior_probs_1(:,class1_idx);
ps_2 = posterior_probs_2(:,class1_idx);

ps_joint_sh = posterior_probs_j_sh(:,class1_idx);
ps_1_sh = posterior_probs_1_sh(:,class1_idx);
ps_2_sh = posterior_probs_2_sh(:,class1_idx);

s=s_all(test_idxs,:);
r=r_all(test_idxs,:);
r_sh=r_all_sh(test_idxs,:);
r_2d = r_all_2d(test_idxs,:);

%% Define sp and consistency

sp = ps_joint>.5;
sp_sh = ps_joint_sh>.5;

correct_s_decoded = sp==s;
correct_s_decoded_shuffle = sp_sh==s;

[~,m1,m2,region]=conical_readout(ps_1,ps_2,alpha_readout);
consistent = region == 1 | region == 2;
inconsistent = region == 0;

bad_trials = (region==1 & sp==1) | (region==2 & sp==0);

[~,~,~,region]=conical_readout(ps_1_sh,ps_2_sh,alpha_readout);
consistent_shuffle = region == 1 | region == 2;
inconsistent_shuffle = region == 0;

bad_trials_sh = (region==1 & sp_sh==1) | (region==2 & sp_sh==0);

%% Perform equalization of stimulus information 

% between consistent and inconsistent trials

switch equalization_type
    case 'joint_posterior'
        n_bins = 30;
        [consistent, inconsistent] = subsample_at_fixed_posterior(ps_joint, consistent, n_bins);
        [consistent_shuffle, inconsistent_shuffle] = subsample_at_fixed_posterior(ps_joint_sh, consistent_shuffle, n_bins);
    case 'corr_con_decoding'
        [consistent, inconsistent]= subsample_by_correctness_and_consistency(correct_s_decoded,consistent);
        [consistent_shuffle, inconsistent_shuffle]= subsample_by_correctness_and_consistency(correct_s_decoded_shuffle,consistent_shuffle);
    case 'correct&error_trials'
        [consistent, inconsistent]= subsample_corr_err_trials_to_equalize_S_decoding_perf(correct_s_decoded, consistent, ~consistent);
        [consistent_shuffle, inconsistent_shuffle]= subsample_corr_err_trials_to_equalize_S_decoding_perf(correct_s_decoded_shuffle,consistent_shuffle,~consistent_shuffle);
end

selected = consistent | inconsistent;
selected_shuffle = consistent_shuffle | inconsistent_shuffle;

% between correlated and uncorrelated trials

if strcmp(equalization_type,'corr&shuffled_trials')
    [selected] = subsample_trials_to_eqShuffledInfo(correct_s_decoded,correct_s_decoded_shuffle);
    selected_shuffle = selected;
end


%% Simulate choices

choice=define_choice_from_sp_and_consistency(sp,consistent,readout_strength,readout_strength_params);
choice_sh=define_choice_from_sp_and_consistency(sp_sh,consistent_shuffle,readout_strength,readout_strength_params);

avg_readout_strength = nnz((sp == choice) & selected)/nnz(selected);
choice_subopt=define_choice_from_sp_and_consistency(sp,ones(size(consistent)),'custom',[avg_readout_strength avg_readout_strength .5 .5]);
choice_subopt_sh=define_choice_from_sp_and_consistency(sp_sh,ones(size(consistent_shuffle)),'custom',[avg_readout_strength avg_readout_strength .5 .5]);

correct_behavior = choice==s;
correct_behavior_sh = choice_sh==s;

correct_behavior_subopt = choice_subopt==s;
correct_behavior_subopt_sh = choice_subopt_sh==s;

%% Compute angles

signal_dir = compute_signal_correlation_direction(r(s==0,:),r(s==1,:)); % points towards class whose label is 1!
[noise_dir_s0, noise_dir_s1]  = compute_noise_correlation_direction(r(s==0,:),r(s==1,:));
consistency_dir = sign(signal_dir); % this is not normalized!
gamma_s0 = md_angle(signal_dir,noise_dir_s0);
gamma_s1 = md_angle(signal_dir,noise_dir_s1);
gamma = abs(circ_mean([gamma_s0;gamma_s1]));
theta = md_angle(signal_dir,decoding_dir);
% signal_dir_angle = atan2(signal_dir(N+1:end),signal_dir(1:N));
% decoding_dir_angle = atan2(decoding_dir(2),decoding_dir(1));
% theta = wrapToPi(decoding_dir_angle - signal_dir_angle);
newtheta = md_angle(consistency_dir,decoding_dir);

signal_dir_sh = compute_signal_correlation_direction(r_sh(s==0,:),r_sh(s==1,:));
noise_dir_sh = compute_noise_correlation_direction(r_sh(s==0,:),r_sh(s==1,:));
consistency_dir_sh = sign(signal_dir_sh); % this is not normalized!
gamma_sh = md_angle(signal_dir_sh,noise_dir_sh);
theta_sh = md_angle(signal_dir_sh,decoding_dir_sh);
% signal_dir_angle_sh = atan2(signal_dir_sh(2),signal_dir_sh(1));
% decoding_dir_angle_sh = atan2(decoding_dir_sh(2),decoding_dir_sh(1));
% theta_sh = wrapToPi(decoding_dir_angle_sh - signal_dir_angle_sh);
newtheta_sh = md_angle(consistency_dir_sh,decoding_dir_sh);

signal_dir_2d = compute_signal_correlation_direction(r_2d(s==0,:),r_2d(s==1,:)); % points towards class whose label is 1!
[noise_dir_s0_2d, noise_dir_s1_2d]  = compute_noise_correlation_direction(r_2d(s==0,:),r_2d(s==1,:));
consistency_dir_2d = sign(signal_dir_2d); % this is not normalized!
gamma_s0_2d = md_angle(signal_dir_2d,noise_dir_s0_2d);
gamma_s1_2d = md_angle(signal_dir_2d,noise_dir_s1_2d);
gamma_2d = abs(circ_mean([gamma_s0_2d;gamma_s1_2d]));
%% Compute deltaI_diag (Averbeck et al.2006)
% here it is computed as fraction of correctly decoded stimuli

% decode using shuffled decision boudnary
sp_corr = feval(f_db,r(:,1:N),r(:,N+1:end))>0;
sp_diag = feval(f_db_sh,r(:,1:N),r(:,N+1:end))>0;
I_corr = nnz((sp_corr == s) & selected)/nnz(selected);
I_diag = nnz((sp_diag == s) & selected)/nnz(selected);
deltaI_diag = I_corr - I_diag;

%% Compute correlation coefficient in correct/incorrect behavior trials
% corr_coeff_s0_c0 = corrcoef(r(s==0 & choice==0 & selected,1),r(s==0 & choice==0 & selected,2));
% corr_coeff_s0_c1 = corrcoef(r(s==0 & choice==1 & selected,1),r(s==0 & choice==1 & selected,2));
% corr_coeff_s1_c0 = corrcoef(r(s==1 & choice==0 & selected,1),r(s==1 & choice==0 & selected,2));
% corr_coeff_s1_c1 = corrcoef(r(s==1 & choice==1 & selected,1),r(s==1 & choice==1 & selected,2));

% out.corr_coeff_BP_correct = nnz(s==0 & correct_behavior & selected)/nnz(correct_behavior & selected) * corr_coeff_s0_c0 + nnz(s==1 & correct_behavior & selected)/nnz(correct_behavior & selected) * corr_coeff_s1_c1;
% out.corr_coeff_BP_error   = nnz(s==0 & not(correct_behavior) & selected)/nnz(not(correct_behavior) & selected) * corr_coeff_s0_c1 + nnz(s==1 & not(correct_behavior) & selected)/nnz(not(correct_behavior) & selected) * corr_coeff_s1_c0;

[pair_corr_coeff_BP_correct,pair_corr_coeff_BP_error,pop_corr_coeff_BP_correct,pop_corr_coeff_BP_error, n_all_selected_trials] = compute_correlation_coefficient_correct_error_N(s(selected),choice(selected),r(selected,1:N),r(selected,N+1:end));

out.pair_corr_coeff_BP_correct = pair_corr_coeff_BP_correct;
out.pair_corr_coeff_BP_error = pair_corr_coeff_BP_error;
out.pop_corr_coeff_BP_correct = pop_corr_coeff_BP_correct;
out.pop_corr_coeff_BP_error = pop_corr_coeff_BP_error;

[pair_corr_coeff_BP_correct_subopt,pair_corr_coeff_BP_error_subopt,pop_corr_coeff_BP_correct_subopt,pop_corr_coeff_BP_error_subopt, n_all_selected_trials] = compute_correlation_coefficient_correct_error_N(s(selected),choice_subopt(selected),r(selected,1:N),r(selected,N+1:end));

out.pair_corr_coeff_BP_correct_subopt = pair_corr_coeff_BP_correct_subopt;
out.pair_corr_coeff_BP_error_subopt = pair_corr_coeff_BP_error_subopt;
out.pop_corr_coeff_BP_correct_subopt = pop_corr_coeff_BP_correct_subopt;
out.pop_corr_coeff_BP_error_subopt = pop_corr_coeff_BP_error_subopt;

%% OUTPUTS

%% Decoding and behavioral performance

out.fraction_correct_S_decoded = nnz(correct_s_decoded&selected)/nnz(selected);
out.fraction_correct_S_decoded_shuffle = nnz(correct_s_decoded_shuffle&selected_shuffle)/nnz(selected_shuffle);

out.fraction_bad_trials = nnz(bad_trials&selected)/nnz(selected);
out.fraction_bad_trials_shuffle = nnz(bad_trials_sh&selected_shuffle)/nnz(selected_shuffle);

out.fraction_consistent = nnz(consistent&selected)/nnz(selected);
out.fraction_consistent_shuffle = nnz(consistent_shuffle&selected_shuffle)/nnz(selected_shuffle);



out.fraction_correct_behavior = nnz(correct_behavior&selected)/nnz(selected);
out.fraction_correct_behavior_subopt = nnz(correct_behavior_subopt&selected)/nnz(selected);
out.fraction_correct_behavior_shuffle = nnz(correct_behavior_sh&selected_shuffle)/nnz(selected_shuffle);
out.fraction_correct_behavior_subopt_shuffle = nnz(correct_behavior_subopt_sh&selected_shuffle)/nnz(selected_shuffle);

out.fraction_consistent_BP_correct        = nnz(correct_behavior       &consistent&selected)/nnz(correct_behavior       &selected);
out.fraction_consistent_BP_correct_subopt = nnz(correct_behavior_subopt&consistent&selected)/nnz(correct_behavior_subopt&selected);
out.fraction_consistent_BP_error          = nnz(not(correct_behavior)       &consistent&selected)/nnz(not(correct_behavior)       &selected);
out.fraction_consistent_BP_error_subopt   = nnz(not(correct_behavior_subopt)&consistent&selected)/nnz(not(correct_behavior_subopt)&selected);

out.gamma = gamma;
out.gamma_2d = gamma_2d;
out.gamma_sh = gamma_sh;
out.theta = theta;
out.theta_sh = theta_sh;
out.newtheta = newtheta;
out.newtheta_sh = newtheta_sh;

out.deltaI_diag = deltaI_diag;
out.I_diag = I_diag;

out.avg_readout_strength = avg_readout_strength;

% out.BP_correct = nnz(correct_behavior&correct_s_decoded&selected)/nnz(correct_s_decoded&selected);
% out.BP_error = nnz(correct_behavior&not(correct_s_decoded)&selected)/nnz(not(correct_s_decoded)&selected);
% out.BP_correct_subopt = nnz(correct_behavior_subopt&correct_s_decoded&selected)/nnz(correct_s_decoded&selected);
% out.BP_error_subopt = nnz(correct_behavior_subopt&not(correct_s_decoded)&selected)/nnz(not(correct_s_decoded)&selected);
% 
% out.BP_correct_con = nnz(correct_behavior&correct_s_decoded&consistent&selected)/nnz(correct_s_decoded&consistent&selected);
% out.BP_correct_inc = nnz(correct_behavior&correct_s_decoded&not(consistent)&selected)/nnz(correct_s_decoded&not(consistent)&selected);
% out.BP_error_con = nnz(correct_behavior&not(correct_s_decoded)&consistent&selected)/nnz(not(correct_s_decoded)&consistent&selected);
% out.BP_error_inc = nnz(correct_behavior&not(correct_s_decoded)&not(consistent)&selected)/nnz(not(correct_s_decoded)&not(consistent)&selected);
% out.BP_correct_con_subopt = nnz(correct_behavior_subopt&correct_s_decoded&consistent&selected)/nnz(correct_s_decoded&consistent&selected);
% out.BP_correct_inc_subopt = nnz(correct_behavior_subopt&correct_s_decoded&not(consistent)&selected)/nnz(correct_s_decoded&not(consistent)&selected);
% out.BP_error_con_subopt = nnz(correct_behavior_subopt&not(correct_s_decoded)&consistent&selected)/nnz(not(correct_s_decoded)&consistent&selected);
% out.BP_error_inc_subopt = nnz(correct_behavior_subopt&not(correct_s_decoded)&not(consistent)&selected)/nnz(not(correct_s_decoded)&not(consistent)&selected);

out.BP_correct = nnz(correct_behavior&correct_s_decoded&selected)/nnz(selected);
out.BP_error = nnz(correct_behavior&not(correct_s_decoded)&selected)/nnz(selected);
out.BP_correct_subopt = nnz(correct_behavior_subopt&correct_s_decoded&selected)/nnz(selected);
out.BP_error_subopt = nnz(correct_behavior_subopt&not(correct_s_decoded)&selected)/nnz(selected);

out.BP_correct_con = nnz(correct_behavior&correct_s_decoded&consistent&selected)/nnz(selected);
out.BP_correct_inc = nnz(correct_behavior&correct_s_decoded&not(consistent)&selected)/nnz(selected);
out.BP_error_con = nnz(correct_behavior&not(correct_s_decoded)&consistent&selected)/nnz(selected);
out.BP_error_inc = nnz(correct_behavior&not(correct_s_decoded)&not(consistent)&selected)/nnz(selected);
out.BP_correct_con_subopt = nnz(correct_behavior_subopt&correct_s_decoded&consistent&selected)/nnz(selected);
out.BP_correct_inc_subopt = nnz(correct_behavior_subopt&correct_s_decoded&not(consistent)&selected)/nnz(selected);
out.BP_error_con_subopt = nnz(correct_behavior_subopt&not(correct_s_decoded)&consistent&selected)/nnz(selected);
out.BP_error_inc_subopt = nnz(correct_behavior_subopt&not(correct_s_decoded)&not(consistent)&selected)/nnz(selected);

%% BP decomposition terms

est_out = estimate_B_perf_decomposition_terms( s(selected), choice(selected), choice_sh(selected), sp(selected), sp_sh(selected), consistent(selected), consistent_shuffle(selected), 'correct', false);

out.p_reg_s1 = est_out.p_reg_s1;
out.p_reg_s2 = est_out.p_reg_s2;
out.p_reg_all = est_out.p_reg_all;

out.p_c_reg_s1 = est_out.p_c_reg_s1;
out.p_c_reg_s2 = est_out.p_c_reg_s2;
out.p_c_reg_all = est_out.p_c_reg_all;

out.p_creg_s1 = est_out.p_creg_s1;
out.p_creg_s2 = est_out.p_creg_s2;
out.p_creg_all = est_out.p_creg_all;

out.est_B_perf_s1 = est_out.B_perf_s1;
out.est_B_perf_s2 = est_out.B_perf_s2;
out.est_B_perf_all = est_out.B_perf_all;

% shuffled values
out.p_reg_s1_shuffle = est_out.p_reg_s1_sh;
out.p_reg_s2_shuffle = est_out.p_reg_s2_sh;
out.p_reg_all_shuffle = est_out.p_reg_all_sh;

out.p_c_reg_s1_shuffle = est_out.p_c_reg_s1_sh;
out.p_c_reg_s2_shuffle = est_out.p_c_reg_s2_sh;
out.p_c_reg_all_shuffle = est_out.p_c_reg_all_sh;

out.p_creg_s1_shuffle = est_out.p_creg_s1_sh;
out.p_creg_s2_shuffle = est_out.p_creg_s2_sh;
out.p_creg_all_shuffle = est_out.p_creg_all_sh;

out.est_B_perf_s1_shuffle = est_out.B_perf_s1_sh;
out.est_B_perf_s2_shuffle = est_out.B_perf_s2_sh;
out.est_B_perf_all_shuffle = est_out.B_perf_all_sh;

est_out_subopt = estimate_B_perf_decomposition_terms( s(selected), choice_subopt(selected), choice_subopt_sh(selected), sp(selected), sp_sh(selected), consistent(selected), consistent_shuffle(selected), 'correct', false);

out.p_reg_s1_subopt = est_out_subopt.p_reg_s1;
out.p_reg_s2_subopt = est_out_subopt.p_reg_s2;
out.p_reg_all_subopt = est_out_subopt.p_reg_all;

out.p_c_reg_s1_subopt = est_out_subopt.p_c_reg_s1;
out.p_c_reg_s2_subopt = est_out_subopt.p_c_reg_s2;
out.p_c_reg_all_subopt = est_out_subopt.p_c_reg_all;

out.p_creg_s1_subopt = est_out_subopt.p_creg_s1;
out.p_creg_s2_subopt = est_out_subopt.p_creg_s2;
out.p_creg_all_subopt = est_out_subopt.p_creg_all;

out.est_B_perf_s1_subopt = est_out_subopt.B_perf_s1;
out.est_B_perf_s2_subopt = est_out_subopt.B_perf_s2;
out.est_B_perf_all_subopt = est_out_subopt.B_perf_all;

%

out.p_reg_s1_subopt_sh = est_out_subopt.p_reg_s1_sh;
out.p_reg_s2_subopt_sh = est_out_subopt.p_reg_s2_sh;
out.p_reg_all_subopt_sh = est_out_subopt.p_reg_all_sh;

out.p_c_reg_s1_subopt_sh = est_out_subopt.p_c_reg_s1_sh;
out.p_c_reg_s2_subopt_sh = est_out_subopt.p_c_reg_s2_sh;
out.p_c_reg_all_subopt_sh = est_out_subopt.p_c_reg_all_sh;

out.p_creg_s1_subopt_sh = est_out_subopt.p_creg_s1_sh;
out.p_creg_s2_subopt_sh = est_out_subopt.p_creg_s2_sh;
out.p_creg_all_subopt_sh = est_out_subopt.p_creg_all_sh;

out.est_B_perf_s1_subopt_sh = est_out_subopt.B_perf_s1_sh;
out.est_B_perf_s2_subopt_sh = est_out_subopt.B_perf_s2_sh;
out.est_B_perf_all_subopt_sh = est_out_subopt.B_perf_all_sh;

if display_figures
    
    b2_x=[-5 5];
    b2_y=m2*b2_x;
    
    if m1 == -inf
        b1_x = [0 0];
        b1_y = [-5 5];
    else
        b1_x=[-5 5];
        b1_y=m1*b1_x;
    end
    
    %% Visualize trials in response space
    
    figure(11); clf;
    set(gcf,'Position',[200 200 900 400])
    
    subplot(1,2,1);
    scatter(r(s==sp & s==0 & selected,1),r(s==sp & s==0 & selected,2),'g','filled'); hold on;
    scatter(r(s~=sp & s==0 & selected,1),r(s~=sp & s==0 & selected,2),'g');
    scatter(r(s==sp & s==1 & selected,1),r(s==sp & s==1 & selected,2),'b','filled');
    scatter(r(s~=sp & s==1 & selected,1),r(s~=sp & s==1 & selected,2),'b');
    fimplicit(f_db,'Color','k','LineStyle','--','Linewidth',2);
    xlim([-1.5 2.5]);
    ylim([-1.5 2.5]);
    xlabel('r1 [a.u.]');
    ylabel('r2 [a.u.]');
    title('Correlated responses');
    axis square
    legend('s=0 sp=0','s=0 sp=1','s=1 sp=1','s=1 sp=0','Location','southeast');
    
    subplot(1,2,2);
    scatter(r_sh(s==sp_sh & s==0 & selected_shuffle,1),r_sh(s==sp_sh & s==0 & selected_shuffle,2),'g','filled'); hold on;
    scatter(r_sh(s~=sp_sh & s==0 & selected_shuffle,1),r_sh(s~=sp_sh & s==0 & selected_shuffle,2),'g');
    scatter(r_sh(s==sp_sh & s==1 & selected_shuffle,1),r_sh(s==sp_sh & s==1 & selected_shuffle,2),'b','filled');
    scatter(r_sh(s~=sp_sh & s==1 & selected_shuffle,1),r_sh(s~=sp_sh & s==1 & selected_shuffle,2),'b');
    fimplicit(f_db_sh,'Color','k','LineStyle','--','Linewidth',2);
    xlim([-1.5 2.5]);
    ylim([-1.5 2.5]);
    xlabel('r1 [a.u.]');
    ylabel('r2 [a.u.]');
    title('Shuffled responses');
    axis square
    
    %% Visualize trials in posterior space
    
    figure(12); clf;
    set(gcf,'Position',[200 200 900 400])
    
    subplot(1,2,1);
    scatter(ps_1(s==sp & s==0 & selected),ps_2(s==sp & s==0 & selected),'g','filled'); hold on;
    scatter(ps_1(s~=sp & s==0 & selected),ps_2(s~=sp & s==0 & selected),'g');
    scatter(ps_1(s==sp & s==1 & selected),ps_2(s==sp & s==1 & selected),'b','filled');
    scatter(ps_1(s~=sp & s==1 & selected),ps_2(s~=sp & s==1 & selected),'b');
    line(b1_x+.5,b1_y+.5,'Color','red','LineStyle','--','Linewidth',2);
    line(b2_x+.5,b2_y+.5,'Color','red','LineStyle','--','Linewidth',2);
    xlim([0 1]);
    ylim([0 1]);
    xlabel('p(s|r1)');
    ylabel('p(s|r2)');
    title('Correlated responses');
    axis square
    legend('s=0 sp=0','s=0 sp=1','s=1 sp=1','s=1 sp=0','Location','southeast');
    
    subplot(1,2,2);
    scatter(ps_1_sh(s==sp_sh & s==0 & selected_shuffle),ps_2_sh(s==sp_sh & s==0 & selected_shuffle),'g','filled');  hold on;
    scatter(ps_1_sh(s~=sp_sh & s==0 & selected_shuffle),ps_2_sh(s~=sp_sh & s==0 & selected_shuffle),'g');
    scatter(ps_1_sh(s==sp_sh & s==1 & selected_shuffle),ps_2_sh(s==sp_sh & s==1 & selected_shuffle),'b','filled');
    scatter(ps_1_sh(s~=sp_sh & s==1 & selected_shuffle),ps_2_sh(s~=sp_sh & s==1 & selected_shuffle),'b');
    line(b1_x+.5,b1_y+.5,'Color','red','LineStyle','--','Linewidth',2);
    line(b2_x+.5,b2_y+.5,'Color','red','LineStyle','--','Linewidth',2);
    xlim([0 1]);
    ylim([0 1]);
    xlabel('p(s|r1)');
    ylabel('p(s|r2)');
    title('Shuffled responses');
    axis square
    
    %% Visualize trials by consistency
    figure(13); clf;
    set(gcf,'Position',[200 200 900 400])
    
    subplot(1,2,1);
    hold off
    scatter(ps_1(consistent & selected),ps_2(consistent & selected),'ro','MarkerEdgeColor',[0 0 0]); hold on
    scatter(ps_1(inconsistent & selected),ps_2(inconsistent & selected),'ko','MarkerEdgeColor',[.5 .5 .5])
    line(b1_x+.5,b1_y+.5,'Color','red','LineStyle','--','Linewidth',2);
    line(b2_x+.5,b2_y+.5,'Color','red','LineStyle','--','Linewidth',2);
    xlim([0 1])
    ylim([0 1])
    legend({'consistent','inconsistent'},'Location','southeast');
    title(sprintf('Correlated reasponses'));
    xlabel('p(s|r1)');
    ylabel('p(s|r2)');
    axis square
    
    
    subplot(1,2,2);
    hold off
    scatter(ps_1_sh(consistent_shuffle & selected_shuffle),ps_2_sh(consistent_shuffle & selected_shuffle),'ro','MarkerEdgeColor',[0 0 0]); hold on
    scatter(ps_1_sh(inconsistent_shuffle & selected_shuffle),ps_2_sh(inconsistent_shuffle & selected_shuffle),'ko','MarkerEdgeColor',[.5 .5 .5])
    line(b1_x+.5,b1_y+.5,'Color','red','LineStyle','--','Linewidth',2);
    line(b2_x+.5,b2_y+.5,'Color','red','LineStyle','--','Linewidth',2);
    xlim([0 1])
    ylim([0 1])
    legend({'consistent','inconsistent'},'FontSize',8);
    title(sprintf('Shuffled responses'));
    xlabel('p(s|r1)');
    ylabel('p(s|r2)');
    axis square
    
end


end