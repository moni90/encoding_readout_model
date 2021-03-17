function run_simu_suppl_model_figures_N_v7(N, n_simu,n_samples,decoder_type, dist_stim, stdev, gamma_range, rho_range, consistency_thr, equalization_type, id_string,n_workers)

% logfile = './log_file.txt';
% initialize_log(logfile);

%rho = .9;
%stdev = .3;
%gamma_range = [0,pi/20,pi/10];
% n_workers = 20;

readout_strength = 'custom';
% rho_within = 0.5;
% gamma_range = [0:.01:.25]; % in pi units
consistency_strength_range = [0:.025:.25]; % 0:ideal consistency readout, .25: Bayesian readout
% rho_range = [.1:.1:.9];

% n_samples = 500;
stdev1=stdev;
stdev2=stdev;
% decoder_type='analytical';
fitting_fraction=.5;
% readout_strength='ideal';
% readout_strength_params=[];
alpha_readout = pi/2;

behav_perf = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
behav_perf_con = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
behav_perf_inc = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
behav_perf_sh = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
behav_perf_subopt = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
behav_perf_subopt_sh = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
S_dec_perf = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
S_dec_perf_con = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
S_dec_perf_inc = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
S_dec_perf_sh = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
frac_cons = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
frac_cons_sh = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
frac_cons_BP_correct = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
frac_cons_BP_correct_subopt = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
frac_cons_BP_error = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
frac_cons_BP_error_subopt = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
p_reg_all = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
p_reg_all_sh = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
p_c_reg_all = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
p_c_reg_all_sh = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
p_creg_all = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
p_creg_all_sh = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
p_c_reg_all_subopt = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
p_creg_all_subopt = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
theta = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
theta_sh = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
gamma_out = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
gamma_out_sh = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
avg_readout_strength = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
corr_coeff_BP_correct = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
corr_coeff_BP_error = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
corr_coeff_BP_correct_subopt = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));
corr_coeff_BP_error_subopt = cell(numel(consistency_strength_range),size(gamma_range,1),numel(rho_range));

success = 0;
trial = 1;
while success == 0
    try
        if isempty(gcp('nocreate'))
            parpool(n_workers);
        end
        success = 1;
    catch
        disp(['trial ' num2str(trial) ' failed. Try again parallel initialization']);
        trial = trial+1;
    end
end
    
for k = 1:numel(rho_range)
    
    rho = rho_range(k);
    
    for i = 1:numel(consistency_strength_range)
        
        consistency_strength = consistency_strength_range(i);
        readout_strength_params = [1-consistency_strength 1-consistency_strength .5+consistency_strength .5+consistency_strength];
        
        parfor j = 1:numel(gamma_range)
            
            gamma = pi*gamma_range(j);
            
            fprintf('Iteration %d/%d\n',...
                (k-1)*numel(consistency_strength_range)*numel(gamma_range)+(i-1)*numel(gamma_range)+j,...
                numel(rho_range)*numel(consistency_strength_range)*numel(gamma_range));
            log_string = ['Iteration ' num2str((k-1)*numel(consistency_strength_range)*numel(gamma_range)+(i-1)*numel(gamma_range)+j) ...
                '/' num2str(numel(rho_range)*numel(consistency_strength_range)*numel(gamma_range)) '\n'];
%             write_log(logfile,log_string);
            
            
            this_rho = rho;

            this_behav_perf = zeros(n_simu,1);
            this_behav_perf_sh = zeros(n_simu,1);
            
            this_behav_perf_subopt = zeros(n_simu,1);
            this_behav_perf_subopt_sh = zeros(n_simu,1);
            
            this_S_dec_perf = zeros(n_simu,1);
            this_S_dec_perf_sh = zeros(n_simu,1);
            
            this_frac_cons = zeros(n_simu,1);
            this_frac_cons_sh = zeros(n_simu,1);
            
            this_frac_cons_BP_correct = zeros(n_simu,1);
            this_frac_cons_BP_correct_subopt = zeros(n_simu,1);
            this_frac_cons_BP_error = zeros(n_simu,1);
            this_frac_cons_BP_error_subopt = zeros(n_simu,1);
            
            this_p_reg_all = zeros(n_simu,4);
            this_p_reg_all_sh = zeros(n_simu,4);
            this_p_c_reg_all = zeros(n_simu,4);
            this_p_c_reg_all_subopt = zeros(n_simu,4);
            this_p_c_reg_all_sh = zeros(n_simu,4);
            this_p_creg_all = zeros(n_simu,4);
            this_p_creg_all_subopt = zeros(n_simu,4);
            this_p_creg_all_sh = zeros(n_simu,4);
            
            this_theta = zeros(n_simu,1);
            this_theta_sh = zeros(n_simu,1);
            
            this_gamma = zeros(n_simu,1);
            this_gamma_2d = zeros(n_simu,1);
            this_gamma_sh = zeros(n_simu,1);
            
            this_avg_readout_strength = zeros(n_simu,1);
            
            this_corr_coeff_BP_correct = zeros(n_simu,1);
            this_corr_coeff_BP_error = zeros(n_simu,1);
            
            this_corr_coeff_BP_correct_subopt = zeros(n_simu,1);
            this_corr_coeff_BP_error_subopt = zeros(n_simu,1);

            for ns=1:n_simu
                
                sign_dir = generate_random_dir(ones(N,1), gamma, 1);
                sign_dir = sign_dir/norm(sign_dir);
                mu2 = dist_stim*sign_dir' + .5;
                mu1 = -dist_stim*sign_dir' + .5;
                
                mu_s0 = mu1;
                mu_s1 = mu2;

                stdev1_s0 = stdev1;
                stdev2_s0 = stdev2;
                stdev1_s1 = stdev1;
                stdev2_s1 = stdev2;

                rho_s0 = this_rho;
                rho_s1 = this_rho;
                
                %                 out=simu_consistency_readout(n_samples,mu1,mu2,stdev1,stdev2,this_rho,decoder_type,fitting_fraction,...
                %                     alpha_readout,readout_strength,readout_strength_params,equalization_type,false);
                
                rho_within = this_rho; %same RHO for all covariance matrix
                
                out=simu_consistency_readout_N_v7(N,n_samples,mu_s0,mu_s1,stdev1_s0,stdev2_s0,stdev1_s1,stdev2_s1,rho_s0,rho_s1,rho_within,decoder_type,fitting_fraction,...
                    alpha_readout,readout_strength,readout_strength_params,consistency_thr,equalization_type,false);
                
                this_behav_perf(ns)=out.fraction_correct_behavior;
                %                 behav_perf_con(ns)=out.fraction_correct_behavior_con;
                %                 behav_perf_inc(ns)=out.fraction_correct_behavior_inc;
                this_behav_perf_sh(ns)=out.fraction_correct_behavior_shuffle;
                
                this_behav_perf_subopt(ns)=out.fraction_correct_behavior_subopt;
                this_behav_perf_subopt_sh(ns)=out.fraction_correct_behavior_subopt_shuffle;
                
                this_S_dec_perf(ns)=out.fraction_correct_S_decoded;
                %                 S_dec_perf_con(ns)=out.fraction_correct_S_decoded_con;
                %                 S_dec_perf_inc(ns)=out.fraction_correct_S_decoded_inc;
                this_S_dec_perf_sh(ns)=out.fraction_correct_S_decoded_shuffle;
                
                this_frac_cons(ns)=out.fraction_consistent;
                this_frac_cons_sh(ns)=out.fraction_consistent_shuffle;
                
                this_frac_cons_BP_correct(ns)=out.fraction_consistent_BP_correct;
                this_frac_cons_BP_correct_subopt(ns)=out.fraction_consistent_BP_correct_subopt;
                this_frac_cons_BP_error(ns)=out.fraction_consistent_BP_error;
                this_frac_cons_BP_error_subopt(ns)=out.fraction_consistent_BP_error_subopt;
                
                this_p_reg_all(ns,:)=out.p_reg_all;
                this_p_reg_all_sh(ns,:)=out.p_reg_all_shuffle;
                this_p_c_reg_all(ns,:)=out.p_c_reg_all;
                this_p_c_reg_all_subopt(ns,:)=out.p_c_reg_all_subopt;
                this_p_c_reg_all_sh(ns,:)=out.p_c_reg_all_shuffle;
                this_p_creg_all(ns,:)=out.p_creg_all;
                this_p_creg_all_subopt(ns,:)=out.p_creg_all_subopt;
                this_p_creg_all_sh(ns,:)=out.p_creg_all_shuffle;
                
                this_theta(ns)=out.theta;
                this_theta_sh(ns)=out.theta_sh;
                
                this_gamma(ns)=out.gamma;
                this_gamma_2d(ns)=out.gamma_2d;
                this_gamma_sh(ns)=out.gamma_sh;
                
                this_avg_readout_strength(ns)=out.avg_readout_strength;
                
                this_corr_coeff_BP_correct(ns) = out.corr_coeff_BP_correct;
                this_corr_coeff_BP_error(ns) = out.corr_coeff_BP_error;
                
                this_corr_coeff_BP_correct_subopt(ns) = out.corr_coeff_BP_correct_subopt;
                this_corr_coeff_BP_error_subopt(ns) = out.corr_coeff_BP_error_subopt;
                
            end
            
            
            behav_perf{i,j,k}=this_behav_perf;
            %                 behav_perf_con{i,j,k}(ns)=out.fraction_correct_behavior_con;
            %                 behav_perf_inc{i,j,k}(ns)=out.fraction_correct_behavior_inc;
            behav_perf_sh{i,j,k}=this_behav_perf_sh;
            
            behav_perf_subopt{i,j,k}=this_behav_perf_subopt;
            behav_perf_subopt_sh{i,j,k}=this_behav_perf_subopt_sh;
            
            S_dec_perf{i,j,k}=this_S_dec_perf;
            %                 S_dec_perf_con{i,j,k}(ns)=out.fraction_correct_S_decoded_con;
            %                 S_dec_perf_inc{i,j,k}(ns)=out.fraction_correct_S_decoded_inc;
            S_dec_perf_sh{i,j,k}=this_S_dec_perf_sh;
            
            frac_cons{i,j,k}=this_frac_cons;
            frac_cons_sh{i,j,k}=this_frac_cons_sh;
            
            frac_cons_BP_correct{i,j,k}=this_frac_cons_BP_correct;
            frac_cons_BP_correct_subopt{i,j,k}=this_frac_cons_BP_correct_subopt;
            frac_cons_BP_error{i,j,k}=this_frac_cons_BP_error;
            frac_cons_BP_error_subopt{i,j,k}=this_frac_cons_BP_error_subopt;
            
            p_reg_all{i,j,k}=this_p_reg_all;
            p_reg_all_sh{i,j,k}=this_p_reg_all_sh;
            p_c_reg_all{i,j,k}=this_p_c_reg_all;
            p_c_reg_all_subopt{i,j,k}=this_p_c_reg_all_subopt;
            p_c_reg_all_sh{i,j,k}=this_p_c_reg_all_sh;
            p_creg_all{i,j,k}=this_p_creg_all;
            p_creg_all_subopt{i,j,k}=this_p_creg_all_subopt;
            p_creg_all_sh{i,j,k}=this_p_creg_all_sh;
            
            theta{i,j,k}=this_theta;
            theta_sh{i,j,k}=this_theta_sh;
            
            gamma_out{i,j,k}=this_gamma;
            gamma_2d_out{i,j,k}=this_gamma_2d;
            gamma_out_sh{i,j,k}=this_gamma_sh;
            
            avg_readout_strength{i,j,k}=this_avg_readout_strength;
            
            corr_coeff_BP_correct{i,j,k}=this_corr_coeff_BP_correct;
            corr_coeff_BP_error{i,j,k}=this_corr_coeff_BP_error;
            
            corr_coeff_BP_correct_subopt{i,j,k}=this_corr_coeff_BP_correct_subopt;
            corr_coeff_BP_error_subopt{i,j,k}=this_corr_coeff_BP_error_subopt;
            
        end
    end
   
end

delete(gcp)

simu_out.behav_perf = behav_perf;
simu_out.behav_perf_con = behav_perf_con;
simu_out.behav_perf_inc = behav_perf_inc;
simu_out.behav_perf_sh = behav_perf_sh;

simu_out.behav_perf_subopt = behav_perf_subopt;
simu_out.behav_perf_subopt_sh = behav_perf_subopt_sh;

simu_out.S_dec_perf = S_dec_perf;
simu_out.S_dec_perf_con = S_dec_perf_con;
simu_out.S_dec_perf_inc = S_dec_perf_inc;
simu_out.S_dec_perf_sh = S_dec_perf_sh;

simu_out.frac_cons = frac_cons;
simu_out.frac_cons_sh = frac_cons_sh;

simu_out.frac_cons_BP_correct = frac_cons_BP_correct;
simu_out.frac_cons_BP_error = frac_cons_BP_error;
simu_out.frac_cons_BP_correct_subopt = frac_cons_BP_correct_subopt;
simu_out.frac_cons_BP_error_subopt = frac_cons_BP_error_subopt;

simu_out.p_reg_all = p_reg_all;
simu_out.p_reg_all_sh = p_reg_all_sh;
simu_out.p_c_reg_all = p_c_reg_all;
simu_out.p_c_reg_all_subopt = p_c_reg_all_subopt;
simu_out.p_c_reg_all_sh = p_c_reg_all_sh;
simu_out.p_creg_all = p_creg_all;
simu_out.p_creg_all_subopt = p_creg_all_subopt;
simu_out.p_creg_all_sh = p_creg_all_sh;

simu_out.theta = theta;
simu_out.theta_sh = theta_sh;
simu_out.gamma = gamma_out;
simu_out.gamma_2d = gamma_2d_out;
simu_out.gamma_sh = gamma_out_sh;

simu_out.avg_readout_strength = avg_readout_strength;

simu_out.corr_coeff_BP_correct = corr_coeff_BP_correct;
simu_out.corr_coeff_BP_error = corr_coeff_BP_error;
simu_out.corr_coeff_BP_correct_subopt = corr_coeff_BP_correct_subopt;
simu_out.corr_coeff_BP_error_subopt = corr_coeff_BP_error_subopt;
    
if not(isempty(equalization_type))
    equalization_type_string = ['_' equalization_type];
else
    equalization_type_string = '';
end

filename = sprintf('simulation_output_files/suppl_model_figure_simulation_output_%s_decoder_%dsimu_%dsamples%s%s.mat',decoder_type,n_simu,n_samples,equalization_type_string,id_string);

save([filename],'simu_out','consistency_strength_range','consistency_strength_range','gamma_range','rho_range','n_samples','n_simu','stdev',...
    'stdev1','stdev2','decoder_type','fitting_fraction','-v7.3');

end