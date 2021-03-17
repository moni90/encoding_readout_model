function [pair_corr_coeff_BP_corr, pair_corr_coeff_BP_err, pop_corr_coeff_BP_corr, pop_corr_coeff_BP_err, n_all_selected_trials] = compute_correlation_coefficient_correct_error_N(stimuli,choices,r1,r2)

% computes correlation coefficient separately for correct behavior and
% incorrect behavior trials

% performs trials subsampling to ensure that a different number of trials
% does not influence the comparison

% computes correlation coefficient separately for the two stimuli, then
% averages across the two

all_stimuli = unique(stimuli);

for s_idx = 1:numel(all_stimuli)
    
    this_stimulus = all_stimuli(s_idx);
    
    trials_this_stimulus = stimuli==this_stimulus;
    
    correct_trials_this_stimulus = stimuli==this_stimulus & choices==this_stimulus;
    error_trials_this_stimulus = stimuli==this_stimulus & not(choices==this_stimulus);
    
    n_error = nnz(error_trials_this_stimulus); % nr of error trials
    n_correct = nnz(correct_trials_this_stimulus); % nr of correct trials
%     n_trials_this_stimulus = nnz(trials_this_stimulus);
    n_trials = numel(stimuli);
    
    selected_correct = zeros(n_trials,1);
    selected_error = zeros(n_trials,1);
    
    if min([n_error,n_correct])<=1    
        corr_coeff_BP_corr_tmp(s_idx) = nan;
        corr_coeff_BP_err_tmp(s_idx) = nan;     
    else
        if n_error<=n_correct
            idx_selected_correct = randsample(find(correct_trials_this_stimulus),n_error);
            selected_correct(idx_selected_correct)=1;
            selected_correct=logical(selected_correct);
            selected_error=logical(error_trials_this_stimulus);
        else
            idx_selected_error = randsample(find(error_trials_this_stimulus),n_correct);
            selected_error(idx_selected_error)=1;
            selected_error=logical(selected_error);
            selected_correct=logical(correct_trials_this_stimulus);
        end
%         cmat_corr = corrcoef(r1(selected_correct),r2(selected_correct));
%         corr_coeff_BP_corr_tmp(s_idx) = cmat_corr(1,2);
%         cmat_err = corrcoef(r1(selected_error),r2(selected_error));
%         corr_coeff_BP_err_tmp(s_idx) = cmat_err(1,2);
        
        [~, ~, ~, ~, var_explained_correct] = pca([r1(selected_correct,:) r2(selected_correct,:)], 'NumComponents', 1);
        pop_corr_coeff_BP_corr_tmp(s_idx) = var_explained_correct(1);
        [~, ~, ~, ~, var_explained_error] = pca([r1(selected_error,:) r2(selected_error,:)], 'NumComponents', 1);
        pop_corr_coeff_BP_err_tmp(s_idx) = var_explained_error(1);
        
        pair_corr_coeff_BP_corr_tmp(s_idx) = corr_coeff_pair(r1(selected_correct,:), r2(selected_correct,:));
        pair_corr_coeff_BP_err_tmp(s_idx) = corr_coeff_pair(r1(selected_error,:), r2(selected_error,:));

    end
    
    trial_weight(s_idx) = nnz(selected_correct); % would be the same as nnz(selected_error) because trials are subsampled. THis weight is useful to opportunely combine the correlation coefficients for the two stimuli
    
    % figure(1); clf;
    % scatter(r1(selected_correct),r2(selected_correct),'b');
    % hold on;
    % scatter(r1(selected_error),r2(selected_error),'r');
    % text(.5,.8,sprintf('corr coeff BP corr: %.2f BP err %.2f',corr_coeff_BP_corr_tmp,corr_coeff_BP_err_tmp))
    
    
end

n_all_selected_trials = sum(trial_weight);

pair_corr_coeff_BP_corr = 0;
pair_corr_coeff_BP_err = 0;
pop_corr_coeff_BP_corr = 0;
pop_corr_coeff_BP_err = 0;
% sum over the two stimuli with appropriate weights

for s_idx = 1:numel(all_stimuli)
    
    pair_corr_coeff_BP_corr = pair_corr_coeff_BP_corr + trial_weight(s_idx)./n_all_selected_trials * pair_corr_coeff_BP_corr_tmp(s_idx);
    pair_corr_coeff_BP_err = pair_corr_coeff_BP_err + trial_weight(s_idx)./n_all_selected_trials * pair_corr_coeff_BP_err_tmp(s_idx);
    pop_corr_coeff_BP_corr = pop_corr_coeff_BP_corr + trial_weight(s_idx)./n_all_selected_trials * pop_corr_coeff_BP_corr_tmp(s_idx);
    pop_corr_coeff_BP_err = pop_corr_coeff_BP_err + trial_weight(s_idx)./n_all_selected_trials * pop_corr_coeff_BP_err_tmp(s_idx);
    
end





end

function this_pairCorr_avg = corr_coeff_pair(r1, r2)

if size(r1,2)==1
    [this_pairCorr, pval_corr] = corr(r1, r2);
else
    this_pairCorr=[];
    for c1=1:size(r1,2)-1
        for c2=c1+1:size(r1,2)
            [this_pairCorr_temp, pval_corr_temp] = corr(r1(:,c1), r2(:,c2));
            this_pairCorr = [this_pairCorr; this_pairCorr_temp];
        end
    end
end
this_pairCorr_avg = nanmean(this_pairCorr);

end
