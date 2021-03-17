function out=estimate_B_perf_decomposition_terms(real_stimuli, real_choices, shuffled_choices, decoded_stimuli, decoded_stimuli_sh, consistency, consistency_sh, which_sp, compute_sp_choice_match)

%% real data
switch which_sp   
    case 'correct'
        regions = 2*(1-decoded_stimuli==real_stimuli)+(1-consistency)+1;
        regions_sh = 2*(1-decoded_stimuli_sh==real_stimuli)+(1-consistency_sh)+1;
        % 1: corr decoded, consistent
        % 2: corr decoded, inconsistent
        % 3: incorrect decoded, consistent
        % 4: incorrect decoded, inconsistent        
    case 'left'
        regions = 2*(1-decoded_stimuli==1)+(1-consistency)+1;
        regions_sh = 2*(1-decoded_stimuli_sh==1)+(1-consistency_sh)+1;
        % 1: left, consistent
        % 2: left, inconsistent
        % 3: right, consistent
        % 4: right, inconsistent  
end

n_trials = numel(real_stimuli);
p_s1 = nnz(real_stimuli==0)/n_trials;
p_s2 = nnz(real_stimuli==1)/n_trials;

for this_region=1:4
    
    this_region_trials_s1 = regions == this_region & real_stimuli==0;
    this_region_trials_s2 = regions == this_region & real_stimuli==1;
    
    p_reg_s1(this_region) = nnz(this_region_trials_s1)/nnz(real_stimuli==0);
    p_reg_s2(this_region) = nnz(this_region_trials_s2)/nnz(real_stimuli==1);
    p_reg_all(this_region) = p_s1*p_reg_s1(this_region)+p_s2*p_reg_s2(this_region);
    
    if compute_sp_choice_match
        p_c_reg_s1(this_region) = nnz(real_choices==decoded_stimuli & this_region_trials_s1)/nnz(this_region_trials_s1);
        p_c_reg_s2(this_region) = nnz(real_choices==decoded_stimuli & this_region_trials_s2)/nnz(this_region_trials_s2);
    else
        p_c_reg_s1(this_region) = nnz(real_choices==real_stimuli & this_region_trials_s1)/nnz(this_region_trials_s1);
        p_c_reg_s2(this_region) = nnz(real_choices==real_stimuli & this_region_trials_s2)/nnz(this_region_trials_s2);
    end
    p_c_reg_all(this_region) = p_s1*p_c_reg_s1(this_region)+p_s2*p_c_reg_s2(this_region);
    
    p_creg_s1(this_region) = p_reg_s1(this_region)*p_c_reg_s1(this_region);
    p_creg_s2(this_region) = p_reg_s2(this_region)*p_c_reg_s2(this_region);
    p_creg_all(this_region) = p_s1*p_creg_s1(this_region)+p_s2*p_creg_s2(this_region);
    
end

B_perf_s1 = sum(p_creg_s1);
B_perf_s2 = sum(p_creg_s2);
B_perf_all = sum(p_creg_all);

%% shuffled data



% n_trials = numel(real_stimuli);
% p_s1 = nnz(real_stimuli==0)/n_trials;
% p_s2 = nnz(real_stimuli==1)/n_trials;

for this_region=1:4
    
    this_region_trials_s1 = regions_sh == this_region & real_stimuli==0;
    this_region_trials_s2 = regions_sh == this_region & real_stimuli==1;
    
    p_reg_s1_sh(this_region) = nnz(this_region_trials_s1)/nnz(real_stimuli==0);
    p_reg_s2_sh(this_region) = nnz(this_region_trials_s2)/nnz(real_stimuli==1);
    p_reg_all_sh(this_region) = p_s1*p_reg_s1_sh(this_region)+p_s2*p_reg_s2_sh(this_region);
    
    if compute_sp_choice_match
        p_c_reg_s1_sh(this_region) = nnz(shuffled_choices==decoded_stimuli_sh & this_region_trials_s1)/nnz(this_region_trials_s1);
        p_c_reg_s2_sh(this_region) = nnz(shuffled_choices==decoded_stimuli_sh & this_region_trials_s2)/nnz(this_region_trials_s2);
    else
        p_c_reg_s1_sh(this_region) = nnz(shuffled_choices==real_stimuli & this_region_trials_s1)/nnz(this_region_trials_s1);
        p_c_reg_s2_sh(this_region) = nnz(shuffled_choices==real_stimuli & this_region_trials_s2)/nnz(this_region_trials_s2);
    end
    p_c_reg_all_sh(this_region) = p_s1*p_c_reg_s1_sh(this_region)+p_s2*p_c_reg_s2_sh(this_region);
    
    p_creg_s1_sh(this_region) = p_reg_s1_sh(this_region)*p_c_reg_s1_sh(this_region);
    p_creg_s2_sh(this_region) = p_reg_s2_sh(this_region)*p_c_reg_s2_sh(this_region);
    p_creg_all_sh(this_region) = p_s1*p_creg_s1_sh(this_region)+p_s2*p_creg_s2_sh(this_region);
    
end

B_perf_s1_sh = sum(p_creg_s1_sh);
B_perf_s2_sh = sum(p_creg_s2_sh);
B_perf_all_sh = sum(p_creg_all_sh);

s_labels = {'s1','s2','all'};

% Correct for nans
% 
for this_region=1:4

    [p_reg_s1(this_region),p_reg_s1_sh(this_region)]=set_to_nan(p_reg_s1(this_region),p_reg_s1_sh(this_region));
    [p_reg_s2(this_region),p_reg_s2_sh(this_region)]=set_to_nan(p_reg_s2(this_region),p_reg_s2_sh(this_region));
    [p_reg_all(this_region),p_reg_all_sh(this_region)]=set_to_nan(p_reg_all(this_region),p_reg_all_sh(this_region));
    
    [p_c_reg_s1(this_region),p_c_reg_s1_sh(this_region)]=set_to_nan(p_c_reg_s1(this_region),p_c_reg_s1_sh(this_region));
    [p_c_reg_s2(this_region),p_c_reg_s2_sh(this_region)]=set_to_nan(p_c_reg_s2(this_region),p_c_reg_s2_sh(this_region));
    [p_c_reg_all(this_region),p_c_reg_all_sh(this_region)]=set_to_nan(p_c_reg_all(this_region),p_c_reg_all_sh(this_region));
    
    [p_creg_s1(this_region),p_creg_s1_sh(this_region)]=set_to_nan(p_creg_s1(this_region),p_creg_s1_sh(this_region));
    [p_creg_s2(this_region),p_creg_s2_sh(this_region)]=set_to_nan(p_creg_s2(this_region),p_creg_s2_sh(this_region));
    [p_creg_all(this_region),p_creg_all_sh(this_region)]=set_to_nan(p_creg_all(this_region),p_creg_all_sh(this_region));

end

[B_perf_s1,B_perf_s1_sh]=set_to_nan(B_perf_s1,B_perf_s1_sh);
[B_perf_s2,B_perf_s2_sh]=set_to_nan(B_perf_s2,B_perf_s2_sh);
[B_perf_all,B_perf_s1_all]=set_to_nan(B_perf_all,B_perf_all_sh);


% Save output

for i=1:numel(s_labels)
    
    this_field=sprintf('p_reg_%s',s_labels{i});  
out.(this_field)=eval(this_field);
    this_field=sprintf('p_c_reg_%s',s_labels{i});  
out.(this_field)=eval(this_field);
    this_field=sprintf('p_creg_%s',s_labels{i});  
out.(this_field)=eval(this_field);
    this_field=sprintf('B_perf_%s',s_labels{i});  
out.(this_field)=eval(this_field);

    this_field=sprintf('p_reg_%s_sh',s_labels{i});  
out.(this_field)=eval(this_field);
    this_field=sprintf('p_c_reg_%s_sh',s_labels{i});  
out.(this_field)=eval(this_field);
    this_field=sprintf('p_creg_%s_sh',s_labels{i});  
out.(this_field)=eval(this_field);
    this_field=sprintf('B_perf_%s_sh',s_labels{i});  
out.(this_field)=eval(this_field);

end





