function [p_value] = myPermTest2_refined(sample1,sample2,group,perm_type,n_perms,tail)

% Performs a permutation test to check if the distributions sample1 and
% sample2 are statistically different. Uses the difference in mean value as
% test statistic. permutations are performed either among all trials
% (perm_type='all') or within paired subsets of the trials
% (perm_type='grouped') identified by the grouping variable group.

n_trials = numel(sample1);

% compute mean difference as test statistic
real_t = compute_mean_diff(sample1,sample2);
perm_t = nan(1,n_perms);

for i=1:n_perms
    
    switch perm_type
        
        case 'all'
            
            sample_all = [sample1; sample2];
            perm_indices = randperm(numel(sample_all));
            sample_all_perm = sample_all(perm_indices);
            sample1_perm = sample_all_perm(1:n_trials);
            sample2_perm = sample_all_perm(n_trials+1:end);
            
        case 'grouped'
            
            sample1_perm = nan(size(sample1));
            sample2_perm = nan(size(sample2));
            
            all_groups = unique(group);
            
            for g = 1:numel(all_groups)
                
                this_group = all_groups(g);
                this_group_trials = group==this_group;
                n_trials_this_group = nnz(this_group_trials);
                
                sample_all_this_group = [sample1(this_group_trials); sample2(this_group_trials)];
                perm_indices = randperm(numel(sample_all_this_group));
                sample_all_this_group_perm = sample_all_this_group(perm_indices);
                
                sample1_this_group_perm = sample_all_this_group_perm(1:n_trials_this_group);
                sample2_this_group_perm = sample_all_this_group_perm(n_trials_this_group+1:end);
                
                sample1_perm(this_group_trials) = sample1_this_group_perm;
                sample2_perm(this_group_trials) = sample2_this_group_perm;
                
            end
            
    end
    
    perm_t(i) = compute_mean_diff(sample1_perm,sample2_perm);
    
end

switch tail
    
    case 'two'
        p_value = nnz(abs(perm_t)>abs(real_t))/n_perms;
        % because our null hypothesis is that the permuted distribution is
        % centered around zero
        
    case 'right' % tests the hypothesis that sample1 mean is greater than sample2 mean (mean_diff > 0)
        p_value = nnz(perm_t>real_t)/n_perms;
         
    case 'left' % tests the hypothesis that sample2 mean in greater than sample1 mean (mean_diff < 0)
        p_value = nnz(perm_t<real_t)/n_perms;
        
end

% figure;
% histogram(perm_t,20);
% hold on;
% yl=ylim;
% plot([real_t real_t],[yl(1) yl(2)]);

end


function t = compute_mean_diff(sample1, sample2)

t = nanmean(sample1)-nanmean(sample2);

end
