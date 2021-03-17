function training_trials_idxs=get_training_trials_indices(s,fitting_fraction)

training_trials_idxs = [];
s_values = unique(s);
for this_s = 1:numel(s_values)
    this_s_trials = s==s_values(this_s);
    n_trials_this_s = nnz(this_s_trials);
    n_training_trials_this_s = round(n_trials_this_s*fitting_fraction);
    training_trials_idx_this_s = randsample(find(this_s_trials),n_training_trials_this_s);
    training_trials_idxs=[training_trials_idxs; training_trials_idx_this_s];
end

training_trials_idxs = sort(training_trials_idxs);
