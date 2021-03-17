function w_s=compute_signal_correlation_direction(r_s1,r_s2)

% r_s1 r_s2 are n_observations x n_features vectors of reponses to stimulus 1 and stimulus 2
% signal correlation direction is computed as the difference between the
% mean responses

w_s = mean(r_s2,1)-mean(r_s1,1);
w_s = w_s';

end