function [w_n_s1,w_n_s2]=compute_noise_correlation_direction(r_s1,r_s2)

% r_s1 r_s2 are n_observations x n_features vectors of reponses to stimulus 1 and stimulus 2
% noise correlation direction is computed as the direction of maximum
% variability at fixed stimulus (first principal component)

%% PCA

pca_coeffs_s1=pca(r_s1);
pca_coeffs_s2=pca(r_s2);

w_n_s1 = pca_coeffs_s1(:,1);
w_n_s2 = pca_coeffs_s2(:,1);


end
