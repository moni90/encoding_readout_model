function [sp,p_s, train_idxs, K, L, mu1, mu2, sigma ]=lda_classify(r, s, fitting_fraction, train_idxs)

% L is the normal to the decoding boundary. It is oriented towards Class1.

n_trials = numel(s);

if nargin < 4
    
    train_idxs = get_training_trials_indices(s,fitting_fraction);
    
end

test_idxs = setdiff(1:n_trials,train_idxs)';

Mdl=fitcdiscr(r(train_idxs,:),s(train_idxs,:),'DiscrimType','linear');
% 
% K = -.5*(mu1+mu2)'*sigma_inv*(mu1-mu2);
% L = sigma_inv * (mu1-mu2);

K=Mdl.Coeffs(2,1).Const;
L=Mdl.Coeffs(2,1).Linear;
% L is the normal to the decoding boundary. It is oriented towards Class1.

mu1 = Mdl.Mu(1,:);
mu2 = Mdl.Mu(2,:);
sigma = Mdl.Sigma;

[sp,p_s]=predict(Mdl,r(test_idxs,:));