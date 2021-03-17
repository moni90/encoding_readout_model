function betas = determine_betas_from_readout_strength_params_model2(readout_strength_params)

% logistic model of reference: p(c=1) = 

% readout_strength_params: p(c=sp|region)
% regions:
% 1: sp=0, con=1
% 2: sp=1, con=1
% 3: sp=0, con=0
% 4: sp=1, con=0

% convert from p(c=sp|region) to p(c=1|region)
readout_strength_params(1) = 1-readout_strength_params(1);
readout_strength_params(3) = 1-readout_strength_params(3);

prob_logit = log(readout_strength_params./(1-readout_strength_params));
% prob_logit = logit(readout_strength_params);

betas(1) = (prob_logit(4) + prob_logit(3))/2;
betas(2) = (prob_logit(4) - prob_logit(3))/2;
betas(3) = prob_logit(2) - prob_logit(4);
betas(4) = prob_logit(3) - prob_logit(1); % check this

end