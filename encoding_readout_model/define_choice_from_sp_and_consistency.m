function choice = define_choice_from_sp_and_consistency(sp,con,readout_strength,readout_strength_params)

% choice is defined from predicted stimulus and consistency
% could be implemented using binornd instead odf randsample. output should
% be the same

choice=nan*ones(size(sp));

switch readout_strength
    case 'ideal'
        choice(sp==0 & con==0)=randsample([0,1],nnz(sp==0 & con==0),true,[.5 .5]);
        choice(sp==1 & con==0)=randsample([0,1],nnz(sp==1 & con==0),true,[.5 .5]);
        choice(sp==0 & con==1)=randsample([0,1],nnz(sp==0 & con==1),true,[1 0]);
        choice(sp==1 & con==1)=randsample([0,1],nnz(sp==1 & con==1),true,[0 1]);
    case 'strong'
        choice(sp==0 & con==0)=randsample([0,1],nnz(sp==0 & con==0),true,[.6 .4]);
        choice(sp==1 & con==0)=randsample([0,1],nnz(sp==1 & con==0),true,[.4 .6]);
        choice(sp==0 & con==1)=randsample([0,1],nnz(sp==0 & con==1),true,[.9 .1]);
        choice(sp==1 & con==1)=randsample([0,1],nnz(sp==1 & con==1),true,[.1 .9]);
    case 'custom'
        choice(sp==0 & con==0)=randsample([0,1],nnz(sp==0 & con==0),true,[readout_strength_params(3) 1-readout_strength_params(3)]);
        choice(sp==1 & con==0)=randsample([0,1],nnz(sp==1 & con==0),true,[1-readout_strength_params(4) readout_strength_params(4)]);
        choice(sp==0 & con==1)=randsample([0,1],nnz(sp==0 & con==1),true,[readout_strength_params(1) 1-readout_strength_params(1)]);
        choice(sp==1 & con==1)=randsample([0,1],nnz(sp==1 & con==1),true,[1-readout_strength_params(2) readout_strength_params(2)]);
end

