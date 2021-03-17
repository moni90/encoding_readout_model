function [trials_sub, fraction_discarded, lowInfo] = subsample_trials_to_eqShuffledInfo(correct_trials,correct_trials_sh)

% note that correct trials refers to correct decodings

% Objective: equalizing corr/all - corr_sh/all

removed_trials = false(size(correct_trials));

I = mean(correct_trials);
I_sh = mean(correct_trials_sh);
I_curr = I;
I_sh_curr = I_sh;
correct_trials_curr = correct_trials;
correct_trials_sh_curr = correct_trials_sh;

lowInfo = I<.5 | I_sh<.5;

% first compute which of the 2 fractions is larger
[h_value,h_idx] = max([I I_sh]);

% then subsample trials accordingly

switch h_idx
    
    case 1
        
        while I_curr > I_sh_curr 
            
            trial_to_remove = randsample(find(correct_trials_curr & ~correct_trials_sh_curr),1);
            correct_trials_curr(trial_to_remove) = false;
            correct_trials_sh_curr(trial_to_remove) = true;
            
            removed_trials(trial_to_remove)=true;
            I_curr = mean(correct_trials(~removed_trials));
            I_sh_curr = mean(correct_trials_sh(~removed_trials));
        end
        
    case 2
        
        while I_sh_curr > I_curr
            
            trial_to_remove = randsample(find(correct_trials_sh_curr & ~correct_trials_curr),1);
            correct_trials_sh_curr(trial_to_remove) = false;
            correct_trials_curr(trial_to_remove) = true;
            
            removed_trials(trial_to_remove)=true;
            I_curr = mean(correct_trials(~removed_trials));
            I_sh_curr = mean(correct_trials_sh(~removed_trials));
        end
        
end

trials_sub = ~removed_trials;
fraction_discarded = mean(removed_trials);

%if I_curr<0.3
%    keyboard
%end

end