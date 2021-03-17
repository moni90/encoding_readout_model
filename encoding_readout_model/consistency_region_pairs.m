function [region] = consistency_region_pairs(ps_margin, consistency_thr)

id_keep = find(triu(ones(size(ps_margin,2)),1));
region = zeros(size(ps_margin,1),1);

for id_trial = 1:size(ps_margin,1)
    cons_pair = squareform(pdist(1*[ps_margin(id_trial,:)>0.5;ps_margin(id_trial,:)>0.5]'));
    frac_con = sum(cons_pair(id_keep)==0)/length(id_keep);
    if frac_con >= consistency_thr
        if mean(ps_margin(id_trial,:)>0.5)>0.5
            region(id_trial) = 1;
        elseif mean(ps_margin(id_trial,:)>0.5)<0.5
            region(id_trial) = 2;
        end
    end
end


    