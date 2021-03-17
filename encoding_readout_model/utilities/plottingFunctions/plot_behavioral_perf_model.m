function plot_behavioral_perf_model(data_folder, file_name, rho_idx, gamma_idx, consistency_strength_idx, readout, save_path, save_name, save_data, save_figs)

curr_dir = pwd;
cd(data_folder);
load(file_name);
cd(curr_dir);


%% Figure 1: readout model characterization (example params)
switch readout
    case 'enhanced-by-consistency'
        BP_real = simu_out.behav_perf{consistency_strength_idx,gamma_idx,rho_idx};
        BP_shuffled = simu_out.behav_perf_sh{consistency_strength_idx,gamma_idx,rho_idx};
        p_values=get_p_values(simu_out.behav_perf(consistency_strength_idx,gamma_idx,rho_idx), simu_out.behav_perf_sh(consistency_strength_idx,gamma_idx,rho_idx),'permtest');
    case 'consistency-independent'
        BP_real = simu_out.behav_perf_subopt{consistency_strength_idx,gamma_idx,rho_idx};
        BP_shuffled = simu_out.behav_perf_subopt_sh{consistency_strength_idx,gamma_idx,rho_idx};
        p_values=get_p_values(simu_out.behav_perf_subopt(consistency_strength_idx,gamma_idx,rho_idx), simu_out.behav_perf_subopt_sh(consistency_strength_idx,gamma_idx,rho_idx),'permtest');
end

figure;
axes; hold on;
boxplot([BP_real BP_shuffled],'symbol','');
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
box1 = a(5); box2 = a(6);
set(box1, 'Color', 'r');
set(box2, 'Color', 'b');
text([1.5],nanmean([BP_real; BP_shuffled])+0.002,p_values,'Fontsize',14,'HorizontalAlignment','center');
set(gca,'XTick',[1.5]);
ylabel('Fraction correct');
if save_figs
    saveas(gca,fullfile(save_path, save_name),'png');
    saveas(gca,fullfile(save_path, save_name),'fig');
    saveas(gca,fullfile(save_path, save_name),'svg');
end
if save_data
    save(fullfile(save_path, 'data_for_figures', ['DATA_' save_name '.mat']),...
        'BP_real', 'BP_shuffled','p_values');
    data_table = table(BP_real, BP_shuffled,...
        'VariableNames',{'BP_real', 'BP_shuffled'});
    writetable(data_table,fullfile(save_path, 'data_for_figures', ['DATA_' save_name '.csv']));
end
cd(curr_dir);

end