function plot_correlation_correct_error_model(data_folder, file_name, rho_idx, gamma_idx, consistency_strength_idx, readout, save_path, save_name, save_data, save_figs)

curr_dir = pwd;
cd(data_folder);
load(file_name);
cd(curr_dir);


%% Figure 1: readout model characterization (example params)
switch readout
    case 'enhanced-by-consistency'
        noiseCorr_correct = simu_out.pair_corr_coeff_BP_correct{consistency_strength_idx,gamma_idx,rho_idx};
        noiseCorr_error = simu_out.pair_corr_coeff_BP_error{consistency_strength_idx,gamma_idx,rho_idx};
        p_values=get_p_values(simu_out.pair_corr_coeff_BP_correct(consistency_strength_idx,gamma_idx,rho_idx), simu_out.pair_corr_coeff_BP_error(consistency_strength_idx,gamma_idx,rho_idx),'permtest');
%         noiseCorr_correct = simu_out.pop_corr_coeff_BP_correct{consistency_strength_idx,gamma_idx,rho_idx};
%         noiseCorr_error = simu_out.pop_corr_coeff_BP_error{consistency_strength_idx,gamma_idx,rho_idx};
%         p_values=get_p_values(simu_out.pop_corr_coeff_BP_correct(consistency_strength_idx,gamma_idx,rho_idx), simu_out.pop_corr_coeff_BP_error(consistency_strength_idx,gamma_idx,rho_idx),'permtest');
    
    case 'consistency-independent'
        noiseCorr_correct = simu_out.pair_corr_coeff_BP_correct_subopt{consistency_strength_idx,gamma_idx,rho_idx};
        noiseCorr_error = simu_out.pair_corr_coeff_BP_error_subopt{consistency_strength_idx,gamma_idx,rho_idx};
        p_values=get_p_values(simu_out.pair_corr_coeff_BP_correct_subopt(consistency_strength_idx,gamma_idx,rho_idx), simu_out.pair_corr_coeff_BP_error_subopt(consistency_strength_idx,gamma_idx,rho_idx),'permtest');
%         noiseCorr_correct = simu_out.pop_corr_coeff_BP_correct_subopt{consistency_strength_idx,gamma_idx,rho_idx};
%         noiseCorr_error = simu_out.pop_corr_coeff_BP_error_subopt{consistency_strength_idx,gamma_idx,rho_idx};
%         p_values=get_p_values(simu_out.pop_corr_coeff_BP_correct_subopt(consistency_strength_idx,gamma_idx,rho_idx), simu_out.pop_corr_coeff_BP_error_subopt(consistency_strength_idx,gamma_idx,rho_idx),'permtest');

end

figure;
axes; hold on;
boxplot([noiseCorr_correct noiseCorr_error],'symbol','');
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
box1 = a(5); box2 = a(6);
set(box1, 'Color', 'r');
set(box2, 'Color', 'b');
text([1.5],nanmean([noiseCorr_correct; noiseCorr_error])+0.002,p_values,'Fontsize',14,'HorizontalAlignment','center');
set(gca,'XTick',[1.5]);
ylabel('Pearson correlations');
if save_figs
    saveas(gca,fullfile(save_path, save_name),'png');
    saveas(gca,fullfile(save_path, save_name),'fig');
    saveas(gca,fullfile(save_path, save_name),'svg');
end
if save_data
    save(fullfile(save_path, 'data_for_figures', ['DATA_' save_name '.mat']),...
        'noiseCorr_correct', 'noiseCorr_error','p_values');
    data_table = table(noiseCorr_correct, noiseCorr_error,...
        'VariableNames',{'noiseCorr_correct', 'noiseCorr_error'});
    writetable(data_table,fullfile(save_path, 'data_for_figures', ['DATA_' save_name '.csv']))
end
cd(curr_dir);

end