function plot_consistency_model(data_folder, file_name, rho_idx, gamma_idx, consistency_strength_idx, save_path, save_name, save_data, save_figs)

curr_dir = pwd;
cd(data_folder);
load(file_name);
cd(curr_dir);


%% Figure 1: readout model characterization (example params)

frac_cons_real = simu_out.frac_cons{consistency_strength_idx,gamma_idx,rho_idx};
frac_cons_shuffled = simu_out.frac_cons_sh{consistency_strength_idx,gamma_idx,rho_idx};

figure;
axes; hold on;
boxplot([frac_cons_real frac_cons_shuffled],'symbol','');
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
box1 = a(5); box2 = a(6);
set(box1, 'Color', 'r');
set(box2, 'Color', 'b');
p_values=get_p_values(simu_out.frac_cons(consistency_strength_idx,gamma_idx,rho_idx), simu_out.frac_cons_sh(consistency_strength_idx,gamma_idx,rho_idx),'permtest');
text([1.5],nanmean([frac_cons_real; frac_cons_shuffled])+0.002,p_values,'Fontsize',14,'HorizontalAlignment','center');
set(gca,'XTick',[1.5]);
ylabel('Fraction consistent');
if save_figs
    saveas(gca,fullfile(save_path, save_name),'png');
    saveas(gca,fullfile(save_path, save_name),'fig');
    saveas(gca,fullfile(save_path, save_name),'svg');
end
if save_data
    save(fullfile(save_path, 'data_for_figures', ['DATA_' save_name '.mat']),...
        'frac_cons_real', 'frac_cons_shuffled','p_values');
    data_table = table(frac_cons_real, frac_cons_shuffled,...
        'VariableNames',{'frac_cons_real', 'frac_cons_shuffled'});
    writetable(data_table,fullfile(save_path, 'data_for_figures', ['DATA_' save_name '.csv']));
end
cd(curr_dir);

end