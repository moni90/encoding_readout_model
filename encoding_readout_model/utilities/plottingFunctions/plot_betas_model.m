function plot_betas_model(data_folder, file_name, rho_idx, gamma_idx, consistency_idx, readout, save_path, save_name, save_data, save_figs)

curr_dir = pwd;
cd(data_folder);
load(file_name);
cd('..');


figure(2); clf;
set(gcf,'Units','pixels','Position',[100 100 600 400]);

betas_color = [0 0 0];
line_color = [0 0 0];
off_x = 7;
off_y = 8;
% gamma_idx = 1;

readout_strength_params = [1-consistency_strength_range(consistency_idx) 1-consistency_strength_range(consistency_idx) ...
    .5+consistency_strength_range(consistency_idx) .5+consistency_strength_range(consistency_idx)];

switch readout
    case 'enhanced-by-consistency'
        readout_strength_params_cons = readout_strength_params;
        betas = determine_betas_from_readout_strength_params_model2(readout_strength_params_cons);
    case 'consistency-independent'
        readout_strength_params_eqsub = [1 1 1 1]*mean(simu_out.avg_readout_strength{consistency_idx,gamma_idx,rho_idx});
        betas = determine_betas_from_readout_strength_params_model2(readout_strength_params_eqsub);
end

% cons readout
axes; hold on;

bar_width = .6;

make_my_bar_nan(1, betas(1), betas_color, line_color, bar_width);
make_my_bar_nan(2, betas(2), betas_color, line_color, bar_width);
make_my_bar_nan(3, betas(3), betas_color, line_color, bar_width);
make_my_bar_nan(4, betas(4), betas_color, line_color, bar_width);

set(gca,'XTick',[1 2 3 4],'XTickLabel',{'\beta_{bias}','\beta_{sp}','\beta_{int1}','\beta_{int2}'},'XTickLabelRotation',0)

my_axis(gca,false,off_x,0)

ylim([0 5]);
xlim([0 5]);
if save_figs
    saveas(gca,fullfile(save_path, save_name),'png');
    saveas(gca,fullfile(save_path, save_name),'fig');
    saveas(gca,fullfile(save_path, save_name),'svg');
end
if save_data
    save(fullfile(save_path, 'data_for_figures', ['DATA_' save_name '.mat']),...
        'betas');
    data_table = table(betas(1),betas(2),betas(3),betas(4),...
        'VariableNames',{'beta_0', 'beta_s', 'beta_i1', 'beta_i2'});
    writetable(data_table,fullfile(save_path, 'data_for_figures', ['DATA_' save_name '.csv']));
end
cd(curr_dir);
