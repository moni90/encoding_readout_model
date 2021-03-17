function heatmap_consistency(data_folder, N, file_name, consistency_strength_idx, save_path, save_name, save_data, save_figs)

plot_Ni = 1;

curr_dir = pwd;
cd(data_folder);
load(file_name);
cd(curr_dir);

if consistency_strength_idx<length(consistency_strength_range)
   readoutStr = '_consEnhanced' ;
else
   readoutStr = '_consIndependent' ;
end

% if plot_Ni == 1
%     Ni_data_nobias = [22.61-0.10 22.61+0.10]/100;
%     Ni_data_nobias_morcos = [19.54-0.12 19.54+0.12]/100;
% else
%     Ni_data_nobias = [0.03-0.001 0.03+0.001];
%     Ni_data_nobias_morcos = [0.012-0.001 0.012+0.001];
% end
% gamma_data = [1.0492-0.0041 1.0492+0.0041]/pi; %bias corrected
% gamma_runyan_m = 1.0492;
% gamma_data_morcos = [0.9226-0.0067 0.9226+0.0067]/pi;  %bias corrected
% gamma_morcos_m = 0.9226;

if plot_Ni
    rho_range = (1+rho_range*(2*N-1))/(2*N);
end

selected_gamma_range = 1:length(gamma_range);

%% Density plots at varying gamma rho and cons strength - correlated vs shuffled
frac_cons = simu_out.frac_cons(:,:,:);
frac_cons_sh = simu_out.frac_cons_sh(:,:,:);

figure;
delta_frac_cons_readout = cellfun(@minus,squeeze(frac_cons(consistency_strength_idx,:,:))',squeeze(frac_cons_sh(consistency_strength_idx,:,:))','UniformOutput',false);
delta_frac_cons_readout_mean = cellfun(@mean,delta_frac_cons_readout);
% delta_frac_cons_readout_mean = cellfun(@mean,squeeze(frac_cons(consistency_strength_idx,:,:))');
clim=max(abs(delta_frac_cons_readout_mean(:)));
imagesc(gamma_range,rho_range,delta_frac_cons_readout_mean)
% colormap(gca,myColorMapRedBlue_log);
set(gca,'YDir','reverse')
% caxis([-clim clim]);
axis square; axis tight;
title('Frac. consistent real - shuffled');
xlabel('\gamma');
if plot_Ni
    ylabel('P(var expl by 1st PC)');
else
    ylabel('\rho');
end
colorbar;
box off
set(gca,'Fontsize',9);
set(gca,'XTick',[0 1/6 1/4 1/3],'XTickLabels',{'0','\pi/6','\pi/4','\pi/3'});
% hold on;
% rectangle('Position',[gamma_data(1) Ni_data_nobias(1) diff(gamma_data) diff(Ni_data_nobias)],...
%     'EdgeColor','m','FaceColor','none');
% rectangle('Position',[gamma_data_morcos(1) Ni_data_nobias_morcos(1) diff(gamma_data_morcos) diff(Ni_data_nobias_morcos)],...
%     'EdgeColor','m','FaceColor','none');
if save_figs
    saveas(gca,[save_path save_name readoutStr],'png');
    saveas(gca,[save_path save_name readoutStr],'fig');
    saveas(gca,[save_path save_name readoutStr],'svg');
end
if save_data
    save(fullfile(save_path, 'data_for_figures', ['DATA_' save_name '.mat']),...
        'gamma_range','rho_range','delta_frac_cons_readout_mean');
    writematrix(delta_frac_cons_readout_mean,fullfile(save_path, 'data_for_figures', ['DATA_' save_name '.csv']));
end
cd(curr_dir);

end