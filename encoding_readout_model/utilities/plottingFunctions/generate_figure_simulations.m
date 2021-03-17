function generate_figure_simulations(which_panel, data_folder, file_name, N, consistency_idx, save_path, save_data, save_figs)

switch which_panel
    case 1
        save_name = ['N' num2str(N) '_fig3_a_heatmap_S_decoder'];
        heatmap_decoder_performance(data_folder, N, file_name, consistency_idx, save_path, save_name, save_data, save_figs);
    case 2
        save_name = ['N' num2str(N) '_fig3_b_heatmap_consistency'];
        heatmap_consistency(data_folder, N, file_name, consistency_idx, save_path, save_name, save_data, save_figs);
    case 3
        consistency_idx = 0;
        save_name = ['N' num2str(N) '_fig3_c_heatmap_Pearson'];
        heatmap_correlation(data_folder, N, file_name, consistency_idx, save_path, save_name, save_data, save_figs);
    case 4
        consistency_idx = 0;
        save_name = ['N' num2str(N) '_fig3_d_heatmap_popCorr'];
        heatmap_pop_correlation(data_folder, N, file_name, consistency_idx, save_path, save_name, save_data, save_figs);
    case 5
        consistency_idx = 0;
        save_name = ['N' num2str(N) '_fig3_e_heatmap_BP'];
        heatmap_BP_performance(data_folder, N, file_name, consistency_idx, save_path, save_name, save_data, save_figs);
    case 6
        save_name = ['N' num2str(N) '_fig3_c_heatmap_Pearson'];
        heatmap_correlation(data_folder, N, file_name, consistency_idx, save_path, save_name, save_data, save_figs);
    case 7
        save_name = ['N' num2str(N) '_fig3_d_heatmap_popCorr'];
        heatmap_pop_correlation(data_folder, N, file_name, consistency_idx, save_path, save_name, save_data, save_figs);
    case 8
        save_name = ['N' num2str(N) '_fig3_e_heatmap_BP'];
        heatmap_BP_performance(data_folder, N, file_name, consistency_idx, save_path, save_name, save_data, save_figs);
    case 9
        save_name = ['N' num2str(N) '_fig3_i_heatmap_BP_readout'];
        heatmap_BP_readout(data_folder, N, file_name, consistency_idx, save_path, save_name, save_data, save_figs);     
end

