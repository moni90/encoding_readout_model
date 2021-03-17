function generate_figure_2(which_panel, data_folder, file_name, rho_idx, gamma_idx, consistency_idx, save_path, save_data, save_figs)
%function to generate panels in Fig1
        
switch which_panel
    case 4
        save_name = 'fig2_d_S_decoder';
        plot_decoder_performance_model(data_folder, file_name, rho_idx, gamma_idx, consistency_idx, save_path, save_name, save_data, save_figs);
    case 5
        save_name = 'fig2_e_frac_consistent';
        plot_consistency_model(data_folder, file_name, rho_idx, gamma_idx, consistency_idx, save_path, save_name, save_data, save_figs);
    case 6
        readout = 'consistency-independent';
        save_name = ['fig2_f_betas_' readout];
        plot_betas_model(data_folder, file_name, rho_idx, gamma_idx, consistency_idx, readout, save_path, save_name, save_data, save_figs);
    case 7
        readout = 'consistency-independent';
        save_name = ['fig2_g_noise_correlations_' readout];
        plot_correlation_correct_error_model(data_folder, file_name, rho_idx, gamma_idx, consistency_idx, readout, save_path, save_name, save_data, save_figs);
    case 8
        readout = 'consistency-independent';
        save_name = ['fig2_h_BP_' readout];
        plot_behavioral_perf_model(data_folder, file_name, rho_idx, gamma_idx, consistency_idx, readout, save_path, save_name, save_data, save_figs);
    case 9
        readout = 'enhanced-by-consistency';
        save_name = ['fig2_i_betas_' readout];
        plot_betas_model(data_folder, file_name, rho_idx, gamma_idx, consistency_idx, readout, save_path, save_name, save_data, save_figs);
    case 10
        readout = 'enhanced-by-consistency';
        save_name = ['fig2_j_noise_correlations_' readout];
        plot_correlation_correct_error_model(data_folder, file_name, rho_idx, gamma_idx, consistency_idx, readout, save_path, save_name, save_data, save_figs);
    case 11
        readout = 'enhanced-by-consistency';
        save_name = ['fig2_k_BP_' readout];
        plot_behavioral_perf_model(data_folder, file_name, rho_idx, gamma_idx, consistency_idx, readout, save_path, save_name, save_data, save_figs);
end

