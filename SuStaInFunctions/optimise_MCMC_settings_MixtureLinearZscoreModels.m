function [seq_sigma_opt,f_sigma_opt] = optimise_MCMC_settings_MixtureLinearZscoreModels(data,...
        min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
        stage_zscore,stage_biomarker_index,...
        seq_init,f_init,likelihood_flag)
% Optimise the perturbation size for the MCMC algorithm
    
n_iterations_MCMC_optimisation = 1e4;
n_passes_optimisation = 3;

seq_sigma_currentpass = 1;
f_sigma_currentpass = 0.01;

N_S = size(seq_init,1);
    
for i = 1:n_passes_optimisation
    [~,~,~,...
        samples_sequence_currentpass,samples_f_currentpass] = ...
        perform_MCMC_MixtureLinearZscoreModels(data,...
        min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
        stage_zscore,stage_biomarker_index,...
        seq_init,f_init,n_iterations_MCMC_optimisation,seq_sigma_currentpass,f_sigma_currentpass,likelihood_flag);
    
    samples_position_currentpass = zeros(size(samples_sequence_currentpass));
    for s = 1:N_S
        for sample = 1:n_iterations_MCMC_optimisation
            temp_seq = samples_sequence_currentpass(s,:,sample);
            temp_inv(temp_seq) = 1:size(samples_sequence_currentpass,2);
            samples_position_currentpass(s,:,sample) = temp_inv;
        end
    end
    
    seq_sigma_currentpass = std(samples_position_currentpass,[],3);
    seq_sigma_currentpass(seq_sigma_currentpass<0.01) = 0.01;
    f_sigma_currentpass = std(samples_f_currentpass,[],2);    
end

seq_sigma_opt = seq_sigma_currentpass;
f_sigma_opt = f_sigma_currentpass;

end

