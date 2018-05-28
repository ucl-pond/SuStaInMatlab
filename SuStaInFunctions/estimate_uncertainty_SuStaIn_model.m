function [ml_sequence,ml_f,ml_likelihood,...
    samples_sequence,samples_f,samples_likelihood] = estimate_uncertainty_SuStaIn_model(data,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
    stage_zscore,stage_biomarker_index,...
    seq_init,f_init,N_iterations_MCMC,likelihood_flag)
% Estimate the uncertainty in the subtype progression patterns and
% proportion of individuals belonging to the SuStaIn model
%
%INPUTS:
% data - !important! needs to be (positive) z-scores! 
%   dim: number of subjects x number of biomarkers
% min_biomarker_zscore - a minimum z-score for each biomarker (usually zero
% for all markers)
%   dim: 1 x number of biomarkers
% max_biomarker_zscore - a maximum z-score for each biomarker - reached at
% the final stage of the linear z-score model
%   dim: 1 x number of biomarkers
% std_biomarker_zscore - the standard devation of each biomarker z-score
% (should be 1 for all markers)
%   dim: 1 x number of biomarkers
% stage_zscore and stage_biomarker_index give the different z-score stages
% for the linear z-score model, i.e. the index of the different z-scores
% for each biomarker
% stage_zscore - the different z-scores of the model
%   dim: 1 x number of z-score stages
% stage_biomarker_index - the index of the biomarker that the corresponding
% entry of stage_zscore is referring to - !important! ensure biomarkers are
% indexed s.t. they correspond to columns 1 to number of biomarkers in your
% data
%   dim: 1 x number of z-score stages
% seq_init - the ordering of the stages for each subtype to initialise the
% MCMC estimation from
%   dim: number of subtypes x number of z-score stages
% f_init - the proportion of individuals belonging to each subtype to
% intialise the MCMC esimation from
%   dim: number of subtypes x 1
% N_iterations_MCMC - the number of MCMC samples to take
% likelihood_flag - whether to use an exact method of inference - when set
% to 'Exact', the exact method is used, the approximate method is used for
% all other settings
%
%OUTPUTS:
% ml_sequence - the most probable ordering of the stages for each subtype
% found across MCMC samples
% ml_f - the most probable proportion of individuals belonging to each
% subtype found across MCMC samples
% ml_likelihood - the likelihood of the most probable SuStaIn model found
% across MCMC samples
% samples_sequence - samples of the ordering of the stages for each subtype
% obtained from MCMC sampling
% samples_f - samples of the proportion of individuals belonging to each
% subtype obtained from MCMC sampling
% samples_likeilhood - samples of the likelihood of each SuStaIn model
% sampled by the MCMC sampling

% Perform a few initial passes where the perturbation sizes of the MCMC
% unertainty estimation are tuned
[seq_sigma_opt,f_sigma_opt] = optimise_MCMC_settings_MixtureLinearZscoreModels(data,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
    stage_zscore,stage_biomarker_index,...
    seq_init,f_init,likelihood_flag);

% Run the full MCMC algorithm to estimate the uncertainty
[ml_sequence,ml_f,ml_likelihood,...
    samples_sequence,samples_f,samples_likelihood] = ...
    perform_MCMC_MixtureLinearZscoreModels(data,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
    stage_zscore,stage_biomarker_index,...
    seq_init,f_init,N_iterations_MCMC,seq_sigma_opt,f_sigma_opt,likelihood_flag);

end

