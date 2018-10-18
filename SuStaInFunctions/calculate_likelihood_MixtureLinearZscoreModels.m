function [loglike,total_prob_subj,total_prob_stage,total_prob_cluster,p_perm_k] = ...
    calculate_likelihood_MixtureLinearZscoreModels(data,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
    stage_zscore,stage_biomarker_index,S,f,likelihood_flag)
% Computes the likelihood of a mixture of linear z-score models using either
% an approximate method (faster, default setting) or an exact method
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
% S - the current ordering of the z-score stages for each subtype
%   dim: number of subtypes x number of z-score stages
% f - the current proportion of individuals belonging to each subtype
%   dim: number of subtypes x 1
% likelihood_flag - whether to use an exact method of inference - when set
% to 'Exact', the exact method is used, the approximate method is used for
% all other settings
%
%OUTPUTS:
% loglike - the log-likelihood of the current model
% total_prob_subj - the total probability of the current SuStaIn model for
% each subject
% total_prob_stage - the total probability of each stage in the current
% SuStaIn model
% total_prob_cluster - the total probability of each subtype in the current
% SuStaIn model
% p_perm_k - the probability of each subjects data at each stage of each
% subtype in the current SuStaIn model

M = size(data,1);
N_S = size(S,1);
N = size(stage_zscore,2);

f_val_mat = repmat(f,[1 N+1 M]);
f_val_mat = permute(f_val_mat,[3 2 1]);

p_perm_k = zeros(M,N+1,N_S);
for s = 1:N_S
    if (strcmp(likelihood_flag,'Exact'))
        p_perm_k(:,:,s) = ...
            calculate_likelihood_stage_LinearZscoreModel(data,...
            min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
            stage_zscore,stage_biomarker_index,S(s,:));
    else
        p_perm_k(:,:,s) = ...
            calculate_likelihood_stage_LinearZscoreModel_Approx(data,...
            min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
            stage_zscore,stage_biomarker_index,S(s,:));
    end
end

total_prob_cluster = squeeze(sum(p_perm_k.*f_val_mat,2));
total_prob_stage = sum(p_perm_k.*f_val_mat,3);
total_prob_subj = sum(total_prob_stage,2);

loglike = sum(log(total_prob_subj+1e-250));

end

