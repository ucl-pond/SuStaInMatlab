function [ml_sequence,ml_f,ml_likelihood,...
    ml_sequence_mat,ml_f_mat,ml_likelihood_mat] = ...
    find_ML_MixtureLinearZscoreModels(data,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
    stage_zscore,stage_biomarker_index,...
    seq_init,f_init,N_startpoints,likelihood_flag)
% Fit a mixture of linear z-score models
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
% seq_init - intial ordering of the stages for each subtype
% f_init - initial proprtion of individuals belonging to each subtype
% N_startpoints - the number of start points for the fitting
% likelihood_flag - whether to use an exact method of inference - when set
% to 'Exact', the exact method is used, the approximate method is used for
% all other settings
%
%OUTPUTS:
% ml_sequence - the ordering of the stages for each subtype for the next
% SuStaIn model in the hierarchy
% ml_f - the most probable proportion of individuals belonging to each
% subtype for the next SuStaIn model in the hierarchy
% ml_likelihood - the likelihood of the most probable SuStaIn model for the
% next SuStaIn model in the hierarchy
% previous outputs _mat - same as before but for each start point

N_S = size(seq_init,1);

terminate = 0;
startpoint = 0;

ml_sequence_mat = zeros(N_S,size(stage_zscore,2),N_startpoints);
ml_f_mat = zeros(N_S,N_startpoints);
ml_likelihood_mat = zeros(N_startpoints,1);
while (terminate==0)
    %tic
    startpoint = startpoint+1;
    fprintf(' ++ startpoint %d\n',startpoint)
    
    [this_ml_sequence,this_ml_f,this_ml_likelihood] = ...
        perform_EM_MixtureLinearZscoreModels(data,...
        min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
        stage_zscore,stage_biomarker_index,seq_init,f_init,likelihood_flag);
    ml_sequence_mat(:,:,startpoint) = this_ml_sequence;
    ml_f_mat(:,startpoint) = this_ml_f;
    ml_likelihood_mat(startpoint) = this_ml_likelihood;
    
    if (startpoint>=N_startpoints)
        terminate = 1;
    end
    %toc
end

ix = find(ml_likelihood_mat==max(ml_likelihood_mat));
ix = ix(1);
ml_sequence = ml_sequence_mat(:,:,ix);
ml_f = ml_f_mat(:,ix);
ml_likelihood = ml_likelihood_mat(ix);

end

