function [ml_sequence,ml_f,ml_likelihood,...
    ml_sequence_mat,ml_f_mat,ml_likelihood_mat] = ...
    find_ML_Mixture2LinearZscoreModels(data,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
    stage_zscore,stage_biomarker_index,N_startpoints,likelihood_flag)
% Fit a mixture of two linear z-score models
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
% N_startpoints - the number of start points for the fitting
% likelihood_flag - whether to use an exact method of inference - when set
% to 'Exact', the exact method is used, the approximate method is used for
% all other settings
%
%OUTPUTS:
% ml_sequence - the ordering of the stages for each subtype
% ml_f - the most probable proportion of individuals belonging to each
% subtype
% ml_likelihood - the likelihood of the most probable SuStaIn model
% previous outputs _mat - same as before but for each start point

N_S = 2;

terminate = 0;
startpoint = 0;

ml_sequence_mat = zeros(N_S,size(stage_zscore,2),N_startpoints);
ml_f_mat = zeros(N_S,N_startpoints);
ml_likelihood_mat = zeros(N_startpoints,1);
while (terminate==0)
    %tic
    startpoint = startpoint+1;
    fprintf(' ++ startpoint %d\n',startpoint)
    
    % randomly initialise individuals as belonging to one of the two
    % subtypes (clusters)
    min_N_cluster = 0;
    while (min_N_cluster==0)
        cluster_assignment = ceil(N_S*rand(size(data,1),1));
        temp_N_cluster = zeros(1,N_S);
        for s = 1:N_S
            temp_N_cluster = sum(cluster_assignment==s);
        end
        min_N_cluster = min(temp_N_cluster);
    end
    % initialise the stages of the two linear z-score models by fitting a
    % single linear z-score model to each of the two sets of individuals
    seq_init = zeros(N_S,size(stage_zscore,2));
    for s = 1:N_S
        temp_data = data(cluster_assignment==s,:);
        temp_seq_init = initialise_sequence_LinearZscoreModel(stage_zscore,stage_biomarker_index);
        seq_init(s,:) = perform_EM_MixtureLinearZscoreModels(temp_data,...
            min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
            stage_zscore,stage_biomarker_index,temp_seq_init,1,likelihood_flag);
    end
    f_init = 1/N_S*ones(N_S,1);
    % optimise the mixture of two linear z-score models from the
    % initialisation
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

