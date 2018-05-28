function [ml_sequence,ml_f,ml_likelihood,...
    ml_sequence_mat,ml_f_mat,ml_likelihood_mat] = ...
    estimate_ML_SuStaIn_model_Nplus1_Clusters(data,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
    stage_zscore,stage_biomarker_index,...
    ml_sequence_prev,ml_f_prev,...
    N_startpoints,likelihood_flag)
% Given the previous SuStaIn model, estimate the next model in the
% hierarchy (i.e. number of subtypes goes from N to N+1)
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
% ml_sequence_prev - the ordering of the stages for each subtype from the
% previous SuStaIn model
%   dim: number of subtypes x number of z-score stages
% ml_f_prev - the proportion of individuals belonging to each subtype from
% the previous SuStaIn model
%   dim: number of subtypes x 1
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

N_S = size(ml_sequence_prev,1)+1;

if (N_S==1)
    % If the number of subtypes is 1, fit a single linear z-score model
    fprintf('Finding ML solution to 1 cluster problem\n')
    [ml_sequence,ml_f,ml_likelihood,...
        ml_sequence_mat,ml_f_mat,ml_likelihood_mat] = ...
        find_ML_LinearZscoreModel(data,...
        min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
        stage_zscore,stage_biomarker_index,N_startpoints,likelihood_flag);
    fprintf('Overall ML likelihood is %f\n',ml_likelihood)
else
    % If the number of subtypes is greater than 1, go through each subtype
    % in turn and try splitting into two subtypes
    [~,~,~,p_sequence] = ...
        calculate_likelihood_MixtureLinearZscoreModels(data,...
        min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
        stage_zscore,stage_biomarker_index,ml_sequence_prev,ml_f_prev,likelihood_flag);
    p_sequence_norm = p_sequence./repmat(sum(p_sequence,2),[1 N_S-1]);
    
    % Assign individuals to a subtype (cluster) based on the previous model
    ml_cluster_subj = zeros(size(data,1),1);
    for m = 1:size(data,1)
        [~,ix] = max(p_sequence_norm(m,:));
        ml_cluster_subj(m) = ix(ceil(rand*length(ix)));
    end
    
    ml_likelihood = -Inf;
    for ix_cluster_split = 1:(N_S-1)
        this_N_cluster = sum(ml_cluster_subj==ix_cluster_split);
        if (this_N_cluster>1)
            % Take the data from the individuals belonging to a particular
            % cluster and fit a two subtype model
            fprintf('Splitting cluster %d of %d\n',ix_cluster_split,N_S-1)
            data_split = data(ml_cluster_subj==ix_cluster_split,:);
            fprintf(' + Resolving 2 cluster problem\n')
            this_ml_sequence_split = ...
                find_ML_Mixture2LinearZscoreModels(data_split,...
                min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
                stage_zscore,stage_biomarker_index,N_startpoints,likelihood_flag);
            % Use the two subtype model combined with the other subtypes to
            % inititialise the fitting of the next SuStaIn model in the
            % hierarchy
            this_seq_init = ml_sequence_prev;
            this_seq_init(ix_cluster_split,:) = this_ml_sequence_split(1,:);
            this_seq_init = [this_seq_init; this_ml_sequence_split(2,:)];
            this_f_init = 1/N_S*ones(N_S,1);
            fprintf(' + Finding ML solution from hierarchical initialisation\n')
            [this_ml_sequence,this_ml_f,this_ml_likelihood,...
                this_ml_sequence_mat,this_ml_f_mat,this_ml_likelihood_mat] = ...
                find_ML_MixtureLinearZscoreModels(data,...
                min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
                stage_zscore,stage_biomarker_index,this_seq_init,this_f_init,N_startpoints,likelihood_flag);
            
            % Choose the most probable SuStaIn model from the different
            % possible SuStaIn models initialised by splitting each subtype
            % in turn
            if (this_ml_likelihood>ml_likelihood)
                ml_likelihood = this_ml_likelihood;
                ml_sequence = this_ml_sequence;
                ml_f = this_ml_f;
                ml_likelihood_mat = this_ml_likelihood_mat;
                ml_sequence_mat = this_ml_sequence_mat;
                ml_f_mat = this_ml_f_mat;
            end
            fprintf('- ML likelihood is %f\n',this_ml_likelihood)
        else
            fprintf('Cluster %d of %d too small for subdivision\n',ix_cluster_split,N_S-1)
        end
    end
    fprintf('Overall ML likelihood is %f\n',ml_likelihood)
end

end