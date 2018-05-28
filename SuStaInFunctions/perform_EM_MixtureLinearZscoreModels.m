function [ml_sequence,ml_f,ml_likelihood,...
        samples_sequence,samples_f,samples_likelihood] = ...
    perform_EM_MixtureLinearZscoreModels(data,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
    stage_zscore,stage_biomarker_index,current_sequence,current_f,likelihood_flag)
% Perform an E-M procedure to estimate parameters of SuStaIn model

MaxIter = 100;

N = size(stage_zscore,2);
N_S = size(current_sequence,1);

current_likelihood = ...
    calculate_likelihood_MixtureLinearZscoreModels(data,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
    stage_zscore,stage_biomarker_index,current_sequence,current_f,likelihood_flag);

terminate = 0;
iteration = 0;

samples_sequence = NaN*ones(MaxIter,N,N_S);
samples_f = NaN*ones(MaxIter,N_S);
samples_likelihood = NaN*ones(MaxIter,1);

samples_sequence(1,:,:) = current_sequence';
samples_f(1,:) = current_f;
samples_likelihood(1) = current_likelihood;
while (terminate==0)
    iteration = iteration + 1;
    % fprintf('++ iteration %d\n',iteration)
    [candidate_sequence,candidate_f,candidate_likelihood] = ...
        optimise_parameters_MixtureLinearZscoreModels(data,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
    stage_zscore,stage_biomarker_index,current_sequence,current_f,likelihood_flag);
    
    HAS_converged = abs((candidate_likelihood-current_likelihood)/max(candidate_likelihood,current_likelihood))<1e-6;
    if HAS_converged
        % fprintf('EM converged in %d iterations\n',iteration)
        terminate = 1;
    else
        if (candidate_likelihood>current_likelihood)
            current_sequence = candidate_sequence;
            current_f = candidate_f;
            current_likelihood = candidate_likelihood;
        end
    end
    samples_sequence(iteration,:,:) = current_sequence';
    samples_f(iteration,:) = current_f;
    samples_likelihood(iteration) = current_likelihood;
    
    if (iteration>=MaxIter)
        terminate = 1;
    end
end

ml_sequence = current_sequence;
ml_f = current_f;
ml_likelihood = current_likelihood;

end