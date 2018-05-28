function [ml_sequence,ml_f,ml_likelihood,...
    samples_sequence,samples_f,samples_likelihood] = ...
    perform_MCMC_MixtureLinearZscoreModels(data,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
    stage_zscore,stage_biomarker_index,...
    seq_init,f_init,n_iterations,seq_sigma,f_sigma,likelihood_flag)
% Take MCMC samples of the uncertainty in the SuStaIn model parameters

N = size(stage_zscore,2);
N_S = size(seq_init,1);

samples_sequence = zeros(N_S,N,n_iterations);
samples_f = zeros(N_S,n_iterations);
samples_likelihood = zeros(n_iterations,1);

samples_sequence(:,:,1) = seq_init;
samples_f(:,1) = f_init;

%tic
for i = 1:n_iterations
    if (mod(i,n_iterations/100)==1)
        %toc
        fprintf('Iteration %d of %d, %d%% complete\n',i,n_iterations,(i-1)/n_iterations*100)
        %tic
    end
    if (i>1)
        seq_order = randperm(N_S);
        for s = seq_order
            move_event_from = ceil(N*rand);
            
            current_sequence = samples_sequence(s,:,i-1);
            current_location(current_sequence) = 1:N;
            
            selected_event = current_sequence(move_event_from);
            
            this_stage_zscore = stage_zscore(selected_event);
            selected_biomarker = stage_biomarker_index(selected_event);
            possible_zscores_biomarker = stage_zscore(stage_biomarker_index==selected_biomarker);
            
            min_zscore_bound = max(possible_zscores_biomarker(possible_zscores_biomarker<this_stage_zscore));
            max_zscore_bound = min(possible_zscores_biomarker(possible_zscores_biomarker>this_stage_zscore));
            events = 1:N;
            
            if (~isempty(min_zscore_bound))
                min_zscore_bound_event = events(((stage_zscore==min_zscore_bound)+(stage_biomarker_index==selected_biomarker))==2);
                move_event_to_lower_bound = current_location(min_zscore_bound_event)+1;
            else
                move_event_to_lower_bound = 1;
            end
            if (~isempty(max_zscore_bound))
                max_zscore_bound_event = events(((stage_zscore==max_zscore_bound)+(stage_biomarker_index==selected_biomarker))==2);
                move_event_to_upper_bound = current_location(max_zscore_bound_event)-1;
            else
                move_event_to_upper_bound = N;
            end
            
            possible_positions = move_event_to_lower_bound:move_event_to_upper_bound;
            distance = possible_positions-move_event_from;
            
            if ((size(seq_sigma,1)==1)+(size(seq_sigma,2)==1)==2)
                this_seq_sigma = seq_sigma;
            else
                this_seq_sigma = seq_sigma(s,selected_event);
            end
            
            weight = normpdf(distance,0,this_seq_sigma);
            weight = weight./sum(weight);
            index = randsample(1:length(possible_positions),1,true,weight);
            
            move_event_to = possible_positions(index);
            
            current_sequence(move_event_from) = [];
            new_sequence = [current_sequence(1:move_event_to-1) selected_event current_sequence(move_event_to:N-1)];
            
            samples_sequence(s,:,i) = new_sequence;
        end
        new_f = samples_f(:,i-1) + f_sigma.*randn(N_S,1);
        new_f = abs(new_f)/sum(abs(new_f));
        samples_f(:,i) = new_f;
    end
    
    S = samples_sequence(:,:,i);
    f = samples_f(:,i);
    
    likelihood_sample = ...
        calculate_likelihood_MixtureLinearZscoreModels(data,...
        min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
        stage_zscore,stage_biomarker_index,S,f,likelihood_flag);
    
    samples_likelihood(i) = likelihood_sample;
    
    if(i>1)
        ratio = exp(samples_likelihood(i)-samples_likelihood(i-1));
        if (ratio<rand)
            samples_likelihood(i) = samples_likelihood(i-1);
            samples_sequence(:,:,i) = samples_sequence(:,:,i-1);
            samples_f(:,i) = samples_f(:,i-1);
        end
    end
    
end
%toc

perm_index = find(samples_likelihood==max(samples_likelihood));
perm_index = perm_index(1);
ml_likelihood = max(samples_likelihood);
ml_sequence = samples_sequence(:,:,perm_index);
ml_f = samples_f(:,perm_index);

end