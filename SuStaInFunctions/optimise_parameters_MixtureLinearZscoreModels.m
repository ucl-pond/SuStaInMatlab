function [S_opt,f_opt,likelihood_opt] = optimise_parameters_MixtureLinearZscoreModels(data,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
    stage_zscore,stage_biomarker_index,S_init,f_init,likelihood_flag)
% Optimise the parameters of the SuStaIn model

M = size(data,1);
N_S = size(S_init,1);
N = size(stage_zscore,2);

S_opt = S_init;
f_opt = f_init;

f_val_mat = repmat(f_opt,[1 N+1 M]);
f_val_mat = permute(f_val_mat,[3 2 1]);

p_perm_k = zeros(M,N+1,N_S);
for s = 1:N_S
    if (strcmp(likelihood_flag,'Exact'))
        p_perm_k(:,:,s) = ...
            calculate_likelihood_stage_LinearZscoreModel_Approx(data,...
            min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
            stage_zscore,stage_biomarker_index,S_opt(s,:));
    else
        p_perm_k(:,:,s) = ...
            calculate_likelihood_stage_LinearZscoreModel(data,...
            min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
            stage_zscore,stage_biomarker_index,S_opt(s,:));
    end
end

p_perm_k_weighted = p_perm_k.*f_val_mat;
p_perm_k_norm = p_perm_k_weighted./repmat(sum(sum(p_perm_k_weighted,2),3),[1 N+1 N_S]);
f_opt = squeeze(sum(sum(p_perm_k_norm)))/sum(sum(sum(p_perm_k_norm)));

f_val_mat = repmat(f_opt,[1 N+1 M]);
f_val_mat = permute(f_val_mat,[3 2 1]);

order_seq = randperm(N_S);

for s = order_seq
    % tic
    
    order_bio = randperm(N);
    
    for i = order_bio
        current_sequence = S_opt(s,:);
        current_location(current_sequence) = 1:N;
        
        selected_event = i;
        
        move_event_from = current_location(selected_event);
        
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
        possible_sequences = zeros(length(possible_positions),N);
        possible_likelihood = zeros(length(possible_positions),1);
        possible_p_perm_k = zeros(M,N+1,length(possible_positions));
        for index = 1:length(possible_positions)
            current_sequence = S_opt(s,:);
            
            move_event_to = possible_positions(index);
            
            current_sequence(move_event_from) = [];
            new_sequence = [current_sequence(1:move_event_to-1) selected_event current_sequence(move_event_to:N-1)];
            
            possible_sequences(index,:) = new_sequence;
            
            if (strcmp(likelihood_flag,'Exact'))
                possible_p_perm_k(:,:,index) = ...
                    calculate_likelihood_stage_LinearZscoreModel_Approx(data,...
                    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
                    stage_zscore,stage_biomarker_index,new_sequence);
            else
                possible_p_perm_k(:,:,index) = ...
                    calculate_likelihood_stage_LinearZscoreModel(data,...
                    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
                    stage_zscore,stage_biomarker_index,new_sequence);
            end
            
            p_perm_k(:,:,s) = possible_p_perm_k(:,:,index);
            total_prob_stage = sum(p_perm_k.*f_val_mat,3);
            total_prob_subj = sum(total_prob_stage,2);
            possible_likelihood(index) = sum(log(total_prob_subj+1e-250));
        end
        
        max_likelihood = max(possible_likelihood);
        this_S = possible_sequences(possible_likelihood==max_likelihood,:);
        this_S = this_S(1,:);
        S_opt(s,:) = this_S;
        this_p_perm_k = possible_p_perm_k(:,:,possible_likelihood==max_likelihood);
        p_perm_k(:,:,s) = this_p_perm_k(:,:,1);
    end
    % toc
    S_opt(s,:) = this_S;
end

p_perm_k_weighted = p_perm_k.*f_val_mat;
p_perm_k_norm = p_perm_k_weighted./repmat(sum(sum(p_perm_k_weighted,2),3),[1 N+1 N_S]);
f_opt = squeeze(sum(sum(p_perm_k_norm)))/sum(sum(sum(p_perm_k_norm)));

f_val_mat = repmat(f_opt,[1 N+1 M]);
f_val_mat = permute(f_val_mat,[3 2 1]);

total_prob_stage = sum(p_perm_k.*f_val_mat,3);
total_prob_subj = sum(total_prob_stage,2);

likelihood_opt = sum(log(total_prob_subj+1e-250));


end

