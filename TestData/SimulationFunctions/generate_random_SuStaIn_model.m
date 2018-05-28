function [S] = generate_random_SuStaIn_model(stage_zscore,stage_biomarker_index,N_S)

N = size(stage_zscore,2);

S = zeros(N_S,N);
for s = 1:N_S
    for i = 1:N
        IS_min_stage_zscore = false(1,N);
        possible_biomarkers = unique(stage_biomarker_index);
        for j = 1:length(possible_biomarkers)
            IS_unselected = false(1,N);
            IS_unselected(setdiff(1:N,S(s,1:i-1))) = true;
            this_biomarkers = ((stage_biomarker_index==possible_biomarkers(j))+(IS_unselected==1))==2;
            this_min_stage_zscore = min(stage_zscore(this_biomarkers));
            if (~isempty(this_min_stage_zscore))
                IS_min_stage_zscore((this_biomarkers+(stage_zscore==this_min_stage_zscore))==2) = true;
            end
        end
        events = 1:N;
        possible_events = events(IS_min_stage_zscore);
        this_index = ceil(rand*length(possible_events));
        S(s,i) = possible_events(this_index);
    end
end


end

