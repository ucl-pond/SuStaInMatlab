function [data,data_denoised,stage_value] = generate_data_SuStaIn(subtypes,stages,gt_ordering,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,stage_zscore,stage_biomarker_index)

N = size(stage_biomarker_index,2);

N_S = size(gt_ordering,1);

possible_biomarkers = unique(stage_biomarker_index);
B = size(possible_biomarkers,2);

stage_value = zeros(B,N+2,N_S);
for s = 1:N_S
    S = gt_ordering(s,:);
    S_inv(S) = 1:N;
    for i = 1:B
        b = possible_biomarkers(i);
        event_location = [0 S_inv(stage_biomarker_index==b) N+1];
        event_value = [min_biomarker_zscore(i) stage_zscore(stage_biomarker_index==b) max_biomarker_zscore(i)];
        for j = 1:length(event_location)-1
            index = event_location(j)+1:event_location(j+1)+1;
            stage_value(i,index,s) = linspace(event_value(j),event_value(j+1),length(index)); 
        end
    end
end

M = size(stages,1);
data_denoised = zeros(M,B);
for m = 1:M
    data_denoised(m,:) = stage_value(:,stages(m)+1,subtypes(m));
end

data = data_denoised + randn(M,B).*repmat(std_biomarker_zscore,[M 1]);

end

