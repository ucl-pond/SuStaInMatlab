function [S] = initialise_sequence_LinearZscoreModel(stage_zscore,stage_biomarker_index)
% Randomly initialises a linear z-score model ensuring that the biomarkers
% are monotonically increasing
%
%INPUTS:
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
%
%OUTPUTS:
% S - a random linear z-score model under the condition that each biomarker
% is monotonically increasing

N = size(stage_zscore,2);

S = zeros(1,N);
for i = 1:N
    IS_min_zscore_event = false(1,N);
    possible_biomarkers = unique(stage_biomarker_index);
    for j = 1:length(possible_biomarkers)
        IS_unselected = false(1,N);
        IS_unselected(setdiff(1:N,S(1:i-1))) = true;
        this_biomarkers = ((stage_biomarker_index==possible_biomarkers(j))+(IS_unselected==1))==2;
        this_min_zscore_event = min(stage_zscore(this_biomarkers));
        if (~isempty(this_min_zscore_event))
            IS_min_zscore_event((this_biomarkers+(stage_zscore==this_min_zscore_event))==2) = true;
        end
    end
    events = 1:N;
    possible_events = events(IS_min_zscore_event);
    this_index = ceil(rand*length(possible_events));
    S(i) = possible_events(this_index);
end


end

