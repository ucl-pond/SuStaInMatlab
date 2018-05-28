function [p_perm_k] = calculate_likelihood_stage_LinearZscoreModel_Approx(data,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
    stage_zscore,stage_biomarker_index,S)
% Computes the likelihood of a single linear z-score model using an
% approximation method (faster)
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
% S - the current ordering of the z-score stages for a particular subtype
%   dim: 1 x number of z-score stages
%
%OUTPUTS:
% p_perm_k - the probability of each subjects data at each stage of a
% particular subtype in the SuStaIn model

N = size(stage_biomarker_index,2);

S_inv(S) = 1:N;

possible_biomarkers = unique(stage_biomarker_index);
B = size(possible_biomarkers,2);

point_value = zeros(B,N+2);
for i = 1:B
    b = possible_biomarkers(i);
    event_location = [0 S_inv(stage_biomarker_index==b) N+1];
    event_value = [min_biomarker_zscore(i) stage_zscore(stage_biomarker_index==b) max_biomarker_zscore(i)];
    for j = 1:(length(event_location)-1)
        point_value(i,1+(event_location(j):event_location(j+1))) = linspace(event_value(j),event_value(j+1),1+event_location(j+1)-event_location(j));
    end    
end
stage_value = 0.5*point_value(:,1:end-1)+0.5*point_value(:,2:end);

% Different method of same computation to check correct:
% stage_value_check = zeros(B,N+1);
% for i = 1:B
%     b = possible_biomarkers(i);
%     event_location = [0 S_inv(stage_biomarker_index==b) N+1];
%     event_value = [min_biomarker_zscore(i) stage_zscore(stage_biomarker_index==b) max_biomarker_zscore(i)];
%     
%     for j = 1:(length(event_location)-1)
%         this_y1 = event_value(j);
%         this_y2 = event_value(j+1);
%         this_x1 = event_location(j);
%         this_x2 = event_location(j+1);
%         this_gradient = (this_y2-this_y1)/(this_x2-this_x1);
%         this_minval = this_y1+this_gradient/2;
%         this_maxval = this_y2-this_gradient/2;
%         temp_values = this_minval:this_gradient:this_maxval;
%         if (abs(this_minval-this_maxval)<1e-6)
%             temp_values = this_minval;
%         end
%         temp_ix = (this_x1+1):(this_x2);
%         if (length(temp_ix)~=length(temp_values))
%             fprintf('Error!\n')
%         end
%         stage_value_check(i,temp_ix) = temp_values;
%     end
% end

M = size(data,1);
p_perm_k = zeros(M,N+1);
for j = 1:(N+1)
    p_perm_k(:,j) = 1/(N+1)*prod(normpdf(data,repmat(stage_value(:,j)',[M 1]),repmat(std_biomarker_zscore,[M 1])),2);
end

end

