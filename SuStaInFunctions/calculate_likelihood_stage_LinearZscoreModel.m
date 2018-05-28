function [p_perm_k] = calculate_likelihood_stage_LinearZscoreModel(data,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
    stage_zscore,stage_biomarker_index,S)
% Computes the likelihood of a single linear z-score model using an exact
% method (slower)
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

tau_val = linspace(0,1,N+2);

point_value = zeros(B,N+2);
for i = 1:B
    b = possible_biomarkers(i);
    event_location = [0 S_inv(stage_biomarker_index==b) N+1];
    event_value = [min_biomarker_zscore(i) stage_zscore(stage_biomarker_index==b) max_biomarker_zscore(i)];
    for j = 1:(length(event_location)-1)
        point_value(i,1+(event_location(j):event_location(j+1))) = linspace(event_value(j),event_value(j+1),1+event_location(j+1)-event_location(j));
    end    
end

stage_initial_value = point_value(:,1:end-1);
stage_final_value = point_value(:,2:end);

stage_initial_tau = tau_val(1:end-1);
stage_final_tau = tau_val(2:end);

stage_a = (stage_final_value-stage_initial_value)./repmat(stage_final_tau-stage_initial_tau,[B 1]);
stage_b = stage_initial_value - stage_a.*repmat(stage_initial_tau,[B 1]);
stage_std = repmat(std_biomarker_zscore',[1 N+1]);

M = size(data,1);

iterative_mean = (repmat(data(:,1),[1 N+1])-repmat(stage_b(1,:),[M 1]))./repmat(stage_a(1,:),[M 1]);
iterative_std = repmat(stage_std(1,:),[M 1])./repmat(stage_a(1,:),[M 1]);
iterative_kappa = ones(M,N+1);
for b = 2:B
    mu1 = iterative_mean;
    mu2 = (repmat(data(:,b),[1 N+1])-repmat(stage_b(b,:),[M 1]))./repmat(stage_a(b,:),[M 1]);
    std1 = iterative_std;
    std2 = repmat(stage_std(b,:),[M 1])./repmat(stage_a(b,:),[M 1]);
    cov1 = std1.^2;
    cov2 = std2.^2;
    munew = (cov1.^-1+cov2.^-1).^-1.*(cov1.^-1.*mu1+cov2.^-1.*mu2);
    covnew = (cov1.^-1+cov2.^-1).^-1;
    kappaval = normpdf(mu1,mu2,sqrt(cov1+cov2));
    iterative_mean = munew;
    iterative_std = sqrt(covnew);
    iterative_kappa = iterative_kappa.*kappaval;
end

iterative_const = repmat(prod(1./abs(stage_a),1),[M 1]);

cdf_diff_val = normcdf(repmat(stage_final_tau,[M 1]),iterative_mean,iterative_std) - ...
    normcdf(repmat(stage_initial_tau,[M 1]),iterative_mean,iterative_std);

if (size(iterative_const,1)~=size(iterative_kappa,1))
    fprintf('Error!\n')
end
if (size(iterative_const,1)~=size(cdf_diff_val,1))
    fprintf('Error!\n')
end
if (size(iterative_const,2)~=size(iterative_kappa,2))
    fprintf('Error!\n')
end
if (size(iterative_const,2)~=size(cdf_diff_val,2))
    fprintf('Error!\n')
end

p_perm_k = iterative_const.*iterative_kappa.*cdf_diff_val;

end

