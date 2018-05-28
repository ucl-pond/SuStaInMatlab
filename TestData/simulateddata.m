addpath('SimulationFunctions')
addpath('../SuStaInFunctions')

N = 10;
M = 500;
N_S_gt = 1;
Z_vals = [1*ones(N,1) 2*ones(N,1) 3*ones(N,1)];
IX_vals = repmat(1:N,[3 1])';
Z_max = 5*ones(N,1);

stage_zscore = Z_vals(:)';
stage_biomarker_index = IX_vals(:)';

min_biomarker_zscore = zeros(1,N);
max_biomarker_zscore = Z_max';
std_biomarker_zscore = ones(1,N);

SuStaInLabels = cell(1,N);
SuStaInStageLabels = cell(size(stage_zscore));
for i = 1:N
    SuStaInLabels{i} = strcat('Biomarker ',num2str(i));
end
for i = 1:length(stage_zscore)
    SuStaInStageLabels{i} = strcat('B',num2str(stage_biomarker_index(i)),' - Z',num2str(stage_zscore(i)));
end

gt_f = 1 + 0.5*(0:N_S_gt-1);
gt_f = fliplr(gt_f/sum(gt_f));

gt_sequence = generate_random_SuStaIn_model(stage_zscore,stage_biomarker_index,N_S_gt);

N_k_gt = size(stage_zscore,2)+1;

subtypes = randsample(1:N_S_gt,M,true,gt_f)';
stages = ceil(rand(M,1)*(N_k_gt+1))-1;

[data,data_denoised,stage_value] = generate_data_SuStaIn(subtypes,stages,gt_sequence,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,stage_zscore,stage_biomarker_index);

save('simulateddataResults/simulateddata_SuStaIn_groundtruth.mat')

N_startpoints = 25;
N_S_max = 3;
N_iterations_MCMC = 1e6;
likelihood_flag = 'Exact';
output_folder = 'simulateddataResults';
dataset_name = 'simulateddata';
run_SuStaIn_algorithm(data,...
    min_biomarker_zscore,max_biomarker_zscore,std_biomarker_zscore,...
    stage_zscore,stage_biomarker_index,N_startpoints,N_S_max,N_iterations_MCMC,...
    likelihood_flag,output_folder,dataset_name);

