addpath('../SuStaInFunctions')

load('simulateddataResults/simulateddata_SuStaIn_groundtruth.mat')
load('simulateddataResults/simulateddata_SuStaIn_input.mat')

N = size(data,2);
SuStaInLabels = cell(1,N);
for i = 1:N
    SuStaInLabels{i} = strcat('Biomarker ',num2str(i));
end

plot_order = 1:size(SuStaInLabels,2);
colour_mat = [1 0 0; 1 0 1; 0 0 1];

h = figure;
for c = 1:2
    load(strcat('simulateddataResults/simulateddata_MCMC_',num2str(c),'Seq.mat'))
    figure(h); hold on; plot(1:N_iterations_MCMC,samples_likelihood)
    
    plot_SuStaIn_model(samples_sequence,samples_f,...
        SuStaInLabels,stage_zscore,stage_biomarker_index,...
        colour_mat,plot_order);
end

