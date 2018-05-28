function [] = plot_SuStaIn_model(samples_sequence,samples_f,...
    biomarker_labels,stage_zscore,stage_biomarker_index,...
    colour_mat,plot_order)

temp_mean_f = mean(samples_f,2);
[vals,ix] = sort(temp_mean_f,'descend');
vals = round(vals*100)/100;

N_S = size(samples_sequence,1);

if (N_S<=3)
    figure('Position',[0 0 400*N_S 300]);
else
    figure('Position',[0 0 400*ceil(N_S/2) 350*floor(N_S/2)]);
end
for i = 1:N_S
    this_samples_sequence = squeeze(samples_sequence(ix(i),:,:))';
    if (N_S<=3)
        subplot(1,N_S,i);
    else
        subplot(floor(N_S/2),ceil(N_S/2),i);
    end
    
    markers = unique(stage_biomarker_index);
    
    N_z = size(colour_mat,1);
    
    N = size(this_samples_sequence,2);
    
    confus_matrix = zeros(N,N);
    for j = 1:N
        confus_matrix(j,:) = sum(this_samples_sequence==j);
    end
    
    confus_matrix = confus_matrix/max(size(this_samples_sequence));
    
    zvalues = unique(stage_zscore);
    
    N_bio = size(markers,2);
    confus_matrix_z = zeros(N_bio,N,N_z);
    for z = 1:N_z
        confus_matrix_z(stage_biomarker_index(stage_zscore==zvalues(z)),:,z) = confus_matrix(stage_zscore==zvalues(z),:);
    end
    
    confus_matrix_c = ones(N_bio,N,3);
    for z = 1:N_z
        this_confus_matrix = confus_matrix_z(:,:,z);
        this_colour = colour_mat(z,:);
        alter_level = this_colour==0;
        this_colour_matrix = zeros(N_bio,N,3);
        this_colour_matrix(:,:,alter_level) = repmat(this_confus_matrix(markers,:),[1 1 sum(alter_level)]);
        confus_matrix_c = confus_matrix_c-this_colour_matrix;
    end
    
    image(confus_matrix_c(plot_order,:,:));
    set(gca,'XTick',1:N);
    set(gca, 'XTickLabel',(1:N));
    set(gca,'YTick',1:N_bio);
    set(gca, 'YTickLabel',biomarker_labels(plot_order));
    set(gca,'FontSize',12);
    xlabel('Event Position');
    
    title(strcat('Group ',num2str(i),', f=',num2str(vals(i))))
end



end

