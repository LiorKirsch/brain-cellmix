function aucs = drawROC(cell_mix, mouse_cell_types, corr_matrix)

    addpath('/home/lab/lior/Projects/general use functions/')
    figure('Name','roc curves');
    hold on;
    
    num_cell_types = length(cell_mix.cell_types);
    aucs = nan(length(mouse_cell_types),1);
    for i =1:num_cell_types
        switch cell_mix.cell_types{i}
            case 'Astrocytes'
                labels = mouse_cell_types.is_astro;
            case 'Neurons'
                labels = mouse_cell_types.is_neuron;
            case 'Oligodendrocytes'
                labels = mouse_cell_types.is_oligo;
            otherwise
                labels = nan;
        end
        if all(labels==0)
            labels = nan;
        end
        fprintf('roc for %s\n',  cell_mix.cell_types{i} );
        scores = corr_matrix(:,i);
        if isnan(labels)
            aucs(i) = nan;
        else
            [auc,fpr,tpr] = auc_uri(labels,scores,true);
            legend_string{i} = sprintf('%s, auc=%g', cell_mix.cell_types{i}, auc);
            aucs(i) = auc;
        end
    end
    
    legend_string( isnan(aucs) ) = [];
    hold off;
    legend(legend_string ,'Location','southeast');
    legend('boxoff');

end