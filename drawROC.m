function drawROC(cell_mix, mouse_cell_types, corr_matrix)

    addpath('/home/lab/lior/Projects/general use functions/')
    figure('Name','roc curves');
    hold on;
    
    
    fprintf('roc for astrocytes (%s)\n',  cell_mix.cell_types{1} );
    labels = mouse_cell_types.is_astro;
    scores = corr_matrix(:,1);
    [auc,fpr,tpr] = auc_uri(labels,scores,true);
    legend_string{1} = sprintf('%s, auc=%g', cell_mix.cell_types{1}, auc);
%     drawRocCurveWithStd({tpr} , {fpr} ,auc);
    
    fprintf('roc for neruons (%s)\n',  cell_mix.cell_types{2} );
    labels = mouse_cell_types.is_neuron;
    scores = corr_matrix(:,2);
    [auc,fpr,tpr] = auc_uri(labels,scores,true);
    legend_string{2} = sprintf('%s, auc=%g', cell_mix.cell_types{2}, auc);
    
    fprintf('roc for Oligodendrocytes (%s)\n',  cell_mix.cell_types{3} );
    labels = mouse_cell_types.is_oligo;
    scores = corr_matrix(:,3);
    [auc,fpr,tpr] = auc_uri(labels,scores,true);
    legend_string{3} = sprintf('%s, auc=%g', cell_mix.cell_types{3}, auc);
    
    hold off;
    legend(legend_string ,'Location','southeast');
    legend('boxoff');


end