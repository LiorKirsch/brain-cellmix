function check_corrleation_within_class()


    addpath('~/Projects/general use functions/');
    addpath('markers profiles/');
    mouse_cell_types = load('mouse_cell_type_profiles.mat');
   
    [mouse_cell_types,sep_ind_neuron,sep_ind_astro,~] = order_sample_by_type(mouse_cell_types);
    
    corr_cell_types = corr(mouse_cell_types.expression);

    draw_imagesc_with_seperation(corr_cell_types, sep_ind_neuron, sep_ind_astro);
    title('cell type correlation');
    
    
    is_neuron = logical(double(mouse_cell_types.sample2type) * double(mouse_cell_types.is_neuron));
    is_astro = logical(double(mouse_cell_types.sample2type) * double(mouse_cell_types.is_astro));
    is_oligo = logical(double(mouse_cell_types.sample2type) * double(mouse_cell_types.is_oligo));
    disp('AUC for cell types');
    getMeanAUC(corr_cell_types, is_neuron, is_astro, is_oligo);

    
    mouse_cell_types = limit_to_samples_from_cortex(mouse_cell_types);
    [mouse_cell_types,sep_ind_neuron,sep_ind_astro,~] = order_sample_by_type(mouse_cell_types);
    
    corr_cell_types = corr(mouse_cell_types.expression);

    draw_imagesc_with_seperation(corr_cell_types, sep_ind_neuron, sep_ind_astro);
    title('cell type correlation - limited to cortex');
    
    
    is_neuron = logical(double(mouse_cell_types.sample2type) * double(mouse_cell_types.is_neuron));
    is_astro = logical(double(mouse_cell_types.sample2type) * double(mouse_cell_types.is_astro));
    is_oligo = logical(double(mouse_cell_types.sample2type) * double(mouse_cell_types.is_oligo));
    disp('AUC for cell types - limited to cortex');
    getMeanAUC(corr_cell_types, is_neuron, is_astro, is_oligo);
end

function getMeanAUC(corr_cell_types, is_neuron, is_astro, is_oligo)

    neuron_ind = find(is_neuron);
    astro_ind = find(is_astro);
    oligo_ind = find(is_oligo);
    
    neuron_auc = nan(length(neuron_ind),1);
    astro_auc = nan(length(astro_ind),1);
    oligo_auc = nan(length(oligo_ind),1);
    
    for i =1:length(neuron_ind)
        current_ind = neuron_ind(i);
        current_corr = corr_cell_types(:,current_ind);
        current_is_neuron = is_neuron;
        
        current_corr(current_ind) = [];
        current_is_neuron(current_ind) = [];
        neuron_auc(i) = auc_uri(current_is_neuron,current_corr,false);
    end
    
    for i =1:length(astro_ind)
        current_ind = astro_ind(i);
        current_corr = corr_cell_types(:,current_ind);
        current_is_astro = is_astro;
        
        current_corr(current_ind) = [];
        current_is_astro(current_ind) = [];
        astro_auc(i) = auc_uri(current_is_astro,current_corr,false);
    end
    
    for i =1:length(oligo_ind)
        current_ind = oligo_ind(i);
        current_corr = corr_cell_types(:,current_ind);
        current_is_oligo = is_oligo;
         
        current_corr(current_ind) = [];
        current_is_oligo(current_ind) = [];
        oligo_auc(i) = auc_uri(current_is_oligo,current_corr,false);
    end
       
       
    fprintf('neurons auc %g (+/- %g)\n', mean(neuron_auc), std(neuron_auc));
    fprintf('astro auc %g (+/- %g)\n', mean(astro_auc), std(astro_auc));
    fprintf('oligo auc %g (+/- %g)\n', mean(oligo_auc), std(oligo_auc));

        
        
        
        
        
        
end


function draw_imagesc_with_seperation(matrix_to_draw, sep_index1, sep_index2)
     
    figure; hold on;
     imagesc(matrix_to_draw);  
     colormap(jet); 
     colorbar;
     
    x = [0.5 ,size(matrix_to_draw,2)+0.5];
    
    y = [sep_index1 ,sep_index1];
    plot(x,y,'Color','w','LineStyle','-');
    plot(x,y,'Color','k','LineStyle',':');
    
    plot(y,x,'Color','w','LineStyle','-');
    plot(y,x,'Color','k','LineStyle',':');
    
    y = [sep_index2 sep_index2];
    plot(x,y,'Color','w','LineStyle','-');
    plot(x,y,'Color','k','LineStyle',':');
    
    plot(y,x,'Color','w','LineStyle','-');
    plot(y,x,'Color','k','LineStyle',':');
   
    hold off;
    
    axis ij;
    ylim([1 size(matrix_to_draw,1)]);
    
%      matrix_to_draw2 = [matrix_to_draw(1:sep_index1,:); -10*ones(3,3); matrix_to_draw( (sep_index1+1):(sep_index2 ) ,:); -10*ones(3,3); matrix_to_draw((sep_index2+1): end,:) ];
%      imagesc(matrix_to_draw2);
%      
end