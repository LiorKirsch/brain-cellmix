function compare_deconv_to_okaty_cell_type(cell_mix, gene_info, species,limit_to_cortical_cell_types)


    addpath('markers profiles/');
    mouse_cell_types = load('mouse_cell_type_profiles.mat');
    % [found_translation, entrez_translation] = translateAffyToEntrez(mouse_cell_types.gene_affi_id, 'entrez', false,gene_info.entrez_ids);
    % mouse_cell_types.gene_entrez = entrez_translation;
    % mouse_cell_types.gene_symbol = mouse_cell_types.gene_symbol(found_translation);
    % 
    % [found_translation, symbol_translation] = translateAffyToEntrez(mouse_cell_types.gene_affi_id, 'symbol', false,gene_info.gene_symbols);
    % mouse_cell_types.gene_symbol = symbol_translation;
    % 
    % mouse_cell_types.expression = mouse_cell_types.expression(found_translation,:);
    % mouse_cell_types.gene_affi_id = mouse_cell_types.gene_affi_id(found_translation);

    [mouse_cell_types,sep_ind_neuron,sep_ind_astro,~] = order_sample_by_type(mouse_cell_types);

    if exist('limit_to_cortical_cell_types','var');
        if (limit_to_cortical_cell_types)
            mouse_cell_types = limit_to_samples_from_cortex(mouse_cell_types);
        end
    end

    switch species
        case 'mouse'
            [reorder_mouse_cell_type, reorder_zapala] = reorderUsingId(mouse_cell_types.all_symbols, gene_info.gene_symbols);
            reorder_mouse_cell_type = mouse_cell_types.refer_to_index(reorder_mouse_cell_type);
            reorder_mouse_cell_type_expresion = mouse_cell_types.expression(reorder_mouse_cell_type,:);
            reorder_cellmix_expresion = cell_mix.celltype_profile(reorder_zapala,:);
        case {'monkey', 'human'}
            addpath('/cortex/code/matlab/homologous_gene_mapping/');
            switch species
                case 'monkey'
                    [gene_to_group_mouse, gene_to_group_primate, homologous_group_id] =  gene_to_homolog_group('mouse_laboratory','rhesus_macaque', mouse_cell_types.gene_symbol, 'symbol',gene_info.entrez_ids,'entrez_gene_ID');
                case 'human'
                    [gene_to_group_mouse, gene_to_group_primate, homologous_group_id] =  gene_to_homolog_group('mouse_laboratory','human', mouse_cell_types.gene_symbol, 'symbol',gene_info.entrez_ids,'entrez_gene_ID');
            end
            
            groups_with_1_to_1 = sum(gene_to_group_mouse,1) == 1  & sum(gene_to_group_primate,1) == 1;
            gene_to_group_mouse = gene_to_group_mouse(:,groups_with_1_to_1);
            gene_to_group_primate = gene_to_group_primate(:,groups_with_1_to_1);
            gene_to_group_mouse = (1:size(gene_to_group_mouse,1)) * gene_to_group_mouse ;
            gene_to_group_primate = (1:size(gene_to_group_primate,1)) * gene_to_group_primate ;

            reorder_mouse_cell_type_expresion = mouse_cell_types.expression(gene_to_group_mouse,:);
            reorder_cellmix_expresion = cell_mix.celltype_profile(gene_to_group_primate,:);
    end

    corr_matrix = corr(reorder_mouse_cell_type_expresion, reorder_cellmix_expresion, 'type','spearman');
    draw_imagesc_with_seperation(corr_matrix,sep_ind_neuron,sep_ind_astro);

    ax = gca;
    ax.XTick = 1:length(cell_mix.cell_types);
    ax.XTickLabel = cell_mix.cell_types;
    % ax.YTick = 1:length(mouse_cell_types.cellTypesDescription);
    % ax.YTickLabel = mouse_cell_types.cellTypesDescription;

    sample_cell_type_id =double(mouse_cell_types.sample2type) * ( ( 1:size(mouse_cell_types.sample2type,2))');
    mouse_cell_types.is_neuron = mouse_cell_types.is_neuron(sample_cell_type_id);
    mouse_cell_types.is_astro = mouse_cell_types.is_astro(sample_cell_type_id);
    mouse_cell_types.is_oligo = mouse_cell_types.is_oligo(sample_cell_type_id);
    drawROC(cell_mix, mouse_cell_types, corr_matrix);

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
    
    y = [sep_index2 sep_index2];
    plot(x,y,'Color','w','LineStyle','-');
    plot(x,y,'Color','k','LineStyle',':');
    hold off;
    
    axis ij;
    ylim([1 size(matrix_to_draw,1)]);
    
%      matrix_to_draw2 = [matrix_to_draw(1:sep_index1,:); -10*ones(3,3); matrix_to_draw( (sep_index1+1):(sep_index2 ) ,:); -10*ones(3,3); matrix_to_draw((sep_index2+1): end,:) ];
%      imagesc(matrix_to_draw2);
%      
end


function [reorder_A, reorder_B] = reorderUsingId(IdA, IdB)
    intersection_ids = intersect(IdA, IdB);
    [~, reorder_A] = ismember(intersection_ids, IdA);
    [~, reorder_B] = ismember(intersection_ids, IdB);
end