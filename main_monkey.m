function main_monkey()
% Do cellmix for macaque data (bernard's, blueprint_micro, blueprint_macro)
addpath('~/Projects/load datasets/');

%===  Bernard's data  ===
% dataset_name = 'bernardMonkey';
% dataset_name = 'bernardAllenDataMonkey';
dataset_name = 'blueprint_monkey_micro';
% dataset_name = 'blueprint_monkey_macro';
[expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info, ~] = load_expression_and_regions(dataset_name, 'gross regions');

[expression, gene_info] = limit_to_genes_with_entrez(expression, gene_info);

[expression,gross_region_vec, samples2subjects] = remove_sample_with_zero_expression(expression,gross_region_vec, samples2subjects);

% create_data_for_deconv(expression, gene_info,sprintf('deconv input/deconv_%s',dataset_name));

cellmix_file = sprintf('cellmix results/cellmix_%s.mat', dataset_name);
cell_mix = load(cellmix_file);
% cell_mix = load('cellmix_bernard_deconf.mat');


% drawHist(cell_mix, 'Neurons');
% drawHist(cell_mix, 'Astrocytes');
% drawHist(cell_mix, 'Oligodendrocytes');

only_v1_samples = find(strcmp('V1', gross_structures_info(:,1) ) );
only_v1_samples = gross_region_vec ==   only_v1_samples;
showCellTypeDistribution(cell_mix.proportions(: , only_v1_samples), cell_mix.cell_types);
title('V1 cell-type porportions');


addpath('markers profiles/');
mouse_cell_types = load('mouse_cell_type_profiles.mat');

% mouse_cell_types = limit_to_samples_from_cortex(mouse_cell_types);

compare_deconv_to_okaty_cell_type(cell_mix, gene_info, 'monkey');
title(regexprep(dataset_name,'_' , ' '));


end

function showCellTypeDistribution(proportions, cell_types)

    figure('Name','Cell type proportions');
    hold on;
    mean_cellmix = mean(proportions,2);
    std_cellmix = std(proportions,0,2);
    bar(1:length(cell_types), mean_cellmix);
    errorbar(1:length(cell_types), mean_cellmix, std_cellmix, '.' );
    ax = gca;
    ax.XTick = 1:length(cell_types);
    ax.XTickLabel = cell_types;
    hold off;
end
function drawHist(cell_mix, cell_type)
figure;
hist(cell_mix.proportions(ismember(cell_mix.cell_types,cell_type),:) );
xlabel( sprintf('distribution of proportions for %s', cell_type));
end


function mouse_cell_types = limit_data_by_cell_type_filter(mouse_cell_types, cell_type_filter)
    mouse_cell_types.cell_type_description = mouse_cell_types.cell_type_description(cell_type_filter);
    mouse_cell_types.anatomical_region = mouse_cell_types.anatomical_region(cell_type_filter);
    mouse_cell_types.age_postnatal_day = mouse_cell_types.age_postnatal_day(cell_type_filter);
    mouse_cell_types.method_for_single_cell = mouse_cell_types.method_for_single_cell(cell_type_filter);
    mouse_cell_types.RNA_isolation_method = mouse_cell_types.RNA_isolation_method(cell_type_filter);
    mouse_cell_types.RNA_amplification_and_labeling_method = mouse_cell_types.RNA_amplification_and_labeling_method(cell_type_filter);
    mouse_cell_types.RNA_input_amount_to_microarray = mouse_cell_types.RNA_input_amount_to_microarray(cell_type_filter);
    mouse_cell_types.microarray_platform = mouse_cell_types.microarray_platform(cell_type_filter);
    mouse_cell_types.reference = mouse_cell_types.reference(cell_type_filter);
    mouse_cell_types.cell_type_id = mouse_cell_types.cell_type_id(cell_type_filter);
    
    mouse_cell_types.is_neuron = mouse_cell_types.is_neuron(cell_type_filter);
    mouse_cell_types.is_astro = mouse_cell_types.is_astro(cell_type_filter);
    mouse_cell_types.is_oligo = mouse_cell_types.is_oligo(cell_type_filter);
    mouse_cell_types.is_cortex_or_hippocampus = mouse_cell_types.is_cortex_or_hippocampus(cell_type_filter);
    
    mouse_cell_types.sample2type = mouse_cell_types.sample2type(:, cell_type_filter);
end
function mouse_cell_types = limit_data_by_sample_filter(mouse_cell_types, sample_filter)    
    mouse_cell_types.expression = mouse_cell_types.expression(:,sample_filter);
    mouse_cell_types.sample2type = mouse_cell_types.sample2type(sample_filter, :);
    mouse_cell_types.samples_id = mouse_cell_types.samples_id(sample_filter);
end

function [expression,gross_region_vec, samples2subjects] = remove_sample_with_zero_expression(expression,gross_region_vec, samples2subjects)

    samples_with_zero_exp = all(expression == 0,2);
    
    if any(samples_with_zero_exp)
        fprintf('removing sample with zero expression:');
        disp(find(samples_with_zero_exp));
        expression = expression(~samples_with_zero_exp,:);
        gross_region_vec = gross_region_vec(~samples_with_zero_exp,:);
        samples2subjects = samples2subjects(~samples_with_zero_exp,:);
    end
    
end
function cell_type_mix_corr()
 
cell_mix_nmf = load('output_cell_mix.mat');
cell_mix_deconf = load('output_cell_mix_deconf.mat');

% %=================== correlations across genes ==================
% corr_matirx = corr(cell_mix_nmf.celltype_profile, cell_mix_deconf.celltype_profile, 'type','spearman');

% %========== correlations proportions in brain samples ===========
corr_matirx = corr(cell_mix_nmf.proportions', cell_mix_deconf.proportions', 'type','spearman');

figure; imagesc(corr_matirx); 
colormap(jet); colorbar;
ax = gca;
ax.XTick = 1:length(cell_mix_deconf.cell_types);
ax.YTick = 1:length(cell_mix_nmf.cell_types);
ax.XTickLabel = cell_mix_deconf.cell_types;
ax.YTickLabel = cell_mix_nmf.cell_types;

end

function [expression, gene_info] = limit_to_genes_with_entrez(expression, gene_info)

%=== limit to probes which have entrez ID ===
has_entrez = ~strcmp(gene_info.entrez_ids ,'');
expression = expression(:, has_entrez);
gene_info.entrez_ids = gene_info.entrez_ids(has_entrez);
gene_info.gene_full_name = gene_info.gene_full_name(has_entrez);
gene_info.gene_symbols = gene_info.gene_symbols(has_entrez);
gene_info.probe_id = gene_info.probe_id(has_entrez);
end

function create_data_for_deconv(expression, gene_info, matfile_name)
% ===== QUESTION ==========
% should I use all probes or just those with a gene symbol ?

    gene_symbols = gene_info.gene_symbols;

%     % make the minimum value zero and not -458
%     expression = expression - min(min(expression));
%     

    mapping = load('/cortex/data/microarray/primate/Bernard 2012 GSE31613/human2MacaqueGeneMapping.mat');
    monkey_markers = {'SSBP2','CELF4','SYT1','GFAP','AQP4','FGFR3','SLC1A2','GJB6','LOC720908','SOX10','MAG','MOG'};
    human_markers = {'STMN2','(CELF4)','SYT1','GFAP','AQP4','FGFR3','SLC1A2','GJB6','(MBP)','SOX10','MAG','MOG'};
    valid_markers = ismember(human_markers, gene_symbols);
    disp( human_markers( ~valid_markers ) );
    

    expression = expression';
    save(sprintf('%s.mat',matfile_name),'expression','gene_symbols');

end

