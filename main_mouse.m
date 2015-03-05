function main_mouse()
addpath('/home/lab/lior/Projects/load datasets/');
addpath('wet celltype compare/');
% Do cellmix for mouse data (Zapala's, Akahoshi's, Allen energy)

algorithem = 'DECONF'; %'NMF'; 'NMFforb'; 'meanProfile'; 'DSA'; 'DECONF';

% marker_file_name = 'markers sets/brain_markers_mouse.txt';
marker_file_name = 'markers sets/brain_markers_mouse_3types.txt';

%===  Zapala data  ===
% dataset_name = 'zapala';
% input_mat_file_name = 'deconv input/deconv_for_zapala.mat';
% [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info, ~] = load_expression_and_regions('zapalaMouse', []);
% gross_structures_info(strcmp('Bed nucleus of the stria terminalis',gross_structures_info))  = {'BNST'};
% % create_data_for_deconv(expression, gene_info,input_mat_file_name);

%==== Zapala cortex and hippo ===
dataset_name = 'zapala cortex';
regions_to_keep = {'Olfactory Bulbs';'Isocortex';'Motor Cortex';'Entorhinal Cortex';'Perirhinal Cortex';'Hippocampus';'Dentate Gyrus';'Hippocampus CA1';'Hippocampus CA3';};
input_mat_file_name = 'deconv input/deconv_for_zapala_cortex.mat';
[expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info, ~] = load_expression_and_regions('zapalaMouse', regions_to_keep);
create_data_for_deconv(expression, gene_info,input_mat_file_name);

%===  Akahoshi data  ===
% dataset_name = 'akahoshi';
% input_mat_file_name = 'deconv input/deconv_for_akahoshi.mat';
% [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info, ~] = load_expression_and_regions('akahoshiMouse', []);
% create_data_for_deconv(expression, gene_info,input_mat_file_name);


%===  Allen mouse cortex data  ===
% dataset_name = 'allen_cortex';
% input_mat_file_name = 'deconv input/deconv_for_allen_cortex.mat';
% [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info, ~] = load_expression_and_regions('mouse', 'cortex');
% % create_data_for_deconv(expression, gene_info,input_mat_file_name);

%===  Su cortex data  ===
% dataset_name = 'Su 2004';
% input_mat_file_name = 'deconv input/su2004.mat';
% [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info, ~] = load_expression_and_regions('suMouse', []);
% % marker_file_name = 'markers sets/brain_markers_for_su.txt';
% % create_data_for_deconv(expression, gene_info,input_mat_file_name);


check_marker_appear_in_data(marker_file_name, input_mat_file_name);

output_mat_file_name = sprintf('cellmix results/cellmix_%s_%s.mat', dataset_name,algorithem);

% load celltype expression from Okaty PLoS One

% do deconvolution using the cell marker with the two methods and compare
% the expression profile with those collect by Okaty.

run_R_script('do_cellmix',input_mat_file_name, output_mat_file_name,marker_file_name,algorithem);
cell_mix = load(output_mat_file_name);


[aucs,median_within, median_outside] = compare_deconv_to_okaty_cell_type(cell_mix, gene_info, 'mouse');
title(dataset_name);
% drawHist(cell_mix, 'Neurons');
% drawHist(cell_mix, 'Astrocytes');
% drawHist(cell_mix, 'Oligodendrocytes');

draw_proprtions_for_regions(cell_mix.proportions, cell_mix.cell_types, gross_structures_info, gross_region_vec); ylim([0,0.8]);
show_proportions(cell_mix.proportions, cell_mix.cell_types);
end


function drawHist(cell_mix, cell_type)
figure;
hist(cell_mix.proportions(ismember(cell_mix.cell_types,cell_type),:) );
xlabel( sprintf('distribution of proportions for %s', cell_type));

end

function cell_type_mix_corr(nmf_file, deconf_file)
 
cell_mix_nmf = load(nmf_file);
cell_mix_deconf = load(deconf_file);

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


function create_data_for_deconv(expression, gene_info, matfile_name)

    gene_symbols = gene_info.gene_symbols;
    
%     % replace {'Brunol4'} with Celf4
%     brunol4_ind = find(strcmp(gene_symbols, 'Brunol4' ));
%     gene_symbols{brunol4_ind} = 'Celf4';

    % replace {'Mtap2'} with Map2
    selected_ind = find(strcmp(gene_symbols, 'A730034C02 /// Mtap2' ));
    if ~isempty(selected_ind)
        gene_symbols{selected_ind} = 'Map2';
    end
    
    if any(any( expression< 0))
    % make the minimum value zero and not -458
    expression = expression - min(min(expression));
    end
    
    expression = expression';
    save(matfile_name,'expression','gene_symbols');

end