function main_mouse_markers_subsets()
addpath('/home/lab/lior/Projects/load datasets/');
% Do cellmix for mouse data (Zapala's, Akahoshi's, Allen energy)

markers_file = 'markers sets/brain_markers_mouse.txt';
%===  Zapala data  ===
dataset_name = 'zapala';
input_mat_file_name = 'deconv input/deconv_for_zapala.mat';
output_mat_file_name = sprintf('cellmix results/cellmix_%s_nmf.mat', dataset_name);
[expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info, ~] = load_expression_and_regions('zapalaMouse', []);
% create_data_for_deconv(expression, gene_info,input_mat_file_name);
check_different_markers_subsets(input_mat_file_name, output_mat_file_name, markers_file);


%===  Akahoshi data  ===
% dataset_name = 'akahoshi';
% [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info, ~] = load_expression_and_regions('akahoshiMouse', []);
% create_data_for_deconv(expression, gene_info,'deconv_input/deconv_for_akahoshi');


%===  Allen mouse cortex data  ===
% data = 'allen_cortex';
% [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info, ~] = load_expression_and_regions('mouse', 'cortex');
% create_data_for_deconv(expression, gene_info,'deconv_input/deconv_for_allen_cortex');


% load celltype expression from Okaty PLoS One

% do deconvolution using the cell marker with the two methods and compare
% the expression profile with those collect by Okaty.

cell_mix = load(output_mat_file_name);
% cell_mix = load(sprintf('cellmix_%s_Deconf.mat', dataset_name));


compare_deconv_to_okaty_cell_type(cell_mix, gene_info, 'mouse');
title(dataset_name);
% drawHist(cell_mix, 'Neurons');
% drawHist(cell_mix, 'Astrocytes');
% drawHist(cell_mix, 'Oligodendrocytes');

end

function create_data_for_deconv(expression, gene_info, matfile_name)

    gene_symbols = gene_info.gene_symbols;
    
    % replace {'Brunol4'} with Celf4
    brunol4_ind = find(strcmp(gene_symbols, 'Brunol4' ));
    gene_symbols(brunol4_ind) = 'Celf4';

    % make the minimum value zero and not -458
    expression = expression - min(min(expression));
    
    expression = expression';
    save(sprintf('%s.mat',matfile_name),'expression','gene_symbols');

end

function check_different_markers_subsets(input_mat_file_name, output_mat_file_name,marker_file_name)
    C = textscan(fopen(marker_file_name),'%q %q','HeaderLines',1,'Delimiter',' ');
    gene_names = C{1};
    marker_type = C{2};
    
    [~,file_name,ext] = fileparts(marker_file_name);
    
    file_name = sprintf('markers subsets/%s%s', file_name, ext);
    for i=1:length(gene_names)
        current_marker_file_name = [file_name, num2str(i)];
        current_output_file = fopen(current_marker_file_name,'w');
        current_gene_names = gene_names;
        current_marker_type = marker_type;
        current_gene_names(i) = [];
        current_marker_type(i) = [];
        
        cellfun(@(x,y) fprintf(current_output_file, '"%s" "%s"\n',x,y) , current_gene_names, current_marker_type)
        
        run_R_script('do_cellmix',input_mat_file_name, output_mat_file_name,current_marker_file_name);
        
    end
        
    

end