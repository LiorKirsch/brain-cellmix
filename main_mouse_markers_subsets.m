function main_mouse_markers_subsets()
addpath('/home/lab/lior/Projects/load datasets/');
% Do cellmix for mouse data (Zapala's, Akahoshi's, Allen energy)

markers_file = 'markers sets/brain_markers_mouse.txt';
%===  Zapala data  ===
dataset_name = 'zapala';
algorithem = 'NMFforb'; %'NMF'; 'NMFforb'; 'meanProfile'; 'DSA'; 'DECONF';
input_mat_file_name = 'deconv input/deconv_for_zapala.mat';
output_mat_file_name = sprintf('cellmix results/cellmix_%s_%s.mat', dataset_name,algorithem);
[expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info, ~] = load_expression_and_regions('zapalaMouse', []);
% create_data_for_deconv(expression, gene_info,input_mat_file_name);


%===  Akahoshi data  ===
% dataset_name = 'akahoshi';
% [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info, ~] = load_expression_and_regions('akahoshiMouse', []);
% create_data_for_deconv(expression, gene_info,'deconv_input/deconv_for_akahoshi');


%===  Allen mouse cortex data  ===
% data = 'allen_cortex';
% [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info, ~] = load_expression_and_regions('mouse', 'cortex');
% create_data_for_deconv(expression, gene_info,'deconv_input/deconv_for_allen_cortex');


[output_files,gene_names] = check_different_markers_subsets(input_mat_file_name, output_mat_file_name, markers_file, algorithem, true);
num_classes = 4;
num_classes_with_labels= 3;

aucs = nan(length(gene_names),num_classes);
median_within = nan(length(gene_names),num_classes_with_labels);
median_outside = nan(length(gene_names),num_classes_with_labels);

for i=1:length(output_files);
    try
    cell_mix = load(output_files{i});
    [aucs(i,:),median_within(i,:),median_outside(i,:)] = compare_deconv_to_okaty_cell_type(cell_mix, gene_info, 'mouse');
    end
    title(sprintf('%s - without %s',dataset_name,gene_names{i}));
end


cell_mix = load(output_mat_file_name);
[aucs_all,median_within_all,median_outside_all] = compare_deconv_to_okaty_cell_type(cell_mix, gene_info, 'mouse');
aucs = [aucs_all; aucs];
median_within = [median_within_all; median_within];
median_outside = [median_outside_all; median_outside];
gene_names = ['no marker missing (NMF)' ;gene_names];
title('all markers (NMF)');

output_mat_file_name = sprintf('cellmix results/cellmix_%s_%s.mat', dataset_name,'DECONF');
cell_mix = load(output_mat_file_name);
[aucs_all,median_within_all,median_outside_all] = compare_deconv_to_okaty_cell_type(cell_mix, gene_info, 'mouse');
aucs = [aucs_all; aucs];
median_within = [median_within_all; median_within];
median_outside = [median_outside_all; median_outside];
gene_names = ['no marker missing (DECONF)' ;gene_names];
title('all markers (DECONF)');

do_bar_fig(aucs,gene_names,cell_mix.cell_types);
title('auc without the marker');

do_bar_fig(median_within,gene_names,cell_mix.cell_types([1,3,4]));
title('median correlation within class');

do_bar_fig(median_outside,gene_names,cell_mix.cell_types([1,3,4]));
title('median correlation with samples from a different class');


end

function do_bar_fig(results_for_bar,gene_names,legendStrings)
figure;
bar(results_for_bar);
legend(legendStrings)
ax = gca;
ax.XTick = 1:length(gene_names);
ax.XTickLabel = gene_names;
ax.XTickLabelRotation = 45;
title('auc without the marker');
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

function [output_files,gene_names] = check_different_markers_subsets(input_mat_file_name, output_mat_file_name,marker_file_name,algorithem,skip_existing)
    C = textscan(fopen(marker_file_name),'%q %q','HeaderLines',1,'Delimiter',' ');
    gene_names = C{1};
    marker_type = C{2};
    
    [~,marker_file_name,ext] = fileparts(marker_file_name);
    marker_file_name = sprintf('markers subsets/%s%s', marker_file_name, ext);
    
    [~,output_file_name,ext] = fileparts(output_mat_file_name);
    output_file_name = sprintf('markers subsets/%s', output_file_name);
    output_files = {};
    for i=1:length(gene_names)
        fprintf('======= (%d/%d) Checking without %s ========\n',i,length(gene_names), gene_names{i} );
        current_marker_file_name = [marker_file_name, num2str(i)];
        current_output_file_name = sprintf('%s%d%s',output_file_name, i,ext);
        current_marker_file = fopen(current_marker_file_name,'w');
        
        current_gene_names = gene_names;
        current_marker_type = marker_type;
        current_gene_names(i) = [];
        current_marker_type(i) = [];
        
        fprintf(current_marker_file,'"GENE NAME" "Class name"\n');
        cellfun(@(x,y) fprintf(current_marker_file, '"%s" "%s"\n',x,y) , current_gene_names, current_marker_type);
        
        if (exist(current_output_file_name, 'file') & skip_existing)
            fprintf('found existing file, skipping ...\n');
        else
            run_R_script('do_cellmix',input_mat_file_name, current_output_file_name,current_marker_file_name,algorithem);
        end
        output_files{i} = current_output_file_name;
    end
        
    

end