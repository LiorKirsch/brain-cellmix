load('/cortex/data/microarray/human/Giger2010/GSE12293-GPL570-Giger-2008.mat');
expression = gene_expression_data;
neuron_samples_mask = ismember(sample_sources, {'neurons from dorsolateral prefrontal cortex'});
expression = expression(:,neuron_samples_mask);
num_samples = sum(neuron_samples_mask);
sample2type = ones(num_samples,1);
cell_type_id = {'prefrontal cortex neurons'};
cell_type_description = {'neurons from dorsolateral prefrontal cortex'};
anatomical_region = {'dorsolateral prefrontal cortex'};
method_for_single_cell = {'laser capture microdissection'};
microarray_platform = {'GPL570'};
reference = {'Giger 2010'};
is_oligo = false;
is_astro = false;
is_neuron = true;
is_cortex_or_hippocampus = true;
% gene_symbol = gene_symbol;
% entrez = entrez;
gene_affi_id = gpl_probe_id;

age_postnatal_day = nan;
method_for_single_cell = nan;
RNA_isolation_method = nan;
RNA_amplification_and_labeling_method = nan;
RNA_input_amount_to_microarray = nan;
cell_type_id = {'neuron'};
samples_id = 1:num_samples;

all_symbols = gene_symbol;
refer_to_index = 1:length(all_symbols);

save('giger_cell_profiles.mat','expression','cell_type_id','cell_type_id',...
'cell_type_description' ,'anatomical_region','method_for_single_cell','microarray_platform', 'reference',...
    'is_oligo','is_astro','is_neuron','is_cortex_or_hippocampus','gene_symbol','entrez',...
    'gene_affi_id','sample2type','age_postnatal_day','method_for_single_cell','RNA_isolation_method',...
    'RNA_amplification_and_labeling_method','RNA_input_amount_to_microarray','cell_type_id','samples_id',...
    'all_symbols','refer_to_index');