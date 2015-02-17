expression = csvread('samples_expression.csv');
z = textscan( fopen('cell_types_metadata.csv'), '%q %q %q %q %q %q %q %q %q %q %q %q %q %q','delimiter',',','HeaderLines',1);

cell_type_to_samples = z{1};
sample2type = [];
samples_id = {};
cell_type_id = {};
for i = 1:length(cell_type_to_samples)
   C = strsplit(cell_type_to_samples{i},':');
   cell_type_id{i} =C{1};
   samples_text = C{2};
   C = strsplit(samples_text,',');
   samples_id = [samples_id, C];
   sample2type = [sample2type; ones(length(C),1) * i ];
end
samples_id = samples_id';
samples_id = strtrim(samples_id);
cell_type_id = cell_type_id';

sample2type = sparse(1:length(sample2type),sample2type, ones(size(sample2type)), length(sample2type), length(cell_type_id));
sample2type = logical(full(sample2type));

cell_type_description = z{2};
anatomical_region = z{3};
age_postnatal_day = z{4};
method_for_single_cell = z{5};
RNA_isolation_method = z{6};
RNA_amplification_and_labeling_method = z{7};
RNA_input_amount_to_microarray = z{8};
microarray_platform = z{9};
reference = z{10};
is_oligo = strcmp(z{11},'TRUE');
is_astro = strcmp(z{12},'TRUE');
is_neuron = strcmp(z{13},'TRUE');
is_cortex_or_hippocampus = strcmp(z{14},'TRUE');

z = textscan( fopen('row_gene_metadata.csv'), '%q %q','delimiter',',','HeaderLines',1);
gene_symbol = z{1};
gene_affi_id = z{2};

z = textscan( fopen('column_samples_id.csv'), '%q','delimiter',',','HeaderLines',1);
samples_id_from_expression_file = z{1};

assert( all(strcmp(samples_id_from_expression_file, samples_id)), 'the order of the samples in the two files should be the same');
clear('z','C','samples_text','samples_id_from_expression_file','i','cell_type_to_samples');

all_symbols = [];
refer_to_index = [];
for i =1:length(gene_symbol)
    parsed_symbols = strsplit(gene_symbol{i},' /// ');
    all_symbols = [all_symbols; parsed_symbols'];
    refer_to_index = [refer_to_index; i*ones(length(parsed_symbols),1)];
end
clear('parsed_symbols');

save('mouse_cell_type_profiles.mat');