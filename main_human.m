function main_human()
% Do cellmix for human data (Kang's, Allen 6 human)

% create_adult_only_mat();
% create_adult_cortex_only_mat();
% cell_type_mix_corr();



%============ Load dataset ==================
% dataset_name = 'kang adults';
% kang_data = load('data/kang_samples_adults.mat');
% cell_mix = load('cellmix_results/cell_mix_kang_nmf.mat');
% only_cortex_celltypes = false;


dataset_name = 'kang adults cortex';
kang_data = load('data/kang_cortex_samples_adults.mat');
cell_mix = load('cellmix results/cellmix_kang_cortex_nmf.mat');
only_cortex_celltypes = true;

% drawHist(cell_mix, 'Neurons');
% drawHist(cell_mix, 'Astrocytes');
% drawHist(cell_mix, 'Oligodendrocytes');

gene_info.gene_symbols = kang_data.gene_names;
gene_info.entrez_ids = arrayfun(@(x) sprintf('%d', x),kang_data.gene_entrez,'UniformOutput', false);

compare_deconv_to_okaty_cell_type(cell_mix, gene_info, 'human',only_cortex_celltypes);
title(dataset_name);

end

function drawHist(cell_mix, cell_type)
figure;
hist(cell_mix.proportions(ismember(cell_mix.cell_types,cell_type),:) );
xlabel( sprintf('distribution of proportions for %s', cell_type));

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

function create_adult_only_mat()

kang_data = load('/cortex/data/microarray/human/Kang2011/kang_samples_with_ontology.mat');
 only_adult_samples_mask = kang_data.samples2periods(:,11:end);
 only_adult_samples_mask = any(only_adult_samples_mask,2);
 
 kang_data.samples2periods(~only_adult_samples_mask,:) = [];
 kang_data.samples2regions(~only_adult_samples_mask,:) = [];
 kang_data.subject_id(~only_adult_samples_mask,:) = [];
 kang_data.hemisphere(~only_adult_samples_mask,:) = [];
 kang_data.data(:, ~only_adult_samples_mask) = []; 
 save('data/kang_samples_adults.mat','-struct','kang_data');
 
end




function create_adult_cortex_only_mat()

 kang_data = load('/cortex/data/microarray/human/Kang2011/kang_samples_with_ontology.mat');
 only_adult_samples_mask = kang_data.samples2periods(:,11:end);
 only_adult_samples_mask = any(only_adult_samples_mask,2);
 
 cortex_only_regions = {'primary auditory (A1) cortex','dorsolateral prefrontal cortex','posterior inferior parietal cortex','inferior temporal cortex','primary motor (M1) cortex','medial prefrontal cortex','orbital prefrontal cortex','primary somatosensory (S1) cortex','superior temporal cortex','primary visual (V1) cortex','ventrolateral prefrontal cortex'};
 cortex_regions = ismember(kang_data.region_names, cortex_only_regions);
 only_cortex_samples_mask = kang_data.samples2regions(:,cortex_regions);
 only_cortex_samples_mask = any(only_cortex_samples_mask,2);
 
 samples_mask = only_cortex_samples_mask & only_adult_samples_mask;
 
 kang_data.samples2periods(~samples_mask,:) = [];
 kang_data.samples2regions(~samples_mask,:) = [];
 kang_data.subject_id(~samples_mask,:) = [];
 kang_data.hemisphere(~samples_mask,:) = [];
 kang_data.data(:, ~samples_mask) = []; 
 save('data/kang_cortex_samples_adults.mat','-struct','kang_data');
 
end
