function draw_proprtions_for_regions(proportions, cell_types, gross_strcutures_info, gross_region_vec)
num_regions = length(gross_strcutures_info);
[num_types,num_samples] = size(proportions);
assert(length(cell_types) == num_types,'should be the same number of cell types');
assert(length(gross_region_vec) == num_samples,'should be the same number of cell types');

region_mean_val = nan(num_regions, num_types);
region_std_val = nan(num_regions, num_types);
for i = 1:num_regions
   selected_samples = gross_region_vec ==i;
   region_mean_val(i,:) = mean(proportions(:,selected_samples),2);
   region_std_val(i,:) = std(proportions(:,selected_samples),1,2);
end

figure;
hold on;
bar(region_mean_val);
errorbar(region_mean_val, region_std_val,'.');
hold off;
legend(cell_types);
ax = gca;
ax.XTick = 1:num_regions;
ax.XTickLabel = gross_strcutures_info;
ax.XTickLabelRotation = 45;
end