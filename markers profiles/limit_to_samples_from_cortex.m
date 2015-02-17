function mouse_cell_types = limit_to_samples_from_cortex(mouse_cell_types)
    cell_type_filter = mouse_cell_types.is_cortex_or_hippocampus;
    samples_filter = double(mouse_cell_types.sample2type) * double(cell_type_filter);
    samples_filter = logical(samples_filter);
    
    mouse_cell_types = limit_data_by_cell_type_filter(mouse_cell_types, cell_type_filter);
    mouse_cell_types = limit_data_by_sample_filter(mouse_cell_types, samples_filter);    
end