function [mouse_cell_types,ind_neuron,ind_astro,ind_oligo] = order_sample_by_type(mouse_cell_types,limit_to_cortical_cell_types)
    [~,cell_type_ind] = max(double([mouse_cell_types.is_neuron, mouse_cell_types.is_astro, mouse_cell_types.is_oligo]),[],2);
    [~,reorder_ind] = sort(cell_type_ind);
    
    [~,samples_cell_type_ind] = max(double([double(mouse_cell_types.sample2type) * double(mouse_cell_types.is_neuron), double(mouse_cell_types.sample2type) * double(mouse_cell_types.is_astro), double(mouse_cell_types.sample2type) * double(mouse_cell_types.is_oligo)]),[],2);
    [~,samples_reorder] = sort(samples_cell_type_ind);
%     samples_reorder = double(mouse_cell_types.sample2type) * reorder_ind;
%     [~, samples_reorder] = sort(samples_reorder);
    
    mouse_cell_types = limit_data_by_cell_type_filter(mouse_cell_types, reorder_ind);
    mouse_cell_types = limit_data_by_sample_filter(mouse_cell_types, samples_reorder);    
    
    if (limit_to_cortical_cell_types)
       mouse_cell_types = limit_to_samples_from_cortex(mouse_cell_types);
    end
        
    ind_neuron = max(find(double(mouse_cell_types.sample2type) * double(mouse_cell_types.is_neuron)));
    ind_astro = max(find(double(mouse_cell_types.sample2type) * double(mouse_cell_types.is_astro)));
    ind_oligo = max(find(double(mouse_cell_types.sample2type) * double(mouse_cell_types.is_oligo)));
    fprintf('neurons -> %d, astro -> %d, oligo -> %d\n',ind_neuron, ind_astro, ind_oligo);
end

