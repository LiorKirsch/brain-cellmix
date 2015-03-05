function check_marker_appear_in_data(marker_file_name, data_file_name)

[gene_name, gene_group] = textread(marker_file_name,'%q %q','headerlines',1);
matVar = whos(matfile(data_file_name));
for i=1:length(matVar)
    assert( ismember(matVar(i).name ,{'gene_symbols','expression'}),'mat file should contain an "expression" and "gene_symbols" variables');
end
load(data_file_name, 'gene_symbols');

is_member = ismember(gene_name, gene_symbols);
if ~all(is_member)
    disp('these genes do not appear in data');
    disp(gene_name( ~is_member ));
end

end