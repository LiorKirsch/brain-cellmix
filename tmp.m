algorithem = 'DECONF'; %'NMF'; 'NMFforb'; 'meanProfile'; 'DSA'; 'DECONF';
% 
% input_mat_file_name = 'deconv input/ginger.mat';
% marker_file_name = 'markers sets/ginger_brain_markers_human.txt';
% output_mat_file_name = 'cellmix results.mat';



input_mat_file_name = 'deconv input/su2004.mat';
marker_file_name = 'markers sets/brain_markers_for_su.txt';
output_mat_file_name = 'cellmix results.mat';

check_marker_appear_in_data(marker_file_name, input_mat_file_name);
run_R_script('do_cellmix',input_mat_file_name, output_mat_file_name,marker_file_name,algorithem);