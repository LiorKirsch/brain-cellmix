function run_R_script(r_function_name,input_MAT_file, output_MAT_file,marker_file)
    script_name = 'tmp_r_script.r';
    create_R_script(script_name, r_function_name, input_MAT_file, output_MAT_file,marker_file);
    system(sprintf('Rscript %s', script_name ));
end


function create_R_script(script_name, r_function_name,input_MAT_file, output_MAT_file,marker_file)

fid = fopen(script_name,'w');
fprintf(fid,'setwd("%s")\n', pwd());
fprintf(fid,'source("cellmix_func.R")\n');
fprintf(fid,'%s("%s","%s","%s","NMF") \n' ,r_function_name,input_MAT_file,output_MAT_file,marker_file);
fclose(fid);

end