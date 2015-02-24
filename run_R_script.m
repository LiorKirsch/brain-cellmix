function run_R_script(r_function_name,input_MAT_file, output_MAT_file,marker_file,alg)
    script_name = 'tmp_r_script.r';
    create_R_script(script_name, r_function_name, input_MAT_file, output_MAT_file,marker_file,alg);
    system(sprintf('Rscript %s', script_name ));
end


function create_R_script(script_name, r_function_name,input_MAT_file, output_MAT_file,marker_file,alg)

fid = fopen(script_name,'w');
fprintf(fid,'setwd("%s")\n', pwd());
fprintf(fid,'source("cellmix_func.R")\n');
fprintf(fid,'%s("%s","%s","%s","%s") \n' ,r_function_name,input_MAT_file,output_MAT_file,marker_file,alg);
fclose(fid);

end