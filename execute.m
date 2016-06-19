


%---------------------------default parameters-----------------------------------------------
num_module=10;
mutation_rate=0.03;
module_size=5;
%--------------------------run------------------------------------------------------
[module,p] = FSME_Cluster(ppi_name,mutation_name,num_sample,num_module,mutation_rate,module_size);
[module,p]=remove_repeat(module,p);
clearvars -except module p;
%-------------------------------------------------------------------------------------------