%----------execute the CE ----
% ppi_name='D:\figure\新方向\newcode\network\data\ppi.filterd.txt';
% mutation_name='D:\figure\新方向\data\gbm_mutation2.txt';
% num_sample=145;
% mutation_name='gbm_somatic_matrix.txt';
% num_sample=550;
% ppi_name='new_ppi.txt';
mutation_name='brca_somatic_mutations_matrix.txt';
num_sample=1078;
ppi_name='new_ppi_brca.txt';

num_module=10;
mutation_rate=0.03;
module_size=5;
[module,p] = mm(ppi_name,mutation_name,num_sample,num_module,mutation_rate,module_size);
[module,p]=remove_repeat(module,p);
clearvars -except module p;
