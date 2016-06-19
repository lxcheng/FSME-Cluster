function [module,p] = FSME_Cluster(ppi_name,mutation_name,num_sample,num_module,mutation_rate,module_size)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
fid = fopen(ppi_name);
pair=textscan(fid,'%s %s %f');
fclose(fid);
%---------judge whether there exis any dupcate gene pair----------------------------
for line=1:length(pair{1})-1
    if(strcmp(pair{1}{line},pair{1}{line+1})==1 && strcmp(pair{2}{line},pair{2}{line+1}==1))
        fprintf('%d',line);
    end
end
for line=1:length(pair{1})-1
    first=pair{1}{line};
    second=pair{2}{line};
    if((strcmp(pair{1}{line+1},first)==1 && strcmp(pair{2}{line+1},second)==1) ||(strcmp(pair{1}{line+1},second)==1 && strcmp(pair{2}{line+1},first)==1))
       fprintf('%d %s %s',line,first,second);
    end
end
%-----------------------------create ppi network-------------------------------------%
gene=unique([pair{1};pair{2}]);
s_gene=length(gene);
value=1:s_gene;
map=containers.Map(gene,value);
ppi=zeros(s_gene,s_gene);
for i=1:length(pair{1})
    first=pair{1}{i};
    second=pair{2}{i};
    l1=map(first);
    l2=map(second);
    ppi(l1,l2)=1;
    ppi(l2,l1)=1;
end
for i=1:s_gene
    ppi(i,i)=1;
end
%------------------------------------load gbm mutation data----------
fm=repmat('%f',1,num_sample);
fm=['%s',fm];
fid=fopen(mutation_name);
gbm=textscan(fid,fm);
fclose(fid);
%--------------------------------------------
p_g=[];
[~,gene_list,list]=intersect(gene,gbm{1});
gene=gene(gene_list);
mutation=zeros(length(gene),length(gbm)-1);
for j=1 :size(mutation,2)
    mutation(:,j)=gbm{1,j+1}(list);
end

%------------------------------------------------run the detection--
sample=mutation;
[module,p] = module_detection(ppi,sample,num_module,gene,mutation_rate,module_size);

end

