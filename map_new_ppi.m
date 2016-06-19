%----create new ppi  for GBM data----
ppi_all_name='D:\figure\ÐÂ·½Ïò\data\Multinet.interactions.txt';
fid_ppi=fopen(ppi_all_name);
ppi_all=textscan(fid_ppi,'%s%s','HeaderLines',1);
fclose(fid_ppi);
%-------run GBM data ------
name='brca_somatic_mutations_matrix.txt';
fid=fopen(name);
num_sample=1078;
format='%s';
for i=1:num_sample
    format=[format,'%*d'];
end
data=textscan(fid,format);
fclose(fid);
gene1=data{1};
gene2=ppi_all{1};
gene3=ppi_all{2};

%------------map the ppi--------%
share_gene1=intersect(gene1,gene2);
share_gene2=intersect(gene1,gene3);
share_gene=union(share_gene1,share_gene2);

value=1:length(share_gene);
map=containers.Map(share_gene,value);

%------------output the ppi-----
new_name='new_ppi_brca.txt';
fid2=fopen(new_name,'w');

vv=0.5;
for j=1:length(gene2)
    t1=gene2{j};
    t2=gene3{j};
    if(map.isKey(t1)==1 && map.isKey(t2)==1)
        fprintf(fid2,'%s\t%s\t%f\n',t1,t2,vv);
    end
end
fclose(fid2);



















