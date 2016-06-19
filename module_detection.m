function [module,p] = module_detection(ppi,sample,num_module,gene,mutation_rate,module_size)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
mutation=sample;
weight_ppi=ppi;
for i=1:size(ppi,2)
    for j=1:size(ppi,2)
        v1=ppi(i,:);
        m1=mutation(i,:);
        w1=sum(mutation(i,:))/size(mutation,2);
       v2=ppi(j,:);
       w2=sum(mutation(j,:))/size(mutation,2);
       m2=mutation(j,:);
       exclusive=sum(m1 | m2)-sum(m1 & m2);
        if i==j
            weight_ppi(i,j)=inf;
            continue;
        end   
        if(sum(v1)==0 || sum(v2)==0)
              weight_ppi(i,j)=0;
        else
            weight_ppi(i,j)=dot(v1,v2)./(norm(v1)*norm(v2))*exclusive;
        end
    end
end
% load('weight_ppi_brca.mat');
% load('weight_ppi.mat');
weight_relation=weight_ppi;
global label;
global candidate;
candidate={};
temp_ppi=ppi;
temp_label=zeros(1,size(ppi,2));
temp_gene=gene;
temp_weight_relation=weight_relation;
% is_connected(temp_ppi);
% temp_gene_connected=temp_gene(label==1);
% temp_ppi_connected=temp_ppi(label==1,label==1);
% temp_weight_ppi_connected=temp_weight_relation(label==1,label==1);
% temp_mutation_connected=mutation(label==1,:);
%------------------------------filter the low aberation gene----
num_sample=size(sample,2);
thres=num_sample*mutation_rate;
list=sum(mutation,2)>=thres;
temp_ppi_connected=temp_ppi(list,list);
temp_weight_ppi_connected=temp_weight_relation(list,list);
temp_gene_connected=temp_gene(list);
temp_mutation_connected=mutation(list,:);
%-----------------------------random set the weight---------------%
module=cell(num_module,1);
p=zeros(num_module,1);
value=1:length(temp_gene_connected);
map2=containers.Map(temp_gene_connected,value);
num_gene=length(temp_gene_connected);
for i=1:num_module
    find_module(temp_ppi_connected,temp_weight_ppi_connected,temp_gene_connected,module_size);
    temp_rs=rank(temp_ppi_connected,temp_mutation_connected,temp_gene_connected);
%     module{i}=temp_rs{1};
%     temp_rs=candidate;
    [rs,p2] = permutation_test(temp_rs,gene,mutation);
    if(p2(1)==-1)
        break;
    end
    module{i}=rs;
    p(i)=p2;
    %   add some tioajian-----
    if(length(module{i})==1)
        lll=map2(module{i}{1});
        temp_weight_ppi_connected(lll,:)=0;
        temp_weight_ppi_connected(:,lll)=0;
        continue;
    end
    if (i>1 && length(intersect(module{i},module{i-1}))==min(length(module{i}),length(module{i-1})))
        l1=map2(module{i}{1});
        l2=map2(module{i}{2});
        temp_weight_ppi_connected(l1,l2)=0;
        temp_weight_ppi_connected(l2,l1)=0; 
        temp_weight_ppi_connected(l1,:)=0;
        temp_weight_ppi_connected(:,l1)=0;
        temp_weight_ppi_connected(l2,:)=0;
        temp_weight_ppi_connected(:,l2)=0;
    end
    %--------
    list=[];
for j=1:length(module{i})
    can_gene=module{i}{j};
    list=[list;map2(can_gene)];
end
   num_rand=length(list)*(length(list)-1)/2;
   num=0;
   w=zeros(1,num_rand);
   while(1)
       if num>num_rand
           break;
       end
       m1=unidrnd(num_gene);
       m2=unidrnd(num_gene);
       in_len=length(intersect([m1,m2],list));
       if(m1~=m2 && in_len<2)
           num=num+1;
           w(num)=temp_weight_ppi_connected(m1,m2);
       end
   end
   num=1;
   for i=1:length(list)
       for j=i+1:length(list)
           temp_weight_ppi_connected(list(i),list(j))=w(num);
           temp_weight_ppi_connected(list(j),list(i))=w(num);
           num=num+1;
       end
   end
end
end

