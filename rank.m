function [rs] = rank(ppi,mutation,gene)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
global candidate;
len=length(gene);
value_set=1:len;
gene_map=containers.Map(gene,value_set);
rs={};
temp_rs={};
score=zeros(length(candidate),1);
for i=1:length(candidate)
   can=candidate{i};
   temp_score=0;
   v=[];
   
   for j=1:length(can)
       list1=gene_map(can{j});
       v1=mutation(list1,:);
       v(j,:)=v1;
   end
   s_v=sum(v,2);
   [~,index]=sort(s_v,'descend');
   score(i)=sum(v(index(1),:));
   v1=v(index(1),:);
   ii=1;
   temp_rs{i}{ii}=can{index(1)};
   ii=ii+1;
   for kk=2:length(index)
       v2=v(index(kk),:);
%        if(sum(v2)/size(mutation,2)<=0.05)
%            continue;
%        end
       v3=v1+v2;
       c=sum(v3>0);
       o=sum(v3==2);
       if(score(i)<c-o)
           v1=v3;
           score(i)=c-o;
           temp_rs{i}{ii}=can{index(kk)};
           ii=ii+1;
       end
   end
end
[~,index]=sort(score,'descend');
rs=temp_rs(index);

end

