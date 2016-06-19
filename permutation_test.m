function [rs,p] = permutation_test(temp_rs,gene,mutation)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
l=length(temp_rs);
len=length(gene);
value_set=1:len;
gene_map=containers.Map(gene,value_set);
value_map=containers.Map(value_set,gene);
%------------parameters settings-----------
num_m=0;
thre=0.05;
N=10000;

%------------------------------------------
for i=1:l
    if(isempty(temp_rs{i}))
        break;
    else
        list=zeros(1,length(temp_rs{i}));
        for j=1:length(temp_rs{i})
            temp_gene=temp_rs{i}{j};
            list(j)=gene_map(temp_gene);
        end
        num=size(mutation,2);
        list=sort(list);
        temp_mutation=mutation(list,:);
        num_mutation=sum(temp_mutation,2);
        temp2=sum(temp_mutation);
        score=sum(temp2>0)-sum(temp2>1);
        p=test(score,num,num_mutation',N);
        if(p<thre)
            num_m=num_m+1;
            p_list(num_m)=p;
            rs{num_m}=temp_rs{i};
        else
            [p,list]=delete_test(temp_mutation,list,N,thre);
            if(p<thre)
                num_m=num_m+1;
                p_list(num_m)=p;
                temp_temp={};
                for k=1:length(list)
                    temp_temp{k}=value_map(list(k));
                end
                rs{num_m}=temp_temp;
            end
        end   
    end
end
if(num_m>0)
    [p_list,index]=sort(p_list);
    p=p_list(1);
    new_rs=rs{index(1)};
    rs=new_rs;
else
    p=-1;
    rs={};

end
end



function p=test(score,num,num_mutation,N)
    len=length(num_mutation);
    p=0;
    for k=1:N
        list=zeros(len,num);
        for i=1:length(num_mutation)
            rand_num=randperm(num,num_mutation(i));
            list(i,rand_num)=1;
        end
        m=sum(list);
        temp_score=sum(m>0)-sum(m>1);
        if(score<temp_score)
            p=p+1;
        elseif(score==temp_score)
            p=p+0.5;
        end
    end
    p=double(p)/(1+N);
end
    
        
function [p,new_list]=delete_test(mutation,list,N,thre)
   len=size(mutation,1);
   temp_p=ones(1,len);
   list2=1:len;
   for i=1:len
       temp_list=setdiff(list2,i);
       temp_mutation=mutation(temp_list,:);
       temp2=sum(temp_mutation);
       score=sum(temp2>0)-sum(temp2>1);
       num=size(temp_mutation,2);
       num_mutation=sum(temp_mutation,2);
       p1=test(score,num,num_mutation',N);
       temp_p(i)=p1;
   end
   [p2,index]=sort(temp_p);
    k=index(1);
    temp3=setdiff(list2,k);
   if(p2(1)>thre && length(p2)>=3)
       temp_mutation=mutation(temp3,:);
       delete_test(temp_mutation,list(temp3),N,thre);
   end
   p=p2(1);
   new_list=list(temp3); 
end











