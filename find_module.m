function  find_module(ppi,weight_relation,gene,module_size)
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明
%----------------------------------calculate similarity-----------
%similar=cal_similar(ppi);
% weight=similar.*weight_relation;
weight=weight_relation;
global candidate;
candidate={};
%----------------remove the low value----------------------------------
list=weight==0;
weight(list)=inf;
ppi(list)=0;

%-----------------temp_rs-------------------
global temp_rs;
global in;
global label;
temp_rs=cell(size(ppi,2),1);
ll=size(ppi,2);
for i=1:size(ppi,2)
    temp_rs{i}=cell(3,1);
end
in=0;
while(1)
    if(is_connected(ppi)==1)
        gene1=gene(label==1);
        weight1=weight(label==1,label==1);
        ppi1=ppi(label==1,label==1);
        in=in+1;
        temp_rs{in}{1}=gene1;
        temp_rs{in}{2}=weight1;
        temp_rs{in}{3}=ppi1;
%         gene=gene2;
%         weight=weight2;
%         ppi=ppi2;
        break;
    else
         gene1=gene(label==1);
         gene2=gene(label==0);
         weight1=weight(label==1,label==1);
         weight2=weight(label==0,label==0);
         ppi1=ppi(label==1,label==1);
         ppi2=ppi(label==0,label==0);
         in=in+1;
         temp_rs{in}{1}=gene1;
         temp_rs{in}{2}=weight1;
         temp_rs{in}{3}=ppi1;
         gene=gene2;
         weight=weight2;
         ppi=ppi2;
         
    end
end
%-------------------------------delete node----------------------------

for j=1:ll
    if(isempty(temp_rs{j}{1})==1)
        break;
    end
    gene=temp_rs{j}{1};
    weight=temp_rs{j}{2};
    ppi=temp_rs{j}{3};
    delete_node(weight,ppi,gene,module_size);
end

end

function delete_node(weight,ppi,gene,module_size)                                                                       
    global candidate;
    global label;
    global in;
    global temp_rs;
%     
%             if (is_connected(ppi)==0)
%                 gene1=gene(label==1);
%                 gene2=gene(label==0);
%                 weight1=weight(label==1,label==1);
%                 weight2=weight(label==0,label==0);
%                 ppi1=ppi(label==1,label==1);
%                 ppi2=ppi(label==0,label==0);
%                 if(length(gene1)>=2)
%                     delete_node(weight1,ppi1,gene1);
%                 end
%                 if(length(gene2)>=2)
%                     delete_node(weight2,ppi2,gene2);
%                 end
%             else
while(1)
    len=size(ppi,2);
    if(len<2)
        break;
    end
    if(len>=2 && len<=module_size)
        if(isempty(candidate))
            candidate={gene};
            break;
        else
         candidate=[candidate;{gene}];
            break
        end
    end
    [value,index]=min(weight);
    [~,index2]=min(value);
     y=index2;
     x=index(index2);
     if(ppi(x,y)==0)
          weight(y,x)=inf;
          weight(x,y)=inf;
      else
          ppi(x,y)=0;
          ppi(y,x)=0;
          weight(x,y)=inf;
          weight(y,x)=inf;
          if(is_connected(ppi)==0)
              if(sum(label)==1)
                  gene=gene(label==0);
                  weight=weight(label==0,label==0);
                  ppi=ppi(label==0,label==0);
              elseif(sum(label==0)==1)
                  if(length(label)==1)
                      break;
                  end
                  gene=gene(label==1);
                  weight=weight(label==1,label==1);
                  ppi=ppi(label==1,label==1);
              else
                  gene1=gene(label==1);
                  gene2=gene(label==0);
                  weight1=weight(label==1,label==1);
                  weight2=weight(label==0,label==0);
                  ppi1=ppi(label==1,label==1);
                  ppi2=ppi(label==0,label==0);
                  in=in+1;
                  temp_rs{in}{1}=gene2;
                  temp_rs{in}{2}=weight2;
                  temp_rs{in}{3}=ppi2;
                  gene=gene1;
                  weight=weight1;
                  ppi=ppi1;
              end
          end
     end
end
end

function [similar]=cal_similar(ppi)
   similar=zeros(size(ppi,1),size(ppi,1));
   for i=1:size(ppi,1)
       for j=i:size(ppi,1)
           v1=ppi(i,:);
           v2=ppi(j,:);
           ww=dot(v1,v2)/(norm(v1)*norm(v2));
           similar(i,j)=ww;
           similar(j,i)=ww;
       end
   end
%     w_change=reshape(similar,1,[]);
%  w_mean=mean(w_change);
%  w_std=std(w_change);
%  similar=(similar-w_mean)./w_std;
%  min_min=min(min(similar));
%  max_max=max(max(similar));
% for i=1:size(similar,1)
%     for j=1:size(similar,2)
%         if(similar(i,j)==min_min)
%             similar(i,j)=0.001;
%         else
%             similar(i,j)=(similar(i,j)-min_min)/(max_max-min_min);
%             if(similar(i,j)<0.001)
%                 similar(i,j)=0.001;
%             end
%         end
%     end
% end
end

%



function [similar]=cal_similar2(ppi)
   similar=zeros(size(ppi,1),size(ppi,1));
   for i=1:size(ppi,1)
       for j=i:size(ppi,1)
           v1=ppi(i,:);
           v2=ppi(j,:);
           ww=sum((v1==1) & (v2==1));
           similar(i,j)=ww;
           similar(j,i)=ww;
       end
   end
   
end



