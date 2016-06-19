function [rs] = is_connected(ppi)
global label;
%UNTITLED7 此处显示有关此函数的摘要
%   此处显示详细说明
len=size(ppi,2);
label=zeros(1,len);
v=1;
depth_search(v,ppi);
if(sum(label)<len)
    rs=0;
else
    rs=1;
end
end

function depth_search(v,ppi)
global label;
label(v)=1;
for i=1:size(ppi,2)
    if(label(i)==0 && ppi(v,i)==1)
        depth_search(i,ppi);
    end
end
end
