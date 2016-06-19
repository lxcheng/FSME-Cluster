function [new_module,new_p]=remove_repeat(module,p)
% remove the repeat modules and sort the answers
%   module: the functional modules; DataType: cell;
%   p: the conrespongding p_value ;DataType: สýื้ฃป
len=length(module);
r=1;
repeat=[];
for i=1:len
    if(isempty(module{i})==1)
        repeat(r)=i;
        r=r+1;
        continue;
    end
    for k=i+1:len
        if(length(module{i})==length(module{k}) && isempty(setdiff(module{i},module{k}))==1)
            repeat(r)=k;
            r=r+1;
        end
    end
end
if(isempty(repeat)==0)
    repeat=unique(repeat);
end
list=setdiff(1:len,repeat);
new_p=p(list);
temp_module=cell(length(list),1);
new_module=temp_module;
for i=1:length(list)
    temp_module{i}=module{list(i)};
end
[new_p,index]=sort(new_p);
for i=1:length(index)
    new_module{i}=temp_module{index(i)};
end

end

