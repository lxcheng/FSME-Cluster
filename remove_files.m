%remove files from RME2 to RME3
noise_list=[5 7 9 11];
% p_list=[5 7 9 11 13];
p_list=[15];
source_path='data/rs/RME2';
des_path='data/rs/RME3';
if(exist(des_path,'dir')==0)
    mkdir(des_path);
end
for i=1:length(noise_list)
    for j=1:length(p_list)
        source_name=[source_path,'/n',num2str(noise_list(i)),'/p',num2str(p_list(j)),'/rs_summary',num2str(p_list(j)),'.txt'];
        des_name=[des_path,'/rs_summary','n',num2str(noise_list(i)),'_p',num2str(p_list(j)),'.txt'];
        copyfile(source_name,des_name);
    end
end

        
        
        
        
        
        
        
        
        
        
        