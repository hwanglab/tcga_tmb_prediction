function main_prediction_proportion

tmb_path='E:\Hongming\projects\tcga-bladder-mutationburden\tmb_blca\';
feat_path='E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\step11)_tmb_map_process\features\';

mat_files=dir(strcat(tmb_path,'*.mat'));


pid_hl=cell(length(mat_files),3);

t=0.5;
for i=1:length(mat_files)
    temp_file=mat_files(i).name;
    load(strcat(tmb_path,temp_file));
    
    tmb_level=[];
    for j=1:size(top_left_tumor,1)
        rs=top_left_tumor(j,1);
        cs=top_left_tumor(j,2);
        re=bottom_right_tumor(j,1)-1;
        ce=bottom_right_tumor(j,2)-1;
        %unique(tmb_map(rs:re,cs:ce))
        tmb_level=[tmb_level,mean(mean(tmb_map(rs:re,cs:ce)))];
    end
    
    bb=(tmb_level>t);
    %pl_value=sum(bb)/sum(~bb);
    
    pid_hl{i,1}=temp_file;
    pid_hl{i,2}=sum(bb);
    pid_hl{i,3}=sum(~bb);
end


save(strcat(feat_path,'feat_hl','.mat'),'pid_hl');

