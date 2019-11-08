function DL_data_generation

load(strcat('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\','blca_MutBurdens.mat')); % load blca using TCGA values with percentiles
image_source='E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\10)image_patches_AP\';
train_data_destination={'E:\Hongming\projects\tcga-bladder-mutationburden\deep_learning_data\data01\train\low\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\deep_learning_data\data01\train\mid\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\deep_learning_data\data01\train\high\'};
validation_data_destination={'E:\Hongming\projects\tcga-bladder-mutationburden\deep_learning_data\data01\validation\low\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\deep_learning_data\data01\validation\mid\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\deep_learning_data\data01\validation\high\'};

rng('default');
cvIndices=crossvalind('Kfold',ones(1,length(blca_MutBurdens)),5);


for i=1:1
    valid_id=(cvIndices==i);
    train_id=~valid_id;
    
    for k=1:length(blca_MutBurdens)
        patientID=blca_MutBurdens{k,1};
        label=blca_MutBurdens{k,3};
        
        temp=dir(strcat(image_source,patientID,'*'));
        if length(temp)>100
            disp(patientID);
        end
        
        if valid_id(k)==1
            switch label
                case 'Low'
                    xu_copy_images(image_source,validation_data_destination{1},temp);
                case 'Mid'
                    xu_copy_images(image_source,validation_data_destination{2},temp);
                case 'High'
                    xu_copy_images(image_source,validation_data_destination{3},temp);
                otherwise
                    disp('impossible');
            end
        else
            switch label
                case 'Low'
                    xu_copy_images(image_source,train_data_destination{1},temp);
                case 'Mid'
                    xu_copy_images(image_source,train_data_destination{2},temp);
                case 'High'
                    xu_copy_images(image_source,train_data_destination{3},temp);
                otherwise
                    disp('impossible');
            end
        end

    end
end