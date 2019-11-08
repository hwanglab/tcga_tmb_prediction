%% divide 1024x1024 tiles to 512x512 tiles for training DL model

%%1) for training pathes
path_train={'Y:\DL_TMB_prediction\data_bladder_20X_norm\training\low\','Y:\DL_TMB_prediction\data_bladder_20X_norm\training\high\'};
des_train={'Y:\DL_TMB_prediction\data_bladder_20X_norm\fold01\training\low\','Y:\DL_TMB_prediction\data_bladder_20X_norm\fold01\training\high\'};

for ip=1:length(path_train)
    path_train_temp=path_train{ip};
    des_train_temp=des_train{ip};
    
    imgs=dir(strcat(path_train_temp,'*.png'));
    
    for ig=1:length(imgs)
        I=imread(strcat(path_train_temp,imgs(ig).name));
        disp(imgs(ig).name);
        for i=1:2
            for j=1:2
                RGB=I((i-1)*512+1:i*512,(j-1)*512+1:j*512,:);
                imgnn=strcat(imgs(ig).name(1:end-3),num2str(i),num2str(j),'.png');
                imwrite(RGB,strcat(des_train_temp,imgnn));
            end
        end
    end
end

%%2) for validation patches
path_validation={'Y:\DL_TMB_prediction\data_bladder_20X_norm\validation\low\','Y:\DL_TMB_prediction\data_bladder_20X_norm\validation\high\'};
des_validation={'Y:\DL_TMB_prediction\data_bladder_20X_norm\fold01\validation\low\','Y:\DL_TMB_prediction\data_bladder_20X_norm\fold01\validation\high\'};

for ip=1:length(path_validation)
    path_validation_temp=path_validation{ip};
    des_validation_temp=des_validation{ip};
    
    imgs=dir(strcat(path_validation_temp,'*.png'));
    for ig=1:length(imgs)
        I=imread(strcat(path_validation_temp,imgs(ig).name));
        disp(imgs(ig).name);
        for i=1:2
            for j=1:2
                RGB=I((i-1)*512+1:i*512,(j-1)*512+1:j*512,:);
                imgnn=strcat(imgs(ig).name(1:end-3),num2str(i),num2str(j),'.png');
                imwrite(RGB,strcat(des_validation_temp,imgnn));
            end
        end
    end
end