% read and process excel spreadsheet file
sheet=2;
[~,~,raw_patientID]=xlsread('Table_S1.2017_08_05.xlsx',sheet,'A1:A413');
% [~,~,raw_KMT2D]=xlsread('Table_S1.2017_08_05.xlsx',sheet,'BK1:BK413');
% [~,~,raw_KMT2C]=xlsread('Table_S1.2017_08_05.xlsx',sheet,'CB1:CB413');


Datasetpath='E:\bladder_project_DX\dataset\';
Unusedpath='E:\bladder_project_DX\unused\';

start_path=fullfile('Z:\Datasets\TCGA_BLCA_WSI\');
topLevelFolder=uigetdir(start_path);
if topLevelFolder==0
    return;
end

% get list of all subfolders
allSubFolders=genpath(topLevelFolder);
remain=allSubFolders;
listOfFolderNames={};
while true
    [singleSubFolder,remain]=strtok(remain,';');
    if isempty(singleSubFolder)
        break;
    end
    listOfFolderNames=[listOfFolderNames singleSubFolder];
end
numberOfFolders=length(listOfFolderNames);

id=1;
patientNum=0;
unusedNum=0;
for k=1:numberOfFolders
    thisFolder=listOfFolderNames{k};
    fprintf('Processing folder %s\n', thisFolder);
    
    % get filenames of all .svs files
    filePattern=sprintf('%s/*.svs',thisFolder);
    baseFileNames=dir(filePattern);
    numberOfFiles=length(baseFileNames);
    
    if numberOfFiles==1
        fprintf('the %d subfile with .svs\n',id);
        id=id+1;
        
        fullFileName=fullfile(thisFolder,baseFileNames.name);
        filename=baseFileNames.name;
        patientID=filename(1:12);
        
        copyfile(fullFileName,strcat(Datasetpath,baseFileNames.name));
        patientNum=patientNum+1;
        
        %         for h=1:length(raw_patientID)
        %             temp_ID=raw_patientID(h);
        %         end
        
        
    elseif  numberOfFiles>1
        fprintf('the %d subfile with .svs\n',id);
        id=id+1;
        for f=1:numberOfFiles
            fullFileName=fullfile(thisFolder,baseFileNames(f).name);
            filename=baseFileNames(f).name;
            patientID=filename(1:12);
            
            if f==1
                copyfile(fullFileName,strcat(Datasetpath,baseFileNames(f).name));
                patientNum=patientNum+1;
            else
                
                copyfile(fullFileName,strcat(Unusedpath,baseFileNames(f).name));
                unusedNum=unusedNum+1;
            end
            
        end
    else
        disp('no svs file under this folder!!!');
    end
end
