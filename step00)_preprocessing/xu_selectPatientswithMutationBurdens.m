% read and process excel spreadsheet file
%sheet=2;
%[~,~,raw_patientID]=xlsread('Table_S1.2017_08_05.xlsx',sheet,'A1:A413');
% [~,~,raw_KMT2D]=xlsread('Table_S1.2017_08_05.xlsx',sheet,'BK1:BK413');
% [~,~,raw_KMT2C]=xlsread('Table_S1.2017_08_05.xlsx',sheet,'CB1:CB413');

load(strcat('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\','blca_MutBurdens.mat'));


%Datasetpath={'E:\blca_mutationBurden\low\','E:\blca_mutationBurden\mid\','E:\blca_mutationBurden\high\'};
Datasetpath={'E:\blca_mutationBurden\tumor2\low\','E:\blca_mutationBurden\tumor2\mid\','E:\blca_mutationBurden\tumor2\high\'};
start_path=fullfile('E:\blca_mutationBurden\tumor_detection\');
%start_path=fullfile('Z:\Datasets\TCGA_BLCA_WSI\');
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
for k=1:numberOfFolders
    thisFolder=listOfFolderNames{k};
    fprintf('Processing folder %s\n', thisFolder);
    
    % get filenames of all .svs files
    filePattern=sprintf('%s/*.svs',thisFolder);
    baseFileNames=dir(filePattern);
    numberOfFiles=length(baseFileNames);
    
    
    fprintf('the %d subfile with .svs\n',id);
    id=id+1;
    for f=1:numberOfFiles
        fullFileName=fullfile(thisFolder,baseFileNames(f).name);
        filename=baseFileNames(f).name;
        patientID=filename(1:12);
        
        indicator=0;
        for pi=1:size(blca_MutBurdens,1)
            pID=blca_MutBurdens{pi,1};
            
            if strcmp(patientID,pID)
                mu=blca_MutBurdens{pi,3};
                switch mu
                    case 'Low'
                        copyfile(fullFileName,strcat(Datasetpath{1},baseFileNames(f).name));
                    case 'Mid'
                        copyfile(fullFileName,strcat(Datasetpath{2},baseFileNames(f).name));
                    case 'High'
                        copyfile(fullFileName,strcat(Datasetpath{3},baseFileNames(f).name));
                    otherwise
                        disp('impossible!!!!');
                end
                indicator=1;
                break;
            end
        end 
        if indicator==0
            disp('impossible error!!!!');
        end
        
    end
    
end
