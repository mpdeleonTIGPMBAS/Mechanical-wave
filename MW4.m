% (Main) Kymograph analysis
%==== Analysis of Kymographs =====
clearvars  % clear pre-assigned variables

%%
% Defining kymos to be analyzed
%--------------------------------------------------------------------------
filePath = './';  % the file path for the folder to be analyzed

folder_ = dir(filePath);     % get "strucutre" within the given folder
LisFold = [folder_.isdir];   % folder_.isdir = 1 for directory, =0 for files
numFolders = length(LisFold(LisFold));  % number of folders
numFiles = length(LisFold(~LisFold));   % number of files (!! there are folders "." & "..")
LisFold_name = {folder_(LisFold).name};  % get list of folders within the given folder
LisFile_name = {folder_(~LisFold).name}; % get list of files...

% show folders of given suffix
% ===========================
% folder_suffix = 'backup';   % folder of interest
% for i=1:length(LisFold_name)  
%     if contains(LisFold_name{i},folder_suffix)  % check if the folder name contains 'folder_suffix'      
%     %:------ analysis of specific files inside the specified folder ------
%        folder2_ = dir(LisFold_name{i});
%        LisFold2 = [folder2_.isdir];
%        LisFold2_name = {folder2_(LisFold2).name}; 
%        LisFile2_name = {folder2_(~LisFold2).name};
%        
%        file2_suffix = 'Kymograph';
%        idx_file2 = 1;
%        for ii=1:length(LisFile2_name)
%            if contains(LisFile2_name{ii},file2_suffix)
%                quest = ['Analyze ' LisFile2_name{ii} ' ?'];
%                set(groot,'defaultUicontrolFontSize', 16);  % change font size
%                answer = questdlg(quest,'Analysis','Yes','No','Yes');
%                if strcmp(answer,'Yes')
%                    ImgAnalysis(LisFile2_name{ii}, idx_file2);  % perform Kymo analysis
%                    idx_file2 = idx_file2 + 1;
%                end               
%            end
%        end
%     %:------ analysis of specific files inside the specified folder ------
%     end
% end
% ===========================

% show kymos of given suffix
% ===========================
% prepare the list of kymographs with different amputation levels; get 'Lis0fFiles'
file_suffix = ["25%","50%","75%"];  % keyword of the kymos to be analyzed
Lis0fFiles = [];

for i=1:length(file_suffix)    % over all the amputation levels
    list_ = [];
    for ii=1:length(LisFile_name)  % search for the kymos with the given keyword
        if contains(LisFile_name{ii},file_suffix(i)) & contains(LisFile_name{ii},".tif")
            list_{end+1} = convertCharsToStrings(LisFile_name{ii}); % list of files of interest   
        end
    end
    
    if length(list_) < 0.5  % empty result
        Lis0fFiles{i} = [];
    else
        Lis0fFiles{i} = [list_{:}];   % Lis0fFiles{i} =[file1,file2,...]; i = {25%, 50%, 75%}
    end 
end
%--------------------------------------------------------------------------
%%

%%
% Analyze {25%,50%,75%} Kymos sequentially
%--------------------------------------------------------------------------
choice = {'25%','50%','75%'};
set(groot,'defaultUicontrolFontSize', 14);
[choice_indx, tf] = listdlg('PromptString','Select an amputation level:',...
                     'SelectionMode','multiple','ListString',choice,'ListSize',[350,150]);
                 
dlg_title = {'25% ','50% ','75% '};
for LevTick = 1:length(choice_indx)
    ampuLevel= choice_indx(LevTick); % assign the amputation level   
    idx_file = 1;  % for counting how many times have the kymos been analyzed
    for i=1:length(Lis0fFiles{ampuLevel})  % over all the kymos with given amputation levels         
        file_suffix = Lis0fFiles{ampuLevel}(i);  % kymos to be analyzed  
        for ii=1:length(LisFile_name)
            if contains(LisFile_name{ii},file_suffix)  
               quest = ['Analyze ' LisFile_name{ii} ' ?'];
               set(groot,'defaultUicontrolFontSize', 16);  % change font size
               answer = questdlg(quest,[dlg_title{ampuLevel} 'Analysis'],'Yes','No','Yes');
               if strcmp(answer,'Yes')
                   ImgAnalysis(LisFile_name{ii}, idx_file);  % perform Kymo analysis
                   idx_file = idx_file+1;
               end
            end
        end         
    end
end
%--------------------------------------------------------------------------
%%
% Completion message
%--------------------------------------------------------------------------
waitfor(msgbox('Analysis Completed, Congrat!'));
%--------------------------------------------------------------------------
%%