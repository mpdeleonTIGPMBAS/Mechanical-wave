% (Main) Analyze files/directories sequentially
clear

wdPath = './';  % working directory path
folder_ = dir(wdPath);     % get "strucutre" within the given folder
LisFold = [folder_.isdir];   % folder_.isdir = 1 for directory, =0 for files
LisFold_name = {folder_(LisFold).name};  % get list of folders within the given folder; (!! there are folders "." & "..")
LisFile_name = {folder_(~LisFold).name}; % get list of files...

for i=3:length(LisFold_name)  % go through all the directories in the current folder, except for "." & ".."
  
    display(['Now analyze:  ' LisFold_name{i}]);
    wdPath2 = ['./' LisFold_name{i}];  % get path of the file in the working directory
    
    folder2_ = dir(wdPath2);
    LisFold2 = [folder2_.isdir];
    LisFold2_name = {folder2_(LisFold2).name};
    LisFile2_name = {folder2_(~LisFold2).name};
    
    %(1)Image sequence
    imgFolder = [];
    directory_suffix = 'Seg-outline_origsize';  % directory_suffix = 'images';
    for ii=1:length(LisFold2_name)
        if contains(LisFold2_name{ii},directory_suffix)
            imgFolder= ['./' LisFold_name{i} '/' LisFold2_name{ii} '/'];
        end
    end    
    %(2)Nucleus position
    dataPosFolder = [];
    directory_suffix = 'ROI_excel';   % directory_suffix = 'csv';
    for ii=1:length(LisFold2_name)
        if contains(LisFold2_name{ii},directory_suffix)
            dataPosFolder= ['./' LisFold_name{i} '/' LisFold2_name{ii} '/'];
        end
    end
    %(3)CMZ
    fname_CMZ = [];
    file_suffix = 'CMZ_detection';
    for ii=1:length(LisFile2_name)
        if contains(LisFile2_name{ii},file_suffix)
            fname_CMZ= ['./' LisFold_name{i} '/' LisFile2_name{ii}];
        end
    end
    %(4)Kymograph
    fname_kymo = [];
    file_suffix = 'img';
    for ii=1:length(LisFile2_name)
        if contains(LisFile2_name{ii},file_suffix)
            fname_kymo= ['./' LisFold_name{i} '/' LisFile2_name{ii}];
        end
    end
    
    
    % Call function for Triangulation analysis
    spacing(imgFolder,dataPosFolder,fname_CMZ,fname_kymo);
    
    % stop temporarily
    if i<length(LisFold_name)
        nexDir = LisFold_name{i+1};
        mydlg = warndlg(['Continue to analyze  ' num2str(nexDir) '?'], 'A Warning Dialog');
        waitfor(mydlg);
        close all;
    end
end