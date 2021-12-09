% Collect data for plotting

%%
% Analyze the folder structure that stores the CSV files
%==========================================================================
filePath = './'; % the file path for the folder to be analyzed

folder_ = dir(filePath);     % get "strucutre" within the given folder
LisFold = [folder_.isdir];   % folder_.isdir = 1 for directory, =0 for files
numFolders = length(LisFold(LisFold));  % number of folders (!! there are folders "." & "..")
numFiles = length(LisFold(~LisFold));   % number of files 
LisFold_name = {folder_(LisFold).name};  % get list of folders within the given folder
LisFile_name = {folder_(~LisFold).name}; % get list of files...
%==========================================================================

%%
% Get list for data to be collected
%==========================================================================
file_suffix = ["50%","75%"];  % keywords of CSV files to be collected
Lis0fCSV = [];

for i=1:length(file_suffix)    % over all the amputation levels
    list_ = [];
    for ii=1:length(LisFile_name)  % search for the kymos with the given keyword
        if contains(LisFile_name{ii},file_suffix(i)) & contains(LisFile_name{ii},".csv")
            list_{end+1} = convertCharsToStrings(LisFile_name{ii}); % list of files of interest   
        end
    end
    
    if length(list_) < 0.5  % empty result
        Lis0fCSV{i} = [];
    else
        Lis0fCSV{i} = [list_{:}];   % Lis0fFiles{i} =[file1,file2,...]; i = {25%, 50%, 75%}
    end 
end
%==========================================================================

%%
% Collect CSV files, unit conversion, statistical analyses, data output
%==========================================================================
% :: Data preparation
T_new = [];
for i=1:length(file_suffix)            % go over [25%, 50%, 75%]
    for jj=1:length(Lis0fCSV{i})       % go over CSV within the amputation level
        file_ = Lis0fCSV{i}(jj);       % the CSV file
        
        file_info = split(file_, '-'); % get details of the CSV file   
        AmpLev = char(file_info(1));   % amputation level
        ExpIdx = char(file_info(2));   % experiment number
        kymo_info = split(file_info(3),'_');
        RayIdx = char(kymo_info(2));         % ray number
        
        T = readtable(file_);
        numRows = height(T);           % number of rows
        
        % clear the NaN data
        TF = ismissing(T,{NaN});       % get matrix for recording NaN element
        NA_colIndx = any(TF,1);        % column index of the NaN column
        T(:, NA_colIndx) = [];         % delete NaN column
        emptyCol = cell(numRows,1);    % empty content
        T = addvars(T,emptyCol,'Before','Comments'); % replace NaN column by empty column
        
        Ray = repmat({RayIdx},numRows,1);          % repeat copies of array
        T = addvars(T,Ray,'Before','FrontSpeed');
        Experiment = repmat({ExpIdx},numRows,1);  
        T = addvars(T,Experiment,'Before','Ray');
        Amputation = repmat({AmpLev},numRows,1);     
        T = addvars(T,Amputation,'Before','Experiment');
        
        T_new = [T_new; T];            % concatenate all the CSV files 
    end
end
%:: Save Data (pixel-based raw data)
fptr = './kymoData_All.csv';
writetable(T_new, fptr);  % rewrite the file


%:: Units transformation
T_new_trans = T_new;
unitTime = (1/20)*2;       % pixel-to-min conversion (X is scaled by 20-fold)
unitLeng = 1/0.88;         % pixel-to-um conversion  (Y is scaled by 1-fold)

%:> List of Speed [L/T]
List_LT = {'FrontSpeed','backSpeed','FrontSpeed2','cutEdgeSpeed'};
%:> List of Length [L]
List_L = {'FrontSpeed_intercept','backSpeed_intercept','FrontSpeed2_intercept','maxCMZLen','maxCMZDist','waveLen','cutEdgeDist'};
%:> List of Time [T]
List_T = {'t_maxCMZlen','t_CMZvanish'};

% unit conversion
for i=1:length(List_LT)
    T_new_trans{:, List_LT{i}} = T_new_trans{:, List_LT{i}}.*(unitLeng/unitTime);
end
for i=1:length(List_L)
    T_new_trans{:, List_L{i}} = T_new_trans{:, List_L{i}}.*(unitLeng);
end
for i=1:length(List_T)
    T_new_trans{:, List_T{i}} = T_new_trans{:, List_T{i}}.*(unitTime);
end
%:: Save Data (real-unit data)
fptr = './kymoData_All_unitTrans.csv';
writetable(T_new_trans, fptr);  % rewrite the file


%:: Statistics (mean & s.d.)
% >> do mean and standard deviation -> insert after the given column

%==========================================================================