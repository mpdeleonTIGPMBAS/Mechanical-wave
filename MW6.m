%%
%:: Statistics 
%>> extract median of repeated measurments
file_ = './kymoData_All_unitTrans.csv';  % the data file path
T_ = readtable(file_);
%:: Choose the ROI data indexed by idxCol (i.e., good kymos)
row_idx = find(T_.idxCol==1);
T = T_(row_idx,:);

% get list of unique IDs
T.Ray = num2str(T.Ray);  % double to string
T.Ray = cellstr(T.Ray);  % string to cell
[idList_,idIdx_] = unique(T(:,{'Amputation','Experiment','Ray'}), 'rows');
% sorting
[idIdx, id_idex] = sort(idIdx_); % ensure that the list are in order from 25% to 75%
idList = idList_(id_idex, :);


col_header = {'FrontSpeed','FrontSpeed_intercept','backSpeed','backSpeed_intercept','FrontSpeed2','FrontSpeed2_intercept','maxCMZLen','maxCMZDist','waveLen','t_maxCMZlen','t_CMZvanish','cutEdgeSpeed','cutEdgeDist'}; 
col_start = 4; % start from the column "FrontSpeed"
col_end = 16;  % end at the column "cutEdgeDist"

T_rev = [];
idIdx(end+1) =  height(T)+1;  % to sovle the boundary issue
for i=1:height(idList)
    row_start = idIdx(i);
    row_end = idIdx(i+1)-1;

    T_median = median(T{row_start:row_end, col_start:col_end}, 1);  % median
    T_mean = mean(T{row_start:row_end, col_start:col_end}, 1);      % mean
    T_std  = std(T{row_start:row_end, col_start:col_end}, 0, 1);       % standard deviation
    % combine
    T_temp = reshape([T_median;T_mean;T_std], 1,[]);
    
    T_rev = [T_rev; T_temp];   %append the table
end

% header transformation
col_header = reshape([col_header;col_header;col_header], 1,[]);
%-------
col_header_median = cell(1,length(T_median));
col_header_median(:) = {'median'};
col_header_mean = cell(1,length(T_mean));
col_header_mean(:) = {'average'};
col_header_std = cell(1,length(T_std));
col_header_std(:) = {'std'};
col_header2 =  reshape([col_header_median;col_header_mean;col_header_std], 1,[]);
%-------
col_header_new = [];
for i=1:length(col_header)
    newCol_header = join([col_header(i),col_header2(i)],"_");
    col_header_new{i} = newCol_header{1};
end

T_rev = array2table(T_rev, 'VariableNames',col_header_new);
T_new = horzcat(idList, T_rev);  % join two tables horizontally

file_ = './kymoData_All_unitTrans_statistics.csv';
writetable(T_new,file_);  % output the file
%==========================================================================