%% Analysis of tailfin thickness
clear

%% Data import
fname = 'M296_Reslice.tif';  % put the tiff name that you want to analyze 
info = imfinfo(fname);  % info of the import figure; channel number = length(info) 
img = []; 
img_b = [];

for i=1:length(info)
    channel = i;    % selected channel
    [img_] = imread(fname,channel);  % read image of the specific channel

    scalFacX = 1;   % scaling factor
    scalFacY = 1;
    [nRows, nCols, ~] = size(img_); 
    img{i} = imresize(img_, [round(nRows*scalFacY) round(nCols*scalFacX)],'nearest');  % adjust image size if necessary

    % Normalization & conversion into double format
    img{i} = mat2gray(img{i});
    img_sm{i} = imgaussfilt(img{i},2);  % Gaussian smoothing filter
    
    % Binarization
    Tc = graythresh(img_sm{i});       % automatic decision of pixel threshold
    img_b{i} = double(img_sm{i}>Tc);  % binary image
end


%% Manual detection of measurement range
%:: create folders to save the analyzed results
fdname_ = regexp(fname,'.tif','split');  % split filename
fdname = [fdname_{1} '_analyzed'];
fpath = sprintf('./%s/',fdname);         % file path of the folder
if ~exist(fpath, 'dir')
    mkdir(fpath);                        % create the folder -- for saving figs
end
% create new sub-folders to save the current analysis results
folder_ = dir(fpath);
LisFold = [folder_.isdir];
numFolders = length(LisFold(LisFold))-2;  % skip './' & '../'
fpath = [fpath 'Measure' num2str(numFolders+1) '/'];
mkdir(fpath);

%:: Mannual correction
dy = 5;  % y-resolution in unit of pixel
yTop = [];
yBot = [];
for i=1:length(info)  % go through each channel
    %:: rule-based selection of top & bottom lines for the specific channel
    for yy=1:dy:nRows
        if find(img{i}(yy,:)>0.5)   % Region Of Interest (ROI)
            yTop{i} = yy;           % top line
            break;
        end
    end
    for yy=nRows:-dy:1
        if find(img{i}(yy,:)>0.5)   % Region Of Interest (ROI)
            yBot{i} = yy;           % bottom line
            break;
        end
    end
        
    % visualization + manual correction
    %----------------------------------
    answer = 'No';
    numMod_(1) = 0.0;   % initialization
    numMod_(2) = 0.0;
    yTop0 = yTop{i};    % keep the original lines
    yBot0 = yBot{i};
    while strcmp(answer,'No')   % continue to correct, if not satisfied
        close(gcf);   % close image
        figure('visible','on');
        set(gcf, 'Position',  [50, 150, 1200, 500]);
        subplot(1,2,1);
        imshow(img{i})  % should overlay with the original figure
        hold on;
        plot([1,nCols], yTop{i}.*ones(1,2), 'm-', 'Linewidth', 1.5);    % plot the top line
        plot([1,nCols], yBot{i}.*ones(1,2), 'g-', 'Linewidth', 1.5);    % plot the Bottom line

        subplot(1,2,2);
        imshow(img_b{i})
        hold on;
        plot([1,nCols], yTop{i}.*ones(1,2), 'm-', 'Linewidth', 1.5);    % plot the top line
        plot([1,nCols], yBot{i}.*ones(1,2), 'g-', 'Linewidth', 1.5);    % plot the Bottom line
        hold off;

        quest = sprintf('Is ROI OK?');
        set(groot,'defaultUicontrolFontSize', 16);  % change font size
        answer = questdlg(quest,'ROI','Yes','No','Yes');

        if strcmp(answer,'No')  % correction
            prompt = {['\fontsize{12} Top line correction Y1:' sprintf('\n') ' (1+Y1) * Ytop  |  ' sprintf('Ytop0=%.3f',yTop0)], ...
                      ['\fontsize{12} Bot line correction Y2:' sprintf('\n') ' (1+Y2) * Ybot  |  ' sprintf('Ybot0=%.3f',yBot0)]};
            dlgtitle = sprintf('X_{pixels} = %d ; \t Y_{pixels} = %d',nCols,nRows);
            dims = [2 80];  % dialog box size
            definput = {num2str(numMod_(1)),num2str(numMod_(2))};  % default values
            opts.Interpreter = 'tex';  % need to switch to "tex" mode, in order to change font size
            answer2 = inputdlg(prompt,dlgtitle,dims,definput,opts);  
                
            yTop{i}= yTop0*(1+str2num(answer2{1})); % slope correction
            yBot{i}= yBot0*(1+str2num(answer2{2})); % intercept correction
            
            numMod_(1) = str2num(answer2{1});  % save the previous input
            numMod_(2) = str2num(answer2{2});
        end
    end
    saveas(gca,fullfile([fpath 'img_ROI_C' num2str(i) '.jpeg'])); % show fitting line with images
    close(gcf);   % close image
    %----------------------------------
end    


%% Detection of ROI at each y-location
% range for profile scanning
yTop_ = max(yTop{:});  % upper boundary of the scanning
yBot_ = min(yBot{:});  % bottom boundary...

locsXR = [];  % X locations of the Right boundaries; locsXR{channel,point_index}
locsXL = [];  % ... of the Left boundaries
locsY = [];
nDY = 50;   % Y-resolution (# of y-pieces)
ySpace = linspace(yTop_,yBot_,nDY);  % y = [yTop_, yTop_+nDY, yTop_+2nDY, ..., yBot_]
for j = 1:length(ySpace(1:end-1)) 
     y1 = ySpace(j+0);
     y2 = ySpace(j+1);    
     for ii=1:length(info)  % go through each channel
         col_ListR = [];   % initialization
         col_ListL = [];
         col_temp  = [];
         for yy = round(y1):round(y2) % y-scanning
             col_temp = find(img_b{ii}(yy,:)==1); % find ROI at the specific y 
             if length(col_temp)>1.5   % signals from the two boundaries
                 col_ListR(end+1) = max(col_temp); % find the right-most point
                 col_ListL(end+1) = min(col_temp); % find the left-most point
             end
         end
         if length(col_ListR)>0.5  % nonzero data
             locsY{ii,j} = round(0.5*(y1+y2));
             locsXR{ii,j} = prctile(col_ListR,50);  % choose 'median' as the representative point
             locsXL{ii,j} = prctile(col_ListL,5);   % "1"- left-most, "100"- right-most, "50"- middle
         end
     end   
end


%% Visualization & Data output
figure('visible','on');
set(gcf, 'Position',  [50, 150, 1200, 500]);
for i=1:length(info)
    subplot(1,2,i);
    hold on;
    imshow(img_b{i})  % channel-i image
    plot([locsXR{i,:}],[locsY{i,:}], 'or', 'Linewidth', 1.5);    % plot the fitting points
    plot([locsXL{i,:}],[locsY{i,:}], '^g', 'Linewidth', 1.5);    % plot the fitting points
    title(['C' num2str(i)]);
end
hold off;
saveas(gca,fullfile([fpath 'img_ROI_fitting.jpeg'])); % show fitting line with images
savefig([fpath 'img_ROI_fitting2.fig']);
%close(gcf);   % close image

%%
%:: Data pre-processing -> collect columns of non-empty data
fpath = ['./' fdname_{1} '.xls'];
dataXR_C1 = []; dataXR_C2 = [];
dataXL_C1 = []; dataXL_C2 = [];
dataY = [];

[~,col_empty_] = find(cellfun(@isempty,locsXR));  % find index of the empty data
col_empty = unique(col_empty_);  % find the column index for which the empty data exist
for i=1:length(locsXR)
    if ~any(col_empty==i)  % column of non-empty data
        dataXL_C1(end+1) = [locsXL{1,i}];
        dataXR_C1(end+1) = [locsXR{1,i}];
        dataXL_C2(end+1) = [locsXL{2,i}];
        dataXR_C2(end+1) = [locsXR{2,i}];
        dataY(end+1) = [locsY{1,i}];
    end
end

% Output into Excel file (!! need to include unit conversion)
col_header = {'PD-axis(Y)','C1_XL','C1_XR','C2_XL','C2_XR'}; 
xlswrite(fpath,col_header,'Sheet1','A1');             %Write column header
xlswrite(fpath,[dataY' dataXL_C1' dataXR_C1' dataXL_C2' dataXR_C2'],'Sheet1','A2');        %Write data

