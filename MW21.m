% Analysis of PIV data
clear;

% :: Data Import
% =========================================================================
load('./PIVlab.mat');  % load the parameter values -- ** save the filename as "time point", may help to do it automatically
numFrames = length(velocity_magnitude);   % # of PIV images
% =========================================================================

%:: NEW
CMZy_ave = [];
%::

for t=1:numFrames
    %:: Key Parameters
    % =========================================================================
    %-- Xcomponent of velocity: u_filtered
    %-- Ycomponent of velocity: v_filtered
    %-- speed: velocity_magnitude
    %-- location: x & y in units of meter(?)

    [numRows,numCols] = size(velocity_magnitude{t}); % # of rows & columns in the velocity heat map 
    % =========================================================================

    %:: Analysis
    % =========================================================================
    % (1) Define the Region of Interest (ROI) <-- should be modified as CMZ region
    ROI_X1{t} = [1:numCols];  % X-range of the first wave
%    ROI_Y1{t} = [1:numRows];  % Y-range of the first wave
    % --- CMZ location ---
    imgTemp = imread('nanocor_0.tif');      % kymo graphs for the CMZ detection
    [imgHeight,~] = size(imgTemp); 
    imgCutInfo = importdata('cutEdge.txt');  % import the data of cutting info
    yCut = min(imgCutInfo(:,2));  % y-location of the cutting plane
    
    locCMZ = csvread('CMZ_detection.csv',1,0); % csvread(csv_file,row_offset,col_offset) 
    dX = (locCMZ(end,1)-locCMZ(1,1)) / numFrames;  % Xspacing in the CMZ kymos for individual PIV image  
    XL = round(1+dX*(t-1));  % Xrange of the PIV image to be analyzed
    XR = min(round(XL+dX), length(locCMZ));
    
    CMZy1 = mean(locCMZ(XL:XR,2))-yCut;  % bottom of the CMZ region
    CMZy2 = mean(locCMZ(XL:XR,3))-yCut;  % top of the CMZ region 
    imgHeight = imgHeight-yCut;          % adjust height due to cutting
    Y1 = ceil( CMZy1/imgHeight*numRows );  % CMZ bottom in the PIV coordinate
    Y2 = ceil( CMZy2/imgHeight*numRows );  % CMZ top ...
    ROI_Y1{t} = [Y2:Y1];  % Y-range of CMZ in the first wave
    
    %:: NEW
    CMZy_ave(t) = 0.5*(CMZy1+CMZy2)+yCut; % average location of CMZ in the Kymo coordinate
    %::
    % ---

    % (2) Find the max Velocity & its location within the ROI
    pct   = 10;  % top 10% of velocities in the ROI
    VList = reshape(velocity_magnitude{t}(ROI_Y1{t},ROI_X1{t}),1,[]);  % flatten matrix into a row vector
    VList_sort = sort(VList,'descend');  % sorting
    VList_ROI  = VList_sort( 1:ceil(length(VList)*pct*0.01) );   % find the Top ?% in the sorted list; data in the ROI

    Vmax{t} = max(VList_ROI);
    [Vmax_locY_ind,Vmax_locX_ind] = find(velocity_magnitude{t}==Vmax{t});
    Vmax_locY{t} = y{t}(Vmax_locY_ind,1);  % X location of max V in real units
    Vmax_locX{t} = x{t}(1,Vmax_locX_ind);  % Y location ...
    % (optional) the average velocity in the ROI
    Vave{t} = mean(VList_ROI);
    % (optional) display the results
    fprintf('For T=%03d \n',t);
    fprintf('The max velocity within ROI is %e \n',Vmax{t});
    fprintf('its location is (%e, %e)\n\n',Vmax_locX{t},Vmax_locY{t});
    fprintf('The average velocity within ROI is %e \n\n',Vave{t});
    
    % (3) Find the Mean Velocity within the CMZ
    % (optional) Plotting histogram of volocities within CMZ
    % ---------------------------------------------------
    figure('visible','off');
    nbins = 25;
    imgHis = histogram(VList_sort,nbins);
    xlabel('Velocity');
    ylabel('Count');
    % --- save into a folder
    fpath = 'histogram_velocities';  % folder name
    mkdir(fpath);
    FIGpath = sprintf([fpath '\\%03d'],t);  % figure path
    saveas(imgHis,FIGpath,'jpeg');
    % ---------------------------------------------------
    
    % identification of Region 0f Interest (ROI) "without outliers"
    pctRov = 10;  % percentage to be removed from the velocities distribution; !! 0 < pctRov < 50
    numRov = ceil(length(VList)*pctRov*0.01);  % number of velocities to be removed
    
    VList_ROI_Rov  = VList_sort( (numRov+1):(length(VList)-numRov) );   % Rmove the top & bottom ? percentage in the sorted list; data in the ROI
    Vave_Rov{t} = mean(VList_ROI_Rov);   % mean velocity in the ROI without outliers
    % =========================================================================
    
    
    % (optional) for check purpose: CMZ location
    % =========================================================================
    figure('visible','off');
    fname = sprintf('PIVlab_out_%03d.jpg',t);
    img = imread(fname);
    [imgRows,imgCols,~] = size(img);
    % --- becareful the location
    reY1 = CMZy1/imgHeight*imgRows;
    reY2 = CMZy2/imgHeight*imgRows;
    % ---   
    img_rev =  insertShape(img,'Line',[1,reY1,imgCols,reY1],'LineWidth',2,'Color','red'); % bottom line
    img_rev2 = insertShape(img_rev,'Line',[1,reY2,imgCols,reY2],'LineWidth',2,'Color','green');  % top line
    imgTT = imshow(img_rev2);
    
    % --- save into a folder
    fpath = 'check_CMZLocation';  % folder name
    mkdir(fpath);
    FIGpath = sprintf([fpath '\\%03d'],t);  % figure path
    saveas(imgTT,FIGpath,'jpeg');
    % =========================================================================    
    
    % (optional) for check purpose: max Velocity location
    % =========================================================================
    figure('visible','off');
    fname = sprintf('PIVlab_out_%03d.jpg',t);
    img = imread(fname);
    [imgRows,imgCols,~] = size(img);
    % --- becareful the location
    reX = Vmax_locX_ind / numCols;
    reY = Vmax_locY_ind / numRows;
    pos = [round(reX*imgCols),round(reY*imgRows)];  % location of the max V
    % ---
    img_rev = insertMarker(img,pos,'x','color','magenta','size',20);
    imgJJ = imshow(img_rev); 
    
     % --- save into a folder
    fpath = 'check_vMaxLocation';  % folder name
    mkdir(fpath);
    FIGpath = sprintf([fpath '\\%03d'],t);  % figure path
    saveas(imgJJ,FIGpath,'jpeg');   
    % =========================================================================
    
    %::(Optional) Replot the heat map
    % =========================================================================
%     surf(x{t},y{t},velocity_magnitude{t});  % plot the heat map
%     view(2);  % Project onto XY plane
%     set(gca,'Ydir','reverse');  % reverse the Y-axis
%     figure;
%     imshow(img_rev);
    % =========================================================================
end

%%
%:: Data Output
% =========================================================================
header = {'Time', 'Max_V', 'Loc_X', 'Loc_Y', 'Ave_V', 'Ave_V_cutted', 'AP_axis'};
dataOut = table([1:numFrames]', Vmax', Vmax_locX', Vmax_locY', Vave', Vave_Rov', CMZy_ave');
dataOut.Properties.VariableNames = header;  % adding header to the table data
writetable(dataOut,'PIV_analyzed.xlsx','Sheet',1);  % output data into excel
% =========================================================================