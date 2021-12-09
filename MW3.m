function ImgAnalysis(kymoFile,idx_file)

% Environment setting
%%
%--------------------------------------------------------------------------
fdname_ = regexp(kymoFile,'.tif','split');  % split filename into XXX and .tif
fdname = [fdname_{1} '_analyzed'];  % name of the folder for data output
fpath = sprintf('./%s/',fdname);    % file path of the folder
if ~exist(fpath, 'dir')
    mkdir(fpath);                   % create the folder -- for saving figs
end
% create new sub-folders to save the current analysis results
folder_ = dir(fpath);
LisFold = [folder_.isdir];
numFolders = length(LisFold(LisFold))-2;  % skip './' & '../'
fpath = [fpath 'Measure' num2str(numFolders+1) '/'];
mkdir(fpath);
%--------------------------------------------------------------------------
%%

% Data Import
%%
%--------------------------------------------------------------------------
[img] = imread(kymoFile);
% imshow(img) % show image
% histogram(img) % show pixel distribution
%--------------------------------------------------------------------------
%%

% Data Pre-processing
%%
%--------------------------------------------------------------------------
% Smooth - removing background noise
img_sm = imgaussfilt(img,2);  % Gaussian smoothing filter
imwrite(img_sm, [sprintf('%s',fpath) 'img_sm.jpeg']);  % save fig

% Binarization
% (method 1)
pixelMax = max(max(img_sm));   % max pixel
pixelThre = round(0.75*pixelMax );   % pixel threshold  0.75*pixelMax
img_b = img_sm;   % copy image, to keep the original one
img_b(img_b > pixelThre) = 255;     % binarized image
img_b(img_b < pixelThre) = 0;
imwrite(img_b, [sprintf('%s',fpath) 'img_b.jpeg']);  % save fig
%(method 2)
% img_sm = mat2gray(img_sm); % tranform into double-format & range [0-1]
% pixelThre = graythresh(img_sm(round(end/2):end,:));  % ignore the upper-half to avoid 3-parts
% img_b = im2bw(img_sm, pixelThre); 
% imwrite(img_b, [sprintf('%s',fpath) 'img_b.jpeg']);  % save fig

% erosion & dilation
se = strel('rectangle',[5 10]);     % structual element
img_b = imopen(img_b,se);   % remove small noise
%img_b = imerode(img_b,se); %(option)
se = strel('rectangle',[5 10]);
img_b = imclose(img_b,se); % fill holes inside object
%img_b = imdilate(img_b,se); %(option)
img_b(img_b > 0) = 255;     % binarized image; again
imwrite(img_b, [sprintf('%s',fpath) 'img_b_rev1.jpeg']);  % save fig

% segmentation
[Lobj,Numobj] = bwlabel(img_b,8);  % finding 'objects' in the image, Lobj=0 is the background

% ------------ analyze the size of each object ----------------------
list_sizeObj = [];
for i=1:Numobj
    list_sizeObj(end+1) = length(find(Lobj==i));  % list of object size 
end
% -------------------------------------------------------------------
% ----- determine the objects of interest based on size -------------
%sizeObj_c = prctile(list_sizeObj,80);       % threshold of the object size
sizeObj_c = 0.8*max(list_sizeObj);

list_Obj = find(list_sizeObj >= sizeObj_c);  % list of objects whose size > size_c
% -------------------------------------------------------------------
% ---- collect the pixel location of the objects of interest --------
pixeLoc = [];
for i=1:length(list_Obj)
    ind_Obj = list_Obj(i);  % the object index
    pixeLoc = cat(1, pixeLoc, find(Lobj==ind_Obj)); % object indexed by "ind_Obj"
end
% -------------------------------------------------------------------
img_b = 0*img_b;  % initialization
img_b(pixeLoc) = 255;  % the selected region
imwrite(img_b, [sprintf('%s',fpath) 'img_b_rev2.jpeg']);  % save fig

%(option) overlay two images
%---------------------------
% c = imfuse(img,img_b);
% imshow(c)
%---------------------------
%(option) image histogram
% imhist(img);
%--------------------------------------------------------------------------
%%

% Analysis
%--------------------------------------------------------------------------
x_tot = length(img_b(1,:));   % total # of horizontal pixels
y_tot = length(img_b(:,1));   % total # of vertical pixels 
DX = round(x_tot/80);   % range of horizontal pixels to be analyzed
DY = round(y_tot/80);   % range of vertical pixels to be analyzed

% extract Y location of high-intensity points in each X
for x=1:x_tot
    indY_CMZ{x} = find(img_b(:,x)==255); 
end

% find the max X of the CMZ zone, i.e., CMZ_xMax 
for i=x_tot:-1:1   
    CMZLen_temp = DY;  % critical size for detection of existence of CMZ
    if length(indY_CMZ{i})>CMZLen_temp
        CMZ_xMax = i;
        break;  % escape loop
    end
    CMZ_xMax = i;
end
%%
%::: (1) Front Wave :::
% extract points of the front wave
% clearing - step1
CMZ_fx = []; % initialization; x-location
CMZ_fy = []; % y-location
for x=1:x_tot
    CMZ_fx(end+1) = x; 
    CMZ_fy(end+1) = prctile(indY_CMZ{x},98);  % extract the 98% position
end
figure('visible','off'); % turn off figure
imshow(img_b) 
hold on;
plot(CMZ_fx,CMZ_fy,'ro'); % plot the points to be fitted
hold off;
saveas(gca,fullfile([fpath 'img_fit_frontWave1_step1.jpeg'])); % show front points to be fitted with clearing-step1 only
% clearing - step2
CMZy1 = prctile(CMZ_fy,10);  % 10% point of the y-location distribution
CMZy2 = prctile(CMZ_fy,90);  % 90% point of ...
ind_y = find( (CMZy1<CMZ_fy) & (CMZ_fy<CMZy2) );  % show the y-location with CMZy1~CMZy2
figure('visible','off'); % turn off figure
imshow(img_b)
hold on;
plot(CMZ_fx(ind_y), CMZ_fy(ind_y),'ro'); % plot the points to be fitted
hold off;
saveas(gca,fullfile([fpath 'img_fit_frontWave1_step2.jpeg'])); % show front points to be fitted with clearing-step2
% Least-squared Linear fitting
coefit = polyfit(CMZ_fx(ind_y), CMZ_fy(ind_y), 1);
% visualization + manual correction
%----------------------------------
answer = 'No';
coefit_(1) = coefit(1);  % backup the original fitting value
coefit_(2) = coefit(2);
numMod_(1) = 0.0;   % initialization
numMod_(2) = 0.0;
while strcmp(answer,'No')   % continue to correct, if not satisfied
    close(gcf);   % close image
    figure('visible','on');
    set(gcf, 'Position',  [50, 150, 1200, 500]);
    subplot(1,2,1);
    imshow(img) % should overlay with the original figure
    hold on;
    plot(CMZ_fx(ind_y), CMZ_fy(ind_y),'ro','MarkerSize',2); % plot the points to be fitted
    plot(CMZ_fx, coefit(1)*CMZ_fx+coefit(2), 'g-', 'Linewidth', 1.5);    % plot new fitting line
    plot(CMZ_fx, coefit_(1)*CMZ_fx+coefit_(2), 'b--', 'Linewidth', 1.5); % plot original fitting line
    subplot(1,2,2);
    imshow(img_b) % should overlay with the original figure
    hold on;
    plot(CMZ_fx(ind_y), CMZ_fy(ind_y),'ro','MarkerSize',2); % plot the points to be fitted
    plot(CMZ_fx, coefit(1)*CMZ_fx+coefit(2), 'g-', 'Linewidth', 1.5);    % plot new fitting line
    plot(CMZ_fx, coefit_(1)*CMZ_fx+coefit_(2), 'b--', 'Linewidth', 1.5); % plot original fitting line
    hold off;
    
    quest = sprintf('Is frontWave1 OK?');
    set(groot,'defaultUicontrolFontSize', 16);  % change font size
    answer = questdlg(quest,'Fitting','Yes','No','Yes');
 
    if strcmp(answer,'No')  % correction
        prompt = {['\fontsize{12} Slope correction X:' sprintf('\n') ' (1+X) * Slope0  |  ' sprintf('Slope0=%.3f',coefit_(1))], ...
                  ['\fontsize{12} Intercept correction Y:' sprintf('\n') ' (1+Y) * Intercept0  |  ' sprintf('Intercept0=%.3f',coefit_(2))]};
        dlgtitle = sprintf('X_{pixels} = %d ; \t Y_{pixels} = %d',x_tot,y_tot);
        dims = [2 80];  % dialog box size
        definput = {num2str(numMod_(1)),num2str(numMod_(2))};  % default values
        opts.Interpreter = 'tex';  % need to switch to "tex" mode, in order to change font size
        answer2 = inputdlg(prompt,dlgtitle,dims,definput,opts);       
        coefit(1)= coefit_(1)*(1+str2num(answer2{1})); % slope correction
        coefit(2)= coefit_(2)*(1+str2num(answer2{2})); % intercept correction
        numMod_(1) = str2num(answer2{1});  % save the previous input
        numMod_(2) = str2num(answer2{2});
    end
end
saveas(gca,fullfile([fpath 'img_fit_frontWave1.jpeg'])); % show fitting line with images
close(gcf);   % close image
%----------------------------------
% Wave front speed
frontSpeed = coefit(1);
frontSpeed_intercept = coefit(2);
%%
%::: (2) Wave back of the 1st wave ::: 
% check whether 2nd wave exists
yHigh = mean(prctile([indY_CMZ{1:DX}],95)); % upper bound of Y exploration range for the 2nd wave
%----- (choose one) method 1
% Ytemp = []; Ytemp_ = []; 
% Ytemp_ = mean2(img_sm(1:DY,1:DX));  % find the average of background pixel
% for i=1:DX
%     Ytemp(end+1)= max( find(img_sm(:,i)<(1.1*Ytemp_)) ); % find max pixel position < (1.1)*background
% end
% Y_t0 = mean(Ytemp);   % Y-location @ t=0
%----- 
%----- (choose one) method 2
%(old version)
%Y_t0 = mean(prctile([indY_CMZ{1:DX}],3));
%(new version)
x_ind = find(~cellfun(@isempty, indY_CMZ));  % find x location of nonzero CMZ
Y_t0 = mean(prctile([indY_CMZ{x_ind(1)}],3));
%-----
yLow  = Y_t0;  % lower bound of ... 
CMZLen_temp = 1.5*DY;  % critical size for detection of CMZ zone existence 
Y_indCMZ_ = [];
ind_temp = [];
ind_2wave = [];   % index for detection of 2nd wave existence
for i=(CMZ_xMax-10*DX):2:CMZ_xMax   % explore DX-range from the X_end of CMZ
    yIndx_ = find(yLow<[indY_CMZ{i}] & [indY_CMZ{i}]<yHigh); % index for pixels in the ROI
    if length(yIndx_)>CMZLen_temp  
        ind_temp(end+1) = 1.0;  % CMZ zone exists
    else 
        ind_temp(end+1) = 0.0;  % CMZ doesn't exist
    end
end
if mean(ind_temp)>0.5  % majority vote
    ind_2wave = 1.0;   % 2nd wave, YES!
else
    ind_2wave = 0.0;   % 2nd wave, NO!
end

%:: Manual check ::
figure('visible','on');
set(gcf, 'Position', [50, 150, 1200, 500]);
subplot(1,2,1);
imshow(img); 
subplot(1,2,2);
imshow(img_b);
switch ind_2wave
    case 1.0
        quest = sprintf('2-wave analysis, OK?');
    case 0.0
        quest = sprintf('1-wave analysis, OK?');
end
set(groot,'defaultUicontrolFontSize', 16);  % change font size
answer = questdlg(quest,'1- or 2-wave','Yes','No','Yes');
if strcmp(answer,'No')  % compare two strings
    ind_2wave = abs(ind_2wave-1.0); % switch the detection result
end
close(gcf);   % close image
%::::::::::::::::::

% choose 1-wave or 2-wave analysis
switch ind_2wave
    case 1.0    
        display('2-wave analysis');
        %----------------------------------------------- 
        %::: (3-1) Wave back of the 1st wave @ 2-wave case :::
        diffY_indCMZ = [];   % !!valid only for existence of 2-waves
        for i=1:x_tot  % save CMZ pixel location between 1- and 2-waves
            diffY_indCMZ{i} = [indY_CMZ{i}(2:end)-indY_CMZ{i}(1:end-1)]; % given X, calculate DY between two successive CMZ pixel
        end 
        diffYC = round(y_tot/20);  % minimal DY for detection betwee 1- and 2-waves
        indTemp = [];              
        indY_CMZ_w12 = [];         % save X-location for appearance of boundary between 1-2 waves
        for i=1:x_tot
            indTemp = find( diffY_indCMZ{i}>diffYC ); % record pixels when pixel spacing is large enough
            listTemp = [];
            for j=1:length(indTemp)
                listTemp(end+1) = indY_CMZ{i}(indTemp(j)+1); % pixel location of one of the pixel pair
            end
            indY_CMZ_w12{i} = listTemp;
        end
        figure('visible','off');
        imshow(img_b) 
        hold on;
        for i=1:x_tot
            plot(i*ones(length(indY_CMZ_w12{i})), [indY_CMZ_w12{i}],'ro','LineWidth', 1.5); % plot the line fitted
        end
        hold off;
        saveas(gca,fullfile([fpath 'img_fit_backWave1_step1.jpeg'])); % save fig
        % clearing - step2
        yloc1 = prctile([indY_CMZ_w12{:}],15); % y-location of the 15% percentile
        yloc2 = prctile([indY_CMZ_w12{:}],75); % y-location of the 75% percentile
        xList = [];
        yList = [];
        for i=1:x_tot
            temp = find( (yloc1<indY_CMZ_w12{i})&(indY_CMZ_w12{i}<yloc2) );
            if temp
                xList(end+1) = i;
                yList(end+1) = mean(indY_CMZ_w12{i}(temp));  % take average over all points at fixed X
            end
        end        
        figure('visible','off');
        imshow(img_b) 
        hold on;
        for i=1:x_tot
            plot(xList, yList,'ro','LineWidth', 1.5); 
        end
        hold off;
        saveas(gca,fullfile([fpath 'img_fit_backWave1_step2.jpeg'])); % save fig
        % Least-squared Linear fitting
        coefit = polyfit(xList, yList, 1);
        % visualization + manual correction
        %----------------------------------
        answer = 'No';
        coefit_(1) = coefit(1);  % backup the original values
        coefit_(2) = coefit(2); 
        numMod_(1) = 0.0;   % initialization
        numMod_(2) = 0.0;
        while strcmp(answer,'No')   % continue to correct, if not satisfied
            close(gcf);   % close image           
            figure('visible','on');
            set(gcf, 'Position',  [50, 150, 1200, 500]);
            subplot(1,2,1);
            imshow(img) % should overlay with the original figure
            hold on;
            plot(xList, yList,'ro'); % plot the points to be fitted
            plot(CMZ_fx, coefit(1)*CMZ_fx+coefit(2), 'g-', 'Linewidth', 1.5); % plot new fitting line
            plot(CMZ_fx, coefit_(1)*CMZ_fx+coefit_(2), 'b--', 'Linewidth', 1.5); % plot original fitting line
            subplot(1,2,2);
            imshow(img_b) % should overlay with the original figure
            hold on;
            plot(xList, yList,'ro'); % plot the points to be fitted
            plot(CMZ_fx, coefit(1)*CMZ_fx+coefit(2), 'g-', 'Linewidth', 1.5); % plot new fitting line
            plot(CMZ_fx, coefit_(1)*CMZ_fx+coefit_(2), 'b--', 'Linewidth', 1.5); % plot original fitting line
            hold off;
            
            quest = sprintf('Is backWave1 OK?');
            set(groot,'defaultUicontrolFontSize', 16);  % change font size
            answer = questdlg(quest,'Analysis','Yes','No','Yes');
            
            if strcmp(answer,'No')  % correction
                prompt = {['\fontsize{12} Slope correction X:' sprintf('\n') ' (1+X) * Slope0  |  ' sprintf('Slope0=%.3f',coefit_(1))], ...
                          ['\fontsize{12} Intercept correction Y:' sprintf('\n') ' (1+Y) * Intercept0  |  ' sprintf('Intercept0=%.3f',coefit_(2))]};
                dlgtitle = sprintf('X_{pixels} = %d ; \t Y_{pixels} = %d',x_tot,y_tot);
                dims = [2 80];  % dialog box size
                definput = {num2str(numMod_(1)),num2str(numMod_(2))};  % default values
                opts.Interpreter = 'tex';  % need to switch to "tex" mode, in order to change font size
                answer2 = inputdlg(prompt,dlgtitle,dims,definput,opts);       
                coefit(1)= coefit_(1)*(1+str2num(answer2{1})); % slope correction
                coefit(2)= coefit_(2)*(1+str2num(answer2{2})); % intercept correction
                numMod_(1) = str2num(answer2{1});  % save the previous input
                numMod_(2) = str2num(answer2{2});
            end 
        end
        saveas(gca,fullfile([fpath 'img_fit_backWave1.jpeg'])); % show fitting line with images
        close(gcf);   % close image
        % Wave back speed
        backSpeed = coefit(1); 
        backSpeed_intercept = coefit(2); 
        %::: (3-2) Wave front of the 2nd wave :::
        for i=1:x_tot
            indTemp = find( diffY_indCMZ{i}>diffYC ); % record pixels when pixel spacing is large enough
            listTemp = [];
            for j=1:length(indTemp)
                listTemp(end+1) = indY_CMZ{i}(indTemp(j)); % pixel location of the one of the pixel pair
            end
            indY_CMZ_w12{i} = listTemp;
        end
        figure('visible','off');
        imshow(img_b) 
        hold on;
        for i=1:x_tot
            plot(i*ones(length(indY_CMZ_w12{i})), [indY_CMZ_w12{i}],'ro','LineWidth', 1.5); % plot the line fitted
        end
        hold off;
        saveas(gca,fullfile([fpath 'img_fit_frontWave2_step1.jpeg']));
        % clearing - step2
        yloc1 = prctile([indY_CMZ_w12{:}],15); % y-location of the 15% percentile
        yloc2 = prctile([indY_CMZ_w12{:}],75); % y-location of the 75% percentile
        xList = [];
        yList = [];
        for i=1:x_tot
            temp = find( (yloc1<indY_CMZ_w12{i})&(indY_CMZ_w12{i}<yloc2) );
            if temp
                xList(end+1) = i;
                yList(end+1) = mean(indY_CMZ_w12{i}(temp));  % take average over all points at fixed X
            end
        end
        figure('visible','off');
        imshow(img_b) 
        hold on;
        for i=1:x_tot
            plot(xList, yList,'ro','LineWidth', 1.5); 
        end
        hold off;
        saveas(gca,fullfile([fpath 'img_fit_frontWave2_step2.jpeg']));
        % Least-squared Linear fitting
        coefit = polyfit(xList, yList, 1);
        % visualization + manual correction
        %----------------------------------
        answer = 'No';
        coefit_(1) = coefit(1);  % backup the original fitting value
        coefit_(2) = coefit(2);  
        numMod_(1) = 0.0;   % initialization
        numMod_(2) = 0.0;
        while strcmp(answer,'No')   % continue to correct, if not satisfied
            close(gcf);   % close image 
            figure('visible','on');
            set(gcf, 'Position',  [50, 150, 1200, 500]);
            subplot(1,2,1);
            imshow(img) % should overlay with the original figure
            hold on;
            plot(xList, yList,'ro'); % plot the points to be fitted
            plot(CMZ_fx, coefit(1)*CMZ_fx+coefit(2), 'g-', 'Linewidth', 1.5); % plot new fitting line
            plot(CMZ_fx, coefit_(1)*CMZ_fx+coefit_(2), 'b--', 'Linewidth', 1.5); % plot original fitting line
            subplot(1,2,2);
            imshow(img_b) % should overlay with the original figure
            hold on;
            plot(xList, yList,'ro'); % plot the points to be fitted
            plot(CMZ_fx, coefit(1)*CMZ_fx+coefit(2), 'g-', 'Linewidth', 1.5); % plot new fitting line
            plot(CMZ_fx, coefit_(1)*CMZ_fx+coefit_(2), 'b--', 'Linewidth', 1.5); % plot original fitting line
            hold off;
            
            quest = sprintf('Is frontWave2 OK?');
            set(groot,'defaultUicontrolFontSize', 16);  % change font size
            answer = questdlg(quest,'Fitting','Yes','No','Yes');            
            
            if strcmp(answer,'No')  % correction
                prompt = {['\fontsize{12} Slope correction X:' sprintf('\n') ' (1+X) * Slope0  |  ' sprintf('Slope0=%.3f',coefit_(1))], ...
                          ['\fontsize{12} Intercept correction Y:' sprintf('\n') ' (1+Y) * Intercept0  |  ' sprintf('Intercept0=%.3f',coefit_(2))]};
                dlgtitle = sprintf('X_{pixels} = %d ; \t Y_{pixels} = %d',x_tot,y_tot);
                dims = [2 80];  % dialog box size
                definput = {num2str(numMod_(1)),num2str(numMod_(2))};  % default values
                opts.Interpreter = 'tex';  % need to switch to "tex" mode, in order to change font size
                answer2 = inputdlg(prompt,dlgtitle,dims,definput,opts);       
                coefit(1)= coefit_(1)*(1+str2num(answer2{1})); % slope correction
                coefit(2)= coefit_(2)*(1+str2num(answer2{2})); % intercept correction
                numMod_(1) = str2num(answer2{1});  % save the previous input
                numMod_(2) = str2num(answer2{2});              
            end
        end
        saveas(gca,fullfile([fpath 'img_fit_frontWave2.jpeg'])); % show fitting line with images
        close(gcf);   % close image
        % Wave back speed
        frontSpeed2 = coefit(1);  
        frontSpeed2_intercept = coefit(2); 
        %-----------------------------------------------
    case 0.0
        display('1-wave analysis');
        %-----------------------------------------------
        %::: (3-1) Wave back of the 1st wave @ 1-wave case :::
        % extract points of CMZ close to the top
        CMZ_fx = []; % initialization; x-location
        CMZ_fy = []; % y-location
        x_start = round(x_tot/2); % X-range for scanning
        x_end = x_tot;   %CMZ_xMax;
        for x=x_start:x_end   
            CMZ_fx(end+1) = x; 
            CMZ_fy(end+1) = prctile(indY_CMZ{x},3);  % extract the 3% position
        end
        figure('visible','off');
        imshow(img) 
        hold on;
        plot(CMZ_fx,CMZ_fy,'ro'); % plot the points to be fitted
        hold off;
        saveas(gca,fullfile([fpath 'img_fit_backWave1_step1.jpeg'])); % show points to be fitted               
        % clearing - step2 
        yloc1 = prctile(CMZ_fy,10); % y-location of the 10% percentile
        yloc2 = prctile(CMZ_fy,90); % y-location of the 90% percentile
        xList = [];
        yList = [];       
        temp = find( (yloc1<CMZ_fy)&(CMZ_fy<yloc2) ); % index of [15%~75%] data points
        xList = CMZ_fx(temp);
        yList = CMZ_fy(temp);  % keep [15%~75%] data points only       
        figure('visible','off');
        imshow(img_b) 
        hold on;
        plot(xList, yList,'ro','LineWidth', 1.5); 
        hold off;
        saveas(gca,fullfile([fpath 'img_fit_backWave1_step2.jpeg']));
        % Least-squared Linear fitting
        coefit = polyfit(xList, yList, 1);
        % visualization + manual correction
        %----------------------------------
        answer = 'No';
        coefit_(1) = coefit(1);  % backup the original fitting value
        coefit_(2) = coefit(2); 
        numMod_(1) = 0.0;   % initialization
        numMod_(2) = 0.0;
        while strcmp(answer,'No')   % continue to correct, if not satisfied
            close(gcf);   % close image        
            CMZ_fx = linspace(1,x_tot,20); % for plotting purpose
            figure('visible','on');
            set(gcf, 'Position',  [50, 150, 1200, 500]);
            subplot(1,2,1);
            imshow(img) % should overlay with the original figure
            hold on;
            plot(xList, yList,'ro','MarkerSize',2); % plot the points to be fitted
            plot(CMZ_fx, coefit(1)*CMZ_fx+coefit(2), 'g-', 'Linewidth', 1.5); % plot new fitting line
            plot(CMZ_fx, coefit_(1)*CMZ_fx+coefit_(2), 'b--', 'Linewidth', 1.5); % plot original fitting line
            subplot(1,2,2);
            imshow(img_b) % should overlay with the original figure
            hold on;
            plot(xList, yList,'ro','MarkerSize',2); % plot the points to be fitted
            plot(CMZ_fx, coefit(1)*CMZ_fx+coefit(2), 'g-', 'Linewidth', 1.5); % plot new fitting line
            plot(CMZ_fx, coefit_(1)*CMZ_fx+coefit_(2), 'b--', 'Linewidth', 1.5); % plot original fitting line
            hold off;
        
            quest = sprintf('Is backWave1 OK?');
            set(groot,'defaultUicontrolFontSize', 16);  % change font size
            answer = questdlg(quest,'Fitting','Yes','No','Yes'); 
            
            if strcmp(answer,'No')  % correction
                prompt = {['\fontsize{12} Slope correction X:' sprintf('\n') ' (1+X) * Slope0  |  ' sprintf('Slope0=%.3f',coefit_(1))], ...
                          ['\fontsize{12} Intercept correction Y:' sprintf('\n') ' (1+Y) * Intercept0  |  ' sprintf('Intercept0=%.3f',coefit_(2))]};
                dlgtitle = sprintf('X_{pixels} = %d ; \t Y_{pixels} = %d',x_tot,y_tot);
                dims = [2 80];  % dialog box size
                definput = {num2str(numMod_(1)),num2str(numMod_(2))};  % default values
                opts.Interpreter = 'tex';  % need to switch to "tex" mode, in order to change font size
                answer2 = inputdlg(prompt,dlgtitle,dims,definput,opts);       
                coefit(1)= coefit_(1)*(1+str2num(answer2{1})); % slope correction
                coefit(2)= coefit_(2)*(1+str2num(answer2{2})); % intercept correction
                numMod_(1) = str2num(answer2{1});  % save the previous input
                numMod_(2) = str2num(answer2{2});                
            end
        end
        saveas(gca,fullfile([fpath 'img_fit_backWave1.jpeg'])); % show fitting line with images
        close(gcf);   % close image
        % Wave back speed
        backSpeed = coefit(1); 
        backSpeed_intercept = coefit(2);         
        %-----------------------------------------------
end
%%
%::: (3) Time points- interception of fitting lines :::
%>> point #1
slope1 = 0;  % properties of Line1
intercept1 = Y_t0;
slope2 = backSpeed;  % properties of Line2
intercept2 = backSpeed_intercept;
% interception for detection of max CMZ length
t_maxCMZlen = (intercept2-intercept1)/(slope1-slope2); % will be used for detection of max CMZ
% visualization + manual correction
%----------------------------------
answer = 'No';
Y_t0_ = Y_t0; % backup the original values
numMod_(1) = 0.0;   % initialization
while strcmp(answer,'No')   % continue to correct, if not satisfied
    close(gcf);   % close image
    figure('visible','on');
    set(gcf, 'Position',  [50, 150, 1200, 500]);
    subplot(1,2,1);
    imshow(img) % should overlay with the original figure
    hold on;
    plot(CMZ_fx, slope1*CMZ_fx+intercept1, 'b-', 'Linewidth', 1.5); % plot line 1
    plot(CMZ_fx, slope2*CMZ_fx+intercept2, 'b-', 'Linewidth', 1.5); % plot line 2
    plot(t_maxCMZlen, slope1*t_maxCMZlen+intercept1, 'ro','Linewidth', 1.5,'MarkerSize',10);   % new intercept point
    subplot(1,2,2);
    imshow(img_b) % should overlay with the original figure
    hold on;
    plot(CMZ_fx, slope1*CMZ_fx+intercept1, 'b-', 'Linewidth', 1.5); % plot line 1
    plot(CMZ_fx, slope2*CMZ_fx+intercept2, 'b-', 'Linewidth', 1.5); % plot line 2
    plot(t_maxCMZlen, slope1*t_maxCMZlen+intercept1, 'ro','Linewidth', 1.5,'MarkerSize',10);   % new intercept point
    hold off;
    
    quest = sprintf('Is T_maxCMZ OK?');
    set(groot,'defaultUicontrolFontSize', 16);  % change font size
    answer = questdlg(quest,'Analysis','Yes','No','Yes');
     
    if strcmp(answer,'No')  % correction
        prompt = {['\fontsize{12} Y position correction on the top line:' sprintf('\n') ' (1+Y)*Y0  |  ' sprintf('Y0=%.3f',Y_t0)]};
        dlgtitle = sprintf('X_{pixels} = %d ; \t Y_{pixels} = %d',x_tot,y_tot);
        dims = [2 80];  % dialog box size
        definput = {num2str(numMod_(1))};  % default values
        opts.Interpreter = 'tex';  % need to switch to "tex" mode, in order to change font size
        answer2 = inputdlg(prompt,dlgtitle,dims,definput,opts);     
        Y_t0 = Y_t0_*(1+str2num(answer2{1})); % Y position correction
        intercept1 = Y_t0;
        t_maxCMZlen = (intercept2-intercept1)/(slope1-slope2);
        numMod_(1) = str2num(answer2{1});  % save the previous input
    end
end
saveas(gca,fullfile([fpath 'img_TmaxCMZlen.jpeg'])); 
close(gcf);
%>> point #2
slope1 = frontSpeed;  % properties of Line1
intercept1 = frontSpeed_intercept;
slope2 = backSpeed;  % properties of Line2
intercept2 = backSpeed_intercept;
% interception for detection of CMZ vanishing
t_CMZvanish = (intercept2-intercept1)/(slope1-slope2); 
% visualization + manual correction
%----------------------------------
answer = 'No';
t_CMZvanish_ = t_CMZvanish;  % backup the original values
numMod_(1) = 0.0;   % initialization
while strcmp(answer,'No')   % continue to correct, if not satisfied
    close(gcf);   % close image
    figure('visible','on');
    set(gcf, 'Position',  [50, 150, 1200, 500]);
    subplot(1,2,1);
    imshow(img) % should overlay with the original figure
    hold on;
    plot(CMZ_fx, slope1*CMZ_fx+intercept1, 'b-', 'Linewidth', 1.5); % plot line 1
    plot(CMZ_fx, slope2*CMZ_fx+intercept2, 'b-', 'Linewidth', 1.5); % plot line 2
    plot(t_CMZvanish, slope1*t_CMZvanish+intercept1, 'go','Linewidth', 1.5,'MarkerSize',10); % new intercept point
    plot(t_CMZvanish_, slope1*t_CMZvanish_+intercept1, 'ro','Linewidth', 1.5,'MarkerSize',10); % original intercept point
    subplot(1,2,2);
    imshow(img_b) % should overlay with the original figure
    hold on;
    plot(CMZ_fx, slope1*CMZ_fx+intercept1, 'b-', 'Linewidth', 1.5); % plot line 1
    plot(CMZ_fx, slope2*CMZ_fx+intercept2, 'b-', 'Linewidth', 1.5); % plot line 2
    plot(t_CMZvanish, slope1*t_CMZvanish+intercept1, 'go','Linewidth', 1.5,'MarkerSize',10); % new intercept point
    plot(t_CMZvanish_, slope1*t_CMZvanish_+intercept1, 'ro','Linewidth', 1.5,'MarkerSize',10);
    hold off;
    
    quest = sprintf('Is T_vanish OK?');
    set(groot,'defaultUicontrolFontSize', 16);  % change font size
    answer = questdlg(quest,'Analysis','Yes','No','Yes');
    
    if strcmp(answer,'No')  % correction
        prompt = {['\fontsize{12} X position correction:' sprintf('\n') ' (1+X)*X0  |  ' sprintf('X0=%.3f',t_CMZvanish_)]};
        dlgtitle = sprintf('X_{pixels} = %d ; \t Y_{pixels} = %d',x_tot,y_tot);
        dims = [2 80];  % dialog box size
        definput = {num2str(numMod_(1))};  % default values
        opts.Interpreter = 'tex';  % need to switch to "tex" mode, in order to change font size
        answer2 = inputdlg(prompt,dlgtitle,dims,definput,opts);       
        t_CMZvanish= t_CMZvanish_*(1+str2num(answer2{1})); % X position correction
        numMod_(1) = str2num(answer2{1});  % save the previous input
    end 
end
saveas(gca,fullfile([fpath 'img_TCMZvanish.jpeg'])); 
close(gcf);
%%
%::: (4) Max CMZ distance traveled :::
%%
Y_tf = frontSpeed*t_CMZvanish+frontSpeed_intercept;  % Y @ t=CMZvanish
% visualization + manual correction
%----------------------------------
answer = 'No';
Y_t0_ = Y_t0;  % backup the original values
Y_tf_ = Y_tf;
numMod_(1) = 0.0;   % initialization
numMod_(2) = 0.0;
while strcmp(answer,'No')   % continue to correct, if not satisfied
    close(gcf);   % close image
    figure('visible','on');
    set(gcf, 'Position',  [50, 150, 1200, 500]);
    subplot(1,2,1);
    imshow(img) 
    hold on;
    plot([2*DX,2*DX],[Y_t0, Y_tf],'ro-','LineWidth', 1.5); % plot the line fitted
    plot([1,x_tot],Y_t0*ones(2),'g-','LineWidth',1.5); % plot the new top horizontal line
    plot([1,x_tot],Y_tf*ones(2),'g-','LineWidth',1.5); % plot the new bottom horizontal line
    plot([1,x_tot],Y_t0_*ones(2),'b--','LineWidth',1.5); % plot the original top horizontal line
    plot([1,x_tot],Y_tf_*ones(2),'b--','LineWidth',1.5); % plot the original bottom horizontal line
    subplot(1,2,2);
    imshow(img_b) 
    hold on;
    plot([2*DX,2*DX],[Y_t0, Y_tf],'ro-','LineWidth', 1.5); % plot the line fitted
    plot([1,x_tot],Y_t0*ones(2),'g-','LineWidth',1.5); % plot the new top horizontal line
    plot([1,x_tot],Y_tf*ones(2),'g-','LineWidth',1.5); % plot the new bottom horizontal line
    plot([1,x_tot],Y_t0_*ones(2),'b--','LineWidth',1.5); % plot the original top horizontal line
    plot([1,x_tot],Y_tf_*ones(2),'b--','LineWidth',1.5); % plot the original bottom horizontal line
    hold off;    

    quest = sprintf('Is maxCMZDist OK?');
    set(groot,'defaultUicontrolFontSize', 16);  % change font size
    answer = questdlg(quest,'Analysis','Yes','No','Yes');
    
    if strcmp(answer,'No')  % correction
        prompt = {['\fontsize{12} Top point correction Y1:' sprintf('\n') ' (1+Y1)*Y10  |  ' sprintf('Y10=%.3f',Y_t0_)], ...
                  ['\fontsize{12} Bottom point correction Y2:' sprintf('\n') ' (1+Y2)*Y20  |  ' sprintf('Y20=%.3f',Y_tf_)]};
        dlgtitle = sprintf('X_{pixels} = %d ; \t Y_{pixels} = %d',x_tot,y_tot);
        dims = [2 80];  % dialog box size
        definput = {num2str(numMod_(1)),num2str(numMod_(2))};  % default values
        opts.Interpreter = 'tex';  % need to switch to "tex" mode, in order to change font size
        answer2 = inputdlg(prompt,dlgtitle,dims,definput,opts);       
        Y_t0= Y_t0_*(1+str2num(answer2{1})); % slope correction
        Y_tf= Y_tf_*(1+str2num(answer2{2})); % intercept correction
        numMod_(1) = str2num(answer2{1});  % save the previous input
        numMod_(2) = str2num(answer2{2});
    end 
end
saveas(gca,fullfile([fpath 'img_fit_maxCMZDist.jpeg'])); % save fig
close(gcf);   % close image
%----------------------------------
% Max CMZ distance traveled
maxCMZDist = Y_tf-Y_t0;
%%
%::: (5) Max CMZ length :::
% extract points of CMZ close to the top
CMZ_fx2 = []; % initialization; x-location
CMZ_fy2 = []; % y-location of 1st backwave
CMZ_fy  = []; % y-location of 1st frontwave
for x=1:x_tot
    CMZ_fx2(end+1) = x; 
    CMZ_fy2(end+1) = prctile(indY_CMZ{x},3);  % extract the 3% position
    CMZ_fy(end+1) = prctile(indY_CMZ{x},95);  % extract the 95% position
end
figure('visible','off');
imshow(img_b) 
hold on;
plot(CMZ_fx2,CMZ_fy2,'ro'); % plot the points to be fitted
hold off;
saveas(gca,fullfile([fpath 'img_maxCMZ_step1.jpeg'])); % show points to be fitted
Xmid = min(abs(t_maxCMZlen),x_tot); % X-location below which, CMZ length is calculated
CMZLen = CMZ_fy(1:Xmid)-CMZ_fy2(1:Xmid);
maxCMZLen = max(CMZLen);  % max CMZ length
indX_CMZLen = find(CMZLen==max(CMZLen));   % X-location of the max CMZ length
% visualization + manual correction
%----------------------------------
answer = 'No';
maxCMZLen_ = maxCMZLen;  % backup the original values
scale = 0.0;   % initialization
shiftX = 0.0;
shiftY = 0.0;
while strcmp(answer,'No')   % continue to correct, if not satisfied
    Xtop_ = (1+shiftX)*indX_CMZLen(end);
    Xbot_ = (1+shiftX)*indX_CMZLen(end);
    Ycen_ = (1+shiftY)*0.5*(CMZ_fy(indX_CMZLen(end))+CMZ_fy2(indX_CMZLen(end))); % center point
    Ytop_ = Ycen_ - (1+scale)*(0.5*maxCMZLen_);
    Ybot_ = Ycen_ + (1+scale)*(0.5*maxCMZLen_);

    close(gcf);   % close image
    figure('visible','on');
    set(gcf, 'Position',  [50, 150, 1200, 500]);
    subplot(1,2,1);
    imshow(img) 
    hold on;
    plot([Xtop_,Xbot_], [Ytop_,Ybot_],'g-','Linewidth', 1.5); % plot new line fitted
    plot([indX_CMZLen(end),indX_CMZLen(end)], [CMZ_fy(indX_CMZLen(end)),CMZ_fy2(indX_CMZLen(end))],'r--','Linewidth', 1.5); % plot original line fitted
    subplot(1,2,2);
    imshow(img_b) 
    hold on;
    plot([Xtop_,Xbot_], [Ytop_,Ybot_],'g-','Linewidth', 1.5); % plot new line fitted
    plot([indX_CMZLen(end),indX_CMZLen(end)], [CMZ_fy(indX_CMZLen(end)),CMZ_fy2(indX_CMZLen(end))],'r--','Linewidth', 1.5);
    hold off;
    
    quest = sprintf('Is maxCMZLen OK?');
    set(groot,'defaultUicontrolFontSize', 16);  % change font size
    answer = questdlg(quest,'Analysis','Yes','No','Yes');
    
    if strcmp(answer,'No')  % correction
        prompt = {['\fontsize{12} Scale correction :' sprintf('\n') ' (1+Scale) * Length0  |  ' sprintf('Length0=%.2f',maxCMZLen_)], ...
                  ['\fontsize{12} X Shift correction on the line center:' sprintf('\n') ' (1+ShiftX) * Xposition0  |  ' sprintf('Xposition0=%d',indX_CMZLen(end))], ...
                  ['\fontsize{12} Y Shift correction on the line center:' sprintf('\n') ' (1+ShiftY) * Yposition0  |  ' sprintf('Yposition0=%.2f',0.5*(CMZ_fy(indX_CMZLen(end))+CMZ_fy2(indX_CMZLen(end))))]};
        dlgtitle = sprintf('X_{pixels} = %d ; \t Y_{pixels} = %d',x_tot,y_tot);
        dims = [2 80];  % dialog box size
        definput = {num2str(scale),num2str(shiftX),num2str(shiftY)};  % default values
        opts.Interpreter = 'tex';  % need to switch to "tex" mode, in order to change font size
        answer2 = inputdlg(prompt,dlgtitle,dims,definput,opts);      
        scale = str2num(answer2{1});   % CMZ length correction
        shiftX = str2num(answer2{2});  % X position correction
        shiftY = str2num(answer2{3});  % Y position correction
        maxCMZLen= maxCMZLen_*(1+scale);       
    end        
end
saveas(gca,fullfile([fpath 'img_maxCMZ.jpeg'])); % save fig
close(gcf);
% Max CMZ Length
maxCMZLen;
%%
%::: (6) CutEdge expansion rate :::
% average pixel values on the background 
pixVal_ground(1) = mean2(img_sm(1:DY,:));
pixVal_ground(2) = mean2(img_sm(end-DY:end,:));
pixVal_ref = mean(pixVal_ground);  % critical pixel value for determining the cut edge
% finding the location of the upper boundary
y_mid = round(y_tot/2);  % focus only on the upper-half image
Xrange = [1:DX:x_tot];   % X-exploration range of the fitting points
indY_cut = [];
for x = Xrange
    for y = 1:y_mid  
        if img_sm(y,x) > pixVal_ref
            indY_cut(end+1) = y;  % save Y loc for given x 
            break;
        end
    end  
end
% Least-squared Linear fitting
x_mid = round(length(Xrange)/2);
coefit = polyfit(Xrange(1:x_mid), indY_cut(1:x_mid), 1);      % 1st fitting curve
coefit2 = polyfit(Xrange(x_mid:end), indY_cut(x_mid:end), 1); % 2nd fitting curve
% visualization + manual correction
%----------------------------------
answer = 'No';
coefit_(1) = coefit(1);  % backup the original fitting value
coefit_(2) = coefit(2);
coefit2_(1) = coefit2(1);  
coefit2_(2) = coefit2(2);
numMod_(1) = 0.0;        % initialization
numMod_(2) = 0.0;
numMod_(3) = 0.0;        
numMod_(4) = 0.0;
while strcmp(answer,'No')   % continue to correct, if not satisfied
    close(gcf);   % close image
    figure('visible','on');
    set(gcf, 'Position',  [50, 150, 1200, 500]);
    subplot(1,2,1);
    imshow(img_sm(1:y_mid,:)) % should overlay with the original figure
    hold on;
    plot(Xrange, indY_cut,'ro','MarkerSize',2); % plot the points to be fitted
    plot(Xrange, coefit(1)*Xrange+coefit(2), 'g-', 'Linewidth', 1.5);    % plot 1st new fitting line
    plot(Xrange, coefit2(1)*Xrange+coefit2(2), 'y-', 'Linewidth', 1.5); % plot 2nd new fitting line
    subplot(1,2,2);
    imshow(img_b(1:y_mid,:))
    hold on;
    plot(Xrange, indY_cut,'ro','MarkerSize',2); % plot the points to be fitted
    plot(Xrange, coefit(1)*Xrange+coefit(2), 'g-', 'Linewidth', 1.5);    % plot 1st new fitting line
    plot(Xrange, coefit2(1)*Xrange+coefit2(2), 'y-', 'Linewidth', 1.5); % plot 2nd new fitting line
    hold off;
    
    quest = sprintf('Is cutEdge OK?');
    set(groot,'defaultUicontrolFontSize', 16);  % change font size
    answer = questdlg(quest,'Fitting','Yes','No','Yes');
 
    if strcmp(answer,'No')  % correction
        prompt = {['\fontsize{12} (Left) Slope correction X:' sprintf('\n') ' (1+X) * Slope0  |  ' sprintf('Slope0=%.3f',coefit_(1))], ...
                  ['\fontsize{12} (Left) Intercept correction Y:' sprintf('\n') ' (1+Y) * Intercept0  |  ' sprintf('Intercept0=%.3f',coefit_(2))], ...
                  ['\fontsize{12} (Right) Slope correction X:' sprintf('\n') ' (1+X) * Slope0  |  ' sprintf('Slope0=%.3f',coefit2_(1))], ...
                  ['\fontsize{12} (Right) Intercept correction Y:' sprintf('\n') ' (1+Y) * Intercept0  |  ' sprintf('Intercept0=%.3f',coefit2_(2))]};
        dlgtitle = sprintf('X_{pixels} = %d ; \t Y_{pixels} = %d',x_tot,y_tot);
        dims = [2 80];  % dialog box size
        definput = {num2str(numMod_(1)),num2str(numMod_(2)),num2str(numMod_(3)),num2str(numMod_(4))};  % default values
        opts.Interpreter = 'tex';  % need to switch to "tex" mode, in order to change font size
        answer2 = inputdlg(prompt,dlgtitle,dims,definput,opts);       
        coefit(1)= coefit_(1)*(1+str2num(answer2{1})); % slope correction
        coefit(2)= coefit_(2)*(1+str2num(answer2{2})); % intercept correction
        coefit2(1)= coefit2_(1)*(1+str2num(answer2{3})); % slope correction
        coefit2(2)= coefit2_(2)*(1+str2num(answer2{4})); % intercept correction
        numMod_(1) = str2num(answer2{1});  % save the previous input
        numMod_(2) = str2num(answer2{2});
        numMod_(3) = str2num(answer2{3});  
        numMod_(4) = str2num(answer2{4});
    end
end
saveas(gca,fullfile([fpath 'img_fit_cutEdge.jpeg'])); % show fitting line with images
close(gcf);   % close image
%----------------------------------
% cutEdge expansion speed & distance
cutEdgeSpeed = coefit(1);   % expansion speed
% intersection for detection of expansion distance
x_cross = (coefit2(2)-coefit(2))/(coefit(1)-coefit2(1));
y1 = coefit(1)*1 + coefit(2);
y2 = coefit(1)*x_cross + coefit(2);
cutEdgeDist = abs(y1-y2);   % expansion distance
%%
%::: (7) Intensity profile :::
xLoc = [round(t_maxCMZlen), ...
        round(t_maxCMZlen+0.3*(min(x_tot,t_CMZvanish)-t_maxCMZlen)),...
        round(t_maxCMZlen+0.6*(min(x_tot,t_CMZvanish)-t_maxCMZlen))]; % horizontal position of ROI
tiLabel = {['Profile @ X=' num2str(xLoc(1)) '(t\_maxCMZlen)'], ...
           ['Profile @ X=' num2str(xLoc(2))], ...
           ['Profile @ X=' num2str(xLoc(3))]};              
for i=1:length(xLoc)  
    y1Loc = frontSpeed*xLoc(i)+frontSpeed_intercept;   % vertical position of ROI
    y2Loc = backSpeed*xLoc(i)+backSpeed_intercept;     % vertical position of ROI
    figure('visible','off');
    plot([1:y_tot],img(:,xLoc(i)),'k-','LineWidth',1);
    xlim([1,y_tot]);
    xlabel('Tail     <----->     Head','FontSize',14);
    ylabel('Pixel intensity','FontSize',14);
    hold on;
    yIntenMax = max(img(:,xLoc(i)));
    p1 = plot([y1Loc,y1Loc],[0,yIntenMax],'r--','LineWidth',1);  % wave-front bound
    p2 = plot([y2Loc,y2Loc],[0,yIntenMax],'b--','LineWidth',1);  % wave-back bound
    legend([p1 p2],'Front point','Back point')
    title(tiLabel{i});
    hold off;
    filename = sprintf('intensityProfile_%d.jpeg',i);
    saveas(gca,fullfile([fpath filename])); % save fig
end  
%%
%::: (8) Wave length :::  
switch ind_2wave
    case 1.0 
        xTemp = linspace(t_maxCMZlen, min(x_tot,t_CMZvanish), 5);  % select x_points to be analyzed
        xTemp = xTemp(2:4); % select the middle points
        waveLen_ = []; y_wave1 = []; y_wave2 =[]; 
        for i=1:length(xTemp)
            y_wave1(i) = frontSpeed*(xTemp(i))+frontSpeed_intercept;
            y_wave2(i) = frontSpeed2*(xTemp(i))+frontSpeed2_intercept;
            waveLen_(end+1) = y_wave1(i) - y_wave2(i);
        end
        waveLen = mean(waveLen_);  % Wave length between 1- and 2-wave 
        figure('visible','off');
        imshow(img);
        hold on;
        for i=1:length(xTemp)
            plot([xTemp(i),xTemp(i)],[y_wave1(i), y_wave2(i)],'g-','LineWidth',1);
        end
        saveas(gca,fullfile([fpath 'waveLength.jpeg'])); % save fig
        hold off;
    case 0.0
        waveLen = 0.0;  % no wave length for 1-wave
end
%%
% Comment message for the analysis
prompt = {'\fontsize{12} Enter comments on the analysis:'};
dlgtitle = 'Analysis Note';
dims = [5 70];
definput = {''};
opts.Interpreter = 'tex'; 
comments = inputdlg(prompt,dlgtitle,dims,definput,opts);
% reshape multiple-line comments
numLine = length(comments{1}(:,1));  % # of lines of the input comments
for i=1:numLine
    if ~isspace(comments{1}(i,end))
        comments{1}(i,end+1) = ' ';  % adding blank space to the end of the sentence if the final is NOT a space
    end
end
comments = reshape([comments{:}'],1,[]);  
%%
%--------------------------------------------------------------------------

% Data output
%-------------------------------------------------------------------------- 
%::: data output ::: 
file_ = ['./' fdname_{1} '_Measure.csv'];  % file path
col_header = {'FrontSpeed','FrontSpeed_intercept','backSpeed','backSpeed_intercept','FrontSpeed2','FrontSpeed2_intercept','maxCMZLen','maxCMZDist','waveLen','t_maxCMZlen','t_CMZvanish','cutEdgeSpeed','cutEdgeDist'};
switch ind_2wave
    case 1.0  % 2-wave case
        data = {frontSpeed, frontSpeed_intercept, backSpeed, backSpeed_intercept, frontSpeed2, frontSpeed2_intercept, maxCMZLen, maxCMZDist, waveLen, t_maxCMZlen, t_CMZvanish, cutEdgeSpeed, cutEdgeDist};
    case 0.0  % 1-wave case
        data = {frontSpeed, frontSpeed_intercept, backSpeed, backSpeed_intercept, 0.0, 0.0, maxCMZLen, maxCMZDist, 0.0, t_maxCMZlen, t_CMZvanish, cutEdgeSpeed, cutEdgeDist};
end
% write header to files
if ~isfile(file_)
    fid = fopen(file_,'w');  %write header to file
    fprintf(fid,'%s,',col_header{:});
    fprintf(fid,'%s, %s',[],'Comments'); 
    fprintf(fid,'\n');
    fclose(fid);
end
% append the data to the files     
fid = fopen(file_,'a');  
fprintf(fid,'%f,',data{:}); 
fprintf(fid,'%s, %s',[],comments); 
fprintf(fid,'\n');
fclose(fid);

% convert comments into numerical values
% ::::::::::
T = readtable(file_);
idxCol = [];
numCols = width(T);   % number of cols
numRows = height(T);  % number of rows
keyword = 'Good';
col_header_temp = col_header;

col_header_temp{end+1} = 'emptyCol';
col_header_temp{end+1} = 'Comments';
col_header_temp_len = length(col_header_temp);
for i=(col_header_temp_len+1):numCols
    col_header_temp{end+1} = 'temp';
end

T.Properties.VariableNames = col_header_temp;
for i= 1:numRows
    if contains(T.Comments{i}, keyword)  % check whether comments contain the keyword 
        idxCol(end+1) = 1;
    else
        idxCol(end+1) = 0;
    end
end
idxCol = idxCol';  % transpose into columns

%:: delete NaN colums
TF = ismissing(T,{NaN}); % get matrix for recording NaN element
NA_colIndx = any(TF,1);  % column index of the NaN column
T(:, NA_colIndx) = [];   % delete column
%:: delete columns after 'Comments'
tidx = find(string(T.Properties.VariableNames) == 'Comments');
T(:, tidx+1:end) = [];   % delete column

emptyCol = cell(numRows,1);  % empty content
T = addvars(T,emptyCol,'Before','Comments');   % write the empty column
T = addvars(T,idxCol,'After','Comments');      % write the reference index

writetable(T,file_);  % rewrite the file
% ::::::::::
%--------------------------------------------------------------------------
%%
% Detection of the front and back points of the 1st CMZ wave
%--------------------------------------------------------------------------
x_pos = [];
y_front = [];
y_back = [];
for i= 1:t_maxCMZlen
    x_pos(end+1) = i;  % x-location
    y_front(end+1) = frontSpeed * i + frontSpeed_intercept; % Y-location of the front point 
    y_back(end+1) = Y_t0; % Y-location of the back point
end
for i= (t_maxCMZlen+1):min(x_tot,t_CMZvanish)
    x_pos(end+1) = i;      
    y_front(end+1) = frontSpeed * i + frontSpeed_intercept;  
    y_back(end+1) = backSpeed * i + backSpeed_intercept; 
end 
% visualization for check purpose
figure('visible','off');
imshow(img) 
hold on;
plot(x_pos, y_front,'ro','MarkerSize',3);
plot(x_pos, y_back,'bo','MarkerSize',3);
hold off;
saveas(gca,fullfile([fpath 'img_CMZ_Detection.jpeg']));
% output into a CSV file
file_ = [fpath '/CMZ_detection.csv']; 
col_header = {'X_Pos', 'Y_front', 'Y_back'};
%-----
fid = fopen(file_,'w');  %write header to file
fprintf(fid,'%s,',col_header{:});
fprintf(fid,'\n');
fclose(fid);
%-----
myData = [x_pos', y_front', y_back']; %transform into a matrix
dlmwrite(file_,myData,'-append');
%--------------------------------------------------------------------------
%%
end
