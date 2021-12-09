function spacing(imgFolder,dataPosFolder,fname_CMZ,fname_kymo)

    %% Data import
    %>> import image sequence
    %:: =======================================================================
    d=dir([imgFolder '*.tif']);
    for i=1:numel(d)  % go over all the imgs in the designated folder
        fprintf(d(i).name + "\n");
        im=imread([imgFolder d(i).name]);
        % Normalization & conversion into double format
        img{i}=mat2gray(im);   
    end
    %:: =======================================================================

    %%
    %>> Import info of nucelus position at each time point 
    %:: =======================================================================
    %(1) import nucelus position at each time frame
    d=dir([dataPosFolder '*.csv']);  % d=dir([dataPosFolder '*.csv']);

    for t=1:numel(d)  % go over all the files in the designated folder
        fprintf(d(t).name + "\n");
        data_raw = readtable([dataPosFolder d(t).name]);
        % Extract nucelus position {X,Y} at each time point
        dataPos{t} = [data_raw.X, data_raw.Y];  % [X, Y]   
    end
    %:: =======================================================================

    %%
    %:: Import info of CMZ at each time points
    %:: =======================================================================         
    data_CMZ = xlsread(fname_CMZ);        
    img_kymo = imread(fname_kymo);    % the smoothed kymo for check purpose

    %!! connect image piece of kymos to slices at a particular time point
    [kymo_row, kymo_col] = size(img_kymo);     % width of kymo
    slice_col = kymo_col / (length(img)-1);    % width of each slice in the kymo

    % calculate the average Y position of CMZ 
    for i=1:length(img)-1   % go over each image (note that a slice is derived from the differential of two imgs)

        x_ROI1 = (i-1)*slice_col+1;   % x-range of the given slice
        x_ROI2 = x_ROI1+slice_col-1;

        data_CMZ_ROI = [];
        data_CMZ_ROI = data_CMZ(data_CMZ(:,1)>x_ROI1 & data_CMZ(:,1)<x_ROI2,:,:); % select data in the region of interest

        Y_front_ave{i} = round(mean(data_CMZ_ROI(:,2)));
        Y_back_ave{i}  = round(mean(data_CMZ_ROI(:,3)));
    end
    %:: =======================================================================

    %%
    %:: Creation a new folder for saving data
    %:: =======================================================================
    fdname_ = regexp(fname_CMZ,'CMZ_detection','split');   % split the data filename
    fpath = [fdname_{1} 'data_analyzed/'];  
    mkdir(fpath)
    %:: =======================================================================
    
    %%
    %:: Crop images & tables according to CMZ positions
    %:: =======================================================================
    % find NaN (caused by CMZ vanishing)
    Y_front_NaN = find(isnan([Y_front_ave{:}]));
    Y_back_NaN  = find(isnan([Y_back_ave{:}]));
    ind_imgMax = min(min(Y_front_NaN), min(Y_back_NaN));

    tMax = max([ind_imgMax, length(img)]); 
    % crop "images" & "cell-nucelus-position tables" with ROI_CMZ
    for t=1:tMax-1
        img_ROI{t} = img{t}(Y_back_ave{t}:Y_front_ave{t},:); 
        img_nonROI{1,t} = img{t}(1:Y_back_ave{t}-1,:);     % upper region of non-CMZ
        img_nonROI{2,t} = img{t}(Y_front_ave{t}+1:end,:);  % lower non-CMZ region
        
        cond1 = (Y_back_ave{t}<dataPos{t}(:,2)) & (dataPos{t}(:,2)<Y_front_ave{t}); 
        dataPos_ROI{t} = dataPos{t}(cond1,:); 
        ind1 = min(find(cond1)); 
        ind2 = max(find(cond1));
        dataPos_nonROI{1,t} = dataPos{t}(1:ind1-1,:);      % upper non-CMZ region
        dataPos_nonROI{2,t} = dataPos{t}(ind2+1:end,:);
        
        %>> save info of CMZ location (used for later analysis)
        t_CMZ(t)   = t;   % time point
        loc_CMZ{t} = [Y_back_ave{t},Y_front_ave{t}];  % upper and lower boundaries of CMZ
    end


    %% Spacing analysis
    %:: Triangulation
    %:: =======================================================================
    for t=1:length(dataPos_ROI)
        x=[]; y=[];
        x = dataPos_ROI{t}(:,1);
        y = dataPos_ROI{t}(:,2);
        DT{t} = delaunay(x,y);  % get connectivity of 'particle-position table'
                                %>> (x,y): point position, DT: row index of points composing the triangle
                                %>> # of rows of DT = # of triangles
    end  
    for i=1:2
        for t=1:length(dataPos_ROI)
            x_non=[]; y_non=[];
            x_non = dataPos_nonROI{i,t}(:,1);
            y_non = dataPos_nonROI{i,t}(:,2);
            DT_non{i,t} = delaunay(x_non,y_non);  
        end    
    end
    %:: =======================================================================

    %%
    %:: Get triangulation properties based on the connectivity list DT
    %:: =======================================================================
    for t=1:length(dataPos_ROI)
        x = dataPos_ROI{t}(:,1);
        y = dataPos_ROI{t}(:,2);
        TR{t} = triangulation(DT{t},[x,y]);
        %>> get list of edges: TR.edges = [vertex1, vertex2]; Row # = edge #
        edges_list{t} = TR{t}.edges;
    end  
    for i=1:2
        for t=1:length(dataPos_ROI)
            x_non = dataPos_nonROI{i,t}(:,1);
            y_non = dataPos_nonROI{i,t}(:,2);
            TR_non{i,t} = triangulation(DT_non{i,t},[x_non,y_non]);
            edges_list_non{i,t} = TR_non{i,t}.edges;
        end
    end
    %:: =======================================================================

    %%
    %:: (Optional) Plotting triangulation with the raw image (CMZ region)
    %:: =======================================================================
    figure('visible','off');
    fname = [fpath 'triangulation_CMZ/'];
    mkdir(fname)
    for t=1:length(dataPos_ROI)
        x = dataPos_ROI{t}(:,1);
        y = dataPos_ROI{t}(:,2);
        imshow(img{t}); hold on,
        triplot(DT{t},x,y,'cyan'),         
        %:: outer boundaries
        edge_outer = freeBoundary(TR{t});  % get vertex # that composes the outer edge
        plot(x(edge_outer),y(edge_outer),'-r','LineWidth',1.5),
        %::
        %:: center points
        plot(dataPos{t}(:,1),dataPos{t}(:,2),'.y', 'MarkerSize',5),
        plot(x,y,'.m', 'MarkerSize',1),        
        %::
        title(['t=' num2str(t)], 'Fontsize',16), hold off
        
        saveas(gcf,[fname num2str(t) '.png']);
        %saveas(gcf,[fname num2str(t) '.tif']);
    end
    %:: =======================================================================
    
    %%
    %:: (Optional) Plotting triangulation with the raw image (non-CMZ region)
    %:: =======================================================================
    figure('visible','off');
    for i=1:2
        fname = [fpath 'triangulation_nonCMZ/' num2str(i) '/'];
        mkdir(fname)      
        for t=1:length(dataPos_ROI)
            x_non = dataPos_nonROI{i,t}(:,1);
            y_non = dataPos_nonROI{i,t}(:,2);
            imshow(img{t}); hold on,
            triplot(DT_non{i,t},x_non,y_non,'cyan'), 
            %:: outer boundaries
            edge_outer = freeBoundary(TR_non{i,t});  % get vertex # that composes the outer edge
            plot(x_non(edge_outer),y_non(edge_outer),'-r','LineWidth',1.5),
            %::
            title(['t=' num2str(t)], 'Fontsize',16), hold off
            
            saveas(gcf,[fname num2str(t) '.png']);
        end
    end
    %:: =======================================================================

    %%    
    %>> calculation of properties according to the edge list (CMZ)
    for t=1:length(dataPos_ROI)
        x = dataPos_ROI{t}(:,1);
        y = dataPos_ROI{t}(:,2);
        for i=1:length(edges_list{t}) % go over all the edges
            vex1 = edges_list{t}(i,1); % vertex 1 that composes the edge
            vex2 = edges_list{t}(i,2); % vertex 2 ...
            vex1_pos = [x(vex1), y(vex1)];  % position of vertex 1
            vex2_pos = [x(vex2), y(vex2)]; 

            edgeLen_list{t}(i) = norm(vex1_pos-vex2_pos);        % list of edge length
            edgeLenX_list{t}(i) = abs(vex1_pos(1)-vex2_pos(1));  % edge length along X
            edgeLenY_list{t}(i) = abs(vex1_pos(2)-vex2_pos(2));  % edge length along Y
        end

        %:: ignore extreme values
        bond_Len  = prctile(edgeLen_list{t},[25,75]); % the 25th & 75th percentiles
        bond_LenX = prctile(edgeLenX_list{t},[25,75]);
        bond_LenY = prctile(edgeLenY_list{t},[25,75]);
        %>> select the specific range 
        cond_Len  = (edgeLen_list{t}>bond_Len(1)) & (edgeLen_list{t}<bond_Len(2));  % points between 25th and 75th
        cond_LenX = (edgeLenX_list{t}>bond_LenX(1)) & (edgeLenX_list{t}<bond_LenX(2));
        cond_LenY = (edgeLenY_list{t}>bond_LenY(1)) & (edgeLenY_list{t}<bond_LenY(2));

        edgeLen_list{t}  = edgeLen_list{t}(cond_Len);  % edge list that satisfies the specific condition
        edgeLenX_list{t} = edgeLenX_list{t}(cond_LenX);
        edgeLenY_list{t} = edgeLenY_list{t}(cond_LenY);
    end
    %:: =======================================================================

    %%    
    %>> calculation of properties according to the edge list (non-CMZ)
    for i=1:2    
        for t=1:length(dataPos_ROI)
            x = dataPos_nonROI{i,t}(:,1);
            y = dataPos_nonROI{i,t}(:,2);
            for j=1:length(edges_list_non{i,t}) 
                vex1 = edges_list_non{i,t}(j,1); 
                vex2 = edges_list_non{i,t}(j,2); 
                vex1_pos = [x(vex1), y(vex1)];  
                vex2_pos = [x(vex2), y(vex2)]; 

                edgeLen_list_non{i,t}(j) = norm(vex1_pos-vex2_pos);        
                edgeLenX_list_non{i,t}(j) = abs(vex1_pos(1)-vex2_pos(1));  
                edgeLenY_list_non{i,t}(j) = abs(vex1_pos(2)-vex2_pos(2));  
            end

            %:: ignore extreme values
            bond_Len  = prctile(edgeLen_list_non{i,t},[25,75]); % the 25th & 75th percentiles
            bond_LenX = prctile(edgeLenX_list_non{i,t},[25,75]);
            bond_LenY = prctile(edgeLenY_list_non{i,t},[25,75]);
            %>> select the specific range 
            cond_Len  = (edgeLen_list_non{i,t}>bond_Len(1)) & (edgeLen_list_non{i,t}<bond_Len(2));  % points between 25th and 75th
            cond_LenX = (edgeLenX_list_non{i,t}>bond_LenX(1)) & (edgeLenX_list_non{i,t}<bond_LenX(2));
            cond_LenY = (edgeLenY_list_non{i,t}>bond_LenY(1)) & (edgeLenY_list_non{i,t}<bond_LenY(2));

            edgeLen_list_non{i,t}  = edgeLen_list_non{i,t}(cond_Len);  % edge list that satisfies the specific condition
            edgeLenX_list_non{i,t} = edgeLenX_list_non{i,t}(cond_LenX);
            edgeLenY_list_non{i,t} = edgeLenY_list_non{i,t}(cond_LenY);
        end
    end
    %:: =======================================================================    
        
    %%
    %>> statisitical analysis according to the lists (CMZ)
    spacing_ave  =[];
    Xspacing_ave =[];
    Yspacing_ave =[];
    for t=1:length(dataPos_ROI)
        %:: Average nucleus-nucleus spacing
        spacing_ave{t} = mean(edgeLen_list{t});
        %:: N-N spacing along the Y-direction
        Xspacing_ave{t} = mean(edgeLenX_list{t});
        %:: N-N spacing along the Y-direction
        Yspacing_ave{t} = mean(edgeLenY_list{t});   
    end
    %:: =======================================================================
      
    %%
    %>> statisitical analysis according to the lists (non-CMZ)
    spacing_ave_non  =[];
    Xspacing_ave_non =[];
    Yspacing_ave_non =[];
    for t=1:length(dataPos_ROI)
        %:: Average nucleus-nucleus spacing
        spacing_ave_non{t} = mean([edgeLen_list_non{1,t}, edgeLen_list_non{2,t}]);
        %:: N-N spacing along the Y-direction
        Xspacing_ave_non{t} = mean([edgeLenX_list_non{1,t}, edgeLenX_list_non{2,t}]);
        %:: N-N spacing along the Y-direction
        Yspacing_ave_non{t} = mean([edgeLenY_list_non{1,t}, edgeLenY_list_non{2,t}]);   
    end
    %:: =======================================================================

    %%
    %:: (Optional) Plotting the statistics of triangulation
    %:: =======================================================================
    figure('visible','off');
    fname = [fpath 'spacing/'];
    mkdir(fname)
    
    t = [1:length(dataPos_ROI)];
    y  = [spacing_ave{:}];  y_non = [spacing_ave_non{:}];
    yX = [Xspacing_ave{:}]; yX_non = [Xspacing_ave_non{:}]; % the X-component
    yY = [Yspacing_ave{:}]; yY_non = [Yspacing_ave_non{:}]; % the Y-component   
    
    subplot(3,1,1); plot(t,y,'-ro'); hold on; plot(t,y_non,'-m^');
    ylabel('Ave Spacing'); legend({'CMZ','non-CMZ'},'Location','northwest');
    subplot(3,1,2); plot(t,yX,'-bo'); hold on; plot(t,yX_non,'-cyan^');
    ylabel('Ave Spacing X'); legend({'CMZ','non-CMZ'},'Location','northwest');
    subplot(3,1,3); plot(t,yY,'-go'); hold on; plot(t,yY_non,'-k^');
    xlabel('Time'); ylabel('Ave Spacing Y'); legend({'CMZ','non-CMZ'},'Location','northwest');
    
    saveas(gcf,[fname 'spacing_analyzed.png']);
    %:: =======================================================================
    
    %%
    %:: Output data into files
    %:: =======================================================================
    fname = [fpath 'spacing/spacing_analyzed.xls'];
    
    col_header  = {'Time','Spacing_ave','XSpacing_ave','YSpacing_ave'};
    xlswrite(fname,col_header,'CMZ','A1');                 %Write column header
    xlswrite(fname,[t',y',yX',yY'],'CMZ','A2');            %Write data   
    xlswrite(fname,col_header,'nonCMZ','A1');                
    xlswrite(fname,[t',y_non',yX_non',yY_non'],'nonCMZ','A2');           
    %:: =======================================================================
    
    %%
    %:: Spacing profile along the PD axis
    %:: =======================================================================
    % divid PD axis
    PDaxis_ = linspace(1,kymo_row,21);
    PDaxis = 0.5*(PDaxis_(1:end-1)+PDaxis_(2:end));
    
    for indCol=1:length(PDaxis)   % go along the PD axis
        col_ind1 = round(PDaxis_(indCol+0));    % top Y boundary
        col_ind2 = round(PDaxis_(indCol+1));    % lower Y boundary
        
        %:: find the nucleus position within the desired region, {col_ind1, col_ind2}
        for t=1:length(dataPos)
            ind_sub = find(dataPos{t}(:,2)>col_ind1 & dataPos{t}(:,2)<col_ind2);
            dataPos_sub{indCol,t} = dataPos{t}(ind_sub,:);  % nucleus position
        end
        
        %:: Triangulation & properties & analysis
        for t=1:length(dataPos)
            [nRow_,nCol_] = size(dataPos_sub{indCol,t});
            %===========================================
            if nRow_>3   % more than 3 cell nucleus composes of a triangle
                %>> Triangulation
                x=[]; y=[];
                x = dataPos_sub{indCol,t}(:,1);
                y = dataPos_sub{indCol,t}(:,2);
                DT_all{indCol,t} = delaunay(x,y); 
                
                %>> Edge properties
                TR{t} = triangulation(DT_all{indCol,t},[x,y]);
                edges_list_all{indCol,t} = TR{t}.edges;
                
                %>> Edge length analysis
                for i=1:length(edges_list_all{indCol,t}) 
                    vex1 = edges_list_all{indCol,t}(i,1); 
                    vex2 = edges_list_all{indCol,t}(i,2); 
                    vex1_pos = [x(vex1), y(vex1)];  
                    vex2_pos = [x(vex2), y(vex2)]; 

                    edgeLen_list_all{indCol,t}(i) = norm(vex1_pos-vex2_pos);        % list of edge length
                    edgeLenX_list_all{indCol,t}(i) = abs(vex1_pos(1)-vex2_pos(1));  % edge length along X
                    edgeLenY_list_all{indCol,t}(i) = abs(vex1_pos(2)-vex2_pos(2));  % edge length along Y
                end
                %:: ignore extreme values
                bond_Len  = prctile(edgeLen_list_all{indCol,t},[25,75]); % the 25th & 75th percentiles
                bond_LenX = prctile(edgeLenX_list_all{indCol,t},[25,75]);
                bond_LenY = prctile(edgeLenY_list_all{indCol,t},[25,75]);
                %>> select the specific range 
                cond_Len  = (edgeLen_list_all{indCol,t}>bond_Len(1)) & (edgeLen_list_all{indCol,t}<bond_Len(2));  % points between 25th and 75th
                cond_LenX = (edgeLenX_list_all{indCol,t}>bond_LenX(1)) & (edgeLenX_list_all{indCol,t}<bond_LenX(2));
                cond_LenY = (edgeLenY_list_all{indCol,t}>bond_LenY(1)) & (edgeLenY_list_all{indCol,t}<bond_LenY(2));
                edgeLen_list_all{indCol,t}  = edgeLen_list_all{indCol,t}(cond_Len);  % edge list that satisfies the specific condition
                edgeLenX_list_all{indCol,t} = edgeLenX_list_all{indCol,t}(cond_LenX);
                edgeLenY_list_all{indCol,t} = edgeLenY_list_all{indCol,t}(cond_LenY); 
                               
                %>> Average nucleus-nucleus spacing
                spacing_ave_all{indCol,t} = mean(edgeLen_list_all{indCol,t});     % spacing_ave_all{PDaxis, Time}
                Xspacing_ave_all{indCol,t} = mean(edgeLenX_list_all{indCol,t});
                Yspacing_ave_all{indCol,t} = mean(edgeLenY_list_all{indCol,t});   
            end
            %===========================================
        end
    end
    %:: =======================================================================
      
    %%
    %:: Plot spacing profile along PD axis at each time point
    %:: =======================================================================  
    fname = [fpath 'spacing_profile/'];
    mkdir(fname) 
    
    %::!! artifically add CMZ location to the final time point
    loc_CMZ{length(dataPos)} = [0,0];
    %::
    
    for t=1:length(dataPos)
        x =[];
        y =[];
        yX=[];
        yY=[];
        for i=1:length(PDaxis)   % select non-empty data
            if ~isempty([spacing_ave_all{i,t}, Xspacing_ave_all{i,t}, Yspacing_ave_all{i,t}])
                x(end+1)  = PDaxis(i);
                y(end+1)  = spacing_ave_all{i,t}; 
                yX(end+1) = Xspacing_ave_all{i,t};
                yY(end+1) = Yspacing_ave_all{i,t};
            end
        end

        fig = figure('visible','off');
        
        subplot(3,1,1); plot(x,y,'-ro'); hold on; 
        YL = get(gca, 'YLim'); plot([loc_CMZ{t}(1) loc_CMZ{t}(1)],[YL(1),YL(2)],'--k'); plot([loc_CMZ{t}(2) loc_CMZ{t}(2)],[YL(1),YL(2)],'--k'); 
        ylabel('Ave Spacing'); hold off;
        subplot(3,1,2); plot(x,yX,'-bo'); hold on; 
        YL = get(gca, 'YLim'); plot([loc_CMZ{t}(1) loc_CMZ{t}(1)],[YL(1),YL(2)],'--k'); plot([loc_CMZ{t}(2) loc_CMZ{t}(2)],[YL(1),YL(2)],'--k'); 
        ylabel('Ave Spacing X'); hold off;
        subplot(3,1,3); plot(x,yY,'-go'); hold on;
        YL = get(gca, 'YLim'); plot([loc_CMZ{t}(1) loc_CMZ{t}(1)],[YL(1),YL(2)],'--k'); plot([loc_CMZ{t}(2) loc_CMZ{t}(2)],[YL(1),YL(2)],'--k'); 
        xlabel('PD axis (Tail -- Body)'); ylabel('Ave Spacing Y'); hold off;
%        suptitle(['t=' num2str(t)]);
    
        saveas(gcf,[fname num2str(t) '.png']); close(fig);
        
        
        %:: output data into Excel
        fname2 = [fname num2str(t) '.xls'];
        col_header  = {'PDaxis','Spacing_ave','XSpacing_ave','YSpacing_ave'};
        xlswrite(fname2,col_header,'sheet1','A1');                 %Write column header
        xlswrite(fname2,[x',y',yX',yY'],'sheet1','A2');            %Write data                            
    end
    %:: =======================================================================
end

 