%:: Plotting Kymograph 

%%
%:: Data import
load('./workspace_wave.mat'); 

%%
%:: Plotting
pic = figure('visible','on');  % make it "off" if you dont wanna show figure
picInfo = pic.Position;  % get default setting of the plot; picInfo = [xLocation, yLocation, width, height]
ASratio = [1, 2];   % the aspect ratio, [width, height], of the plot; width = height = 1 for a squared figure 
set(gcf,'position',[picInfo(1), picInfo(2), picInfo(3)*ASratio(1), picInfo(4)*ASratio(2)]); % picInfo(3) is for width & picInfo(4) is for height of the plot
for t=1:tF:Tmax-1
    indxCMZ{t} = find(abs(v(t,:)) > threV0);        % index of cell pairing within CMZ
    
    for ii=1:length(indxCMZ{t})    % run over all cell pairs of CMZ       
        plot([t], [x(t,indxCMZ{t}(ii))],'ro','MarkerFaceColor','red','MarkerSize',3);  %!! should be in the real space 
        hold on;
    end
    ylim([0 200]);
    %:: (new) reverse yticks
    yticks([0, 50, 100, 150, 200]);
    yticklabels({'200','150','100','50','0'});
    xticks([0, 14, 28, 42, 56, 70]);
end
set(gca,'FontSize',16);  % fontsize for the ticks
set(gca,'LineWidth',5);  % frame tickness
xlabel('Time','interpreter','Latex','FontSize',24);
ylabel('Proximal $\leftrightarrow$ Distal','interpreter','Latex','FontSize',24);
hold off;

%>> for saving figure into a file (Comment it out if you want to save into file)
% fileID = './CMZ_Kymo';
% print(pic, fileID,'-dpng','-r600');  % save into a png file with resolution 600 dpi (you can increase resolution if you want)