% Modeling epithelial migration during re-epithelialization process
% Developed by Dr.Fu-Lai Wen (email: fulai@gs.ncku.edu.tw)
% tested on the versions R2018b & R2021a
%%
%:: Model Parameters ::
%==========================================================================
N = 85;        % # of cell points
x0 = 1.00;     % preferred cell-cell spacing
L0 = 2.00;     % cell-cell spacing @ t=0
Tmax = 66;     % max simulation time step 
dT =  0.5;     % time interval for Euler update

%>>>> frictional gradient >>>>>>> 
etaP         = 30.0;    % viscosity @ the proximal end   
etaD         = 12.0;    % viscosity @ the distal end
eta(1:Tmax,1:N) = ones(Tmax,1) .*linspace(etaP,etaD,N); 

%>>>> motility forces >>>>>>>>>>>>>
F(1:Tmax,1:N)   = 0.00;   
K0              = 3.5;                  % elastic modulus of cell-cell interaction
K(1:Tmax,1:N-1) = K0;
f0              = 7.0*K0*abs(L0-1.0);   % initial motility force 
Td = 50;    
F(1:Td,N)       = -f0/Td.*[0:Td-1]+f0;  % motility force profile of the most top cell
F(Td+1:Tmax,N)  = 0.0;                  % motility force profile of the most top cell
Fc              = 2.5*K0*abs(L0-1.0);   % threshold for activation of motility force 
indF(2:N-1)     = 0.0;                 

%>>>> confinement forces >>>>>>>>
FBC1(1:Tmax,1:N) = 0.00;   
xC1  = (N+2)*L0;                        % location of the wound edge
fC1  = 0.02*f0;                         % resistive force
xC1o = 1.0*L0;                          % characteristic length 
%==========================================================================
%%
%:: Initial conditions
%==========================================================================
x(1,1) = 0.0;   % cell point position, x(Time, cell index)
v(1,1) = 0.0;   % velocity of the most bottom cell 

for i=2:N
    x(1,i) = x(1,i-1)+L0;
    v(1,i) = 0.0;   
    FBC1(1,i) = -fC1*exp(-(xC1-x(1,i))/xC1o);  % resistive force from wounds
end

v(1,N) = ( FBC1(1,N)+F(1,N)-K(1,N-1)*(x(1,N)-x(1,N-1)-x0) )/eta(1,N);   % velocity of the most top cell
%==========================================================================
%%
%:: Evolution ::
%==========================================================================
for t=2:Tmax
        %:: the most bottom cell
        x(t,1) = x(t-1,1);  % fixed boundary condition  
        v(t,1) = 0.0;       % fixed boundary condition
        
    for i=2:N-1  
        F1 = +K(t-1,i+0) * (x(t-1,i+1)-x(t-1,i+0)-x0);    % cell-cell interaction
        F2 = -K(t-1,i-1) * (x(t-1,i+0)-x(t-1,i-1)-x0);    % cell-cell interaction
        Fnet = F1+F2;

        %:: stretch-induced motility force  
        if (F1-F2 > Fc) && (indF(i) < 0.5) 
             F(t-1,i) = f0;        
             indF(i)  = 1.0;
             dF = 0.04*dT*f0;          
        elseif (indF(i) > 0.5) && (F(t-2,i) > dF)
             F(t-1,i) = F(t-2,i) - dF;    % degradation of motility force
        else
             F(t-1,i) = 0.0;
             indF(i) = 0.0;               % reset for re-excitability
        end
                   
        %:: wound-induced resistive force
        FBC1(t-1,i) = -fC1*exp(-(xC1-x(t-1,i))/xC1o);
        
        %:: position-dependent drag coefficient
        eta(t-1,i) = x(t-1,i)/x(1,end)*(etaD-etaP)+ etaP; 
        
        %:: update cell-point position using Euler's forward method         
        x(t,i)   = x(t-1,i)+ dT*(FBC1(t-1,i)+F(t-1,i)+Fnet)/eta(t-1,i);
        v(t-1,i) = (x(t,i)-x(t-1,i))/dT;   
    end
    
        %:: the most top cell
        FBC1(t-1,N) = -fC1*exp(-(xC1-x(t-1,N))/xC1o);
        eta(t-1,N) = x(t-1,N)/x(1,end)*(etaD-etaP)+ etaP;
        x(t,N)   = x(t-1,N) + dT*( FBC1(t-1,N)+F(t-1,N)-K(t-1,N-1)*(x(t-1,N)-x(t-1,N-1)-x0) )/eta(t-1,N); 
        v(t-1,N) = (x(t,N)-x(t-1,N))/dT; 
end
%==========================================================================
%%
%:: (Optional) Data Visualization: In silico Kymograph
%==========================================================================
fileID = sprintf('./KymoCMZ');
pic = figure('visible','off');
fs = 16;    
threV0 = 0.874;   % critical velocity for determination of cells in CMZ

for t=1:Tmax-1
    indxCMZ{t} = find(abs(v(t,:)) > threV0);        
    
    for ii=1:length(indxCMZ{t})    
        plot([t], [x(t,indxCMZ{t}(ii))],'r.','LineWidth',2,'MarkerFaceColor','red','MarkerSize',10);  
        hold on;
    end
    ylim([0 200]);
    %:: reverse Yticks to make consistence with experiments
    yticks([0, 50, 100, 150, 200]);
    yticklabels({'200','150','100','50','0'});
    xticks([0, 14, 28, 42, 56, 70]);
end
xlabel('Simulation Time Step','interpreter','Latex','FontSize',28);
ylabel('Proximal $\leftrightarrow$ Distal','interpreter','Latex','FontSize',28);
title('{\it In silico} Kymograph','interpreter','Latex','FontSize',28);
set(gca,'FontSize',fs);
hold off;
print(pic, fileID,'-dpng','-r600');   % save into a file
%==========================================================================