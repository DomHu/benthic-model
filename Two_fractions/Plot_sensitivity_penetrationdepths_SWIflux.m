%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  PLOT sensitivity of TEA penetration depth and SWI flux O2 and NO3 to e.g. k1 or TOC load 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = Plot_sensitivity_penetrationdepths_SWIflux(Results, ind, range, NameVar, xaxis, str_date)

% % % set date-time
% % str_date = datestr(now,'ddmmyy_HH_MM_SS')

% for i=1:length(k1)
%     if(abs(Results(i,7))>0.1)
%         Results(i,7)=0.0;
%     end
% end
    set(0,'defaultLineLineWidth', 2)
    set(0,'DefaultAxesFontSize',12) % plots 18            

    figure
    % F_O2 and F_NO3 and FSO4
%    subplot(2,1,1)
    hold on          
    box on
    plot(ind, Results(:,2),'bo')   
    plot(ind, Results(:,3),'go')
    plot(ind, Results(:,4),'ro')
    axis ([range(1) range(2) -Inf Inf])
    xlabel (xaxis)
    ylabel(' SWI fluxes [mol cm^{-2} yr^{-1}] ')
    hleg=legend('F_{O_2}','F_{NO_3}','F_{SO_4}');
    set(hleg,'Location','SouthWest')
%	ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%	text(0.5, 1,'\bf Test 4\_1: 500m anoxic (no Mn, Fe)','HorizontalAlignment','center','VerticalAlignment', 'top')
	print('-depsc2', ['./Sensitivity/1_' NameVar '_SWIflux_Sensitivity_' str_date '.eps']);

    figure
    %    subplot(2,1,2)
    % zox zNO3 zSO4
    hold on            
    box on
    plot(ind, Results(:,8),'bo')    
    plot(ind, Results(:,9),'go')        
    plot(ind, Results(:,10),'ro')     
    axis ([range(1) range(2) -Inf Inf])
    xlabel (xaxis)
    ylabel(' Penetration depth [cm] ')
    set(gca,'YDir','reverse');
    hleg=legend('z_{O_2}','z_{NO_3}','z_{SO_4}');
    set(hleg,'Location','SouthEast')
%	ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%	text(0.5, 1,'\bf Test 4\_1: 500m anoxic (no Mn, Fe)','HorizontalAlignment','center','VerticalAlignment', 'top')
	print('-depsc2', ['./Sensitivity/1_' NameVar '_TEA_Sensitivity_' str_date '.eps']);
end

