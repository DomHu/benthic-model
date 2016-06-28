%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  PLOT the Results of Sensitivity study

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = Plot_sensitivity_single(Results, ind, range, NameVar, xaxis)


% for i=1:length(k1)
%     if(abs(Results(i,7))>0.1)
%         Results(i,7)=0.0;
%     end
% end

    figure
    % F_O2
    subplot(3,3,1)
    plot(ind, Results(:,2),'.')            
    hold on            
    axis ([range(1) range(2) -Inf Inf])
%    xlabel ('log(k_1) [yr^{-1}]')
    ylabel('F_{O_2}')

    % F_NO3
    subplot(3,3,2)
    plot(ind, Results(:,3),'.')            
    hold on            
    axis ([range(1) range(2) -Inf Inf])
%    xlabel ('log(k_1) [yr^{-1}]')
    ylabel('F_{NO_3}')

    % F_SO4
    subplot(3,3,3)
    plot(ind, Results(:,4),'.')            
    hold on            
    axis ([range(1) range(2) -Inf Inf])
%    xlabel ('log(k_1) [yr^{-1}]')
    ylabel('F_{SO_4}')

    % F_PO4
    subplot(3,3,4)
    plot(ind, Results(:,7),'.')            
    hold on            
    axis ([range(1) range(2) -Inf Inf])
%    xlabel ('log(k_1) [yr^{-1}]')
    ylabel('F_{PO_4}')

    % F_NH4
    subplot(3,3,5)
    plot(ind, Results(:,5),'.')            
    hold on            
    axis ([range(1) range(2) -Inf Inf])
%    xlabel ('log(k_1) [yr^{-1}]')
    ylabel('F_{NH_4}')

    % F_H2S
    subplot(3,3,6)
    plot(ind, Results(:,6),'.')            
    hold on            
    axis ([range(1) range(2) -Inf Inf])
%    xlabel ('log(k_1) [yr^{-1}]')
    ylabel('F_{H_2S}')
    
    % zox
    subplot(3,3,7)
    plot(ind, Results(:,8),'.')            
    hold on            
    axis ([range(1) range(2) -Inf Inf])
    xlabel (xaxis)
    ylabel('z_{O_2}')

    % zNO3
    subplot(3,3,8)
    plot(ind, Results(:,9),'.')            
    hold on            
    axis ([range(1) range(2) -Inf Inf])
    xlabel (xaxis)
    ylabel('z_{NO_3}')

    % zox
    subplot(3,3,9)
    plot(ind, Results(:,10),'.')            
    hold on            
    axis ([range(1) range(2) -Inf Inf])
    xlabel (xaxis)
    ylabel('z_{SO_4}')
    
%	ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%	text(0.5, 1,'\bf Test 4\_1: 500m anoxic (no Mn, Fe)','HorizontalAlignment','center','VerticalAlignment', 'top')
	print('-dpsc2', ['1_' NameVar '_Sensitivity.ps']);
end