%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  PLOT the Results of Sensitivity study

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = Plot_sensitivity_singleOutput(Params, Results, index,NameVar, yaxis, str_date)
% % % Set date and time
% % str_date = datestr(now,'ddmmyy_HH_MM_SS')


    figure
    % k1
    subplot(3,3,1)
    plot(log10(Params.k1), Results(:,index),'.')            
    hold on            
    axis ([Params.range_k1(1) Params.range_k1(2) -Inf Inf])
    xlabel ('log(k_1) [yr^{-1}]')
    ylabel(yaxis)

    % f1
    subplot(3,3,2)
    plot(Params.f1, Results(:,index),'.')            
    hold on            
    axis ([0 1 -Inf Inf])
    xlabel ('labile fraction')
%    ylabel('F_{O_2}')

    % KNH4
    subplot(3,3,3)
    plot(Params.KNH4, Results(:,index),'.')            
    hold on            
    axis ([Params.range_KNH4(1) Params.range_KNH4(2) -Inf Inf])
    xlabel ('KNH4 [-]')
%    ylabel('F_{O_2}')

    % KPO4ox
    subplot(3,3,4)
    plot(Params.KPO4ox, Results(:,index),'.')            
    hold on            
    axis ([Params.range_KPO4ox(1) Params.range_KPO4ox(2) -Inf Inf])
    xlabel ('KPO4ox [-]')
    ylabel(yaxis)

    % KPO4anox
    subplot(3,3,5)
    plot(Params.KPO4anox, Results(:,index),'.')            
    hold on            
    axis ([Params.range_KPO4anox(1) Params.range_KPO4anox(2) -Inf Inf])
    xlabel ('KPO4anox [-]')
 %   ylabel(yaxis)

    % ksPO4
    subplot(3,3,6)
    plot(Params.ksPO4, Results(:,index),'.')            
    hold on            
    axis ([Params.range_ksPO4(1) Params.range_ksPO4(2) -Inf Inf])
    xlabel ('ksPO4 [yr^{-1}]')
 %   ylabel(yaxis)
    
    % kmPO4
    subplot(3,3,7)
    plot(Params.kmPO4, Results(:,index),'.')            
    hold on            
    axis ([Params.range_kmPO4(1) Params.range_kmPO4(2) -Inf Inf])
    xlabel ('kmPO4 [yr^{-1}]')
    ylabel(yaxis)

    % kaPO4
    subplot(3,3,8)
    plot(log10(Params.kaPO4), Results(:,index),'.')            
    hold on            
    axis ([Params.range_kaPO4(1) Params.range_kaPO4(2) -Inf Inf])
    xlabel ('log(kaPO4) [yr^{-1}]')
 %   ylabel(yaxis)

    % gammaNH4 gammaH2S
    subplot(3,3,9)
    plot(Params.gammaNH4, Results(:,index),'.')
    plot(Params.gammaH2S, Results(:,index),'r.')            
    hold on            
    axis ([Params.range_gammaNH4(1) Params.range_gammaNH4(2) -Inf Inf])
%    xlabel ('gammaNH4 [-]')
    xlabel ('gammaNH4   gammaH2S [-]')
 %   ylabel(yaxis)
    
%	ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%	text(0.5, 1,'\bf Test 4\_1: 500m anoxic (no Mn, Fe)','HorizontalAlignment','center','VerticalAlignment', 'top')
	print('-depsc2', ['./Sensitivity/2_' NameVar '_Sensitivity_' str_date '.eps']);
end