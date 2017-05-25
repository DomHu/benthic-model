function plot_OMEN_Thullner_results()

    % load data
	data.Thullner=xlsread('ThullnerGcubed2009 Fluxes-rates.xlsx','exponential','A2:L8'); 
    data.OMEN_gamma95=xlsread('ThullnerGcubed2009 Fluxes-rates.xlsx','exponential','A31:H37'); 
    data.OMEN_gamma05=xlsread('ThullnerGcubed2009 Fluxes-rates.xlsx','exponential','J31:Q37'); 

    set(0,'defaultLineLineWidth', 2)
    set(0,'DefaultAxesFontSize',12)

    figure;
    % SWI-fluxes
    subplot(1,3,1)
    % O2 and TOU
    plot(-data.OMEN_gamma95(:,2), -data.OMEN_gamma95(:,1)/1000, '-ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    hold on
    plot(-data.OMEN_gamma05(:,2), -data.OMEN_gamma05(:,1)/1000, '--ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    plot(data.Thullner(:,2), -data.Thullner(:,1)/1000, '-ro', 'MarkerEdgeColor','r')    
    plot(data.Thullner(:,7), -data.Thullner(:,1)/1000, '-ro', 'MarkerEdgeColor','r','MarkerFaceColor','r')    
    hold off
    %                ylim([-50 0.0])
    xlabel ({'O_2 flux, TOU';'(\mumol cm^{-2}yr^{-1})'})
    ylabel('SFD (km)')
    hleg=legend('OMEN \gamma=0.95','OMEN \gamma=0.05', 'O_2 flux - Thullner', 'TOU - Thullner');
    set(hleg,'FontSize',4);
    set(hleg,'Location','SouthEast')

    subplot(1,3,2)
    % NO3
    plot(-data.OMEN_gamma95(:,3), -data.OMEN_gamma95(:,1)/1000, '-ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    hold on
    plot(-data.OMEN_gamma05(:,3), -data.OMEN_gamma05(:,1)/1000, '--ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    plot(data.Thullner(:,3), -data.Thullner(:,1)/1000, '-ro', 'MarkerEdgeColor','r')    
    hold off
    xlim([-50 100.0])
   	xlabel ({'NO_3 flux';'(\mumol cm^{-2}yr^{-1})'})
  %  ylabel('SFD (km)')
    
    subplot(1,3,3)
    % NO3
    plot(-data.OMEN_gamma95(:,4), -data.OMEN_gamma95(:,1)/1000, '-ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    hold on
    plot(-data.OMEN_gamma05(:,4), -data.OMEN_gamma05(:,1)/1000, '--ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    plot(data.Thullner(:,4), -data.Thullner(:,1)/1000, '-ro', 'MarkerEdgeColor','r')    
    hold off
%    xlim([-50 100.0])
    xlabel ({'SO_4 flux';'(\mumol cm^{-2}yr^{-1})'})
%    ylabel('SFD (km)')

	print('-depsc2', ['0_OMEN_Thullner_hypsometry_fluxes.eps']);
    
    
   figure;
    % Degradatrion rates
    subplot(1,3,1)
    % O2 and TOU
    plot(-data.OMEN_gamma95(:,5)*10^6, -data.OMEN_gamma95(:,1)/1000, '-ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    hold on
    plot(-data.OMEN_gamma05(:,5)*10^6, -data.OMEN_gamma05(:,1)/1000, '--ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    plot(data.Thullner(:,8), -data.Thullner(:,1)/1000, '-ro', 'MarkerEdgeColor','r')    
    hold off
    %                ylim([-50 0.0])
    xlabel ({'Aerobic';'(\mumol cm^{-2}yr^{-1})'})
    ylabel('SFD (km)')
    hleg=legend('OMEN \gamma=0.95','OMEN \gamma=0.05', 'Thullner results');
    set(hleg,'FontSize',4);
    set(hleg,'Location','SouthEast')

    subplot(1,3,2)
    % NO3
    plot(-data.OMEN_gamma95(:,6)*10^6, -data.OMEN_gamma95(:,1)/1000, '-ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    hold on
    plot(-data.OMEN_gamma05(:,6)*10^6, -data.OMEN_gamma05(:,1)/1000, '--ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    plot(data.Thullner(:,9), -data.Thullner(:,1)/1000, '-ro', 'MarkerEdgeColor','r')    
    hold off
    xlim([-50 100.0])
   	xlabel ({'Denitrification';'(\mumol cm^{-2}yr^{-1})'})
  %  ylabel('SFD (km)')
    
    subplot(1,3,3)
    % NO3
    plot(-data.OMEN_gamma95(:,7)*10^6, -data.OMEN_gamma95(:,1)/1000, '-ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    hold on
    plot(-data.OMEN_gamma05(:,7)*10^6, -data.OMEN_gamma05(:,1)/1000, '--ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    plot(data.Thullner(:,10), -data.Thullner(:,1)/1000, '-ro', 'MarkerEdgeColor','r')    
    hold off
%    xlim([-50 100.0])
    xlabel ({'SO_4 reduction';'(\mumol cm^{-2}yr^{-1})'})
%    ylabel('SFD (km)')

	print('-depsc2', ['0_OMEN_Thullner_hypsometry_rates.eps']);
