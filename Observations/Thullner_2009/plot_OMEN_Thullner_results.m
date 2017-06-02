function plot_OMEN_Thullner_results()

    % load data
	data.Thullner=xlsread('ThullnerGcubed2009 Fluxes-rates.xlsx','exponential','A2:L8'); 
    % Middelburg all in one
    data.MiddelburgO2flux=xlsread('Domink1_Middelburg96_GBCdataFig_1.xlsx','Combined','C5:D108'); 
    data.MiddelburgNO3flux=xlsread('Domink1_Middelburg96_GBCdataFig_1.xlsx','Combined','G5:H74'); 
    data.MiddelburgNH4flux=xlsread('Domink1_Middelburg96_GBCdataFig_1.xlsx','Combined','K5:L19'); 
    % Middelburg per ocean region
    data.MiddelburgO2fluxAtlantic=xlsread('Domink1_Middelburg96_GBCdataFig_1.xlsx','Combined','C5:D48'); 
    data.MiddelburgO2fluxPacific=xlsread('Domink1_Middelburg96_GBCdataFig_1.xlsx','Combined','C49:D90'); 
    data.MiddelburgO2fluxIndArc=xlsread('Domink1_Middelburg96_GBCdataFig_1.xlsx','Combined','C91:D108'); 
    data.MiddelburgNO3fluxAtlantic=xlsread('Domink1_Middelburg96_GBCdataFig_1.xlsx','Combined','G48:H72'); 
    data.MiddelburgNO3fluxPacific=xlsread('Domink1_Middelburg96_GBCdataFig_1.xlsx','Combined','G5:H47'); 
    data.MiddelburgNO3fluxInd=xlsread('Domink1_Middelburg96_GBCdataFig_1.xlsx','Combined','G73:H74'); 

    
    data.OMEN_gamma95=xlsread('ThullnerGcubed2009 Fluxes-rates.xlsx','exponential','A31:H37'); 
    data.OMEN_gamma05=xlsread('ThullnerGcubed2009 Fluxes-rates.xlsx','exponential','J31:Q37'); 

    set(0,'defaultLineLineWidth', 1)
    set(0,'DefaultAxesFontSize',12)

    figure;
    % SWI-fluxes
    subplot(1,3,1)
    % O2 and TOU
%    scatter(data.MiddelburgO2flux(:,2), -data.MiddelburgO2flux(:,1)/1000,'kd')  % all in one symbol
    hold on
    ylim([-5 0.0])
    scatter(data.MiddelburgO2fluxAtlantic(:,2), -data.MiddelburgO2fluxAtlantic(:,1)/1000,'kd')  % Atlantic
    scatter(data.MiddelburgO2fluxPacific(:,2), -data.MiddelburgO2fluxPacific(:,1)/1000,'ks')  % Pacific
    scatter(data.MiddelburgO2fluxIndArc(:,2), -data.MiddelburgO2fluxIndArc(:,1)/1000,'kx')  % Indian and Arctic
    plot(data.Thullner(:,2), -data.Thullner(:,1)/1000, '-ro', 'MarkerEdgeColor','r')    
    plot(data.Thullner(:,7), -data.Thullner(:,1)/1000, '-ro', 'MarkerEdgeColor','r','MarkerFaceColor','r')    
    plot(-data.OMEN_gamma95(:,2), -data.OMEN_gamma95(:,1)/1000, '-ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    plot(-data.OMEN_gamma05(:,2), -data.OMEN_gamma05(:,1)/1000, '-kv', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    box on
    hold off
    set(gca,'YTick',[-5:0])
    xlabel ({'O_2 flux, TOU'}); %;'(\mumol cm^{-2}yr^{-1})'})
    ylabel('SFD (km)')
%    hleg=legend('Obs Atl', 'Obs PAcific', 'Obs Ind+Arc', 'O_2 flux - Thullner', 'TOU - Thullner','OMEN \gamma=0.95','OMEN \gamma=0.05');
%    set(hleg,'FontSize',4);
%    set(hleg,'Location','SouthEast')

    subplot(1,3,2)
    % NO3
%    scatter(data.MiddelburgNO3flux(:,2), -data.MiddelburgNO3flux(:,1)/1000,'kd')  % all in one symbol
    hold on
    ylim([-5 0.0])
    scatter(data.MiddelburgNO3fluxAtlantic(:,2), -data.MiddelburgNO3fluxAtlantic(:,1)/1000,'kd')  % Atlantic
    scatter(data.MiddelburgNO3fluxPacific(:,2), -data.MiddelburgNO3fluxPacific(:,1)/1000,'ks')  % Pacific
    scatter(data.MiddelburgNO3fluxInd(:,2), -data.MiddelburgNO3fluxInd(:,1)/1000,'kx')  % Indian
    plot(data.Thullner(:,3), -data.Thullner(:,1)/1000, '-ro', 'MarkerEdgeColor','r')    
    plot(-data.OMEN_gamma95(:,3), -data.OMEN_gamma95(:,1)/1000, '-ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    plot(-data.OMEN_gamma05(:,3), -data.OMEN_gamma05(:,1)/1000, '-kv', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    xlim([-20 60.0])
    set(gca,'XTick',[-20:20:60])
    set(gca,'YTick',[-5:0])
   	xlabel ({'NO_3 flux'}); %;'(\mumol cm^{-2}yr^{-1})'})
    box on
    hold off
 %  ylabel('SFD (km)')
    
    subplot(1,3,3)
    % NO3
    plot(data.Thullner(:,4), -data.Thullner(:,1)/1000, '-ro', 'MarkerEdgeColor','r')    
    hold on
    ylim([-5 0.0])
    plot(-data.OMEN_gamma95(:,4), -data.OMEN_gamma95(:,1)/1000, '-ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    plot(-data.OMEN_gamma05(:,4), -data.OMEN_gamma05(:,1)/1000, '-kv', 'MarkerEdgeColor','k','MarkerFaceColor','k')
%    xlim([-50 100.0])
    set(gca,'YTick',[-5:0])
    xlabel ({'SO_4 flux'}); %;'(\mumol cm^{-2}yr^{-1})'})
%    ylabel('SFD (km)')
    box on
    hold off

	print('-depsc2', ['0_OMEN_Thullner_hypsometry_fluxes.eps']);
   	print('-dpsc2', ['0_OMEN_Thullner_hypsometry_fluxes.ps']);
 
   figure;
    % Degradatrion rates
    subplot(1,3,1)
    % O2 and TOU
    plot(-data.OMEN_gamma95(:,5)*10^6, -data.OMEN_gamma95(:,1)/1000, '-ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    hold on
    plot(-data.OMEN_gamma05(:,5)*10^6, -data.OMEN_gamma05(:,1)/1000, '-kv', 'MarkerEdgeColor','k','MarkerFaceColor','k')
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
    plot(-data.OMEN_gamma05(:,6)*10^6, -data.OMEN_gamma05(:,1)/1000, '-kv', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    plot(data.Thullner(:,9), -data.Thullner(:,1)/1000, '-ro', 'MarkerEdgeColor','r')    
    hold off
    xlim([-50 100.0])
   	xlabel ({'Denitrification';'(\mumol cm^{-2}yr^{-1})'})
  %  ylabel('SFD (km)')
    
    subplot(1,3,3)
    % NO3
    plot(-data.OMEN_gamma95(:,7)*10^6, -data.OMEN_gamma95(:,1)/1000, '-ko', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    hold on
    plot(-data.OMEN_gamma05(:,7)*10^6, -data.OMEN_gamma05(:,1)/1000, '-kv', 'MarkerEdgeColor','k','MarkerFaceColor','k')
    plot(data.Thullner(:,10), -data.Thullner(:,1)/1000, '-ro', 'MarkerEdgeColor','r')    
    hold off
%    xlim([-50 100.0])
    xlabel ({'SO_4 reduction';'(\mumol cm^{-2}yr^{-1})'})
%    ylabel('SFD (km)')

	print('-dpsc2', ['0_OMEN_Thullner_hypsometry_rates.ps']);
