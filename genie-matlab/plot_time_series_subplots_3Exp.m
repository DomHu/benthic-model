%function [] = plot_time_series(PEXP1,PEXP2,PVAR1,PVAR2,PT1,PT2,PIK,PMASK,PCSCALE,PCMIN,PCMAX,PCN,PDATA,POPT,PNAME)
% plot time-series
clear all;

% plot mean (true) or total (false)
plot_mean = false;

% set experiment 
exp_1 = './cgenie_output/EXAMPLE.worjh2.Archeretal2009.SPIN1';
exp_2 = './cgenie_output/1507_43_Archeretal2009_OMEN.boudreau1997_100';
exp_3 = './cgenie_output/1707_41_Archeretal2009_OMEN.boudreau1997_25_233years';
% %%%% load other data

REF_sed_O2_exp_1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_O2.res'),'ascii');
REF_sed_O2_exp_2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_O2.res'),'ascii');
REF_sed_O2_exp_3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_O2.res'),'ascii');

REF_sed_SO4_exp_1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_SO4.res'),'ascii');
REF_sed_SO4_exp_2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_SO4.res'),'ascii');
REF_sed_SO4_exp_3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_SO4.res'),'ascii');

REF_sed_H2S_exp_1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_H2S.res'),'ascii');
REF_sed_H2S_exp_2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_H2S.res'),'ascii');
REF_sed_H2S_exp_3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_H2S.res'),'ascii');

REF_sed_PO4_exp_1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_PO4.res'),'ascii');
REF_sed_PO4_exp_2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_PO4.res'),'ascii');
REF_sed_PO4_exp_3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_PO4.res'),'ascii');

REF_sed_ALK_exp_1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_ALK.res'),'ascii');
REF_sed_ALK_exp_2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_ALK.res'),'ascii');
REF_sed_ALK_exp_3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_ALK.res'),'ascii');

REF_sed_DIC_exp_1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_DIC.res'),'ascii');
REF_sed_DIC_exp_2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_DIC.res'),'ascii');
REF_sed_DIC_exp_3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_DIC.res'),'ascii');

set(0,'defaultLineLineWidth', 1)
set(0,'DefaultAxesFontSize',8)

figure
grid on
hold on

subplot(3,2,1)
if(plot_mean) % mean (mol/kg)
    plot(REF_sed_O2_exp_1(:,1),REF_sed_O2_exp_1(:,3)*1e+6,'b',REF_sed_O2_exp_2(:,1),REF_sed_O2_exp_2(:,3)*1e+6,'r--',REF_sed_O2_exp_3(:,1),REF_sed_O2_exp_3(:,3)*1e+6,'g:'); 
    ylabel('O_2 (\mumol kg^{-1})');
else
    % total (mol)
    plot(REF_sed_O2_exp_1(:,1),REF_sed_O2_exp_1(:,2),'b',REF_sed_O2_exp_2(:,1),REF_sed_O2_exp_2(:,2),'r--',REF_sed_O2_exp_3(:,1),REF_sed_O2_exp_3(:,2),'g:'); 
    ylabel('O_2 (mol)');
end
xlabel('yrs ');
%hleg=legend('Abiotic - No OMEN', 'Abiotic - with OMEN', 'Biotic - all remin. k=0.1', 'Biotic - smaller k'); 
%ylim([0 300])
%set(hleg,'FontSize',4);
%set(hleg,'Location','best')


subplot(3,2,2)
if(plot_mean)  % mean (mol/kg)
	plot(REF_sed_SO4_exp_1(:,1),REF_sed_SO4_exp_1(:,3)*1e+6,'b',REF_sed_SO4_exp_2(:,1),REF_sed_SO4_exp_2(:,3)*1e+6,'r--',REF_sed_SO4_exp_3(:,1),REF_sed_SO4_exp_3(:,3)*1e+6,'g:'); 
    ylabel('SO_4 (\mumol kg^{-1})');
else    % total (mol)
    plot(REF_sed_SO4_exp_1(:,1),REF_sed_SO4_exp_1(:,2),'b',REF_sed_SO4_exp_2(:,1),REF_sed_SO4_exp_2(:,2),'r--',REF_sed_SO4_exp_3(:,1),REF_sed_SO4_exp_3(:,2),'g:'); 
    ylabel('SO_4 (mol)');
end
xlabel('yrs ');
hleg=legend('Abiotic - with OMEN', 'Biotic - k=0.1, no 2nd redox', 'Biotic - k=1.0, no 2nd redox'); 
set(hleg,'FontSize',4);
set(hleg,'Location','SouthEast');

subplot(3,2,3)
if(plot_mean) % mean (mol/kg)
	plot(REF_sed_H2S_exp_1(:,1),REF_sed_H2S_exp_1(:,3)*1e+6,'b',REF_sed_H2S_exp_2(:,1),REF_sed_H2S_exp_2(:,3)*1e+6,'r--', REF_sed_H2S_exp_3(:,1),REF_sed_H2S_exp_3(:,3)*1e+6,'g:');
    ylabel('H_2S (\mumol kg^{-1})');
else % total (mol)
    plot(REF_sed_H2S_exp_1(:,1),REF_sed_H2S_exp_1(:,2),'b',REF_sed_H2S_exp_2(:,1),REF_sed_H2S_exp_2(:,2),'r--', REF_sed_H2S_exp_3(:,1),REF_sed_H2S_exp_3(:,2),'g:', REF_sed_H2S_exp_3(:,1),REF_sed_H2S_exp_3(:,2),'g:');
    ylabel('H_2S (mol)');
end
xlabel('yrs ');
%hleg=legend('small k', 'higher k'); 
%set(hleg,'Location','SouthEast')

subplot(3,2,4)
if(plot_mean) % mean (mol/kg)
    plot(REF_sed_PO4_exp_1(:,1),REF_sed_PO4_exp_1(:,3)*1e+6,'b',REF_sed_PO4_exp_2(:,1),REF_sed_PO4_exp_2(:,3)*1e+6,'r--', REF_sed_PO4_exp_3(:,1),REF_sed_PO4_exp_3(:,3)*1e+6,'g:'); 
    ylabel('PO_4 (\mumol kg^{-1})');
else % total (mol)
    plot(REF_sed_PO4_exp_1(:,1),REF_sed_PO4_exp_1(:,2),'b',REF_sed_PO4_exp_2(:,1),REF_sed_PO4_exp_2(:,2),'r--', REF_sed_PO4_exp_3(:,1),REF_sed_PO4_exp_3(:,2),'g:');
    ylabel('PO_4 (mol)');
end
xlabel('yrs ');
%ylim([2.15 2.17])
% hleg=legend('small k', 'higher k', 'SPIN no Corg'); 
% set(hleg,'Location','NorthEast')

subplot(3,2,5)
if(plot_mean)  % mean (mol/kg)
    plot(REF_sed_ALK_exp_1(:,1),REF_sed_ALK_exp_1(:,3)*1e+6,'b',REF_sed_ALK_exp_2(:,1),REF_sed_ALK_exp_2(:,3)*1e+6,'r--', REF_sed_ALK_exp_3(:,1),REF_sed_ALK_exp_3(:,3)*1e+6,'g:'); 
    ylabel('ALK (\mumol kg^{-1})');
else
	plot(REF_sed_ALK_exp_1(:,1),REF_sed_ALK_exp_1(:,2),'b',REF_sed_ALK_exp_2(:,1),REF_sed_ALK_exp_2(:,2),'r--', REF_sed_ALK_exp_3(:,1),REF_sed_ALK_exp_3(:,2),'g:');
    ylabel('ALK (mol)');
end
xlabel('yrs ');
%hleg=legend('small k', 'higher k'); 
%set(hleg,'FontSize',10)
%set(hleg,'Location','SouthEast')
%title('oceanic Alkalinity','FontSize',18);

subplot(3,2,6)
if(plot_mean)  % mean (mol/kg)
    plot(REF_sed_DIC_exp_1(:,1),REF_sed_DIC_exp_1(:,3)*1e+6,'b',REF_sed_DIC_exp_2(:,1),REF_sed_DIC_exp_2(:,3)*1e+6,'r--', REF_sed_DIC_exp_3(:,1),REF_sed_DIC_exp_3(:,3)*1e+6,'g:'); 
    ylabel('DIC (\mumol kg^{-1})');
else
	plot(REF_sed_DIC_exp_1(:,1),REF_sed_DIC_exp_1(:,2),'b',REF_sed_DIC_exp_2(:,1),REF_sed_DIC_exp_2(:,2),'r--', REF_sed_DIC_exp_3(:,1),REF_sed_DIC_exp_3(:,2),'g:');
    ylabel('DIC (mol)');
end
xlabel('yrs ');
% hleg=legend('small k', 'higher k', 'SPIN no Corg'); 
% set(hleg,'Location','SouthEast')

% set(gcf,'NextPlot','add');
% axes;
% h = title('Global ocean values (micromol kg**-1)');
% set(gca,'Visible','off');
% set(h,'Visible','on'); 

if(plot_mean)  % mean (mol/kg)
    print('-depsc', 'cgenie_output/0_PLOTS/plots_1002/3_MEAN-time-series_NO_2redox_1002_5Exp');
else
    print('-depsc', 'cgenie_output/1707_CHECK_CRASH');
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% atm O2 and CO2
% %%%% load other data
% REF experiments
if(false)

REF_sed_pCO2_allremin = load(fullfile(exp_1,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
REF_sed_pCO2_smallerk = load(fullfile(exp_2,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
REF_sed_pCO2_NoOMEN = load(fullfile(exp_3,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
REF_sed_pCO2_SPIN = load(fullfile(exp_4,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
REF_sed_pCO2_Exp_5 = load(fullfile(exp_5,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
% REF_sed_pCO2_smallerk_DISS_HACK = load(fullfile(exp_6,'/biogem/biogem_series_atm_pCO2.res'),'ascii');

REF_sed_pO2_allremin = load(fullfile(exp_1,'/biogem/biogem_series_atm_pO2.res'),'ascii');
REF_sed_pO2_smallerk = load(fullfile(exp_2,'/biogem/biogem_series_atm_pO2.res'),'ascii');
REF_sed_pO2_NoOMEN = load(fullfile(exp_3,'/biogem/biogem_series_atm_pO2.res'),'ascii');
REF_sed_pO2_SPIN = load(fullfile(exp_4,'/biogem/biogem_series_atm_pO2.res'),'ascii');
REF_sed_pO2_Exp_5 = load(fullfile(exp_1,'/biogem/biogem_series_atm_pO2.res'),'ascii');
% REF_sed_pO2_smallerk_DISS_HACK = load(fullfile(exp_2,'/biogem/biogem_series_atm_pO2.res'),'ascii');

% 
% % with fake sediments
% SED_NO_BURIAL = load('./cgenie_output/EXAMPLE101.p0093k.PO4Fe_S36x36.SPIN_2P4CO2_REF-SED_O2_rest_nothing_buried_1201/biogem/biogem_series_ocn_ALK.res','ascii');
% SED_NO_BURIAL_10 = load('./cgenie_output/EXAMPLE101.p0093k.PO4Fe_S36x36.SPIN_2P4CO2_REF-SED_O2_rest_nothing_buried_10_1201/biogem/biogem_series_ocn_ALK.res','ascii');
% 
% 
% SED_ALL_BURIAL = load('./cgenie_output/EXAMPLE101.p0093k.PO4Fe_S36x36.SPIN_2P4CO2_REF-SED_O2_rest_all_buried_1201/biogem/biogem_series_ocn_ALK.res','ascii');
% SED_ALL_BURIAL_10 = load('./cgenie_output/EXAMPLE101.p0093k.PO4Fe_S36x36.SPIN_2P4CO2_REF-SED_O2_rest_all_buried_10_1201/biogem/biogem_series_ocn_ALK.res','ascii');

set(0,'defaultLineLineWidth', 1)
set(0,'DefaultAxesFontSize',8)

figure
grid on
hold on
subplot(2, 1,1)
%plot(REF_sed_pCO2(:,1),REF_sed_pCO2(:,3)*1e+6,'b',REF_sed_pCO2_smallerk(:,1),REF_sed_pCO2_smallerk(:,3)*1e+6,'r' ); 
plot(REF_sed_pCO2_SPIN(:,1),REF_sed_pCO2_SPIN(:,3)*1e+6,'k', REF_sed_pCO2_NoOMEN(:,1),REF_sed_pCO2_NoOMEN(:,3)*1e+6,'g:', REF_sed_pCO2_allremin(:,1),REF_sed_pCO2_allremin(:,3)*1e+6,'b-.',REF_sed_pCO2_smallerk(:,1),REF_sed_pCO2_smallerk(:,3)*1e+6,'r-.',REF_sed_pCO2_Exp_5(:,1),REF_sed_pCO2_Exp_5(:,3)*1e+6,'m:'); 
%plot(REF_sed_pCO2_allremin(:,1),REF_sed_pCO2_allremin(:,3)*1e+6,'b',REF_sed_pCO2_smallerk(:,1),REF_sed_pCO2_smallerk(:,3)*1e+6,'r',REF_sed_pCO2_NoOMEN(:,1),REF_sed_pCO2_NoOMEN(:,3)*1e+6,'g:',REF_sed_pCO2_SPIN(:,1),REF_sed_pCO2_SPIN(:,3)*1e+6,'k:'); 
xlabel('time (years) ');
ylabel('global pCO2 (ppm)');
%x = [0,5000,10000,15000,20000,25000];
%set(gca,'XTick',x,'XTickLabel',sprintf('%5.0f|',x))
hleg=legend( 'CLOSED: with CaCO3 - No OMEN', 'OPEN: No CaCO3, No-OMEN', 'No CaCO3, OMEN: all remin','No CaCO3, OMEN: smaller k','No CaCO3, OMEN: all remin - DISS-HACK'); 
set(hleg,'Location','best')
%title('global pCO2 (ppm)','FontSize',18);

subplot(2,1,2)
%plot(REF(:,1),REF(:,3)*1e+6,'k')
%plot(REF_sed_pO2(:,1),REF_sed_pO2(:,3),'b',REF_sed_pO2_smallerk(:,1),REF_sed_pO2_smallerk(:,3),'r' );
plot(REF_sed_pO2_SPIN(:,1),REF_sed_pO2_SPIN(:,3),'k', REF_sed_pO2_NoOMEN(:,1),REF_sed_pO2_NoOMEN(:,3),'g:', REF_sed_pO2_allremin(:,1),REF_sed_pO2_allremin(:,3),'b-.',REF_sed_pO2_smallerk(:,1),REF_sed_pO2_smallerk(:,3),'r-.',REF_sed_pO2_Exp_5(:,1),REF_sed_pO2_Exp_5(:,3),'m:'); 
%plot(REF_sed_pO2_allremin(:,1),REF_sed_pO2_allremin(:,3),'b',REF_sed_pO2_smallerk(:,1),REF_sed_pO2_smallerk(:,3),'r',REF_sed_pO2_NoOMEN(:,1),REF_sed_pO2_NoOMEN(:,3),'g:',REF_sed_pO2_SPIN(:,1),REF_sed_pO2_SPIN(:,3),'k:'); 
xlabel('time (years) ');
ylabel('global pO2 (atm)');
ylim([0.20 0.21])
%ylim([0.207 0.2095])
%x = [0,5000,10000,15000,20000,25000];
%set(gca,'XTick',x,'XTickLabel',sprintf('%5.0f|',x))
% hleg=legend('small k', 'higher k', 'SPIN no Corg'); 
% set(hleg,'Location','NorthEast')

%print('-depsc', fullfile(exp_1,'/1_OUTPUT_PLOTS/0_Atmp-time-series'));
print('-depsc', 'cgenie_output/plots0212/0_Atmp-time-series_0212');

end