%function [] = plot_time_series(PEXP1,PEXP2,PVAR1,PVAR2,PT1,PT2,PIK,PMASK,PCSCALE,PCMIN,PCMAX,PCN,PDATA,POPT,PNAME)
% plot time-series
clear all;

% set experiment 
exp_1 = './cgenie_output/EXAMPLE05._rwlma.PO4_S18x18.OPEN_noCaCO3_Corg_all_remin_2311';
exp_2 = './cgenie_output/EXAMPLE05._rwlma.PO4_S18x18.OPEN_noCaCO3_Corg_ksmaller_2311';
exp_3 = './cgenie_output/EXAMPLE05._rwlma.PO4_S18x18.OPEN_noCaCO3_2311';
exp_4 = './cgenie_output/EXAMPLE._rwlma.PO4_S18x18.SPIN1';
% %%%% load other data

REF_sed_O2_allremin = load(fullfile(exp_1,'/biogem/biogem_series_ocn_O2.res'),'ascii');
REF_sed_O2_smallerk = load(fullfile(exp_2,'/biogem/biogem_series_ocn_O2.res'),'ascii');
REF_sed_O2_NoOMEN = load(fullfile(exp_3,'/biogem/biogem_series_ocn_O2.res'),'ascii');
REF_sed_O2_SPIN = load(fullfile(exp_4,'/biogem/biogem_series_ocn_O2.res'),'ascii');

REF_sed_SO4_allremin = load(fullfile(exp_1,'/biogem/biogem_series_ocn_SO4.res'),'ascii');
REF_sed_SO4_smallerk = load(fullfile(exp_2,'/biogem/biogem_series_ocn_SO4.res'),'ascii');
REF_sed_SO4_NoOMEN = load(fullfile(exp_3,'/biogem/biogem_series_ocn_SO4.res'),'ascii');
REF_sed_SO4_SPIN = load(fullfile(exp_4,'/biogem/biogem_series_ocn_SO4.res'),'ascii');

REF_sed_H2S_allremin = load(fullfile(exp_1,'/biogem/biogem_series_ocn_H2S.res'),'ascii');
REF_sed_H2S_smallerk = load(fullfile(exp_2,'/biogem/biogem_series_ocn_H2S.res'),'ascii');
REF_sed_H2S_NoOMEN = load(fullfile(exp_3,'/biogem/biogem_series_ocn_H2S.res'),'ascii');
REF_sed_H2S_SPIN = load(fullfile(exp_4,'/biogem/biogem_series_ocn_H2S.res'),'ascii');

REF_sed_PO4_allremin = load(fullfile(exp_1,'/biogem/biogem_series_ocn_PO4.res'),'ascii');
REF_sed_PO4_smallerk = load(fullfile(exp_2,'/biogem/biogem_series_ocn_PO4.res'),'ascii');
REF_sed_PO4_NoOMEN = load(fullfile(exp_3,'/biogem/biogem_series_ocn_PO4.res'),'ascii');
REF_sed_PO4_SPIN = load(fullfile(exp_4,'/biogem/biogem_series_ocn_PO4.res'),'ascii');

REF_sed_ALK_allremin = load(fullfile(exp_1,'/biogem/biogem_series_ocn_ALK.res'),'ascii');
REF_sed_ALK_smallerk = load(fullfile(exp_2,'/biogem/biogem_series_ocn_ALK.res'),'ascii');
REF_sed_ALK_NoOMEN = load(fullfile(exp_3,'/biogem/biogem_series_ocn_ALK.res'),'ascii');
REF_sed_ALK_SPIN = load(fullfile(exp_4,'/biogem/biogem_series_ocn_ALK.res'),'ascii');

REF_sed_DIC_allremin = load(fullfile(exp_1,'/biogem/biogem_series_ocn_DIC.res'),'ascii');
REF_sed_DIC_smallerk = load(fullfile(exp_2,'/biogem/biogem_series_ocn_DIC.res'),'ascii');
REF_sed_DIC_NoOMEN = load(fullfile(exp_3,'/biogem/biogem_series_ocn_DIC.res'),'ascii');
REF_sed_DIC_SPIN = load(fullfile(exp_4,'/biogem/biogem_series_ocn_DIC.res'),'ascii');

set(0,'defaultLineLineWidth', 1)
set(0,'DefaultAxesFontSize',8)

figure
grid on
hold on

subplot(3,2,1)
%plot(REF_sed_O2(:,1),REF_sed_O2(:,3)*1e+6,'b',REF_sed_O2_smallerk(:,1),REF_sed_O2_smallerk(:,3)*1e+6,'r' ); 
plot(REF_sed_O2_allremin(:,1),REF_sed_O2_allremin(:,3)*1e+6,'b',REF_sed_O2_smallerk(:,1),REF_sed_O2_smallerk(:,3)*1e+6,'r',REF_sed_O2_NoOMEN(:,1),REF_sed_O2_NoOMEN(:,3)*1e+6,'g--',REF_sed_O2_SPIN(:,1),REF_sed_O2_SPIN(:,3)*1e+6,'k:'); 
%plot(REF_sed_O2_allremin(:,1),REF_sed_O2_allremin(:,3)*1e+6,'b',REF_sed_O2_smallerk(:,1),REF_sed_O2_smallerk(:,3)*1e+6,'r',REF_sed_O2_NoOMEN(:,1),REF_sed_O2_NoOMEN(:,3)*1e+6,'g--',REF_sed_O2_SPIN(:,1),REF_sed_O2_SPIN(:,3)*1e+6,'k:'); 
xlabel('yrs ');
ylabel('O_2 (\mumol kg^{-1})');
hleg=legend('OPEN: No CaCO3, OMEN: all remin','OPEN: No CaCO3, OMEN: smaller k','OPEN: No CaCO3, No-OMEN', 'SPIN: Closed with CaCO3 - No OMEN'); 
ylim([0 300])
set(hleg,'FontSize',4);
set(hleg,'Location','best')


subplot(3,2,2)
%plot(REF_sed_SO4(:,1),REF_sed_SO4(:,3)*1e+6,'b',REF_sed_SO4_smallerk(:,1),REF_sed_SO4_smallerk(:,3)*1e+6,'r' ); 
plot(REF_sed_SO4_allremin(:,1),REF_sed_SO4_allremin(:,3)*1e+6,'b',REF_sed_SO4_smallerk(:,1),REF_sed_SO4_smallerk(:,3)*1e+6,'r',REF_sed_SO4_NoOMEN(:,1),REF_sed_SO4_NoOMEN(:,3)*1e+6,'g--',REF_sed_SO4_NoOMEN(:,1),REF_sed_SO4_NoOMEN(:,3)*1e+6,'g--',REF_sed_SO4_SPIN(:,1),REF_sed_SO4_SPIN(:,3)*1e+6,'k:' ); 
xlabel('yrs ');
ylabel('SO_4 (\mumol kg^{-1})');
hleg=legend('OPEN: No CaCO3, OMEN: all remin','OPEN: No CaCO3, OMEN: smaller k','OPEN: No CaCO3, No-OMEN', 'SPIN: Closed with CaCO3 - No OMEN'); 
set(hleg,'FontSize',4);
set(hleg,'Location','best')

subplot(3,2,3)
%plot(REF_sed_H2S(:,1),REF_sed_H2S(:,3)*1e+6,'b', REF_sed_H2S_smallerk(:,1),REF_sed_H2S_smallerk(:,3)*1e+6,'r');
plot(REF_sed_H2S_allremin(:,1),REF_sed_H2S_allremin(:,3)*1e+6,'b',REF_sed_H2S_smallerk(:,1),REF_sed_H2S_smallerk(:,3)*1e+6,'r', REF_sed_H2S_NoOMEN(:,1),REF_sed_H2S_NoOMEN(:,3)*1e+6,'g--', REF_sed_H2S_NoOMEN(:,1),REF_sed_H2S_NoOMEN(:,3)*1e+6,'g--',REF_sed_H2S_SPIN(:,1),REF_sed_H2S_SPIN(:,3)*1e+6,'k:');
%plot(REF_sed_H2S_allremin(:,1),REF_sed_H2S_allremin(:,3)*1e+6,'b',REF_sed_H2S_smallerk(:,1),REF_sed_H2S_smallerk(:,3)*1e+6,'r',REF_sed_H2S_NoOMEN(:,1),REF_sed_H2S_NoOMEN(:,3)*1e+6,'g--',REF_sed_H2S_SPIN(:,1),REF_sed_H2S_SPIN(:,3)*1e+6,'k:'); 
xlabel('yrs ');
ylabel('H_2S (\mumol kg^{-1})');
%hleg=legend('small k', 'higher k'); 
%set(hleg,'Location','SouthEast')

subplot(3,2,4)
%plot(REF_sed_PO4(:,1),REF_sed_PO4(:,3)*1e+6,'b', REF_sed_PO4_smallerk(:,1),REF_sed_PO4_smallerk(:,3)*1e+6,'r'); 
plot(REF_sed_PO4_allremin(:,1),REF_sed_PO4_allremin(:,3)*1e+6,'b',REF_sed_PO4_smallerk(:,1),REF_sed_PO4_smallerk(:,3)*1e+6,'r', REF_sed_PO4_NoOMEN(:,1),REF_sed_PO4_NoOMEN(:,3)*1e+6,'g--',REF_sed_PO4_SPIN(:,1),REF_sed_PO4_SPIN(:,3)*1e+6,'k:'); 
%plot(REF_sed_PO4_allremin(:,1),REF_sed_PO4_allremin(:,3)*1e+6,'b',REF_sed_PO4_smallerk(:,1),REF_sed_PO4_smallerk(:,3)*1e+6,'r',REF_sed_PO4_NoOMEN(:,1),REF_sed_PO4_NoOMEN(:,3)*1e+6,'g--',REF_sed_PO4_SPIN(:,1),REF_sed_PO4_SPIN(:,3)*1e+6,'k:'); 
xlabel('yrs ');
ylabel('PO_4 (\mumol kg^{-1})');
%ylim([2.15 2.17])
% hleg=legend('small k', 'higher k', 'SPIN no Corg'); 
% set(hleg,'Location','NorthEast')

subplot(3,2,5)
%plot(REF_sed_ALK(:,1),REF_sed_ALK(:,3)*1e+6,'b', REF_sed_ALK_smallerk(:,1),REF_sed_ALK_smallerk(:,3)*1e+6,'r'); 
plot(REF_sed_ALK_allremin(:,1),REF_sed_ALK_allremin(:,3)*1e+6,'b',REF_sed_ALK_smallerk(:,1),REF_sed_ALK_smallerk(:,3)*1e+6,'r', REF_sed_ALK_NoOMEN(:,1),REF_sed_ALK_NoOMEN(:,3)*1e+6,'g--',REF_sed_ALK_SPIN(:,1),REF_sed_ALK_SPIN(:,3)*1e+6,'k:'); 
%plot(REF_sed_ALK_allremin(:,1),REF_sed_ALK_allremin(:,3)*1e+6,'b',REF_sed_ALK_smallerk(:,1),REF_sed_ALK_smallerk(:,3)*1e+6,'r',REF_sed_ALK_NoOMEN(:,1),REF_sed_ALK_NoOMEN(:,3)*1e+6,'g--',REF_sed_ALK_SPIN(:,1),REF_sed_ALK_SPIN(:,3)*1e+6,'k:'); 
xlabel('yrs ');
ylabel('ALK (\mumol kg^{-1})');
%hleg=legend('small k', 'higher k'); 
%set(hleg,'FontSize',10)
%set(hleg,'Location','SouthEast')
%title('oceanic Alkalinity','FontSize',18);

subplot(3,2,6)
%plot(REF_sed_DIC(:,1),REF_sed_DIC(:,3)*1e+6,'b', REF_sed_DIC_smallerk(:,1),REF_sed_DIC_smallerk(:,3)*1e+6,'r'); 
plot(REF_sed_DIC_allremin(:,1),REF_sed_DIC_allremin(:,3)*1e+6,'b',REF_sed_DIC_smallerk(:,1),REF_sed_DIC_smallerk(:,3)*1e+6,'r', REF_sed_DIC_NoOMEN(:,1),REF_sed_DIC_NoOMEN(:,3)*1e+6,'g--',REF_sed_DIC_SPIN(:,1),REF_sed_DIC_SPIN(:,3)*1e+6,'k:'); 
%plot(REF_sed_DIC_allremin(:,1),REF_sed_DIC_allremin(:,3)*1e+6,'b',REF_sed_DIC_smallerk(:,1),REF_sed_DIC_smallerk(:,3)*1e+6,'r',REF_sed_DIC_NoOMEN(:,1),REF_sed_DIC_NoOMEN(:,3)*1e+6,'g--',REF_sed_DIC_SPIN(:,1),REF_sed_DIC_SPIN(:,3)*1e+6,'k:'); 
xlabel('yrs ');
ylabel('DIC (\mumol kg^{-1})');
% hleg=legend('small k', 'higher k', 'SPIN no Corg'); 
% set(hleg,'Location','SouthEast')

% set(gcf,'NextPlot','add');
% axes;
% h = title('Global ocean values (micromol kg**-1)');
% set(gca,'Visible','off');
% set(h,'Visible','on'); 

% print('-depsc', fullfile(exp_1,'/1_OUTPUT_PLOTS/0_ALL-time-series'));
print('-depsc', 'cgenie_output/plots2311_NEWDissflux/0_ALL-time-series_2311');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% atm O2 and CO2
% %%%% load other data
% REF experiments

REF_sed_pCO2_allremin = load(fullfile(exp_1,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
REF_sed_pCO2_smallerk = load(fullfile(exp_2,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
REF_sed_pCO2_NoOMEN = load(fullfile(exp_3,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
REF_sed_pCO2_SPIN = load(fullfile(exp_4,'/biogem/biogem_series_atm_pCO2.res'),'ascii');

REF_sed_pO2_allremin = load(fullfile(exp_1,'/biogem/biogem_series_atm_pO2.res'),'ascii');
REF_sed_pO2_smallerk = load(fullfile(exp_2,'/biogem/biogem_series_atm_pO2.res'),'ascii');
REF_sed_pO2_NoOMEN = load(fullfile(exp_3,'/biogem/biogem_series_atm_pO2.res'),'ascii');
REF_sed_pO2_SPIN = load(fullfile(exp_4,'/biogem/biogem_series_atm_pO2.res'),'ascii');

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
plot(REF_sed_pCO2_allremin(:,1),REF_sed_pCO2_allremin(:,3)*1e+6,'b',REF_sed_pCO2_smallerk(:,1),REF_sed_pCO2_smallerk(:,3)*1e+6,'r', REF_sed_pCO2_NoOMEN(:,1),REF_sed_pCO2_NoOMEN(:,3)*1e+6,'g--', REF_sed_pCO2_SPIN(:,1),REF_sed_pCO2_SPIN(:,3)*1e+6,'k:' ); 
%plot(REF_sed_pCO2_allremin(:,1),REF_sed_pCO2_allremin(:,3)*1e+6,'b',REF_sed_pCO2_smallerk(:,1),REF_sed_pCO2_smallerk(:,3)*1e+6,'r',REF_sed_pCO2_NoOMEN(:,1),REF_sed_pCO2_NoOMEN(:,3)*1e+6,'g--',REF_sed_pCO2_SPIN(:,1),REF_sed_pCO2_SPIN(:,3)*1e+6,'k:'); 
xlabel('time (years) ');
ylabel('global pCO2 (ppm)');
%x = [0,5000,10000,15000,20000,25000];
%set(gca,'XTick',x,'XTickLabel',sprintf('%5.0f|',x))
hleg=legend('OPEN: No CaCO3, OMEN: all remin','OPEN: No CaCO3, OMEN: smaller k','OPEN: No CaCO3, No-OMEN', 'SPIN: Closed with CaCO3 - No OMEN'); 
set(hleg,'Location','best')
%title('global pCO2 (ppm)','FontSize',18);

subplot(2,1,2)
%plot(REF(:,1),REF(:,3)*1e+6,'k')
%plot(REF_sed_pO2(:,1),REF_sed_pO2(:,3),'b',REF_sed_pO2_smallerk(:,1),REF_sed_pO2_smallerk(:,3),'r' );
plot(REF_sed_pO2_allremin(:,1),REF_sed_pO2_allremin(:,3),'b',REF_sed_pO2_smallerk(:,1),REF_sed_pO2_smallerk(:,3),'r',REF_sed_pO2_NoOMEN(:,1),REF_sed_pO2_NoOMEN(:,3),'g--',REF_sed_pO2_SPIN(:,1),REF_sed_pO2_SPIN(:,3),'k:' );
%plot(REF_sed_pO2_allremin(:,1),REF_sed_pO2_allremin(:,3),'b',REF_sed_pO2_smallerk(:,1),REF_sed_pO2_smallerk(:,3),'r',REF_sed_pO2_NoOMEN(:,1),REF_sed_pO2_NoOMEN(:,3),'g--',REF_sed_pO2_SPIN(:,1),REF_sed_pO2_SPIN(:,3),'k:'); 
xlabel('time (years) ');
ylabel('global pO2 (atm)');
ylim([0.20 0.21])
%ylim([0.207 0.2095])
%x = [0,5000,10000,15000,20000,25000];
%set(gca,'XTick',x,'XTickLabel',sprintf('%5.0f|',x))
% hleg=legend('small k', 'higher k', 'SPIN no Corg'); 
% set(hleg,'Location','NorthEast')

%print('-depsc', fullfile(exp_1,'/1_OUTPUT_PLOTS/0_Atmp-time-series'));
print('-depsc', 'cgenie_output/plots2311_NEWDissflux/0_Atmp-time-series_2311');