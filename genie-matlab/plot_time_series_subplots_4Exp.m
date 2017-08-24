%function [] = plot_time_series(PEXP1,PEXP2,PVAR1,PVAR2,PT1,PT2,PIK,PMASK,PCSCALE,PCMIN,PCMAX,PCN,PDATA,POPT,PNAME)
% plot time-series
clear all;

% plot mean (true) or total (false)
plot_mean = false;

% set experiment 
exp_1 = './cgenie_output/0606_02_EXAMPLE.worjh2.Archeretal2009.SPIN1_fastsinking';                                   
exp_2 = './cgenie_output/2308_01_EXAMPLE.worjh2.Archeretal2009.SPIN1_fastsinking_solidfields_OPEN';
exp_3 = './cgenie_output/2308_03_Archeretal2009_OMEN.inv_k2_0.005_k1_0.015_ord_3_solidfields_OPEN';
exp_4 = './cgenie_output/2308_05_Archeretal2009_OMEN.boudreau1997_k_depthdep_solidfields_OPEN';
% exp_2 = './cgenie_output/2308_07_EXAMPLE.worjh2.Archeretal2009.SPIN1_fastsinking_solidfields_ClosedCaCO3';
% exp_3 = './cgenie_output/2308_09_Archeretal2009_OMEN.inv_k2_0.005_k1_0.015_ord_3_solidfields_ClosedCaCO3';
% exp_4 = './cgenie_output/2308_11_Archeretal2009_OMEN.boudreau1997_k_depthdep_solidfields_ClosedCaCO3';
% %%%% load other data

REF_sed_O2_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_O2.res'),'ascii');
REF_sed_O2_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_O2.res'),'ascii');
REF_sed_O2_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_O2.res'),'ascii');
REF_sed_O2_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_ocn_O2.res'),'ascii');

REF_sed_SO4_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_SO4.res'),'ascii');
REF_sed_SO4_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_SO4.res'),'ascii');
REF_sed_SO4_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_SO4.res'),'ascii');
REF_sed_SO4_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_ocn_SO4.res'),'ascii');

REF_sed_H2S_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_H2S.res'),'ascii');
REF_sed_H2S_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_H2S.res'),'ascii');
REF_sed_H2S_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_H2S.res'),'ascii');
REF_sed_H2S_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_ocn_H2S.res'),'ascii');

REF_sed_PO4_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_PO4.res'),'ascii');
REF_sed_PO4_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_PO4.res'),'ascii');
REF_sed_PO4_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_PO4.res'),'ascii');
REF_sed_PO4_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_ocn_PO4.res'),'ascii');

REF_sed_ALK_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_ALK.res'),'ascii');
REF_sed_ALK_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_ALK.res'),'ascii');
REF_sed_ALK_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_ALK.res'),'ascii');
REF_sed_ALK_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_ocn_ALK.res'),'ascii');

REF_sed_DIC_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_DIC.res'),'ascii');
REF_sed_DIC_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_DIC.res'),'ascii');
REF_sed_DIC_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_DIC.res'),'ascii');
REF_sed_DIC_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_ocn_DIC.res'),'ascii');

set(0,'defaultLineLineWidth', 2)
set(0,'DefaultAxesFontSize', 10)

figure
grid on
hold on

subplot(3,2,1)
if(plot_mean) % mean (mol/kg)
    plot(REF_sed_O2_exp1(:,1),REF_sed_O2_exp1(:,3)*1e+6,'b',REF_sed_O2_exp2(:,1),REF_sed_O2_exp2(:,3)*1e+6,'r--',REF_sed_O2_exp3(:,1),REF_sed_O2_exp3(:,3)*1e+6,'g--',REF_sed_O2_exp4(:,1),REF_sed_O2_exp4(:,3)*1e+6,'k:'); 
    ylabel('O_2 (\mumol kg^{-1})');
else
    % total (mol)
    plot(REF_sed_O2_exp1(:,1),REF_sed_O2_exp1(:,2),'b',REF_sed_O2_exp2(:,1),REF_sed_O2_exp2(:,2),'ro',REF_sed_O2_exp3(:,1),REF_sed_O2_exp3(:,2),'g--',REF_sed_O2_exp4(:,1),REF_sed_O2_exp4(:,2),'k:'); 
    ylabel('O_2 (mol)');
end
xlabel('yrs ');
%hleg=legend('Abiotic - No OMEN', 'Abiotic - with OMEN', 'Biotic - all remin. k=0.1', 'Biotic - smaller k'); 
xlim([0 10000])
%set(hleg,'FontSize',4);
%set(hleg,'Location','best')


subplot(3,2,2)
if(plot_mean)  % mean (mol/kg)
	plot(REF_sed_SO4_exp1(:,1),REF_sed_SO4_exp1(:,3)*1e+6,'b',REF_sed_SO4_exp2(:,1),REF_sed_SO4_exp2(:,3)*1e+6,'r--',REF_sed_SO4_exp3(:,1),REF_sed_SO4_exp3(:,3)*1e+6,'g--',REF_sed_SO4_exp4(:,1),REF_sed_SO4_exp4(:,3)*1e+6,'k:' ); 
    ylabel('SO_4 (\mumol kg^{-1})');
else    % total (mol)
    plot(REF_sed_SO4_exp1(:,1),REF_sed_SO4_exp1(:,2),'b',REF_sed_SO4_exp2(:,1),REF_sed_SO4_exp2(:,2),'ro',REF_sed_SO4_exp3(:,1),REF_sed_SO4_exp3(:,2),'g--',REF_sed_SO4_exp4(:,1),REF_sed_SO4_exp4(:,2),'k:' ); 
    ylabel('SO_4 (mol)');
end
xlabel('yrs ');
xlim([0 10000])
hleg=legend('0606-02 Archer SPIN - No OMEN - fast sinking', '2308-01 Archer SPIN - No OMEN - fast sinking - solid fields', '2308-03 with OMEN - solid fields - invariant k1 = 0.015', '2308-05  with OMEN - solid fields - Boudreau depth'); 
%hleg=legend('0606-02 Archer SPIN - No OMEN - fast sinking', '2308-07 Archer SPIN - No OMEN - fast sinking - solid fields closedCaCO3', '2308-09 with OMEN - solid fields - invariant k1 = 0.015 closedCaCO3', '2308-11  with OMEN - solid fields - Boudreau depth closedCaCO3'); 
set(hleg,'FontSize',4);
set(hleg,'Location','SouthEast');

subplot(3,2,3)
if(plot_mean) % mean (mol/kg)
	plot(REF_sed_H2S_exp1(:,1),REF_sed_H2S_exp1(:,3)*1e+6,'b',REF_sed_H2S_exp2(:,1),REF_sed_H2S_exp2(:,3)*1e+6,'r--', REF_sed_H2S_exp3(:,1),REF_sed_H2S_exp3(:,3)*1e+6,'g--', REF_sed_H2S_exp3(:,1),REF_sed_H2S_exp3(:,3)*1e+6,'g--',REF_sed_H2S_exp4(:,1),REF_sed_H2S_exp4(:,3)*1e+6,'k:');
    ylabel('H_2S (\mumol kg^{-1})');
else % total (mol)
    plot(REF_sed_H2S_exp1(:,1),REF_sed_H2S_exp1(:,2),'b',REF_sed_H2S_exp2(:,1),REF_sed_H2S_exp2(:,2),'ro', REF_sed_H2S_exp3(:,1),REF_sed_H2S_exp3(:,2),'g--', REF_sed_H2S_exp3(:,1),REF_sed_H2S_exp3(:,2),'g--',REF_sed_H2S_exp4(:,1),REF_sed_H2S_exp4(:,2),'k:');
    ylabel('H_2S (mol)');
end
xlabel('yrs ');
xlim([0 10000])
%hleg=legend('small k', 'higher k'); 
%set(hleg,'Location','SouthEast')

subplot(3,2,4)
if(plot_mean) % mean (mol/kg)
    plot(REF_sed_PO4_exp1(:,1),REF_sed_PO4_exp1(:,3)*1e+6,'b',REF_sed_PO4_exp2(:,1),REF_sed_PO4_exp2(:,3)*1e+6,'r--', REF_sed_PO4_exp3(:,1),REF_sed_PO4_exp3(:,3)*1e+6,'g--',REF_sed_PO4_exp4(:,1),REF_sed_PO4_exp4(:,3)*1e+6,'k:'); 
    ylabel('PO_4 (\mumol kg^{-1})');
else % total (mol)
    plot(REF_sed_PO4_exp1(:,1),REF_sed_PO4_exp1(:,2),'b',REF_sed_PO4_exp2(:,1),REF_sed_PO4_exp2(:,2),'ro', REF_sed_PO4_exp3(:,1),REF_sed_PO4_exp3(:,2),'g--', REF_sed_PO4_exp3(:,1),REF_sed_PO4_exp3(:,2),'g--',REF_sed_PO4_exp4(:,1),REF_sed_PO4_exp4(:,2),'k:');
    ylabel('PO_4 (mol)');
end
xlabel('yrs ');
xlim([0 10000])
%ylim([2.15 2.17])
% hleg=legend('small k', 'higher k', 'exp4 no Corg'); 
% set(hleg,'Location','NorthEast')

subplot(3,2,5)
if(plot_mean)  % mean (mol/kg)
    plot(REF_sed_ALK_exp1(:,1),REF_sed_ALK_exp1(:,3)*1e+6,'b',REF_sed_ALK_exp2(:,1),REF_sed_ALK_exp2(:,3)*1e+6,'r--', REF_sed_ALK_exp3(:,1),REF_sed_ALK_exp3(:,3)*1e+6,'g--',REF_sed_ALK_exp4(:,1),REF_sed_ALK_exp4(:,3)*1e+6,'k:'); 
    ylabel('ALK (\mumol kg^{-1})');
else
	plot(REF_sed_ALK_exp1(:,1),REF_sed_ALK_exp1(:,2),'b',REF_sed_ALK_exp2(:,1),REF_sed_ALK_exp2(:,2),'ro', REF_sed_ALK_exp3(:,1),REF_sed_ALK_exp3(:,2),'g--', REF_sed_ALK_exp3(:,1),REF_sed_ALK_exp3(:,2),'g--',REF_sed_ALK_exp4(:,1),REF_sed_ALK_exp4(:,2),'k:');
    ylabel('ALK (mol)');
end
xlabel('yrs ');
xlim([0 10000])
%hleg=legend('small k', 'higher k'); 
%set(hleg,'FontSize',10)
%set(hleg,'Location','SouthEast')
%title('oceanic Alkalinity','FontSize',18);

subplot(3,2,6)
if(plot_mean)  % mean (mol/kg)
    plot(REF_sed_DIC_exp1(:,1),REF_sed_DIC_exp1(:,3)*1e+6,'b',REF_sed_DIC_exp2(:,1),REF_sed_DIC_exp2(:,3)*1e+6,'r--', REF_sed_DIC_exp3(:,1),REF_sed_DIC_exp3(:,3)*1e+6,'g--',REF_sed_DIC_exp4(:,1),REF_sed_DIC_exp4(:,3)*1e+6,'k:'); 
    ylabel('DIC (\mumol kg^{-1})');
else
	plot(REF_sed_DIC_exp1(:,1),REF_sed_DIC_exp1(:,2),'b',REF_sed_DIC_exp2(:,1),REF_sed_DIC_exp2(:,2),'ro', REF_sed_DIC_exp3(:,1),REF_sed_DIC_exp3(:,2),'g--', REF_sed_DIC_exp3(:,1),REF_sed_DIC_exp3(:,2),'g--',REF_sed_DIC_exp4(:,1),REF_sed_DIC_exp4(:,2),'k:');
    ylabel('DIC (mol)');
end
xlabel('yrs ');
xlim([0 10000])
% hleg=legend('small k', 'higher k', 'exp4 no Corg'); 
% set(hleg,'Location','SouthEast')

% set(gcf,'NextPlot','add');
% axes;
% h = title('Global ocean values (micromol kg**-1)');
% set(gca,'Visible','off');
% set(h,'Visible','on'); 

% print('-depsc', fullfile(exp_1,'/1_OUTPUT_PLOTS/0_ALL-time-series'));
if(plot_mean)  % mean (mol/kg)
    print('-depsc', 'cgenie_output/0_PLOTS/0606_Archer_Boudreau/0606_ArcherExp_also_withBoudreau_k');
else
    print('-depsc', 'cgenie_output/000_FOR_GMD_V1608_1908/2308_Timeseries_SPIN_and_OMEN_NoPO4_with_solid_fields');
%    print('-depsc', 'cgenie_output/000_FOR_GMD_V1608_1908/2308_Timeseries_SPIN_and_OMEN_NoPO4_with_solid_fields_ClosedCaCO3');
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% atm O2 and CO2
% %%%% load other data
% REF experiments

if(false)
REF_sed_pCO2_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
REF_sed_pCO2_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
REF_sed_pCO2_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
REF_sed_pCO2_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_atm_pCO2.res'),'ascii');

REF_sed_pO2_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_atm_pO2.res'),'ascii');
REF_sed_pO2_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_atm_pO2.res'),'ascii');
REF_sed_pO2_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_atm_pO2.res'),'ascii');
REF_sed_pO2_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_atm_pO2.res'),'ascii');

% 
% % with fake sediments
% SED_NO_BURIAL = load('./cgenie_output/EXAMPLE101.p0093k.PO4Fe_S36x36.exp4_2P4CO2_REF-SED_O2_rest_nothing_buried_1201/biogem/biogem_series_ocn_ALK.res','ascii');
% SED_NO_BURIAL_10 = load('./cgenie_output/EXAMPLE101.p0093k.PO4Fe_S36x36.exp4_2P4CO2_REF-SED_O2_rest_nothing_buried_10_1201/biogem/biogem_series_ocn_ALK.res','ascii');
% 
% 
% SED_ALL_BURIAL = load('./cgenie_output/EXAMPLE101.p0093k.PO4Fe_S36x36.exp4_2P4CO2_REF-SED_O2_rest_all_buried_1201/biogem/biogem_series_ocn_ALK.res','ascii');
% SED_ALL_BURIAL_10 = load('./cgenie_output/EXAMPLE101.p0093k.PO4Fe_S36x36.exp4_2P4CO2_REF-SED_O2_rest_all_buried_10_1201/biogem/biogem_series_ocn_ALK.res','ascii');

set(0,'defaultLineLineWidth', 1)
set(0,'DefaultAxesFontSize',8)

figure
grid on
hold on
subplot(2, 1,1)
%plot(REF_sed_pCO2(:,1),REF_sed_pCO2(:,3)*1e+6,'b',REF_sed_pCO2_exp2(:,1),REF_sed_pCO2_exp2(:,3)*1e+6,'r' ); 
plot(REF_sed_pCO2_exp1(:,1),REF_sed_pCO2_exp1(:,3)*1e+6,'b',REF_sed_pCO2_exp2(:,1),REF_sed_pCO2_exp2(:,3)*1e+6,'r--', REF_sed_pCO2_exp3(:,1),REF_sed_pCO2_exp3(:,3)*1e+6,'g--', REF_sed_pCO2_exp4(:,1),REF_sed_pCO2_exp4(:,3)*1e+6,'k:' ); 
%plot(REF_sed_pCO2_exp1(:,1),REF_sed_pCO2_exp1(:,3)*1e+6,'b',REF_sed_pCO2_exp2(:,1),REF_sed_pCO2_exp2(:,3)*1e+6,'r',REF_sed_pCO2_exp3(:,1),REF_sed_pCO2_exp3(:,3)*1e+6,'g--',REF_sed_pCO2_exp4(:,1),REF_sed_pCO2_exp4(:,3)*1e+6,'k:'); 
xlabel('time (years) ');
ylabel('global pCO2 (ppm)');
%x = [0,5000,10000,15000,20000,25000];
%set(gca,'XTick',x,'XTickLabel',sprintf('%5.0f|',x))
hleg=legend('Abiotic - No OMEN', 'Abiotic - with OMEN', 'Biotic - all remin. k=0.1', 'Biotic - smaller k'); ylim([0 300])
set(hleg,'Location','best')
%title('global pCO2 (ppm)','FontSize',18);

subplot(2,1,2)
%plot(REF(:,1),REF(:,3)*1e+6,'k')
%plot(REF_sed_pO2(:,1),REF_sed_pO2(:,3),'b',REF_sed_pO2_exp2(:,1),REF_sed_pO2_exp2(:,3),'r' );
plot(REF_sed_pO2_exp1(:,1),REF_sed_pO2_exp1(:,3),'b',REF_sed_pO2_exp2(:,1),REF_sed_pO2_exp2(:,3),'r--',REF_sed_pO2_exp3(:,1),REF_sed_pO2_exp3(:,3),'g--',REF_sed_pO2_exp4(:,1),REF_sed_pO2_exp4(:,3),'k:' );
%plot(REF_sed_pO2_exp1(:,1),REF_sed_pO2_exp1(:,3),'b',REF_sed_pO2_exp2(:,1),REF_sed_pO2_exp2(:,3),'r',REF_sed_pO2_exp3(:,1),REF_sed_pO2_exp3(:,3),'g--',REF_sed_pO2_exp4(:,1),REF_sed_pO2_exp4(:,3),'k:'); 
xlabel('time (years) ');
ylabel('global pO2 (atm)');
ylim([0.20 0.21])
%ylim([0.207 0.2095])
%x = [0,5000,10000,15000,20000,25000];
%set(gca,'XTick',x,'XTickLabel',sprintf('%5.0f|',x))
% hleg=legend('small k', 'higher k', 'exp4 no Corg'); 
% set(hleg,'Location','NorthEast')

%print('-depsc', fullfile(exp_1,'/1_OUTPUT_PLOTS/0_Atmp-time-series'));
print('-depsc', 'cgenie_output/0_PLOTS/plots_0602/0_Atmp-time-series_abiotic_biotic_0602');
end