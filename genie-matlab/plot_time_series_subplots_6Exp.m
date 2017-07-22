%function [] = plot_time_series(PEXP1,PEXP2,PVAR1,PVAR2,PT1,PT2,PIK,PMASK,PCSCALE,PCMIN,PCMAX,PCN,PDATA,POPT,PNAME)
% plot time-series
clear all;

% plot mean (true) or total (false)
plot_mean = false;

% set experiment 
%exp_1 = './cgenie_output/2_Feb2017/1602_09_art_world__nogas_noweath_WITH2ndRedox+ALK_k_0.1_abiotic';
%exp_1 = './cgenie_output/02_OMEN_GENIE_GMD/0903_01_art_world_shelves_nogas_noweath_No2ndRedox_NoALK_k_0.1_abiotic';
% % 16 level ocean
% exp_1 = './cgenie_output/02_OMEN_GENIE_GMD/2103_00_worjh2_nogas_noweather_open_NoOMEN';
% exp_2 = './cgenie_output/02_OMEN_GENIE_GMD/2103_01_worjh2_nogas_noweather_open_k_0.01';  
% exp_3 = './cgenie_output/02_OMEN_GENIE_GMD/2103_02_worjh2_nogas_noweather_open_k_0.05';
% exp_4 = './cgenie_output/02_OMEN_GENIE_GMD/2103_03_worjh2_nogas_noweather_open_k_0.1';
% exp_5 = './cgenie_output/02_OMEN_GENIE_GMD/2103_04_worjh2_nogas_noweather_open_k_1.0';
% exp_6 = './cgenie_output/02_OMEN_GENIE_GMD/2103_05_worjh2_nogas_noweather_open_k_4.0';
% %%%% load other data
% 8 level ocean
exp_1 = './cgenie_output/EXAMPLE.worjh2.Archeretal2009.SPIN1';                                   
exp_2 = './cgenie_output/1707_26_Archeretal2009_OMEN.inv_k2_0.007_k1_0.014_ord_2';
exp_3 = './cgenie_output/1707_27_Archeretal2009_OMEN.inv_k2_0.007_k1_0.021_ord_3';
exp_4 = './cgenie_output/1707_28_Archeretal2009_OMEN.inv_k2_0.007_k1_0.028_ord_4';
exp_5 = './cgenie_output/1707_29_Archeretal2009_OMEN.inv_k2_0.007_k1_0.035_ord_5';
exp_6 = './cgenie_output/1707_30_Archeretal2009_OMEN.inv_k2_0.007_k1_0.07_ord_10';

REF_sed_O2_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_O2.res'),'ascii');
REF_sed_O2_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_O2.res'),'ascii');
REF_sed_O2_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_O2.res'),'ascii');
REF_sed_O2_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_ocn_O2.res'),'ascii');
REF_sed_O2_exp5 = load(fullfile(exp_5,'/biogem/biogem_series_ocn_O2.res'),'ascii');
REF_sed_O2_exp6 = load(fullfile(exp_6,'/biogem/biogem_series_ocn_O2.res'),'ascii');

REF_sed_SO4_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_SO4.res'),'ascii');
REF_sed_SO4_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_SO4.res'),'ascii');
REF_sed_SO4_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_SO4.res'),'ascii');
REF_sed_SO4_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_ocn_SO4.res'),'ascii');
REF_sed_SO4_exp5 = load(fullfile(exp_5,'/biogem/biogem_series_ocn_SO4.res'),'ascii');
REF_sed_SO4_exp6 = load(fullfile(exp_6,'/biogem/biogem_series_ocn_SO4.res'),'ascii');

REF_sed_H2S_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_H2S.res'),'ascii');
REF_sed_H2S_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_H2S.res'),'ascii');
REF_sed_H2S_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_H2S.res'),'ascii');
REF_sed_H2S_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_ocn_H2S.res'),'ascii');
REF_sed_H2S_exp5 = load(fullfile(exp_5,'/biogem/biogem_series_ocn_H2S.res'),'ascii');
REF_sed_H2S_exp6 = load(fullfile(exp_6,'/biogem/biogem_series_ocn_H2S.res'),'ascii');

REF_sed_PO4_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_PO4.res'),'ascii');
REF_sed_PO4_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_PO4.res'),'ascii');
REF_sed_PO4_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_PO4.res'),'ascii');
REF_sed_PO4_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_ocn_PO4.res'),'ascii');
REF_sed_PO4_exp5 = load(fullfile(exp_5,'/biogem/biogem_series_ocn_PO4.res'),'ascii');
REF_sed_PO4_exp6 = load(fullfile(exp_6,'/biogem/biogem_series_ocn_PO4.res'),'ascii');

REF_sed_ALK_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_ALK.res'),'ascii');
REF_sed_ALK_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_ALK.res'),'ascii');
REF_sed_ALK_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_ALK.res'),'ascii');
REF_sed_ALK_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_ocn_ALK.res'),'ascii');
REF_sed_ALK_exp5 = load(fullfile(exp_5,'/biogem/biogem_series_ocn_ALK.res'),'ascii');
REF_sed_ALK_exp6 = load(fullfile(exp_6,'/biogem/biogem_series_ocn_ALK.res'),'ascii');

REF_sed_DIC_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_ocn_DIC.res'),'ascii');
REF_sed_DIC_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_ocn_DIC.res'),'ascii');
REF_sed_DIC_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_ocn_DIC.res'),'ascii');
REF_sed_DIC_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_ocn_DIC.res'),'ascii');
REF_sed_DIC_exp5 = load(fullfile(exp_5,'/biogem/biogem_series_ocn_DIC.res'),'ascii');
REF_sed_DIC_exp6 = load(fullfile(exp_6,'/biogem/biogem_series_ocn_DIC.res'),'ascii');

set(0,'defaultLineLineWidth', 2)
set(0,'DefaultAxesFontSize',10)

figure
grid on
hold on

subplot(3,2,1)
if(plot_mean) % mean (mol/kg)
    plot(REF_sed_O2_exp1(:,1),REF_sed_O2_exp1(:,3)*1e+6,'b', ...
        REF_sed_O2_exp2(:,1),REF_sed_O2_exp2(:,3)*1e+6,'r--', ...
        REF_sed_O2_exp3(:,1),REF_sed_O2_exp3(:,3)*1e+6,'g--', ...
        REF_sed_O2_exp4(:,1),REF_sed_O2_exp4(:,3)*1e+6,'k--', ... % k:
        REF_sed_O2_exp5(:,1),REF_sed_O2_exp5(:,3)*1e+6,'m:', ...  % '---m:'
        REF_sed_O2_exp6(:,1),REF_sed_O2_exp6(:,3)*1e+6,'c:'); 
    ylabel('O_2 (\mumol kg^{-1})');
else
    % total (mol)
    plot(REF_sed_O2_exp1(:,1),REF_sed_O2_exp1(:,2),'b', ...
         REF_sed_O2_exp2(:,1),REF_sed_O2_exp2(:,2),'r--', ...
         REF_sed_O2_exp3(:,1),REF_sed_O2_exp3(:,2),'g--', ...
         REF_sed_O2_exp4(:,1),REF_sed_O2_exp4(:,2),'k--', ...
         REF_sed_O2_exp5(:,1),REF_sed_O2_exp5(:,2),'m:', ...
         REF_sed_O2_exp6(:,1),REF_sed_O2_exp6(:,2),'c:'); 
    ylabel('O_2 (mol)');
end
xlabel('yrs ');
%hleg=legend('Abiotic - No OMEN', 'Abiotic - with OMEN', 'Biotic - all remin. k=0.1', 'Biotic - smaller k'); 
%ylim([0 300])
%set(hleg,'FontSize',4);
%set(hleg,'Location','best')


subplot(3,2,2)
if(plot_mean)  % mean (mol/kg)
	plot(REF_sed_SO4_exp1(:,1),REF_sed_SO4_exp1(:,3)*1e+6,'b',REF_sed_SO4_exp2(:,1),REF_sed_SO4_exp2(:,3)*1e+6,'r--',REF_sed_SO4_exp3(:,1),REF_sed_SO4_exp3(:,3)*1e+6,'g--',REF_sed_SO4_exp4(:,1),REF_sed_SO4_exp4(:,3)*1e+6,'k--', REF_sed_SO4_exp5(:,1),REF_sed_SO4_exp5(:,3)*1e+6,'m:' , REF_sed_SO4_exp6(:,1),REF_sed_SO4_exp6(:,3)*1e+6,'c:'); 
    ylabel('SO_4 (\mumol kg^{-1})');
else    % total (mol)
    plot(REF_sed_SO4_exp1(:,1),REF_sed_SO4_exp1(:,2),'b',REF_sed_SO4_exp2(:,1),REF_sed_SO4_exp2(:,2),'r--',REF_sed_SO4_exp3(:,1),REF_sed_SO4_exp3(:,2),'g--',REF_sed_SO4_exp4(:,1),REF_sed_SO4_exp4(:,2),'k--', REF_sed_SO4_exp5(:,1),REF_sed_SO4_exp5(:,2),'m:', REF_sed_SO4_exp6(:,1),REF_sed_SO4_exp6(:,2),'c:' ); 
    ylabel('SO_4 (mol)');
end
xlabel('yrs ');
hleg=legend('Archer default', 'Invariant - kord 2', 'Invariant - kord 3', 'Invariant - kord 4', 'Invariant - kord 5', 'Invariant - kord 10'); 
set(hleg,'FontSize',6);
set(hleg,'Location','SouthEast');

subplot(3,2,3)
if(plot_mean) % mean (mol/kg)
	plot(REF_sed_H2S_exp1(:,1),REF_sed_H2S_exp1(:,3)*1e+6,'b',REF_sed_H2S_exp2(:,1),REF_sed_H2S_exp2(:,3)*1e+6,'r--', REF_sed_H2S_exp3(:,1),REF_sed_H2S_exp3(:,3)*1e+6,'g--',REF_sed_H2S_exp4(:,1),REF_sed_H2S_exp4(:,3)*1e+6,'k--', REF_sed_H2S_exp5(:,1),REF_sed_H2S_exp5(:,3)*1e+6,'m:', REF_sed_H2S_exp6(:,1),REF_sed_H2S_exp6(:,3)*1e+6,'c:' );
    ylabel('H_2S (\mumol kg^{-1})');
else % total (mol)
    plot(REF_sed_H2S_exp1(:,1),REF_sed_H2S_exp1(:,2),'b',REF_sed_H2S_exp2(:,1),REF_sed_H2S_exp2(:,2),'r--', REF_sed_H2S_exp3(:,1),REF_sed_H2S_exp3(:,2),'g--',REF_sed_H2S_exp4(:,1),REF_sed_H2S_exp4(:,2),'k--', REF_sed_H2S_exp5(:,1),REF_sed_H2S_exp5(:,2),'m:', REF_sed_H2S_exp6(:,1),REF_sed_H2S_exp6(:,2),'c:' );
    ylabel('H_2S (mol)');
end
xlabel('yrs ');
%hleg=legend('small k', 'higher k'); 
%set(hleg,'Location','SouthEast')

subplot(3,2,4)
if(plot_mean) % mean (mol/kg)
    plot(REF_sed_PO4_exp1(:,1),REF_sed_PO4_exp1(:,3)*1e+6,'b',REF_sed_PO4_exp2(:,1),REF_sed_PO4_exp2(:,3)*1e+6,'r--', REF_sed_PO4_exp3(:,1),REF_sed_PO4_exp3(:,3)*1e+6,'g--',REF_sed_PO4_exp4(:,1),REF_sed_PO4_exp4(:,3)*1e+6,'k--', REF_sed_PO4_exp5(:,1),REF_sed_PO4_exp5(:,3)*1e+6,'m:', REF_sed_PO4_exp6(:,1),REF_sed_PO4_exp6(:,3)*1e+6,'c:' ); 
    ylabel('PO_4 (\mumol kg^{-1})');
else % total (mol)
    plot(REF_sed_PO4_exp1(:,1),REF_sed_PO4_exp1(:,2),'b',REF_sed_PO4_exp2(:,1),REF_sed_PO4_exp2(:,2),'r--', REF_sed_PO4_exp3(:,1),REF_sed_PO4_exp3(:,2),'g--', REF_sed_PO4_exp3(:,1),REF_sed_PO4_exp3(:,2),'g--',REF_sed_PO4_exp4(:,1),REF_sed_PO4_exp4(:,2),'k--', REF_sed_PO4_exp5(:,1),REF_sed_PO4_exp5(:,2),'m:', REF_sed_PO4_exp6(:,1),REF_sed_PO4_exp6(:,2),'c:' );
    ylabel('PO_4 (mol)');
end
xlabel('yrs ');
%ylim([2.15 2.17])
% hleg=legend('small k', 'higher k', 'SPIN no Corg'); 
% set(hleg,'Location','NorthEast')

subplot(3,2,5)
if(plot_mean)  % mean (mol/kg)
    plot(REF_sed_ALK_exp1(:,1),REF_sed_ALK_exp1(:,3)*1e+6,'b',REF_sed_ALK_exp2(:,1),REF_sed_ALK_exp2(:,3)*1e+6,'r--', REF_sed_ALK_exp3(:,1),REF_sed_ALK_exp3(:,3)*1e+6,'g--',REF_sed_ALK_exp4(:,1),REF_sed_ALK_exp4(:,3)*1e+6,'k--', REF_sed_ALK_exp5(:,1),REF_sed_ALK_exp5(:,3)*1e+6,'m:', REF_sed_ALK_exp6(:,1),REF_sed_ALK_exp6(:,3)*1e+6,'c:'); 
    ylabel('ALK (\mumol kg^{-1})');
else
	plot(REF_sed_ALK_exp1(:,1),REF_sed_ALK_exp1(:,2),'b',REF_sed_ALK_exp2(:,1),REF_sed_ALK_exp2(:,2),'r--', REF_sed_ALK_exp3(:,1),REF_sed_ALK_exp3(:,2),'g--',REF_sed_ALK_exp4(:,1),REF_sed_ALK_exp4(:,2),'k--', REF_sed_ALK_exp5(:,1),REF_sed_ALK_exp5(:,2),'m:', REF_sed_ALK_exp6(:,1),REF_sed_ALK_exp6(:,2),'c:');
    ylabel('ALK (mol)');
end
xlabel('yrs ');
%hleg=legend('small k', 'higher k'); 
%set(hleg,'FontSize',10)
%set(hleg,'Location','SouthEast')
%title('oceanic Alkalinity','FontSize',18);

subplot(3,2,6)
if(plot_mean)  % mean (mol/kg)
    plot(REF_sed_DIC_exp1(:,1),REF_sed_DIC_exp1(:,3)*1e+6,'b',REF_sed_DIC_exp2(:,1),REF_sed_DIC_exp2(:,3)*1e+6,'r--', REF_sed_DIC_exp3(:,1),REF_sed_DIC_exp3(:,3)*1e+6,'g--',REF_sed_DIC_exp4(:,1),REF_sed_DIC_exp4(:,3)*1e+6,'k--', REF_sed_DIC_exp5(:,1),REF_sed_DIC_exp5(:,3)*1e+6,'m:', REF_sed_DIC_exp6(:,1),REF_sed_DIC_exp6(:,3)*1e+6,'c:'); 
    ylabel('DIC (\mumol kg^{-1})');
else
	plot(REF_sed_DIC_exp1(:,1),REF_sed_DIC_exp1(:,2),'b',REF_sed_DIC_exp2(:,1),REF_sed_DIC_exp2(:,2),'r--', REF_sed_DIC_exp3(:,1),REF_sed_DIC_exp3(:,2),'g--',REF_sed_DIC_exp4(:,1),REF_sed_DIC_exp4(:,2),'k--', REF_sed_DIC_exp5(:,1),REF_sed_DIC_exp5(:,2),'m:', REF_sed_DIC_exp6(:,1),REF_sed_DIC_exp6(:,2),'c:');
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
    print('-depsc', 'cgenie_output/000_FOR_GMD_V1507_1707_20kyr/TIMESERIES/1507_Boudreau_2_mean');
else
    print('-depsc', 'cgenie_output/000_FOR_GMD_V1507_1707_20kyr/TIMESERIES/1707_invariant_007');
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% atm O2 and CO2
% %%%% load other data
% REF experiments
if(false)

REF_sed_pCO2_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
REF_sed_pCO2_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
REF_sed_pCO2_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
REF_sed_pCO2_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
REF_sed_pCO2_exp5 = load(fullfile(exp_5,'/biogem/biogem_series_atm_pCO2.res'),'ascii');
% REF_sed_pCO2_smallerk_DISS_HACK = load(fullfile(exp_6,'/biogem/biogem_series_atm_pCO2.res'),'ascii');

REF_sed_pO2_exp1 = load(fullfile(exp_1,'/biogem/biogem_series_atm_pO2.res'),'ascii');
REF_sed_pO2_exp2 = load(fullfile(exp_2,'/biogem/biogem_series_atm_pO2.res'),'ascii');
REF_sed_pO2_exp3 = load(fullfile(exp_3,'/biogem/biogem_series_atm_pO2.res'),'ascii');
REF_sed_pO2_exp4 = load(fullfile(exp_4,'/biogem/biogem_series_atm_pO2.res'),'ascii');
REF_sed_pO2_exp5 = load(fullfile(exp_1,'/biogem/biogem_series_atm_pO2.res'),'ascii');
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
plot(REF_sed_pCO2_exp1(:,1),REF_sed_pCO2_exp1(:,3)*1e+6,'k', REF_sed_pCO2_exp2(:,1),REF_sed_pCO2_exp2(:,3)*1e+6,'g--', REF_sed_pCO2_exp3(:,1),REF_sed_pCO2_exp3(:,3)*1e+6,'b-.',REF_sed_pCO2_exp4(:,1),REF_sed_pCO2_exp4(:,3)*1e+6,'r-.',REF_sed_pCO2_exp5(:,1),REF_sed_pCO2_exp5(:,3)*1e+6,'m:'); 
%plot(REF_sed_pCO2_exp1(:,1),REF_sed_pCO2_exp1(:,3)*1e+6,'b',REF_sed_pCO2_smallerk(:,1),REF_sed_pCO2_smallerk(:,3)*1e+6,'r',REF_sed_pCO2_NoOMEN(:,1),REF_sed_pCO2_NoOMEN(:,3)*1e+6,'g--',REF_sed_pCO2_SPIN(:,1),REF_sed_pCO2_SPIN(:,3)*1e+6,'k--'); 
xlabel('time (years) ');
ylabel('global pCO2 (ppm)');
%x = [0,5000,10000,15000,20000,25000];
%set(gca,'XTick',x,'XTickLabel',sprintf('%5.0f|',x))
hleg=legend('Abiotic - with OMEN', 'Biotic - k=0.1', 'Biotic - k=1.0', 'Biotic - k=4.0', 'Biotic - k1=0.01, k2=0.001'); 
%title('global pCO2 (ppm)','FontSize',18);

subplot(2,1,2)
%plot(REF(:,1),REF(:,3)*1e+6,'k')
%plot(REF_sed_pO2(:,1),REF_sed_pO2(:,3),'b',REF_sed_pO2_smallerk(:,1),REF_sed_pO2_smallerk(:,3),'r' );
plot(REF_sed_pO2_exp1(:,1),REF_sed_pO2_exp1(:,3),'k', REF_sed_pO2_exp2(:,1),REF_sed_pO2_exp2(:,3),'g--', REF_sed_pO2_exp3(:,1),REF_sed_pO2_exp3(:,3),'b-.',REF_sed_pO2_exp4(:,1),REF_sed_pO2_exp4(:,3),'r-.',REF_sed_pO2_exp5(:,1),REF_sed_pO2_exp5(:,3),'m:'); 
%plot(REF_sed_pO2_exp1(:,1),REF_sed_pO2_exp1(:,3),'b',REF_sed_pO2_smallerk(:,1),REF_sed_pO2_smallerk(:,3),'r',REF_sed_pO2_NoOMEN(:,1),REF_sed_pO2_NoOMEN(:,3),'g--',REF_sed_pO2_SPIN(:,1),REF_sed_pO2_SPIN(:,3),'k--'); 
xlabel('time (years) ');
ylabel('global pO2 (atm)');
ylim([0.20 0.21])
%ylim([0.207 0.2095])
%x = [0,5000,10000,15000,20000,25000];
%set(gca,'XTick',x,'XTickLabel',sprintf('%5.0f|',x))
% hleg=legend('small k', 'higher k', 'SPIN no Corg'); 
% set(hleg,'Location','NorthEast')

%print('-depsc', fullfile(exp_1,'/1_OUTPUT_PLOTS/0_Atmp-time-series'));
print('-depsc', 'cgenie_output/0_PLOTS/plots_1602/1_Atmp-time-series_NO_2redox_diff_k-values_1diff-coeff_1602');

end