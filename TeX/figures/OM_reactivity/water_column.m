clear all 
close all

%__________________________________________________________________________

depth_max = 400;
depth=[100:10:depth_max];
f1=0.945*exp((100-depth)/589);

f2=0.055;

set(0,'defaultLineLineWidth', 4)
set(0,'DefaultAxesFontSize',16)

figure      
hold on
box on

%subplot(2,1,1)
plot(f1, -depth, 'r-', [f2, f2], [-100, -depth_max], 'g-')
xlabel('Fraction of OM (ocean)');
ylabel('Depth (m)');
print('-depsc', 'OM_degradation_400m_ocean');

% calculate fraction at the top of the sediments
sed_f1 = f1(end)/(f1(end)+f2)
sed_f2 = 1- sed_f1;

%subplot(2,1,2)
figure      
hold on
box on

plot([sed_f1, sed_f1], [0, -100], 'r-', [sed_f2, sed_f2], [0, -100], 'g-')
xlabel('Fraction of OM  (sediments)');
ylabel('Depth (cm)');

print('-depsc', 'OM_degradation_400m_sediments');

