clear all 
close all

%__________________________________________________________________________
% RCM parameters that describe OM distribution
% high nu, low a high reactivity, low nu, high a low reactivity
%__________________________________________________________________________
nu=0.125; % deep sea 10.0; fresh 0.125
a=3e-4; % fresh 3e-4, unreactive 1e4
nG=200;  % define size of multiG model- e.g. 200G
emin=-15;% reactivity range lower boundary kmin=10^emin
emax=-log10(a)+2;% reactivity range upper boundary kmin=10^emax


%__________________________________________________________________________
%site specific environmental parameters
%__________________________________________________________________________

POC0=0.3437; 
w0=4.1141e-04;  
db0=0.5900; % bioturb
beta=0.5e-5;
por0=0.45;
zbio=10; 

%__________________________________________________________________________
%calculate ki and Fi for nG-G model
%__________________________________________________________________________

k(1)= 10^(emin);
kk(1)=10^(emin);
F(1) = gammainc(a*10^emin,nu,'lower');
kk(nG)=10^(emax);
k(nG)=10^(emax);
F(nG) = gammainc(a*10^emax,nu,'upper');
%k(nG+1)= 10^emax;
%F(nG+1)= gammainc(a*10^10000,nu,'upper')-gammainc(a*10^emax,nu,'upper');


                for i=2:nG-1
                    ne=emin+(i-1)*(emax-emin)/(nG-1);
                    kk(i)=10^ne;
                    G_inc_0 = gammainc(a*kk(i-1),nu,'upper');
                    G_inc_1 = gammainc(a*kk(i),nu,'upper');
                    F(i) = (G_inc_0 - G_inc_1);
                    k(i)=kk(i-1)+(kk(i)-kk(i-1))/2;
                end 
                
                
figure                
plot(log10(k), F, '.')

aaa=sum(F) %check sum


%__________________________________________________________________________
%Check analytical solution vs rcm solution t-dependent
%__________________________________________________________________________

for tt=1:1000
    t=(tt-1)*10;
 for i=1:nG
    POC_t(tt,i)=F(i)*POC0*exp(-k(i)*t);
 end
 POC_rcm(tt)=POC0*(a/(a+t))^nu;
 sPOC_t(tt)=sum(POC_t(tt,:));
 diff_POC(tt)=(POC_rcm(tt)-sPOC_t(tt))/POC_rcm(tt);
end

figure
subplot(2,2,1)
plot(sPOC_t)
hold on
plot(POC_rcm, 'r')
hold off

%subplot(2,2,2)
%plot(diff_POC)

%__________________________________________________________________________
%Calculate example depth profiles
%__________________________________________________________________________

for x=1:1000
    z=(x-1)/100;
    for i=1:nG
         aa = (w0-(w0^2+4*db0*k(i))^(1/2))/(2*db0);
         bb = (w0+(w0^2+4*db0*k(i))^(1/2))/(2*db0);
         
         
         B = (exp(aa*zbio)*POC0*aa*F(i))/(-exp(aa*zbio)*aa+exp(bb*zbio)*bb);
         A = (F(i)*POC0)-B;
         
        POC(x,i)=A*exp(aa*z) + B*exp(bb*z);
        FF(x,i)=POC(x,i)/POC0;
          
    end
    
% 1G 
%__________________________________________________________________________
%      kk=nu/a;
%      aa = (w0-(w0^2+4*db0*kk)^(1/2))/(2*db0);
%      bb = (w0+(w0^2+4*db0*kk)^(1/2))/(2*db0);
%           
%      AA = -(POC0*exp(bb*z)*bb)/(exp(aa*z)*aa-exp(bb*z)*bb);
%      BB = (exp(aa*z)*POC0*aa)/(exp(aa*z)*aa-exp(bb*z)*bb);
%      
%      POC_1G(x)=AA*exp(aa*z) + BB*exp(bb*z);


% % RCM
% %__________________________________________________________________________
% 
%     depth=z;
%     agercm=depth/w0;%(depth+1/beta*por0*(exp(-beta*depth)-1))/(w0*(1-por0));
%     POCrcm(x)=POC0*(a/(a+agercm))^nu;
%     
POC(isnan(POC))=0;
sPOC(x)=sum(POC(x,:));
%     
% %test age  
% 
% agetest(x)=-a*(exp(log(sPOC(x)/POC0)/nu)-1)/exp(log(sPOC(x)/POC0)/nu);
% agelinear(x)=z/w0;
% POCrcmtest(x)=POC0*(a/(a+agetest(x)))^nu;
    
end



figure
plot(POC(1:1000,1:nG), -(1:1000))

figure
plot(sPOC(1:1000),-(1:1000), 'k-')
hold on
% plot(POCrcm,-(1:1000), 'r-')
% plot(POCrcmtest,-(1:1000), 'g-')
hold off

set(0,'defaultLineLineWidth', 4)
set(0,'DefaultAxesFontSize',16)

figure
plot(log10(k), F, 'k')
hold on
plot(log10(k), FF(2,:), 'k--')
%plot(log10(k), FF(6,:), 'k--')
plot(log10(k), FF(100,:), 'k--')
plot(log10(k), FF(1000,:), 'k--')
hleg=legend('t = 0 years','t < 1 year','t ~ 500 years','t ~ 5000 years')
set(hleg,'Location','NorthWest');
hold off

print('-depsc', 'OM_degradation_rate_distribution_400m_0908');

% figure
% plot(agetest)
% hold on
% %plot(agelinear, 'r')
% hold off






