format longEng
clear all;
%close all

tol_const = 1e-18;

w = 0.2;
C0 = 2.166666667e-3;
k=0.05;
DC1 = 0.10;
DC2 = 0.0;
por = 0.85;

zbio = 10.0;
zinf = 100;

a1=(w-sqrt(w^2+4*DC1*k))/(2*DC1);
b1=(w+sqrt(w^2+4*DC1*k))/(2*DC1);
A1=-(C0*b1*exp(b1*zbio))/(a1*exp(a1*zbio)-b1*exp(b1*zbio));% for C(inf) = 0: C0;

a2= -k/w;

% Int_bio_OLD = k*( (A1*(exp(a1*zinf)/a1 - exp(b1*zinf)/b1)) + C0/b1 * exp(b1*zinf) - (A1*(1/a1 - 1/b1)+ C0/b1)) % for C(inf) = 0:  k* (A1/a1*exp(a1*zinf)-A1/a1)%


Fnonbio =C0*(1-por)*w

% % for C(inf) = 0: C0bio= Fnonbio/((1-por)*(w-DC1*a1));
% %C0bio = (Fnonbio + (1-por)*DC1*A1/a1 - (1-por)*DC1*A1/b1)*1/((1-por)*w-(1-por)*DC1/b1);

%  % MY CALCULATION
% N1=b1*exp(b1*zbio)/(a1*exp(a1*zbio)-b1*exp(b1*zbio));
% N2=DC1/(w-DC1/b1);
% N3=(1/a1-1/b1);
% C0bio = Fnonbio/((1-por)*(w-DC1/b1)) * 1/(1+N1*N2*N3); % MY CALCULATION
% 
% A1=-(C0bio*b1*exp(b1*zbio))/(a1*exp(a1*zbio)-b1*exp(b1*zbio)); % for C(inf) = 0: C0bio;
% Int_bio_new_MINE = k*( (A1*(exp(a1*zinf)/a1 - exp(b1*zinf)/b1)) + C0bio/b1 * exp(b1*zinf) - (A1*(1/a1 - 1/b1)+ C0bio/b1)) % for C(inf) = 0:   k* (A1/a1*exp(a1*zinf)-A1/a1)

 % MAPLE CALCULATION
C0bio_M = (Fnonbio*(-a1*exp(a1*zbio)+b1*exp(b1*zbio)))/(-DC1*b1*a1*exp(b1*zbio) + DC1*b1*a1*exp(a1*zbio) + ...
        DC1*b1*a1*por*exp(b1*zbio) - DC1*b1*a1*por*exp(a1*zbio) - w*a1*exp(a1*zbio) + w*b1*exp(b1*zbio) + ...
        w*por*a1*exp(a1*zbio) - w*por*b1*exp(b1*zbio))
%C0bio_M = C0;    

A1_M=-(C0bio_M*b1*exp(b1*zbio))/(a1*exp(a1*zbio)-b1*exp(b1*zbio)) % for C(inf) = 0: C0bio;
%A2=C0; % when zbio = 0, whole column is non-bioturbated
A2 = (A1_M*(exp(a1*zbio)-exp(b1*zbio))+C0bio_M*exp(b1*zbio))/exp(a2*zbio) % when zbio > 0

% for zbio = 0, all column non-bioturbated Int_nonbio = k*((A2/a2*exp(a2*zinf)) - A2/a2)
Int_nonbio = k*((A2/a2*exp(a2*zinf)) - A2/a2*exp(a2*zbio))

% First two terms are veery large +/-
%for all column bioturbated Int_bio_Maple_ALL = k*( (A1_M*(exp(a1*zinf)/a1 - exp(b1*zinf)/b1)) + C0bio_M/b1 * exp(b1*zinf) - (A1_M*(1/a1 - 1/b1)+ C0bio_M/b1)) % for C(inf) = 0:   k* (A1/a1*exp(a1*zinf)-A1/a1)
Int_bio_Maple = k*( (A1_M*(exp(a1*zbio)/a1 - exp(b1*zbio)/b1)) + C0bio_M/b1 * exp(b1*zbio) - (A1_M*(1/a1 - 1/b1)+ C0bio_M/b1)) % for C(inf) = 0:   k* (A1/a1*exp(a1*zinf)-A1/a1)
first_term = (A1_M*(exp(a1*zbio)/a1 - exp(b1*zbio)/b1))
second_term = C0bio_M/b1 * exp(b1*zbio)

if(first_term == -second_term)
    equal = true
else
    equal = false
end
difference = first_term + second_term

Int_bio_Maple_last_term = k*( - (A1_M*(1/a1 - 1/b1)+ C0bio_M/b1)) % for C(inf) = 0:   k* (A1/a1*exp(a1*zinf)-A1/a1)

Sum_bio_nonbio_ALL = Int_nonbio + Int_bio_Maple
Sum_bio_nonbio_laster_term = Int_nonbio + Int_bio_Maple_last_term


for i=1:1:1001
    z(i) = (i-1)/10;
    if(z<=zbio)
        C1(i)=A1_M*(exp(a1*z(i))-exp(b1*z(i)))+C0bio_M*exp(b1*z(i));
    else                   
        C1(i)=A2*exp(a2*z(i));
    end
%   for bioturbated and unbioturbated separately
%     C(i) = k*(A1_M*exp(a1*z));
%     %C(i) = k*(A1*(exp(a1*z)-exp(b1*z))+C0bio*exp(b1*z))
%     C2(i) = k*(C0bio_M*exp(a2*z));
    
end
figure
plot(C1,-z)
 hold on
% plot(C2, 'r')
% legend('bio new', 'non-bio')
t=xlim;
plot([0,t(1,2)], [-zbio,-zbio], 'k--') 
xlabel ('[TOC] ')
ylabel('Depth (cm)')
 hold off
