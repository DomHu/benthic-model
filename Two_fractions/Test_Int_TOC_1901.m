format longEng
clear all;
close all

w = 0.1;
C0 = 0.1;
k=0.05;
DC1 = 0.10;
DC2 = 0.0;
por = 0.85;

zbio = 100;
zinf = 100;

a1=(w-sqrt(w^2+4*DC1*k))/(2*DC1);
b1=(w+sqrt(w^2+4*DC1*k))/(2*DC1);
A1=-(C0*b1*exp(b1*zbio))/(a1*exp(a1*zbio)-b1*exp(b1*zbio));% for C(inf) = 0: C0;

a2= -k/w;

Int_bio = k*( (A1*(exp(a1*zinf)/a1 - exp(b1*zinf)/b1)) + C0/b1 * exp(b1*zinf) - (A1*(1/a1 - 1/b1)+ C0/b1)) % for C(inf) = 0:  k* (A1/a1*exp(a1*zinf)-A1/a1)%

Int_nonbio = k*((C0/a2*exp(a2*zinf)) - C0/a2)

for i=1:1:1001
    z = (i-1)/10;
    C(i) = k*(A1*exp(a1*z));
    C2(i) = k*(C0*exp(a2*z));
end

figure
plot(C)
hold on
plot(C2, 'r')
legend('bio', 'non-bio')
hold off

for i=1:100000
    DC1v(i)=1/i;
    a2v(i)=a2;
a1v(i)=(w-sqrt(w^2+4*DC1v(i)*k))/(2*DC1v(i));
end

Fnonbio =C0*(1-por)*w;
% for C(inf) = 0: C0bio= Fnonbio/((1-por)*(w-DC1*a1));
%C0bio = (Fnonbio + (1-por)*DC1*A1/a1 - (1-por)*DC1*A1/b1)*1/((1-por)*w-(1-por)*DC1/b1);
N1=b1*exp(b1*zbio)/(a1*exp(a1*zbio)-b1*exp(b1*zbio));
N2=DC1/(w-DC1/b1);
N3=(1/a1-1/b1);
C0bio = Fnonbio/((1-por)*(w-DC1/b1)) * 1/(1+N1*N2*N3);
%C0bio = Fnonbio/((1-por)*(w-DC1/b1)) * 1/(1+((b1*exp(b1*zbio)/(a1*exp(a1*zbio)-b1*exp(b1*zbio))*DC1/(w-DC1/b1)*(1/a1-1/b1))));

%c1=a1;
%a1=b1;
%b1=c1;

A1=-(C0bio*b1*exp(b1*zbio))/(a1*exp(a1*zbio)-b1*exp(b1*zbio)); % for C(inf) = 0: C0bio;
Int_bio_new = k*( (A1*(exp(a1*zinf)/a1 - exp(b1*zinf)/b1)) + C0bio/b1 * exp(b1*zinf) - (A1*(1/a1 - 1/b1)+ C0bio/b1)) % for C(inf) = 0:   k* (A1/a1*exp(a1*zinf)-A1/a1)

for i=1:1:1001
    z = (i-1)/10;
    C(i) = k*(A1*exp(a1*z));
    %C(i) = k*(A1*(exp(a1*z)-exp(b1*z))+C0bio*exp(b1*z))
    C2(i) = k*(C0bio*exp(a2*z));
end
figure
plot(C)
hold on
plot(C2, 'r')
legend('bio new', 'non-bio')
hold off
