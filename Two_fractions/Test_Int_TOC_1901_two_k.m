format longEng
clear all;

w = 0.1;
C0 = 0.1;
k(1)=0.01;
k(2)=0.1;
DC1 = 10;
DC2 = 0.0;

zbio = 10000;
zinf = 10000;

for i=1:2
a1=(w-sqrt(w^2+4*DC1*k))/(2*DC1);
b1=(w+sqrt(w^2+4*DC1*k))/(2*DC1);
A1=-(C0*b1*exp(b1*zbio))/(a1*exp(a1*zbio)-b1*exp(b1*zbio));

a2= -k/w;

Int_bio =k*( (A1*(exp(a1*zinf)/a1 - exp(b1*zinf)/b1)) + C0/b1 * exp(b1*zinf) - (A1*(1/a1 - 1/b1)+ C0/b1))

Int_nonbio = k*((C0/a2*exp(a2*zinf)) - C0/a2)
end

for i=1:1:10001
    z = (i-1)/10;
    C(i) = k*(A1*exp(a1*z)+(C0-A1)*exp(b1*z));
    C2(i) = k*(C0*exp(a2*z));
end

figure
plot(C)
hold on
plot(C2, 'r')
legend('bio', 'non-bio')
hold off