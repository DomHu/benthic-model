# _________________________________________________________________________

# AUTOMATIC CODE GENERATOR (ACG) FOR CONSTRUCTING USER-DEFINED
# BIOGEOCHEMICAL REACTION NETWORKS
# 
# version 1.2
# COPYRIGHT (c) 2001 P.A.G. Regnier 
# All Rights Reserved
# 
# Research Unit on Biogeochemical Systems Dynamics
# Department of Geochemistry, Utrecht University, 
# The Netherlands
# __________________________________________________________________________
# __________________________________________________________________________
# INPUT TYPES
# 
# OOO : Sections that should be modified by the user
# OOO : Sections that should NOT be modified by the user

# OOO : comments 
# OOO : Maple input
# OOO : Maple output (appears only after you have executed the spreadsheet) 
# OOO  : Maple input entries that have to be specified by the user 
# 
# WWW  : Hyperlink to the Knowledge Book
# 
# _________________________________________________________________________ 

#  Maple specific info
restart ;
#with(Spread) :
 precision := double :
#  Summary and governing equations
#  Governing equation 
# The governing equation solved is of the form :
# theta(x)*dC/dt = -d(Fdiff+Fadv)/dx + R,
# where 
# - dX/dp is the partial derivative of X with respect to p [M/L^3/[p]]
# - theta(x) = A(x)*por(x) for dissolved, and A(x)*(1-por(x)) for solid species [L^2]
# - Fdiff = -(D(x)*theta(x)*dC/dx) [M/L^2/T]
# - Fadv  = (v*theta(x)*C)[M/L^2/T]
# - A is the cross section area [L^2], por is porosity [-]
# - v is the advection velocity [L/T], the sum of 
# the water flow velocity (vwat = water flow q in [L^3/T]/theta(x), solutes only) and 
# an advective velocity w acting upon both solids and solutes (e.g. movement due to fixed reference frame)
# - D is the effective diffusion coefficient [L^2/T], the sum of 
# the molecular diffusion coefficient (Dmol, solutes only), 
# the bioturbation coefficient (Db) and 
# the dispersion coefficient (Disp = aL*|vwat|, solutes only, aL in [L]). 
# - R is the sum of the reaction terms [M/L^3/T] as specified in the reaction network established below
#  Caveats 
# - for real numbers you should add a point after the number.
# - make sure all units match. 
# e.g. rates and dXdt, or flux boundary conditions and concentrations
# for the conversion of solid to solute units, 
# one may define temporary variables that can be used in the rate laws below:
# s_dens := 2.5; # solid density in [g/cm_solid^3]
# sd := 1000. * s_dens * (1. - por(j)) / por(j); 
# the factor 1.d03 converts cm^3 to liter, e.g. [g/cm^3] to [g/l]. 
# note that you need to refer to porosity exactly as por(j).
#  Physics - Parameters
# The list of length nparphys & nparphys2 is given by phys_name & phys_name2; the values collected in phys_val & phys_val2.
#  Spatial and temporal domain size
# tot_time: total length of simulation (T, e.g. years)
# depth_max: spatial extent, e.g. total depth of simulation (L, e.g. cm)
tot_time := 500. ; # Dom was 2000.
depth_max := 100.0 ;


#  Temperature and salinity
# T_C: temperature (Celsius)
# S : salinity (PSU)
# T and S are used to calculate the molecular diffusion coefficients
# T: absolute temperature (Kelvin)
T_C := 8.1; #Dom was 10.014.3;#siteK:15.0; #siteA:16.8;
T := T_C + 273.15 ;
S := 35.00; # Dom was 35.00 Jens 38
Depth:= 500.0; #77.0; #site1: 18.3;
P := Depth*9.81*1027./100000.;
R := 83.145;  
#  Transport coefficients
# vaL: longitudinal dispersivity (L). The dispersion coefficient is calculated as Disp = aL*|vwat| 
# viq: depth dependency of flow. 0: constant, else changing with depth (i.e. distance)
# vq0: water flow (either constant or else at x=0, L^3/T)
# viw: depth dependency of w. 0: constant, else changing with depth (i.e. distance)
# vw0: advection velocity working on solids and solutes (either constant or else at x=0, L/T)
# viDb: depth dependency of Db. 0: constant, else changing with depth (i.e. distance)
# vDb0: bioturbation coefficient working on solids and solutes (either constant or else at x=0, L^2/T)
# 
# - user specified profiles of q(x), w(x) or Db(x) must be defined in the fortran subroutines advcoeff.f and diffcoeff.f, respectively.
# - the molecular diffusion coefficients for the species used are specified further below
val := 0.0 ;
viq := 0 ;
vq0 := 0. ;
viw := 0 ;
vw0 := 267.0e-3;#3.3*10.0^(-0.87478367-0.00043512*Depth); Jens 30.0
viDb := 1;
vDb0 := 5.2*(10.0^(0.7624-0.0003972*Depth)); # Jens 5.0
#  Porosity profile and cross section area
# vipor: depth dependency of porosity. 0: constant, else changing with depth (i.e. distance)
# vpor0: porosity value (either constant or else at x=0)
# viarea: depth dependency of cross section area. 0: constant, else changing with depth (i.e. distance)
# varea0: cross section area (either constant or else at x=0, L^2)
# - user specified profiles of por(x) or area(x) must be defined in the fortran subroutine porarea.f.
vipor := 1 ;
vpor0 := 0.85; # Dom was 0.85 #Jens 0.73
viarea := 0 ;
varea0 := 1.0 ;
# 
#  Grid and discretization
# Dt: time step of numerical integration (T). 
# note that there are options defined in drivervalues.f that allow automatical selection of timestep
# nnodes: number of nodes in the spatial domain. 
# for a regular grid, the grid spacing of the concentration profile is then depth_max/(nnodes-1)
# vigridv: type of grid. 
# 0 is for regular, evenly spaced grid, 
# else the user needs to specify the grid in the fortran subroutine gridsetup.f
Dt := 0.0005 ;
nnodes := 122;
vigrid := 1;
#  Reaction Network - Size and Variables
#  Size of reaction network
# nsolids : number of solid species 
# ndissolved : number of dissolved species 
# ncompo : total number of species 
# nreactions : total number of reactions (including equilibrium rxns) 
# neqrxns : number of equilibrium reaction
nsolids := 13;
ndissolved := 19;
ncompo := nsolids + ndissolved;
nreactions := 43;
neqrxns := 8 ;
#  List of variables
# variables: list of variables to model. example: 
# listsolids: species number which is a SOLID species. 
# note: all other variables are temporary and are NOT parsed to the ACG
# 
# Example: 
# variables:=[O2, so4, MnOx, FeOx, hco3, co3, hplus, hs];
# listsolids: = [3,4];
variables := [G1, o2, no3, mno2, feoh3, so4, ch4, nh4, po4, mn, fe, h2s, hs, ch4g, h2co3, hco3, co3, boh4, boh3, hplus, caco3, ca, snh4, spo4, fes, feco3, s0, fes2, sfe, sfp, mnco3, G2] ;
listsolids := [1, 4, 5, 21, 23, 24, 25, 26, 28, 29, 30, 31, 32] ;
#  Biogeochemistry - Rate laws
# Definition of kinetic rate laws
# rate.i : array of rates. 
# - For equilibrium rate expression, a kinetic rate MUST be specified as well. It will be overwritten in the equilibrium section below, but you need it as space holder and stoichiometry. Furthermore, the steadystate module uses detailed balancing method with fast kinetics. Therefore, in the example below, kf will have to be defined as a rather large number and kb = kf*Keq
# note: all other variables are temporary and are NOT parsed to the ACG
# - conditional statements: if a rate law depends on a conditional statement you need to make use of the subroutine switches.f. Example: dissolution (Rd) is only to take place at undersaturation, thus  Rd= f(saturation). If saturation > 1, Rd>0, else Rd=0. This canbe implemented as Rd:=k*H1*saturation, where H1 is toggled between 0 and 1. Rather than giving the condition here in maple, for now you need to do this in "switches.f", where you program the conditions for H1, e.g. H1=0, If (A*B/K>1) then H1 = 1.
#  
# example: 
# rate1 := 1000.*O2*hs; # rate law for 2O2 + HS -> SO4 + Hplus
# rate2 := kf*hplus*co3 - kb*hco3; # kinetic rate law for HCO3 = CO3 + Hplus (equilibrium)
# 
#  Primary redox reactions WWW
#primary redox

fo2   :=(sw02+(1.0-sw02)*o2/kmo2);
fno3  :=sw03*(1.0-fo2)*(sw04+(1.0-sw04)*no3/kmno3);
fmno2 :=sw05*(1.0-fo2-fno3)*(sw06+(1.0-sw06)*mno2/kmmno2);
ffeoh3:=sw07*(1.0-fo2-fno3-fmno2)*(sw08+(1.0-sw08)*feoh3/kmfeoh3);
fso4  :=sw09*(1.0-fo2-fno3-fmno2-ffeoh3)*(sw10+(1.0-sw10)*so4/kmso4);
fch4  :=sw11*(1.0-fo2-fno3-fmno2-ffeoh3-fso4);

rate1 := k1*G1;
rate37:= k2*G2;
rate2 := k1*G1*fo2;
rate38:= k2*G2*fo2;
rate3 := k1*G1*fno3;
rate39:= k2*G2*fno3;
rate4 := k1*G1*fmno2;
rate40:= k2*G2*fmno2;
rate5 := k1*G1*ffeoh3;
rate41:= k2*G2*ffeoh3;
rate6 := k1*G1*fso4;
rate42:= k2*G2*fso4;
rate7 := k1*G1*fch4;
rate43:= k2*G2*fch4;


#secondary redox
rate8 := knit*nh4*o2;
rate9 := kmnox*mn*o2;
rate10:= kfemno2*fe*mno2;
rate11:= kfeo2*fe*o2; 
rate12:= kh2so2*(h2s+hs)*o2;
rate13:= kh2smno2*mno2*(hs+h2s);
rate14:= kh2sfeoh3*feoh3*(hs+h2s);
rate15:= kaom*ch4*so4;
rate16:= kaomo2*ch4*o2;
rate17:= kfeso2*fes*o2;

#mineral
rate18:= kmnco3precip*sw12*((mn*co3/KsMnCO3)-1.0);
rate19:= kfesprecip*sw13*(((fe*hs)/(hplus*KsFeS))-1.0);
rate20:= kfesdiss*(1.0-sw13)*fes*(1.0-((fe*hs)/(hplus*KsFeS)));
rate21:= kfeco3precip*sw14*((fe*co3/KsFeCO3)-1.0);
rate22:= kpyr*fes*(h2s+hs); 
rate23:= kfess0*fes*s0;
rate24:= sw15*kcaldiss *caco3* ((1.0 - (ca*(co3)/kspcal)));
rate25:= sw16*kapa*(po4-po4_eq);

#misc
rate26:= (1.0-sw17)*kdis*sw18*(ch4g)*(ch4eq-(ch4));
rate27:= sw17*kgas*(ch4-ch4eq);
rate28:= kdi*sw19*s0;

#eq
rate29:= kf1*hplus*hco3-kb1*h2co3;
rate30:= kf2*hplus*co3-kb2*hco3;
rate31:= kf3*hplus*hs-kb3*h2s ;
rate32:= kf4*hplus*boh4-kb4*boh3;
rate33:= kf5*po4*por(j) - kb5*spo4*(1. - por(j));
rate34:= kf6*nh4*por(j) - kb6*snh4*(1. - por(j));
rate35:= kf7*feoh3*po4- kb7*sfp;
rate36:= kf8*fe*por(j)  - kb8*sfe*(1. - por(j));

#  Biogeochemistry - Stoichiometry

# Stoichiometry of the biogeochemical reactions
# d.sp.dt : rates of change of sp due to the sum of biogeochemical reactions
# note that rateX must be referred to as rX  
#  
# example:
# dO2dt := -2*r1;
# dhco3dt = -r2; 
SD := (1.0 - por(j)) / por(j);
SD1:= 1/por(j);
x  := 106;
y  :=  16; 
z  :=   1;

# 
dG1dt := -r1;
dG2dt := -r37;
do2dt   := -(x+2*y)/x*SD*(r2+r38)-2.0*r8-2.0*r16-2.0*r12-0.25*r11-2.0*SD*r17-0.5*r9;
dno3dt  := -(4.0*x+3.0*y)/(5.0*x)*SD*(r3+r39)+r8;
dmno2dt := -2.0*(r4+r40)+r9/SD-r10-r13;
dfeoh3dt:= -4*(r5+r41)-2.0*r14+1.0*r11/SD+2.0*r10;
dso4dt  := -0.5*SD*(r6+r42)-r15+0.5*r12+SD*r17;#+r28;
dch4dt  :=  0.5*SD*(r7+r43)-r15-r16+r26-r27;
dch4gdt := -r26+r27; 
dnh4dt  := y/x*SD*(r1-r3+r37-r39)-r8-r34*SD1;
dpo4dt  := z/x*SD*(r1+r37)-r25-r33*SD1-r35*SD1;
dmndt   := 2.0*SD*(r4+r40)-r9+r10*SD+r13*SD-r18;
dfedt   := 4.0*SD*(r5+r41)+2.0*SD*r14-1.0*r11+SD*r17-1.0*r19-1.0*r21+SD*r20-r36*SD1-2.0*r10*SD;
dh2sdt  := 0.5*SD*(r6+r42)+r15-r12+r31-1.0*SD*r14-1.0*r19-1.0*SD*r22+3.0*r28+SD*r20-r13*SD;
dhsdt   :=-r31;
dh2co3dt:=((x+y+2.0*z)/x)*SD*(r2+r38)+((x-3.0*y+10.0*z)/(5.0*x))*SD*(r3+r39)-(3.0*x+y-2.0*z)*SD*(r4+r40)/x-((y-2.0*z)/x)*SD*(r6+r42)+((x-2.0*y+4.0*z)/(2.0*x))*SD*(r7+r43)-(7.0*x+y-2.0*z)/x*SD*(r5+r41)+2.0*r8-r15+r16+2.0*r12-4.0*SD*r14+r29+2.0*r11+2.0*r19+1.0*r21+2.0*r28-2.0*SD*r20+2.0*r9+2.0*r10*SD-2.0*r13*SD+r18;
dhco3dt:=-((y+2.0*z)/x)*SD*(r2+r38)+((4.0*x+3.0*y-10.0*z)/(5.0*x))*SD*(r3+r39)+(4.0*x+y-2.0*z)*SD*(r4+r40)/x+((x+y-2.0*z)/x)*SD*(r6+r42)+((y-2.0*z)/x)*SD*(r7+r43)+(8.0*x+y-2.0*z)/x*SD*(r5+r41)-2.0*r8+2.0*r15-2.0*r12+4.0*SD*r14-r29+r30-2.0*r11-2.0*r19-2.0*r21-2.0*r28+2.0*SD*r20-2.0*r9-2.0*r10*SD+2.0*r13*SD-2.0*r18*SD;
dco3dt:=-r30+SD*r24;
dboh4dt:=-r32;
dboh3dt:=r32;
dhplusdt:=-r29-r30-r31-r32;
dcaco3dt:=-r24;
dcadt:=SD*r24;
dsnh4dt:=r34/(1.0-por(j));
dspo4dt:=r33/(1.0-por(j));
dfesdt:=(1.0*r19/SD-1.0*r17-1.0*r22-1.0*r23-r20);
dfeco3dt:=(1.0*r21/SD);
dmnco3dt:=(1.0*r18/SD);
dfes2dt:=1.0*r22+1.0*r23;
ds0dt:=1.0*SD*r14-4.0*r28-1.0*SD*r23;
dsfedt:=r36/(1.0-por(j));
dsfpdt:=r35/(1.0-por(j));



#  Biogeochemistry - Equilibria
# Specification of equilibrium constraints
# eqrxnsId : set of kinetic reactions which are overuled by a thermodynamic constraint 
# equilibriumseqns[i] : Equilibrium constraint for reaction i  
# 
# example:
# eqrxnID := [r2,rX];
# equilibriumeqns[1] := hplus*co3 - Keq*hco3;
# equilibriumeqns[2] := ...;
eqrxnId := [r29, r30, r31, r32, r33, r34, r35, r36] ;
equilibriumeqns[1] := hplus*hco3-keq1*h2co3 ;
equilibriumeqns[2] := hplus*co3-keq2*hco3  ;
equilibriumeqns[3] := hplus*hs-keq3*h2s ;
equilibriumeqns[4] := hplus*boh4-keq4*boh3;
equilibriumeqns[5] := kspo4*po4*por(j) - spo4*(1. - por(j)) ;
equilibriumeqns[6] := ksnh4*nh4*por(j) - snh4*(1. - por(j)) ;
equilibriumeqns[7] := ksfp*feoh3*po4 - sfp ; #(assume forward reaction is dictated by solid phase)
equilibriumeqns[8] := ksfe*fe*por(j) - sfe*(1. - por(j)) ;

#  Biogeochemistry - Parameters
# Values of rates constants and parameters 
# In this section, all parameters defined in section 'Rate laws' should be defined.
# nparam: number of parameters to define 
# The list is given by bio_name; the values collected in bio_val.
# note that for double precision, 10 should be written as 10.
# 
# example:
# nparam:=4;
# bio_name:=[kmo2hs,kf,kb,Keq];
# vkf :=1.0*10^(5);
# vKeq:=1.0*10^(-10.4);
# vkb :=vkf*vKeq;
# bio_val:=[1000.,vkf,vkb,vKeq];
# 
# 
nparam := 79;
# calculate the dissociation constants, etc, based on temp and salinity (baseline P = 0) from Millero (1995)
lnKC1_0 := 2.83655 - 2307.1266/ T - 1.5529413* ln(T) + (- (0.20760841 + 4.0484/ T)* sqrt(S)) + (0.08468345* S - 0.00654208* S* sqrt(S)) + (ln(1 - 0.001005* S));
lnKC2_0 :=-9.226508 - 3351.6106/ T - 0.2005743* ln(T)+(-0.106901773 - 23.9722/ T)* sqrt(S)+(0.1130822* S - 0.00846934* S^1.5 + ln(1 - 0.001005 * S));
lnKS_0 := 225.838 - 13275.3/T - 34.6435*ln(T) + 0.3449*S^0.5 - 0.0274*S ;
lnKB_0:=(-8966.90-2890.53*sqrt(S)-77.942*S+1.728*S^(3.0/2.0)-0.0996*S*S)/T+(148.0248+137.1942*sqrt(S)+1.62142*S)+((-24.4344-25.085*sqrt(S)-0.2474*S)*ln(T))+ (0.053105*sqrt(S)*T);

# Now apply conversion for P (Millero 1995).
# units are M

deltav:=  -25.5 + 0.1271*T_C + 0.0*T_C*T_C;
deltak:= (-3.08e-3 + 0.0877e-3*T_C + 0.0*T_C*T_C);  
lnkpok0:= -(deltav/(R*T))*P + (0.5*deltak/(R*T))*P*P;
val_kkeq1 := exp(lnKC1_0)*exp(lnkpok0);


deltav:=  -15.82  -0.0219*T_C;
deltak:=  (1.13e-3 -0.1475e-3*T_C);  
lnkpok0:= -(deltav/(R*T))*P + (0.5*deltak/(R*T))*P*P;
val_kkeq2:= exp(lnKC2_0)*exp(lnkpok0);

val_kkeq3 := exp(lnKS_0 + ((14.8  - 0.002*T_C + 0.0004   *T_C^2)*P + 0.5* (2.89 + 0.054*T_C) *P^2/1.01345)/R/T);

deltav:=  -29.48 + 0.1622*T_C -2.608e-3*T_C*T_C;
deltak:= (-2.84e-3+ 0.0794e-3*T_C);  
lnkpok0:= -(deltav/(R*T))*P + (0.5*deltak/(R*T))*P*P;
val_kkeq4:= exp(lnKB_0)*exp(lnkpok0);

# calculate the solubility product for caco3
readlib(log10) :
tmp1 := -171.9065-0.077993*T+2839.319/T+71.595*log10(T);
tmp2 := +(-0.77712+0.0028426*T+178.34/T)*sqrt(S);
tmp3 := -0.07711*S+0.0041249*S^1.5;
log10Kspc := tmp1 + tmp2 + tmp3;
ksca := 10.0^(log10Kspc);

deltav:=   -48.7600 + 0.5304*T_C;
deltak:= (-0.0118 +  3.6920e-4*T_C);  
lnkpok0:= -(deltav/(R*T))*P + (0.5*deltak/(R*T))*P*P;

#________________________________________________________________________________________________


#units mol, l, yr
val_xbiot       :=10.0; # Dom was 5.0 jens 3
val_k1          :=0.174; # Dom was 0.1; Jens 2.0
val_k2          :=0.001; # Dom was 0.001;Jens 0.12
val_kmo2        :=8.0e-9; # Dom was 8.0e-9; Jens 0.8e-6
val_kmno3       :=5.0e-9;
val_kmmno2      :=5.0e-6;
val_kmfeoh3     :=1.25e-5;
val_kmso4       :=100.0e-9;
val_knit        :=1.0e7;
val_kmnox       :=2.0e12;
val_kfemno2     :=2.0e11;
val_kfeo2       :=0.1e12;
val_kh2so2      :=0.1e10;
val_kh2smno2    :=1.0e7;
val_kh2sfeoh3   :=1.0e7;
val_kaom        :=5.0e6;
val_kaomo2      :=1.0e13;
val_kfeso2      :=1.0e9;
val_kmnco3precip:=0.1e-7;
val_KsMnCO3     :=3.2e-15;
val_kfesprecip  :=5.0e-9;
val_kfesdiss    :=1.0e-3;
val_KsFeS       :=6.3e-6/0.77;
val_kfeco3precip:=1.0e-9;
val_KsFeCO3     :=4.0e-15;
val_kpyr        :=60.0e6;
val_kfess0      :=2.5e4;
val_kcaldiss    :=0.1;
val_kapa        :=1.0;
val_po4_eq      :=3.7e-9;
val_kdis        :=0.0e10;
val_ch4eq       :=10000.0;
val_kgas        :=0.0e-3;
val_kdi         :=0.001;
val_h2sstar     :=10.0e-6;
val_ksnh4       :=1.6;
val_kspo4       :=1.3; 
val_ksfp        :=100.0;
val_ksfe        :=400.0;

val_keq1 :=val_kkeq1*1e-3;
val_keq2 := val_kkeq2*1e-3;
val_keq3 := val_kkeq3*1e-3;
val_keq4 := val_kkeq4*1e-3;
val_kf1 := 1.0;
val_kb1 := val_keq1;
val_kf2 := 1.0 ;
val_kb2 := val_keq2;
val_kf3 := 1.0;
val_kb3 := val_keq3;
val_kf4 := 1.0;
val_kb4 := val_keq4;
val_kf5 := val_kspo4;
val_kb5 := 1.0;
val_kf6 := val_ksnh4;
val_kb6 := 1.0;
val_kf7 := val_ksfp;
val_kb7 := 1.0;
val_kf8 := val_ksfe;
val_kb8 := 1.0;
val_kspcal := ksca * exp(lnkpok0)*1e-6;


bio_name := [sw01, sw02, sw03, sw04, sw05, sw06, sw07, sw08, sw09, sw10, sw11, sw12, sw13, sw14, sw15, sw16, sw17, sw18, sw19, xbiot, k1, k2, kmo2, kmno3, kmmno2, kmfeoh3, kmso4, knit, kmnox, kfemno2, kfeo2, kh2so2, kh2smno2, kh2sfeoh3, kaom, kaomo2, kfeso2, kmnco3precip, KsMnCO3, kfesprecip, kfesdiss, KsFeS, kfeco3precip, KsFeCO3, kpyr, kfess0, kcaldiss, kapa, po4_eq, kdis, h2sstar, ch4eq, kgas, kdi, ksnh4, kspo4, ksfp, ksfe, kf1, kb1, kf2, kb2, kf3, kb3, kf4, kb4, kf5, kb5, kf6, kb6, kf7, kb7, kf8, kb8, kspcal, keq1, keq2, keq3, keq4]; 
bio_val :=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, val_xbiot, val_k1, val_k2, val_kmo2, val_kmno3, val_kmmno2, val_kmfeoh3, val_kmso4, val_knit, val_kmnox, val_kfemno2, val_kfeo2, val_kh2so2, val_kh2smno2, val_kh2sfeoh3, val_kaom, val_kaomo2, val_kfeso2, val_kmnco3precip, val_KsMnCO3, val_kfesprecip, val_kfesdiss, val_KsFeS, val_kfeco3precip, val_KsFeCO3, val_kpyr, val_kfess0, val_kcaldiss, val_kapa, val_po4_eq, val_kdis, val_h2sstar, val_ch4eq, val_kgas, val_kdi, val_ksnh4, val_kspo4, val_ksfp, val_ksfe, val_kf1, val_kb1, val_kf2, val_kb2, val_kf3, val_kb3, val_kf4, val_kb4, val_kf5, val_kb5, val_kf6, val_kb6, val_kf7, val_kb7, val_kf8, val_kb8, val_kspcal, val_keq1, val_keq2, val_keq3, val_keq4]; 




# Switches
# Switches can be used in the rate equations. Specify in nswitches, how many switches are in use, name them and define the switch expressions. The switch names must also appear in bio_name and must be assigned a dummy value there. The switch equals 1 if the switch expression is >0, 0 otherwise. To reference the coordinates in the domain, use x_pos, y_pos and z_pos.
nswitches := 19;
switchlist := [sw01, sw02, sw03, sw04, sw05, sw06, sw07, sw08, sw09, sw10, sw11, sw12, sw13, sw14, sw15, sw16, sw17, sw18, sw19] ;
switchcrit := [(xbiot-x_pos), (o2-kmo2), -(sw02-1.0), sw03*(no3-kmno3), -sw03*(sw04-1.0), sw05*(mno2-kmmno2),-sw05*(sw06-1.0), sw07*(feoh3-kmfeoh3), -sw07*(sw08-1.0), sw09*(so4-kmso4), -sw09*(sw10-1.0), ((mn*co3/KsMnCO3)-1.0), ((fe*hs)/(KsFeS*hplus)-1.0), ((fe*co3/KsFeCO3)-1.0), (1.0-(ca*co3/kspcal)), (po4-po4_eq), (ch4-ch4eq), (ch4-0.0), (h2sstar-(h2s+hs))] ;
# 
#  Transport - Molecular Diffusion
# Spec ification of the molecular diffusion coefficients
# diffdata: molecular diffusion coefficient at 0 degree celsius (cm^2/yr)
# alphadata: temperature dependence of the diffusion coefficient (1/K) 
# the in situ molecular diffusion coefficient, corrected for tortuosity is calculated as:
# D(T,sal) = [(0.95-0.001*sal)* D(T=0,sal=0)*(1 + muc*T[C])]/(1-ln(por^2))  
# example:
# diffdata := [100., 304.,0.,0.,100.];
# alphadata:= [0.006, 0.04,0.,0.,.0.05];
# [ch2o, o2, no3, mno2, feoh3, so4, ch4, nh4, po4, mn, fe, h2s, hs, ch4g, h2co3, hco3, co3, boh4, boh3, hplus, caco3, ca, snh4, spo4, fes, feco3, s0, fes2, sfe, sfp, mnco3] 
diffdata := [0.0, 380.449545, 394.5878727, 0.0, 0.0, 173.9205889, 263.9351889, 395.8731752, 112.35777, 123.3890416, 136.2420668, 331.6080494, 392.0172677, 5000.0, 320.0403267, 217.2161254, 176.0864448, 96.29573485, 110.0522684, 600.0, 0.0, 150.3803945, 0.0, 0.0, 0.0, 0.0, 173.9205889, 0.0, 0.0, 0.0, 0.0, 0.0];
alphadata := [0.0, 0.06, 0.038, 0.0, 0.0, 0.045, 0.0520, 0.041, 0.054, 0.05, 0.044, 0.06, 0.031, 0.0, 0.06, 0.048, 0.047, 0.048, 0.048, 0.06, 0.0, 0.045, 0.0, 0.0, 0.0, 0.0, 0.045, 0.0, 0.0, 0.0, 0.0, 0.0] ; 

#  Transport - Boundary Conditions
# Specification of upper and lower boundary conditions for each species.
# There are 3 options 
# 0. kmnown concentration (Dirichlet, M/L^3)
# 1. kmnown concentration gradient (Neumann, M/L^3*L)
# 2. kmnown total (diffusive and advective) flux (Robin, M/L^2/T)
# technical note: option 1 and 2 involve ghost points outside the domain. If the mixing parameters vary with depth, one needs to assign a mixing intensities at the ghost points. By default this is done by linear extrapolation. To ovrwrite this, the user has to edit gridsetup.f and advdiffcoeff.f (both at the bottom; explanations are given there)
# type_up: array defining the type of condition for each species at the upper boundary (0, 1 or 2) 
# bnddata_up: array containing the values specified at the upper boundary.
# type_down: array defining the type of condition for each species at the lower boundary 
# bnddata_down: array containing the values specified at the upper boundary.
# 
# example:
# type_up := [0,0,2];
# bnddata_up := [1.4,0.001,0.001];
# type_down := [1,1,1];
# bnddata_down := [0.,0.,0.];
# [ch2o, o2, no3, mno2, feoh3, so4, ch4, nh4, po4, mn, fe, h2s, hs, ch4g, h2co3, hco3, co3, boh4, boh3, hplus, caco3, ca, snh4, spo4, fes, feco3, s0, fes2, sfe, sfp, mnco3]  ;
type_up := [0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ;
bnddata_up := [1.0417e-04, 300.0e-9, 20.0e-9, 3.25e-6, 12.0e-6, 28000.0e-9, 0.0, 0.0, 0.06e-8, 0.0, 0.0, 0.0, 0.0, 0.0, 2.687126616977748e-08, 2.207192533546868e-06, 1.159362002833545e-07, 5.129572516684913e-08, 3.737042748331509e-07, 1.2303e-11, 0.35, 11.2e-6, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0417e-04] ; 
type_down := [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ;
bnddata_down := [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ;
#  Transport - Identify species that are not transported
# nrnotransp: number of species that are not getting transported at all
# notransp: list of species that are not getting transported at all 
# (i.e. for which the above defined dispersion and advection velocities are not used!)
# 
# example:
# nrnotransp:=3;
# notransp := [1,3,4];
nrnotransp := 1;
listnotransp := [20] ;
#  Initial conditions
#  vic: options to select initial guesses
# - 1: read from a file "initialconc.txt", which contains the columns 
#      z, conc(species1), conc(species2)....
# - 2: fixed concentration, defined in the array iniconc
# - 3: individual files for each species. listinput contains the species number (lenght: number of species), filenames are given in "file_in_names".inp, containing in column 1 the conc, in column 2 depth. 
# 
# example (only relevant input data is given for each option, assuming 17 species):
# NOTE that you have to provide something for all ncomp species in the arrays iniconc, listinput and file_in_names, EVEN IF YOU DON'T USE IT WITH THE OPTION YOU SELECTED
# vic:=1;
# vic:=2;iniconc:=[0.,2.,...]; 
# vic:=3;listinput:=[1,2,...,17];file_in_names:=[o2,no3,...,sp17];
vic := 2 ;
iniconc := [0.5*30.19744e-6, 0.0, 0.0, 0.0, 0.0, 28000.0e-9, 0.0, 0.0, 1.0e-9, 0.0, 0.0, 0.0, 0.0, 0.0, 27.111e-9, 2143.30e-9, 79.582e-9, 50.298e-9, 374.70e-9, 1.5849e-8*1e-3, 20.0e-5, 1e-5, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5*30.19744e-6] ;
listinput := [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32] ;
file_in_names := [dummy1,dummy2,dummy3,dummy4,dummy5,dummy6,dummy7,dummy8, dummy9, dummy10, dummy11, dummy12, dummy13, dummy14, dummy15, dummy16, dummy17, dummy18, dummy19, dummy20, dummy21, dummy22, dummy23, dummy24, dummy25, dummy26, dummy27, dummy28, dummy29, dummy30, dummy31, dummy32] ;
#  Output
# noutput: number of species to be printed
# nroutput: number of rates to be printed
# listoutput: species number to print
# listroutput: rate number to print
# file_names: Respective file name for each of the species to print
# file_rnames: Respective file name for each of the rates to print
# time_iniout: First time (in years) for which a printout is requested
# time_intvout: time interval (in years) at which the printing is performed, starting from time_iniout
# 
# example:
# noutput:=4;
# listoutput:=[1,2,3,5];
# file_names:=[o2, so4, MnOx, hco3];
# time_iniout:=10.;
# time_intvout:= 100.;
# 
noutput := 32;
nroutput := 43;
listoutput := [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32] ;
listroutput := [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43] ;
file_names := [g1, zzo2, zno3, mno2, feoh3, zso4, zch4, znh4, zpo4, zzmn, zzfe, zh2s, zzhs, ch4g, h2co, hco3, zco3, boh4, boh3, zzph, caco, zzca, snh4, spo4, zfes, feco, zzs0, fes2, zsfe, zfp, mnco, g2] ;
file_rnames := [xrate1, xrate2, xrate3, xrate4, xrate5, xrate6, xrate7, xrate8, xrate9, xrate10, xrate11, xrate12, xrate13, xrate14, xrate15, xrate16, xrate17, xrate18, xrate19, xrate20, xrate21, xrate22, xrate23, xrate24, xrate25, xrate26, xrate27, xrate28, xrate29, xrate30, xrate31, xrate32, xrate33, xrate34, xrate35, xrate36, xrate37, xrate38, xrate39, xrate40, xrate41, xrate42, xrate43] ;
time_iniout := 0.01 ;
time_intvout := 20.0;#tot_time/10.0;
#  Optimization
# there are several optimization options available. here you can specify what needs to be optimized and where the data is stored. to identify what kind of algorithm you want to use please select the appropriate options in drivervalues.f
# nopt_v: number of parameters to be optimized
# ntopt_v: number of time points where measurements are available
# nparam_opt: total number of parameters. can include the nparam but also the physical parameters.
# maxxmeas_v: maximum number of depth points at any given time measured (used to make array sizes)
# maxspmeas_v: maximum number of species measured at any given time (used to make array sizes)
# opt_name: names of the parameters. they have to match the names given above, so best you make a copy paste!
# idpar_v: identify the parameters to be optimized from the parameter list opt_name
# filemeas_name: name of the files containing the measured data
# if nopt_v is set to 0 then the rest of the input doesn't matter
# example:
# nopt_v := 2; # number of parameters to be optimized
# ntopt_v := 3; # number of timepoints with measurements
# nparam_opt := nparam; # total number of parameters, set equal to all parameters except physcial ones
# maxxmeas_v:=20; # maximum 20 points in a profile at any given time
# maxspmeas_v :=2; # maximum 2 species measured at one timepoint
# opt_name := bioname; # (note that this does not include the physical parameters! if you want them to be adapted you need to specify them explicitly)
# idpar_v:= [1,3]; # optimize parameters 1 and 3 in the above list
# filemeas_name:= [meas1.dat, meas2.dat, meas3.dat]; # filenames with measurements at timepoints
nopt_v := 0 ;
ntopt_v := 0 ;
nparam_opt := 0 ;
maxxmeas_v := 0;
maxspmeas_v := 0;
opt_name := [];
idpar_v := [];
filemeas_name := [];
#  Maple specific info
# dir_f: directory where the FORTRAN routines and Maple spread.m files are parsed
# format Mac: "Macinthosh HD:UU:...:code"
# format PC: "C:\\maple\\...\\code"
# WAS: dir_f := "C:\\Dokumente und Einstellungen\\centler\\Desktop\\Labor\\Simulations": 
# currentdir(dir_f):
# save "spread.m" ;
dir_f := "/home/dh12423/Documents/5_BRNS/BRNS_full_for_OMEN_15062016":
currentdir(dir_f) :
parse(sprintf("save %q,\"spread.m\";",anames()), statement) ;

"now execute processor - make sure the directories are set correctly";
# 
# 
# 
# 
# 
# 
# 
# 
# ACG
