!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!     OMEN-SED model
!!!     Authors: Dominik Hülse, Sandra Arndt, Stuart Daines
!!!     Date: June 2015

!!!     modelled stuff: (labile, refractory) TOC, O2, Sulfate (SO4), Hydrogen sulfide (H2S),
!!!                     Nitrate (NO3), Ammonium (NH4), todo: Phosphate (PO4), Methane (CH4 ) implicitly

!!!     gets bottom water concentrations from an Ocean/Earth System Model (e.g. GENIE) calculates burial flux of TOC, sediment water interfaces fluxes of other species
!!!     and gives this back to the other model

!!!     zbrent:     root-finding algorithm combining the bisection method, the secant method and inverse quadratic interpolation.
!!!                 It has the reliability of bisection but it can be as quick as some of the less reliable methods.
!!!                 The algorithm tries to use the potentially fast-converging secant method or inverse quadratic interpolation if possible,
!!!                 but it falls back to the more robust bisection method if necessary.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE benthic
    implicit none

    !sediment characteristics
    real rho_sed                            ! sediment density (g/cm3)
    real wdepth                            ! water depth (m)
    real w                                      ! burial velocity  (cm/yr)
    real z0, zox                                ! surface
    real zbio                                ! bioturbation depth (cm)

    real zinf                               !Inifinity (cm)
    !zinf = 1000
    !zlow=100
    real Dbio                                 !bioturbation coefficient (cm2/yr)
    real por                                !porosity (-)
    real tort                              !tortuosity (-)
    real irrigationFactor                   !irrigation factor (-)
    real dispFactor                             !dispersion factor (-)

    !stoichiometric factors
    real SD                                     !volume factor solid->dissolved phase
    real OC                                     !O2/C (mol/mol)
    real NC1                                    !N/C first TOC fraction (mol/mol)
    real NC2                                    !N/C second TOC fraction (mol/mol)
    real PC1                                    !P/C first TOC fraction (mol/mol)
    real PC2                                    !P/C second TOC fraction (mol/mol)
    real SO4C                                   !SO4/C (mol/mol)
    real DICC1                                  !DIC/C until zSO4 (mol/mol)
    real DICC2                                  !DIC/C below zSO$ (mol/mol)
    real MC                                     !CH4/C (mol/mol)
    real gamma                                  !fraction of NH4 that is oxidised in oxic layer
    real gammaH2S                               !fraction of H2S that is oxidised in oxic layer
    real gammaCH4                               !fraction of CH4 that is oxidised at SO4
    real satSO4                                 ! SO4 saturation
    real NO3CR                                  ! NO3 consumed by Denitrification
    real ALKROX;                                 ! Aerobic degradation
    real ALKRNIT;                                ! Nitrification
    real ALKRDEN;                                ! Denitrification
    real ALKRSUL;                                ! Sulfato reduction
    real ALKRH2S;                                ! H2S oxydation (CHECK THIS VALUE!!!)
    real ALKRMET;                                ! Methanogenesis
    real ALKRAOM;                                ! AOM

    real zoxgf                            ! cm, rolloff NH4, H2S oxidation for small zox depth

    !bottom water concentrations
    real T                           !temperature (degree C)
    real C01                         !TOC concentration at SWI (wt!) -> (mol/cm3 bulk phase)
    real C02                         !TOC concentration at SWI (wt!) -> (mol/cm3 bulk phase)
    real O20                                               !O2  concentration at SWI (mol/cm3)
    real SO40                                             !SO4 concentration at SWI (mol/cm3)
    real H2S0                                               !H2S concentration at SWI (mol/cm3)
    real DIC0                                             !DIC concentration at SWI (mol/cm3)
    real ALK0                                             !ALK concentration at SWI (mol/cm3)
    real NO30                                               !NO3 concentration at SWI (mol/cm3)
    real NH40                                               !NH4 concentration at SWI (mol/cm3)
    real PO40                                                  !PO4 concentration at SWI (mol/cm3)
    real Mflux0                       ! flux of M to the sediment (mol/(cm2*yr))
    real S0                                                     !Salinity at SWI


    ! ORGANIC MATTER
    real DC1                                           !TOC diffusion coefficient (cm2/yr)
    real C, C1, C2
    real k1                                             !TOC degradation rate constnat (1/yr)
    real k2                                              !TOC degradation rate constant (1/yr)

    ! O2
    real qdispO2                          !O2 diffusion coefficient in water (cm2/yr)
    real adispO2                          !O2 linear coefficient for temperature dependence (cm2/yr/oC)
    real DO21                             !O2 diffusion coefficient in bioturbated layer (cm2/yr)
    real DO22                             !O2 diffusion coefficient in non-bioturbated layer (cm2/yr)
    real r_zxf                            !roll off oxidation at low zox

    real aa11, bb11, aa21, A11, A21, aa12, bb12, aa22, A12, A22
    real ls_a, ls_b, ls_c, ls_d, ls_e, ls_f

    ! Nitrate (NO3)
    real qdispNO3                 ! NO3 diffusion coefficient in water (cm2/yr)
    real adispNO3                 ! NO3 linear coefficient for temperature dependence (cm2/yr/oC)
    real DN1                     ! NO3 diffusion coefficient in bioturbated layer (cm2/yr)
    real DN2                      ! NO3 diffusion coefficient in non-bioturbated layer (cm2/yr)
    real zno3
    real KNH4                     ! Adsorption coefficient (same in ocix and anoxic layer) (-)

    ! Sulfate (SO4)
    real qdispSO4                 ! SO4 diffusion coefficient in water (cm2/yr)
    real adispSO4                 ! SO4 linear coefficient for temperature dependence (cm2/yr/oC)
    real DSO41                    ! SO4 diffusion coefficient in bioturbated layer (cm2/yr)
    real DSO42                    ! SO4 diffusion coefficient in non-bioturbated layer (cm2/yr)
    real zso4

    ! Ammonium (NH4)
    real qdispNH4                 ! NH4 diffusion coefficient in water (cm2/yr)
    real adispNH4                 ! NH4 linear coefficient for temperature dependence (cm2/yr/oC)
    real DNH41                    ! NH4 diffusion coefficient in bioturbated layer (cm2/yr)
    real DNH42                    ! NH4 diffusion coefficient in non-bioturbated layer (cm2/yr)

    ! Hydrogen sulfide (H2S)
    real qdispH2S                 ! H2S diffusion coefficient in water (cm2/yr)
    real adispH2S                 ! H2S linear coefficient for temperature dependence (cm2/yr/oC)
    real DH2S1                    ! H2S diffusion coefficient in bioturbated layer (cm2/yr)
    real DH2S2                    ! H2S diffusion coefficient in non-bioturbated layer (cm2/yr)

    ! Phosphate (PO4)
    real qdispPO4                 ! PO4 diffusion coefficient in water (cm2/yr)
    real adispPO4                 ! PO4 linear coefficient for temperature dependence (cm2/yr/oC)
    real DPO41                    ! PO4 diffusion coefficient in bioturbated layer (cm2/yr)
    real DPO42                    ! PO4 diffusion coefficient in non-bioturbated layer (cm2/yr)
    real KPO4_ox                    ! Adsorption coefficient in oxic layer (-)
    real KPO4_anox                    ! Adsorption coefficient in anoxic layer (-)
    real ksPO4                    ! Rate constant for kinetic P sorption (1/yr)
    real kmPO4                    ! Rate constant for Fe-bound P release upon Fe oxide reduction
    real kaPO4                    ! Rate constant for authigenic P formation (1/yr)
    real PO4s                     ! Equilibrium concentration for P sorption (mol/cm3)
    real PO4a                     ! Equilibrium concentration for authigenic P formation (mol/cm3)
    real Minf                     ! asymptotic concentration for Fe-bound P (mol/cm3)

    ! DIC
    real qdispDIC                   ! DIC diffusion coefficient in water (cm2/yr)
    real adispDIC                   ! DIC linear coefficient for temperature dependence (cm2/yr/oC)
    real DDIC1                      ! DIC diffusion coefficient in bioturbated layer (cm2/yr)
    real DDIC2                      ! DIC diffusion coefficient in non-bioturbated layer (cm2/yr)

    ! Alkalinity
    real qdispALK                   ! ALK diffusion coefficient in water (cm2/yr)
    real adispALK                   ! ALK linear coefficient for temperature dependence (cm2/yr/oC)
    real DALK1                      ! ALK diffusion coefficient in bioturbated layer (cm2/yr)
    real DALK2                      ! ALK diffusion coefficient in non-bioturbated layer (cm2/yr)

CONTAINS

    SUBROUTINE initialize()
        !   __________________________________________________________
        !
        !   initalize
        !   __________________________________________________________

        !   MOST OF THE FOLLOWING VALUES WILL BE PASSED DOWN FROM GENIE
        ! *****************************************************************

        print*, ' '
        print*, '----------- start initialize --------------'
        print*, ' '

        rho_sed=2.6                               ! sediment density (g/cm3)
        wdepth=600.0                                ! water depth (m)
        z0  = 0.0                                  ! surface
        zox = 0.0
        zbio=10.0                                  ! bioturbation depth (cm)

        zinf=100.0                                  !Inifinity (cm)
        Dbio = 5.2*(10.0**(0.7624-0.0003972*wdepth))  !bioturbation coefficient (cm2/yr) - after Middelburg at al. 1997
        por=0.85                                !porosity (-) defined as: porewater_vol./(solid_sed_vol.+porewater_vol.)
        tort=3.0                               !tortuosity (-)
        irrigationFactor=1.0

        gamma=0.9                      !fraction of NH4 that is oxidised in oxic layer
        gammaH2S=0.95                           !fraction of H2S that is oxidised in oxic layer
        gammaCH4=0.99                   !fraction of CH4 that is oxidised at SO4
        satSO4=0.0                      ! SO4 saturation

        zoxgf = 0.1                         ! cm, rolloff NH4, H2S oxidation for small zox depth

        !bottom water concentrations
        T=8.0                                                     ! temperature (degree C)
        C01=1.0*1e-2/12*rho_sed                                ! TOC concentration at SWI (wt!) -> (mol/cm3 bulk phase)
        C02=1.0*1e-2/12*rho_sed                                ! TOC concentration at SWI (wt!) -> (mol/cm3 bulk phase)
        O20=300.0e-9         !was  300.0e-9                                     ! O2  concentration at SWI (mol/cm3)
        NO30=40.0e-9                                               ! NO3 concentration at SWI (mol/cm3)
        NH40=0.0e-9                                                ! NH4 concentration at SWI (mol/cm3)
        SO40=28000.0e-9      !  28000.0e-9                                     ! SO4 concentration at SWI (mol/cm3)
        H2S0=0.0e-13           !was 0.0e-9                                 ! H2S concentration at SWI (mol/cm3)
        PO40=40.0e-9                                                  ! PO4 concentration at SWI (mol/cm3)
        Mflux0=365*0.2e-10                                          ! flux of M to the sediment (mol/(cm2*yr))
        DIC0 = 2000.0e-9                                             ! DIC concentration at SWI (mol/cm^3)
        ALK0=2400.0e-9                                             !ALK concentration at SWI (mol/cm^3)

        w=10.0**(-0.87478367-0.00043512*wdepth)*3.3             ! sedimentation rate, cm/yr - after Middelburg at al. 1997
        dispFactor=por**(tort-1.0)*irrigationFactor             ! dispersion factor (-) - Ausbreitung - type of mixing that accompanies hydrodynamic                                    ! flows -> ~builds out paths of flow
        SD=(1-por)/por                                          ! volume factor solid->dissolved phase
        OC=1.0*SD                                               ! O2/C (mol/mol)
        NC1=0.1509*SD                                           ! N/C first TOC fraction (mol/mol)
        NC2=0.13333*SD                                          ! N/C second TOC fraction (mol/mol)
        PC1=0.0094*SD                                          ! P/C first TOC fraction (mol/mol)
        PC2=0.0094*SD                                          ! P/C second TOC fraction (mol/mol)
        SO4C=0.5*SD                                             ! SO4/C (mol/mol)
        DICC1=1.0*SD                                           ! DIC/C until zSO4 (mol/mol)
        DICC2=0.5*SD                                           ! DIC/C below zSO4 (mol/mol)
        MC=0.5*SD                                              ! CH4/C (mol/mol)
        NO3CR=(94.4/106)*SD                                    ! NO3 consumed by Denitrification
        ALKROX=15.0/106;                                        ! Aerobic degradation
        ALKRNIT=-2.0;                                           ! Nitrification
        ALKRDEN=93.4/106;                                       ! Denitrification
        ALKRSUL=15.0/106;                                       ! Sulfato reduction
        ALKRH2S=-1.0;                                           ! H2S oxydation (CHECK THIS VALUE!!!)
        ALKRMET=14.0/106;                                       ! Methanogenesis
        ALKRAOM=2.0;                                            ! AOM

        ! ORGANIC MATTER
        DC1 = Dbio
        k1=0.01
        k2=0.0001

        ! O2
        qdispO2=348.62172
        adispO2=14.08608

        DO21=(qdispO2+adispO2*T)*dispFactor+Dbio
        DO22=(qdispO2+adispO2*T)*dispFactor

        r_zxf=0.0

        ! Nitrate (NO3)
        qdispNO3=308.42208
        adispNO3=12.2640
        DN1=(qdispNO3+adispNO3*T)*dispFactor+Dbio
        DN2=(qdispNO3+adispNO3*T)*dispFactor
        zno3 = 0.0
        KNH4 = 1.3                                      !Adsorption coefficient (same in ocix and anoxic layer) (-)

        ! Sulfate (SO4)
        qdispSO4=309.0528                               ! SO4 diffusion coefficient in water (cm2/yr)
        adispSO4=12.2640                                ! SO4 linear coefficient for temperature dependence (cm2/yr/oC)
        DSO41=(qdispSO4+adispSO4*T)*dispFactor+Dbio     ! SO4 diffusion coefficient in bioturbated layer (cm2/yr)
        DSO42=(qdispSO4+adispSO4*T)*dispFactor          ! SO4 diffusion coefficient in non-bioturbated layer (cm2/yr)
        zso4 = 0.0

        ! Ammonium (NH4)
        qdispNH4=309.0528
        adispNH4=12.2640
        DNH41=((qdispNH4+adispNH4*T)*dispFactor+Dbio)/(1.0+KNH4)
        DNH42=((qdispNH4+adispNH4*T)*dispFactor)/(1.0+KNH4)

        ! Hydrogen sulfide (H2S)
        qdispH2S=309.0528
        adispH2S=12.2640
        DH2S1=(qdispH2S+adispH2S*T)*dispFactor+Dbio
        DH2S2=(qdispH2S+adispH2S*T)*dispFactor

        ! Phosphate (PO4)
        qdispPO4=112.90764
        adispPO4=5.586252
        DPO41=((qdispPO4+adispPO4*T)*dispFactor+Dbio)               ! PO4 diffusion coefficient in bioturbated layer (cm2/yr)
        DPO42=((qdispPO4+adispPO4*T)*dispFactor);                   ! PO4 diffusion coefficient in non-bioturbated layer (cm2/yr)
        KPO4_ox = 200.0                    ! Adsorption coefficient in oxic layer (-)
        KPO4_anox = 1.3                   ! Adsorption coefficient in anoxic layer (-)
        ksPO4 = 1.0                   ! Rate constant for kinetic P sorption (1/yr)
        kmPO4 = 2.2e-6*24*365                   ! Rate constant for Fe-bound P release upon Fe oxide reduction
        kaPO4 = 10.0                   ! Rate constant for authigenic P formation (1/yr)
        PO4s = 1.0e-9                    ! Equilibrium concentration for P sorption (mol/cm3)
        PO4a = 0.5e-8                    ! Equilibrium concentration for authigenic P formation (mol/cm3)
        Minf = 1.0e-10                    ! asymptotic concentration for Fe-bound P (mol/cm3)

        ! DIC
        qdispDIC=309.0528
        adispDIC=12.2640
        DDIC1=(qdispDIC+adispDIC*T)*dispFactor+Dbio                 ! DIC diffusion coefficient in bioturbated layer (cm2/yr)
        DDIC2=(qdispDIC+adispDIC*T)*dispFactor                      ! DIC diffusion coefficient in non-bioturbated layer (cm2/yr)



    end SUBROUTINE initialize

    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************

    !                           TOC

    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE zTOC()
        !   __________________________________________________________
        !
        !   calculate benthic burial/recycling fluxes (see documentation for details!)
        !   __________________________________________________________

        !   organic matter burial - 2 fractions

        ! local variables
        real loc_POC1_conc_zinf, loc_POC2_conc_zinf, dum_sed_pres_fracC
        !    real aa11, bb11, aa21, A11, A21, aa12, bb12, aa22, A12, A22
        real dC1dz, C1flx, dC2dz, C2flx, Cflx             ! Cflx: Sed input flux to upper boundary, per cm^2 water column
        real F_TOC1, F_TOC2, F_TOC                        ! Flux through lower boundary zinf, per cm^2 water-column


        !    print*,' ------------------ START zTOC ---------------------'

        aa11 = (w-sqrt(w**2+4*DC1*k1))/(2*DC1)
        bb11 = (w+sqrt(w**2+4*DC1*k1))/(2*DC1)
        aa21 = (-k1/w)
        A11 = -(C01*bb11*exp(bb11*zbio))/(aa11*exp(aa11*zbio)-bb11*exp(bb11*zbio))
        A21=(A11*(exp(aa11*zbio)-exp(bb11*zbio))+C01*exp(bb11*zbio))/exp(aa21*zbio)

        aa12 = (w-sqrt(w**2+4*DC1*k2))/(2*DC1)
        bb12 = (w+sqrt(w**2+4*DC1*k2))/(2*DC1)
        aa22 = (-k2/w)
        A12=-(C02*bb12*exp(bb12*zbio))/(aa12*exp(aa12*zbio)-bb12*exp(bb12*zbio))
        A22=(A12*(exp(aa12*zbio)-exp(bb12*zbio))+C02*exp(bb12*zbio))/exp(aa22*zbio)


        ! Cflx: Sed input flux to upper boundary, per cm^2 water column
        if(z0<zbio) then
            C1=A11*(exp(aa11*z0)-exp(bb11*z0))+C01*exp(bb11*z0)
            C2=A12*(exp(aa12*z0)-exp(bb12*z0))+C02*exp(bb12*z0)
        else
            C1=A21*exp(aa21*z0)
            C2=A22*exp(aa22*z0)
        end if
        C = C1 + C2

        print*, 'C = C1 + C2 ', C
        print*, ' '

        if(z0 < zbio) then
            dC1dz =  A11*(aa11*exp(aa11*z0)-bb11*exp(bb11*z0))+C01*bb11*exp(bb11*z0)
            C1flx = - (1-por)*(-DC1*dC1dz + w*C1)
            dC2dz =  A12*(aa12*exp(aa12*z0)-bb12*exp(bb12*z0))+C02*bb12*exp(bb12*z0)
            C2flx = - (1-por)*(-DC1*dC2dz + w*C2)
        else
            C1flx = - (1-por)*w*C1
            C2flx = - (1-por)*w*C2
        end if
        Cflx = C1flx + C2flx

        print*, 'Cflx', char(9), Cflx
        print*, 'C1flx', char(9), C1flx
        print*, 'C2flx', char(9), C2flx


        ! Flux through lower boundary zinf, per cm^2 water-column
        F_TOC1 = -(1-por)*w*A21*exp(aa21*zinf)
        F_TOC2 = -(1-por)*w*A22*exp(aa22*zinf)
        F_TOC = F_TOC1 + F_TOC2

        !    print*, ' '
        print*, 'F_TOC1', char(9), F_TOC1
        print*, 'F_TOC2', char(9), F_TOC2
        print*, 'F_TOC', char(9), F_TOC

        ! Concentration at lower boundary zinf
        if(zinf<zbio) then
            loc_POC1_conc_zinf=A11*(exp(aa11*zinf)-exp(bb11*zinf))+C1*exp(bb11*zinf)
            loc_POC2_conc_zinf=A12*(exp(aa12*zinf)-exp(bb12*zinf))+C2*exp(bb12*zinf)
        else
            loc_POC1_conc_zinf=A21*exp(aa21*zinf)
            loc_POC2_conc_zinf=A22*exp(aa22*zinf)
        end if

        ! DH: need to give back fraction buried of initially deposited (so fraction of the input values to this subroutine)
        !        print*, 'loc_POC1_conc_zinf ', char(9), loc_POC1_conc_zinf
        !        print*, 'loc_POC2_conc_zinf ', char(9), loc_POC2_conc_zinf

        dum_sed_pres_fracC = (loc_POC1_conc_zinf+loc_POC2_conc_zinf)/(C1+C2)
        print*,'Fraction POC-preserved/POC-deposited =' , dum_sed_pres_fracC
        print*,' '

    end SUBROUTINE zTOC


    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------



    FUNCTION calcReac(zU, zL, reac1, reac2)

        real,intent(in):: zU, zL, reac1, reac2
        real calcReac

        ! Integral of reacted organic matter from zU to zL,
        ! multiplied by stoichiometric factors reac1, reac2 (for the two OC phases)

        ! Vector-friendly way of handling 3 cases:
        ! 1) wholly within bioturbated layer:    calcReac_l1(zU,zL)     + (0 =) calcReac_l2(zbio, zbio)
        ! 2) wholly within non-bio     layer:  (0=) calcReac_l1(zbio, zbio) +   calcReac_l2(zU, zL)
        ! 3) crossing zbio                       calcRead_l1(zU,zbio)   +       calcReac_l2(zbio, zL)

        calcReac = calcReac_l1(min(zU,zbio), min(zL,zbio), reac1, reac2) &
        + calcReac_l2(max(zU,zbio), max(zL, zbio), reac1, reac2)

    ! TODO confirm (1-por)*  has been added (to k1 & k2 ?)

    END FUNCTION calcReac

    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    FUNCTION calcReac_l1(zU, zL, reac1, reac2)

        real,intent(in):: zU, zL, reac1, reac2
        real calcReac_l1, reacf1, reacf2

        reacf1 = k1*reac1
        reacf2 = k2*reac2
        calcReac_l1 = -reacf1*(A11*(exp(aa11*zU)*bb11 - exp(bb11*zU)*aa11 - exp(aa11*zL)*bb11 + exp(bb11*zL)*aa11) &
        + C01*exp(bb11*zU)*aa11 - C01*exp(bb11*zL)*aa11)/(aa11*bb11) &
        -reacf2*(A12*(exp(aa12*zU)*bb12 - exp(bb12*zU)*aa12 - exp(aa12*zL)*bb12 + exp(bb12*zL)*aa12) &
        + C02*exp(bb12*zU)*aa12 - C02*exp(bb12*zL)*aa12)/(aa12*bb12)


    END FUNCTION calcReac_l1

    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    FUNCTION calcReac_l2(zU, zL, reac1, reac2)

        real,intent(in):: zU, zL, reac1, reac2
        real calcReac_l2, reacf1, reacf2

        reacf1 = k1*reac1
        reacf2 = k2*reac2

        calcReac_l2 = -reacf1*A21*(exp(aa21*zU) - exp(aa21*zL))/aa21 &
        -reacf2*A22*(exp(aa22*zU) - exp(aa22*zL))/aa22

    !    print*,'calcReac_l2', calcReac_l2

    END FUNCTION calcReac_l2



    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------



    SUBROUTINE calcfg_l12(z, reac1, reac2, ktemp, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ls_D1, ls_D2,&
    ltype, e, dedz, f, dfdz, g, dgdz)
            ! calculate solution basis functions, for layer which may cross bioturbation boundary
            !
            ! reac1, reac2        - mol/mol S released per organic carbon C
            !
            ! General solution for solute S is given by
            !  S(z) = A * e(z) + B * f(z) + g(z)
            !
            ! Where e,f,g are generated by matching solutions across bioturbation boundary (if necessary)
            ! Solution properties (matching etc) are input in ls
            ! On input, ls should contain fields generated by prepfg_l12

        real,intent(in)::       z, reac1, reac2, ktemp
        INTEGER, INTENT(in)::   ltype
        real,intent(inout)::    e, dedz, f, dfdz, g, dgdz
        real                    ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ls_D1, ls_D2

        ! local variables
        real e_1, f_1, g_1, dedz_1, dfdz_1, dgdz_1
        !real ls_a, ls_b, ls_c, ls_d, ls_e, ls_f


        select case (ltype)
            case (1)    ! bioturbated
                !                print*, 'calcfg_l12 CASE 1 bioturbated'
                call calcfg_l1(z, reac1, reac2, ls_D1, ktemp, e, dedz, f, dfdz, g, dgdz)
            case (2)    ! not bioturbated
                call calcfg_l2(z, reac1, reac2, ls_D2, ktemp, e, dedz, f, dfdz, g, dgdz)
            case (3)    ! crossing boundary
                !                   print*, 'calcfg_l12 CASE 3 crossing boundary'
                IF(z >= zbio) THEN      ! below bioturbated region
                    call calcfg_l2(z, reac1, reac2, ls_D2, ktemp, e, dedz, f, dfdz, g, dgdz)
                else    ! above bioturbated region
                    !                    print*, 'CASE 3 crossing boundary ELSEEEEE'
                    call calcfg_l1(z, reac1, reac2, ls_D1, ktemp, e_1, dedz_1, f_1, dfdz_1, g_1, dgdz_1)
                    !                    print*, 'calcfg_l1111111111111: z, reac1, reac2, DO21, ktemp, e_1, dedz_1, f_1, dfdz_1,&
                    !                            & g_1, dgdz_1', z, reac1, reac2, DO21, ktemp, e_1, dedz_1, f_1, dfdz_1, g_1, dgdz_1
                    ! Now find 'transformed' basis functions such that in layer 1, O2 = A_2*et + B_2*ft + gt
                    ! (ie layer 1 soln written in terms of layer 2 coeffs A_2, B_2)
                    call xformsoln(e_1, f_1, g_1, dedz_1, dfdz_1, dgdz_1, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f,&
                    e, f, g, dedz, dfdz, dgdz)
                end if
            case default
                STOP
        end select

    END SUBROUTINE calcfg_l12

    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------


    SUBROUTINE calcfg_l1(z, reac1, reac2, Dtemp, ktemp, e, dedz, f, dfdz, g, dgdz)
        ! Basis functions for solutes, case z <= zbio
        !
        ! reac1, reac2        - mol./mol S released per organic carbon C
        !
        ! General solution for solute S is given by
        !  S(z) = A * e(z) + B * f(z) + g(z)

        real z, reac1, reac2, Dtemp, ktemp                                                          ! in from SUBROUTINE before
        real,INTENT(inout)::  e, dedz, f, dfdz, g, dgdz             ! out

        ! local variables
        real b1, pfac, PhiI1, PhiII1, PhiIII1, PhiI2, PhiII2, PhiIII2
        real ea11z, eb11z, ea12z, eb12z


        e = 1.0
        dedz = 0.0
        b1=w/Dtemp
        f=exp(z*b1)
        dfdz = b1*exp(z*b1)

        pfac = 1                    ! in fact, already has (1-por)/por

        PhiI1 = -pfac*k1*(reac1)*A11/(Dtemp*aa11**2-w*aa11-ktemp)
        PhiII1  = pfac*k1*(reac1)*A11/(Dtemp*bb11**2-w*bb11-ktemp)
        PhiIII1 =-pfac*k1*(reac1)*C01/(Dtemp*bb11**2-w*bb11-ktemp)
        PhiI2   =-pfac*k2*(reac2)*A12/(Dtemp*aa12**2-w*aa12-ktemp)
        PhiII2  = pfac*k2*(reac2)*A12/(Dtemp*bb12**2-w*bb12-ktemp)
        PhiIII2 =-pfac*k2*(reac2)*C02/(Dtemp*bb12**2-w*bb12-ktemp)


        ea11z = exp(aa11*z)
        eb11z = exp(bb11*z)
        ea12z = exp(aa12*z)
        eb12z = exp(bb12*z)

        g =  PhiI1*ea11z + PhiII1*eb11z + PhiIII1*eb11z + &
        PhiI2*ea12z + PhiII2*eb12z + PhiIII2*eb12z


        dgdz = PhiI1*aa11*ea11z + PhiII1*bb11*eb11z + PhiIII1*bb11*eb11z + &
        PhiI2*aa12*ea12z + PhiII2*bb12*eb12z + PhiIII2*bb12*eb12z

    !            print*, 'INPUT calcfg_l1', z, reac1, reac2, Dtemp, ktemp
    !            print*, 'IN  calcfg_l1 g', g
    !            print*, 'IN  calcfg_l1 dgdz', dgdz

    end SUBROUTINE calcfg_l1

    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE calcfg_l2(z, reac1, reac2, Dtemp, ktemp, e, dedz, f, dfdz, g, dgdz)
        ! Basis functions for solutes, case z > zbio
        ! reac1, reac2        - mol/mol S released per organic carbon C
        !
        ! General solution for solute S is given by
        !  S(z) = A * e(z) + B * f(z) + g(z)

        real z, reac1, reac2, Dtemp, ktemp                                                          ! in from SUBROUTINE before
        real,INTENT(inout)::  e, dedz, f, dfdz, g, dgdz             ! out

        ! local variables
        real b2, pfac, PhiI1, PhiI2
        real ea11z, eb11z, ea12z, eb12z

        e = 1.0
        dedz = 0.0
        b2 = w/Dtemp
        f=exp(z*b2)
        dfdz = b2*exp(z*b2)

        !pfac=1./por;   ! assume org matter already *(1-por)
        pfac = 1            !in fact, already has (1-por)/por



        PhiI1 = -pfac*k1*(reac1)*A21/(Dtemp*aa21**2-w*aa21-ktemp)
        PhiI2 = -pfac*k2*(reac2)*A22/(Dtemp*aa22**2-w*aa22-ktemp)
!                print*, 'IN  calcfg_l2 z, PhiI1', z, PhiI1
!                print*, ' pfac, k1, reac1, A21, Dtemp, aa21, w, aa21, ktemp ',  pfac, k1, reac1, A21, Dtemp, aa21, w, aa21, ktemp
!                print*, ' A22', A22
        g = PhiI1*exp(aa21*z) + PhiI2*exp(aa22*z)
        dgdz = PhiI1*aa21*exp(aa21*z) + PhiI2*aa22*exp(aa22*z)


    end SUBROUTINE calcfg_l2


    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE prepfg_l12(reac1, reac2, ktemp, zU, zL, D1, D2, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ltype)
        real reac1, reac2, ktemp, zU, zL, D1, D2
        real e_zbio_l1, dedz_zbio_l1, f_zbio_l1, dfdz_zbio_l1, g_zbio_l1, dgdz_zbio_l1
        real e_zbio_l2, dedz_zbio_l2, f_zbio_l2, dfdz_zbio_l2, g_zbio_l2, dgdz_zbio_l2
        real ls_a, ls_b, ls_c, ls_d, ls_e, ls_f

        INTEGER,INTENT(inout):: ltype

        if(zL <= zbio)then  ! wholly within bioturbated layer
            ltype = 1
        elseif(zU >= zbio) then ! wholly within non-bioturbated layer
            ltype = 2
        else             ! crossing boundary - sort out solution matching at zbio
!            print*, 'IN  CROSS BOUNDARY CASE: zL ', zL
            ltype = 3
            call calcfg_l1(zbio, reac1, reac2, D1, ktemp, e_zbio_l1, dedz_zbio_l1, f_zbio_l1, dfdz_zbio_l1, g_zbio_l1, dgdz_zbio_l1)
            call calcfg_l2(zbio, reac1, reac2, D2, ktemp, e_zbio_l2, dedz_zbio_l2, f_zbio_l2, dfdz_zbio_l2, g_zbio_l2, dgdz_zbio_l2)
!             print*, ' '
!            print*, 'e_zbio_l2, dedz_zbio_l2, f_zbio_l2, dfdz_zbio_l2, g_zbio_l2, dgdz_zbio_l2 ', &
!                        e_zbio_l2, dedz_zbio_l2, f_zbio_l2, dfdz_zbio_l2, g_zbio_l2, dgdz_zbio_l2

            ! match solutions at zbio - continuous concentration and flux
            call matchsoln(e_zbio_l1, f_zbio_l1, g_zbio_l1, D1*dedz_zbio_l1, D1*dfdz_zbio_l1, D1*dgdz_zbio_l1, &
            e_zbio_l2, f_zbio_l2, g_zbio_l2, D2*dedz_zbio_l2, D2*dfdz_zbio_l2, D2*dgdz_zbio_l2, &
            0.0, 0.0, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f)
        !           print*, 'in prepfg_l12 AFTER matchsoln:  ls_a, ls_b, ls_c, ls_d, ls_e, ls_f', ls_a, ls_b, ls_c, ls_d, ls_e, ls_f
        end if

    END SUBROUTINE prepfg_l12

    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************

    !                               Utility subroutines for PO4/Fe-bound P

    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE calcfg_l12_PO4_M(z, reac1P, reac2P, dum_ktempP, dum_QtempP, &
    dum_D1P, dum_D2P, dum_alphaP, dum_mat_C, dum_vec_D, dum_ltype, &
    dum_ktempM, dum_QtempM, dum_D1M, dum_D2M, dum_alphaM, e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, &
    p_P, dpdz_P, q_P, dqdz_P, e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M)
            ! calculate solution basis functions, for layer which may cross bioturbation boundary
            !
            ! reac1, reac2        - mol/mol S released per organic carbon C
            !
            ! General solution for solute S is given by
            !  S(z) = A * e(z) + B * f(z) + g(z)
            !
            ! Where e,f,g are generated by matching solutions across bioturbation boundary (if necessary)
            ! Solution properties (matching etc) are input in ls
            ! On input, ls should contain fields generated by prepfg_l12

        real,intent(in):: z, reac1P, reac2P, dum_ktempP, dum_QtempP, dum_alphaP, dum_D1P, dum_D2P
        real,intent(in):: dum_ktempM, dum_QtempM, dum_alphaM, dum_D1M, dum_D2M
        INTEGER, INTENT(in):: dum_ltype
        real, dimension (4, 4), INTENT(in) :: dum_mat_C
        real, dimension (1:4), INTENT(in) ::  dum_vec_D

        ! ODE solutions (E, F, P, Q) and the particulat integral (G) and their derivatives
        real,intent(inout):: e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P
        real,intent(inout):: e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M

        ! local variables
        real loc_a1_P, loc_b1_P, loc_a1_M, loc_b1_M, loc_a2_M, loc_a2_P, loc_b2_P
        real, dimension (1:6) :: loc_Phi1_P, loc_Phi2_P
        real loc_e_P_1, loc_dedz_P_1, loc_f_P_1, loc_dfdz_P_1, loc_g_P_1, loc_dgdz_P_1
        real loc_p_P_1, loc_dpdz_P_1, loc_q_P_1, loc_dqdz_P_1
        real loc_e_M_1, loc_dedz_M_1, loc_f_M_1, loc_dfdz_M_1, loc_g_M_1, loc_dgdz_M_1
        real loc_p_M_1, loc_dpdz_M_1, loc_q_M_1, loc_dqdz_M_1
        ! save the ODE solutions in vectors to make calculation easier (DH?: however, is this faster?)
        real, dimension (1:4) :: loc_EFPQ_P, loc_dEFPQdz_P, loc_EFPQ_M, loc_dEFPQdz_M
        ! the transformed ODE solutions coming from xformsoln_PO4_M
        real, dimension (1:4) :: loc_EFPQ_P_t, loc_dEFPQdz_P_t, loc_EFPQ_M_t, loc_dEFPQdz_M_t

        loc_Phi1_P = (/ 0, 0, 0, 0, 0, 0 /)
        loc_Phi2_P = (/ 0, 0, 0, 0, 0, 0 /)

        select case (dum_ltype)
            case (1)    ! bioturbated
                if(dum_alphaP==0)then   ! oxic layer -> call PO4 first
                    call calcfg_l1_PO4(z, reac1P, reac2P, dum_D1P, dum_ktempP, dum_QtempP, 0.0, 0.0, dum_alphaP, &
                    e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, loc_a1_P, loc_b1_P, loc_Phi1_P)
                    call calcfg_l1_M(z, dum_D1M, dum_ktempM, dum_QtempM, loc_a1_P, loc_b1_P, loc_Phi1_P, dum_alphaM, &
                    e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M, loc_a1_M, loc_b1_M)
                else        ! anoxic layer -> call M first
                    call calcfg_l1_M(z, dum_D1M, dum_ktempM, dum_QtempM, 0.0, 0.0, loc_Phi1_P, dum_alphaM, &
                    e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M, loc_a1_M, loc_b1_M)
                    call calcfg_l1_PO4(z, reac1P, reac2P, dum_D1P, dum_ktempP, dum_QtempP, loc_a1_M, loc_b1_M, dum_alphaP, &
                    e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, loc_a1_P, loc_b1_P, loc_Phi1_P)
                end if

            case (2)    ! not bioturbated
                if(dum_alphaP==0)then   ! oxic layer -> call PO4 first
                    call calcfg_l2_PO4(z, reac1P, reac2P, dum_D2P, dum_ktempP, dum_QtempP, 0.0, dum_alphaP, &
                    e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, loc_a2_P, loc_b2_P, loc_Phi2_P)
                    call calcfg_l2_M(z, dum_ktempM, dum_QtempM, loc_a2_P, loc_b2_P, loc_Phi2_P, dum_alphaM, &
                    e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M, loc_a2_M)
                else    ! anoxic layer -> call M first
                    call calcfg_l2_M(z, dum_ktempM, dum_QtempM, 0.0, 0.0, loc_Phi2_P, dum_alphaM, &
                    e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M, loc_a2_M)
                    call calcfg_l2_PO4(z, reac1P, reac2P, dum_D2P, dum_ktempP, dum_QtempP, loc_a2_M, dum_alphaP, &
                    e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, loc_a2_P, loc_b2_P, loc_Phi2_P)
                end if

            case (3)    ! crossing boundary
                IF(z > zbio) THEN      ! not bioturbated region
                    if(dum_alphaP==0)then   ! oxic layer -> call PO4 first NOTE: BUT DECIDE VIA ALPHA_M NOT WITH <= ZOX!!! DOESN't WORK FOR BOUNDARY ZOX
                        call calcfg_l2_PO4(z, reac1P, reac2P, dum_D2P, dum_ktempP, dum_QtempP, 0.0, dum_alphaP, &
                        e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, loc_a2_P, loc_b2_P, loc_Phi2_P)
                        call calcfg_l2_M(z, dum_ktempM, dum_QtempM, loc_a2_P, loc_b2_P, loc_Phi2_P, dum_alphaM, &
                        e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M, loc_a2_M)
                    else ! anoxic layer -> call M first
                        call calcfg_l2_M(z, dum_ktempM, dum_QtempM, 0.0, 0.0, loc_Phi2_P, dum_alphaM, &
                        e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M, loc_a2_M)
                        call calcfg_l2_PO4(z, reac1P, reac2P, dum_D2P, dum_ktempP, dum_QtempP, loc_a2_M, dum_alphaP, &
                        e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, loc_a2_P, loc_b2_P, loc_Phi2_P)
                    end if ! (dum_alphaP==0)
                ELSE    ! bioturbated region z <= zbio
                    if(dum_alphaP==0)then   ! oxic layer -> call PO4 first
                        ! CASE 1 & 2: LAYER 1: have 4 int. const.
                        call calcfg_l1_PO4(z, reac1P, reac2P, dum_D1P, dum_ktempP, dum_QtempP, 0.0, 0.0, dum_alphaP, &
                        loc_e_P_1, loc_dedz_P_1, loc_f_P_1, loc_dfdz_P_1, loc_g_P_1, loc_dgdz_P_1, loc_p_P_1, loc_dpdz_P_1, &
                        loc_q_P_1, loc_dqdz_P_1, loc_a1_P, loc_b1_P, loc_Phi1_P)
                        call calcfg_l1_M(z, dum_D1M, dum_ktempM, dum_QtempM, loc_a1_P, loc_b1_P, loc_Phi1_P, dum_alphaM, &
                        loc_e_M_1, loc_dedz_M_1, loc_f_M_1, loc_dfdz_M_1, loc_g_M_1, loc_dgdz_M_1, loc_p_M_1, loc_dpdz_M_1, &
                        loc_q_M_1, loc_dqdz_M_1, loc_a1_M, loc_b1_M)
                        ! DH: FOR CASE 2: DON'T HAVE D FROM LAYER 2
                    else    ! anoxic layer -> call M first
                        ! DH: CASE 1: LAYER 2: have 4 int. const.
                        call calcfg_l1_M(z, dum_D1M, dum_ktempM, dum_QtempM, 0.0, 0.0, loc_Phi1_P, dum_alphaM, &
                        loc_e_M_1, loc_dedz_M_1, loc_f_M_1, loc_dfdz_M_1, loc_g_M_1, loc_dgdz_M_1, loc_p_M_1, loc_dpdz_M_1, &
                        loc_q_M_1, loc_dqdz_M_1, loc_a1_M, loc_b1_M)
                        call calcfg_l1_PO4(z, reac1P, reac2P, dum_D1P, dum_ktempP, dum_QtempP, loc_a1_M, loc_b1_M, dum_alphaP, &
                        loc_e_P_1, loc_dedz_P_1, loc_f_P_1, loc_dfdz_P_1, loc_g_P_1, loc_dgdz_P_1, loc_p_P_1, loc_dpdz_P_1, &
                        loc_q_P_1, loc_dqdz_P_1, loc_a1_P, loc_b1_P, loc_Phi1_P)
                    end if  ! (dum_alphaP==0)

                    ! Now find 'transformed' basis functions such that in layer 1,
                    ! O2 = A_2*et + B_2*ft + gt  (ie layer 1 solution written in terms of layer 2 coeffs A_2, B_2)

                    loc_EFPQ_P = (/ loc_e_P_1, loc_f_P_1, loc_p_P_1, loc_q_P_1 /)
                    loc_dEFPQdz_P = (/ loc_dedz_P_1, loc_dfdz_P_1, loc_dpdz_P_1, loc_dqdz_P_1 /)
                    loc_EFPQ_M = (/ loc_e_M_1, loc_f_M_1, loc_p_M_1, loc_q_M_1 /)
                    loc_dEFPQdz_M = (/loc_dedz_M_1, loc_dfdz_M_1, loc_dpdz_M_1, loc_dqdz_M_1/)

                    call xformsoln_PO4_M(loc_EFPQ_P, loc_EFPQ_M, loc_dEFPQdz_P, loc_dEFPQdz_M, &
                    loc_g_P_1, loc_g_M_1,loc_dgdz_P_1, loc_dgdz_M_1, dum_mat_C, dum_vec_D, &
                    loc_EFPQ_P_t, g_P, loc_dEFPQdz_P_t, dgdz_P, loc_EFPQ_M_t, &
                    g_M, loc_dEFPQdz_M_t, dgdz_M)

                    ! WHEN lTYPE=3 - DEAL WITH ONE VARIABLE SHORT FROM LAYER BELOW
                    ! FOR CASE 1 & 2: deal with missing values from layer below
                    if(zox .LE. zbio)then   ! CASE 1: no F from layer below - DH: no e_M & f_M as well !? but 0 anyway at the moment
                        e_P = loc_EFPQ_P_t(1)
                        f_P = loc_EFPQ_P_t(2)
                        p_P = loc_EFPQ_P_t(3)
                        q_P = loc_q_P_1

                        dedz_P = loc_dEFPQdz_P_t(1)
                        dfdz_P = loc_dEFPQdz_P_t(2)
                        dpdz_P = loc_dEFPQdz_P_t(3)
                        dqdz_P = loc_dqdz_P_1

                        e_M = loc_EFPQ_M_t(1)
                        f_M = loc_EFPQ_M_t(2)
                        p_M = loc_EFPQ_M_t(3)
                        q_M = loc_q_M_1

                        dedz_M = loc_dEFPQdz_M_t(1)
                        dfdz_M = loc_dEFPQdz_M_t(2)
                        dpdz_M = loc_dEFPQdz_M_t(3)
                        dqdz_M = loc_dqdz_M_1

                    else            ! CASE 2: no Q from layer below - DH: no p_P as well !?
                        e_P = loc_EFPQ_P_t(1)
                        f_P = loc_EFPQ_P_t(2)
                        p_P = loc_p_P_1             ! DH:  was = EFPQ_P_t(3);
                        q_P = loc_q_P_1

                        dedz_P = loc_dEFPQdz_P_t(1)
                        dfdz_P = loc_dEFPQdz_P_t(2)
                        dpdz_P = loc_dpdz_P_1       !DH: was = dEFPQ_P_t(3);
                        dqdz_P = loc_dqdz_P_1

                        e_M = loc_EFPQ_M_t(1)
                        f_M = loc_EFPQ_M_t(2)
                        p_M = loc_EFPQ_M_t(3)
                        q_M = loc_q_M_1

                        dedz_M = loc_dEFPQdz_M_t(1)
                        dfdz_M = loc_dEFPQdz_M_t(2)
                        dpdz_M = loc_dEFPQdz_M_t(3)
                        dqdz_M = loc_dqdz_M_1

                    end if ! (zox .LE. zbio)

                END IF  ! (z > zbio)
            case default
            print*, ' unrecognized ltype in  calcfg_l2_PO4 ', dum_ltype
                STOP
        end select

!        print*, ' '
!        print*, 'IN  calcfg_l12_PO4 --------'
!        print*, ' z = ', z
!        print*,'e_P, dedz_P ', char(9), e_P, dedz_P
!        print*,'f_P, dfdz_P', char(9), f_P, dfdz_P
!        print*,'g_P, dgdz_P', char(9), g_P, dgdz_P
!        print*,'p_P, dpdz_P', char(9), p_P, dpdz_P
!        print*,' q_P, dqdz_P ', char(9), q_P, dqdz_P
!        print*, ' '
!        print*,'e_M, dedz_M ', char(9), e_M, dedz_M
!        print*,'f_M, dfdz_M', char(9), f_M, dfdz_M
!        print*,'g_M, dgdz_M', char(9), g_M, dgdz_M
!        print*,'p_M, dpdz_M', char(9), p_M, dpdz_M
!        print*,' q_M, dqdz_M ', char(9), q_M, dqdz_M


    END SUBROUTINE calcfg_l12_PO4_M


    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------



    SUBROUTINE calcfg_l1_PO4(z, reac1, reac2, Dtemp, ktemp, Qtemp, a1_M, b1_M, alpha, e, dedz, f, dfdz, g, dgdz, &
    p, dpdz, q, dqdz, a1, b1, Phi1)
        ! Basis functions for solutes, case z <= zbio
        !
        ! reac1, reac2        - mol./mol S released per organic carbon C
        ! depend1,   depend2 coming from other species
        !
        ! General solution for solute S is given by
        !  S(z) = A * e(z) + B * f(z) + g(z)
        ! and for dependent species
        !       S(z) = A .* e(z) + B .* f(z) + C .* p(z) +  D.* q(z) + g(z)

        real z, reac1, reac2, Dtemp, ktemp, Qtemp, a1_M, b1_M, alpha              ! in from SUBROUTINE before
        real,INTENT(inout)::  e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, a1, b1 ! out
        real, dimension (1:6), intent (inout) :: Phi1                               ! out

        ! local variables
        real pfac
        real ea11z, eb11z, ea12z, eb12z


        a1=(w-sqrt(w**2+4.0*Dtemp*ktemp))/(2.0*Dtemp)
        e=exp(z*a1)
        dedz = a1*exp(z*a1)

        b1=(w+sqrt(w**2+4.0*Dtemp*ktemp))/(2.0*Dtemp)
        f=exp(z*b1)
        dfdz = b1*exp(z*b1);

        pfac = 1                    ! in fact, already has (1-por)/por

        ! NOW to OM reaction terms!
        ! save all Phis in one variable to pass back
        Phi1(1) = -pfac*k1*(reac1)*A11/(Dtemp*aa11**2-w*aa11-ktemp)
        Phi1(2)  = pfac*k1*(reac1)*A11/(Dtemp*bb11**2-w*bb11-ktemp)
        Phi1(3) =-pfac*k1*(reac1)*C01/(Dtemp*bb11**2-w*bb11-ktemp)
        Phi1(4)   =-pfac*k2*(reac2)*A12/(Dtemp*aa12**2-w*aa12-ktemp)
        Phi1(5)  = pfac*k2*(reac2)*A12/(Dtemp*bb12**2-w*bb12-ktemp)
        Phi1(6) =-pfac*k2*(reac2)*C02/(Dtemp*bb12**2-w*bb12-ktemp)


        ea11z = exp(aa11*z);
        eb11z = exp(bb11*z);
        ea12z = exp(aa12*z);
        eb12z = exp(bb12*z);

        if(ktemp==0)then        !CHECK: actually no need as ktemp always <> 0
            g =  Phi1(1)*ea11z + Phi1(2)*eb11z + Phi1(3)*eb11z + &
            Phi1(4)*ea12z + Phi1(5)*eb12z + Phi1(6)*eb12z
        else
            g =  Phi1(1)*ea11z + Phi1(2)*eb11z + Phi1(3)*eb11z + &
            Phi1(4)*ea12z + Phi1(5)*eb12z + Phi1(6)*eb12z + Qtemp/ktemp   ! here problem if ktemp=0
        end if

        dgdz = Phi1(1)*aa11*ea11z + Phi1(2)*bb11*eb11z + Phi1(3)*bb11*eb11z + &
        Phi1(4)*aa12*ea12z + Phi1(5)*bb12*eb12z + Phi1(6)*bb12*eb12z

        if(alpha == 0)then      ! was z<=res.zox PO4 is independent of M (no info in alpha)
            p = 0
            dpdz = 0
            q = 0
            dqdz = 0
        else                    ! PO4 is dependent on M
            p = -alpha/(Dtemp*a1_M**2-w*a1_M-ktemp)*exp(z*a1_M)
            dpdz = a1_M*p
            q = -alpha/(Dtemp*b1_M**2-w*b1_M-ktemp)*exp(z*b1_M)
            dqdz = b1_M*q
        end if

    !            print*, 'INPUT calcfg_l1', z, reac1, reac2, Dtemp, ktemp
    !            print*, 'IN  calcfg_l1 g', g
    !            print*, 'IN  calcfg_l1 dgdz', dgdz

    end SUBROUTINE calcfg_l1_PO4

    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE calcfg_l2_PO4(z, reac1, reac2, Dtemp, ktemp, Qtemp, a2_M, alpha, &
    e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, a2, b2, Phi2)
        ! Basis functions for solutes, case z > zbio
        ! reac1, reac2        - mol/mol S released per organic carbon C
        !
        ! General solution for solute S is given by
        !  S(z) = A * e(z) + B * f(z) + g(z)

        ! in from SUBROUTINE before
        real z, reac1, reac2, Dtemp, ktemp, Qtemp , a2_M, alpha
        ! out
        real,INTENT(inout)::  e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, a2, b2
        real, dimension (1:2), intent (inout) :: Phi2

        ! local variables
        real pfac
        real ea11z, eb11z, ea12z, eb12z


        a2=(w-sqrt(w**2+4.0*Dtemp*ktemp))/(2.0*Dtemp)
        e=exp(z*a2)
        dedz = a2*exp(z*a2)

        b2=(w+sqrt(w**2+4.0*Dtemp*ktemp))/(2.0*Dtemp)
        f=exp(z*b2)
        dfdz = b2*exp(z*b2)

        !pfac=1./por;   ! assume org matter already *(1-por)
        pfac = 1            !in fact, already has (1-por)/por

        Phi2(1) = -pfac*k1*(reac1)*A21/(Dtemp*aa21**2-w*aa21-ktemp)
        Phi2(2) = -pfac*k2*(reac2)*A22/(Dtemp*aa22**2-w*aa22-ktemp)


        if(ktemp==0)then            ! CHECK: think no need for this as always ktemp <> 0
            g = Phi2(1)*exp(aa21*z) + Phi2(2)*exp(aa22*z)
        else
            g = Phi2(1)*exp(aa21*z) + Phi2(2)*exp(aa22*z) + Qtemp/ktemp
        end if
        dgdz = Phi2(1)*aa21*exp(aa21*z) + Phi2(2)*aa22*exp(aa22*z)

        if(alpha==0)then            ! was z<=res.zox PO4 is independent of M (no info in alpha)
            p = 0
            dpdz = 0
            q = 0
            dqdz = 0
        else                        ! PO4 is dependent on M
            p = -alpha/(Dtemp*a2_M**2-w*a2_M-ktemp)*exp(a2_M*z)
            dpdz = a2_M*p
            q=0
            dqdz=0
        end if

    !            print*, 'IN  calcfg_l2 g', g
    !            print*, 'IN  calcfg_l1 dgdz', dgdz

    end SUBROUTINE calcfg_l2_PO4


    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE calcfg_l1_M(z, Dtemp, ktemp, Qtemp, a1_P, b1_P, Phi1_P, alpha, e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, c1, d1)
        ! Basis functions for solutes, case z <= zbio
        !
        ! reac1, reac2        - mol./mol S released per organic carbon C
        ! depend1,   depend2 coming from other species
        !
        ! General solution for solute S is given by
        !  S(z) = A * e(z) + B * f(z) + g(z)
        ! and for dependent species
        !       S(z) = A .* e(z) + B .* f(z) + C .* p(z) +  D.* q(z) + g(z)

        real z, Dtemp, ktemp, Qtemp, a1_P, b1_P, alpha                                ! in from SUBROUTINE before
        real, dimension (1:6), intent (in) :: Phi1_P                                            ! in from SUBROUTINE before
        real,INTENT(inout)::  e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, c1, d1             ! out

        ! local variables


        c1=(w-sqrt(w**2+4.0*Dtemp*ktemp))/(2.0*Dtemp)
        p=exp(z*c1)
        dpdz = c1*exp(z*c1)

        d1=(w+sqrt(w**2+4.0*Dtemp*ktemp))/(2.0*Dtemp)
        q=exp(z*d1)
        dqdz = d1*exp(z*d1);

        if(alpha .NE. 0)then      ! oxic layer: was z<=res.zox BUT problems at boundary. M is dependent on PO4
            c1=0
            d1=0
            e = -alpha/(Dtemp*a1_P**2-w*a1_P-ktemp)*exp(z*a1_P)
            dedz = a1_P*e
            f = -alpha/(Dtemp*b1_P**2-w*b1_P-ktemp)*exp(z*b1_P)
            dfdz = b1_P*f
            g = -alpha*(Phi1_P(1)/(Dtemp*aa11**2-w*aa11-ktemp)*exp(z*aa11) + Phi1_P(2)/(Dtemp*bb11**2-w*bb11-ktemp)*exp(z*bb11) + &
            Phi1_P(3)/(Dtemp*bb11**2-w*bb11-ktemp)*exp(z*bb11) + &
            Phi1_P(4)/(Dtemp*aa12**2-w*aa12-ktemp)*exp(z*aa12) + Phi1_P(5)/(Dtemp*bb12**2-w*bb12-ktemp)*exp(z*bb12) + &
            Phi1_P(6)/(Dtemp*bb12**2-w*bb12-ktemp)*exp(z*bb12))
            dgdz = -alpha*(Phi1_P(1)/(Dtemp*aa11**2-w*aa11-ktemp)*exp(z*aa11)*aa11 + Phi1_P(2)/(Dtemp*bb11**2-w*bb11-ktemp) &
            *exp(z*bb11)*bb11 + Phi1_P(3)/(Dtemp*bb11**2-w*bb11-ktemp)*exp(z*bb11)*bb11 + &
            Phi1_P(4)/(Dtemp*aa12**2-w*aa12-ktemp)*exp(z*aa12)*aa12 + Phi1_P(5)/(Dtemp*bb12**2-w*bb12-ktemp)*exp(z*bb12)*bb12 + &
            Phi1_P(6)/(Dtemp*bb12**2-w*bb12-ktemp)*exp(z*bb12)*bb12)

        else                    ! anoxic layer: M is independent of PO4 (no value in alpha!)
            g = Qtemp/ktemp
            dgdz = 0
            e = 0
            dedz = 0
            f = 0
            dfdz = 0
        end if

    !            print*, 'INPUT calcfg_l1_M', z, reac1, reac2, Dtemp, ktemp
    !            print*, 'IN  calcfg_l1 g', g
    !            print*, 'IN  calcfg_l1 dgdz', dgdz

    end SUBROUTINE calcfg_l1_M


    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE calcfg_l2_M(z, ktemp, Qtemp, a2_P, b2_P, Phi2_P, alpha, e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, c2)
        ! Basis functions for solutes, case z > zbio
        ! reac1, reac2        - mol/mol S released per organic carbon C
        !
        ! General solution for solute S is given by
        !  S(z) = A * e(z) + B * f(z) + g(z)

        real z, ktemp, Qtemp, a2_P, b2_P, alpha                                           ! in from SUBROUTINE before
        real,INTENT(inout)::  e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, c2             ! out
        real, dimension (1:2), intent (in) :: Phi2_P                                        ! out

        ! local variables

        c2=0

        if(alpha .NE. 0)then            ! M is dependent of PO4, was z<=res.zox
            e=alpha/(w*a2_P)*exp(z*a2_P)
            dedz = a2_P*e
            f=alpha/(w*b2_P)*exp(z*b2_P)
            dfdz = b2_P*f
            p = 1                       ! DH CHECK/TODO: integration constant just C
            dpdz = 0
            q=0
            dqdz=0
            g = alpha/(w)*(Phi2_P(1)/(aa21)*exp(aa21*z) + &
            Phi2_P(2)/(aa22)*exp(aa22*z))
            dgdz = alpha/(w)*(Phi2_P(1)*exp(aa21*z) + &
            Phi2_P(2)*exp(aa22*z))
        else                        ! M is independent of PO4 - z > res.zox
            c2=-ktemp/w
            p=exp(c2*z)
            dpdz = c2*exp(c2*z)
            q=0
            dqdz=0
            g = Qtemp/ktemp
            dgdz = 0
            e=0
            dedz=0
            f=0
            dfdz=0

        end if

    end SUBROUTINE calcfg_l2_M

    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE prepfg_l12_PO4_M(reac1, reac2, ktempP, QtempP, zU, zL, D1P, D2P, alphaP, &
    ktempM, QtempM, D1M, D2M, alphaM, loc_mat_C, loc_vec_D, ltype)
        real reac1, reac2, ktempP, QtempP, zU, zL, D1P, D2P, alphaP, ktempM, QtempM, D1M, D2M, alphaM
        real, dimension (4, 4), intent (inout) :: loc_mat_C
        real, dimension (1:4), intent (inout) ::  loc_vec_D
        INTEGER,INTENT(inout):: ltype

        ! local variables
        real e_zbio_l1_P, dedz_zbio_l1_P, f_zbio_l1_P, dfdz_zbio_l1_P, g_zbio_l1_P, dgdz_zbio_l1_P
        real p_zbio_l1_P, dpdz_zbio_l1_P, q_zbio_l1_P, dqdz_zbio_l1_P, a1_P, b1_P
        real e_zbio_l2_P, dedz_zbio_l2_P, f_zbio_l2_P, dfdz_zbio_l2_P, g_zbio_l2_P, dgdz_zbio_l2_P
        real p_zbio_l2_P, dpdz_zbio_l2_P, q_zbio_l2_P, dqdz_zbio_l2_P, a2_P, b2_P
        real e_zbio_l1_M, dedz_zbio_l1_M, f_zbio_l1_M, dfdz_zbio_l1_M, g_zbio_l1_M, dgdz_zbio_l1_M
        real p_zbio_l1_M, dpdz_zbio_l1_M, q_zbio_l1_M, dqdz_zbio_l1_M, a1_M, b1_M
        real e_zbio_l2_M, dedz_zbio_l2_M, f_zbio_l2_M, dfdz_zbio_l2_M, g_zbio_l2_M, dgdz_zbio_l2_M
        real p_zbio_l2_M, dpdz_zbio_l2_M, q_zbio_l2_M, dqdz_zbio_l2_M, a2_M
        real Vb, Fb
        real, dimension (1:6) :: Phi1_P, Phi2_P
        real, dimension (1:4,1:4) :: mat_X, mat_Y
        real, dimension (1:4) ::  vec_Z
        integer :: loc_dim, loc_i, loc_j

        Phi1_P = (/ 0, 0, 0, 0, 0, 0 /)
        Phi2_P = (/ 0, 0, 0, 0, 0, 0 /)

        if(zL <= zbio)then  ! wholly within bioturbated layer
            ltype = 1
        elseif(zU >= zbio) then ! wholly within non-bioturbated layer
            ltype = 2
        else             ! crossing boundary - sort out solution matching at zbio
            !            print*, 'IN  CROSS BOUNDARY CASE '
            ltype = 3

            if(zL <= zox)then       ! oxic layer -> call PO4 first
                call calcfg_l1_PO4(zbio, reac1, reac2, D1P, ktempP, QtempP, 0.0, 0.0, alphaP, &
                e_zbio_l1_P, dedz_zbio_l1_P, f_zbio_l1_P, dfdz_zbio_l1_P, g_zbio_l1_P, &
                dgdz_zbio_l1_P, p_zbio_l1_P, dpdz_zbio_l1_P, q_zbio_l1_P, dqdz_zbio_l1_P, a1_P, b1_P, Phi1_P)
                call calcfg_l2_PO4(zbio, reac1, reac2, D2P, ktempP, QtempP, 0.0, alphaP, &
                e_zbio_l2_P, dedz_zbio_l2_P, f_zbio_l2_P, dfdz_zbio_l2_P, g_zbio_l2_P, &
                dgdz_zbio_l2_P, p_zbio_l2_P, dpdz_zbio_l2_P, q_zbio_l2_P, dqdz_zbio_l2_P, a2_P, b2_P, Phi2_P)
                call calcfg_l1_M(zbio, D1M, ktempM, QtempM, a1_P, b1_P, Phi1_P, alphaM, &
                e_zbio_l1_M, dedz_zbio_l1_M, f_zbio_l1_M, dfdz_zbio_l1_M, g_zbio_l1_M, dgdz_zbio_l1_M, &
                p_zbio_l1_M, dpdz_zbio_l1_M, q_zbio_l1_M, dqdz_zbio_l1_M, a1_M, b1_M)
                call calcfg_l2_M(zbio, ktempM, QtempM, a2_P, b2_P, Phi2_P, alphaM, &
                e_zbio_l2_M, dedz_zbio_l2_M, f_zbio_l2_M, dfdz_zbio_l2_M, g_zbio_l2_M, dgdz_zbio_l2_M, &
                p_zbio_l2_M, dpdz_zbio_l2_M, q_zbio_l2_M, dqdz_zbio_l2_M, a2_M)

            else                ! anoxic layer -> call M first
                call calcfg_l1_M(zbio, D1M, ktempM, QtempM, 0.0, 0.0, Phi1_P, alphaM, &
                e_zbio_l1_M, dedz_zbio_l1_M, f_zbio_l1_M, dfdz_zbio_l1_M, g_zbio_l1_M, dgdz_zbio_l1_M, &
                p_zbio_l1_M, dpdz_zbio_l1_M, q_zbio_l1_M, dqdz_zbio_l1_M, a1_M, b1_M)
                call calcfg_l2_M(zbio, ktempM, QtempM, 0.0, 0.0, Phi2_P, alphaM, &
                e_zbio_l2_M, dedz_zbio_l2_M, f_zbio_l2_M, dfdz_zbio_l2_M, g_zbio_l2_M, dgdz_zbio_l2_M, &
                p_zbio_l2_M, dpdz_zbio_l2_M, q_zbio_l2_M, dqdz_zbio_l2_M, a2_M)
                call calcfg_l1_PO4(zbio, reac1, reac2, D1P, ktempP, QtempP, a1_M, b1_M, alphaP, &
                e_zbio_l1_P, dedz_zbio_l1_P, f_zbio_l1_P, dfdz_zbio_l1_P, g_zbio_l1_P, dgdz_zbio_l1_P, &
                p_zbio_l1_P, dpdz_zbio_l1_P, q_zbio_l1_P, dqdz_zbio_l1_P, a1_P, b1_P, Phi1_P)
                call calcfg_l2_PO4(zbio, reac1, reac2, D2P, ktempP, QtempP, a2_M, alphaP, &
                e_zbio_l2_P, dedz_zbio_l2_P, f_zbio_l2_P, dfdz_zbio_l2_P, g_zbio_l2_P, dgdz_zbio_l2_P, &
                p_zbio_l2_P, dpdz_zbio_l2_P, q_zbio_l2_P, dqdz_zbio_l2_P, a2_P, b2_P, Phi2_P)
            end if

            ! match solutions at zbio - continuous concentration and flux
            ! organize the data in matrices, and use the intrinsic fortran fct.
            ! DH: Maybe more efficient when written out !?

            !  |x1        |   | A_l |      | y1        | | A_r|    |z1|    always PO4 continuity
            !  |    .     |   | B_l |      |    .      | | B_r|    |z2|    always PO4 flux
            !  |      .   |   | C_l |   =  |      .    | | C_r|  + |z3|    always M continuity
            !  |       x16|   | D_l |      |        y16| | D_r|    |z4|    SD always M _diffusive_ flux  = 0 (cf org C)

            ! discontinuity constants
            Vb = 0
            Fb = 0

            ! weird FORTRAN matrices makes the transpose necessary
            ! matrix mat_X
            mat_X = transpose(reshape((/ e_zbio_l1_P, f_zbio_l1_P, p_zbio_l1_P, q_zbio_l1_P, &
            D1P*dedz_zbio_l1_P, D1P*dfdz_zbio_l1_P, D1P*dpdz_zbio_l1_P, D1P*dqdz_zbio_l1_P, &
            e_zbio_l1_M, f_zbio_l1_M, p_zbio_l1_M, q_zbio_l1_M, &
            D1M*dedz_zbio_l1_M, D1M*dfdz_zbio_l1_M, D1M*dpdz_zbio_l1_M, D1M*dqdz_zbio_l1_M/), shape(mat_X)))
            !            data (mat_X(1,loc_i), loc_i=1,4) /  e_zbio_l1_P, f_zbio_l1_P, p_zbio_l1_P, q_zbio_l1_P /
            !            data (mat_X(2,loc_i), loc_i=1,4) /  D1P*dedz_zbio_l1_P, D1P*dfdz_zbio_l1_P, D1P*dpdz_zbio_l1_P, D1P*dqdz_zbio_l1_P /
            !            data (mat_X(3,loc_i), loc_i=1,4) /  e_zbio_l1_M, f_zbio_l1_M, p_zbio_l1_M, q_zbio_l1_M /
            !            data (mat_X(4,loc_i), loc_i=1,4) /  D1M*dedz_zbio_l1_M, D1M*dfdz_zbio_l1_M, D1M*dpdz_zbio_l1_M, D1M*dqdz_zbio_l1_M /

            ! matrix mat_Y
            mat_Y = transpose(reshape((/ e_zbio_l2_P, f_zbio_l2_P, p_zbio_l2_P, q_zbio_l2_P, &
            D2P*dedz_zbio_l2_P, D2P*dfdz_zbio_l2_P, D2P*dpdz_zbio_l2_P, D2P*dqdz_zbio_l2_P, &
            e_zbio_l2_M, f_zbio_l2_M, p_zbio_l2_M, q_zbio_l2_M, &
            D2M*dedz_zbio_l2_M, D2M*dfdz_zbio_l2_M, D2M*dpdz_zbio_l2_M, D2M*dqdz_zbio_l2_M/), shape(mat_Y)))
            !            data (mat_Y(1,loc_i), loc_i=1,4) /  e_zbio_l2_P, f_zbio_l2_P, p_zbio_l2_P, q_zbio_l2_P /
            !            data (mat_Y(2,loc_i), loc_i=1,4) /  D2P*dedz_zbio_l2_P, D2P*dfdz_zbio_l2_P, D2P*dpdz_zbio_l2_P, D2P*dqdz_zbio_l2_P /
            !            data (mat_Y(3,loc_i), loc_i=1,4) /  e_zbio_l2_M, f_zbio_l2_M, p_zbio_l2_M, q_zbio_l2_M /
            !            data (mat_Y(4,loc_i), loc_i=1,4) /  D2M*dedz_zbio_l2_M, D2M*dfdz_zbio_l2_M, D2M*dpdz_zbio_l2_M, D2M*dqdz_zbio_l2_M /

            vec_Z = (/ g_zbio_l2_P-g_zbio_l1_P + Vb, &
            D2P*dgdz_zbio_l2_P - D1P*dgdz_zbio_l1_P + Fb - w*Vb, &
            g_zbio_l2_M-g_zbio_l1_M + Vb, &
            D2M*dgdz_zbio_l2_M - D1M*dgdz_zbio_l1_M + Fb - w*Vb /)


!            201         format (6f12.6)
!                        print*,'mat_X '
!                        do loc_i=1,4
!                            write (*,201) (mat_X(loc_i,loc_j),loc_j=1,4)
!                        end do
!                        print*,'mat_Y '
!                        do loc_i=1,4
!                            write (*,201) (mat_Y(loc_i,loc_j),loc_j=1,4)
!                        end do
!                        print*,'vec_Z ', char(9), vec_Z
!                        print*,' '

            loc_dim = 4
            call matchsoln_PO4_M(mat_X, mat_Y, vec_Z, loc_dim, loc_mat_C, loc_vec_D)
        end if

    END SUBROUTINE prepfg_l12_PO4_M


    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************

    !                               Oxygen

    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE zO2()

        real flxzox, conczox, flxswi, fun0, fun1, zL, tol
        integer bctype
        real FO2

!            print*, ''
!            print*, '---------------------- START zO2 ------------------------ '

        zox = 1e-10
        bctype = 1
        flxzox = 0.0
        conczox = 0.0
        flxswi = 0.0

        !    call zO2_calcbc(zox, bctype, flxzox, conczox, flxswi, r_zxf)

        !    fun0 = flxzox + calcFO2(zox)
        !    print*,'!!!!!!!!!!!! fun0', fun0

        fun0 = FUN_zO2(zox)
        !    print*,' fun0', fun0

        !    print*,' '
        !    print*,'Try zero flux at zinf and see if we have any O2 left'
        ! Try zero flux at zinf and see if we have any O2 left
        call zO2_calcbc(zinf, 2, flxzox, conczox, flxswi, r_zxf);
!            print*,'flxzox', flxzox
!            print*,'conczox', conczox
!            print*,'flxswi', flxswi

        if (fun0 .ge. 0)then   ! eg zero oxygen at swi
            zox = 0
            bctype = 1
        elseif (conczox .ge. 0)then      ! still O2 at zinf -> zox = zinf
            zox = zinf
            bctype = 2
        else                        ! search zox in the interval
            bctype = 1
            zL=1e-10
            tol=1e-16
            zox = zbrent(FUN_zO2, zL, zinf, tol)
        !        print*,'$$$$$$$$$$$$$4   CALCULATE zox = ', zox
        !        zox = 2.7134
        end if

        call zO2_calcbc(zox, bctype, flxzox, conczox, flxswi, r_zxf)
        !    print*,' '
        !    print*,'---------- FINAL RESULTS zO2 --------- '
        print*,'zox ', char(9), zox
        print*,'r_zxf', char(9), r_zxf
        print*,''
        print*,'flxzox', char(9), flxzox
        print*,'conczox', char(9), conczox
        print*,'flxswiO2', char(9), flxswi


    END SUBROUTINE zO2

    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE zO2_calcbc(zox, bctype, flxzox, conczox, flxswi,r_zxf)

        !   Solve O2
        real, intent(in)::zox
        integer, intent(in)::bctype
        real, intent(inout)::flxzox, conczox, flxswi, r_zxf
        !   local variables

        !    real qdispO2                          !O2 diffusion coefficient in water (cm2/yr)
        !    real adispO2                          !O2 linear coefficient for temperature dependence (cm2/yr/oC)
        !    real DO21                             !O2 diffusion coefficient in bioturbated layer (cm2/yr)
        !    real DO22                             !O2 diffusion coefficient in non-bioturbated layer (cm2/yr)
        real reac1, reac2, ls !, z0, zox
        integer ltype
        real ls_a, ls_b, ls_c, ls_d, ls_e, ls_f
        real e_0, dedz_0, f_0, dfdz_0, g_0, dgdz_0
        real e_zox, dedz_zox, f_zox, dfdz_zox, g_zox, dgdz_zox
        !    real bctype1_AO2, bctype1_BO2, bctype2_AO2, bctype2_BO2

        real rO2_AO2, rO2_BO2, Dzox, Dswi

        !    real FO2

        !    qdispO2=348.5750
        !    adispO2=14.0890

        !    DO21=(qdispO2+adispO2*T)*dispFactor+Dbio
        !    DO22=(qdispO2+adispO2*T)*dispFactor


        !   reactive terms: OM degradation (-) and nitrification (-)
        reac1=-OC-2*gamma*NC1
        reac2=-OC-2*gamma*NC2

!            print*, ''
!            print*, '------- START zO2_calcbc ----- zox:', zox

        ! calculate solution for given zox

        !    print*, 'Preparation: sort out solution-matching across bioturbation boundary (if necessary)'
        ! Preparation: sort out solution-matching across bioturbation boundary (if necessary)
        call prepfg_l12(reac1, reac2, 0.0, z0, zox, DO21, DO22, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ltype)
!        print*, ''
!        print*, ' result prepfg_l12: ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ltype', ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ltype
!        print*, 'e_0, dedz_0, f_0, dfdz_0, g_0, dgdz_0:', e_0, dedz_0, f_0, dfdz_0, g_0, dgdz_0

        ! basis functions at upper boundary
        call calcfg_l12(z0, reac1, reac2, 0.0, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, DO21, DO22, ltype, &
        e_0, dedz_0, f_0, dfdz_0, g_0, dgdz_0)


        ! ... and lower boundary
        call calcfg_l12(zox, reac1, reac2, 0.0, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, DO21, DO22, ltype, e_zox, dedz_zox,&
        f_zox, dfdz_zox, g_zox, dgdz_zox)

        ! Solve for AO2, BO2 given boundary conditions (expressed in terms of transformed soln)

        IF(bctype==1)THEN
            ! Case 1 zero concentration at zox
            ! AO2*e_zox   +   BO2*f_zox  + g_zox = 0;
            ! AO2*e_0     +   BO2*f_0     + g_0  = swi.O20;

            ! | e_zox f_zox |  |AO2|   = | -g_zox         |
            ! | e_0   f_0   |  |BO2|     | swi.O20 - gt_0 |

            call solve2eqn(e_zox, f_zox, e_0, f_0, -g_zox, O20 - g_0, rO2_AO2, rO2_BO2)
        ELSE
            ! Case  2 zero flux at zox
            ! AO2*dedz_zox +  BO2*dfz_zox + dgz_zox = 0;
            ! AO2*e_0     +   BO2*f_0     + g_0     = swi.O20;
                         ! a            b       c   d       e           f
            call solve2eqn(dedz_zox, dfdz_zox, e_0, f_0, -dgdz_zox, O20 - g_0, rO2_AO2, rO2_BO2)
        END IF

        IF(zox < zbio)THEN
            Dzox = DO21
        ELSE
            Dzox = DO22
        END IF

        flxzox =  Dzox*(rO2_AO2*dedz_zox + rO2_BO2*dfdz_zox + dgdz_zox)         ! no por factor as this is per cm^2 pore area
        conczox = rO2_AO2*e_zox + rO2_BO2 * f_zox + g_zox


        ! TODO: ASK STUART or ANDY
        if(flxzox /= flxzox)then !check for NaN if then give value as in matlab.....
!            print*,' '
!            print*,' '
!            print*,'------ zO2_calcbc --------- flxzox is INFFFFFFFFFFFFFFFF at zox', flxzox, zox
!            print*,' '
!            print*,' '
            flxzox = -1.2e6
        end if

        if(0 < zbio)then
            Dswi = DO21
        else
            Dswi = DO22
        end if

        flxswi = por*(Dswi*(rO2_AO2*dedz_0+rO2_BO2*dfdz_0 + dgdz_0) - w*O20)   ! por fac so this is per cm^2 water column

        r_zxf = zox/(zoxgf + zox)   ! roll off oxidation at low zox
    !    FO2 = calcFO2()

    !    print*,'END zO2_calcbc: flxzox', flxzox
    !    print*,'conczox', conczox
    !    print*,'rO2_AO2, e_zox, rO2_BO2, f_zox, g_zox', rO2_AO2, e_zox, rO2_BO2, f_zox, g_zox

    end SUBROUTINE zO2_calcbc

    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    FUNCTION calcFO2(z)

        real calcFO2, z, tmpreac1, tmpreac2

        ! Oxydation of reduced species at zox (NEED A RATIO for ODU! and add NH4
        ! adsporption!
        !tmpreac1=0.5*OC+2*gamma*NC1;
        !tmpreac2=0.5*OC+2*gamma*NC2;

        !    print*,' '
        !    print*,'..... START calcFO2'
        tmpreac1=gammaH2S*1.0*SO4C+2.0*gamma*NC1;
        tmpreac2=gammaH2S*1.0*SO4C+2.0*gamma*NC2;

        !FLUX of NH4 and Reduced species from ZOX to ZINF


        calcFO2 = z/(zoxgf + z) * calcReac(z, zinf, tmpreac1, tmpreac2)

    !    print*,'calcFO2', calcFO2

    END FUNCTION calcFO2

    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************

    !                       Nitrate

    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE benthic_zNO3()

        real flxzinf, conczinf, flxswi, fun0, flxzno3, conczno3, flxswiNO3, zL, tol
        integer bctype
        real FNO3

        !    print*, ''
        !    print*, '------------------------------------------------------------------'
        !    print*, '---------------------- START zNO3 ------------------------------- '

        IF(zox == zinf)THEN
            zno3 = zinf
            bctype = 2
        ELSE

            bctype = 2
            ! Try zero flux at zinf and see if we have any NO3 left
            !    print*,'Try zero flux at zinf and see if we have any NO3 left'
            call zNO3_calcbc(zinf, bctype, flxzinf, conczinf, flxswi)

            !    print*,'RESULTS Try zero flux at zinf zNO3_calcbc flxzinf, conczinf, flxswi', flxzinf, conczinf, flxswi

            IF (conczinf > 0.0) THEN
                zno3 = zinf
                bctype = 2;
            ELSE
                zL=1e-10
                tol=1e-16
                zno3 = zbrent(FUN_zNO3, max(zL,zox), zinf, tol)
                !        print*,'$$$$ calculated zno3 =', zno3
                !        zno3 = 7.4319         ! use qualifier d0 fuer double: 7.4319d0
                bctype = 1;
            END IF
        END IF
        call zNO3_calcbc(zno3, bctype, flxzno3, conczno3, flxswiNO3)

        !    print*,' ---- FINAL RESULTS zNO3: ----'
        print*,'zno3', char(9), zno3
        print*,''
        print*,'flxzno3', char(9), flxzno3
        print*,'conczno3', char(9), conczno3
        print*,'flxswiNO3', char(9), flxswiNO3

    !            r.flxzno3 = flxzno3;  % should be zero
    !            r.conczno3 = conczno3; % may be non-zero if eg fully oxic sediment
    !            r.flxswiNO3 = flxswiNO3;
    END SUBROUTINE benthic_zNO3


    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE zNO3_calcbc(zNO3, bctype, flxzno3, conczno3, flxswi)

        real zNO3, flxzNO3, conczNO3, flxswi
        integer bctype, ltype1, ltype2
        real ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1
        real ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2
        real e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3
        real e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox
        real e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox
        real zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
        real e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00
        real e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0
        real bctype1_A2, bctype1_B2, bctype2_A2, bctype2_B2
        real rNO3_A2, rNO3_B2, rNO3_A1, rNO3_B1

        real FNH4

        ! Calculate trial solution for given zno3, matching boundary conditions from layer-by-layer solutions


        ! Preparation: for each layer, sort out solution - matching across bioturbation boundary (if necessary)
        ! layer 1: 0 < z < zox, nitrification
        !    prepfg_l12(reac1, reac2,  ktemp ,     zU , zL , D1,  D2)
        call prepfg_l12(gamma*NC1, gamma*NC2,  0.0,  0.0, zox, DN1, DN2, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, ltype1)
        !    print*, ''
        !    print*, '1. prepfg_l12 RESULTS: ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, ltype1', ls_a1, ls_b1, ls_c1, &
        !           & ls_d1, ls_e1, ls_f1, ltype1
        ! layer 2: zox < z < zno3, denitrification
        call prepfg_l12(-NO3CR, -NO3CR, 0.0,  zox, zNO3, DN1, DN2, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, ltype2)
        !    print*, ''
        !    print*, '2. prepfg_l12 RESULTS: ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, ltype2', ls_a2, ls_b2, ls_c2, ls_d2, &
        !            & ls_e2, ls_f2, ltype2


        ! Work up from the bottom, matching solutions at boundaries
        ! Basis functions at bottom of layer 2 zno3
        call calcfg_l12(zno3, -NO3CR, -NO3CR, 0.0, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DN1, DN2, ltype2, &
        e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3)

        ! Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from NH4 -> NO3 source)

        ! basis functions at bottom of layer 1
        call calcfg_l12(zox, gamma*NC1, gamma*NC2, 0.0, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DN1, DN2, ltype1, &
        e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox)

        ! basis functions at top of layer 2
        call calcfg_l12(zox, -NO3CR, -NO3CR, 0.0, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DN1, DN2, ltype2, &
        e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox)

        ! flux of NH4 to zox  TODO NH4 production by denitrification?

        FNH4 = calcReac(zno3, zinf, NC1/(1.0+KNH4), NC2/(1.0+KNH4))   ! MULTIPLY BY 1/POR ????
        !    print*, 'FNH4 = ', FNH4

        ! match solutions at zox - continuous concentration, flux discontinuity from H2S ox

        !    call matchsoln(e_zbio_l1, f_zbio_l1, g_zbio_l1, D1*dedz_zbio_l1, D1*dfdz_zbio_l1, D1*dgdz_zbio_l1, &
        !                         & e_zbio_l2, f_zbio_l2, g_zbio_l2, D2*dedz_zbio_l2, D2*dfdz_zbio_l2, D2*dgdz_zbio_l2, &
        !                                0.0, 0.0, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f)
        IF(zox .le. zbio)then
            call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
            e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
            0.0, -r_zxf*gamma*FNH4/DN1, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
        !        print*, ''
        !        print*,'zox<= zbio: RESULTS:  zox_a, zox_b, zox_c, zox_d, zox_e, zox_f', zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
        ELSE
            call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
            e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
            0.0, -r_zxf*gamma*FNH4/DN2, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
        !        print*, ''
        !        print*,'zox> zbio: RESULTS:  zox_a, zox_b, zox_c, zox_d, zox_e, zox_f', zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
        END IF

        ! Solution at swi, top of layer 1
        call calcfg_l12(0.0, gamma*NC1, gamma*NC2, 0.0, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DN1, DN2, ltype1, &
        e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00)
        !    print*, ''
        !    print*,'top of layer 1; calcfg_l12: RESULTS:  e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00', e1_00, &
        !            & dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00

        ! transform to use coeffs from l2
        call xformsoln(e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00, zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f, &
        e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0)
        !    print*, ''
        !    print*,'transform to use coeffs from l2; xformsoln: RESULTS:  e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0', e1_0, &
        !            & f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0

        ! Solve for ANO3, BNO3 given boundary conditions (expressed in terms of transformed basis fns, layer 2 A, B)

        ! Case 1 zero concentration at zno3
        ! ANO3*e2_zno3   +  BNO3*f2_zno3  + g2_zno3 = 0;
        ! ANO3*e1_0     +   BNO3*f1_0     + g1_0  = swi.NO30;

        ! | e2_zno3 f2_zno3 |  |ANO3|   = | -g2_zno3       |
        ! | e1_0     f1_0   |  |BNO3|     | swi.NO30 - g1_0 |

        call solve2eqn(e2_zno3, f2_zno3, e1_0, f1_0, -g2_zno3, NO30 - g1_0, bctype1_A2, bctype1_B2)

        ! Case  2 zero flux at zno3
        ! ANO3*de2dz_zno3   +  BNO3*dfdz2_zno3  + dgdz2_zno3 = 0;
        ! ANO3*e1_0         +   BNO3*f1_0       + g1_0       = swi.NO30;

        call solve2eqn(dedz2_zno3, dfdz2_zno3, e1_0, f1_0, -dgdz2_zno3, NO30 - g1_0, bctype2_A2, bctype2_B2)

        ! Choose type of solution requested
        IF(bctype==1) THEN
            rNO3_A2 = bctype1_A2
            rNO3_B2 = bctype1_B2
        ELSE
            rNO3_A2 = bctype2_A2
            rNO3_B2 = bctype2_B2
        END IF

        ! calculate flux at zno3
        IF(zno3 .le. zbio) THEN
            flxzno3 = DN1*(rNO3_A2*dedz2_zno3+rNO3_B2*dfdz2_zno3 + dgdz2_zno3)      ! includes 1/por ie flux per (cm^2 pore area)
        ELSE
            flxzno3 = DN2*(rNO3_A2*dedz2_zno3+rNO3_B2*dfdz2_zno3 + dgdz2_zno3)      ! includes 1/por ie flux per (cm^2 pore area)
        END IF

        conczno3 = rNO3_A2*e2_zno3+rNO3_B2*f2_zno3 + g2_zno3
        ! flux at swi - DO include por so this is per cm^2 water column area
        flxswi = por*(DN1*(rNO3_A2*dedz1_0+rNO3_B2*dfdz1_0 + dgdz1_0) - w*NO30)              ! NB: use A2, B2 as these are _xformed_ layer 1 basis functions

        !    print*,'flxzno3', flxzno3
        !    print*,'conczno3', conczno3
        !    print*,'flxswi', flxswi

        ! save coeffs for layer 1
        rNO3_A1 = zox_a*rNO3_A2 + zox_b*rNO3_B2 + zox_e
        rNO3_B1 = zox_c*rNO3_A2 + zox_d*rNO3_B2 + zox_f

    !    print*,'rNO3_A1, rNO3_B1', rNO3_A1, rNO3_B1

    END SUBROUTINE zNO3_calcbc


    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************

    !                           Sulfate

    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE benthic_zSO4()
        !    real, intent(in)::zox, zno3
        !    real, intent(inout)::zso4

        ! local variables
        real flxzso4, conczso4, flxswiso4, zL, tol
        integer bctype


        !    print*, ''
        !    print*, '------------------------------------------------------------------'
        !    print*, '---------------------- START zSO4 ------------------------------- '


        ! Iteratively solve for zso4

        IF(zno3 == zinf)THEN
            zso4 = zinf
            bctype = 2
        ELSE

            ! try zero flux at zinf and see if we have any SO4 left
            !    print*, ''
            !    print*, '-----try zero flux at zinf and see if we have any SO4 left------'
            bctype = 2

            call zSO4_calcbc(zinf, bctype, flxzso4, conczso4, flxswiso4)

            IF(conczso4 .ge.0)THEN
                zso4 = zinf
                bctype = 2
            ELSE
                bctype = 1
                zL=1e-10
                tol=1e-16
                zso4 = zbrent(FUN_zSO4, max(zno3,zL), zinf, tol)
            !        print*,'$$$$$$$$$$$$$4   CALCULATE zso4 = ', zso4
            END IF
        !    print*,'bctype, zso4 ', bctype, zso4
        END IF
        call zSO4_calcbc(zso4, bctype, flxzso4, conczso4, flxswiso4)

        !    print*,' '
        !    print*,'-------------------------------- FINAL RESULTS zSO4 --------------------------------'
        print*,'zso4', char(9), zso4
        print*,' '
        print*,'flxzso4', char(9), flxzso4
        print*,'conczso4', char(9), conczso4
        print*,'flxswiso4', char(9), flxswiso4

    END SUBROUTINE benthic_zSO4

    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE zSO4_calcbc(zso4, bctype, flxzso4, conczso4, flxswi)

        real, intent(in)::zso4
        integer, intent(in)::bctype
        real, intent(inout)::flxzso4, conczso4, flxswi !, r_zxf

        ! local variable
        integer ltype1, ltype2, ltype3
        real reac1_so4, reac2_so4
        real ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1
        real ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2
        real ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3
        real e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4
        real e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3
        real e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3
        real zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f
        real e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox
        real e2_zox0, f2_zox0, g2_zox0, dedz2_zox0, dfdz2_zox0, dgdz2_zox0
        real e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox
        real zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
        real e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0
        real e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00
        real bctype1_A3, bctype1_B3, bctype2_A3, bctype2_B3
        real rSO4_A1, rSO4_B1, rSO4_A2, rSO4_B2, rSO4_A3, rSO4_B3

        real FH2S

        !    print*, ' '
        !    print*, '---------------------- START zSO4_calcbc ------------------------------- '

        reac1_so4=-SO4C
        reac2_so4=-SO4C

        ! Calculate trial solution for given zso4, matching boundary conditions from layer-by-layer solutions


        ! Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
        ! layer 1: 0 < z < zox, passive diffn
        !      ls =      prepfg_l12( bsd, swi, r, reac1,     reac2,     ktemp, zU, zL, D1,        D2)

        call prepfg_l12(0.0, 0.0, 0.0, 0.0, zox, DSO41, DSO42, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, ltype1)

        ! layer 2: zox < z < zno3, passive diffn
        call prepfg_l12(0.0, 0.0, 0.0, zox, zno3, DSO41, DSO42, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, ltype2)

        ! layer 3: zno3 < z < zso4, SO4 consumption by OM oxidation
                !rSO4.ls3 = r.zTOC.prepfg_l12(bsd, swi, r, obj.reac1, obj.reac2, 0, r.zno3, zso4, obj.DSO41, obj.DSO42)
        call prepfg_l12(reac1_so4, reac2_so4, 0.0, zno3, zso4, DSO41, DSO42, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, ltype3)


        ! Work up from the bottom, matching solutions at boundaries
        ! Basis functions at bottom of layer 3 zso4

        call calcfg_l12(zso4, reac1_so4, reac2_so4, 0.0, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, DSO41, DSO42, ltype3, &
        e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4)

        ! Match at zno3, layer 2 - layer 3 (continuity and flux)
        ! basis functions at bottom of layer 2
        call calcfg_l12(zno3, 0.0, 0.0, 0.0, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DSO41, DSO42, ltype2, &
        e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3)

        ! ... and top of layer 3
        call calcfg_l12(zno3, reac1_so4, reac2_so4, 0.0, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, DSO41, DSO42, ltype3, &
        e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3)

        ! match solutions at zno3 - continuous concentration and flux
        call matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, &
        e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, &
        0.0, 0.0, zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f)


        ! Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from H2S source)
        ! flux of H2S to oxic interface (Source of SO4)
        ! NB: include methane region as AOM will produce sulphide as well..

        FH2S = calcReac(zno3, zso4, SO4C, SO4C) & ! MULTIPLY BY 1/POR ????
        + gammaCH4*calcReac(zso4, zinf, MC, MC)

        !    print*,' '
        !    print*,'FH2S ', FH2S

        ! basis functions at bottom of layer 1
        call calcfg_l12(zox, 0.0, 0.0, 0.0, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DSO41, DSO42, ltype1, &
        e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox)

        ! basis functions at top of layer 2
        call calcfg_l12(zox, 0.0, 0.0, 0.0, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DSO41, DSO42, ltype2, &
        e2_zox0, dedz2_zox0, f2_zox0, dfdz2_zox0, g2_zox0, dgdz2_zox0)

        ! transform to use coeffs from l3
        call xformsoln(e2_zox0, f2_zox0, g2_zox0, dedz2_zox0, dfdz2_zox0, dgdz2_zox0, &
        zno3_a , zno3_b , zno3_c , zno3_d , zno3_e ,zno3_f, &
        e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox )

        ! match solutions at zox - continuous concentration, flux discontinuity from H2S ox
        IF(zox .le. zbio)THEN
            call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
            e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
            0.0, -r_zxf*gammaH2S*FH2S/DSO41, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
        ELSE
            call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
            e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
            0.0, -r_zxf*gammaH2S*FH2S/DSO42, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
        END IF

        ! Solution at swi, top of layer 1
        call calcfg_l12(0.0, 0.0, 0.0, 0.0, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DSO41, DSO42, ltype1, &
        e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00)

        ! transform to use coeffs from l3
        call xformsoln(e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00, zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f, &
        e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0)

        ! Find solutions for two possible types of lower bc
        !  case 1  zero concentration at zso4
        ! Solve for ASO4, BSO4 given boundary conditions (expressed in terms of transformed basis fns, layer 3 A, B)
        ! ASO4*e3_zso4   +  BSO4*f3_zso4  + g3_zso4 = 0;
        ! ASO4*e1_0     +   BSO4*f1_0     + g1_0  = swi.SO40

        ! | e3_zso4 f3_zso4 |  |ASO4|   = | -g3_zso4       |
        ! | e1_0     f1_0   |  |BSO4|     | swi.SO40 - g1_0 |

        call solve2eqn(e3_zso4, f3_zso4, e1_0, f1_0, -g3_zso4, SO40 - g1_0, bctype1_A3, bctype1_B3)

        ! case  2 zero flux at zso4
        ! ASO4*de3dz_zso4   +  BSO4*dfdz3_zso4  + dgdz3_zso4 = 0;
        ! ASO4*e1_0         +   BSO4*f1_0       + g1_0       = swi.SO40;
        call solve2eqn(dedz3_zso4, dfdz3_zso4, e1_0, f1_0, -dgdz3_zso4, SO40 - g1_0, bctype2_A3, bctype2_B3)

        ! Choose type of solution requested
        !rSO4.A3 = (bctype==1).*bctype1_A3 + (bctype==2).*bctype2_A3;
        !        rSO4.B3 = (bctype==1).*bctype1_B3 + (bctype==2).*bctype2_B3;
        IF(bctype==1)THEN
            rSO4_A3=bctype1_A3
            rSO4_B3=bctype1_B3
        ELSE
            rSO4_A3=bctype2_A3
            rSO4_B3=bctype2_B3
        END IF

        ! calculate conc and flux at zso4
        conczso4 = rSO4_A3*e3_zso4+rSO4_B3*f3_zso4 + g3_zso4
        !D = (zso4 <= zbio).*obj.DSO41 + (zso4 > zbio).*obj.DSO42
        IF(zso4 .le. zbio)THEN
            flxzso4 = DSO41*(rSO4_A3*dedz3_zso4+rSO4_B3*dfdz3_zso4 + dgdz3_zso4)        ! includes 1/por ie flux per (cm^2 pore area)
        ELSE
            flxzso4 = DSO42*(rSO4_A3*dedz3_zso4+rSO4_B3*dfdz3_zso4 + dgdz3_zso4)
        END IF

        ! flux at swi - DO include por so this is per cm^2 water column area
        flxswi = por*(DSO41*(rSO4_A3*dedz1_0+rSO4_B3*dfdz1_0 + dgdz1_0) - w*SO40)   ! NB: use A3, B3 as these are _xformed_ layer 1 basis functions

        !    print*,' '
        !    print*,'RESULTS zso4_calcbc_: conczso4, flxzso4, flxswi', conczso4, flxzso4, flxswi

        ! save coeffs for layers 2 and 1
        rSO4_A2 = zno3_a*rSO4_A3 + zno3_b*rSO4_B3 + zno3_e
        rSO4_B2 = zno3_c*rSO4_A3 + zno3_d*rSO4_B3 + zno3_f

        rSO4_A1 = zox_a*rSO4_A3 + zox_b*rSO4_B3 + zox_e
        rSO4_B1 = zox_c*rSO4_A3 + zox_d*rSO4_B3 + zox_f

    !    print*,'rSO4_A3, rSO4_B3, rSO4_A2, rSO4_B2, rSO4_A1, rSO4_B1', rSO4_A3, rSO4_B3, rSO4_A2, rSO4_B2, rSO4_A1, rSO4_B1
    !    print*,' '

    END SUBROUTINE zSO4_calcbc

    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    FUNCTION calcFSO4(z)

        real calcFSO4, z, tmpreac1, tmpreac2

        ! Calculate SO4 consumption below zso4, by organic matter and indirectly via methane oxidation

        !    print*,' '
        !    print*,'..... START calcFSO4'

        tmpreac1    = MC*gammaCH4
        tmpreac2    = MC*gammaCH4

        calcFSO4 = calcReac(z, zinf, tmpreac1, tmpreac2)
    ! TODO confirm (1-por)*  has been added (to k1 & k2 ?)
    !    print*,'=============== IN calcFSO4 =====', calcFSO4

    END FUNCTION calcFSO4

    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************

    !                       Ammonium

    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE benthic_zNH4()

        !    real, intent(in)::zox, zno3

        ! local variables
        integer ltype1, ltype2, ltype3
        real ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1
        real ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2
        real ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3
        real e3_zinf, dedz3_zinf, f3_zinf, dfdz3_zinf, g3_zinf, dgdz3_zinf
        real e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3
        real e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3
        real zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f
        real e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox
        real e2_zox0, dedz2_zox0, f2_zox0, dfdz2_zox0, g2_zox0, dgdz2_zox0
        real e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox
        real zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
        real e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00
        real e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0
        real rNH4_A3, rNH4_B3
        real rNH4_A2, rNH4_B2
        real rNH4_A1, rNH4_B1

        real flxswiNH4
        real FNH4

        !    print*, ''
        !    print*, '------------------------------------------------------------------'
        !    print*, '---------------------- START zNH4 ------------------------------- '

        ! Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
        ! layer 1: 0 < z < zox, NH4 prod (remaining after oxidation)
        !      ls =      prepfg_l12( bsd, swi, r, reac1,     reac2,     ktemp, zU, zL, D1,        D2)
        call prepfg_l12((1-gamma)*NC1/(1.0+KNH4),(1-gamma)*NC2/(1.0+KNH4),0.0, 0.0, zox, DNH41, DNH42, &
        ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, ltype1)

        ! layer 2: zox < z < zno3, passive diffn TODO NH4 from denitrification?
        call prepfg_l12(0.0, 0.0, 0.0, zox, zno3, DNH41, DNH42, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, ltype2)

        ! layer 3: zno3 < z < zinf, NH4 production
        call prepfg_l12(NC1/(1.0+KNH4), NC2/(1.0+KNH4), 0.0, zno3, zinf, DNH41, DNH42, &
        ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, ltype3)

        ! Work up from the bottom, matching solutions at boundaries
        ! Basis functions at bottom of layer 3 zinf
        call calcfg_l12(zinf, NC1/(1.0+KNH4), NC2/(1.0+KNH4), 0.0, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, &
        DNH41, DNH42, ltype3, e3_zinf, dedz3_zinf, f3_zinf, dfdz3_zinf, g3_zinf, dgdz3_zinf)

        ! Match at zno3, layer 2 - layer 3 (continuity and flux)
        ! basis functions at bottom of layer 2
        call calcfg_l12(zno3, 0.0, 0.0, 0.0, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DNH41, DNH42, ltype2, &
        e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3)

        ! ... and top of layer 3
        call calcfg_l12(zno3, NC1/(1.0+KNH4), NC2/(1.0+KNH4), 0.0, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, &
        DNH41, DNH42, ltype3, e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3)

        ! match solutions at zno3 - continuous concentration and flux
        call matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, &
        e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, &
        0.0, 0.0, zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f)

        ! Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from NH4 sink)
        ! flux of NH4 to oxic interface  TODO NH4 prod by denitrification?
        FNH4 = calcReac(zno3, zinf, NC1/(1.0+KNH4), NC2/(1.0+KNH4))   ! MULTIPLY BY 1/POR ????

        ! basis functions at bottom of layer 1
        call calcfg_l12(zox, (1-gamma)*NC1/(1.0+KNH4),(1-gamma)*NC2/(1.0+KNH4) , 0.0, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, &
        DNH41, DNH42, ltype1, e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox)
        ! basis functions at top of layer 2
        call calcfg_l12(zox, 0.0, 0.0, 0.0, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DNH41, DNH42, ltype2, &
        e2_zox0, dedz2_zox0, f2_zox0, dfdz2_zox0, g2_zox0, dgdz2_zox0)

        ! transform to use coeffs from l3
        call xformsoln(e2_zox0, f2_zox0, g2_zox0, dedz2_zox0, dfdz2_zox0, dgdz2_zox0, &
        zno3_a , zno3_b , zno3_c , zno3_d , zno3_e ,zno3_f, &
        e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox)

        ! match solutions at zox - continuous concentration, flux discontinuity from NH4 ox
        IF(zox .le. zbio)THEN
            call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
            e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
            0.0, r_zxf*gamma*FNH4/DNH41, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
        ELSE
            call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
            e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
            0.0, r_zxf*gamma*FNH4/DNH42, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)

        END IF

        ! Solution at swi, top of layer 1
        call calcfg_l12(0.0, (1-gamma)*NC1/(1.0+KNH4), (1-gamma)*NC2/(1.0+KNH4), 0.0, &
        ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DNH41, DNH42, ltype1, &
        e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00)

        ! transform to use coeffs from l3
        call xformsoln(e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00, &
        zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f, &
        e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0)

        ! Solve for ANH4, BNH4 given boundary conditions (expressed in terms of transformed basis fns, layer 3 A, B)
        ! ANH4*dedz3_zinf   +  BNH4*dfdz3_zinf  + dgdz3_zinf = 0;
        ! ANH4*e1_0     +   BNH4*f1_0     + g1_0  = swi.NH40;

        ! | dedz3_zinf dfdz3_zinf |  |ANH4|   = | -dgdz3_zinf       |
        ! | e1_0     f1_0         |  |BNH4|     | swi.NH40 - g1_0 |

        call solve2eqn(dedz3_zinf, dfdz3_zinf, e1_0, f1_0, -dgdz3_zinf, NH40 - g1_0, rNH4_A3, rNH4_B3)

        ! flux at swi - DO include por so this is per cm^2 water column area
        flxswiNH4 = por*(DNH41*(rNH4_A3*dedz1_0+rNH4_B3*dfdz1_0 + dgdz1_0) - w*NH40)   ! NB: use A3, B3 as these are _xformed_ layer 1 basis functions

        ! save coeffs for layers 2 and 1
        rNH4_A2 = zno3_a*rNH4_A3 + zno3_b*rNH4_B3 + zno3_e
        rNH4_B2 = zno3_c*rNH4_A3 + zno3_d*rNH4_B3 + zno3_f

        rNH4_A1 = zox_a*rNH4_A3 + zox_b*rNH4_B3 + zox_e
        rNH4_B1 = zox_c*rNH4_A3 + zox_d*rNH4_B3 + zox_f


        !    print*,'INPUT1: por, DNH41, rNH4_A3, dedz1_0, rNH4_B3, dfdz1_0 , dgdz1_0', &
        !            & por, DNH41, rNH4_A3, dedz1_0, rNH4_B3, dfdz1_0 , dgdz1_0
        !    print*,'INPUT2: zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f', &
        !            & zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f
        !    print*,' '
        !    print*,' ---------------- RESULTS: benthic_zNH4 ---------------- '
        print*,'flxswiNH4', char(9), flxswiNH4
        print*,' '
    !    print*,'rNH4_A3, rNH4_B3, rNH4_A2, rNH4_B2, rNH4_A1, rNH4_B1', &
    !          &  rNH4_A3, rNH4_B3, rNH4_A2, rNH4_B2, rNH4_A1, rNH4_B1


    END SUBROUTINE benthic_zNH4

    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************

    !                           Hydrogen Sulfide

    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE benthic_zH2S()

        !    real, intent(in)::zox, zso4

        ! local variables
        real reac1_h2s, reac2_h2s                 ! reactive terms: OM degradation
        integer ltype1, ltype2, ltype3, ltype4
        real ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1
        real ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2
        real ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3
        real ls_a4, ls_b4, ls_c4, ls_d4, ls_e4, ls_f4
        real e4_zinf, dedz4_zinf, f4_zinf, dfdz4_zinf, g4_zinf, dgdz4_zinf
        real e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4
        real e4_zso4, dedz4_zso4, f4_zso4, dfdz4_zso4, g4_zso4, dgdz4_zso4
        real zso4_a, zso4_b, zso4_c, zso4_d, zso4_e, zso4_f
        real e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3
        real e3_zno30, dedz3_zno30, f3_zno30, dfdz3_zno30, g3_zno30, dgdz3_zno30
        real e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3
        real zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f
        real e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox
        real e2_zox0, dedz2_zox0, f2_zox0, dfdz2_zox0, g2_zox0, dgdz2_zox0
        real e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox
        real zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
        real e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00
        real e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0

        real rH2S_A4, rH2S_B4
        real rH2S_A3, rH2S_B3
        real rH2S_A2, rH2S_B2
        real rH2S_A1, rH2S_B1

        real zso4FH2S, zoxFH2S
        real flxswiH2S

        reac1_h2s=SO4C
        reac2_h2s=SO4C

        !    print*, ''
        !    print*, '------------------------------------------------------------------'
        !    print*, '---------------------- START zH2S ------------------------------- '
        !    print*, ''

        ! Calculate H2S

        ! Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
        ! layer 1: 0 < z < zox, passive diffn
        !      ls =      prepfg_l12( bsd, swi, r, reac1,     reac2,     ktemp, zU, zL, D1,        D2)
        call prepfg_l12(0.0, 0.0, 0.0, 0.0, zox, DH2S1, DH2S2, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, ltype1)

        ! layer 2: zox < z < zno3, passive diffn
        call prepfg_l12(0.0, 0.0, 0.0, zox, zno3, DH2S1, DH2S2, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, ltype2)

        ! layer 3: zno3 < z < zso4, H2S consumption by OM oxidation
        call prepfg_l12(reac1_h2s, reac2_h2s, 0.0, zno3, zso4, DH2S1, DH2S2, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, ltype3)

        ! layer 4: zso4 < z < zinf, passive diffn
        call prepfg_l12(0.0, 0.0, 0.0, zso4, zinf, DH2S1, DH2S2, ls_a4, ls_b4, ls_c4, ls_d4, ls_e4, ls_f4, ltype4)

        ! Work up from the bottom, matching solutions at boundaries
        ! Basis functions at bottom of layer 4 zinf
        call calcfg_l12(zinf, 0.0, 0.0, 0.0, ls_a4, ls_b4, ls_c4, ls_d4, ls_e4, ls_f4, DH2S1, DH2S2, ltype4, &
        e4_zinf, dedz4_zinf, f4_zinf, dfdz4_zinf, g4_zinf, dgdz4_zinf)

        ! Match at zso4, layer 3 - layer 4 (continuity and flux with AOM production)
        ! basis functions at bottom of layer 3
        call calcfg_l12(zso4, reac1_h2s, reac2_h2s, 0.0, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, DH2S1, DH2S2, ltype3, &
        e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4)

        ! ... and top of layer 4
        call calcfg_l12(zso4, 0.0,  0.0, 0.0, ls_a4, ls_b4, ls_c4, ls_d4, ls_e4, ls_f4, DH2S1, DH2S2, ltype4, &
        e4_zso4, dedz4_zso4, f4_zso4, dfdz4_zso4, g4_zso4, dgdz4_zso4)

        ! flux of H2S produced by AOM interface (Source of H2S)
        zso4FH2S = calcReac(zso4, zinf, MC, MC) ! MULTIPLY BY 1/POR ????
        !    print*,'flux of H2S produced by AOM interface zso4FH2S = ', zso4FH2S

        ! match solutions at zso4 - continuous concentration and flux
        call matchsoln(e3_zso4, f3_zso4, g3_zso4, dedz3_zso4, dfdz3_zso4, dgdz3_zso4, &
        e4_zso4, f4_zso4, g4_zso4, dedz4_zso4, dfdz4_zso4, dgdz4_zso4, &
        0.0, -gammaCH4*zso4FH2S/DH2S2, zso4_a, zso4_b, zso4_c, zso4_d, zso4_e, zso4_f)

        ! Match at zno3, layer 2 - layer 3 (continuity and flux)
        ! basis functions at bottom of layer 2
        call calcfg_l12(zno3, 0.0, 0.0, 0.0, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DH2S1, DH2S2, ltype2, &
        e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3)

        ! ... and top of layer 3
        call calcfg_l12(zno3, reac1_h2s, reac2_h2s, 0.0, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, DH2S1, DH2S2, ltype3, &
        e3_zno30, dedz3_zno30, f3_zno30, dfdz3_zno30, g3_zno30, dgdz3_zno30)

        ! ... transformed to use coeffs from l4
        call xformsoln(e3_zno30, f3_zno30, g3_zno30, dedz3_zno30, dfdz3_zno30, dgdz3_zno30, &
        zso4_a , zso4_b , zso4_c , zso4_d , zso4_e ,zso4_f, &
        e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3)
        ! match solutions at zno3 - continuous concentration and flux
        call matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, &
        e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, &
        0.0, 0.0, zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f)

        ! Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from H2S source)
        ! flux of H2S to oxic interface (from all sources of H2S below)
        ! NB: include methane region as AOM will produce sulphide as well..
        zoxFH2S = calcReac(zno3, zso4, SO4C, SO4C) &
        +  calcReac(zso4, zinf, MC, MC) ! MULTIPLY BY 1/POR ????
        !    print*,' '
        !    print*,'flux of H2S to oxic interface zoxFH2S = ', zoxFH2S
        !    print*,' '

        ! basis functions at bottom of layer 1
        call calcfg_l12(zox, 0.0, 0.0, 0.0, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DH2S1, DH2S2, ltype1, &
        e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox)

        ! basis functions at top of layer 2
        call calcfg_l12(zox, 0.0, 0.0, 0.0, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DH2S1, DH2S2, ltype2, &
        e2_zox0, dedz2_zox0, f2_zox0, dfdz2_zox0, g2_zox0, dgdz2_zox0)

        !   transform to use coeffs from l4
        call xformsoln(e2_zox0, f2_zox0, g2_zox0, dedz2_zox0, dfdz2_zox0, dgdz2_zox0, &
        zno3_a , zno3_b , zno3_c , zno3_d , zno3_e ,zno3_f, &
        e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox)

        ! match solutions at zox - continuous concentration, flux discontinuity from H2S ox

        IF(zox .le. zbio) THEN
            call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
            e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
            0.0, r_zxf*gammaH2S*zoxFH2S/DH2S1, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
        ELSE
            call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
            e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
            0.0, r_zxf*gammaH2S*zoxFH2S/DH2S2, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
        END IF

        ! Solution at swi, top of layer 1
        call calcfg_l12(0.0, 0.0, 0.0, 0.0, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DH2S1, DH2S2, ltype1, &
        e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00)

        ! transform to use coeffs from l4
        call xformsoln(e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00, &
        zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f, &
        e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0)

        ! Solve for AH2S, BH2S given boundary conditions (expressed in terms of transformed basis fns, layer 4 A, B)
        !  AH2S*dedz4_zinf   +  BH2S*dfz4_zinf  + dgz4_zinf = 0;          % zero flux at zinf
        !  AH2S*e1_0     +   BH2S*f1_0     + g1_0  = swi.H2S0;

        !  | dedz4_zinf dfdz4_zinf |  |AH2S|   = | -dgz4_zinf       |
        !  | e1_0     f1_0         |  |BH2S|     | swi.H2S0 - g1_0 |

        call solve2eqn(dedz4_zinf, dfdz4_zinf, e1_0, f1_0, -dgdz4_zinf, H2S0 - g1_0, rH2S_A4, rH2S_B4)

        ! flux at swi - DO include por so this is per cm^2 water column area
        flxswiH2S = por*(DH2S1*(rH2S_A4*dedz1_0+rH2S_B4*dfdz1_0 + dgdz1_0) - w*H2S0)   ! NB: use A4, B4 as these are _xformed_ layer 1 basis functions

        ! save coeffs for layers 3, 2 and 1
        rH2S_A3 = zso4_a*rH2S_A4 + zso4_b*rH2S_B4 + zso4_e
        rH2S_B3 = zso4_c*rH2S_A4 + zso4_d*rH2S_B4 + zso4_f

        rH2S_A2 = zno3_a*rH2S_A4 + zno3_b*rH2S_B4 + zno3_e
        rH2S_B2 = zno3_c*rH2S_A4 + zno3_d*rH2S_B4 + zno3_f

        rH2S_A1 = zox_a*rH2S_A4 + zox_b*rH2S_B4 + zox_e
        rH2S_B1 = zox_c*rH2S_A4 + zox_d*rH2S_B4 + zox_f

        !    print*,' ---------------- RESULTS: benthic_zH2S ---------------- '
        print*,'flxswiH2S ', char(9), flxswiH2S
    !    print*,'rH2S_A4, rH2S_B4, rH2S_A3, rH2S_B3, rH2S_A2, rH2S_B2, rH2S_A1, rH2S_B1', &
    !          &  rH2S_A4, rH2S_B4, rH2S_A3, rH2S_B3, rH2S_A2, rH2S_B2, rH2S_A1, rH2S_B1


    !   print*,'INPUT1: dedz4_zinf, dfdz4_zinf, e1_0, f1_0, -dgdz4_zinf, H2S0 - g1_0', &
    !           & dedz4_zinf, dfdz4_zinf, e1_0, f1_0, -dgdz4_zinf, H2S0 - g1_0
    !    print*,'INPUT2: zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f', &
    !            & zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f
    !    print*,'RESULTS:  rH2S_A4, rH2S_B4', &
    !            & rH2S_A4, rH2S_B4

    END SUBROUTINE benthic_zH2S

    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************

    !                           Phosphate and Fe-bound P

    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE benthic_zPO4_M()
        !    real, intent(in)::zox, zno3
        !    real, intent(inout)::zso4

        ! local variables
        ! Integration constants
        real rPO4_M_A2, rPO4_M_B2, rPO4_M_C2, rPO4_M_D2
        real rPO4_M_A1, rPO4_M_B1, rPO4_M_C1, rPO4_M_D1

        integer dum_ltype1, dum_ltype2
        real, dimension (4, 4) :: dum_mat_C1, dum_mat_C2
        real, dimension (1:4) ::  dum_vec_D1, dum_vec_D2
        integer loc_i, loc_j

        real reac1_po4_ox, reac2_po4_ox, reac1_po4_anox, reac2_po4_anox                 ! reactive terms: OM degradation

        ! all the base functions for boundary matching
        ! i.e. ODE solutions (E, F, P, Q) and the particulat integral (G) and their derivatives
        ! at zinf
        real e2_zinf_P, dedz2_zinf_P, f2_zinf_P, dfdz2_zinf_P, g2_zinf_P, dgdz2_zinf_P
        real p2_zinf_P, dpdz2_zinf_P, q2_zinf_P, dqdz2_zinf_P
        real e2_zinf_M, dedz2_zinf_M, f2_zinf_M, dfdz2_zinf_M, g2_zinf_M, dgdz2_zinf_M
        real p2_zinf_M, dpdz2_zinf_M, q2_zinf_M, dqdz2_zinf_M
        ! at zox (above)
        real e1_zox_P, dedz1_zox_P, f1_zox_P, dfdz1_zox_P, g1_zox_P, dgdz1_zox_P
        real p1_zox_P, dpdz1_zox_P, q1_zox_P, dqdz1_zox_P
        real e1_zox_M, dedz1_zox_M, f1_zox_M, dfdz1_zox_M, g1_zox_M, dgdz1_zox_M
        real p1_zox_M, dpdz1_zox_M, q1_zox_M, dqdz1_zox_M
        ! at zox (below)
        real e2_zox_P, dedz2_zox_P, f2_zox_P, dfdz2_zox_P, g2_zox_P, dgdz2_zox_P
        real p2_zox_P, dpdz2_zox_P, q2_zox_P, dqdz2_zox_P
        real e2_zox_M, dedz2_zox_M, f2_zox_M, dfdz2_zox_M, g2_zox_M, dgdz2_zox_M
        real p2_zox_M, dpdz2_zox_M, q2_zox_M, dqdz2_zox_M
        ! at SWI (z0)
        real e1_z0_P, dedz1_z0_P, f1_z0_P, dfdz1_z0_P, g1_z0_P, dgdz1_z0_P
        real p1_z0_P, dpdz1_z0_P, q1_z0_P, dqdz1_z0_P
        real e1_z0_M, dedz1_z0_M, f1_z0_M, dfdz1_z0_M, g1_z0_M, dgdz1_z0_M
        real p1_z0_M, dpdz1_z0_M, q1_z0_M, dqdz1_z0_M
        ! the final g's (other are saved in matrices, e.g. loc_EFPQ_P_t, ...)
        real g_P, dgdz_P
        real g_M, dgdz_M


        real loc_Vb, loc_Fb                                   ! discontinuity constants
        real, dimension (1:4,1:4) :: loc_mat_X_4x4, loc_mat_Y_4x4
        real, dimension (1:4) ::  loc_vec_Z_4
        real, dimension (1:3,1:3) :: loc_mat_X_3x3, loc_mat_Y_3x3
        real, dimension (1:3) ::  loc_vec_Z_3

        ! Matrix/Vector calculated by matchsoln_PO4_M
        real, dimension (4, 4) :: loc_mat_C_4x4
        real, dimension (1:4) ::  loc_vec_D_4
        real, dimension (3, 3) :: loc_mat_C_3x3
        real, dimension (1:3) ::  loc_vec_D_3
        integer loc_dim

        ! save the ODE solutions in vectors to make calculation easier (DH?: however, is this faster?)
        real, dimension (1:4) :: loc_EFPQ_P, loc_dEFPQdz_P, loc_EFPQ_M, loc_dEFPQdz_M
        ! the transformed ODE solutions coming from xformsoln_PO4_M
        real, dimension (1:4) :: loc_EFPQ_P_t, loc_dEFPQdz_P_t, loc_EFPQ_M_t, loc_dEFPQdz_M_t

        ! the final SWI fluxes
        real loc_flxswipo4, loc_flxswiM
        ! calculated integration constants for Layer 1 & 2 - in case we want to calculate profiles later - calculated first as vector to save code
        real, dimension (1:4) :: loc_Layer1_IC, loc_Layer2_IC

        reac1_po4_ox = 1/(1+KPO4_ox)*PC1
        reac2_po4_ox = 1/(1+KPO4_ox)*PC2
        reac1_po4_anox = 1/(1+KPO4_anox)*PC1
        reac2_po4_anox = 1/(1+KPO4_anox)*PC2

        loc_Vb = 0.0
        loc_Fb = 0.0

        ! Initialize loc_mat_C_4x4 & loc_vec_D_4 with zeros as, so I don't need to add them later manually ...
        loc_mat_C_4x4 = 0.0
        loc_vec_D_4 = 0.0

        !    print*, ''
        !    print*, '------------------------------------------------------------------'
        !    print*, '---------------------- START zPO4_M ------------------------------- '

        !   Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary

        !   layer 1: 0 < z < zox, OM degradation (-) Sorption to sediment Fe-oxides (ktemp)

        call prepfg_l12_PO4_M(reac1_po4_ox, reac2_po4_ox, ksPO4/(1+KPO4_ox), PO4s*ksPO4/(1+KPO4_ox), 0.0, zox, &
                                DPO41/(1+KPO4_ox), DPO42/(1+KPO4_ox), 0.0, 0.0, 0.0, Dbio, 0.0, (1/SD)*ksPO4, &
                                dum_mat_C1, dum_vec_D1, dum_ltype1)

        !   layer 2: zox < z < zinf,
        !   OM degradation (-) authigenic P formation (ktemp) (+) P desorption due to Fe-bound P release upon Fe oxide reduction
        !       rPO4_M.ls2 = r.zTOC.prepfg_l12_PO4_M(bsd, swi, r, obj.reac1_anox, obj.reac2_anox, obj.kaPO4/(1+obj.KPO4_anox), obj.PO4a*obj.kaPO4/(1+obj.KPO4_anox), r.zox, bsd.zinf, obj.DPO41/(1+obj.KPO4_anox), obj.DPO42/(1+obj.KPO4_anox), bsd.SD*obj.kmPO4/(1+obj.KPO4_anox), ...
        !                                           obj.kmPO4, obj.kmPO4.*obj.Minf, bsd.Dbio, 0, 0);
        call prepfg_l12_PO4_M(reac1_po4_anox, reac2_po4_anox, kaPO4/(1+KPO4_anox), PO4a*kaPO4/(1+KPO4_anox), zox, zinf, &
                                DPO41/(1+KPO4_anox), DPO42/(1+KPO4_anox), SD*kmPO4/(1+KPO4_anox), kmPO4, kmPO4*Minf, Dbio, &
                                0.0, 0.0, dum_mat_C2, dum_vec_D2, dum_ltype2)

            ! Work up from the bottom, matching solutions at boundaries
            ! Basis functions at bottom of layer 2 zinf
        call calcfg_l12_PO4_M(zinf, reac1_po4_anox, reac2_po4_anox, kaPO4/(1+KPO4_anox), PO4a*kaPO4/(1+KPO4_anox), &
                                DPO41/(1+KPO4_anox), DPO42/(1+KPO4_anox), SD*kmPO4/(1+KPO4_anox), dum_mat_C2, dum_vec_D2, &
                                dum_ltype2, kmPO4, kmPO4*Minf, Dbio, 0.0, 0.0, e2_zinf_P, dedz2_zinf_P, f2_zinf_P, &
                                dfdz2_zinf_P, g2_zinf_P, dgdz2_zinf_P, p2_zinf_P, dpdz2_zinf_P, q2_zinf_P, dqdz2_zinf_P, &
                                e2_zinf_M, dedz2_zinf_M, f2_zinf_M, dfdz2_zinf_M, g2_zinf_M, dgdz2_zinf_M, &
                                p2_zinf_M, dpdz2_zinf_M, q2_zinf_M, dqdz2_zinf_M)

        ! Match at zox, layer 1 - layer 2 (continuity and flux)
        ! basis functions at bottom of layer 1
        call calcfg_l12_PO4_M(zox, reac1_po4_ox, reac2_po4_ox, ksPO4/(1+KPO4_ox), PO4s*ksPO4/(1+KPO4_ox), &
                                DPO41/(1+KPO4_ox), DPO42/(1+KPO4_ox), 0.0, dum_mat_C1, dum_vec_D1, &
                                dum_ltype1, 0.0, 0.0, Dbio, 0.0, (1/SD)*ksPO4, e1_zox_P, dedz1_zox_P, f1_zox_P, &
                                dfdz1_zox_P, g1_zox_P, dgdz1_zox_P, p1_zox_P, dpdz1_zox_P, q1_zox_P, dqdz1_zox_P, &
                                e1_zox_M, dedz1_zox_M, f1_zox_M, dfdz1_zox_M, g1_zox_M, dgdz1_zox_M, &
                                p1_zox_M, dpdz1_zox_M, q1_zox_M, dqdz1_zox_M)

        !  and top of layer 2
         call calcfg_l12_PO4_M(zox, reac1_po4_anox, reac2_po4_anox, kaPO4/(1+KPO4_anox), PO4a*kaPO4/(1+KPO4_anox), &
                                DPO41/(1+KPO4_anox), DPO42/(1+KPO4_anox), SD*kmPO4/(1+KPO4_anox), dum_mat_C2, dum_vec_D2, &
                                dum_ltype2, kmPO4, kmPO4*Minf, Dbio, 0.0, 0.0, e2_zox_P, dedz2_zox_P, f2_zox_P, &
                                dfdz2_zox_P, g2_zox_P, dgdz2_zox_P, p2_zox_P, dpdz2_zox_P, q2_zox_P, dqdz2_zox_P, &
                                e2_zox_M, dedz2_zox_M, f2_zox_M, dfdz2_zox_M, g2_zox_M, dgdz2_zox_M, &
                                p2_zox_M, dpdz2_zox_M, q2_zox_M, dqdz2_zox_M)

        ! match solutions at zox - continuous concentration and flux
        ! organize the data in matrices and let the intrinsic fortran function do the calculation
        ! DH: Maybe this could be done more efficiently !?
        !  |x1        |   | A_l |      | y1        | | A_r|    |z1|    always PO4 continuity
        !  |    .     |   | B_l |      |    .      | | B_r|    |z2|    always PO4 flux
        !  |      .   |   | C_l |   =  |      .    | | C_r|  + |z3|    always M continuity
        !  |       x16|   | D_l |      |        y16| | D_r|    |z4|    SD M flux only in bioturbated case, otherwise not an independent constraint

        if(zox < zbio)then  ! 1. CASE: 4 int const. in each layer
            ! weird FORTRAN matrices makes the transpose necessary
            loc_mat_X_4x4 = transpose(reshape((/ e1_zox_P, f1_zox_P, p1_zox_P, q1_zox_P, &
            dedz1_zox_P, dfdz1_zox_P, dpdz1_zox_P, dqdz1_zox_P, &
            e1_zox_M, f1_zox_M, p1_zox_M, q1_zox_M, &
            dedz1_zox_M, dfdz1_zox_M, dpdz1_zox_M, dqdz1_zox_M /), shape(loc_mat_X_4x4)))

            loc_mat_Y_4x4 = transpose(reshape((/ e2_zox_P, f2_zox_P, p2_zox_P, q2_zox_P, &
            dedz2_zox_P, dfdz2_zox_P, dpdz2_zox_P, dqdz2_zox_P, &
            e2_zox_M, f2_zox_M, p2_zox_M, q2_zox_M, &
            dedz2_zox_M, dfdz2_zox_M, dpdz2_zox_M, dqdz2_zox_M /), shape(loc_mat_Y_4x4)))

            loc_vec_Z_4 = (/ g2_zox_P-g1_zox_P + loc_Vb, &
            dgdz2_zox_P - dgdz1_zox_P + loc_Fb - w*loc_Vb, &
            g2_zox_M-g1_zox_M + loc_Vb, &
            dgdz2_zox_M - dgdz1_zox_M + loc_Fb - w*loc_Vb /)

            loc_dim = 4
            call matchsoln_PO4_M(loc_mat_X_4x4, loc_mat_Y_4x4, loc_vec_Z_4, loc_dim, loc_mat_C_4x4, loc_vec_D_4)


        else    ! 2. CASE: 3 int const. in each layer
                ! DH: zox non-bioturbated
                ! DH: this should generate 3x3 matrices as no M flux boundary condition (and then 4x4 C with zeros)

            loc_mat_X_3x3 = transpose(reshape((/ e1_zox_P, f1_zox_P, p1_zox_P, &
            dedz1_zox_P, dfdz1_zox_P, dpdz1_zox_P, &
            e1_zox_M, f1_zox_M, p1_zox_M /), shape(loc_mat_X_3x3)))

            loc_mat_Y_3x3 = transpose(reshape((/ e2_zox_P, f2_zox_P, p2_zox_P, &
            dedz2_zox_P, dfdz2_zox_P, dpdz2_zox_P, &
            e2_zox_M, f2_zox_M, p2_zox_M /), shape(loc_mat_Y_3x3)))

            loc_vec_Z_3 = (/ g2_zox_P-g1_zox_P + loc_Vb, &
            dgdz2_zox_P - dgdz1_zox_P + loc_Fb - w*loc_Vb, &
            g2_zox_M-g1_zox_M + loc_Vb /)
            loc_dim = 3

            call matchsoln_PO4_M(loc_mat_X_3x3, loc_mat_Y_3x3, loc_vec_Z_3, loc_dim, loc_mat_C_3x3, loc_vec_D_3)

! integrate the 3x3 matrix in the 4x4 matrix
            do loc_i = 1,3
                do loc_j = 1,3
                    loc_mat_C_4x4(loc_i, loc_j) = loc_mat_C_3x3(loc_i, loc_j)
                end do
                loc_vec_D_4(loc_i)  = loc_vec_D_3(loc_i)
            end do

        end if !(zox < zbio)

!201     format (6f12.6)
!        print*,'final loc_mat_C_4x4 '
!        do loc_i=1,4
!            write (*,201) (loc_mat_C_4x4(loc_i,loc_j),loc_j=1,4)
!        end do
!        print*,'loc_vec_D_4 ', char(9), loc_vec_D_4

        ! Solution at SWI, top of layer 1
        call calcfg_l12_PO4_M(z0, reac1_po4_ox, reac2_po4_ox, ksPO4/(1+KPO4_ox), PO4s*ksPO4/(1+KPO4_ox), &
                                DPO41/(1+KPO4_ox), DPO42/(1+KPO4_ox), 0.0, dum_mat_C1, dum_vec_D1, &
                                dum_ltype1, 0.0, 0.0, Dbio, 0.0, (1/SD)*ksPO4, e1_z0_P, dedz1_z0_P, f1_z0_P, &
                                dfdz1_z0_P, g1_z0_P, dgdz1_z0_P, p1_z0_P, dpdz1_z0_P, q1_z0_P, dqdz1_z0_P, &
                                e1_z0_M, dedz1_z0_M, f1_z0_M, dfdz1_z0_M, g1_z0_M, dgdz1_z0_M, &
                                p1_z0_M, dpdz1_z0_M, q1_z0_M, dqdz1_z0_M)

        ! transform to use coeffs from l2
        ! Now find 'transformed' basis functions such that in layer 1 (here for O2, PO4 is a bit more complex)
        ! O2 = A_2*et + B_2*ft + gt  (ie layer 1 soln written in terms of layer 2 coeffs A_2, B_2)
        loc_EFPQ_P = (/ e1_z0_P, f1_z0_P, p1_z0_P, q1_z0_P /)
        loc_dEFPQdz_P = (/ dedz1_z0_P, dfdz1_z0_P, dpdz1_z0_P, dqdz1_z0_P /)
        loc_EFPQ_M = (/ e1_z0_M, f1_z0_M, p1_z0_M, q1_z0_M /)
        loc_dEFPQdz_M = (/dedz1_z0_M, dfdz1_z0_M, dpdz1_z0_M, dqdz1_z0_M /)

        call xformsoln_PO4_M(loc_EFPQ_P, loc_EFPQ_M, loc_dEFPQdz_P, loc_dEFPQdz_M, &
                g1_z0_P, g1_z0_M, dgdz1_z0_P, dgdz1_z0_M, loc_mat_C_4x4, loc_vec_D_4, &
                loc_EFPQ_P_t, g_P, loc_dEFPQdz_P_t, dgdz_P, loc_EFPQ_M_t, &
                g_M, loc_dEFPQdz_M_t, dgdz_M)

        ! SD assume D2 == 0 (as q, dqdz2_zinf = 0 ) and solve for 3 unknowns
        loc_mat_X_3x3 =  transpose(reshape((/ dedz2_zinf_P, dfdz2_zinf_P, dpdz2_zinf_P, &
            loc_EFPQ_P_t(1), loc_EFPQ_P_t(2), loc_EFPQ_P_t(3), &
            w*loc_EFPQ_M_t(1), w*loc_EFPQ_M_t(2), w*loc_EFPQ_M_t(3) /), shape(loc_mat_X_3x3)))
        loc_vec_Z_3= (/ -dgdz2_zinf_P, &
                            PO40 - g_P, &
                           Mflux0 - w*g_M  /)

        ! calculate the integration conctants for Layer 2
        ! just need it once, so actually no need for subroutine, but maybe for later
        loc_dim = 3
        call solve2eqn_PO4_M(loc_mat_X_3x3, loc_vec_Z_3, rPO4_M_A2, rPO4_M_B2, rPO4_M_C2, loc_dim)
        rPO4_M_D2 = 0.0
        ! save IC in a vector for a later calculation
        loc_Layer2_IC = (/ rPO4_M_A2, rPO4_M_B2, rPO4_M_C2, rPO4_M_D2 /)

        ! CALCULATE FINAL SWI fluxes and save the coefficients for
        ! DH: use A2, B2, C2, D2 as these are _xformed_ layer 1 basis functions
        loc_flxswipo4 = por*(DPO41/(1+KPO4_ox)*(rPO4_M_A2*loc_dEFPQdz_P_t(1)+rPO4_M_B2*loc_dEFPQdz_P_t(2) &
                    + rPO4_M_C2*loc_dEFPQdz_P_t(3)+rPO4_M_D2*loc_dEFPQdz_P_t(4) + dgdz_P) - w*PO40)
        ! Does actually not exist, as it is a solid, just calculate for debugging
        loc_flxswiM = por*Dbio*(rPO4_M_A2*loc_dEFPQdz_M_t(1)+rPO4_M_B2*loc_dEFPQdz_M_t(2) + &
                        rPO4_M_C2*loc_dEFPQdz_M_t(3) + rPO4_M_D2*loc_dEFPQdz_M_t(4) + dgdz_M)

        ! save coeffs for layer 1 - in case I want to calculate a profile later
        loc_Layer1_IC = matmul(loc_mat_C_4x4, loc_Layer2_IC) + loc_vec_D_4

        print*,' '
        print*,'loc_flxswipo4 ', char(9), loc_flxswipo4
        print*,'loc_flxswiM ', char(9), loc_flxswiM
!        print*,'loc_Layer1_IC ', char(9), loc_Layer1_IC
!        print*,'loc_Layer2_IC ', char(9), loc_Layer2_IC

    END SUBROUTINE benthic_zPO4_M


    !------------------------------------------------------------------------------------
    !   *****************************************************************
    !   *****************************************************************

    !                   Dissolved Inorganic Carbon (DIC)

    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE benthic_zDIC()


        ! local variables
        real reac1_dic, reac2_dic                 ! reactive terms: OM degradation
        integer ltype1, ltype2, ltype3, ltype4
        real ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1
        real ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2
!        real e4_zinf, dedz4_zinf, f4_zinf, dfdz4_zinf, g4_zinf, dgdz4_zinf
!        real e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4
!        real e4_zso4, dedz4_zso4, f4_zso4, dfdz4_zso4, g4_zso4, dgdz4_zso4
!        real zso4_a, zso4_b, zso4_c, zso4_d, zso4_e, zso4_f
!        real e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3
!        real e3_zno30, dedz3_zno30, f3_zno30, dfdz3_zno30, g3_zno30, dgdz3_zno30
!        real e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3
!        real zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f
!        real e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox
!        real e2_zox0, dedz2_zox0, f2_zox0, dfdz2_zox0, g2_zox0, dgdz2_zox0
!        real e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox
!        real zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
!        real e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00
!        real e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0
!
!        real rH2S_A4, rH2S_B4
!        real rH2S_A3, rH2S_B3
!        real rH2S_A2, rH2S_B2
!        real rH2S_A1, rH2S_B1
!
!        real zso4FH2S, zoxFH2S
!        real flxswiH2S

        reac1_dic = DICC1                 ! DIC/C until zSO4 (mol/mol)
        reac2_dic = DICC2                 !DIC/C below zSO4 (mol/mol)

        !    print*, ''
        !    print*, '------------------------------------------------------------------'
        !    print*, '---------------------- START zDIC ------------------------------- '
        !    print*, ''

        ! Calculate DIC

        ! Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
        ! layer 1: 0 < z < zso4, DIC produced my OM degradation
        !    prepfg_l12(reac1,      reac2,  ktemp, zU,  zL,     D1, D2, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ltype)
        call prepfg_l12(reac1_dic, reac1_dic, 0.0, 0.0, zso4, DDIC1, DDIC2, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, ltype1)

        ! layer 2: zso4 < z < zinf, DIC production by OM degradation (Methanogenesis) -> different production rate
        call prepfg_l12(reac2_dic, reac2_dic, 0.0, zso4, zinf, DDIC1, DDIC2, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, ltype2)


        ! Work up from the bottom, matching solutions at boundaries
        ! Basis functions at bottom of layer 2 zinf
        !   calcfg_l12(  z,     reac1,      reac2, ktemp, ls_a, ls_b,  ls_c, ls_d,   ls_e, ls_f,  ls_D1, ls_D2, ltype, e, dedz, f, dfdz, g, dgdz)
        call calcfg_l12(zinf, reac2_dic, reac2_dic, 0.0, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DDIC1, DDIC2, ltype2, &
        e2_zinf, dedz2_zinf, f2_zinf, dfdz2_zinf, g2_zinf, dgdz2_zinf)
! DH 17.11.2016 TODO: Go on from here !!!!
        ! Match at zso4, layer 3 - layer 4 (continuity and flux with AOM production)
        ! basis functions at bottom of layer 3
        call calcfg_l12(zso4, reac1_h2s, reac2_h2s, 0.0, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, DH2S1, DH2S2, ltype3, &
        e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4)

        ! ... and top of layer 4
        call calcfg_l12(zso4, 0.0,  0.0, 0.0, ls_a4, ls_b4, ls_c4, ls_d4, ls_e4, ls_f4, DH2S1, DH2S2, ltype4, &
        e4_zso4, dedz4_zso4, f4_zso4, dfdz4_zso4, g4_zso4, dgdz4_zso4)

        ! flux of H2S produced by AOM interface (Source of H2S)
        zso4FH2S = calcReac(zso4, zinf, MC, MC) ! MULTIPLY BY 1/POR ????
        !    print*,'flux of H2S produced by AOM interface zso4FH2S = ', zso4FH2S

        ! match solutions at zso4 - continuous concentration and flux
        call matchsoln(e3_zso4, f3_zso4, g3_zso4, dedz3_zso4, dfdz3_zso4, dgdz3_zso4, &
        e4_zso4, f4_zso4, g4_zso4, dedz4_zso4, dfdz4_zso4, dgdz4_zso4, &
        0.0, -gammaCH4*zso4FH2S/DH2S2, zso4_a, zso4_b, zso4_c, zso4_d, zso4_e, zso4_f)

        ! Match at zno3, layer 2 - layer 3 (continuity and flux)
        ! basis functions at bottom of layer 2
        call calcfg_l12(zno3, 0.0, 0.0, 0.0, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DH2S1, DH2S2, ltype2, &
        e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3)

        ! ... and top of layer 3
        call calcfg_l12(zno3, reac1_h2s, reac2_h2s, 0.0, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, DH2S1, DH2S2, ltype3, &
        e3_zno30, dedz3_zno30, f3_zno30, dfdz3_zno30, g3_zno30, dgdz3_zno30)

        ! ... transformed to use coeffs from l4
        call xformsoln(e3_zno30, f3_zno30, g3_zno30, dedz3_zno30, dfdz3_zno30, dgdz3_zno30, &
        zso4_a , zso4_b , zso4_c , zso4_d , zso4_e ,zso4_f, &
        e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3)
        ! match solutions at zno3 - continuous concentration and flux
        call matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, &
        e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, &
        0.0, 0.0, zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f)

        ! Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from H2S source)
        ! flux of H2S to oxic interface (from all sources of H2S below)
        ! NB: include methane region as AOM will produce sulphide as well..
        zoxFH2S = calcReac(zno3, zso4, SO4C, SO4C) &
        +  calcReac(zso4, zinf, MC, MC) ! MULTIPLY BY 1/POR ????
        !    print*,' '
        !    print*,'flux of H2S to oxic interface zoxFH2S = ', zoxFH2S
        !    print*,' '

        ! basis functions at bottom of layer 1
        call calcfg_l12(zox, 0.0, 0.0, 0.0, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DH2S1, DH2S2, ltype1, &
        e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox)

        ! basis functions at top of layer 2
        call calcfg_l12(zox, 0.0, 0.0, 0.0, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DH2S1, DH2S2, ltype2, &
        e2_zox0, dedz2_zox0, f2_zox0, dfdz2_zox0, g2_zox0, dgdz2_zox0)

        !   transform to use coeffs from l4
        call xformsoln(e2_zox0, f2_zox0, g2_zox0, dedz2_zox0, dfdz2_zox0, dgdz2_zox0, &
        zno3_a , zno3_b , zno3_c , zno3_d , zno3_e ,zno3_f, &
        e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox)

        ! match solutions at zox - continuous concentration, flux discontinuity from H2S ox

        IF(zox .le. zbio) THEN
            call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
            e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
            0.0, r_zxf*gammaH2S*zoxFH2S/DH2S1, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
        ELSE
            call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
            e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
            0.0, r_zxf*gammaH2S*zoxFH2S/DH2S2, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
        END IF

        ! Solution at swi, top of layer 1
        call calcfg_l12(0.0, 0.0, 0.0, 0.0, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DH2S1, DH2S2, ltype1, &
        e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00)

        ! transform to use coeffs from l4
        call xformsoln(e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00, &
        zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f, &
        e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0)

        ! Solve for AH2S, BH2S given boundary conditions (expressed in terms of transformed basis fns, layer 4 A, B)
        !  AH2S*dedz4_zinf   +  BH2S*dfz4_zinf  + dgz4_zinf = 0;          % zero flux at zinf
        !  AH2S*e1_0     +   BH2S*f1_0     + g1_0  = swi.H2S0;

        !  | dedz4_zinf dfdz4_zinf |  |AH2S|   = | -dgz4_zinf       |
        !  | e1_0     f1_0         |  |BH2S|     | swi.H2S0 - g1_0 |

        call solve2eqn(dedz4_zinf, dfdz4_zinf, e1_0, f1_0, -dgdz4_zinf, H2S0 - g1_0, rH2S_A4, rH2S_B4)

        ! flux at swi - DO include por so this is per cm^2 water column area
        flxswiH2S = por*(DH2S1*(rH2S_A4*dedz1_0+rH2S_B4*dfdz1_0 + dgdz1_0) - w*H2S0)   ! NB: use A4, B4 as these are _xformed_ layer 1 basis functions

        ! save coeffs for layers 3, 2 and 1
        rH2S_A3 = zso4_a*rH2S_A4 + zso4_b*rH2S_B4 + zso4_e
        rH2S_B3 = zso4_c*rH2S_A4 + zso4_d*rH2S_B4 + zso4_f

        rH2S_A2 = zno3_a*rH2S_A4 + zno3_b*rH2S_B4 + zno3_e
        rH2S_B2 = zno3_c*rH2S_A4 + zno3_d*rH2S_B4 + zno3_f

        rH2S_A1 = zox_a*rH2S_A4 + zox_b*rH2S_B4 + zox_e
        rH2S_B1 = zox_c*rH2S_A4 + zox_d*rH2S_B4 + zox_f

        !    print*,' ---------------- RESULTS: benthic_zH2S ---------------- '
        print*,'flxswiH2S ', char(9), flxswiH2S
    !    print*,'rH2S_A4, rH2S_B4, rH2S_A3, rH2S_B3, rH2S_A2, rH2S_B2, rH2S_A1, rH2S_B1', &
    !          &  rH2S_A4, rH2S_B4, rH2S_A3, rH2S_B3, rH2S_A2, rH2S_B2, rH2S_A1, rH2S_B1


    !   print*,'INPUT1: dedz4_zinf, dfdz4_zinf, e1_0, f1_0, -dgdz4_zinf, H2S0 - g1_0', &
    !           & dedz4_zinf, dfdz4_zinf, e1_0, f1_0, -dgdz4_zinf, H2S0 - g1_0
    !    print*,'INPUT2: zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f', &
    !            & zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f
    !    print*,'RESULTS:  rH2S_A4, rH2S_B4', &
    !            & rH2S_A4, rH2S_B4

    END SUBROUTINE benthic_zDIC


    ! *******************************************************************
    !   *****************************************************************

    !                       UTILITY FUNCTIONS

    !   *****************************************************************
    ! *******************************************************************


    SUBROUTINE matchsoln(E_l, F_l, G_l, dEdx_l, dFdx_l, dGdx_l, &
    E_r, F_r, G_r, dEdx_r, dFdx_r, dGdx_r, &
    Vb, Db, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f)

        real, INTENT(in):: E_l, F_l, G_l, dEdx_l, dFdx_l, dGdx_l
        real, INTENT(in):: E_r, F_r, G_r, dEdx_r, dFdx_r, dGdx_r, Vb, Db
        real, INTENT(inout):: ls_a, ls_b, ls_c, ls_d, ls_e, ls_f
        real:: alden, blden

        ! Match two solutions at a boundary:
        ! 'left' solution   y_l(x) = A_l*E_l(x) + B_l*F_l(x) + G_l(x)
        ! 'right' solution  y_r(x) = A_r*E_r(x) + B_r*F_l(x) + G_r(x)
        !
        ! (Dis)continuity conditions at boundary:
        !                   y_r(xb)    = y_l(xb)     + Vb
        !                   dydx_r(xb) = dydx_l(xb)  + Db
        !
        ! Find a,b,c,d,e,f such that:
        !         | A_l |   =  | a  b | | A_r|  + |e|
        !         | B_l |      | c  d | | B_r|    |f|

        alden = dFdx_l*E_l - F_l*dEdx_l
        ls_a     = (dFdx_l*E_r - F_l*dEdx_r)/alden
        ls_b     = (dFdx_l*F_r - F_l*dFdx_r)/alden
        ls_e     = (F_l*(dGdx_l - dGdx_r + Db) + dFdx_l*(-G_l + G_r - Vb))/alden

        blden = dEdx_l*F_l - E_l*dFdx_l
        ls_c     = (dEdx_l*E_r - E_l*dEdx_r)/blden
        ls_d     = (dEdx_l*F_r - E_l*dFdx_r)/blden;
        ls_f     = (E_l*(dGdx_l - dGdx_r + Db) + dEdx_l*(-G_l+G_r - Vb))/blden

    END SUBROUTINE matchsoln


    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE matchsoln_PO4_M(dum_mat_X, dum_mat_Y, dum_vec_Z, dum_dim, loc_mat_C, loc_vec_D)

        integer dum_dim                      ! dimension of the matrices
        real, dimension (1:dum_dim,1:dum_dim), intent (in) :: dum_mat_X, dum_mat_Y
        real, dimension (1:dum_dim), intent (in) ::  dum_vec_Z
        real, dimension (1:dum_dim,1:dum_dim), intent (inout) :: loc_mat_C
        real, dimension (1:dum_dim), intent (inout) ::  loc_vec_D

        ! local variables
        real, dimension (1:dum_dim,1:dum_dim) :: loc_mat_X, loc_inv_mat_X
        !        integer :: loc_i, loc_j

        ! Match four solutions at a boundary:
        !  for PO4
        ! 'left' solution   y_l(z) = A_l*E_l(z) + B_l*F_l(z) + C_l*P_l(z) + D_l*Q_l(z) + G_l(z)
        ! 'right' solution  y_r(z) = A_r*E_r(z) + B_r*F_r(z) + C_r*P_r(z) + D_r*Q_r(z) + G_r(z)
        !
        !  and the same for M
        !
        ! (Dis)continuity conditions at boundary:
        !                   y_r(xb)    = y_l(xb)     + Vb
        !                   dydx_r(xb) = dydx_l(xb)  + Db

        !         | A_l |         | A_r|
        !         | B_l |         | B_r|
        !     X   | C_l |   =  Y  | C_r|  + Z
        !         | D_l |         | D_r|

        !
        ! Find C and D such that:
        !         | A_l |         | A_r|
        !         | B_l |         | B_r|
        !         | C_l |   =  C  | C_r|  +  D
        !         | D_l |         | D_r|

        ! save matrix locally, as the original matrix loc_mat_X(4,4) will be destroyed during the calculation
        loc_mat_X = dum_mat_X
        ! calculate loc_mat_X^{-1}
        call inverse(loc_mat_X,loc_inv_mat_X,dum_dim)
        loc_mat_C = matmul(loc_inv_mat_X, dum_mat_Y)
        loc_vec_D = matmul(loc_inv_mat_X, dum_vec_Z)


    END SUBROUTINE matchsoln_PO4_M


    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------


    SUBROUTINE xformsoln(E, F, G, dEdx, dFdx, dGdx, ls_a , ls_b , ls_c , ls_d , ls_e ,ls_f, Et, Ft, Gt, dEtdx, dFtdx, dGtdx)

        real, INTENT(in):: E, F, G, dEdx, dFdx, dGdx, ls_a , ls_b , ls_c , ls_d , ls_e ,ls_f
        real, INTENT(inout):: Et, Ft, Gt, dEtdx, dFtdx, dGtdx

        ! Find 'transformed' soln such that in layer l,
        !    y_l = A_r*et + B_r*ft + gt
        ! (ie l soln written in terms of r solution coefficents A_r, B_r)

        Et      = ls_a*E    + ls_c*F
        dEtdx   = ls_a*dEdx + ls_c*dFdx
        Ft      = ls_b*E    + ls_d*F
        dFtdx   = ls_b*dEdx + ls_d*dFdx
        Gt      = G       + ls_e*E      + ls_f*F
        dGtdx   = dGdx    + ls_e*dEdx   + ls_f*dFdx

    !            print*, 'IN xformsoln E, F, G, dEdx, dFdx, dGdx, ls_a , ls_b , ls_c , ls_d , ls_e ,ls_f', &
    !                    & E, F, G, dEdx, dFdx, dGdx, ls_a , ls_b , ls_c , ls_d , ls_e ,ls_f
    !            print*, 'OUT xformsoln Et, Ft, Gt, dEtdx, dFtdx, dGtdx', Et, Ft, Gt, dEtdx, dFtdx, dGtdx

    END SUBROUTINE xformsoln

    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE xformsoln_PO4_M(dum_EFPQ_P, dum_EFPQ_M, dum_dEFPQdz_P, dum_dEFPQdz_M, &
    dum_g_P, dum_g_M, dum_dgdz_P, dum_dgdz_M, dum_mat_C, dum_vec_D, &
    loc_EFPQ_P_t, loc_G_P_t, loc_dEFPQ_P_t, loc_dG_P_t, loc_EFPQ_M_t, loc_G_M_t, &
    loc_dEFPQ_M_t, loc_dG_M_t)

        real, dimension (1:4), intent(in) :: dum_EFPQ_P, dum_dEFPQdz_P, dum_EFPQ_M, dum_dEFPQdz_M
        real,INTENT(IN):: dum_g_P, dum_g_M, dum_dgdz_P, dum_dgdz_M
        real, dimension (4, 4), INTENT(in) :: dum_mat_C
        real, dimension (1:4), INTENT(in) ::  dum_vec_D

        ! output variables
        real, dimension (1:4), intent(inout) :: loc_EFPQ_P_t, loc_dEFPQ_P_t, loc_EFPQ_M_t, loc_dEFPQ_M_t
        real,intent(inout):: loc_G_P_t, loc_dG_P_t, loc_G_M_t, loc_dG_M_t

        ! Find 'transformed' soln such that in layer l,
        !    y_l = A_r*et + B_r*ft + gt
        ! (ie l soln written in terms of r solution coefficents A_r, B_r)
        !
        ! here save values in matrices - as this saves a lot of code

        loc_EFPQ_P_t = matmul(transpose(dum_mat_C), dum_EFPQ_P)     ! DH TODO: check multiplication with transpose, especially if vector*vector
        loc_dEFPQ_P_t = matmul(transpose(dum_mat_C), dum_dEFPQdz_P)
        loc_G_P_t = dot_product(dum_vec_D, dum_EFPQ_P)+dum_g_P
        loc_dG_P_t = dot_product(dum_vec_D,dum_dEFPQdz_P)+dum_dgdz_P

        loc_EFPQ_M_t = matmul(transpose(dum_mat_C), dum_EFPQ_M)
        loc_dEFPQ_M_t = matmul(transpose(dum_mat_C), dum_dEFPQdz_M)
        loc_G_M_t = dot_product(dum_vec_D, dum_EFPQ_M) + dum_g_M
        loc_dG_M_t = dot_product(dum_vec_D,dum_dEFPQdz_M) + dum_dgdz_M

!        print*,' '
!        print*, 'IN xformsoln_PO4_M '
!        print*, 'loc_EFPQ_P_t', loc_EFPQ_P_t
!        print*, 'loc_dEFPQ_P_t', loc_dEFPQ_P_t
!        print*, 'loc_G_P_t', loc_G_P_t
!        print*, 'loc_dG_P_t', loc_dG_P_t
!        print*,' '
!        print*, 'loc_EFPQ_M_t', loc_EFPQ_M_t
!        print*, 'loc_dEFPQ_M_t', loc_dEFPQ_M_t
!        print*, 'loc_G_M_t', loc_G_M_t
!        print*, 'loc_dG_M_t', loc_dG_M_t

    END SUBROUTINE xformsoln_PO4_M

    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE solve2eqn(a, b, c, d, e, f, x, y)

        ! Find soln of
        ! | a    b |  |x|   = | e |
        ! | c    d |  |y|     | f |

        real,INTENT(IN)::a, b, c, d, e, f
        real,INTENT(OUT)::x, y

        ! local variable
        real::det

        det = a*d-b*c
        x    =  (e*d-b*f)/det
        y    =  (a*f-e*c)/det

    !    print*,'a, b, c, d, e, f', a, b, c, d, e, f
    !    print*,'det', det
    !    print*,'x', x
    !    print*,'y', y

    END SUBROUTINE solve2eqn

    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE solve2eqn_PO4_M(dum_mat_X, dum_vec_Y, dum_A, dum_B, dum_C, dum_dim)

        ! Find solution of
        ! | x1  .  . x4|  |A|     | y1 |
        ! |     .      |  |B|     | y2 |
        ! |       .    |  |C|   = | y3 |
        ! | .       x16|  |D|     | y4 |
        integer, intent(in) :: dum_dim                             ! dim of input matrix to inverse
        real, dimension (1:dum_dim,1:dum_dim), intent(in) :: dum_mat_X
        real, dimension (1:dum_dim), intent(in) ::  dum_vec_Y
        real,INTENT(INOUT) :: dum_A, dum_B, dum_C

        ! local variable
        real, dimension (1:dum_dim,1:dum_dim) :: loc_mat_X, loc_inv_mat_X
        real, dimension (1:dum_dim) :: loc_vec_Z
        ! save matrix locally, as the original matrix dum_mat_X(4,4) will be destroyed during the calculation
        loc_mat_X = dum_mat_X
        ! calculate loc_mat_X^{-1}
        call inverse(loc_mat_X,loc_inv_mat_X,dum_dim)
        loc_vec_Z = matmul(loc_inv_mat_X, dum_vec_Y)

        dum_A = loc_vec_Z(1)
        dum_B = loc_vec_Z(2)
        dum_C = loc_vec_Z(3)

!        print*,' IN SOLVE2EQN_PO4_M '
!        print*,'dum_A ', dum_A
!        print*,'dum_B', dum_B
!        print*,'dum_C', dum_C

    END SUBROUTINE solve2eqn_PO4_M


    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    FUNCTION FUN_zO2(z)

        real FUN_zO2, z, flxzox, conczox, flxswi, r_zxf

        !    print*,' '
        !    print*,'..... START FUN_zO2'

        call zO2_calcbc(z, 1, flxzox, conczox, flxswi, r_zxf)

        FUN_zO2 = flxzox + calcFO2(z)

    !    print*,'FUN_zO2, flxzox, calcFO2(z)', FUN_zO2, flxzox, calcFO2(z)
    !    print*,' '
    END FUNCTION FUN_zO2

    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    FUNCTION FUN_zNO3(z)

        real FUN_zNO3, z, flxzno3, conczno3, flxswi, r_zxf

        !    print*,' '
        !    print*,'..... START FUN_zNO3'

        call zNO3_calcbc(z, 1, flxzno3, conczno3, flxswi)

        FUN_zNO3 = -flxzno3

    !    print*,'FUN_zNO3', FUN_zNO3
    !    print*,' '
    END FUNCTION FUN_zNO3

    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    FUNCTION FUN_zSO4(z)

        real FUN_zSO4, z, flxzso4, conczso4, flxswi, r_zxf

        !    print*,' '
        !    print*,'..... START FUN_zSO4'

        call zSO4_calcbc(z, 1, flxzso4, conczso4, flxswi)

        FUN_zSO4 = -flxzso4 - calcFSO4(z)

    !    print*,'FUN_zSO4, flxzso4, calcFSO4(z)', FUN_zSO4, flxzso4, calcFSO4(z)
    !    print*,' '
    END FUNCTION FUN_zSO4


    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------


    subroutine inverse(a,c,n)
        !============================================================
        ! Inverse matrix
        ! Method: Based on Doolittle LU factorization for Ax=b
        ! Alex G. December 2009
        !-----------------------------------------------------------
        ! input ...
        ! a(n,n) - array of coefficients for matrix A
        ! n      - dimension
        ! output ...
        ! c(n,n) - inverse matrix of A
        ! comments ...
        ! the original matrix a(n,n) will be destroyed
        ! during the calculation
        !===========================================================
        implicit none
        integer n
        real, dimension (1:n,1:n) :: a, c
        real, dimension (1:n,1:n) :: L, U
        real, dimension (1:n) :: b, d, x
        real coeff
!        double precision a(n,n), c(n,n)
!        double precision L(n,n), U(n,n), b(n), d(n), x(n)
!        double precision coeff
        integer i, j, k

        ! step 0: initialization for matrices L and U and b
        ! Fortran 90/95 aloows such operations on matrices
        L=0.0
        U=0.0
        b=0.0

        ! step 1: forward elimination
        do k=1, n-1
            do i=k+1,n
                coeff=a(i,k)/a(k,k)
                L(i,k) = coeff
                do j=k+1,n
                    a(i,j) = a(i,j)-coeff*a(k,j)
                end do
            end do
        end do

        ! Step 2: prepare L and U matrices
        ! L matrix is a matrix of the elimination coefficient
        ! + the diagonal elements are 1.0
        do i=1,n
            L(i,i) = 1.0
        end do
        ! U matrix is the upper triangular part of A
        do j=1,n
            do i=1,j
                U(i,j) = a(i,j)
            end do
        end do

        ! Step 3: compute columns of the inverse matrix C
        do k=1,n
            b(k)=1.0
            d(1) = b(1)
            ! Step 3a: Solve Ld=b using the forward substitution
            do i=2,n
                d(i)=b(i)
                do j=1,i-1
                    d(i) = d(i) - L(i,j)*d(j)
                end do
            end do
            ! Step 3b: Solve Ux=d using the back substitution
            x(n)=d(n)/U(n,n)
            do i = n-1,1,-1
                x(i) = d(i)
                do j=n,i+1,-1
                    x(i)=x(i)-U(i,j)*x(j)
                end do
                x(i) = x(i)/u(i,i)
            end do
            ! Step 3c: fill the solutions x(n) into column k of C
            do i=1,n
                c(i,k) = x(i)
            end do
            b(k)=0.0
        end do
    end subroutine inverse

    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    !!!! TODO: better put outside module, in kind of collection of auxiliary functions

    FUNCTION zbrent(func,x1,x2,tol)

        ! calculate root of func in the interval [x1,x2]

        INTEGER ITMAX
        real zbrent,tol,x1,x2,func,EPS
        EXTERNAL func
        PARAMETER (ITMAX=100,EPS=3.e-8)
        INTEGER iter
        real a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm

        !    print*,' '
        !    print*,'++++++++++++ START zbrent ++++++++++++++++ '

        a=x1
        b=x2
        fa=func(a)
        fb=func(b)
        !was      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))pause
        !was     *'root must be bracketed for zbrent'
        IF((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))THEN
            print*,'root must be bracketed for zbrent'
            STOP
        ELSE
            c=b
            fc=fb
            do 11 iter=1,ITMAX
                if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
                    c=a
                    fc=fa
                    d=b-a
                    e=d
                endif
                if(abs(fc).lt.abs(fb)) then
                    a=b
                    b=c
                    c=a
                    fa=fb
                    fb=fc
                    fc=fa
                endif
                tol1=2.*EPS*abs(b)+0.5*tol
                xm=.5*(c-b)
                if(abs(xm).le.tol1 .or. fb.eq.0.)then
                    zbrent=b
                    return
                endif
                if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
                    s=fb/fa
                    if(a.eq.c) then
                        p=2.*xm*s
                        q=1.-s
                    else
                        q=fa/fc
                        r=fb/fc
                        p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
                        q=(q-1.)*(r-1.)*(s-1.)
                    endif
                    if(p.gt.0.) q=-q
                    p=abs(p)
                    if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
                        e=d
                        d=p/q
                    else
                        d=xm
                        e=d
                    endif
                else
                    d=xm
                    e=d
                endif
                a=b
                fa=fb
                if(abs(d) .gt. tol1) then
                    b=b+d
                else
                    b=b+sign(tol1,xm)
                endif
                fb=func(b)
11          continue
        END IF
        ! was      pause 'zbrent exceeding maximum iterations'
        print*,'zbrent exceeding maximum iterations'
        zbrent=b
        STOP

        return

    END FUNCTION zbrent


end module
