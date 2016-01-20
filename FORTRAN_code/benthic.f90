!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!     Sediment model - to be named (not called SEDGEM ;-) )
!!!     Authors: Sandra Arndt, Dominik Hülse, Stuart Daines
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
        real*8 rho_sed                            ! sediment density (g/cm3)
        real*8 wdepth                            ! water depth (m)
        real*8 w                                      ! burial velocity  (cm/yr)
        real*8 z0, zox                                ! surface
        real*8 zbio                                ! bioturbation depth (cm)

        real*8 zinf                               !Inifinity (cm)
        !zinf = 1000
        !zlow=100
        real*8 Dbio                                 !bioturbation coefficient (cm2/yr)
        real*8 por                                !porosity (-)
        real*8 tort                              !tortuosity (-)
        real*8 irrigationFactor                   !irrigation factor (-)
        real*8 dispFactor                             !dispersion factor (-)

        !stoichiometric factors
        real*8 SD                                     !volume factor solid->dissolved phase
        real*8 OC                                     !O2/C (mol/mol)
        real*8 NC1                                    !N/C first TOC fraction (mol/mol)
        real*8 NC2                                    !N/C second TOC fraction (mol/mol)
        real*8 PC1                                    !P/C first TOC fraction (mol/mol)
        real*8 PC2                                    !P/C second TOC fraction (mol/mol)
        real*8 SO4C                                   !SO4/C (mol/mol)
        real*8 DICC1                                  !DIC/C until zSO4 (mol/mol)
        real*8 DICC2                                  !DIC/C below zSO$ (mol/mol)
        real*8 MC                                     !CH4/C (mol/mol)
        real*8 gamma                                !fraction of NH4 that is oxidised in oxic layer
        real*8 gammaH2S                           !fraction of H2S that is oxidised in oxic layer
        real*8 gammaCH4                           !fraction of CH4 that is oxidised at SO4
        real*8 satSO4                               ! SO4 saturation
        real*8 NO3CR                                  ! NO3 consumed by Denitrification

        real*8 zoxgf                            ! cm, rolloff NH4, H2S oxidation for small zox depth

        !bottom water concentrations
        real*8 T                                                     !temperature (degree C)
        real*8 C01                                         !TOC concentration at SWI (wt!) -> (mol/cm3 bulk phase)
        real*8 C02                                        !TOC concentration at SWI (wt!) -> (mol/cm3 bulk phase)
        real*8 O20                                               !O2  concentration at SWI (mol/cm3)
        real*8 SO40                                             !SO4 concentration at SWI (mol/cm3)
            !SO40 = 2800e-9
        real*8 H2S0                                               !H2S concentration at SWI (mol/cm3)
        real*8 DIC0                                             !DIC concentration at SWI (mol/cm3)
        real*8 ALK0                                             !ALK concentration at SWI (mol/cm3)
        real*8 NO30                                               !NO3 concentration at SWI (mol/cm3)
        real*8 NH40                                               !NH4 concentration at SWI (mol/cm3)
        real*8 PO40                                                  !PO4 concentration at SWI (mol/cm3)
        real*8 S0                                                     !Salinity at SWI


        ! ORGANIC MATTER
        real*8 DC1                                           !TOC diffusion coefficient (cm2/yr)
        real*8 C, C1, C2
        real*8 k1                                             !TOC degradation rate constnat (1/yr)
        real*8 k2                                              !TOC degradation rate constant (1/yr)

        ! O2
        real*8 qdispO2                          !O2 diffusion coefficient in water (cm2/yr)
        real*8 adispO2                          !O2 linear coefficient for temperature dependence (cm2/yr/oC)
        real*8 DO21                             !O2 diffusion coefficient in bioturbated layer (cm2/yr)
        real*8 DO22                             !O2 diffusion coefficient in non-bioturbated layer (cm2/yr)
        real*8 r_zxf                            !roll off oxidation at low zox

        real*8 aa11, bb11, aa21, A11, A21, aa12, bb12, aa22, A12, A22
        real*8 ls_a, ls_b, ls_c, ls_d, ls_e, ls_f

        ! Nitrate (NO3)
        real*8 qdispNO3                 ! NO3 diffusion coefficient in water (cm2/yr)
        real*8 adispNO3                 ! NO3 linear coefficient for temperature dependence (cm2/yr/oC)
        real *8 DN1                     ! NO3 diffusion coefficient in bioturbated layer (cm2/yr)
        real*8 DN2                      ! NO3 diffusion coefficient in non-bioturbated layer (cm2/yr)
        real*8 zno3

        ! Sulfate (SO4)
        real*8 qdispSO4                 ! SO4 diffusion coefficient in water (cm2/yr)
        real*8 adispSO4                 ! SO4 linear coefficient for temperature dependence (cm2/yr/oC)
        real*8 DSO41                    ! SO4 diffusion coefficient in bioturbated layer (cm2/yr)
        real*8 DSO42                    ! SO4 diffusion coefficient in non-bioturbated layer (cm2/yr)
        real*8 zso4

        ! Ammonium (NH4)
        real*8 qdispNH4                 ! NH4 diffusion coefficient in water (cm2/yr)
        real*8 adispNH4                 ! NH4 linear coefficient for temperature dependence (cm2/yr/oC)
        real*8 DNH41                    ! NH4 diffusion coefficient in bioturbated layer (cm2/yr)
        real*8 DNH42                    ! NH4 diffusion coefficient in non-bioturbated layer (cm2/yr)

        ! Hydrogen sulfide (H2S)
        real*8 qdispH2S                 ! H2S diffusion coefficient in water (cm2/yr)
        real*8 adispH2S                 ! H2S linear coefficient for temperature dependence (cm2/yr/oC)
        real*8 DH2S1                    ! H2S diffusion coefficient in bioturbated layer (cm2/yr)
        real*8 DH2S2                    ! H2S diffusion coefficient in non-bioturbated layer (cm2/yr)

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

        rho_sed=2.5                           ! sediment density (g/cm3)
        wdepth=600                            ! water depth (m)
        z0  = 0.0                               ! surface
        zox = 0.0
        zbio=5.0                                ! bioturbation depth (cm)

        zinf=100.0                               !Inifinity (cm)
        !zinf = 1000
        !zlow=100
        Dbio=3                                 !bioturbation coefficient (cm2/yr)
        por=0.8                                !porosity (-)
        tort=3.0                               !tortuosity (-)
        irrigationFactor=1.0

        gamma=1                                !fraction of NH4 that is oxidised in oxic layer
        gammaH2S=0.8                           !fraction of H2S that is oxidised in oxic layer
        gammaCH4=0.9                           !fraction of CH4 that is oxidised at SO4
        satSO4=0                               ! SO4 saturation

        zoxgf = 0.1                            ! cm, rolloff NH4, H2S oxidation for small zox depth

        !bottom water concentrations
        T=20.0                                                     ! temperature (degree C)
        C01=0.2*1e-2/12*rho_sed                                ! TOC concentration at SWI (wt!) -> (mol/cm3 bulk phase)
        C02=0.2*1e-2/12*rho_sed                                ! TOC concentration at SWI (wt!) -> (mol/cm3 bulk phase)
        O20=300.0e-9         !was  300.0e-9                                     ! O2  concentration at SWI (mol/cm3)
        SO40=28000.0e-9      !  28000.0e-9                                     ! SO4 concentration at SWI (mol/cm3)
        H2S0=0.0e-9           !was 0.0e-9                                 ! H2S concentration at SWI (mol/cm3)
        DIC0=2000.0e-9                                             ! DIC concentration at SWI (mol/cm3)
        ALK0=2400.0e-9                                             ! ALK concentration at SWI (mol/cm3)
        NO30=30.0e-9                                               ! NO3 concentration at SWI (mol/cm3)
        NH40=1.0e-9                                                ! NH4 concentration at SWI (mol/cm3)
        PO40=1e-9                                                  ! PO4 concentration at SWI (mol/cm3)
        S0=35                                                      ! Salinity at SWI


        w=10.0**(-0.87478367-0.00043512*wdepth)*3.3             ! sedimentation rate, cm/yr
        dispFactor=por**(tort-1.0)*irrigationFactor             !dispersion factor (-) - Ausbreitung - type of mixing that accompanies hydrodynamic flows -> ~builds out paths of flow
        SD=(1-por)/por                                          !volume factor solid->dissolved phase
        OC=1.0*SD                                               ! O2/C (mol/mol)
        NC1=0.1509*SD                                           ! N/C first TOC fraction (mol/mol)
        NC2=0.13333*SD                                          ! N/C second TOC fraction (mol/mol)
        SO4C=0.5*SD                                             ! SO4/C (mol/mol)
        PC1=0.0094*SD                                          ! P/C first TOC fraction (mol/mol)
        PC2=0.0094*SD                                          ! P/C second TOC fraction (mol/mol)
        DICC1=1.0*SD                                           ! DIC/C until zSO4 (mol/mol)
        DICC2=0.5*SD                                           ! DIC/C below zSO$ (mol/mol)
        MC=0.5*SD                                              ! CH4/C (mol/mol)
        NO3CR=(94.4/106)*SD                                    ! NO3 consumed by Denitrification

        ! ORGANIC MATTER
        DC1 = Dbio
        k1=0.01
        k2=0.001

        ! O2
        qdispO2=348.5750
        adispO2=14.0890

        DO21=(qdispO2+adispO2*T)*dispFactor+Dbio
        DO22=(qdispO2+adispO2*T)*dispFactor

        r_zxf=0.0

        ! Nitrate (NO3)
        qdispNO3=308.4221
        adispNO3=12.2640
        DN1=(qdispNO3+adispNO3*T)*dispFactor+Dbio
        DN2=(qdispNO3+adispNO3*T)*dispFactor
        zno3 = 0.0

        ! Sulfate (SO4)
        qdispSO4=309.0528                               ! SO4 diffusion coefficient in water (cm2/yr)
        adispSO4=12.2640                                ! SO4 linear coefficient for temperature dependence (cm2/yr/oC)
        DSO41=(qdispSO4+adispSO4*T)*dispFactor+Dbio     ! SO4 diffusion coefficient in bioturbated layer (cm2/yr)
        DSO42=(qdispSO4+adispSO4*T)*dispFactor          ! SO4 diffusion coefficient in non-bioturbated layer (cm2/yr)
        zso4 = 0.0

        ! Ammonium (NH4)
        qdispNH4=309.0528
        adispNH4=12.2640
        DNH41=(qdispNH4+adispNH4*T)*dispFactor+Dbio
        DNH42=(qdispNH4+adispNH4*T)*dispFactor

        ! Hydrogen sulfide (H2S)
        qdispH2S=309.0528
        adispH2S=12.2640
        DH2S1=(qdispH2S+adispH2S*T)*dispFactor+Dbio
        DH2S2=(qdispH2S+adispH2S*T)*dispFactor

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

!   organic matter burial

! local variables

!    real*8 aa11, bb11, aa21, A11, A21, aa12, bb12, aa22, A12, A22
    real*8 dC1dz, C1flx, dC2dz, C2flx, Cflx             ! Cflx: Sed input flux to upper boundary, per cm^2 water column
    real*8 F_TOC1, F_TOC2, F_TOC                        ! Flux through lower boundary zinf, per cm^2 water-column


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

    end SUBROUTINE zTOC


! ****************************************************************************************************************************** !
!   *****************************************************************
!   *****************************************************************
!------------------------------------------------------------------------------------



    FUNCTION calcReac(zU, zL, reac1, reac2)

        real*8,intent(in):: zU, zL, reac1, reac2
        real*8 calcReac

    ! Integral of reacted organic matter from zU to zL,
    ! multiplied by stoichiometric factors reac1, reac2 (for the two OC phases)

    ! Vector-friendly way of handling 3 cases:
    ! 1) wholly within bioturbated layer:    calcReac_l1(zU,zL)     + (0 =) calcReac_l2(bsd.zbio, bsd.zbio)
    ! 2) wholly within non-bio     layer:  (0=) calcReac_l1(zbio, zbio) +   calcReac_l2(zU, zL)
    ! 3) crossing zbio                       calcRead_l1(zU,zbio)   +       calcReac_l2(zbio, zL)

    calcReac = calcReac_l1(min(zU,zbio), min(zL,zbio), reac1, reac2) &
           & + calcReac_l2(max(zU,zbio), max(zL, zbio), reac1, reac2)

   ! TODO confirm (1-bsd.por)*  has been added (to k1 & k2 ?)

    END FUNCTION calcReac

! ****************************************************************************************************************************** !
!   *****************************************************************
!   *****************************************************************
!------------------------------------------------------------------------------------

    FUNCTION calcReac_l1(zU, zL, reac1, reac2)

    real*8,intent(in):: zU, zL, reac1, reac2
    real*8 calcReac_l1, reacf1, reacf2

    reacf1 = k1*reac1
    reacf2 = k2*reac2
    calcReac_l1 = -reacf1*(A11*(exp(aa11*zU)*bb11 - exp(bb11*zU)*aa11 - exp(aa11*zL)*bb11 + exp(bb11*zL)*aa11) &
                & + C01*exp(bb11*zU)*aa11 - C01*exp(bb11*zL)*aa11)/(aa11*bb11) &
                & -reacf2*(A12*(exp(aa12*zU)*bb12 - exp(bb12*zU)*aa12 - exp(aa12*zL)*bb12 + exp(bb12*zL)*aa12) &
                & + C02*exp(bb12*zU)*aa12 - C02*exp(bb12*zL)*aa12)/(aa12*bb12)


    END FUNCTION calcReac_l1

! ****************************************************************************************************************************** !
!   *****************************************************************
!   *****************************************************************
!------------------------------------------------------------------------------------

    FUNCTION calcReac_l2(zU, zL, reac1, reac2)

    real*8,intent(in):: zU, zL, reac1, reac2
    real*8 calcReac_l2, reacf1, reacf2

    reacf1 = k1*reac1
    reacf2 = k2*reac2

    calcReac_l2 = -reacf1*A21*(exp(aa21*zU) - exp(aa21*zL))/aa21 &
                & -reacf2*A22*(exp(aa22*zU) - exp(aa22*zL))/aa22

!    print*,'calcReac_l2', calcReac_l2

    END FUNCTION calcReac_l2



! ****************************************************************************************************************************** !
!   *****************************************************************
!   *****************************************************************
!------------------------------------------------------------------------------------



    SUBROUTINE calcfg_l12(z, reac1, reac2, ktemp, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ls_D1, ls_D2,&
                            & ltype, e, dedz, f, dfdz, g, dgdz)
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

        real*8,intent(in)::       z, reac1, reac2, ktemp
        INTEGER, INTENT(in)::   ltype
        real*8,intent(inout)::    e, dedz, f, dfdz, g, dgdz
        real*8                    ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ls_D1, ls_D2

        ! local variables
        real*8 e_1, f_1, g_1, dedz_1, dfdz_1, dgdz_1
        !real*8 ls_a, ls_b, ls_c, ls_d, ls_e, ls_f


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
                                        &  e, f, g, dedz, dfdz, dgdz)
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

            real*8 z, reac1, reac2, Dtemp, ktemp                                                          ! in from SUBROUTINE before
            real*8,INTENT(inout)::  e, dedz, f, dfdz, g, dgdz             ! out

            ! local variables
            real*8 b1, pfac, PhiI1, PhiII1, PhiIII1, PhiI2, PhiII2, PhiIII2
            real*8 ea11z, eb11z, ea12z, eb12z


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
                & PhiI2*ea12z + PhiII2*eb12z + PhiIII2*eb12z


            dgdz = PhiI1*aa11*ea11z + PhiII1*bb11*eb11z + PhiIII1*bb11*eb11z + &
                  & PhiI2*aa12*ea12z + PhiII2*bb12*eb12z + PhiIII2*bb12*eb12z

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

            real*8 z, reac1, reac2, Dtemp, ktemp                                                          ! in from SUBROUTINE before
            real*8,INTENT(inout)::  e, dedz, f, dfdz, g, dgdz             ! out

            ! local variables
            real*8 b2, pfac, PhiI1, PhiI2
            real*8 ea11z, eb11z, ea12z, eb12z

            e = 1.0
            dedz = 0.0
            b2 = w/Dtemp
            f=exp(z*b2)
            dfdz = b2*exp(z*b2)

            !pfac=1./bsd.por;   ! assume org matter already *(1-bsd.por)
            pfac = 1            !in fact, already has (1-por)/por



            PhiI1 = -pfac*k1*(reac1)*A21/(Dtemp*aa21**2-w*aa21-ktemp)
            PhiI2 = -pfac*k2*(reac2)*A22/(Dtemp*aa22**2-w*aa22-ktemp)

            g = PhiI1*exp(aa21*z) + PhiI2*exp(aa22*z)
            dgdz = PhiI1*aa21*exp(aa21*z) + PhiI2*aa22*exp(aa22*z)

!            print*, 'IN  calcfg_l2 g', g
!            print*, 'IN  calcfg_l1 dgdz', dgdz

        end SUBROUTINE calcfg_l2


!------------------------------------------------------------------------------------
!   *****************************************************************
!   *****************************************************************
!------------------------------------------------------------------------------------

    SUBROUTINE prepfg_l12(reac1, reac2, ktemp, zU, zL, D1, D2, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ltype)
        real*8 reac1, reac2, ktemp, zU, zL, D1, D2
        real*8 e_zbio_l1, dedz_zbio_l1, f_zbio_l1, dfdz_zbio_l1, g_zbio_l1, dgdz_zbio_l1
        real*8 e_zbio_l2, dedz_zbio_l2, f_zbio_l2, dfdz_zbio_l2, g_zbio_l2, dgdz_zbio_l2
        real*8 ls_a, ls_b, ls_c, ls_d, ls_e, ls_f

        INTEGER,INTENT(inout):: ltype

        if(zL <= zbio)then  ! wholly within bioturbated layer
            ltype = 1
        elseif(zU >= zbio) then ! wholly within non-bioturbated layer
            ltype = 2
        else             ! crossing boundary - sort out solution matching at zbio
!            print*, 'IN  CROSS BOUNDARY CASE '
            ltype = 3
            call calcfg_l1(zbio, reac1, reac2, D1, ktemp, e_zbio_l1, dedz_zbio_l1, f_zbio_l1, dfdz_zbio_l1, g_zbio_l1, dgdz_zbio_l1)
            call calcfg_l2(zbio, reac1, reac2, D2, ktemp, e_zbio_l2, dedz_zbio_l2, f_zbio_l2, dfdz_zbio_l2, g_zbio_l2, dgdz_zbio_l2)

            ! match solutions at zbio - continuous concentration and flux
            call matchsoln(e_zbio_l1, f_zbio_l1, g_zbio_l1, D1*dedz_zbio_l1, D1*dfdz_zbio_l1, D1*dgdz_zbio_l1, &
                         & e_zbio_l2, f_zbio_l2, g_zbio_l2, D2*dedz_zbio_l2, D2*dfdz_zbio_l2, D2*dgdz_zbio_l2, &
                                0.0D00, 0.0D00, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f)
!           print*, 'in prepfg_l12 AFTER matchsoln:  ls_a, ls_b, ls_c, ls_d, ls_e, ls_f', ls_a, ls_b, ls_c, ls_d, ls_e, ls_f
        end if

    END SUBROUTINE prepfg_l12



!------------------------------------------------------------------------------------
!   *****************************************************************
!   *****************************************************************

!                               Oxygen

!   *****************************************************************
!   *****************************************************************
!------------------------------------------------------------------------------------

    SUBROUTINE zO2()

    real*8 flxzox, conczox, flxswi, fun0, fun1, zL, tol
    integer bctype
    real*8 FO2

!    print*, ''
!    print*, ''
!    print*, '---------------------- START zO2 ------------------------ '

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
!    print*,'flxzox', flxzox
!    print*,'conczox', conczox
!    print*,'flxswi', flxswi

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
    real*8, intent(in)::zox
    integer, intent(in)::bctype
    real*8, intent(inout)::flxzox, conczox, flxswi, r_zxf
    !   local variables

!    real*8 qdispO2                          !O2 diffusion coefficient in water (cm2/yr)
!    real*8 adispO2                          !O2 linear coefficient for temperature dependence (cm2/yr/oC)
!    real*8 DO21                             !O2 diffusion coefficient in bioturbated layer (cm2/yr)
!    real*8 DO22                             !O2 diffusion coefficient in non-bioturbated layer (cm2/yr)
    real*8 reac1, reac2, ls !, z0, zox
    integer ltype
    real*8 ls_a, ls_b, ls_c, ls_d, ls_e, ls_f
    real*8 e_0, dedz_0, f_0, dfdz_0, g_0, dgdz_0
    real*8 e_zox, dedz_zox, f_zox, dfdz_zox, g_zox, dgdz_zox
!    real*8 bctype1_AO2, bctype1_BO2, bctype2_AO2, bctype2_BO2

    real*8 rO2_AO2, rO2_BO2, Dzox, Dswi

!    real*8 FO2

!    qdispO2=348.5750
!    adispO2=14.0890

!    DO21=(qdispO2+adispO2*T)*dispFactor+Dbio
!    DO22=(qdispO2+adispO2*T)*dispFactor


!   reactive terms: OM degradation (-) and nitrification (-)
    reac1=-OC-2*gamma*NC1
    reac2=-OC-2*gamma*NC2

!    print*, ''
!    print*, '------- START zO2_calcbc ----- zox:', zox

    ! calculate solution for given zox

!    print*, 'Preparation: sort out solution-matching across bioturbation boundary (if necessary)'
    ! Preparation: sort out solution-matching across bioturbation boundary (if necessary)
    call prepfg_l12(reac1, reac2, 0.0D00, z0, zox, DO21, DO22, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ltype)

    ! basis functions at upper boundary
    call calcfg_l12(z0, reac1, reac2, 0.0D00, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, DO21, DO22, ltype, &
            &e_0, dedz_0, f_0, dfdz_0, g_0, dgdz_0)

    ! ... and lower boundary
    call calcfg_l12(zox, reac1, reac2, 0.0D00, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, DO21, DO22, ltype, e_zox, dedz_zox,&
            & f_zox, dfdz_zox, g_zox, dgdz_zox)

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
        print*,' '
        print*,' '
        print*,'------ zO2_calcbc --------- flxzox is INFFFFFFFFFFFFFFFF', flxzox
        print*,' '
        print*,' '
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

    real*8 calcFO2, z, tmpreac1, tmpreac2

    ! Oxydation of reduced species at zox (NEED A RATIO for ODU! and add NH4
    ! adsporption!
    !tmpreac1=0.5*bsd.OC+2*bsd.gamma*bsd.NC1;
    !tmpreac2=0.5*bsd.OC+2*bsd.gamma*bsd.NC2;

!    print*,' '
!    print*,'..... START calcFO2'

    tmpreac1=OC+2*gamma*NC1
    tmpreac2=OC+2*gamma*NC2
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

    real*8 flxzinf, conczinf, flxswi, fun0, flxzno3, conczno3, flxswiNO3, zL, tol
    integer bctype
    real*8 FNO3

!    print*, ''
!    print*, '------------------------------------------------------------------'
!    print*, '---------------------- START zNO3 ------------------------------- '

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

    real*8 zNO3, flxzNO3, conczNO3, flxswi
    integer bctype, ltype1, ltype2
    real*8 ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1
    real*8 ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2
    real*8 e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3
    real*8 e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox
    real*8 e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox
    real*8 zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
    real*8 e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00
    real*8 e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0
    real*8 bctype1_A2, bctype1_B2, bctype2_A2, bctype2_B2
    real*8 rNO3_A2, rNO3_B2, rNO3_A1, rNO3_B1

    real*8 FNH4

    ! Calculate trial solution for given zno3, matching boundary conditions from layer-by-layer solutions


    ! Preparation: for each layer, sort out solution - matching across bioturbation boundary (if necessary)
    ! layer 1: 0 < z < zox, nitrification
    !    prepfg_l12(reac1, reac2,  ktemp ,     zU , zL , D1,  D2)
    call prepfg_l12(gamma*NC1, gamma*NC2,  0.0D00,  0.0D00, zox, DN1, DN2, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, ltype1)
!    print*, ''
!    print*, '1. prepfg_l12 RESULTS: ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, ltype1', ls_a1, ls_b1, ls_c1, &
!           & ls_d1, ls_e1, ls_f1, ltype1
   ! layer 2: zox < z < zno3, denitrification
    call prepfg_l12(-NO3CR, -NO3CR, 0.0D00,  zox, zNO3, DN1, DN2, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, ltype2)
!    print*, ''
!    print*, '2. prepfg_l12 RESULTS: ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, ltype2', ls_a2, ls_b2, ls_c2, ls_d2, &
!            & ls_e2, ls_f2, ltype2


    ! Work up from the bottom, matching solutions at boundaries
    ! Basis functions at bottom of layer 2 zno3
    call calcfg_l12(zno3, -NO3CR, -NO3CR, 0.0D00, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DN1, DN2, ltype2, &
            &e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3)

    ! Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from NH4 -> NO3 source)

    ! basis functions at bottom of layer 1
    call calcfg_l12(zox, gamma*NC1, gamma*NC2, 0.0D00, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DN1, DN2, ltype1, &
            & e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox)

    ! basis functions at top of layer 2

    call calcfg_l12(zox, -NO3CR, -NO3CR, 0.0D00, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DN1, DN2, ltype2, &
            & e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox)

    ! flux of NH4 to zox  TODO NH4 production by denitrification?


    FNH4 = calcReac(zno3, zinf, NC1, NC1)   ! MULTIPLY BY 1/POR ????
!    print*, 'FNH4 = ', FNH4

    ! match solutions at zox - continuous concentration, flux discontinuity from H2S ox

!    call matchsoln(e_zbio_l1, f_zbio_l1, g_zbio_l1, D1*dedz_zbio_l1, D1*dfdz_zbio_l1, D1*dgdz_zbio_l1, &
!                         & e_zbio_l2, f_zbio_l2, g_zbio_l2, D2*dedz_zbio_l2, D2*dfdz_zbio_l2, D2*dgdz_zbio_l2, &
!                                0.0D00, 0.0D00, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f)
    IF(zox .le. zbio)then
        call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                           & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                           & 0.0D00, -r_zxf*gamma*FNH4/DN1, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
!        print*, ''
!        print*,'zox<= zbio: RESULTS:  zox_a, zox_b, zox_c, zox_d, zox_e, zox_f', zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
    ELSE
        call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                           & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                           & 0.0D00, -r_zxf*gamma*FNH4/DN2, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
!        print*, ''
!        print*,'zox> zbio: RESULTS:  zox_a, zox_b, zox_c, zox_d, zox_e, zox_f', zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
    END IF

    ! Solution at swi, top of layer 1
    call calcfg_l12(0.0D00, gamma*NC1, gamma*NC2, 0.0D00, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DN1, DN2, ltype1, &
            & e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00)
!    print*, ''
!    print*,'top of layer 1; calcfg_l12: RESULTS:  e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00', e1_00, &
!            & dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00

    ! transform to use coeffs from l2
    call xformsoln(e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00, zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f, &
                                        & e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0)
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
    flxswi = por*DN1*(rNO3_A2*dedz1_0+rNO3_B2*dfdz1_0 + dgdz1_0)              ! NB: use A2, B2 as these are _xformed_ layer 1 basis functions

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
!    real*8, intent(in)::zox, zno3
!    real*8, intent(inout)::zso4

    ! local variables
    real*8 flxzso4, conczso4, flxswiso4, zL, tol
    integer bctype


!    print*, ''
!    print*, '------------------------------------------------------------------'
!    print*, '---------------------- START zSO4 ------------------------------- '


    ! Iteratively solve for zso4


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

    real*8, intent(in)::zso4
    integer, intent(in)::bctype
    real*8, intent(inout)::flxzso4, conczso4, flxswi !, r_zxf

    ! local variable
    integer ltype1, ltype2, ltype3
    real*8 reac1_so4, reac2_so4
    real*8 ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1
    real*8 ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2
    real*8 ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3
    real*8 e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4
    real*8 e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3
    real*8 e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3
    real*8 zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f
    real*8 e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox
    real*8 e2_zox0, f2_zox0, g2_zox0, dedz2_zox0, dfdz2_zox0, dgdz2_zox0
    real*8 e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox
    real*8 zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
    real*8 e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0
    real*8 e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00
    real*8 bctype1_A3, bctype1_B3, bctype2_A3, bctype2_B3
    real*8 rSO4_A1, rSO4_B1, rSO4_A2, rSO4_B2, rSO4_A3, rSO4_B3

    real*8 FH2S

!    print*, ' '
!    print*, '---------------------- START zSO4_calcbc ------------------------------- '

    reac1_so4=-SO4C
    reac2_so4=-SO4C

    ! Calculate trial solution for given zso4, matching boundary conditions from layer-by-layer solutions


    ! Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
    ! layer 1: 0 < z < zox, passive diffn
    !      ls =      prepfg_l12( bsd, swi, r, reac1,     reac2,     ktemp, zU, zL, D1,        D2)

    call prepfg_l12(0.0D00, 0.0D00, 0.0D00, 0.0D00, zox, DSO41, DSO42, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, ltype1)

    ! layer 2: zox < z < zno3, passive diffn
    call prepfg_l12(0.0D00, 0.0D00, 0.0D00, zox, zno3, DSO41, DSO42, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, ltype2)

    ! layer 3: zno3 < z < zso4, SO4 consumption by OM oxidation
            !rSO4.ls3 = r.zTOC.prepfg_l12(bsd, swi, r, obj.reac1, obj.reac2, 0, r.zno3, zso4, obj.DSO41, obj.DSO42)
    call prepfg_l12(reac1_so4, reac2_so4, 0.0D00, zno3, zso4, DSO41, DSO42, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, ltype3)


    ! Work up from the bottom, matching solutions at boundaries
    ! Basis functions at bottom of layer 3 zso4

    call calcfg_l12(zso4, reac1_so4, reac2_so4, 0.0D00, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, DSO41, DSO42, ltype3, &
            &e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4)

    ! Match at zno3, layer 2 - layer 3 (continuity and flux)
    ! basis functions at bottom of layer 2
    call calcfg_l12(zno3, 0.0D00, 0.0D00, 0.0D00, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DSO41, DSO42, ltype2, &
            &e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3)

    ! ... and top of layer 3
    call calcfg_l12(zno3, reac1_so4, reac2_so4, 0.0D00, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, DSO41, DSO42, ltype3, &
            &e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3)

    ! match solutions at zno3 - continuous concentration and flux
!    [zno3.a, zno3.b, zno3.c, zno3.d, zno3.e, zno3.f] = benthic_utils.matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, ...
!                                                                e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, ...
!                                                                0, 0);
    call matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, &
                           & e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, &
                           & 0.0D00, 0.0D00, zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f)


    ! Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from H2S source)
    ! flux of H2S to oxic interface (Source of SO4)
    ! NB: include methane region as AOM will produce sulphide as well..

    FH2S = calcReac(zno3, zso4, SO4C, SO4C) & ! MULTIPLY BY 1/POR ????
                  & + gammaCH4*calcReac(zso4, zinf, SO4C, SO4C)
!    print*,' '
!    print*,'FH2S ', FH2S

    ! basis functions at bottom of layer 1
    call calcfg_l12(zox, 0.0D00, 0.0D00, 0.0D00, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DSO41, DSO42, ltype1, &
            & e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox)

    ! basis functions at top of layer 2
    call calcfg_l12(zox, 0.0D00, 0.0D00, 0.0D00, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DSO41, DSO42, ltype2, &
            e2_zox0, dedz2_zox0, f2_zox0, dfdz2_zox0, g2_zox0, dgdz2_zox0)

    ! transform to use coeffs from l3
    call xformsoln(e2_zox0, f2_zox0, g2_zox0, dedz2_zox0, dfdz2_zox0, dgdz2_zox0, &
            & zno3_a , zno3_b , zno3_c , zno3_d , zno3_e ,zno3_f, &
            & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox )

    ! match solutions at zox - continuous concentration, flux discontinuity from H2S ox
    IF(zox .le. zbio)THEN
        call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                    & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                    & 0.0D00, -r_zxf*FH2S/DSO41, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
    ELSE
        call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                    & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                    & 0.0D00, -r_zxf*FH2S/DSO42, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
    END IF

    ! Solution at swi, top of layer 1
    call calcfg_l12(0.0D00, 0.0D00, 0.0D00, 0.0D00, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DSO41, DSO42, ltype1, &
            & e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00)

    ! transform to use coeffs from l3
    call xformsoln(e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00, zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f, &
            & e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0)

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
    !D = (zso4 <= bsd.zbio).*obj.DSO41 + (zso4 > bsd.zbio).*obj.DSO42
    IF(zso4 .le. zbio)THEN
        flxzso4 = DSO41*(rSO4_A3*dedz3_zso4+rSO4_B3*dfdz3_zso4 + dgdz3_zso4)        ! includes 1/por ie flux per (cm^2 pore area)
    ELSE
        flxzso4 = DSO42*(rSO4_A3*dedz3_zso4+rSO4_B3*dfdz3_zso4 + dgdz3_zso4)
    END IF

    ! flux at swi - DO include por so this is per cm^2 water column area
    flxswi = por*DSO41*(rSO4_A3*dedz1_0+rSO4_B3*dfdz1_0 + dgdz1_0)   ! NB: use A3, B3 as these are _xformed_ layer 1 basis functions

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

    real*8 calcFSO4, z, tmpreac1, tmpreac2

    ! Calculate SO4 consumption below zso4, by organic matter and indirectly via methane oxidation

!    print*,' '
!    print*,'..... START calcFSO4'

    tmpreac1    = SO4C*gammaCH4
    tmpreac2    = SO4C*gammaCH4

    calcFSO4 = calcReac(z, zinf, tmpreac1, tmpreac2)
    ! TODO confirm (1-bsd.por)*  has been added (to k1 & k2 ?)
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

!    real*8, intent(in)::zox, zno3

    ! local variables
    integer ltype1, ltype2, ltype3
    real*8 ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1
    real*8 ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2
    real*8 ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3
    real*8 e3_zinf, dedz3_zinf, f3_zinf, dfdz3_zinf, g3_zinf, dgdz3_zinf
    real*8 e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3
    real*8 e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3
    real*8 zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f
    real*8 e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox
    real*8 e2_zox0, dedz2_zox0, f2_zox0, dfdz2_zox0, g2_zox0, dgdz2_zox0
    real*8 e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox
    real*8 zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
    real*8 e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00
    real*8 e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0
    real*8 rNH4_A3, rNH4_B3
    real*8 rNH4_A2, rNH4_B2
    real*8 rNH4_A1, rNH4_B1

    real*8 flxswiNH4
    real*8 FNH4

!    print*, ''
!    print*, '------------------------------------------------------------------'
!    print*, '---------------------- START zNH4 ------------------------------- '

    ! Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
    ! layer 1: 0 < z < zox, NH4 prod (remaining after oxidation)
    !      ls =      prepfg_l12( bsd, swi, r, reac1,     reac2,     ktemp, zU, zL, D1,        D2)
    call prepfg_l12((1-gamma),(1-gamma),0.0D00, 0.0D00, zox, DNH41, DNH42, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, ltype1)

    ! layer 2: zox < z < zno3, passive diffn TODO NH4 from denitrification?
    call prepfg_l12(0.0D00, 0.0D00, 0.0D00, zox, zno3, DNH41, DNH42, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, ltype2)

    ! layer 3: zno3 < z < zinf, NH4 production
    call prepfg_l12(NC1, NC2, 0.0D00, zno3, zinf, DNH41, DNH42, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, ltype3)

    ! Work up from the bottom, matching solutions at boundaries
    ! Basis functions at bottom of layer 3 zinf
    call calcfg_l12(zinf, NC1, NC2, 0.0D00, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, DNH41, DNH42, ltype3, &
            & e3_zinf, dedz3_zinf, f3_zinf, dfdz3_zinf, g3_zinf, dgdz3_zinf)

    ! Match at zno3, layer 2 - layer 3 (continuity and flux)
    ! basis functions at bottom of layer 2
    call calcfg_l12(zno3, 0.0D00, 0.0D00, 0.0D00, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DNH41, DNH42, ltype2, &
                & e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3)

    ! ... and top of layer 3
    call calcfg_l12(zno3, NC1, NC2, 0.0D00, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, DNH41, DNH42, ltype3, &
            & e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3)

    ! match solutions at zno3 - continuous concentration and flux
    call matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, &
                    e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, &
                    0.0D00, 0.0D00, zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f)

    ! Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from NH4 sink)
    ! flux of NH4 to oxic interface  TODO NH4 prod by denitrification?
    FNH4 = calcReac(zno3, zinf, NC1, NC2)   ! MULTIPLY BY 1/POR ????

    ! basis functions at bottom of layer 1
    call calcfg_l12(zox, (1-gamma),(1-gamma) , 0.0D00, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DNH41, DNH42, ltype1, &
               & e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox)
    ! basis functions at top of layer 2
    call calcfg_l12(zox, 0.0D00, 0.0D00, 0.0D00, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DNH41, DNH42, ltype2, &
                    & e2_zox0, dedz2_zox0, f2_zox0, dfdz2_zox0, g2_zox0, dgdz2_zox0)

    ! transform to use coeffs from l3
    call xformsoln(e2_zox0, f2_zox0, g2_zox0, dedz2_zox0, dfdz2_zox0, dgdz2_zox0, &
                    zno3_a , zno3_b , zno3_c , zno3_d , zno3_e ,zno3_f, &
                    & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox)

    ! match solutions at zox - continuous concentration, flux discontinuity from NH4 ox
    IF(zox .le. zbio)THEN
        call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                      e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                     0.0D00, r_zxf*gamma*FNH4/DNH41, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
    ELSE
        call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                       e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                      0.0D00, r_zxf*gamma*FNH4/DNH42, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)

    END IF

    ! Solution at swi, top of layer 1
    call calcfg_l12(0.0D00, 0.0D00, 0.0D00, 0.0D00, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DNH41, DNH42, ltype1, &
                    & e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00)

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
    flxswiNH4 = por*DNH41*(rNH4_A3*dedz1_0+rNH4_B3*dfdz1_0 + dgdz1_0)   ! NB: use A3, B3 as these are _xformed_ layer 1 basis functions

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

!    real*8, intent(in)::zox, zso4

    ! local variables
    real*8 reac1_h2s, reac2_h2s                 ! reactive terms: OM degradation
    integer ltype1, ltype2, ltype3, ltype4
    real*8 ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1
    real*8 ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2
    real*8 ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3
    real*8 ls_a4, ls_b4, ls_c4, ls_d4, ls_e4, ls_f4
    real*8 e4_zinf, dedz4_zinf, f4_zinf, dfdz4_zinf, g4_zinf, dgdz4_zinf
    real*8 e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4
    real*8 e4_zso4, dedz4_zso4, f4_zso4, dfdz4_zso4, g4_zso4, dgdz4_zso4
    real*8 zso4_a, zso4_b, zso4_c, zso4_d, zso4_e, zso4_f
    real*8 e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3
    real*8 e3_zno30, dedz3_zno30, f3_zno30, dfdz3_zno30, g3_zno30, dgdz3_zno30
    real*8 e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3
    real*8 zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f
    real*8 e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox
    real*8 e2_zox0, dedz2_zox0, f2_zox0, dfdz2_zox0, g2_zox0, dgdz2_zox0
    real*8 e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox
    real*8 zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
    real*8 e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00
    real*8 e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0

    real*8 rH2S_A4, rH2S_B4
    real*8 rH2S_A3, rH2S_B3
    real*8 rH2S_A2, rH2S_B2
    real*8 rH2S_A1, rH2S_B1

    real*8 zso4FH2S, zoxFH2S
    real*8 flxswiH2S

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
    call prepfg_l12(0.0D00, 0.0D00, 0.0D00, 0.0D00, zox, DH2S1, DH2S2, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, ltype1)

    ! layer 2: zox < z < zno3, passive diffn
    call prepfg_l12(0.0D00, 0.0D00, 0.0D00, zox, zno3, DH2S1, DH2S2, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, ltype2)

    ! layer 3: zno3 < z < zso4, H2S consumption by OM oxidation
    call prepfg_l12(reac1_h2s, reac2_h2s, 0.0D00, zno3, zso4, DH2S1, DH2S2, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, ltype3)

    ! layer 4: zso4 < z < zinf, passive diffn
    call prepfg_l12(0.0D00, 0.0D00, 0.0D00, zso4, zinf, DH2S1, DH2S2, ls_a4, ls_b4, ls_c4, ls_d4, ls_e4, ls_f4, ltype4)

    ! Work up from the bottom, matching solutions at boundaries
    ! Basis functions at bottom of layer 4 zinf
    call calcfg_l12(zso4, 0.0D00, 0.0D00, 0.0D00, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, DH2S1, DH2S2, ltype3, &
                & e4_zinf, dedz4_zinf, f4_zinf, dfdz4_zinf, g4_zinf, dgdz4_zinf)

    ! Match at zso4, layer 3 - layer 4 (continuity and flux with AOM production)
    ! basis functions at bottom of layer 3
    call calcfg_l12(zso4, reac1_h2s, reac2_h2s, 0.0D00, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, DH2S1, DH2S2, ltype3, &
            & e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4)

    ! ... and top of layer 4
    call calcfg_l12(zso4, 0.0D00,  0.0D00, 0.0D00, ls_a4, ls_b4, ls_c4, ls_d4, ls_e4, ls_f4, DH2S1, DH2S2, ltype4, &
            & e4_zso4, dedz4_zso4, f4_zso4, dfdz4_zso4, g4_zso4, dgdz4_zso4)

    ! flux of H2S produced by AOM interface (Source of H2S)
    zso4FH2S = calcReac(zso4, zinf, SO4C, SO4C) ! MULTIPLY BY 1/POR ????
!    print*,'flux of H2S produced by AOM interface zso4FH2S = ', zso4FH2S

    ! match solutions at zso4 - continuous concentration and flux
    call matchsoln(e3_zso4, f3_zso4, g3_zso4, dedz3_zso4, dfdz3_zso4, dgdz3_zso4, &
                 & e4_zso4, f4_zso4, g4_zso4, dedz4_zso4, dfdz4_zso4, dgdz4_zso4, &
                 & 0.0D00, -zso4FH2S/DH2S2, zso4_a, zso4_b, zso4_c, zso4_d, zso4_e, zso4_f)

    ! Match at zno3, layer 2 - layer 3 (continuity and flux)
    ! basis functions at bottom of layer 2
    call calcfg_l12(zno3, 0.0D00, 0.0D00, 0.0D00, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DH2S1, DH2S2, ltype2, &
                    & e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3)

    ! ... and top of layer 3
    call calcfg_l12(zno3, reac1_h2s, reac2_h2s, 0.0D00, ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3, DH2S1, DH2S2, ltype3, &
                    & e3_zno30, dedz3_zno30, f3_zno30, dfdz3_zno30, g3_zno30, dgdz3_zno30)

    ! ... transformed to use coeffs from l4
    call xformsoln(e3_zno30, f3_zno30, g3_zno30, dedz3_zno30, dfdz3_zno30, dgdz3_zno30, &
                    zso4_a , zso4_b , zso4_c , zso4_d , zso4_e ,zso4_f, &
                    e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3)
    ! match solutions at zno3 - continuous concentration and flux
    call matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, &
                 & e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, &
                 & 0.0D00, 0.0D00, zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f)

    ! Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from H2S source)
    ! flux of H2S to oxic interface (from all sources of H2S below)
    ! NB: include methane region as AOM will produce sulphide as well..
    zoxFH2S = calcReac(zno3, zinf, SO4C, SO4C)  ! MULTIPLY BY 1/POR ????
!    print*,' '
!    print*,'flux of H2S to oxic interface zoxFH2S = ', zoxFH2S
!    print*,' '

    ! basis functions at bottom of layer 1
    call calcfg_l12(zox, 0.0D00, 0.0D00, 0.0D00, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DH2S1, DH2S2, ltype1, &
                    & e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox)

    ! basis functions at top of layer 2
    call calcfg_l12(zox, 0.0D00, 0.0D00, 0.0D00, ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2, DH2S1, DH2S2, ltype2, &
                    & e2_zox0, dedz2_zox0, f2_zox0, dfdz2_zox0, g2_zox0, dgdz2_zox0)

    !   transform to use coeffs from l4
    call xformsoln(e2_zox0, f2_zox0, g2_zox0, dedz2_zox0, dfdz2_zox0, dgdz2_zox0, &
                 & zno3_a , zno3_b , zno3_c , zno3_d , zno3_e ,zno3_f, &
                 & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox)

    ! match solutions at zox - continuous concentration, flux discontinuity from H2S ox

    IF(zox .le. zbio) THEN
        call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                     & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                     & 0.0D00, r_zxf*zoxFH2S/DH2S1, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
    ELSE
        call matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                     & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                     & 0.0D00, r_zxf*zoxFH2S/DH2S2, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
    END IF

    ! Solution at swi, top of layer 1
    call calcfg_l12(0.0D00, 0.0D00, 0.0D00, 0.0D00, ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1, DH2S1, DH2S2, ltype1, &
                  & e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00)

    ! transform to use coeffs from l4
    call xformsoln(e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00, &
                 & zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f, &
                 & e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0)

    ! Solve for AH2S, BH2S given boundary conditions (expressed in terms of transformed basis fns, layer 4 A, B)
    !  AH2S*dedz4_zinf   +  BH2S*dfz4_zinf  + dgz4_zinf = 0;          % zero flux at zinf
    !  AH2S*e1_0     +   BH2S*f1_0     + g1_0  = swi.H2S0;

    !  | dedz4_zinf dfdz4_zinf |  |AH2S|   = | -dgz4_zinf       |
    !  | e1_0     f1_0         |  |BH2S|     | swi.H2S0 - g1_0 |

    call solve2eqn(dedz4_zinf, dfdz4_zinf, e1_0, f1_0, -dgdz4_zinf, H2S0 - g1_0, rH2S_A4, rH2S_B4)

    ! flux at swi - DO include por so this is per cm^2 water column area
    flxswiH2S = por*DH2S1*(rH2S_A4*dedz1_0+rH2S_B4*dfdz1_0 + dgdz1_0)   ! NB: use A4, B4 as these are _xformed_ layer 1 basis functions

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



! *******************************************************************
!   *****************************************************************

!                       UTILITY FUNCTIONS

!   *****************************************************************
! *******************************************************************


    SUBROUTINE matchsoln(E_l, F_l, G_l, dEdx_l, dFdx_l, dGdx_l, &
                       & E_r, F_r, G_r, dEdx_r, dFdx_r, dGdx_r, &
                       & Vb, Db, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f)

        real*8, INTENT(in):: E_l, F_l, G_l, dEdx_l, dFdx_l, dGdx_l
        real*8, INTENT(in):: E_r, F_r, G_r, dEdx_r, dFdx_r, dGdx_r, Vb, Db
        real*8, INTENT(inout):: ls_a, ls_b, ls_c, ls_d, ls_e, ls_f
        real*8:: alden, blden

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


    SUBROUTINE xformsoln(E, F, G, dEdx, dFdx, dGdx, ls_a , ls_b , ls_c , ls_d , ls_e ,ls_f, Et, Ft, Gt, dEtdx, dFtdx, dGtdx)

            real*8, INTENT(in):: E, F, G, dEdx, dFdx, dGdx, ls_a , ls_b , ls_c , ls_d , ls_e ,ls_f
            real*8, INTENT(inout):: Et, Ft, Gt, dEtdx, dFtdx, dGtdx

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

    SUBROUTINE solve2eqn(a, b, c, d, e, f, x, y)

    ! Find soln of
    ! | a    b |  |x|   = | e |
    ! | c    d |  |y|     | f |

    real*8,INTENT(IN)::a, b, c, d, e, f
    real*8,INTENT(OUT)::x, y

    ! local variable
    real*8::det

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

    FUNCTION FUN_zO2(z)

    real*8 FUN_zO2, z, flxzox, conczox, flxswi, r_zxf

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

    real*8 FUN_zNO3, z, flxzno3, conczno3, flxswi, r_zxf

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

    real*8 FUN_zSO4, z, flxzso4, conczso4, flxswi, r_zxf

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

!!!! TODO: better put outside module, in kind of collection of auxiliary functions

    FUNCTION zbrent(func,x1,x2,tol)

    ! calculate root of func in the interval [x1,x2]

      INTEGER ITMAX
      REAL*8 zbrent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL*8 a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm

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
      11      continue
    END IF
! was      pause 'zbrent exceeding maximum iterations'
        print*,'zbrent exceeding maximum iterations'
        zbrent=b
        STOP

    return

    END FUNCTION zbrent


end module
