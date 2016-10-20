program main
   use benthic

    implicit none

    call initialize()
    call zTOC()
    call zO2()
    call benthic_zNO3()
    call benthic_zSO4()
    call benthic_zNH4()
    call benthic_zH2S()
    call benthic_zPO4_M()
end program main
