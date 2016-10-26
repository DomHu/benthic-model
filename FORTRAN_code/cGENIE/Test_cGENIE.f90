    SUBROUTINE sub_matchsoln(E_l, F_l, G_l, dEdx_l, dFdx_l, dGdx_l, &
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

    END SUBROUTINE sub_matchsoln


    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------


    SUBROUTINE sub_xformsoln(E, F, G, dEdx, dFdx, dGdx, ls_a , ls_b , ls_c , ls_d , ls_e ,ls_f, Et, Ft, Gt, dEtdx, dFtdx, dGtdx)

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

    !            print*, 'IN sub_xformsoln E, F, G, dEdx, dFdx, dGdx, ls_a , ls_b , ls_c , ls_d , ls_e ,ls_f', &
    !                    & E, F, G, dEdx, dFdx, dGdx, ls_a , ls_b , ls_c , ls_d , ls_e ,ls_f
    !            print*, 'OUT sub_xformsoln Et, Ft, Gt, dEtdx, dFtdx, dGtdx', Et, Ft, Gt, dEtdx, dFtdx, dGtdx

    END SUBROUTINE sub_xformsoln

    ! ****************************************************************************************************************************** !
    !   *****************************************************************
    !   *****************************************************************
    !------------------------------------------------------------------------------------

    SUBROUTINE sub_solve2eqn(a, b, c, d, e, f, x, y)

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

    END SUBROUTINE sub_solve2eqn


