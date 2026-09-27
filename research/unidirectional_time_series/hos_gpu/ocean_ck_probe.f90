! SPDX-License-Identifier: GPL-3.0-or-later
! Fixed attempted steps of the original HOS-Ocean Cash-Karp routine.
SUBROUTINE export_ck_reference(case_id)
USE Runge_Kutta, ONLY: RK_parameters,fill_butcher_array,RK_adapt_2var_3D_in_mo_lin
IMPLICIT NONE
INTEGER,INTENT(IN) :: case_id
INTEGER :: trial,unit
REAL(RP) :: start,step,estimate
TYPE(RK_parameters) :: method
COMPLEX(CP),ALLOCATABLE :: e(:,:),p(:,:),derivative(:,:)
CHARACTER(LEN=80) :: filename
ALLOCATE(e(m1o2p1,m2),p(m1o2p1,m2),derivative(m1o2p1,m2))
method%A=0.0_rp;method%b=0.0_rp;method%c=0.0_rp;method%e=0.0_rp
CALL fill_butcher_array(method)
i_rlx=0
DO trial=1,3
    start=0.0_rp;step=.2_rp
    IF(trial==1)step=.01_rp
    IF(trial==3)start=.37_rp
    e=a_eta;p=a_phis
    CALL RK_adapt_2var_3D_in_mo_lin(0,method,start,step,p,e,derivative,estimate,1.0_rp,1.0_rp)
    WRITE(filename,'(A,I1,A,I1,A)') 'ck_case',case_id,'_trial',trial,'.bin'
    OPEN(newunit=unit,file=TRIM(filename),access='stream',form='unformatted',status='new')
    WRITE(unit)INT(n1,4),INT(n2,4),INT(M,4),INT(case_id,4),INT(trial,4)
    WRITE(unit)xlen_star,ylen_star,depth_star,g_star,start,step,estimate
    WRITE(unit)a_eta(1:n1o2p1,1:n2),a_phis(1:n1o2p1,1:n2),e(1:n1o2p1,1:n2),p(1:n1o2p1,1:n2)
    CLOSE(unit)
    WRITE(*,*) 'CASH_KARP_EXPORTED ',TRIM(filename),estimate
ENDDO
END SUBROUTINE export_ck_reference
