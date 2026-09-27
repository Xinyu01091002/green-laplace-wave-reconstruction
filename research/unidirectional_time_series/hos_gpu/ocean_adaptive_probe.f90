! SPDX-License-Identifier: GPL-3.0-or-later
! Main controller is inserted verbatim from the production main, with logging only.
SUBROUTINE export_adaptive_reference(case_id)
IMPLICIT NONE
INTEGER,INTENT(IN)::case_id
INTEGER::out_unit,trace_unit,attempts,trace_accept,noutputs
REAL(RP)::trace_before,trace_proposal,trace_actual
LOGICAL::is_actual
CHARACTER(LEN=80)::filename
INQUIRE(file='actual_initial.flag',exist=is_actual)
IF(.NOT.is_actual) THEN
    T_stop_star=.57_rp;dt_out=.13_rp;dt=.2_rp;dt_lin=.5_rp
ELSE
    T_stop_star=.3_rp/T
    ! Keep the original dt_out, dt and dt_lin computed by the production main.
ENDIF
toler=1.0e-12_rp;phis_scale=1.0_rp;eta_scale=1.0_rp
RK_param%A=0.0_rp;RK_param%b=0.0_rp;RK_param%c=0.0_rp;RK_param%e=0.0_rp
CALL fill_butcher_array(RK_param)
time_cur=0.0_rp;attempts=0;ibrk=0;i_rlx=0
! RHS export has already evaluated initial gradients for the original slope check.
noutputs=CEILING(T_stop_star/dt_out)
WRITE(filename,'(A,I1,A)')'adaptive_case',case_id,'.bin'
OPEN(newunit=out_unit,file=TRIM(filename),access='stream',form='unformatted',status='new')
WRITE(out_unit)INT(n1,4),INT(n2,4),INT(M,4),INT(case_id,4),INT(noutputs,4)
WRITE(out_unit)xlen_star,ylen_star,depth_star,g_star,T_stop_star,dt_out,dt,dt_lin,toler,L,T
WRITE(out_unit)a_eta(1:n1o2p1,1:n2),a_phis(1:n1o2p1,1:n2)
WRITE(filename,'(A,I1,A)')'adaptive_case',case_id,'.trace'
OPEN(newunit=trace_unit,file=TRIM(filename),status='new')
WRITE(trace_unit,'(A)')'attempt t_before h_proposed h_actual error t_after accepted h_next'
DO WHILE(time_cur<T_stop_star)
!__ORIGINAL_CONTROLLER__
    WRITE(out_unit)time_cur,a_eta(1:n1o2p1,1:n2),a_phis(1:n1o2p1,1:n2)
ENDDO
CLOSE(out_unit);CLOSE(trace_unit)
WRITE(*,*)'ADAPTIVE_EXPORTED ',case_id,time_cur,attempts
END SUBROUTINE export_adaptive_reference
