! SPDX-License-Identifier: GPL-3.0-or-later
! Read-only test harness appended to a COPY of the actual HOS-Ocean MPI main.
! Existing numerical module objects are linked unchanged. No time integration.
SUBROUTINE export_gpu_reference()
USE resol_HOS, ONLY: solveHOS_lin
IMPLICIT NONE
INTEGER :: test_id, ix, iy, mx, my, unit, count
INTEGER :: first_case,last_case
LOGICAL :: actual_initial,cash_karp_check,adaptive_check
REAL(RP) :: angle, amplitude_test, kk
REAL(RP), ALLOCATABLE :: ee(:,:), pp(:,:), en(:,:), pn(:,:), ef(:,:), pf(:,:)
COMPLEX(CP), ALLOCATABLE :: de(:,:), dp(:,:)
CHARACTER(LEN=64) :: filename
IF (nb_procs /= 1) STOP 'This reference exporter requires MPI with one rank'
ALLOCATE(ee(Nd1,Nd2),pp(Nd1,Nd2),en(n1,n2),pn(n1,n2),ef(n1,n2),pf(n1,n2))
ALLOCATE(de(m1o2p1,m2),dp(m1o2p1,m2))
! Explicit flat-bed, no-current, no-breaking, no-absorption test scope.
cuxst=0.0_rp;dUdX=0.0_rp;v_eddy=0.0_rp;filt_eta=0.0_rp;filt_phis=0.0_rp
diff_term_eta=0.0_rp;diff_term_phis=0.0_rp;nu=0.0_rp;gamma_d=0.0_rp
i_abs=0;Ta=0.0_rp;zbeta=.FALSE.
INQUIRE(file='actual_initial.flag',exist=actual_initial)
INQUIRE(file='cash_karp.flag',exist=cash_karp_check)
INQUIRE(file='adaptive.flag',exist=adaptive_check)
first_case=1;last_case=4
IF(actual_initial) THEN
    first_case=5;last_case=5
ENDIF
DO test_id=first_case,last_case
  IF (.NOT.actual_initial) THEN
    a_eta=0.0_cp;a_phis=0.0_cp
    a_eta(1,1)=0.001_rp;a_phis(1,1)=0.003_rp
    IF (test_id==1 .OR. test_id==4) THEN
        a_eta(2,3)=CMPLX(.006_rp,.001_rp,CP);a_phis(2,3)=CMPLX(-.002_rp,.004_rp,CP)
        a_eta(4,n2)=CMPLX(.002_rp,-.001_rp,CP);a_phis(4,n2)=CMPLX(.003_rp,.001_rp,CP)
    ELSEIF (test_id==2) THEN
        a_eta(n1/2-1,n2/2-1)=CMPLX(.0002_rp,.0001_rp,CP)
        a_phis(n1/2-1,n2/2-1)=CMPLX(-.0001_rp,.0002_rp,CP)
        a_eta(n1/4,n2-n2/4+2)=CMPLX(.0003_rp,-.0001_rp,CP)
        a_phis(n1/4,n2-n2/4+2)=CMPLX(.0001_rp,.0002_rp,CP)
    ELSEIF (test_id==3) THEN
        count=(n1/2-1)*(n2-1)
        amplitude_test=.002_rp/SQRT(REAL(count,RP))
        DO mx=1,n1/2-1
            DO my=-n2/2+1,n2/2-1
                ix=mx+1;iy=MODULO(my,n2)+1
                angle=SIN(17.13_rp*mx+9.71_rp*my)*100.0_rp
                a_eta(ix,iy)=amplitude_test*CMPLX(COS(angle),SIN(angle),CP)
                a_phis(ix,iy)=amplitude_test*CMPLX(SIN(angle),-COS(angle),CP)
            ENDDO
        ENDDO
    ENDIF
    IF(test_id==4) THEN
        a_eta(n1/2+1,1)=.0002_rp;a_phis(n1/2+1,1)=.0001_rp
        a_eta(2,n2/2+1)=CMPLX(.0001_rp,.0002_rp,CP)
        a_phis(2,n2/2+1)=CMPLX(.0002_rp,-.0001_rp,CP)
        a_eta(1,n2/2+1)=.0001_rp;a_phis(1,n2/2+1)=.0002_rp
    ENDIF
  ENDIF
    CALL fourier_2_space(a_eta,eta);CALL fourier_2_space(a_phis,phis)
    temp_C_Nd=extend_C(a_eta);CALL fourier_2_space_big(temp_C_Nd,ee)
    temp_C_Nd=extend_C(a_phis);CALL fourier_2_space_big(temp_C_Nd,pp)
    CALL solveHOS_lin(a_phis,dp,a_eta,de,0.0_rp)
    ! Remove its separately included mean gravity term to define pure N(eta,psi).
    dp(1,1)=dp(1,1)+g_star*a_eta(1,1)
    CALL fourier_2_space(de,en);CALL fourier_2_space(dp,pn)
    DO iy=1,n2
        my=iy-1;IF(my>n2/2)my=my-n2
        DO ix=1,n1o2p1
            kk=SQRT((TWOPI*(ix-1)/xlen_star)**2+(TWOPI*my/ylen_star)**2)
            de(ix,iy)=de(ix,iy)+kk*TANH(kk*depth_star)*a_phis(ix,iy)
            dp(ix,iy)=dp(ix,iy)-g_star*a_eta(ix,iy)
        ENDDO
    ENDDO
    CALL fourier_2_space(de,ef);CALL fourier_2_space(dp,pf)
    WRITE(filename,'(A,I1,A)') 'ocean_case',test_id,'.bin'
    OPEN(newunit=unit,file=TRIM(filename),access='stream',form='unformatted',status='new')
    WRITE(unit) INT(n1,4),INT(n2,4),INT(Nd1,4),INT(Nd2,4),INT(M,4),INT(test_id,4)
    WRITE(unit) xlen_star,ylen_star,depth_star,g_star
    WRITE(unit) ee,pp,eta(1:n1,1:n2),phis(1:n1,1:n2),en,pn,ef,pf
    CLOSE(unit)
    WRITE(*,*) 'EXPORTED ',TRIM(filename),n1,n2,Nd1,Nd2,M,RP,CP
    IF(cash_karp_check)CALL export_ck_reference(test_id)
    IF(adaptive_check .AND. test_id/=3)CALL export_adaptive_reference(test_id)
ENDDO
END SUBROUTINE export_gpu_reference
