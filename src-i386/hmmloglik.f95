subroutine hmmloglik(a,b,init,obs,m,p,n,miss,loglik)

    implicit none

    integer, intent(in) ::  m,p,n
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(n) :: obs,miss
    double precision, intent(inout) :: loglik
    double precision sumtmp
    double precision, dimension(m) :: tmp,tmp2
    integer i

    if(miss(1)==0) then !check this!
        tmp = init*b(:,obs(1))
    else
        tmp = init
    end if


    sumtmp = sum(tmp)
    loglik =log(sumtmp)
    tmp = tmp/sumtmp
    do i = 2,n
        call dgemv('t', m, m, 1.0d0, a, m, tmp, 1,0.0d0, tmp2,1)
        if(miss(i)==0) then !check this!
            tmp = tmp2*b(:,obs(i))
        else
            tmp = tmp2
        end if
        sumtmp = sum(tmp)
        loglik =loglik+log(sumtmp)
        tmp = tmp/sumtmp
    end do
!    implicit none
!
!    integer, intent(in) ::  m,p,n
!    integer, intent(inout) ::  error
!    double precision, intent(in), dimension(m,m) :: a
!    double precision, intent(in), dimension(m,p) :: b
!    double precision, intent(in), dimension(m) :: pi
!    integer, intent(in), dimension(n) :: obs,miss
!    double precision, intent(inout) :: loglik
!    double precision, intent(in) :: tol
!    double precision, dimension(m) :: alpha,tmp
!    double precision s
!    integer i
!
!    loglik = 0.0d0
!    ! Unscaled forward probability:
!    alpha = pi*b(:,obs(1))
!    s = sum(alpha)
!    if(s< tol) then
!        error = 1
!    else
!        loglik = loglik + log(s)
!        ! Scaled forward probability:
!        alpha = alpha/s
!
!        ! Recursion for scaled forward probabilities
!        do i = 1, n-1
!            if(miss(i+1)==0) then !check this!
!                tmp = alpha*b(:,obs(i+1))
!            else
!                tmp = alpha
!            end if
!            call dgemv('t', m, m, 1.0d0, a, m, tmp, 1,0.0d0, alpha,1)
!            s = sum(alpha)
!            if(s< tol) then
!                error = i+1
!            else
!                if(miss(i+1)==0) then !check this!
!                    loglik = loglik + log(s) !Rabiner eq 91 and 103
!                end if
!                ! Scaled forward probability:
!                alpha = alpha/s
!            end if
!        end do
!    end if

end subroutine hmmloglik

subroutine mvhmmloglik(a,b,init,obs,m,p,n,miss,k,loglik)

    implicit none

    integer, intent(in) ::  m,p,n,k
    integer i
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(k,n) :: obs,miss
    double precision, intent(inout) :: loglik
    double precision tmp

    loglik = 0.0d0

    do i = 1,k
        call hmmloglik(a,b,init,obs(i,:),m,p,n,miss(i,:),tmp)
        loglik=loglik+tmp
    end do

end subroutine mvhmmloglik
