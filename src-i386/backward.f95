subroutine backward(a,b,init,obs,m,p,n,miss,beta)

    implicit none

    integer, intent(in) ::  m,p,n
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(n) :: obs,miss
    double precision, intent(inout), dimension(m,n) :: beta
    double precision sumtmp,logscale
    double precision, dimension(m) :: tmp,tmp2
    integer i

    beta(:,n) = 0.0d0
    tmp = 1.0d0/dble(m)
    logscale =log(dble(m))
    do i = n-1,1,-1

        if(miss(i)==0) then !check this!
            tmp2 = tmp*b(:,obs(i+1))
        else
            tmp2 = tmp
        end if
        call dgemv('n', m, m, 1.0d0, a, m, tmp2, 1,0.0d0, tmp,1)
        beta(:,i) = log(tmp)+logscale
        sumtmp = sum(tmp)
        tmp = tmp/sumtmp
        logscale =logscale+log(sumtmp)

    end do


end subroutine backward


subroutine mvbackward(a,b,init,obs,m,p,n,miss,k,beta)

    implicit none

    integer, intent(in) ::  m,p,n,k
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(k,n) :: obs,miss
    double precision, intent(inout), dimension(m,n,k) :: beta
    integer i

    do i = 1,k
        call backward(a,b,init,obs(i,:),m,p,n,miss(i,:),beta(:,:,i))
    end do
end subroutine mvbackward
