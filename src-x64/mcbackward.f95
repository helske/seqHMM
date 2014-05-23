subroutine mcbackward(a,b,obs,m,p,n,miss,beta,r)

    implicit none

    integer, intent(in) ::  m,p,n,r
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p,r) :: b
    integer, intent(in), dimension(n,r) :: obs,miss
    double precision, intent(inout), dimension(m,n) :: beta
    double precision sumtmp,logscale
    double precision, dimension(m) :: tmp,tmp2
    integer i,j

    beta(:,n) = 0.0d0
    tmp = 1.0d0/dble(m)
    logscale =log(dble(m))
    do i = (n-1),1,-1
        tmp2 = tmp
        do j=1,r
            if(miss(i+1,j)==0) then !check this!
                tmp2 = tmp2*b(:,obs(i+1,j),j)
            end if
        end do
        call dgemv('n', m, m, 1.0d0, a, m, tmp2, 1,0.0d0, tmp,1)
        beta(:,i) = log(tmp)+logscale
        sumtmp = sum(tmp)
        tmp = tmp/sumtmp
        logscale =logscale+log(sumtmp)

    end do


end subroutine mcbackward


subroutine mvmcbackward(a,b,obs,m,p,n,miss,k,beta,r)

    implicit none

    integer, intent(in) ::  m,p,n,k,r
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p,r) :: b
    integer, intent(in), dimension(k,n,r) :: obs,miss
    double precision, intent(inout), dimension(m,n,k) :: beta
    integer i

    do i = 1,k
        call mcbackward(a,b,obs(i,:,:),m,p,n,miss(i,:,:),beta(:,:,i),r)
    end do
end subroutine mvmcbackward
