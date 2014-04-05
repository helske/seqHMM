subroutine mchmmloglik(a,b,init,obs,m,p,n,miss,loglik,r)

    implicit none

    integer, intent(in) ::  m,p,n,r
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p,r) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(n,r) :: obs,miss
    double precision, intent(inout) :: loglik
    double precision sumtmp,logscale
    double precision, dimension(m) :: tmp,tmp2
    integer i,j

    tmp = init
        do j=1, r
        if(miss(1,j)==0) then !check this!
            tmp=tmp*b(:,obs(1,j),j)
        end if
    end do

    sumtmp = sum(tmp)
    loglik =log(sumtmp)
    tmp = tmp/sumtmp
    do i = 2,n
        call dgemv('t', m, m, 1.0d0, a, m, tmp, 1,0.0d0, tmp2,1)
        tmp = tmp2
        do j=1,r
            if(miss(i,j)==0) then !check this!
                tmp = tmp*b(:,obs(i,j),j)

            end if
        end do
        sumtmp = sum(tmp)
        loglik =loglik+log(sumtmp)
        tmp = tmp/sumtmp
    end do

end subroutine mchmmloglik

subroutine mvmchmmloglik(a,b,init,obs,m,p,n,miss,k,loglik,r)

    implicit none

    integer, intent(in) ::  m,p,n,k,r
    integer i
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p,r) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(k,n,r) :: obs,miss
    double precision, intent(inout) :: loglik
    double precision tmp

    loglik = 0.0d0
    do i = 1,k
        call mchmmloglik(a,b,init,obs(i,:,:),m,p,n,miss(i,:,:),tmp,r)
        loglik=loglik+tmp
    end do

end subroutine mvmchmmloglik
