subroutine posterior(a,b,init,obs,m,p,n,miss,post)

    implicit none

    integer, intent(in) ::  m,p,n
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(n) :: obs,miss
    double precision, intent(inout), dimension(m,n) :: post
    double precision, dimension(m,n) :: beta,alpha
    double precision sumtmp,logscale,loglik
    double precision, dimension(m) :: tmp,tmp2
    integer i

    !forward

    if(miss(1)==0) then !check this!
        tmp = init*b(:,obs(1))
    else
        tmp = init
    end if
    alpha(:,1) = log(tmp)
    sumtmp = sum(tmp)
    logscale =log(sumtmp)
    tmp = tmp/sumtmp
    do i = 2,n
        call dgemv('t', m, m, 1.0d0, a, m, tmp, 1,0.0d0, tmp2,1)
        if(miss(i)==0) then !check this!
            tmp = tmp2*b(:,obs(i))
        else
            tmp = tmp2
        end if
        sumtmp = sum(tmp)
        logscale =logscale+log(sumtmp)
        tmp = tmp/sumtmp
        alpha(:,i) = log(tmp)+logscale
    end do

    !loglik
    loglik = logscale

    ! backward
    beta(:,n) = 0.0d0
    tmp = 1.0d0/dble(m)
    logscale =log(dble(m))
    do i = (n-1),1,-1

        if(miss(i+1)==0) then !check this!
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

   do i=1, n
    post(:,i) = alpha(:,i) + beta(:,i)
   end do
   post = post - loglik

end subroutine posterior

subroutine mvposterior(a,b,init,obs,m,p,n,miss,k,post)

    implicit none

    integer, intent(in) ::  m,p,n,k
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(k,n) :: obs,miss
    double precision, intent(inout), dimension(m,n,k) :: post
    integer i
    do i=1,k
        call posterior(a,b,init,obs(i,:),m,p,n,miss(i,:),post(:,:,i))
    end do
end subroutine mvposterior
