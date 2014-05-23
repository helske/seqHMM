subroutine mcposterior(a,b,init,obs,m,p,n,miss,post,r)

    implicit none

    integer, intent(in) ::  m,p,n,r
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p,r) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(n,r) :: obs,miss
    double precision, intent(inout), dimension(m,n) :: post
    double precision, dimension(m,n) :: beta,alpha
    double precision sumtmp,logscale,loglik
    double precision, dimension(m) :: tmp,tmp2
    integer i,j

    !forward

    tmp = init
    do j=1, r
        if(miss(1,j)==0) then !check this!
            tmp=tmp*b(:,obs(1,j),j)
        end if
    end do
    alpha(:,1) = log(tmp)
    sumtmp = sum(tmp)
    logscale =log(sumtmp)
    tmp = tmp/sumtmp
    do i = 2,n
        call dgemv('t', m, m, 1.0d0, a, m, tmp, 1,0.0d0, tmp2,1)
        tmp = tmp2
        do j = 1,r
            if(miss(i,j)==0) then !check this!
                tmp = tmp*b(:,obs(i,j),j)
            end if
        end do
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
        tmp2 = tmp
        do j = 1,r
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

    do i=1, n
        post(:,i) = alpha(:,i) + beta(:,i)
    end do
    post = post - loglik

end subroutine mcposterior

subroutine mvmcposterior(a,b,init,obs,m,p,n,miss,k,post,r)

    implicit none

    integer, intent(in) ::  m,p,n,k,r
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p,r) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(k,n,r) :: obs,miss
    double precision, intent(inout), dimension(m,n,k) :: post
    integer i
    do i=1,k
        call mcposterior(a,b,init,obs(i,:,:),m,p,n,miss(i,:,:),post(:,:,i),r)
    end do
end subroutine mvmcposterior
