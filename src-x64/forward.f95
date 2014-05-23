subroutine forward(a,b,init,obs,m,p,n,miss,alpha)

    implicit none

    integer, intent(in) ::  m,p,n
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(n) :: obs,miss
    double precision, intent(inout), dimension(m,n) :: alpha
    double precision sumtmp,logscale
    double precision, dimension(m) :: tmp,tmp2
    integer i

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
end subroutine forward

subroutine mvforward(a,b,init,obs,m,p,n,miss,k,alpha)

    implicit none

    integer, intent(in) ::  m,p,n,k
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(k,n) :: obs,miss
    double precision, intent(inout), dimension(m,n,k) :: alpha
    integer i

    do i = 1,k
        call forward(a,b,init,obs(i,:),m,p,n,miss(i,:),alpha(:,:,i))
    end do
end subroutine mvforward
