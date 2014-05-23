subroutine viterbi(a,b,init,obs,m,p,n,miss,q,logp)

    implicit none

    integer, intent(in) ::  m,p,n
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(n) :: obs,miss
    double precision, dimension(m,n) :: delta
    integer, dimension(m,n) :: phi
    integer,  intent(inout),dimension(n) :: q
    double precision, intent(inout) :: logp
    integer t,j

    if(miss(1)==0) then !check this!
        delta(:,1) = init+b(:,obs(1))
    else
        delta(:,1) = init
    end if

    phi(:,1) = 0
    do t=2,n
        do j=1,m
            phi(j,t) = maxloc(delta(:,t-1)+a(:,j),dim=1)
            if(miss(t)==0) then !check this!
                delta(j,t) = delta(phi(j,t),t-1)+a(phi(j,t),j)+b(j,obs(t))
                !delta(j,t) = maxval(delta(:,t-1)+a(:,j))+b(j,obs(t))
            else
                delta(j,t) = delta(phi(j,t),t-1)+a(phi(j,t),j)
            end if
        end do
    end do
    q(n) = maxloc(delta(:,n),dim=1)

    do t = (n-1), 1,-1
        q(t) = phi(q(t+1),t+1)
    end do
    logp = maxval(delta(:,n))

end subroutine viterbi

subroutine mvviterbi(a,b,init,obs,m,p,n,miss,k,q,logp)

    implicit none

    integer, intent(in) ::  m,p,n,k
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(k,n) :: obs,miss
    integer, intent(inout), dimension(k,n) :: q
    double precision, dimension(k), intent(inout) :: logp
    integer i
    do i=1,k
        call viterbi(a,b,init,obs(i,:),m,p,n,miss(i,:),q(i,:),logp(i))
    end do
end subroutine mvviterbi
