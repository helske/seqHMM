subroutine mcviterbi(a,b,init,obs,m,p,n,miss,q,r,logp)

    implicit none

    integer, intent(in) ::  m,p,n,r
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p,r) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(n,r) :: obs,miss
    double precision, dimension(m,n) :: delta
    integer, dimension(m,n) :: phi
    integer,  intent(inout),dimension(n) :: q
    double precision, intent(inout) :: logp
    integer t,j,i

    delta(:,1) = init
    do i= 1,r
        if(miss(1,i)==0) then
            delta(:,1) =  delta(:,1)+b(:,obs(1,i),i)
        end if
    end do

    phi(:,1) = 0
    do t=2,n
        do j=1,m
            phi(j,t) = maxloc(delta(:,t-1)+a(:,j),dim=1)
            delta(j,t) = delta(phi(j,t),t-1)+a(phi(j,t),j)

            do i= 1,r
                if(miss(t,i)==0) then !check this!
                    delta(j,t) = delta(j,t)+b(j,obs(t,i),i)
                end if
            end do
        end do
    end do
    q(n) = maxloc(delta(:,n),dim=1)

    do t = (n-1), 1,-1
        q(t) = phi(q(t+1),t+1)
    end do
    logp = maxval(delta(:,n))

end subroutine mcviterbi

subroutine mvmcviterbi(a,b,init,obs,m,p,n,miss,k,q,r,logp)

    implicit none

    integer, intent(in) ::  m,p,n,k,r
    double precision, intent(in), dimension(m,m) :: a
    double precision, intent(in), dimension(m,p,r) :: b
    double precision, intent(in), dimension(m) :: init
    integer, intent(in), dimension(k,n,r) :: obs,miss
    integer, intent(inout), dimension(k,n) :: q
    double precision, dimension(k), intent(inout) :: logp
    integer i
    do i=1,k
        call mcviterbi(a,b,init,obs(i,:,:),m,p,n,miss(i,:,:),q(i,:),r,logp(i))
    end do
end subroutine mvmcviterbi
