module subprog  ! runge-kutta
implicit none
contains
 subroutine fn(E,y,l,j,h,V,W,isospin)
real(8)::d(4,2),x(201),jl,T ! y=y(,1) dy/dx=y(,2)
integer::n,i,k
integer,intent(in)::isospin
real(8),intent(in)::E,l,j,h,V(201,2),W(201,2)
real(8),intent(out):: y(201,2)
T = (939.d0*2.d0/197.d0**2.d0)
do i=1,201
   do n=1,2
     y(i,n)=0.d0
  end do
end do
do i=1,4
   do n=1,2
      d(i,n)=0.d0
   end do
enddo
jl = (j*(j+1.d0)-l*(l+1.d0)-3.d0/4.d0)
y(1,1)= 1.d0
y(1,2)=y(1,1)/h
k = isospin
do i=1,200
   x(i)=h*dble(i)

   d(1,2)=h*T*(l*(l+1.d0)/(x(i)**2.d0*T)-E+V(i,k)&
        +jl*0.5d0*W(i,k))*y(i,1)

   d(1,1)=h*y(i,2)

   d(2,2)=h*T*(l*(l+1.d0)/((x(i)+h*0.5d0)**2.d0*T)-E+V(i,k)&
        +jl*0.5d0*W(i,k))*(y(i,1)+d(1,1)/2.d0)

   d(2,1)=h*(y(i,2)+d(1,2)/2.d0)

   d(3,2)=h*T*(l*(l+1.d0)/((x(i)+h*0.5d0)**2.d0*T)-E+V(i,k)&
        +jl*0.5d0*W(i,k))*(y(i,1)+d(2,1)/2.d0)

   d(3,1)=h*(y(i,2)+d(2,2)/2.d0)

   d(4,2)=h*T*(l*(l+1.d0)/((x(i)+h)**2.d0*T)-E+V(i,k)&
        +jl*0.5d0*W(i,k))*(y(i,1)+d(3,1))

   d(4,1)=h*(y(i,2)+d(3,2))

   do n = 1,2
      y(i+1,n)=y(i,n)+(d(1,n)+2.0d0*d(2,n)+2.0d0*d(3,n)+d(4,n))/6.d0
   end do
end do
end subroutine fn
end module subprog

module node     !count node
implicit none
contains
subroutine gn (y,P)
  real(8),intent(in)::y(201,2)
  integer,intent(out)::P
  integer::i
  P =0
   do i = 2,200
      if (y(i,1)*y(i+1,1)<0.d0) then
         P = P+1
       endif
   enddo
end subroutine gn
end module node

module normalize ! normalize wavefunction
implicit none
contains
  subroutine hn(y,S,A,h)
real(8),intent(in)::h
real(8),intent(out)::S,A
real(8),intent(inout)::y(201,2)
integer::i
S  = 1.d0
A  = 1.d0
do i = 2,200
   S = S + h*y(i,1)**2.d0
end do
   S = S + h*y(1,1)**2.d0 +h*y(201,1)**2.d0
A = S**(0.5d0)
do i = 1,201
   y(i,1)=y(i,1)/A
end do
end subroutine hn
end module normalize

program main
  use subprog
  use node
  use normalize
  implicit none
 real(8) :: EL, EM, ESm, h,S,A,l,R,r0,B,a0,Z,e,j,NA,ea,pi
 real(8) :: y(201,2),f(201),df(201),x(201),V(201,2),W(201,2)
 real(8) ::SM(201),es(42),sa
 integer :: i,k,ni,PL,PM,PS,li,n,ji,isospin,NAi,Zi,Ai,m,sf,state(3,42),testni,testli,testji,statenum,iso
 EM = 0.d0; PL = 0; PM = 0; PS = 0; li = 0; ji = 0; j  = 0.d0 !substitute IV
 e  = 197.d0/137.d0 !実際はe**2
 r0 = 127.d-2
 a0 = 67.d-2
 pi = atan(1.d0)*4.d0
 testni=0;testli=0;testji=0
 !write(*,*)"nutron",":proton",":n0~3",",:l0~6",",:Ji1,2"
 read(*,*)NAi,Zi !中性子数 陽子数の順に読み込む。
 read(*,*)iso !testni,testli,testji,iso
  NA =dble(NAi)
  Z = dble(Zi)
  A =NA+Z
  R = r0*A**(1.d0/3.d0)
  h  = 0.1d0
  do i = 1,201
    x(i)=h*dble(i)
    f(i)=1.d0/(1+exp((x(i)-R)/a0))
    df(i)= -exp((x(i)-R)/a0)/(1+exp((x(i)-R)/a0))**2.d0/a0
    do n =1,2
        y(i,n)=0.d0
     end do
  end do
  do i = 1,201
     SM(i) = 0.d0
if (x(i) < R) then
  V(i,1) = -(51.d0+33.d0*(NA-Z)/A)*f(i) + (Z-1.d0)*e*(3.d0-(x(i)/R)**2)/(2.d0*R)
  else
   V(i,1) = -(51.d0+33.d0*(NA-Z)/A)*f(i) + (Z-1.d0)*e/x(i)+(l*(l+1.d0))
  endif
     V(i,2) = -(51.d0-33.d0*(NA-Z)/A)*f(i)
     W(i,1) =  (22.d0+14.d0*(NA-Z)/A)*r0**2.d0/x(i)*df(i)
     W(i,2) =  (22.d0-14.d0*(NA-Z)/A)*r0**2.d0/x(i)*df(i)
  end do

  sf = 0     !end subatitute Initial value
do isospin  = 1,2  !1=proton,2=nutron
 do li = 0,6
   l = dble(li)
   do n  = 1,3
      do ji = 1,2
        if (ji == 1) then
           j = l + 0.5d0
        else
           j = l - 0.5d0
        end if
           EL =  100.d0
           ESm = -50.d0
          do i=1,100   !bisection method
            EM = (EL + ESm)/2.d0
            call fn(EL,y,l,j,h,V,W,isospin)
            call gn(y,PL)
            call fn(EM,y,l,j,h,V,W,isospin)
            call gn(y,PM)
           if((dble(PL) -dble(n) +0.5d0)*(dble(PM) -dble(n) +0.5d0) > 0.d0) then
              EL = EM
            else
              ESm = EM
            end if
         end do     !end bisectoin method
   if( (n ==testni .and. li == testli).and.ji == testji.and.isospin ==iso ) then
         call fn(em,y,l,j,h,V,W,2)
         call hn(y,S,A,h)
         do i = 1,201
            !  write(*,*) x(i), y(i,1)
          enddo
      endif
          if (isospin == iso) then
             sf = sf+1
             es(sf) = EM
             state(1,sf) = n
             state(2,sf) = li
             state(3,sf) = ji
          endif
       enddo
    enddo
 enddo
end do

  isospin = 2
  do sf = 1,41   !sort energy
     do k = sf+1,42
        if (es(sf) > es(k)) then
           ea = es(sf)
           es(sf) = es(k)
           es(k)= ea
           do i = 1,3
              sa = state(i,k)
              state(i,k) = state(i,sf)
              state(i,sf)= sa
           enddo
        endif
     enddo
  enddo  !end sort
  statenum = 1
  m = 0
  do while (statenum < Nai+1 )
     m = m + 1
     em = es (m)
     l =dble(state(2,m)) !li
     if  (state(3,m)==1) then !ji
        j = l+0.5d0
     else
        j = l-0.5d0
     endif
     if (j >= 0.d0) then
        call fn(em,y,l,j,h,V,W,iso)
        call hn(y,S,A,h)
        do k = 1,nint(2.d0*j+1.d0)
           ! write(*,*)"E=",em,"n=",state(1,m),"L=",nint(l),"J=",state(3,m),statenum
           if( (state(1,m) ==testni .and. state(2,m) ==
testli).and.(state(3,m) == testji) ) then
              do i = 1,201
                 !  write(*,*) x(i),y(i,1)
              enddo
           endif
           do i = 1,201
              SM(i) = SM(i)+y(i,1)**2/(x(i)**2*pi*4.d0)
           enddo
           statenum = statenum+1
        enddo
     endif
  enddo
  !do i = 1,42
   ! write(*,*) es(i),state(1,i),state(2,i),state(3,i)
  !end do

 do i = 1,201
     write(*,*)x(i),SM(i)
 enddo
end program main