module subprog  ! a=2000*hの範囲内の波関数をmainに返す。
implicit none 
contains 
 subroutine fn(E,y,l,j,h,V,W,isospin,T) 
real(8)::d(4,2),x(2001) ! y=y(,1) dy/dx=y(,2)
integer::n,i,k
integer,intent(in)::isospin
real(8),intent(in)::E,l,j,h,T,V(2001,2),W(2001,2)
real(8),intent(out):: y(2001,2)
do i=1,2001
   do n=1,2
     y(i,n)=0.d0 
  end do
end do
do i=1,4
   do n=1,2
      d(i,n)=0.d0
   end do
enddo
y(1,1)=0.d0
y(1,2)=1.d0
k = isospin
do i=1,2000                    !l=1 j = 1.5 で特異
   x(i)=h*dble(i) 
   d(1,2)=h*T*(l*(l+1.d0)/(x(i)**2.d0*T)-E+V(i,k)&
        +(j*(j+1.d0)-l*(l+1.d0)-3.d0/4.d0)*0.5d0*W(i,k))*y(i,1)
   d(1,1)=h*y(i,2)
   d(2,2)=h*T*(l*(l+1.d0)/((x(i)+h*0.5d0)**2.d0*T)-E+V(i,k)&
        +(j*(j+1.d0)-l*(l+1.d0)-3.d0/4.d0)*0.5d0*W(i,k))*(y(i,1)+d(1,1)/2.d0)
   d(2,1)=h*(y(i,2)+d(1,2)/2.d0)
   d(3,2)=h*T*(l*(l+1.d0)/((x(i)+h*05d0)**2.d0*T)-E+V(i,k)&
        +(j*(j+1.d0)-l*(l+1.d0)-3.d0/4.d0)*0.5d0*W(i,k))*(y(i,1)+d(2,1)/2.d0)
   d(3,1)=h*(y(i,2)+d(2,2)/2.d0)
   d(4,2)=h*T*(l*(l+1.d0)/((x(i)+h)**2.d0*T)-E+V(i,k)&
        +(j*(j+1.d0)-l*(l+1.d0)-3.d0/4.d0)*0.5d0*W(i,k))*(y(i,1)+d(3,1))
   d(4,1)=h*(y(i,2)+d(3,2))
   do n = 1,2
      y(i+1,n)=y(i,n)+(d(1,n)+2.0d0*d(2,n)+2.0d0*d(3,n)+d(4,n))/6.d0
   end do
end do
end subroutine fn
end module subprog

module node !  波動関数の解がiとi+1の間にあるか二分法によって判定する。
implicit none 
contains 
subroutine gn (y,P)
  real(8),intent(in)::y(2001,2)
  integer,intent(out)::P
  integer::i
  P =0
   do i = 2,2000
      if (y(i,1)*y(i+1,1)<0.d0) then
         P = P+1
       endif  
   enddo
end subroutine gn
end module node 

module normalize ! 波動関数を規格化する
implicit none
contains
subroutine hn(y,S,A,h)
real(8),intent(in)::h
real(8),intent(out)::S,A  
real(8),intent(inout)::y(2001,2)  
integer::i
do i = 2,2000
   S = S + h*y(i,1)**2.d0
end do
   S = S + h*y(1,1)**2 +h*y(2001,1)**2
A = S**(1/2.d0)
do i = 1,2001
   y(i,1)=y(i,1)/A
end do
end subroutine hn
end module normalize

program main  !N励起状態のEを計算する
  use subprog
  use node
  use normalize
  implicit none
 real(8) :: EL, EM, ES, h,S,A,l,R,r0,B,a0,Z,e,j,T,NA,sn,st
 real(8)::y(2001,2),f(2001),df(2001),x(2001),V(2001,2),W(2001,2)
 real(8) ::EN(3,2,7),yt(2001)
 integer :: i,k,ni,PL,PM,PS,li,n,ji,isospin,NAi,Zi,Ai,m,lt,g,dn
 EM = 0.d0
  h  = 3.d-3
  PL = 0
  PM = 0
  PS = 0
  li = 0
  ji = 0
  j  = 0.d0
  e  = 197.d0/137.d0
  r0 = 127.d-2
  a0 = 67.d-2
  T = (939.d0*2.d0/197.d0**2.d0)
  read(*,*)NAi,Zi
  NA =dble(NAi)
  Z = dble(Zi)
  A =NA+Z
 ! n = ni-li !n(node)
  R  = r0*A**(1.d0/3.d0)
  do n = 1,3 
     do ji = 1,2 
        do li = 1,7
     EN(n,ji,li)= 0.d0
        end do
     enddo
  enddo

  do i = 1,2001
   x(i)=h*dble(i)
    f(i)=1.d0/(1+exp((x(i)-R)/a0))
    df(i)= -exp((x(i)-R)/a0)/(1+exp((x(i)-R)/a0))**2.d0/a0
    do n =1,2
        y(i,n)=0.d0
     end do
  end do
do i = 1,2001
if (x(i) < R) then
   V(i,1) = -(51.d0+33.d0*(NA-Z)/A)*f(i) + (Z-1.d0)*e**2.d0*(3.d0-(x(i)/R)**2)/(2.d0*R)
     else
    V(i,1) = -(51.d0+33.d0*(NA-Z)/A)*f(i) + (Z-1.d0)*e**2.d0/x(i)+(l*(l+1.d0))
     endif
     V(i,2)  =-(51.d0-33.d0*(NA-Z)/A)*f(i)
     W(i,1) = (22.d0+14.d0*(NA-Z)/A)*r0**2.d0/x(i)*df(i)
     W(i,2) = (22.d0-14.d0*(NA-Z)/A)*r0**2.d0/x(i)*df(i)
  end do
  
  do i = 1,2001
   write (*,*) V(i,1) + W (i,1)
  enddo

  li = 1
  n =1
  isospin =1
  ji = 1
 do li = 1,7
   lt = li-1
   l = dble(lt)
  do n  = 1,3
    do ji = 1,2
      if (ji == 1) then
         j = l + 0.5d0
       else
         j =abs(l - 0.5d0)
      end if
    do isospin  = 1,2  !1=proton,2=nutron   
        EL =  1.d3
        ES = -5.d2
    do i=1,100
      EM = (EL + ES)/2.d0       
      call fn(EL,y,l,j,h,V,W,isospin,T)
      call gn(y,PL)
      call fn(EM,y,l,j,h,V,W,isospin,T)
      call gn(y,PM)
      if((dble(PL) -dble(n) +0.5d0)*(dble(PM) -dble(n) +0.5d0) > 0.d0) then
        EL = EM
      else
         ES = EM
      end if
      EN(n,ji,li) = EM
      ! write(*,*) EM
    end do
    if (isospin == 2 ) then
    if (EN(n,ji,li)<10.d0 ) then
      if(li== 1) then 
         sn = 1.d0
      else 
         sn = 2*j+1.d0
      endif
         st = st +sn
     call fn(EM,y,l,j,h,V,W,isospin,T)
     call hn(y,S,A,h)
     do dn =1,2001
    yt(dn) = yt(dn) + y(dn,1)**2.d0*sn/(4.d0*3.14d0*x(dn)**2.d0)
      end do
   !      write(*,*) EN(n,ji,li),",",n,",",j,",",lt,",",sn,",",st,","!!
        endif
    endif
 ! do i = 1,2001
   ! write(*,*) x(i),yt(i)
  ! end do  
           enddo
        enddo
      enddo  
   end do       
   do dn =1,2001
   ! write(*,*)x(dn),yt(dn)
   ! sn = sn + h*yt(dn)
    end do
  
end program main 
