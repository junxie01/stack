subroutine taperf(sf,fre,dom,nk,npow)
integer *4 nmax
parameter(nmax=4000000)
real*4    fre(4)
real*8    dom,f
real*8    d1,d2,dpi,ss,s(nmax)
integer   i,j,nk,npow
complex   sf(nmax)
dpi = datan(1.0d0)*4.0d0
s=0d0
do i = 1,nk
   f = dble(i-1)*dom
   if(f.le.dble(fre(1))) then
      cycle
   else if(f.le.dble(fre(2))) then
      d1 = dpi/(dble(fre(2))-dble(fre(1)))
      ss = 1.0d0
      do j = 1,npow
         ss = ss*(1.d0-dcos(d1*(dble(fre(1))-f)))/2.0d0
      enddo
      s(i) = ss
   else if(f.le.dble(fre(3))) then
      s(i) = 1.0d0
   else if(f.le.dble(fre(4))) then
      d2 = dpi/(dble(fre(4))-dble(fre(3)))
      ss = 1.0d0
      do j = 1,npow
         ss = ss*(1.0d0+dcos(d2*(fre(3)-f)))/2.0d0
      enddo
      s(i) = ss
   endif
enddo
sf(1:nk) = sf(1:nk)*sngl(s(1:nk))
return
end subroutine
