subroutine filter(sig,fre,dt,npts)
implicit none
integer*4 nmax
parameter(nmax=4000000)
integer*4 k,nk,npts,ns,npow
real*4    fre(4),dt,sig(nmax),sigout(nmax)
real*8    dom
complex   czero,s(nmax)
czero = (0.0,0.0)
ns=2;npow=1
do while(ns.lt.npts)
   ns=ns*2
   npow=npow+1
enddo
dom = dble(1.0/dt/ns)
s=czero
nk = ns/2+1
s(1:npts) = cmplx(sig(1:npts),0)
call clogc(npow,s,1,dt)
s(nk+1:ns) = czero
call taperf(s,fre,dom,nk,npow)
call clogc(npow,s,-1,dt)
sig(1:npts)=real(s(1:npts))
return
end subroutine
