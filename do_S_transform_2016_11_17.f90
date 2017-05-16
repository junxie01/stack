program main
use sacio
type (sac_head):: sachead
integer,parameter :: nsmax=4000000
complex cdum(nsmax)
complex csig(nsmax)
integer npow,nn,n21,ns,nf,ddf
integer nerr,i,j,nsamp
integer nf1,mf,mn,fi,fb,fe
real dt,rdum,rr
real rmf,rns,f1,f2,df
real argp,argn,rpp,pp
real pi,fact,cyc,freq
real,allocatable,dimension(:) :: sig
character(180) :: sacfile,par,output
if(iargc().ne.1)stop 'Usage: do_S_tran sacfile'
call getarg(1,sacfile)
!call getarg(2,output)
call read_sachead(sacfile,sachead,nerr)
nsamp=sachead%npts
allocate(sig(nsamp))
!write(*,*)'no. of points is ',sachead%npts
call read_sac(sacfile,sig,sachead,nerr)
if(nerr.eq.-1)stop 'read sac file problem'
write(*,*)'Read file ',trim(sacfile)
cyc=4.0;npow=1;nn=2
dt=sachead%delta
do while(nn.lt.nsamp)
   npow=npow+1
   nn=nn*2
enddo
n21=nn/2
pi=4.0*atan(1.0)
fact=cyc/2.0
csig=cmplx(0.0d0,0d0)
csig(1:nsamp)=cmplx(sig(1:nsamp),0.0)
call clogc(npow,csig,1,dt)    ! do fft
pp=-2.0*pi*pi*fact*fact
f1=0.030;f2=0.045
df=0.001
df=1.0/nn/dt
fb=int(f1*nn*dt)
fe=int(f2*nn*dt)
ddf=int(df*nn*dt)
nf=int((f2-f1)/df)+1
open(15,file='a.out')
do fi=1, nf               ! loop over frequency
   nf1=(fi-1)*ddf+fb
   cdum=cmplx(0.0,0.0)
   rpp=pp/float(nf1)/float(nf1)
   rpp=-pi/nf1**2
   do mf=1,nn
      rmf=float(mf)
      rns=float(nn-mf+1)
      argp=rpp*rmf*rmf
      argn=rpp*rns*rns
      if (argp.gt.-25.or.argn.gt.-25) then
         mn=mf+nf1
         if(mn.gt.nn)mn=mn-nn
         rr=exp(argp)+exp(argn)
         rr=exp(argp)
         cdum(mf)=csig(mn)*rr
         !write(*,*) rr
      endif
   enddo
   call clogc(npow,cdum,-1,dt)
   !ctf(1:nsamp,nf)=cdum(1:nsamp)
   freq=1.0/nn/dt*((fi-1)*ddf+fb) 
   do j=1,nsamp
      write(15,*)dt*(j-1),freq,real(cdum(j))
   enddo
enddo
close(15)
deallocate(sig)
end program
