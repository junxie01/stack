program main
use sacio
type (sac_head):: sachead
integer nsmax,nfmax
parameter(nsmax=8192,nfmax=4096)
integer npow,nn,n21,ns,nf
integer nerr,i,j,nsamp
complex ctf(nsmax,nfmax),cdum
complex ctf2(nsmax,nfmax),ctfps(nsmax,nfmax)
real cyc,dt,rdum,freq
real,allocatable,dimension(:) :: sig
character(180) :: sacfile,par,output
if(iargc().ne.1)stop 'Usage: do_S_tran sacfile'
call getarg(1,sacfile)
!call getarg(2,output)
call read_sachead(sacfile,sachead,nerr)
nsamp=sachead%npts
allocate(sig(nsamp))
write(*,*)'no. of points is ',sachead%npts
call read_sac(sacfile,sachead,sig,nerr)
if(nerr.eq.-1)stop 'read sac file problem'
cyc=4.0
npow=1
dt=sachead%delta
nn=2
do while(nn.lt.nsamp)
   npow=npow+1
   nn=nn*2
enddo
n21=nn/2
call S_trans(sig,nsamp,nn,n21,npow,cyc,dt,ctf)

!open(15,file=output)
!do i=2,n21
!   freq=1.0/nn/dt*(i-1) 
!   do j=1,nsamp
!      write(15,*)freq,dt*j-1,real(ctf(j,i))
!   enddo
!enddo
!close(15)
deallocate(sig)
end program
