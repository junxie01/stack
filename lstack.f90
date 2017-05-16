program stack
use sacio
implicit none
type(sac_head) :: sachead,sacheado
integer,parameter:: nsmax=4000000,maxtr=50000
character(18)par
character(180)filein,fileout,output,list
logical ext
integer do_n
integer nsamp,nerr,icor,nsmp,ntime,itr
integer narg,npts,nn,n21,n1,ns,nf,i,nfre,npow
real dt0
real t1,t2
real beg,dt,depmax
real sig(nsmax),rsig(nsmax)
real sig0(nsmax),hilb(nsmax)

narg=iargc()
call cpu_time(t1)
if(narg.lt.3 )stop "Usage: stack sac_file_list output_file_name do_norm"
call getarg(1,list)
call getarg(2,fileout)
call getarg(3,par)
read(par,'(bn,i20.0)')do_n

inquire(file=list,exist=ext)
if(.not.ext)stop 'SAC file list not found'
icor=0;sig=0.0;sig0=0.0
rsig=0.0
open(21,file=list)
write(*,*)'maxtrace is ',maxtr
do itr=1,maxtr
   read(21,'(a180)',end=80,err=80)filein
   call read_sachead(trim(filein),sachead,nerr)
   write(*,'(" Read ",i5.5, " file ",1a)')itr,trim(filein)
   if(nerr.eq.-1)cycle
   if(icor==0)dt0=sachead%delta
   write(*,*)'dt=',dt0,'depmax=',sachead%depmax
   if(sachead%delta/=dt0)cycle
   if(sachead%depmax.ne.sachead%depmax)cycle ! data null
   if(abs(sachead%depmax).lt.1e-6)cycle      ! data null
   nsmp=sachead%npts 
   call read_sac(trim(filein),sig0,sachead,nerr)
   beg=sachead%b
   if(nerr.eq.-1)cycle   ! data read problem
   depmax=maxval(sig0(1:nsmp))
   write(*,*)'dt=',dt0,'depmax=',depmax
   if(abs(depmax).lt.1e-16)cycle
   if(depmax.ne.depmax)cycle
   if(do_n.ne.1)depmax=1.0
   rsig(1:nsmp)=sig0(1:nsmp)/depmax+rsig(1:nsmp) ! linear stack
   !write(*,*)depmax
   icor=icor+1
enddo   ! end loop over traces
80 close(21)

write(*,*) 'Number of the trace is: ',icor
if(icor.eq.0)stop 'No ncf found!'
call initial_sachead(sacheado)
sacheado%stla=sachead%stla
sacheado%stlo=sachead%stlo
sacheado%evla=sachead%evla
sacheado%evlo=sachead%evlo
sacheado%npts=nsmp
sacheado%b=beg
sacheado%unused13=icor
sacheado%delta=dt0
sacheado%o=0
output='tl_'//trim(fileout)
call write_sac(output,rsig,sacheado,nerr)
if(nerr.eq.-1)stop 'Error in writing linear stack file'
call cpu_time(t2)
write(*,*)"Total running time is ",t2-t1, "s"
end program
