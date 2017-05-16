program stack
use sacio
implicit none
type(sac_head) :: sachead
integer,parameter:: nsmax=8192,nfmax=4096,maxtr=50000
!parameter (nsmax=16384,nfmax=8192,maxtr=10000)
character(18)par
character(180)filein,fileout,output,list
logical ext
integer do_n
integer do_linear,do_pws,do_tf_pws
integer nsamp,nerr,icor,nsmp,ntime,itr
integer narg,npts,nn,n21,n1,ns,nf,i,nfre,npow
real rdum
real theta
real cyc,wu
real twpi,t1,t2
real beg,dt,depmax
real sig(nsmax),rsig(nsmax)
real sig0(nsmax),hilb(nsmax)
complex pw_sig(nsmax)
complex ctf(nsmax,nfmax)
complex carray1(nsmax),cdum,carg
complex ctf2(nsmax,nfmax),ctfps(nsmax,nfmax)
interface
    function diftc_fast_f(carray1,nn,n21,carg)
    integer nsmax,nfmax,maxtr,nn,n21
    parameter (nsmax=8192,nfmax=4096,maxtr=50000)
    complex carray1(nsmax),carg
    end function
end interface 

twpi=8.0*atan(1.0)
cyc=4.0;wu=2.0
narg=iargc()
call cpu_time(t1)
if(narg.lt.6 )then
   stop "Usage: stack sac_file_list output_file_name do_linear do_pws &
   do_tf_pws do_normalization"
endif
call getarg(1,list)
call getarg(2,fileout)
call getarg(3,par)
read(par,'(bn,i20.0)')do_linear
call getarg(4,par)
read(par,'(bn,i20.0)')do_pws
call getarg(5,par)
read(par,'(bn,i20.0)')do_tf_pws
call getarg(6,par)
read(par,'(bn,i20.0)')do_n

inquire(file=list,exist=ext)
if(.not.ext)stop 'SAC file list not found'
if(do_linear.eq.0.and.do_pws.eq.0.and.do_tf_pws.eq.0)stop 'Program will do nothing because of your parameter input'
if(do_linear.eq.1)write(*,*)'Do linear stacking'
if(do_pws.eq.1)write(*,*)'Do pws stacking'
if(do_tf_pws.eq.1)write(*,*)'Do tf_pws stacking'
icor=0;sig=0.0;sig0=0.0
rsig=0.0;pw_sig=0.0
ctf=cmplx(0.0,0.0)
ctf2=cmplx(0.0,0.0)
ctfps=cmplx(0.0,0.0)
open(20,file=list,status="old")
do itr=1,maxtr
   read(20,'(a180)',end=80,err=80)filein
   call read_sachead(trim(filein),sachead,nerr)
   write(*,'(" Read ",i5.5, " file ",1a)')itr,trim(filein)
   if(nerr.eq.-1)cycle
   if(sachead%depmax.ne.sachead%depmax)cycle ! data null
   nsmp=sachead%npts 
   call read_sac(trim(filein),sig0,sachead,nerr)
   beg=sachead%b
   if(nerr.eq.-1)cycle   ! data read problem
   icor=icor+1
   sig=sig0
   if(nsmp.gt.nsmax)then ! if the number of samples is too big
      n1=int(sachead%dist/3.0)+2000
      nsmp=2*n1+1
      if(nsmp.gt.nsmax)n1=(nsmax-1)/2
      nsmp=2*n1+1
      n21=(nsmp-1)/2
      sig(1:nsmp)=sig0(n21-n1:n21+n1)
      beg=-real(n1)*sachead%delta
   endif
   depmax=0.0
   if(do_n.eq.1)then
      do i=1,nsmp
         depmax=max(depmax,sig(i))
      enddo
   else
      depmax=1.0
   endif
   rsig(1:nsmp)=sig(1:nsmp)/depmax+rsig(1:nsmp) ! linear stack
   if(do_tf_pws.eq.1)then
      if(icor.eq.1)then
         npow=1
         nn=2
         do while(nn.lt.nsmp)
            npow=npow+1
            nn=nn*2
         enddo
         n21=nn/2
      endif
      call S_trans(sig,nsmp,nn,n21,npow,cyc,dt,ctf)
      do ns=1,nn
         do nf=1,n21
            cdum=ctf(ns,nf)
            ctf2(ns,nf)=ctf2(ns,nf)+cdum
            rdum=cabs(cdum)
            if(rdum.lt.0.00001) then
               cdum=cmplx(0.,0.)
            else
               cdum=cdum/(cabs(cdum))
            endif
            ctfps(ns,nf)=ctfps(ns,nf)+cdum
         enddo
      enddo
   endif ! do tf PWS
   if(do_pws.eq.1)then
      call hilbert(sig,nsmp,dt,hilb)              ! do hilbert transform
      do i=1,nsmp
         theta=atan2(hilb(i),sig(i))
         pw_sig(i)=pw_sig(i)+cmplx(cos(theta),sin(theta))
      enddo
   endif
enddo   ! end loop over traces
80 close(20)

write(*,*) 'Number of the trace is: ',itr-1
if(icor.eq.0)stop 'No ncf found!'
sachead%npts=nsmp
sachead%b=beg
sachead%unused13=icor
if(do_tf_pws.eq.1)then
   output='tf_pws_'//trim(fileout)
   do ns=1,nn
      do nf=1,n21
         cdum=ctf2(ns,nf)/real(icor)
         rdum=cabs(ctfps(ns,nf))/real(icor)
         ctf(ns,nf)=cdum*(rdum**wu)
      enddo
   enddo
   cdum=cmplx(0.0,-twpi/float(nn))
   rdum=float(nn)*twpi*2.
   do ntime=1,nsmp
      carg=cdum*(ntime-1)
      do nfre=1,n21
         carray1(nfre)=ctf(ntime,nfre)*rdum
      enddo
      sig(ntime)=diftc_fast_f(carray1,nn,n21,carg)
   enddo
   call write_sac(output,sig,sachead,nerr)
   if(nerr.eq.-1) stop 'Error in writing linear stack file'
endif
if(do_linear.eq.1)then
   output='tl_'//trim(fileout)
   call write_sac(output,rsig,sachead,nerr)
   if(nerr.eq.-1)stop 'Error in writing linear stack file'
endif
if(do_pws.eq.1)then
   output='pws_'//trim(fileout)
   do i=1,nsmp 
      rsig(i)=(cabs(pw_sig(i))/real(icor))**wu*rsig(i)
   enddo
   call write_sac(output,rsig,sachead,nerr)
   if(nerr.eq.-1)stop 'Error in writing pws file'
endif
call cpu_time(t2)
write(*,*)"Total running time is ",t2-t1, "s"
end program
