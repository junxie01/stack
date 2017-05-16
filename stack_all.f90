program stack_all
use sacio
implicit none
type(sac_head) :: sachead1,sachead2
integer nsmax,nfmax,maxtr,nstmax
parameter (nsmax=8192,nfmax=4096,maxtr=50000,nstmax=20000)
!parameter (nsmax=16384,nfmax=8192,maxtr=10000)
character(18) :: par,year_day
character(180):: sac_ncc,sac_pcc,input,command
character(180):: filein,fileout,output,list,dirinn,dirout
character(8) ::  sta1(nstmax),sta2(nstmax)
character(2) ::  net1(nstmax),net2(nstmax),nd
character(3) ::  com(3),comm
logical ext
integer do_n,nlen,year_b,year_e,nh,nnh,nseg,do_fil
integer day_b,day_e,is,iy,jday,ih,ddhour,multpt,iseg
integer narg,npts,nn,n21,n1,ns,nf,i,nfre,npow,id,ic1,ic2
integer dotl,dopws,dotfpws,doncc,dopcc,dhour,cn,nerr1,nerr2
integer nsamp,nerr,icor,nsmp,ntime,itr,endday,begday,nptseg,nst
real rdum
real theta
real cyc,wu
real twpi,t1,t2
real sig01(nsmax),dseg
real sig11(nsmax),sig(nsmax)
real beg,dt,depmax1,depmax2
real sig02(nsmax),hilb(nsmax),fre(4)
real sig12(nsmax),rsig1(nsmax),rsig2(nsmax)
complex pw_sig1(nsmax)
complex pw_sig2(nsmax)
complex ctfps(nsmax,nfmax),cdum,carg
complex ctfps2(nsmax,nfmax),carray1(nsmax)
complex ctf11(nsmax,nfmax),ctf12(nsmax,nfmax)
complex ctf21(nsmax,nfmax),ctf22(nsmax,nfmax)
complex ctfps1(nsmax,nfmax),ctf(nsmax,nfmax)
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
if (iargc().ne.1)then
   write(*,*)'Usage: Stack_all param.dat '
   write(*,*)'param.dat is like:'
   write(*,*)'tobecomputed.list'
   write(*,*)'year_b day_b year_e day_e'
   write(*,*)'cn,com,nlen,do_n'
   write(*,*)'dhour,ddhour,multpt'
   write(*,*)'f1,f2,do_fil'
   write(*,*)'doncc,dopcc'
   write(*,*)'dotl,dopws,dotfpws'
   write(*,*)'/inputdir/'
   write(*,*)'/outputdir/'
   stop
endif
call getarg(1,input)
open(10,file=input)
read(10,*)list
read(10,*)year_b,day_b,year_e,day_e
read(10,*)cn,comm,nlen,do_n
read(10,*)dhour,ddhour,multpt
read(10,*)fre(2),fre(3),do_fil
read(10,*)doncc,dopcc
read(10,*)dotl,dopws,dotfpws
read(10,*)dirinn
read(10,*)dirout
close(10)

inquire(file=list,exist=ext)
if(.not.ext)stop "Station pair list not found!"
!write(*,'("mkdir ",1a)') dirout
write(command,'("mkdir ",1a)') trim(dirout)
call system(command)
if(dotl.eq.0.and.dopws.eq.0.and.dotfpws.eq.0)stop 'Program will do nothing because of your parameter input'
if(dotl.eq.1)   write(*,*)'Do linear stacking'
if(dopws.eq.1)  write(*,*)'Do pws stacking'
if(dotfpws.eq.1)write(*,*)'Do tf_pws stacking'
if(cn.eq.3)then
   com(1)=trim(comm)//'Z'
   com(2)=trim(comm)//'N'
   com(3)=trim(comm)//'E'
else
   com(1)=comm
endif
nptseg=2*nlen+1
open(11,file=list)   ! read in station pairs
do i=1,nstmax
   read(11,*,err=13,end=13) net1(i),sta1(i),net2(i),sta2(i)
enddo
13 close(11)
nst=i-1
nh=24/dhour
fre(1)=0.95*fre(2)
fre(4)=1.05*fre(3)
write(*,*)'No of station paris is: ',nst
write(*,*)'No of the points is: ',nptseg
do is=1,nst  ! loop over station pair
   icor=0;
   do ic1=1,cn        ! loop over first component       
      do ic2=1,cn     ! loop over second components 
         ctf12=cmplx(0.0,0.0)
         ctf22=cmplx(0.0,0.0)
         ctfps=cmplx(0.0,0.0)
         sig11=0.0;sig01=0.0
         sig12=0.0;sig02=0.0
         rsig1=0.0;pw_sig1=cmplx(0.0,0.0)
         rsig2=0.0;pw_sig2=cmplx(0.0,0.0)
         do iy=year_b,year_e  ! loop over year
            jday=365
            if(mod(iy,4).eq.0.and.mod(iy,100).ne.0.or.mod(iy,400).eq.0)jday=366
            endday=day_e
            if(iy.ne.year_e)endday=jday
            begday=day_b
            if(iy.ne.year_b)begday=1
            do id=begday,endday ! loop over day
               write(year_day,'(i0,"_",i3.3)')iy,id
               do ih=1,nh            ! loop over hour segment
                  nnh=(ih-1)*dhour 
                  write(nd,'(i2.2)')nnh
                  dseg=(1.0-real(multpt)/100.0)*ddhour ! the left points without overlapping
                  nseg=int((dhour-ddhour)/dseg)+1 ! number of segments
                  do iseg=1,nseg
                     write(sac_ncc,'(1a,"/",1a,"_",1a,"/cc_",1a,"_",1a,"_",1a,"_",1a,"_",1a,"_",1a,"_",1a,"_",1a,".SAC_",i2.2)')&
                     trim(dirinn),trim(sta1(is)),trim(sta2(is)),trim(year_day),trim(nd),trim(net1(is)),trim(sta1(is)),&
                     trim(net2(is)),trim(sta2(is)),trim(com(ic1)),trim(com(ic2)),iseg
                     write(sac_pcc,'(1a,"/",1a,"_",1a,"/pcc_",1a,"_",1a,"_",1a,"_",1a,"_",1a,"_",1a,"_",1a,"_",1a,".SAC_",i2.2)')&
                     trim(dirinn),trim(sta1(is)),trim(sta2(is)),trim(year_day),trim(nd),trim(net1(is)),trim(sta1(is)),&
                     trim(net2(is)),trim(sta2(is)),trim(com(ic1)),trim(com(ic2)),iseg
                     write(*,*)trim(sac_ncc),trim(sac_pcc)
                     nerr1=0;nerr2=0
                     if(doncc.eq.1)call read_sachead(trim(sac_ncc), sachead1,nerr1)
                     if(dopcc.eq.1)call read_sachead(trim(sac_pcc),sachead2,nerr2)
                     write(*,'(" Read file ",1a)')trim(sac_pcc)
                     if(nerr1.eq.-1)cycle
                     if(nerr2.eq.-1)cycle
                     if(sachead1%depmax.ne.sachead1%depmax)cycle ! data null
                     if(sachead2%depmax.ne.sachead2%depmax)cycle ! data null
                     nsmp=sachead1%npts 
                     if(doncc.eq.1)call read_sac(trim(sac_ncc), sig01,sachead1,nerr1)
                     if(dopcc.eq.1)call read_sac(trim(sac_pcc),sig02,sachead2,nerr2)
                     beg=sachead1%b
                     if(nerr1.eq.-1)cycle   ! data read problem
                     if(nerr2.eq.-1)cycle   ! data read problem
                     write(*,*)'hello beg=',beg
                     icor=icor+1
                     if(doncc.eq.1 )sig11=sig01
                     if(dopcc.eq.1)sig12=sig02
                     if(nsmp.gt.nsmax)then ! if the number of samples is too big
                        n1=int(sachead1%dist/3.0)+2000
                        nsmp=2*n1+1
                        if(nsmp.gt.nsmax)n1=(nsmax-1)/2
                        nsmp=2*n1+1
                        n21=(nsmp-1)/2
                        if(doncc.eq.1)sig11(1:nsmp)=sig01(n21-n1:n21+n1)
                        if(dopcc.eq.1)sig12(1:nsmp)=sig02(n21-n1:n21+n1)
                        beg=-real(n1)*sachead1%delta
                     endif
                     depmax1=0.0;depmax2=0.0
                     if(doncc.eq.1.and.do_fil.eq.1)call filter(sig11,fre,sachead1%delta,nsmp)
                     if(dopcc.eq.1.and.do_fil.eq.1)call filter(sig12,fre,sachead1%delta,nsmp)
                     if(do_n.eq.1)then
                        do i=1,nsmp
                           if(doncc.eq.1)depmax1=max(depmax1,sig11(i))
                           if(dopcc.eq.1)depmax2=max(depmax2,sig12(i))
                        enddo
                     else
                        depmax1=1.0
                        depmax2=1.0
                     endif
                     dt=sachead1%delta
                     if(doncc.eq.1)rsig1(1:nsmp)=sig11(1:nsmp)/depmax1+rsig1(1:nsmp) ! linear stack
                     if(dopcc.eq.1)rsig2(1:nsmp)=sig12(1:nsmp)/depmax2+rsig2(1:nsmp) ! linear stack
                     if(dotfpws.eq.1)then
                        if(doncc.eq.1)call do_tfpws(sig11,nsmp,nn,n21,dt,ctf12,ctfps1)
                        if(dopcc.eq.1)call do_tfpws(sig12,nsmp,nn,n21,dt,ctf22,ctfps2)
                        write(*,*)'hello'
                     endif ! do tf PWS
                     if(dopws.eq.1)then
                        if(doncc.eq.1)call do_pws(sig11,nsmp,dt,pw_sig1)
                        if(dopcc.eq.1)call do_pws(sig12,nsmp,dt,pw_sig2)
                     endif
                  enddo  ! end loop over segments
               enddo           ! end loop over hour segments
            enddo              ! end loop over days
         enddo                 ! end loop over years
         write(*,*) 'Number of the trace is: ',icor
         if(icor.eq.0)cycle
         sachead1%npts=nsmp
         sachead1%b=beg
         sachead1%unused13=icor
         sachead2%npts=nsmp
         sachead2%b=beg
         sachead2%unused13=icor
         write(*,*)'hello10'
         if(dotfpws.eq.1)then
            if(doncc.eq.1)then
               write(output,'(1a,"/tf_ncc_",1a,"_",1a,"_",1a,"_",1a,"_",1a,"_",1a,".SAC")')&
               trim(dirout),net1(is),trim(sta1(is)),net2(is),trim(sta2(is)),com(ic1),com(ic2)
               do ns=1,nn
                  do nf=1,n21
                     cdum=ctf12(ns,nf)/real(icor)
                     rdum=cabs(ctfps1(ns,nf))/real(icor)
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
               call write_sac(output,sig,sachead1,nerr)
               if(nerr.eq.-1) stop 'Error in writing linear stack file'
            endif
            if(dopcc.eq.1)then
               write(output,'(1a,"/tf_pcc_",1a,"_",1a,"_",1a,"_",1a,"_",1a,"_",1a,".SAC")')&
               trim(dirout),net1(is),trim(sta1(is)),net2(is),trim(sta2(is)),com(ic1),com(ic2)
               do ns=1,nn
                  do nf=1,n21
                     cdum=ctf22(ns,nf)/real(icor)
                     rdum=cabs(ctfps2(ns,nf))/real(icor)
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
               call write_sac(output,sig,sachead2,nerr)
               if(nerr.eq.-1) stop 'Error in writing linear stack file'
            endif
         endif
         if(dotl.eq.1)then
            if(doncc.eq.1)then
               write(output,'(1a,"/tl_ncc_",1a,"_",1a,"_",1a,"_",1a,"_",1a,"_",1a,".SAC")')&
               trim(dirout),net1(is),trim(sta1(is)),net2(is),trim(sta2(is)),com(ic1),com(ic2)
               call write_sac(output,rsig1,sachead1,nerr)
               if(nerr.eq.-1)stop 'Error in writing linear stack file'
            endif
            if(dopcc.eq.1)then
               write(output,'(1a,"/tl_pcc_",1a,"_",1a,"_",1a,"_",1a,"_",1a,"_",1a,".SAC")')&
               trim(dirout),net1(is),trim(sta1(is)),net2(is),trim(sta2(is)),com(ic1),com(ic2)
               call write_sac(output,rsig2,sachead2,nerr)
               if(nerr.eq.-1)stop 'Error in writing linear stack file'
            endif
         endif
         if(dopws.eq.1)then
            if(doncc.eq.1)then
               write(output,'(1a,"/pw_ncc_",1a,"_",1a,"_",1a,"_",1a,"_",1a,"_",1a,".SAC")')&
               trim(dirout),net1(is),trim(sta1(is)),net2(is),trim(sta2(is)),com(ic1),com(ic2)
               do i=1,nsmp 
                  rsig1(i)=(cabs(pw_sig1(i))/real(icor))**wu*rsig1(i)
               enddo
               call write_sac(output,rsig1,sachead1,nerr)
               if(nerr.eq.-1)stop 'Error in writing pws file'
            endif
            if(dopcc.eq.1)then
               write(output,'(1a,"/pw_pcc_",1a,"_",1a,"_",1a,"_",1a,"_",1a,"_",1a,".SAC")')&
               trim(dirout),net1(is),trim(sta1(is)),net2(is),trim(sta2(is)),com(ic1),com(ic2)
               do i=1,nsmp 
                  rsig2(i)=(cabs(pw_sig2(i))/real(icor))**wu*rsig2(i)
               enddo
               call write_sac(output,rsig2,sachead2,nerr)
               if(nerr.eq.-1)stop 'Error in writing pws file'
            endif
         endif
      enddo     ! end loop over component1
   enddo        ! end loop over component2
enddo  ! end loop over station pairs
call cpu_time(t2)
write(*,*)"Total running time is ",t2-t1, "s"
end program
