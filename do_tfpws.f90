subroutine do_tfpws(sig,nsmp,nn,n21,dt,ctf2,ctfps)
integer nsmax,nfmax
parameter(nsmax=8192,nfmax=4096)
integer npow,nn,n21,nsamp,ns,nf
complex ctf(nsmax,nfmax),cdum
complex ctf2(nsmax,nfmax),ctfps(nsmax,nfmax)
real sig(nsmax),cyc,dt,rdum
cyc=4.0
npow=1
nn=2
do while(nn.lt.nsmp)
   npow=npow+1
   nn=nn*2
enddo
n21=nn/2
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
return
end subroutine
