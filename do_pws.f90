subroutine do_pws(sig,nsmp,dt,pw_sig)
integer nsmax,nfmax
parameter(nsmax=8192,nfmax=4096)
real hilb(nsmax),sig(nsmax),dt,theta
integer nsmp,i
complex pw_sig(nsmax)
call hilbert(sig,nsmp,dt,hilb) ! do hilbert transform
do i=1,nsmp
   theta=atan2(hilb(i),sig(i))
   pw_sig(i)=pw_sig(i)+cmplx(cos(theta),sin(theta))
enddo
end subroutine
