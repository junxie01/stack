REAL FUNCTION diftc_fast_f(y,nn,n21,carg)
!
! Discrete Inverse Fourier Transformation / |f|:
! (Note: constant factors are not included.)
!
!  y        is array with complex spectrum (positive freqs).
!  nsmp     is # of samples in time domain.
!  nsmp2    is the length of data vector y (pos. freqs).
!  ntime    is desired time index.
!  carg     is cmplx(0.0,twpi*ntime/nsmp)
!===============================================================
complex y(1),carg,yy,cdum,cphs1,cphs2
integer nn,n21

cdum=cmplx(0.D0,0.D0)
do nm=1,n21
        yy=y(nm)/float(nm)
        cphs1=yy*cexp(carg*(nm-1))
        cphs2=conjg(yy)*cexp(carg*(nn-nm+1))
        cdum=cdum+cphs1+cphs2
enddo
diftc_fast_f=real(cdum)/float(nn)
end
!---------------------------------------------------------------

