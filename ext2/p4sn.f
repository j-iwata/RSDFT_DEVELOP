C************************* NEW FFT **********************
      subroutine fft3bx(nrx,nry,nrz,ng,rhog,work,
     &                  wsavex,wsavey,wsavez,ifacx,ifacy,ifacz,
     &                  lx1,lx2,ly1,ly2,lz1,lz2)
C***********************************************************
C     3-D Fourier transformation ( G-SPACE--> REAL SPACE)
C                                   (1990-04-12) OSAMU SUGINO
C     INPUT :RHO,NR?,NG,WSAVE?,IFAC?,L??
C     OUTPUT:RHOG
C     WORK  :WORK
C**********************************************************
C nac -> ncc (1995-05-26) Jun Yamauchi
C
      implicit real*8 (a-h,o-z)
      complex*16 wsavex(nrx),wsavey(nry),wsavez(nrz)
      dimension rhog(2,ng),work(2,ng)
      dimension ifacx(30),ifacy(30),ifacz(30)
      dimension lx1(ng),lx2(ng),ly1(ng),ly2(ng),lz1(ng),lz2(ng)
C
C     DO 15 IG=1,NG
C     WORK(1,IG)=RHOG(1,IG)
C     WORK(2,IG)=RHOG(2,IG)
C  15 CONTINUE
      call fftsv1(ng,rhog,work)
      call cfft3b(ng,nrx*nry,nrz,work,rhog,wsavez,ifacz)
C
      call fftxyz(ng,nrx*nry,nrz,work,rhog,lz1,lz2)
      call cfft3b(ng,nrz*nrx,nry,rhog,work,wsavey,ifacy)
C
      call fftxyz(ng,nrz*nrx,nry,rhog,work,ly1,ly2)
      call cfft3b(ng,nry*nrz,nrx,work,rhog,wsavex,ifacx)
C
      call fftxyz(ng,nry*nrz,nrx,work,rhog,lx1,lx2)
      call fftsv2(ng,rhog,work)
C     FAC=1.0D0/DBLE(NG)
      fac=1.0d0
      do 40 i=1,ng
      rhog(1,i)= work(1,i)*fac
      rhog(2,i)= work(2,i)*fac
   40 continue
C
      return
      end
C***********************************************************
      subroutine fft3fx(nrx,nry,nrz,ng,rhog,work,
     &                  wsavex,wsavey,wsavez,ifacx,ifacy,ifacz,
     &                  lx1,lx2,ly1,ly2,lz1,lz2)
C***********************************************************
C     3-D Fourier transformation (REAL-SPACE --> G-SPACE)
C                                   (1990-04-12) OSAMU SUGINO
C     INPUT :RHOG,NR?,NG,WSAVE?,IFAC?,L??
C     OUTPUT:WORK
C     WORK  :NONE
C
      implicit real*8 (a-h,o-z)
      dimension  rhog(2,ng),work(2,ng)
      complex*16 wsavex(nrx),wsavey(nry),wsavez(nrz)
      dimension ifacx(30),ifacy(30),ifacz(30)
      dimension lx1(ng),lx2(ng),ly1(ng),ly2(ng),lz1(ng),lz2(ng)
C
      call fftsv1(ng,rhog,work)
      call cfft3f(ng,nrx*nry,nrz,work,rhog,wsavez,ifacz)
C
      call fftxyz(ng,nrx*nry,nrz,work,rhog,lz1,lz2)
      call cfft3f(ng,nrz*nrx,nry,rhog,work,wsavey,ifacy)
C
      call fftxyz(ng,nrz*nrx,nry,rhog,work,ly1,ly2)
      call cfft3f(ng,nry*nrz,nrx,work,rhog,wsavex,ifacx)
C
      call fftxyz(ng,nry*nrz,nrx,work,rhog,lx1,lx2)
      call fftsv2(ng,rhog,work)
      fac=1.0d0/dble(ng)
C     FAC=1.0D0
      do 40 i=1,ng
      rhog(1,i)= work(1,i)*fac
      rhog(2,i)= work(2,i)*fac
   40 continue
      return
      end
C
      subroutine fftsv1(ng,rhog,work)
      implicit real*8 (a-h,o-z)
      dimension rhog(2,ng),work(ng,2)
      do 10 i=1,ng
      work(i,1)=rhog(1,i)
   10 work(i,2)=rhog(2,i)
      return
      end
      subroutine fftsv2(ng,rhog,work)
      implicit real*8 (a-h,o-z)
      dimension rhog(ng,2),work(2,ng)
      do 10 i=1,ng
      work(1,i)=rhog(i,1)
   10 work(2,i)=rhog(i,2)
      return
      end
C***********************************************************
      subroutine fftxyz(ng,nrxy,nrz,work,rhog,lt1,lt2)
C***********************************************************
C                                   (1990-04-12) OSAMU SUGINO
      implicit real*8 (a-h,o-z)
      dimension work(nrxy,nrz,2),rhog(ng,2),lt1(ng),lt2(ng)
C
      do 20 i=1,ng
      rhog(i,1)=work(lt2(i),lt1(i),1)
      rhog(i,2)=work(lt2(i),lt1(i),2)
   20 continue
      return
      end
C***************************************************************
C*****************FFT3D*****************************************
C SUBROUTINE PACKAGE OF THE 3-DIMENSIONAL FFT
C WRITTEN BY OSAMU SUGINO 1990-03-14 FIRST VERSION
C *************
C THIS IS A VECTORIZED VERSION OF THE NCARL FFTPACK IN U-TOKYO
C ( FTPKCV     VERSION T03                       1985-11-18 )
C ( PART OF VECTORIZED VERSION OF FFTPACK IN NCARL.         )
C *************
C BY THIS ROUTINE THE THIRD VARIABLE IS FOURIER-TRANSFORMED.
C SO, YOU HAVE TO SUCCESSIVELY CALL THIS ROUTINE 3-TIMES FOR THE
C COMPLETE 3-DIMENSIONAL FFT.
C     **********************************************************
C     NOTE:   CPU_TIME=1.6E-3 SEC (FOR SX-2)
C
C ****************USAGE****************************************
C
C CFFT3I:EQUIPMENT ROUTINE
C       N:DIMENSION OF THE THIRD VARIABLE
C       WSAVE:ARRAY[1..2*N] OF REAL*8; WORK
C       IFAC:ARRAY[1..30] OF INTEGER*4; WORK
C       ************************
C       CALL CFFT3I(N,WSAVE,IFAC)
C       ************************
C
C
C CFFT3B:REAL-SPACE --> K-SPACE FFT
C       NRXY:PRODUCT OF THE DIMENSION OF THE FIRST AND THE SECOND
C            VARIABLES
C       N:DIMENSION OF THE THIRD VARIABLE
C       C:ARRAY[1..NRXY,1..N] OF COMPLEX*16; INPUT AND OUTPUT
C       WORK:ARRAY[1..NRXY,1..N] OF COMPLEX*16; WORK
C       WSAVE:ARRAY[1..2*N] OF REAL*8; WORK
C       IFAC:ARRAY[1..30] OF INTEGER*4; WORK
C       ************************************
C       CALL CFFT3B(NRXY,N,C,WORK,WSAVE,IFAC)
C       ************************************
C
C
C CFFT3F:K-SPACE --> REAL-SPACE FFT
C       NRXY:PRODUCT OF THE DIMENSION OF THE FIRST AND THE SECOND
C            VARIABLES
C       N:DIMENSION OF THE THIRD VARIABLE
C       C:ARRAY[1..NRXY,1..N] OF COMPLEX*16; INPUT AND OUTPUT
C       WORK:ARRAY[1..NRXY,1..N] OF COMPLEX*16; WORK
C       WSAVE:ARRAY[1..2*N] OF REAL*8; WORK
C       IFAC:ARRAY[1..30] OF INTEGER*4; WORK
C       ************************************
C       CALL CFFT3F(NRXY,N,C,WORK,WSAVE,IFAC)
C       ************************************
C
C
C FOR FUTURE PROGRAMMING:
C     THIS ROUTINE MAY BECOME FASTER IF YOU CHANGE
C
*     DO 203 I=1,IDO
*     DO 203 IJ=1,NRXY
C
C     TO
C
*     DO 203 IL=1,IDO*NRXY
*     I=(IL-1)/NRXY+1
*     IJ=MOD(IL-1,NRXY)+1
C
C     BUT THIS WAS NOT AN IMPROVEMENT FOR THE SX-2 SUPER-COMPUTER
C     SUSTEM
C****************************************************************
      subroutine prefft(nrx,nry,nrz,ng,wsavex,wsavey,wsavez,
     & ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)
      implicit real*8 (a-h,o-z)
      complex*16 wsavex(nrx),wsavey(nry),wsavez(nrz)
      dimension ifacx(30),ifacy(30),ifacz(30)
      dimension lx1(ng),lx2(ng),ly1(ng),ly2(ng),lz1(ng),lz2(ng)
C
C        MAKE LIST VECTORS
C
      ij=0
      do 961 j=1,nrx*nry
      do 961 i=1,nrz
      ij=ij+1
      lz1(ij)=i
  961 lz2(ij)=j
      ij=0
      do 962 j=1,nry*nrz
      do 962 i=1,nrx
      ij=ij+1
      lx1(ij)=i
  962 lx2(ij)=j
      ij=0
      do 963 j=1,nrz*nrx
      do 963 i=1,nry
      ij=ij+1
      ly1(ij)=i
  963 ly2(ij)=j
      call cfft3i(nrx,wsavex,ifacx)
      call cfft3i(nry,wsavey,ifacy)
      call cfft3i(nrz,wsavez,ifacz)
      return
      end
      subroutine cfft3i(n,wsave,ifac)
      implicit real*8(a-h,o-z)
      dimension wsave(*),ifac(*)
      if (n .eq. 1) return
      call cft3i1 (n,wsave,ifac)
      return
      end
      subroutine cft3i1 (n,wa,ifac)
      implicit real*8(a-h,o-z)
      dimension wa(*),ifac(*),ntryh(4)
      data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/3,4,2,5/
      data tpi/6.28318530717959d0/
      nl = n
      nf = 0
      j = 0
  101 j = j+1
      if (j-4) 102,102,103
  102 ntry = ntryh(j)
      go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr) 101,105,101
  105 nf = nf+1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         ifac(ib+2) = ifac(ib+1)
  106 continue
      ifac(3) = 2
  107 if (nl .ne. 1) go to 104
      ifac(1) = n
      ifac(2) = nf
      argh = tpi/dble(n)
      i2 = 2
      l1 = 1
      do 110 k1=1,nf
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do 109 j=1,ipm
            i1 = i2
            wa(i1-1) = 1.0d0
            wa(i1)   = 0.0d0
            ld = ld+l1
            argld = dble(ld)*argh
            do 108 ifi=1,ido
               arg = dble(ifi)*argld
               wa(2*ifi+i1-1) = cos(arg)
               wa(2*ifi+i1) = sin(arg)
  108       continue
            i2 = i1+ido+ido
            if (ip .le. 5) go to 109
            wa(i1-1) = wa(i2-1)
            wa(i1)   = wa(i2)
  109    continue
         l1 = l2
  110 continue
      return
      end
      subroutine cfft3b(ng,nrxy,n,c,work,wsave,ifac)
      implicit real*8(a-h,o-z)
      dimension c(ng,2),wsave(*),work(ng,2),ifac(*)
      if (n .eq. 1) return
      call cft3b1 (nrxy,n,c(1,1),c(1,2),work(1,1),work(1,2),wsave,ifac)
      return
      end
      subroutine cft3b1 (nrxy,n,cr,ci,chr,chi,wa,ifac)
      implicit real*8(a-h,o-z)
      dimension chi(*),chr(*),cr(*),ci(*),wa(*),ifac(*)
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido
         idl1 = idot*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+idot*2
         ix3 = ix2+idot*2
         if (na .ne. 0) go to 101
         call thdb4v(nrxy,idot,l1,cr,ci,chr,chi,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call thdb4v(nrxy,idot,l1,chr,chi,cr,ci,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call thdb2v(nrxy,idot,l1,cr,ci,chr,chi,wa(iw))
         go to 105
  104    call thdb2v(nrxy,idot,l1,chr,chi,cr,ci,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip .ne. 3) go to 109
         ix2 = iw+idot*2
         if (na .ne. 0) go to 107
         call thdb3v(nrxy,idot,l1,cr,ci,chr,chi,wa(iw),wa(ix2))
         go to 108
  107    call thdb3v(nrxy,idot,l1,chr,chi,cr,ci,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip .ne. 5) go to 112
         ix2 = iw+idot*2
         ix3 = ix2+idot*2
         ix4 = ix3+idot*2
         if (na .ne. 0) go to 110
         call thdb5v(nrxy,idot,l1,cr,ci,chr,chi,
     &               wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call thdb5v(nrxy,idot,l1,chr,chi,cr,ci,
     &               wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na .ne. 0) go to 113
         call thdbg (nrxy,ncc,idot,ip,l1,idl1,cr,ci,cr,ci,cr,ci,
     &               chr,chi,chr,chi,wa(iw))
         go to 114
  113    call thdbg (nrxy,ncc,idot,ip,l1,idl1,chr,chi,chr,chi,chr,chi,
     &               cr,ci,cr,ci,wa(iw))
  114    if (ncc .ne. 0) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*idot*2
  116 continue
      if (na .eq. 0) return
      do 1171 i=1,n*nrxy
      cr(i) = chr(i)
 1171 ci(i) = chi(i)
      return
      end
      subroutine cfft3f(ng,nrxy,n,c,work,wsave,ifac)
      implicit real*8(a-h,o-z)
      dimension c(ng,2),wsave(*),work(ng,2),ifac(*)
      if (n .eq. 1) return
      call cft3f1 (nrxy,n,c(1,1),c(1,2),work(1,1),work(1,2),wsave,ifac)
      return
      end
      subroutine cft3f1 (nrxy,n,cr,ci,chr,chi,wa,ifac)
      implicit real*8(a-h,o-z)
      dimension chr(*),chi(*),cr(*),ci(*),wa(*),ifac(*)
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido
         idl1 = idot*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+idot*2
         ix3 = ix2+idot*2
         if (na .ne. 0) go to 101
         call thdf4v(nrxy,idot,l1,cr,ci,chr,chi,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call thdf4v(nrxy,idot,l1,chr,chi,cr,ci,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call thdf2v(nrxy,idot,l1,cr,ci,chr,chi,wa(iw))
         go to 105
  104    call thdf2v(nrxy,idot,l1,chr,chi,cr,ci,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip .ne. 3) go to 109
         ix2 = iw+idot*2
         if (na .ne. 0) go to 107
         call thdf3v(nrxy,idot,l1,cr,ci,chr,chi,wa(iw),wa(ix2))
         go to 108
  107    call thdf3v(nrxy,idot,l1,chr,chi,cr,ci,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip .ne. 5) go to 112
         ix2 = iw+idot*2
         ix3 = ix2+idot*2
         ix4 = ix3+idot*2
         if (na .ne. 0) go to 110
         call thdf5v(nrxy,idot,l1,cr,ci,chr,chi,
     &               wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call thdf5v(nrxy,idot,l1,chr,chi,cr,ci,
     &               wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na .ne. 0) go to 113
         call thdfg (nrxy,ncc,idot,ip,l1,idl1,cr,ci,cr,ci,cr,ci,
     &               chr,chi,chr,chi,wa(iw))
         go to 114
  113    call thdfg (nrxy,ncc,idot,ip,l1,idl1,chr,chi,chr,chi,chr,chi,
     &               cr,ci,cr,ci,wa(iw))
  114    if (ncc .ne. 0) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*idot*2
  116 continue
      if (na .eq. 0) return
      do 1171 i=1,n*nrxy
      cr(i) = chr(i)
 1171 ci(i) = chi(i)
      return
      end
C-----------------------------------------------------------------------
      subroutine thdb2v(nrxy,ido,l1,ccr,cci,chr,chi,wa1)
      implicit real*8(a-h,o-z)
      dimension ccr(nrxy,ido,2,l1),chr(nrxy,ido,l1,2),wa1(*)
      dimension cci(nrxy,ido,2,l1),chi(nrxy,ido,l1,2)
      do 203 i=1,ido
      do 203 k=1,l1
C THE FOLLOWING DUMMY STATEMENT IS ADDED SO THAT THE LONGEST
C IJ LOOP REMAINS THE INNERMOST LOOP
      if(i.eq.0.or.k.eq.0) stop
      do 203 ij=1,nrxy
          chr(ij,i,k,1) = ccr(ij,i,1,k)+ccr(ij,i,2,k)
          tr2 = ccr(ij,i,1,k)-ccr(ij,i,2,k)
          chi(ij,i,k,1) = cci(ij,i,1,k)+cci(ij,i,2,k)
          ti2 = cci(ij,i,1,k)-cci(ij,i,2,k)
          twr = wa1(2*i-1)*tr2
          twi = wa1(2*i-1)*ti2
          chi(ij,i,k,2) = twi + wa1(2*i)*tr2
          chr(ij,i,k,2) = twr - wa1(2*i)*ti2
  203 continue
      return
      end
      subroutine thdb3v(nrxy,ido,l1,ccr,cci,chr,chi,wa1,wa2)
      implicit real*8(a-h,o-z)
      dimension ccr(nrxy,ido,3,l1),chr(nrxy,ido,l1,3),
     &          cci(nrxy,ido,3,l1),chi(nrxy,ido,l1,3),
     &          wa1(*),wa2(*)
      data taur,taui /-.5d0,.866025403784439d0/
      do 203 i=1,ido
      do 203 k=1,l1
      do 203 ij=1,nrxy
          tr2 = ccr(ij,i,2,k)+ccr(ij,i,3,k)
          cr2 = ccr(ij,i,1,k)+taur*tr2
          chr(ij,i,k,1) = ccr(ij,i,1,k)+tr2
          ti2 = cci(ij,i,2,k)+cci(ij,i,3,k)
          ci2 = cci(ij,i,1,k)+taur*ti2
          chi(ij,i,k,1) = cci(ij,i,1,k)+ti2
          cr3 = taui*(ccr(ij,i,2,k)-ccr(ij,i,3,k))
          ci3 = taui*(cci(ij,i,2,k)-cci(ij,i,3,k))
          dr2 = cr2-ci3
          dr3 = cr2+ci3
          di2 = ci2+cr3
          di3 = ci2-cr3
          dwr2= wa1(2*i-1)*dr2
          dwi2= wa1(2*i-1)*di2
          dwr3= wa2(2*i-1)*dr3
          dwi3= wa2(2*i-1)*di3
          chi(ij,i,k,2) = dwi2+wa1(2*i)*dr2
          chr(ij,i,k,2) = dwr2-wa1(2*i)*di2
          chi(ij,i,k,3) = dwi3+wa2(2*i)*dr3
          chr(ij,i,k,3) = dwr3-wa2(2*i)*di3
  203 continue
      return
      end
      subroutine thdb4v(nrxy,ido,l1,ccr,cci,chr,chi,wa1,wa2,wa3)
      implicit real*8(a-h,o-z)
      dimension ccr(nrxy,ido,4,l1),chr(nrxy,ido,l1,4),
     &          cci(nrxy,ido,4,l1),chi(nrxy,ido,l1,4),
     &                wa1(*)     ,wa2(*)     ,wa3(*)
      do 203 i=1,ido
      do 203 k=1,l1
      do 203 ij=1,nrxy
          ti1 = cci(ij,i,1,k)-cci(ij,i,3,k)
          ti2 = cci(ij,i,1,k)+cci(ij,i,3,k)
          ti3 = cci(ij,i,2,k)+cci(ij,i,4,k)
          tr4 = cci(ij,i,4,k)-cci(ij,i,2,k)
          tr1 = ccr(ij,i,1,k)-ccr(ij,i,3,k)
          tr2 = ccr(ij,i,1,k)+ccr(ij,i,3,k)
          ti4 = ccr(ij,i,2,k)-ccr(ij,i,4,k)
          tr3 = ccr(ij,i,2,k)+ccr(ij,i,4,k)
          chr(ij,i,k,1) = tr2+tr3
          cr3 = tr2-tr3
          chi(ij,i,k,1) = ti2+ti3
          ci3 = ti2-ti3
          cr2 = tr1+tr4
          cr4 = tr1-tr4
          ci2 = ti1+ti4
          ci4 = ti1-ti4
          wr2 = wa1(2*i-1)*cr2
          wi2 = wa1(2*i-1)*ci2
          wr3 = wa2(2*i-1)*cr3
          wi3 = wa2(2*i-1)*ci3
          wr4 = wa3(2*i-1)*cr4
          wi4 = wa3(2*i-1)*ci4
          chr(ij,i,k,2) = wr2-wa1(2*i)*ci2
          chi(ij,i,k,2) = wi2+wa1(2*i)*cr2
          chr(ij,i,k,3) = wr3-wa2(2*i)*ci3
          chi(ij,i,k,3) = wi3+wa2(2*i)*cr3
          chr(ij,i,k,4) = wr4-wa3(2*i)*ci4
          chi(ij,i,k,4) = wi4+wa3(2*i)*cr4
  203 continue
      return
      end
      subroutine thdb5v(nrxy,ido,l1,ccr,cci,chr,chi,wa1,wa2,wa3,wa4)
      implicit real*8(a-h,o-z)
      dimension  ccr(nrxy,ido,5,l1),chr(nrxy,ido,l1,5),
     1           cci(nrxy,ido,5,l1),chi(nrxy,ido,l1,5),
     1           wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
      data tr11,ti11,tr12,ti12 /.309016994374947d0,.951056516295154d0,
     1-.809016994374947d0,.587785252292473d0/
      do 203 i=1,ido
      do 203 k=1,l1
      do 203 ij=1,nrxy
          ti5 = cci(ij,i,2,k)-cci(ij,i,5,k)
          ti2 = cci(ij,i,2,k)+cci(ij,i,5,k)
          ti4 = cci(ij,i,3,k)-cci(ij,i,4,k)
          ti3 = cci(ij,i,3,k)+cci(ij,i,4,k)
          tr5 = ccr(ij,i,2,k)-ccr(ij,i,5,k)
          tr2 = ccr(ij,i,2,k)+ccr(ij,i,5,k)
          tr4 = ccr(ij,i,3,k)-ccr(ij,i,4,k)
          tr3 = ccr(ij,i,3,k)+ccr(ij,i,4,k)
          chr(ij,i,k,1) = ccr(ij,i,1,k)+tr2+tr3
          chi(ij,i,k,1) = cci(ij,i,1,k)+ti2+ti3
          cr2 = ccr(ij,i,1,k)+tr11*tr2+tr12*tr3
          ci2 = cci(ij,i,1,k)+tr11*ti2+tr12*ti3
          cr3 = ccr(ij,i,1,k)+tr12*tr2+tr11*tr3
          ci3 = cci(ij,i,1,k)+tr12*ti2+tr11*ti3
          cr5 = ti11*tr5+ti12*tr4
          ci5 = ti11*ti5+ti12*ti4
          cr4 = ti12*tr5-ti11*tr4
          ci4 = ti12*ti5-ti11*ti4
          dr3 = cr3-ci4
          dr4 = cr3+ci4
          di3 = ci3+cr4
          di4 = ci3-cr4
          dr5 = cr2+ci5
          dr2 = cr2-ci5
          di5 = ci2-cr5
          di2 = ci2+cr5
          chr(ij,i,k,2) = wa1(2*i-1)*dr2-wa1(2*i)*di2
          chi(ij,i,k,2) = wa1(2*i-1)*di2+wa1(2*i)*dr2
          chr(ij,i,k,3) = wa2(2*i-1)*dr3-wa2(2*i)*di3
          chi(ij,i,k,3) = wa2(2*i-1)*di3+wa2(2*i)*dr3
          chr(ij,i,k,4) = wa3(2*i-1)*dr4-wa3(2*i)*di4
          chi(ij,i,k,4) = wa3(2*i-1)*di4+wa3(2*i)*dr4
          chr(ij,i,k,5) = wa4(2*i-1)*dr5-wa4(2*i)*di5
          chi(ij,i,k,5) = wa4(2*i-1)*di5+wa4(2*i)*dr5
  203 continue
      return
      end
      subroutine thdf2v(nrxy,ido,l1,ccr,cci,chr,chi,wa1)
      implicit real*8(a-h,o-z)
      dimension ccr(nrxy,ido,2,l1),chr(nrxy,ido,l1,2),
     &          cci(nrxy,ido,2,l1),chi(nrxy,ido,l1,2),
     1          wa1(*)
      do 203 i=1,ido
      do 203 k=1,l1
C THE FOLLOWING DUMMY STATEMENT IS ADDED SO THAT THE LONGEST
C IJ LOOP REMAINS THE INNERMOST LOOP
      if(i.eq.0.or.k.eq.0) stop
      do 203 ij=1,nrxy
          chr(ij,i,k,1) = ccr(ij,i,1,k)+ccr(ij,i,2,k)
          tr2 = ccr(ij,i,1,k)-ccr(ij,i,2,k)
          chi(ij,i,k,1) = cci(ij,i,1,k)+cci(ij,i,2,k)
          ti2 = cci(ij,i,1,k)-cci(ij,i,2,k)
          chi(ij,i,k,2) = wa1(2*i-1)*ti2-wa1(2*i)*tr2
          chr(ij,i,k,2) = wa1(2*i-1)*tr2+wa1(2*i)*ti2
  203 continue
      return
      end
      subroutine thdf3v(nrxy,ido,l1,ccr,cci,chr,chi,wa1,wa2)
      implicit real*8(a-h,o-z)
      dimension ccr(nrxy,ido,3,l1),chr(nrxy,ido,l1,3),
     &          cci(nrxy,ido,3,l1),chi(nrxy,ido,l1,3),
     1                wa1(*)     ,wa2(*)
      data taur,taui /-.5d0,-.866025403784439d0/
      do 203 i=1,ido
      do 203 k=1,l1
      do 203 ij=1,nrxy
          tr2 = ccr(ij,i,2,k)+ccr(ij,i,3,k)
          cr2 = ccr(ij,i,1,k)+taur*tr2
          chr(ij,i,k,1) = ccr(ij,i,1,k)+tr2
          ti2 = cci(ij,i,2,k)+cci(ij,i,3,k)
          ci2 = cci(ij,i,1,k)+taur*ti2
          chi(ij,i,k,1) = cci(ij,i,1,k)+ti2
          cr3 = taui*(ccr(ij,i,2,k)-ccr(ij,i,3,k))
          ci3 = taui*(cci(ij,i,2,k)-cci(ij,i,3,k))
          dr2 = cr2-ci3
          dr3 = cr2+ci3
          di2 = ci2+cr3
          di3 = ci2-cr3
          chi(ij,i,k,2) = wa1(2*i-1)*di2-wa1(2*i)*dr2
          chr(ij,i,k,2) = wa1(2*i-1)*dr2+wa1(2*i)*di2
          chi(ij,i,k,3) = wa2(2*i-1)*di3-wa2(2*i)*dr3
          chr(ij,i,k,3) = wa2(2*i-1)*dr3+wa2(2*i)*di3
  203 continue
      return
      end
      subroutine thdf4v(nrxy,ido,l1,ccr,cci,chr,chi,wa1,wa2,wa3)
      implicit real*8(a-h,o-z)
      dimension ccr(nrxy,ido,4,l1),chr(nrxy,ido,l1,4),
     &          cci(nrxy,ido,4,l1),chi(nrxy,ido,l1,4),
     1          wa1(*)     ,wa2(*)     ,wa3(*)
      do 203 i=1,ido
      do 203 k=1,l1
      do 203 ij=1,nrxy
          ti1 = cci(ij,i,1,k)-cci(ij,i,3,k)
          ti2 = cci(ij,i,1,k)+cci(ij,i,3,k)
          ti3 = cci(ij,i,2,k)+cci(ij,i,4,k)
          tr4 = cci(ij,i,2,k)-cci(ij,i,4,k)
          tr1 = ccr(ij,i,1,k)-ccr(ij,i,3,k)
          tr2 = ccr(ij,i,1,k)+ccr(ij,i,3,k)
          ti4 = ccr(ij,i,4,k)-ccr(ij,i,2,k)
          tr3 = ccr(ij,i,2,k)+ccr(ij,i,4,k)
          chr(ij,i,k,1) = tr2+tr3
          cr3 = tr2-tr3
          chi(ij,i,k,1) = ti2+ti3
          ci3 = ti2-ti3
          cr2 = tr1+tr4
          cr4 = tr1-tr4
          ci2 = ti1+ti4
          ci4 = ti1-ti4
          chr(ij,i,k,2) = wa1(2*i-1)*cr2+wa1(2*i)*ci2
          chi(ij,i,k,2) = wa1(2*i-1)*ci2-wa1(2*i)*cr2
          chr(ij,i,k,3) = wa2(2*i-1)*cr3+wa2(2*i)*ci3
          chi(ij,i,k,3) = wa2(2*i-1)*ci3-wa2(2*i)*cr3
          chr(ij,i,k,4) = wa3(2*i-1)*cr4+wa3(2*i)*ci4
          chi(ij,i,k,4) = wa3(2*i-1)*ci4-wa3(2*i)*cr4
  203 continue
      return
      end
      subroutine thdf5v(nrxy,ido,l1,ccr,cci,chr,chi,wa1,wa2,wa3,wa4)
      implicit real*8(a-h,o-z)
      dimension ccr(nrxy,ido,5,l1),chr(nrxy,ido,l1,5),
     &          cci(nrxy,ido,5,l1),chi(nrxy,ido,l1,5),
     1          wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
      data tr11,ti11,tr12,ti12 /.309016994374947d0,-.951056516295154d0,
     1-.809016994374947d0,-.587785252292473d0/
      do 203 i=1,ido
      do 203 k=1,l1
      do 203 ij=1,nrxy
          ti5 = cci(ij,i,2,k)-cci(ij,i,5,k)
          ti2 = cci(ij,i,2,k)+cci(ij,i,5,k)
          ti4 = cci(ij,i,3,k)-cci(ij,i,4,k)
          ti3 = cci(ij,i,3,k)+cci(ij,i,4,k)
          tr5 = ccr(ij,i,2,k)-ccr(ij,i,5,k)
          tr2 = ccr(ij,i,2,k)+ccr(ij,i,5,k)
          tr4 = ccr(ij,i,3,k)-ccr(ij,i,4,k)
          tr3 = ccr(ij,i,3,k)+ccr(ij,i,4,k)
          chr(ij,i,k,1) = ccr(ij,i,1,k)+tr2+tr3
          chi(ij,i,k,1) = cci(ij,i,1,k)+ti2+ti3
          cr2 = ccr(ij,i,1,k)+tr11*tr2+tr12*tr3
          ci2 = cci(ij,i,1,k)+tr11*ti2+tr12*ti3
          cr3 = ccr(ij,i,1,k)+tr12*tr2+tr11*tr3
          ci3 = cci(ij,i,1,k)+tr12*ti2+tr11*ti3
          cr5 = ti11*tr5+ti12*tr4
          ci5 = ti11*ti5+ti12*ti4
          cr4 = ti12*tr5-ti11*tr4
          ci4 = ti12*ti5-ti11*ti4
          dr3 = cr3-ci4
          dr4 = cr3+ci4
          di3 = ci3+cr4
          di4 = ci3-cr4
          dr5 = cr2+ci5
          dr2 = cr2-ci5
          di5 = ci2-cr5
          di2 = ci2+cr5
          chr(ij,i,k,2) = wa1(2*i-1)*dr2+wa1(2*i)*di2
          chi(ij,i,k,2) = wa1(2*i-1)*di2-wa1(2*i)*dr2
          chr(ij,i,k,3) = wa2(2*i-1)*dr3+wa2(2*i)*di3
          chi(ij,i,k,3) = wa2(2*i-1)*di3-wa2(2*i)*dr3
          chr(ij,i,k,4) = wa3(2*i-1)*dr4+wa3(2*i)*di4
          chi(ij,i,k,4) = wa3(2*i-1)*di4-wa3(2*i)*dr4
          chr(ij,i,k,5) = wa4(2*i-1)*dr5+wa4(2*i)*di5
          chi(ij,i,k,5) = wa4(2*i-1)*di5-wa4(2*i)*dr5
  203 continue
      return
      end
      subroutine thdbg1(nrxy,ido,ip,l1,ipph,ipp2,ccr,cci,chr,chi)
      implicit real*8(a-h,o-z)
      dimension chr(nrxy,ido,l1,ip),ccr(nrxy,ido,ip,l1)
      dimension chi(nrxy,ido,l1,ip),cci(nrxy,ido,ip,l1)
C
      do 208 j=2,ipph
          jc = ipp2-j
        do 207 i=1,ido
          do 206 k=1,l1
          do 206 ij=1,nrxy
            chr(ij,i,k,j) = ccr(ij,i,j,k)+ccr(ij,i,jc,k)
            chi(ij,i,k,j) = cci(ij,i,j,k)+cci(ij,i,jc,k)
            chr(ij,i,k,jc) = ccr(ij,i,j,k)-ccr(ij,i,jc,k)
            chi(ij,i,k,jc) = cci(ij,i,j,k)-cci(ij,i,jc,k)
  206     continue
  207   continue
  208 continue
      return
      end
      subroutine thdbg2(nrxy,ido,ip,l1,ccr,cci,chr,chi)
      implicit real*8(a-h,o-z)
      dimension chr(nrxy,ido,l1,ip),ccr(nrxy,ido,ip,l1)
      dimension chi(nrxy,ido,l1,ip),cci(nrxy,ido,ip,l1)
C
      do 518 k=1,l1
        do 517 i=1,ido
        do 517 ij=1,nrxy
          chi(ij,i,k,1) = cci(ij,i,1,k)
          chr(ij,i,k,1) = ccr(ij,i,1,k)
  517   continue
  518 continue
      return
      end
      subroutine thdbg3(nrxy,idl1,ip,ido,ipph,ipp2,c2r,c2i,ch2r,ch2i,wa)
      implicit real*8(a-h,o-z)
      dimension c2r(nrxy,idl1,ip),ch2r(nrxy,idl1,ip),wa(*)
      dimension c2i(nrxy,idl1,ip),ch2i(nrxy,idl1,ip)
C
      do 133 ik=1,idl1
          idl = 2-ido*2
        do 132 l=2,ipph
          lc = ipp2-l
          idl = idl+ido*2
          do 132 ij=1,nrxy
          c2r(ij,ik,l) = ch2r(ij,ik,1)+wa(idl-1)*ch2r(ij,ik,2)
          c2i(ij,ik,l) = ch2i(ij,ik,1)+wa(idl-1)*ch2i(ij,ik,2)
          c2r(ij,ik,lc) = wa(idl)*ch2r(ij,ik,ip)
          c2i(ij,ik,lc) = wa(idl)*ch2i(ij,ik,ip)
  132   continue
  133 continue
      return
      end
      subroutine thdbg4(nrxy,idl1,ip,ido,idp,ipph,ipp2,c2r,c2i,ch2r,ch2i
     &                 ,wa)
      implicit real*8(a-h,o-z)
      dimension c2r(nrxy,idl1,ip),ch2r(nrxy,idl1,ip),wa(*)
      dimension c2i(nrxy,idl1,ip),ch2i(nrxy,idl1,ip)
C
      do 138 l=2,ipph
          lc = ipp2-l
          idlj = 2+ido*2*(l-2)
          inc = ido*2*(l-1)
        do 137 j=3,ipph
            jc = ipp2-j
            idlj = idlj+inc
            if (idlj .gt. idp) idlj = idlj-idp
            war = wa(idlj-1)
            wai = wa(idlj)
          do 136 ik=1,idl1
          do 136 ij=1,nrxy
            c2r(ij,ik,l) = c2r(ij,ik,l)+war*ch2r(ij,ik,j)
            c2i(ij,ik,l) = c2i(ij,ik,l)+war*ch2i(ij,ik,j)
            c2r(ij,ik,lc) = c2r(ij,ik,lc)+wai*ch2r(ij,ik,jc)
            c2i(ij,ik,lc) = c2i(ij,ik,lc)+wai*ch2i(ij,ik,jc)
  136     continue
  137   continue
  138 continue
      return
      end
      subroutine thdbg5(nrxy,idl1,ip,ipph,ch2r,ch2i)
      implicit real*8(a-h,o-z)
      dimension ch2r(nrxy,idl1,ip),ch2i(nrxy,idl1,ip)
C
      do 541 j=2,ipph
        do 540 ik=1,idl1
        do 540 ij=1,nrxy
          ch2r(ij,ik,1) = ch2r(ij,ik,1)+ch2r(ij,ik,j)
          ch2i(ij,ik,1) = ch2i(ij,ik,1)+ch2i(ij,ik,j)
  540   continue
  541 continue
      return
      end
      subroutine thdbg6(nrxy,idl1,ip,ipph,ipp2,c2r,c2i,ch2r,ch2i)
      implicit real*8(a-h,o-z)
      dimension c2r(nrxy,idl1,ip),ch2r(nrxy,idl1,ip)
      dimension c2i(nrxy,idl1,ip),ch2i(nrxy,idl1,ip)
C
      do 253 ik=1,idl1
        do 252 j=2,ipph
          jc = ipp2-j
          do 252 ij=1,nrxy
          ch2r(ij,ik,j) = c2r(ij,ik,j)-c2i(ij,ik,jc)
          ch2i(ij,ik,j) = c2i(ij,ik,j)+c2r(ij,ik,jc)
          ch2r(ij,ik,jc) = c2r(ij,ik,j)+c2i(ij,ik,jc)
          ch2i(ij,ik,jc) = c2i(ij,ik,j)-c2r(ij,ik,jc)
  252   continue
  253 continue
      return
      end
      subroutine thdbg7(nrxy,idl1,ip,c2r,c2i,ch2r,ch2i)
      implicit real*8(a-h,o-z)
      dimension c2r(nrxy,idl1,ip),ch2r(nrxy,idl1,ip)
      dimension c2i(nrxy,idl1,ip),ch2i(nrxy,idl1,ip)
C
      do 255 ik=1,idl1
      do 255 ij=1,nrxy
       c2r(ij,ik,1) = ch2r(ij,ik,1)
       c2i(ij,ik,1) = ch2i(ij,ik,1)
  255 continue
      return
      end
      subroutine thdbg8(nrxy,ido,ip,l1,c1r,c1i,chr,chi)
      implicit real*8(a-h,o-z)
      dimension c1r(nrxy,ido,l1,ip),chr(nrxy,ido,l1,ip)
      dimension c1i(nrxy,ido,l1,ip),chi(nrxy,ido,l1,ip)
C
      do 357 j=2,ip
        do 356 k=1,l1
        do 356 ij=1,nrxy
          c1r(ij,1,k,j) = chr(ij,1,k,j)
          c1i(ij,1,k,j) = chi(ij,1,k,j)
  356   continue
  357 continue
      return
      end
      subroutine thdbg9(nrxy,ido,ip,l1,c1r,c1i,chr,chi,wa)
      implicit real*8(a-h,o-z)
      dimension c1r(nrxy,ido,l1,ip),chr(nrxy,ido,l1,ip),wa(*)
      dimension c1i(nrxy,ido,l1,ip),chi(nrxy,ido,l1,ip)
C
        idj = -ido*2
      do 264 j=2,ip
          idj = idj+ido*2
        do 263 i=2,ido
            idij = idj+i*2
          do 262 k=1,l1
          do 262 ij=1,nrxy
            c1r(ij,i,k,j) = wa(idij-1)*chr(ij,i,k,j)
     &                     - wa(idij  )*chi(ij,i,k,j)
            c1i(ij,i,k,j) = wa(idij-1)*chi(ij,i,k,j)
     &                     + wa(idij  )*chr(ij,i,k,j)
  262     continue
  263   continue
  264 continue
      return
      end
      subroutine thdfg1(nrxy,ido,ip,l1,ipph,ipp2,ccr,cci,chr,chi)
      implicit real*8(a-h,o-z)
      dimension chr(nrxy,ido,l1,ip),ccr(nrxy,ido,ip,l1)
      dimension chi(nrxy,ido,l1,ip),cci(nrxy,ido,ip,l1)
C
      do 108 i=1,ido
        do 107 j=2,ipph
            jc = ipp2-j
          do 106 k=1,l1
          do 106 ij=1,nrxy
            chr(ij,i,k,j) = ccr(ij,i,j,k)+ccr(ij,i,jc,k)
            chi(ij,i,k,j) = cci(ij,i,j,k)+cci(ij,i,jc,k)
            chr(ij,i,k,jc) = ccr(ij,i,j,k)-ccr(ij,i,jc,k)
            chi(ij,i,k,jc) = cci(ij,i,j,k)-cci(ij,i,jc,k)
  106     continue
  107   continue
  108 continue
      return
      end
      subroutine thdfg2(nrxy,ido,ip,l1,ccr,cci,chr,chi)
      implicit real*8(a-h,o-z)
      dimension chr(nrxy,ido,l1,ip),ccr(nrxy,ido,ip,l1)
      dimension chi(nrxy,ido,l1,ip),cci(nrxy,ido,ip,l1)
C
      do 518 k=1,l1
        do 517 i=1,ido
        do 517 ij=1,nrxy
          chr(ij,i,k,1) = ccr(ij,i,1,k)
          chi(ij,i,k,1) = cci(ij,i,1,k)
  517   continue
  518 continue
      return
      end
      subroutine thdfg3(nrxy,idl1,ip,ido,ipph,ipp2,c2r,c2i,ch2r,ch2i,wa)
      implicit real*8(a-h,o-z)
      dimension c2r(nrxy,idl1,ip),ch2r(nrxy,idl1,ip),wa(*)
      dimension c2i(nrxy,idl1,ip),ch2i(nrxy,idl1,ip)
C
      do 133 ik=1,idl1
          idl = 2-ido*2
        do 132 l=2,ipph
          lc = ipp2-l
          idl = idl+ido*2
          do 132 ij=1,nrxy
          c2r(ij,ik,l) = ch2r(ij,ik,1)+wa(idl-1)*ch2r(ij,ik,2)
          c2i(ij,ik,l) = ch2i(ij,ik,1)+wa(idl-1)*ch2i(ij,ik,2)
          c2r(ij,ik,lc) = -wa(idl)*ch2r(ij,ik,ip)
          c2i(ij,ik,lc) = -wa(idl)*ch2i(ij,ik,ip)
  132   continue
  133 continue
      return
      end
      subroutine thdfg4(nrxy,idl1,ip,ido,idp,ipph,ipp2,c2r,c2i,ch2r,ch2i
     &                 ,wa)
      implicit real*8(a-h,o-z)
      dimension c2r(nrxy,idl1,ip),ch2r(nrxy,idl1,ip),wa(*)
      dimension c2i(nrxy,idl1,ip),ch2i(nrxy,idl1,ip)
C
      do 138 l=2,ipph
          lc =ipp2-l
          idlj = 2+ido*2*(l-2)
          inc = ido*2*(l-1)
        do 137 j=3,ipph
            jc = ipp2-j
            idlj = idlj+inc
            if (idlj .gt. idp) idlj = idlj-idp
            war = wa(idlj-1)
            wai = wa(idlj)
          do 136 ik=1,idl1
          do 136 ij=1,nrxy
             c2r(ij,ik,l) = c2r(ij,ik,l)+war*ch2r(ij,ik,j)
             c2i(ij,ik,l) = c2i(ij,ik,l)+war*ch2i(ij,ik,j)
             c2r(ij,ik,lc) = c2r(ij,ik,lc)-wai*ch2r(ij,ik,jc)
             c2i(ij,ik,lc) = c2i(ij,ik,lc)-wai*ch2i(ij,ik,jc)
  136     continue
  137   continue
  138 continue
      return
      end
      subroutine thdfg5(nrxy,idl1,ip,ipph,ch2r,ch2i)
      implicit real*8(a-h,o-z)
      dimension ch2r(nrxy,idl1,ip)
      dimension ch2i(nrxy,idl1,ip)
C
      do 541 j=2,ipph
        do 540 ik=1,idl1
        do 540 ij=1,nrxy
          ch2r(ij,ik,1) = ch2r(ij,ik,1)+ch2r(ij,ik,j)
          ch2i(ij,ik,1) = ch2i(ij,ik,1)+ch2i(ij,ik,j)
  540   continue
  541 continue
      return
      end
      subroutine thdfg6(nrxy,idl1,ip,ipph,ipp2,c2r,c2i,ch2r,ch2i)
      implicit real*8(a-h,o-z)
      dimension c2r(nrxy,idl1,ip),ch2r(nrxy,idl1,ip)
      dimension c2i(nrxy,idl1,ip),ch2i(nrxy,idl1,ip)
C
      do 253 ik=1,idl1
        do 252 j=2,ipph
          jc = ipp2-j
          do 252 ij=1,nrxy
          ch2r(ij,ik,j) = c2r(ij,ik,j)-c2i(ij,ik,jc)
          ch2i(ij,ik,j) = c2i(ij,ik,j)+c2r(ij,ik,jc)
          ch2r(ij,ik,jc) = c2r(ij,ik,j)+c2i(ij,ik,jc)
          ch2i(ij,ik,jc) = c2i(ij,ik,j)-c2r(ij,ik,jc)
  252   continue
  253 continue
      return
      end
      subroutine thdfg7(nrxy,idl1,ip,c2r,c2i,ch2r,ch2i)
      implicit real*8(a-h,o-z)
      dimension c2r(nrxy,idl1,ip),ch2r(nrxy,idl1,ip)
      dimension c2i(nrxy,idl1,ip),ch2i(nrxy,idl1,ip)
C
      do 255 ik=1,idl1
      do 255 ij=1,nrxy
       c2r(ij,ik,1) = ch2r(ij,ik,1)
       c2i(ij,ik,1) = ch2i(ij,ik,1)
  255 continue
      return
      end
      subroutine thdfg8(nrxy,ido,ip,l1,c1r,c1i,chr,chi)
      implicit real*8(a-h,o-z)
      dimension c1r(nrxy,ido,l1,ip),chr(nrxy,ido,l1,ip)
      dimension c1i(nrxy,ido,l1,ip),chi(nrxy,ido,l1,ip)
C
      do 357 j=2,ip
        do 356 k=1,l1
        do 356 ij=1,nrxy
          c1r(ij,1,k,j) = chr(ij,1,k,j)
          c1i(ij,1,k,j) = chi(ij,1,k,j)
  356   continue
  357 continue
      return
      end
      subroutine thdfg9(nrxy,ido,l1,ip,c1r,c1i,chr,chi,wa)
      implicit real*8(a-h,o-z)
      dimension chr(nrxy,ido,l1,ip),c1r(nrxy,ido,l1,ip),wa(*)
      dimension chi(nrxy,ido,l1,ip),c1i(nrxy,ido,l1,ip)
C
        idj = -ido*2
      do 264 j=2,ip
          idj = idj+ido*2
        do 263 i=2,ido
            idij = idj+i*2
          do 262 k=1,l1
          do 262 ij=1,nrxy
            c1r(ij,i,k,j) = wa(idij-1)*chr(ij,i,k,j)
     &                     + wa(idij  )*chi(ij,i,k,j)
            c1i(ij,i,k,j) = wa(idij-1)*chi(ij,i,k,j)
     &                     - wa(idij  )*chr(ij,i,k,j)
  262     continue
  263   continue
  264 continue
      return
      end
      subroutine thdbg(nrxy,ncc,ido,ip,l1,idl1,
     & ccr,cci,c1r,c1i,c2r,c2i,chr,chi,ch2r,ch2i,wa)
      implicit real*8(a-h,o-z)
      dimension chr(nrxy,ido,l1,ip),ccr(nrxy,ido,ip,l1),
     1          c1r(nrxy,ido,l1,ip),wa(*),c2r(nrxy,idl1,ip),
     2         ch2r(nrxy,idl1,ip)
      dimension chi(nrxy,ido,l1,ip),cci(nrxy,ido,ip,l1),
     1          c1i(nrxy,ido,l1,ip),      c2i(nrxy,idl1,ip),
     2         ch2i(nrxy,idl1,ip)
C
      ipp2 = ip+2
      ipph = (ip+1)/2
      idp = ip*ido*2
C
        call thdbg1(nrxy,ido,ip,l1,ipph,ipp2,ccr,cci,chr,chi)
        call thdbg2(nrxy,ido,ip,l1,ccr,cci,chr,chi)
        call thdbg3(nrxy,idl1,ip,ido,ipph,ipp2,c2r,c2i,ch2r,ch2i,wa)
        call thdbg4(nrxy,idl1,ip,ido,idp,ipph,ipp2,c2r,c2i,ch2r,ch2i,wa)
        call thdbg5(nrxy,idl1,ip,ipph,ch2r,ch2i)
        call thdbg6(nrxy,idl1,ip,ipph,ipp2,c2r,c2i,ch2r,ch2i)
      ncc = 1
      if (ido*2 .eq. 2) return
      ncc = 0
C
        call thdbg7(nrxy,idl1,ip,c2r,c2i,ch2r,ch2i)
        call thdbg8(nrxy,ido,ip,l1,c1r,c1i,chr,chi)
        call thdbg9(nrxy,ido,ip,l1,c1r,c1i,chr,chi,wa)
      return
      end
      subroutine thdfg(nrxy,ncc,ido,ip,l1,idl1,
     & ccr,cci,c1r,c1i,c2r,c2i,chr,chi,ch2r,ch2i,wa)
      implicit real*8(a-h,o-z)
      dimension chr(nrxy,ido,l1,ip),ccr(nrxy,ido,ip,l1),
     1          c1r(nrxy,ido,l1,ip),wa(*),c2r(nrxy,idl1,ip),
     2          ch2r(nrxy,idl1,ip)
      dimension chi(nrxy,ido,l1,ip),cci(nrxy,ido,ip,l1),
     1          c1i(nrxy,ido,l1,ip),      c2i(nrxy,idl1,ip),
     2          ch2i(nrxy,idl1,ip)
C
      ipp2 = ip+2
      ipph = (ip+1)/2
      idp = ip*ido*2
C
        call thdfg1(nrxy,ido,ip,l1,ipph,ipp2,ccr,cci,chr,chi)
        call thdfg2(nrxy,ido,ip,l1,ccr,cci,chr,chi)
        call thdfg3(nrxy,idl1,ip,ido,ipph,ipp2,c2r,c2i,ch2r,ch2i,wa)
        call thdfg4(nrxy,idl1,ip,ido,idp,ipph,ipp2,c2r,c2i,ch2r,ch2i,wa)
        call thdfg5(nrxy,idl1,ip,ipph,ch2r,ch2i)
        call thdfg6(nrxy,idl1,ip,ipph,ipp2,c2r,c2i,ch2r,ch2i)
      ncc = 1
      if (ido*2 .eq. 2) return
      ncc = 0
C
        call thdfg7(nrxy,idl1,ip,c2r,c2i,ch2r,ch2i)
        call thdfg8(nrxy,ido,ip,l1,c1r,c1i,chr,chi)
        call thdfg9(nrxy,ido,l1,ip,c1r,c1i,chr,chi,wa)
      return
      end
