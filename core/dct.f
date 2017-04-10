
      subroutine map14 (y,x,iel)
C---------------------------------------------------------------
C
C     Map the elemental array X from mesh M1 to mesh M4
C
C---------------------------------------------------------------
       INCLUDE 'SIZE'
       INCLUDE 'GEOM'
       INCLUDE 'IXYZ'
       INCLUDE 'INPUT'
       INCLUDE 'COMPRESS' 
 
       REAL X(LX1,LY1,LZ1)
       REAL Y(LX4,LY4,LZ4)
      
       COMMON /CTMP0/ XA(LX4,LY1,LZ1), XB(LX4,LY4,LZ1)
          
       NYZ1 = NY1*NZ1
       NXY4 = NX4*NY4
       NYZ4 = NY4*NZ4
      
       if(IF3D) then
       CALL MXM (IXM14,NX4,X,NX4,XA,NYZ1)
       DO IZ=1,NZ1    
           CALL MXM (XA(1,1,IZ),NX4,IYTM14,NY1,XB(1,1,IZ),NY4)
       enddo 
       CALL MXM (XB,NXY4,IZTM14,NZ1,Y,NZ4)
       else
          CALL MXM (IXM14,NX4,X,NX1,XA,NYZ4)
          CALL MXM (XA,NX4,IYTM14,NY1,Y,NY4)
       endif

       RETURN
       END

c----------------------------------------------------------------
      subroutine DCT14 (y,x,iel)

       INCLUDE 'SIZE'
       INCLUDE 'GEOM'
       INCLUDE 'IXYZ'
       INCLUDE 'INPUT'
       INCLUDE 'COMPRESS' 
 
       REAL X(LX4,LY4,LZ4)
       REAL Y(LX4,LY4,LZ4)
      
       COMMON /CTMP0/ XA(LX4,LY1,LZ1), XB(LX4,LY4,LZ1)
          
       NYZ1 = NY1*NZ1
       NXY4 = NX4*NY4
       NYZ4 = NY4*NZ4
      

       if(IF3D) then
       CALL MXM (DCT1d,NX4,X,NX4,XA,NYZ1)
       DO IZ=1,NZ1    
           CALL MXM (XA(1,1,IZ),NX4,DCT1dT,NY1,XB(1,1,IZ),NY4)
       enddo 
       CALL MXM (XB,NXY4,DCT1dT,NZ1,Y,NZ4)
       else
           CALL MXM (DCT1d,NX4,X,NX1,XA,NYZ4)
           CALL MXM (XA,NX4,DCT1dT,NY1,Y,NY4)
       endif


       RETURN
       END

c-------------------------------------------------------------

       subroutine map41 (y,x,iel)
C---------------------------------------------------------------
C
C     Map the elemental array X from mesh M2 to mesh M1 (Y)
C
C---------------------------------------------------------------
       INCLUDE 'SIZE'
       INCLUDE 'GEOM'
       INCLUDE 'IXYZ'
       INCLUDE 'INPUT'
       INCLUDE 'COMPRESS'

       REAL X(LX4,LY4,LZ4)
       REAL Y(LX1,LY1,LZ1)
C
       COMMON /CTMP0/ XA(LX1,LY4,LZ4), XB(LX1,LY1,LZ4)

       NYZ4 = NY4*NZ4
       NXY1 = NX1*NY1
      
C     Use the appropriate derivative- and interpolation operator in
C     the y-direction (= radial direction if axisymmetric).
C
       IF (IF3D) THEN
          CALL MXM (IXM41,NX1,X,NX4,XA,NYZ4)
          do IZ=1,NZ4
             CALL MXM (XA(1,1,IZ),NX1,IYTM41,NY4,XB(1,1,IZ),NY1)
          enddo
          CALL MXM (XB,NXY1,IZTM41,NZ4,Y,NZ1)
       ELSE
          CALL MXM (IXM41,NX1,X,NX4,XA,NYZ4)
          CALL MXM (XA,NX1,IYTM41,NY4,Y,NY1)
       ENDIF
       RETURN
       END

c--------------------------------------------------------------
       subroutine IDCT41 (y,x,iel)
C---------------------------------------------------------------
C
C     Take IDCT of x return result y 
C---------------------------------------------------------------
       INCLUDE 'SIZE'
       INCLUDE 'GEOM'
       INCLUDE 'IXYZ'
       INCLUDE 'INPUT'
       INCLUDE 'COMPRESS'

       REAL X(LX4,LY4,LZ4)
       REAL Y(LX1,LY1,LZ1)
C
       COMMON /CTMP0/ XA(LX1,LY4,LZ4), XB(LX1,LY1,LZ4)

       NYZ4 = NY4*NZ4
       NXY1 = NX1*NY1
       
       IF (IF3D) THEN
          CALL MXM (DCT1dT,NX1,X,NX4,XA,NYZ4)
          do IZ=1,NZ4
             CALL MXM (XA(1,1,IZ),NX1,DCT1d,NY4,XB(1,1,IZ),NY1)
          enddo
          CALL MXM (XB,NXY1,DCT1d,NZ4,Y,NZ1)
       ELSE
          CALL MXM (DCT1dT,NX1,X,NX4,XA,NYZ4)
          CALL MXM (XA,NX1,DCT1d,NY4,Y,NY1)
       endif

       RETURN
       END

c----------------------------------------------------------------------
      subroutine truncate(ftrunc,fin,comp,iel)
C----------------------------------------------------------------------
c     truncate in spectral space 
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
     
      integer indd(lx4*ly4*lz4)
      real fin(lx4,ly4,lz4) 
      real f(lx4,ly4,lz4) 
      real ftrunc(lx4,ly4,lz4) 
      real error, temp1, temp2,vol_el, locthres
      integer i,nxyz4,ordind,iel
      real comp,bmlocal(lx1,ly1,lz1)
      logical tt

      nxyz4=nx4*ny4*nz4;
      call copy(f,fin,nxyz4) 
      call copy(ftrunc,fin,nxyz4)   
      call copy(Bmlocal,BM1(:,:,:,iel),nxyz4)

      call col2(f,f,nxyz4)
      !call dssum(f,lx1,ly1,lz1)
      call col2(f,bmlocal,nxyz4)
      !!call absolute(f,nxyz4)
      !! this for l2 norm only, but works for all others simlarly
      ! return f ordered ascending, and array indd of indeces
      call sort(f,indd,nxyz4)   
      vol_el=elvolm1(iel)

      !vlsum(BMlocal,nxyz4)

C      locthres=compthres*sqrt(vol_el)!/sqrt(volm)!*sqrt(real(nelgv))
      locthres=compthres*sqrt(vol_el)!/sqrt(volm)!*sqrt(real(nelgv))
      error=0.0
      i=0.0
      temp1=0.0
      tt=.true.
      do while (tt)
         temp1=temp1+f(i+1,1,1)
         error=sqrt(temp1)!*vol_el!/volm
         if  ((error.lt.locthres) .and. (i+2.le.lx4*ly4*lz4)) then
             i=i+1
             ordind=indd(i) 
             ftrunc(ordind,1,1)=0.0 
         else 
             tt=.false.
         endif
              
      end do   
      !write(*,*) 'compthres', i
      comp=real(i)/real(nxyz4)
      
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine truncatem4(ftrunc,fin,comp,iel)
C----------------------------------------------------------------------
c     truncate in spectral space 
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
     
      integer indd(lx4*ly4*lz4)
      real fin(lx4,ly4,lz4) 
      real f(lx4,ly4,lz4) 
      real ftrunc(lx4,ly4,lz4) 
      real error, temp1, temp2,vol_el, locthres
      integer i,nxyz4,ordind,iel
      real comp, bmlocal(lx4,ly4,lz4)
      logical tt

      nxyz4=nx4*ny4*nz4;
      call copy(ftrunc,fin,nxyz4)   
      call copy(bmlocal,BM4(:,:,:,iel),nxyz4)

      call col3(f,fin,fin,nxyz4)
      call col2(f,bmlocal,nxyz4)
      !!call absolute(f,nxyz4)
      !! this for l2 norm only, but works for all others simlarly
      ! return f ordered ascending, and array indd of indeces
      call sort(f,indd,nxyz4)   
      vol_el=elvolm4(iel)

      locthres=compthres*sqrt(vol_el)
      error=0.0
      i=0.0
      temp1=0.0
      tt=.true.
      do while (tt)
         temp1=temp1+f(i+1,1,1)
         error=sqrt(temp1)
         if  ((error.lt.locthres) .and. (i+2.le.lx4*ly4*lz4)) then
             i=i+1
             ordind=indd(i) 
             ftrunc(ordind,1,1)=0.0 
          else 
             tt=.false.
         endif
              
      end do   
      comp=real(i)/real(nxyz4)
      
      RETURN
      END

c-----------------------------------------------------------------------
      subroutine truncatedlt(ftrunc,fin,comp,iel)
C----------------------------------------------------------------------
c     truncate in spectral space 
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
     
      integer indd(lx1*ly1*lz1)
      real fin(lx1,ly1,lz1) 
      real f(lx1,ly1,lz1) 
      real ftrunc(lx1,ly1,lz1) 
      real error, temp1, temp2,vol_el, locthres,fact
      integer i,nxyz1,ordind,iel
      real comp, bmlocal(lx1,ly1,lz1)
      logical tt

      nxyz1=nx1*ny1*nz1
      call copy(ftrunc,fin,nxyz1)   
      call copy(bmlocal,BM1(:,:,:,iel),nxyz1)

      call rone(indd,nxyz1)
      call col3(f,fin,fin,nxyz1)
      call col2(f,bmlocal,nxyz1)
      !! this for l2 norm only, but works for all others simlarly
      ! return f ordered ascending, and array indd of indeces
      call sort(f,indd,nxyz1)   
      vol_el=elvolm1(iel)
      !write(*,*) 'indd', (f(i,1,1), i=1,nx1*nx1)
      locthres=compthres*sqrt(vol_el)
      error=0.0
      i=0
      temp1=0.0
      tt=.true.
      do while (tt.eqv..true.)!.and.i.lt.nxyz1)
         
         fact=fDLT(indd(i+1),1,1)/bmlocal(indd(i+1),1,1) 
         temp1=temp1+dble(f(i+1,1,1)*fact)
         !write(*,*) 'temp',i+1,indd(i+1),f(i+1,1,1),fact,temp1
         error=sqrt(temp1)
         i=i+1
         if  ((error.lt.locthres) .and. (i+2.le.lx1*ly1*lz1)) then
             i=i+1
             ordind=indd(i) 
             ftrunc(ordind,1,1)=0.0 
          else 
             tt=.false.
         endif   
      end do  
      !call exitt
      !write(*,*), 'comp',i 
      comp=dble(dble(i)/dble(nxyz1))
      
      RETURN
      END
c----------------------------------------------------------------------
      SUBROUTINE DChebT
C----------------------------------------------------------------------
C
C     Compute the one-dimensional DCT matrix
c     extra care needed since this is not done flexibly accors dimensions
C--------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
     
      real fact(lx4)
      integer i,j

      call rone(fact,nx4)
      call cmult(fact,sqrt(2.0/nx4),nx4)
      fact(1)=1/sqrt(real(nx4))

       do j=1,ny4
          do i=1,nx4
          DCT1d(j,i)=cos((2.0*i-1.0)*(j-1.0)*PI/2.0/real(nx4))*fact(j)
          DCT1dT(i,j)=DCT1d(j,i)
         end do
      end do 
      RETURN
      END

c----------------------------------------------------------------------
      SUBROUTINE DLT_build
C----------------------------------------------------------------------
C
C     Compute the one-dimensional DCT matrix
c     extra care needed since this is not done flexibly accors dimensions
C--------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'COMPRESS'
     
      real Lj(lx1), fac(lx1)
      integer i,j


      call rone(fac,nx1) 

      do i=1,nx1
          call legendre_poly(Lj,ZGM1(i,1),nx1)
          do j=1,ny1   
          fac(j)=dble(dble(2.0*(j-1.0)+1.0)/2.0);
          if (j.eq.nx1) fac(j)=dble(dble(nx1-1.0)/2.0);    
          DLT1d(j,i)=Lj(j)*wxm1(i)
          iDLT1d(i,j)=dble(Lj(j)*fac(j))
          end do
      end do 
      call copy(DLT1dT, DLT1d,nx1*nx1) 
      call copy(iDLT1dT, iDLT1d,nx1*nx1) 
      call transpose1(DLT1dT,nx1)
      call transpose1(iDLT1dT,nx1)

      do i=1,nx1
         do j=1,ny1
            fDLT(i,j,1)=fac(i)*fac(j)
            if (if3d) then
            do k=1,nz1
               fDLT(i,j,k)=fac(i)*fac(j)*fac(k) 
            enddo 
            endif
         enddo
      enddo

      !write(*,*) 'DLT1d', (fac(i), i=1,nx1)
      !write(*,*) 'iDLT1dT', (iDLT1d(i,1), i=1,nx1*nx1)
      !call exitt
      RETURN
      END

c----------------------------------------------------------------
      subroutine DLT (y,x,iel)

       INCLUDE 'SIZE'
       INCLUDE 'GEOM'
       INCLUDE 'IXYZ'
       INCLUDE 'INPUT'
       INCLUDE 'COMPRESS' 
 
       REAL X(LX1,LY1,LZ1)
       REAL Y(LX1,LY1,LZ1)
      
       COMMON /CTMP0/ XA(LX1,LY1,LZ1), XB(LX1,LY1,LZ1)
          
       NYZ1 = NY1*NZ1
       NXY4 = NX4*NY4
       NYZ4 = NY4*NZ4
      

       if(IF3D) then
       CALL MXM (DLT1d,NX1,X,NX1,XA,NYZ1)
       DO IZ=1,NZ1    
           CALL MXM (XA(1,1,IZ),NX4,DLT1dT,NY1,XB(1,1,IZ),NY1)
       enddo 
       CALL MXM (XB,NXY1,DLT1dT,NZ1,Y,NZ1)
       else
           CALL MXM (DLT1d,NX1,X,NX1,XA,NYZ1)
           CALL MXM (XA,NX1,DLT1dT,NY1,Y,NY1)
       endif


       RETURN
       END

c--------------------------------------------------------------
       subroutine iDLT (y,x,iel)
C---------------------------------------------------------------
C
C     Take IDCT of x return result y 
C---------------------------------------------------------------
       INCLUDE 'SIZE'
       INCLUDE 'GEOM'
       INCLUDE 'IXYZ'
       INCLUDE 'INPUT'
       INCLUDE 'COMPRESS'

       REAL X(LX1,LY1,LZ1)
       REAL Y(LX1,LY1,LZ1)
C
       COMMON /CTMP0/ XA(LX1,LY1,LZ1), XB(LX1,LY1,LZ1)

       NXY1 = NX1*NY1
       NYZ1 = NY1*NZ1
       
       if(IF3D) then
       CALL MXM (iDLT1d,NX1,X,NX1,XA,NYZ1)
       DO IZ=1,NZ1    
           CALL MXM (XA(1,1,IZ),NX1,iDLT1dT,NY1,XB(1,1,IZ),NY1)
       enddo 
       CALL MXM (XB,NXY1,iDLT1dT,NZ1,Y,NZ1)
       else
           CALL MXM (iDLT1d,NX1,X,NX1,XA,NYZ1)
           CALL MXM (XA,NX1,iDLT1dT,NY1,Y,NY1)
       endif

       RETURN
       END

c-----------------------------------------------------------------------
      subroutine errorsM4(l2norm,maxerr,f1,f2)
c     computes l2 norm error and maxnorm on grid M4

      include 'SIZE'
      include 'TOTAL'
      INCLUDE 'COMPRESS' 

      integer n
      real l2norm, maxerr

      real f1(lx4,ly4,lz4,nelv)
      real f2(lx4,ly4,lz4,nelv)
      real error(lx4,ly4,lz4,nelv)

      n=nx4*ny4*nz4*nelv 
      call sub3(error,f1,f2,n)
      call vsq(error,n)

      l2norm= sqrt(glsc2(bm4,error,n))/volm
      maxerr= sqrt(glmax(error,n))/volm
      
      end subroutine

c-----------------------------------------------------------------------
      subroutine errorsM1(l2norm,maxerr,f1,f2)
c     computes l2 norm error and maxnorm on grid M1

      include 'SIZE'
      include 'TOTAL'
      INCLUDE 'COMPRESS' 

      integer n
      real l2norm, maxerr
      
      real f1(lx1,ly1,lz1,nelv)
      real f2(lx1,ly1,lz1,nelv)
      real error(lx1,ly1,lz1,nelv)
           
      n=nx1*ny1*nz1*nelv 
      
      call sub3(error,f1,f2,n)
      call absolute(error,n)
      maxerr= glmax(error,n)/volm   
      call vsq(error,n)
C      l2norm= sqrt(glsc2(bm1,error,n))/sqrt(volm)
      l2norm= sqrt(glsc2(bm1,error,n))/sqrt(volm)

      !write(*,*) 'volume', volm

      end subroutine

      subroutine genmeshM4
C---------------------------------------------------------------
C
C     Generate grid arrays ZGM4(1,1:3)
C     Generate Grids XM4, YM4, ZM4 and Interpolants
C     Generate Jacobians and mass matrix also, the other geometric factors  
C     emulates routine geom2 in coef.f, leaving out the axissymmteric cases 
C---------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'COMPRESS'
 
      integer n, n4,nxyz4,ii

      nxyz4=nx4*ny4*nz4
      n=nx1*ny1*nz1*nelv
      n4=nx4*ny4*nz4*nelv

      call rzero(ZGM4,lx4*ldim)
      call rone(wzm4,lz4)
c     generate chebyshev grid
      call ZWGLJD (ZGM4(1,1),WXM4,lx4,-0.5,-0.5) 
      call ZWGLJD (ZGM4(1,2),WYM4,ly4,-0.5,-0.5)  
      if(IF3D) then
          call ZWGLJD (ZGM4(1,3),WZM4,lz4,-0.5,-0.5)  
      endif
      
      DO IY=1,NY4
       DO IX=1,NX4
         W3M4(IX,IY,1)=WXM4(IX)*WYM4(IY)
       enddo
      enddo

      if (if3d) then
        DO Ix=1,Nx4
           DO Iy=1,Ny4
             do iz=1,Nz4
               W3M4(IX,IY,IZ)=WXM4(IX)*WYM4(IY)*WZM4(IZ)
             enddo
           enddo
         enddo
      endif
      
     
c     build interpolants, not for 3d yet
c     interpolate from mesh  M4(Chebyshev) to M1(GLL)
      CALL IGLJM (IXM41,IXTM41,ZGM4(1,1),ZGM1(1,1),NX4,NX1,
     $ NX4,NX1,-0.5,-0.5)
      CALL IGLJM (IYM41,IYTM41,ZGM4(1,2),ZGM1(1,2),NY4,NY1,
     $ NY4,NY1,-0.5,-0.5) 
      if (IF3D) then
      CALL IGLJM (IZM41,IZTM41,ZGM4(1,3),ZGM1(1,3),NZ4,NZ1,
     $ NZ4,NZ1,-0.5,-0.5) 
      endif
      
c     interpolate from mesh M1(GLL) to mesh M4(Chebyshev)
      CALL IGLLM (IXM14,IXTM14,ZGM1(1,1),ZGM4(1,1),NX1,NX4,
     $ NX1,NX4) 
      CALL IGLLM (IYM14,IYTM14,ZGM1(1,2),ZGM4(1,2),NY1,NY4,
     $ NY1,NY4) 
      if (IF3D) then
      CALL IGLLM (IZM14,IZTM14,ZGM1(1,3),ZGM4(1,3),NZ1,NZ4,
     $ NZ1,NZ4) 
      endif

c     build geometric factors
       DO IEL=1,NELV

C        Mapping from mesh M1 to mesh M4
C         CALL MAP14 (RXM4(1,1,1,IEL),RXM1(1,1,1,IEL),IEL)
C         CALL MAP14 (RYM4(1,1,1,IEL),RYM1(1,1,1,IEL),IEL)
C         CALL MAP14 (SXM4(1,1,1,IEL),SXM1(1,1,1,IEL),IEL)
C         CALL MAP14 (SYM4(1,1,1,IEL),SYM1(1,1,1,IEL),IEL)
C         IF (NDIM.EQ.3) THEN
C             CALL MAP14 (RZM4(1,1,1,IEL),RZM1(1,1,1,IEL),IEL)
C             CALL MAP14 (SZM4(1,1,1,IEL),SZM1(1,1,1,IEL),IEL)
C             CALL MAP14 (TXM4(1,1,1,IEL),TXM1(1,1,1,IEL),IEL)
C             CALL MAP14 (TYM4(1,1,1,IEL),TYM1(1,1,1,IEL),IEL)
C             CALL MAP14 (TZM4(1,1,1,IEL),TZM1(1,1,1,IEL),IEL)
C          ENDIF
          CALL MAP14 (JACM4(1,1,1,IEL),JACM1(1,1,1,IEL),IEL)
C
          CALL MAP14 (XM4(1,1,1,IEL),XM1(1,1,1,IEL),IEL)
          CALL MAP14 (YM4(1,1,1,IEL),YM1(1,1,1,IEL),IEL)
          CALL MAP14 (ZM4(1,1,1,IEL),ZM1(1,1,1,IEL),IEL)
C
C        Compute the mass matrix on mesh M4
          CALL COL3 (BM4(1,1,1,IEL),W3M4,JACM4(1,1,1,IEL),NXYZ4)
          elvolm1(iel)=vlsum(BM1(1,1,1,iel),nxyz4)
          !call invcol2 (BM1local(1,1,1,iel),vmult(1,1,1,iel),nxyz4)
          elvolm4(iel)=vlsum(BM4(1,1,1,iel),nxyz4)
          !write(*,*) 'volumes'
       end do
      volm=glsum(bm1,n)


      CALL INVERS2(BM4INV,BM4,n4)

      RETURN
      END
c----------------------------------------------------------------------
      subroutine absolute(a,n)
      REAL A(1)
       DO I=1,N
          A(I) = abs(A(I))
       end do    
      return
      END
c----------------------------------------------------------------------
