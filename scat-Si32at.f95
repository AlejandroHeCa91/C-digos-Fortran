subroutine scat(i,n,mtxc)
      implicit none
      integer, parameter :: DP=SELECTED_REAL_KIND(30)
      integer :: i, n, j !n = Numero de frecuencias con las que se probara el sistema
      intent(in) i, n
      real(kind=DP) :: SiMass, fmin, fmax
      real(kind=DP), dimension(:,:) :: mtxa(2,2), mtxb(2,2), mtxc(2,2), mtxHrev(2,2), frec(n+1,1), frec2(n+1,1)
      real(kind=DP), dimension(:,:) :: SiBond(7,1)
      intent(out) mtxc
      !open(1,file='matrizscat')
      SiBond(1,1) = 295.035706d0
      SiBond(2,1) = 275.166748d0
      SiBond(3,1) = 296.208008d0
      SiBond(4,1) = 274.449768d0
      SiBond(5,1) = 275.005768d0
      SiBond(6,1) = 294.765503d0
      SiBond(7,1) = 295.035706d0
      SiMass = 4.6637066E-26

      fmin = 0.0d0
      fmax = 2.0*sqrt(296.208008d0/SiMass)

      frec(i,1) = fmin + (i-1)*((fmax - fmin)/n)
      frec2(i,1) = frec(i,1)*frec(i,1)

      mtxa(1,1) = 1.0d0
      mtxa(1,2) = 0.0d0
      mtxa(2,1) = 0.0d0
      mtxa(2,2) = 1.0d0
      

      do j=1,6
      mtxb(1,1) = (SiBond(j,1) + SiBond(j+1,1) - SiMass*frec2(i,1))/SiBond(j+1,1)
      mtxb(1,2) = -SiBond(j+1,1)/SiBond(j,1)
      mtxb(2,1) = 1.0d0
      mtxb(2,2) = 0.0d0
      mtxa = matmul(mtxa,mtxb) ! matriz del bloque unidad del reservorio
      enddo
      

      mtxc(1,1) = 1.0d0
      mtxc(1,2) = 0.0d0
      mtxc(2,1) = 0.0d0
      mtxc(2,2) = 1.0d0

      do j=1,1
      mtxc = matmul(mtxc,mtxa)!! Esto es para crecer el bloque de la parte scat
      enddo
      !write(1,*) mtxc
      end subroutine scat
