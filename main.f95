program principal
      implicit none
      integer, parameter :: DP=SELECTED_REAL_KIND(30)
     
      !-------------Argumentos que se intercambiara con ls subrutinas
      integer :: i, n=10000 !Numero de frecuencias con las que se probara el sistema
      real(kind=DP), dimension(:,:), allocatable :: mtxHrev, mtxScat, frec, frec2
     
      !-----------------------------------------------------------------
      !-------------Argumentos que solo se utilizan en este programa
      real(kind=DP) :: SiBond = 296.208008d0, SiMass=4.6637066E-26
      real(kind=DP), dimension(:,:), allocatable :: mtxHrevScat, mtxensamble, sinqa, sin2qa, cosqa
      real(kind=DP), dimension(:,:), allocatable :: num, denA, den2A, denB, tau, fragile 
      allocate(mtxHrev(2,2), mtxScat(2,2), frec(n+1,1), frec2(n+1,1))
      allocate(mtxHrevScat(2,2), mtxensamble(2,2),sinqa(n+1,1), sin2qa(n+1,1), cosqa(n+1,1))
      allocate(num(n+1,1), denA(n+1,1), den2A(n+1,1), denB(n+1,1), tau(n+1,1), fragile(n+1,1))
         
      open(1,file='tau-canal6-z.dat')
      !open(2,file='revision.dat')
      !open(3,file='sinocoseno.dat')
      do i=1,n+1   
      call      hrev(i,n,mtxHrev,frec,frec2) !reservorio caliente, region de dispersion y reservorio frio son iguales
      call      scat(i,n,mtxScat)
      
      mtxHrevScat = matmul(mtxHrev,mtxScat)
      mtxensamble = matmul(mtxHrevScat,mtxHrev)

      !-------------Relacion de dispersion-------------------------------------------------------------------
      fragile(i,1) = 4.0*SiBond-SiMass*frec2(i,1)
      
      if(fragile(i,1).lt.0.0)then
              tau(i,1) = 0.0d0
      else


      sinqa(i,1) = sqrt((SiMass*frec2(i,1)*fragile(i,1))/(4.0*(SiBond**2)))
      sin2qa(i,1) = sinqa(i,1)*sinqa(i,1)                                                                  
      cosqa(i,1) = (2.0*SiBond-SiMass*frec2(i,1))/(2.0*SiBond)                                 
      !------------------------------------------------------------------------------------------------------
      
      !-----------------Calculo del Coeficiente de Transmision---------------------------------------------
      num(i,1) = 4.0*sin2qa(i,1)
      denA(i,1) = mtxensamble(1,2) - mtxensamble(2,1) + (mtxensamble(1,1) - mtxensamble(2,2))*cosqa(i,1)
      den2A(i,1) = denA(i,1)*denA(i,1)
      denB(i,1) = (mtxensamble(1,1) + mtxensamble(2,2))*(mtxensamble(1,1) + mtxensamble(2,2))*sin2qa(i,1)
      
     ! if(num(i,1).eq.0.and.den2A(i,1).eq.0)then
     !         tau(i,1) = 0.0d0
     ! else
      
      tau(i,1) = num(i,1)/(den2A(i,1) + denB(i,1))
      
      if(tau(i,1)>=0.9999.and.tau(i,1)<=1.00001)then
                tau(i,1) = 1.0
      endif

      endif
      write(1,*) frec(i,1), frec2(i,1), tau(i,1)
      !write(2,*) num(i,1), den2A(i,1), denB(i,1)
      !write(3,*) sinqa(i,1), cosqa(i,1)
      enddo
      end program principal
        
