program conductance
      implicit none
      integer, parameter :: DP=SELECTED_REAL_KIND(25)
      !-------------------Constantes-------------------------------------
      real(kind=DP) :: hbar=1.054571628E-34, kb=1.380649E-23, pi=3.1416
      !------------------------------------------------------------------
      !------------------ Argumentos-------------------------------------
      integer :: n=10000 !Numero de frecuencias que con las que se prueba el sistema
      real(kind=DP) :: fmax=2.0d0, fmin=0.0d0
      real(kind=DP), dimension(:,:), allocatable :: tau, frec, frec2
      !------------------------------------------------------------------
      integer :: divt=1000, i, j !n=numero de frecuencias 
      real(kind=DP) :: tmin, tmax, h, fa, fb
      real(kind=DP), dimension(:,:), allocatable ::  temp, temp2, a, b, x
      real(kind=DP), dimension(:,:), allocatable :: sum1, sum2, sum3, kappa
      allocate(tau(n+1,1), frec(n+1,1), temp(divt+1,1), a(divt+1,1), b(n+1,1), x(n+1,1))
      allocate(sum1(n,1), sum2(n,1), sum3(n,1), kappa(divt+1,1), temp2(divt+1,1), frec2(n+1,1))
      open(1,file='tau-canal1-x.dat')
      open(2,file='K-canal1-x-1500k.dat')
     ! open(3,file='argumento.dat')
      do j=1,n+1
        read(1,*) frec(j,1), frec2(j,1), tau(j,1)
        enddo
      
      tmin = 1.0d0
      tmax = 1500.0d0
        

      do i=1,divt+1
        temp(i,1) = tmin + (i-1)*((tmax-tmin)/divt)
        temp2(i,1) = temp(i,1)*temp(i,1)
        
        a(i,1) = (hbar*hbar)/(2.0*pi*kb*temp2(i,1))
        

        !---------------------Parte del barrido de frecuencias---------------------------------------
        !do j=1,n+1
        !read(1,*) frec(j,1), frec2(j,1), tau(j,1)
        !enddo
        do j=1,n+1
        if(frec(j,1)==0)then
               x(j,1) = 0.0d0
       else 
               b(j,1) = (hbar*frec(j,1))/(kb*temp(i,1))
               x(j,1) = (frec2(j,1)*exp(b(j,1))*tau(j,1))/((exp(b(j,1))-1)*(exp(b(j,1))-1))
        endif
        enddo 
        !--------------------------------------------------------------------------------------------
        
        h = (fmax-fmin)/(n+1)

        do j=1,n-2,3
        sum1(j,1) = x(j,1)
        enddo

        do j=2,n-1,3
        sum2(j,1) = x(j,1)
        end do

        do j=3,n-3,3
        sum3(j,1) = x(j,1)
        enddo
        !f(a)
        fa = x(1,1)
       ! f(b)
        fb = x(n+1,1)
        kappa(i,1) = (1.0)*((a(i,1)*3.0*h)/8.0)*(fa + 3.0*sum(sum1) + 3.0*sum(sum2) + 2.0*sum(sum3) + fb)
      enddo
      !do j=1,n+1
      !write(3,*) x(j,1)
      !enddo
      do i=1,divt+1
      write(2,*) temp(i,1), kappa(i,1)
      enddo
      end program conductance
