      ! 1dim Kalman Filter
      program main
      parameter (size=120)
      implicit real (a-h,o-z)
      ! F:simulation H:obs operator x:state vector y:obs data R:system noise var Q:obs noise var
      real*8  F,H,R,Q,Ft,Ht,obs_dif,predict_v,S,I,predict_x
      real*8 ans
      real*8  x(121),K(121),v(121),y(121)

      x = 0
      v = 0
      ans = 5
      F = 1
      Ft = 1
      H = 1
      Ht = 1
      R = 1.0E-3
      Q = 1.0E-5
      obs_dif = 0
      predict_x = 0
      predict_v = 0
      S = 0
      K = 0
      I = 1

! data input
      open(17, file='gauss2.dat', status='old')
      do i=1,size
         read (17,*) y(i)                                            
      end do
      close(17)
! Kalman Filter
      do i=1,size
! predict
         predict_x = F * x(i)
         predict_v = F * v(i) * Ft + Q
! update
         obs_dif = y(i) - H * predict_x
         S = predict_v + R
! compute S inverse -> 1dim's S inverse is div        
         K(i) = predict_v / S
         x(i+1) = predict_x + K(i) * obs_dif
         v(i+1) = (1-K(i)) * predict_v
         print *,v(1)
      end do

! print
      print *,"x"
      do i=1,size
         print *, x(i)
      end do
      print *,"v"
      do i=1,size
         print *, v(i)
      end do
      end program main
