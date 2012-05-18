

      subroutine inv_cuts(p,passed)
!---- Take a final state point p and check invariant cuts
      include 'user_cuts.f'
      include 'interface_settings.f' 
      include 'process.f' 
      double precision p(n_final,4) 
      logical passed 
      double precision my_inv_m,my_inv_m3
      logical first
      data first /.true. /
      save first


      passed=.true. 
      
      if(first) then 
         write(6,*) '********* Invariant mass cuts **********'
         write(6,93) '* ',mll_s,' < m(lept,lept) < ',mll_h,' *'
         write(6,93) '* ',mllgam_min,' < m(lept,lept,gamma) < '
     &        ,mllgam_max,' *'
         write(6,93) '*',mlgam_s,' < m(lept,gam) < ',mlgam_h,' *'
         write(6,*) '****************************************'
         first=.false. 
      endif

 93   format(1x,a6,f7.1,a25,f7.1,a3) 

      if(case.eq.'Z_only')then 
!----------- cuts done elsewhere (in masscuts) 
         return 
      elseif(case.eq.'Zgamma') then 
!         write(6,*) my_inv_m(1,2,p),my_inv_m(1,3,p),my_inv_m(2,3,p) 
!         write(6,*) my_inv_m3(1,2,3,p)
!         pause
         if((my_inv_m(1,2,p).lt.mll_s)
     &        .or.(my_inv_m(1,2,p).gt.mll_h)) then 
            passed=.false.
            return 
         endif
         if((my_inv_m3(1,2,3,p).lt.mllgam_min)
     &        .or.(my_inv_m3(1,2,3,p).gt.mllgam_max)) then 
            passed=.false.
            return 
         endif
         if((my_inv_m(2,3,p).lt.mlgam_s)
     &        .or.(my_inv_m(2,3,p).gt.mlgam_h)) then 
            passed=.false.
            return 
         endif
         if((my_inv_m(1,3,p).lt.mlgam_s)
     &        .or.(my_inv_m(1,3,p).gt.mlgam_h)) then 
            passed=.false.
            return 
         endif

      elseif((case.eq.'ZZlept').or.(case.eq.'HZZ_4l')) then 
         return 

      endif

      return 
      end
      
      double precision function my_inv_m(i,j,p) 
      implicit none 
      include 'interface_settings.f'
      integer i,j
      double precision p(n_final,4) 
      integer nu

      my_inv_m=0d0 
      do nu=1,3
         my_inv_m=my_inv_m-(p(i,nu)+p(j,nu))**2
      enddo

      my_inv_m=my_inv_m+(p(i,4)+p(j,4))**2
      my_inv_m=dsqrt(max(my_inv_m,0d0))

      return 
      end 

      double precision function my_inv_m3(i,j,k,p) 
      implicit none 
      include 'interface_settings.f'
      integer i,j,k
      double precision p(n_final,4) 
      integer nu

      my_inv_m3=0d0 
      do nu=1,3
         my_inv_m3=my_inv_m3-(p(i,nu)+p(j,nu)+p(k,nu))**2
      enddo

      my_inv_m3=my_inv_m3+(p(i,4)+p(j,4)+p(k,4))**2
      my_inv_m3=dsqrt(max(my_inv_m3,0d0))

      return 
      end 
