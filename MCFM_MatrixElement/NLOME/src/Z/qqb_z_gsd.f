      subroutine qqb_z_gsd(p,msq)
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  Z^0((p3)+a(p4))+g(p5)
      implicit none 
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'
      integer j,k,nd

      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      double precision msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     . sub15_2(4),sub25_1(4),dummyv(-nf:nf,-nf:nf),dsubv
      external qqb_z,donothing_gvec
      integer iab 
      common/iab/iab
      
      ndmax=2


      if(iab.eq.1) then 
!======= a is spectator 2,5 is naturally singular 
       
         call dips_fb_fi(2,p,2,5,1,sub25_1,dsubv,msq25_1,dummyv,
     .        qqb_z,donothing_gvec)
          call dips_fb_if(1,p,1,5,2,sub15_2,dsubv,msq15_2,dummyv,
     .        qqb_z,donothing_gvec)
        
         do j=-nf,nf
            do k=-nf,nf
               
               do nd=1,ndmax
                  msq(nd,j,k)=0d0
               enddo
               
               if(((j .ne. 0) .and. (k .ne. 0) .and. (j.ne.-k)) 
     .              .or. ((j .eq. 0) .and. (k .eq. 0))) goto 19

               if     (((j .gt. 0) .and. (k .lt. 0))
     .              .or. ((j .lt. 0) .and. (k .gt. 0))) then
                 
                  msq(1,j,k)=2d0*cf*sub15_2(qq)*msq15_2(j,k)
                  msq(2,j,k)=2d0*cf*sub25_1(qq)*msq25_1(j,k)

!               elseif ((j .eq. 0) .and. (k .ne. 0)) then
!                  msq(1,j,k)=2d0*tr*sub15_2(qg)*msq15_2(-k,k)
!                  msq(2,j,k)=0d0
                  
               elseif ((j .ne. 0) .and. (k .eq. 0)) then
                  msq(1,j,k)=0d0
                  msq(2,j,k)=2d0*tr*sub25_1(qg)*msq25_1(j,-j)
                  
              
                  
               endif
               
 19            continue
            enddo
         enddo

      elseif(iab.eq.2) then 
!======= b is spectator 1,5 is natrually singular
       
         call dips_fb_fi(1,p,1,5,2,sub15_2,dsubv,msq15_2,dummyv,
     .        qqb_z,donothing_gvec)
         call dips_fb_if(2,p,2,5,1,sub25_1,dsubv,msq25_1,dummyv,
     .        qqb_z,donothing_gvec)

         do j=-nf,nf
            do k=-nf,nf
               
               do nd=1,ndmax
                  msq(nd,j,k)=0d0
               enddo
               
               if(((j .ne. 0) .and. (k .ne. 0) .and. (j.ne.-k)) 
     .              .or. ((j .eq. 0) .and. (k .eq. 0))) goto 20

               if     (((j .gt. 0) .and. (k .lt. 0))
     .              .or. ((j .lt. 0) .and. (k .gt. 0))) then
                  
                  msq(1,j,k)=2d0*cf*sub15_2(qq)*msq15_2(j,k)
                  msq(2,j,k)=2d0*cf*sub25_1(qq)*msq25_1(j,k)

               elseif ((j .ne. 0) .and. (k .eq. 0)) then
                  msq(1,j,k)=0d0
                  msq(2,j,k)=2d0*tr*sub25_1(qg)*msq25_1(j,-j)

               elseif ((j .eq. 0) .and. (k .ne. 0)) then
                  msq(1,j,k)=2d0*tr*sub15_2(qg)*msq15_2(-k,k)
                  msq(2,j,k)=0d0
                  
               endif
               
 20            continue
            enddo
         enddo
      else
         write(6,*) 'unregconized iab = ',iab 
         stop 
      endif
         
      return      
      end


