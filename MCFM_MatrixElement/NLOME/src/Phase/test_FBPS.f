      progra mtester_FBPS 
      implicit none 
!---- This is the testing area for the FBPS
      integer i,j,k,n
      logical FBPS
      double precision sqrts,psweight,weight
      double precision xa,xb,Et
      double precision r(3),sum(4)
      double precision pte(4,4),ptemp(4,4),p(12,4),pa(4),pb(4),pr(4)
      common/pabr/pa,pb,pr
      common/pswt/psweight
      common/PSsupport/xa,xb,Et

!---- SETUP PARAMETERS 
      sqrts=7000d0  
      n=6
!----- filler routines 
      call i2mcfm_fill_sqrts(sqrts)

      call srand(86456)
      call CMS_Phase(pte,1)
!------ Boost away missing pT 
      call my_boostZZ(pte,ptemp) 
!----- rescale Energies 
      call recastE(ptemp) 
!----- create inital state momentum 
      call fixx1x2(ptemp,p) 
      write(6,*) 'born event'
      do i=1,n
         write(6,*) i,p(i,1),p(i,2),p(i,3),p(i,4)
      enddo
      write(6,*) '---------------------------------------------------'
      weight=0d0
      do k=1,100000
         do j=1,4
            sum(j)=p(1,j)+p(2,j)
         enddo
         do i=1,3
            r(i)=rand()
         enddo
         write(6,*) 'bremsstrahlung event',k
         if (FBPS(-p(1,4),-p(2,4),r)) then
            write(6,*) 'pa',pa(1),pa(2),pa(3),pa(4)
            write(6,*) 'pb',pb(1),pb(2),pb(3),pb(4)
            write(6,*) 'pr',pr(1),pr(2),pr(3),pr(4)
            do j=1,4
               sum(j)=sum(j)+pa(j)+pb(j)+pr(j)
            enddo
            write(6,*) 'sum momenta',sum(1),sum(2),sum(3),sum(4)
            write(6,*) 'event weight',psweight
            write(6,*) 'Et =',Et
            write(6,*) 'xa =',xa
            write(6,*) 'xb =',xb
            weight=weight+psweight
         else
            write(6,*) '********** Event failed ***********'
         endif
         write(6,*)
      enddo
      write(6,*) 'average weight: ',weight/100000d0
      stop 
      end program
