      subroutine qqb_zz_gs(p,msq)
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  H  + Z +g(p7)
c                           |    |
c                           |    --> fermion(p3)+antifermion(p4)
c                           |
c                           ---> b(p5)+bb(p6)
      implicit none 
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'
      integer j,k,nd

      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      double precision msq17_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),
     . sub17_2(4),sub27_1(4),dummyv(-nf:nf,-nf:nf),dsubv
      external qqb_zz,donothing_gvec

      ndmax=2

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,7,2,sub17_2,dsubv,msq17_2,dummyv,
     . qqb_zz,donothing_gvec)
      call dips(2,p,2,7,1,sub27_1,dsubv,msq27_1,dummyv,
     . qqb_zz,donothing_gvec)


      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
      msq(nd,j,k)=0d0
      enddo

      if((j .ne. 0 .and. k .ne. 0 .and. j.ne.-k) 
     .  .or. (j .eq. 0 .and. k .eq. 0)) goto 19

      if     (((j .gt. 0) .and. (k .lt. 0))
     .   .or. ((j .lt. 0) .and. (k .gt. 0))) then

      msq(1,j,k)=2d0*cf*sub17_2(qq)*msq17_2(j,k)
      msq(2,j,k)=2d0*cf*sub27_1(qq)*msq27_1(j,k)

      elseif ((j .ne. 0) .and. (k .eq. 0)) then
      msq(1,j,k)=0d0
      msq(2,j,k)=2d0*tr*sub27_1(qg)*msq27_1(j,-j)

      elseif ((j .eq. 0) .and. (k .ne. 0)) then
      msq(1,j,k)=2d0*tr*sub17_2(qg)*msq17_2(-k,k)
      msq(2,j,k)=0d0

      endif

   19 continue
      enddo
      enddo

  
      return      
      end


