

      subroutine fixx1x2(pin,pfull) 
      implicit none 
      include 'constants.f'
      include 'interface_settings.f' 
      double precision pin(n_final,4),pfull(mxpart,4),xx(2),yrap,pt
      double precision y_ar(n_final),pt_ar(n_final)
      double precision sqrts,rtson2 
      double precision x1,x2
      integer i,nu
      common/energy/sqrts 
      common/x1x2/x1,x2

      do i=1,mxpart
         do nu=1,4
            pfull(i,nu)=0d0 
         enddo
      enddo

      
     
      rtson2=0.5d0*sqrts

   

!------ copy accross known components 
      do i=1,n_final 
         do nu=1,4
            pfull(i+2,nu)=pin(i,nu) 
         enddo
      enddo

    

      xx(1)=0d0
      xx(2)=0d0
!-----pts/sqrts and rapidity of leptons
      do i=1,n_final 
         pt_ar(i)=pt(i+2,pfull)/rtson2
         y_ar(i)=yrap(i+2,pfull)
         xx(1)=xx(1)+pt_ar(i)*exp(y_ar(i))        
         xx(2)=xx(2)+pt_ar(i)*exp(-y_ar(i))
      enddo
    
   
      xx(1)=xx(1)*half
      xx(2)=xx(2)*half
      
      pfull(1,4)=-0.5d0*xx(1)*sqrts
      pfull(1,1)=0d0 
      pfull(1,2)=0d0 
      pfull(1,3)=-0.5d0*xx(1)*sqrts
       
      pfull(2,4)=-0.5d0*xx(2)*sqrts
      pfull(2,1)=0d0 
      pfull(2,2)=0d0 
      pfull(2,3)=+0.5d0*xx(2)*sqrts
      
      x1=xx(1) 
      x2=xx(2) 

!      write(6,*) x1,x2 
!      pause

!      write(6,*) pfull(1,1),pfull(1,2),pfull(1,3),pfull(1,4)      
!      write(6,*) pfull(2,1),pfull(2,2),pfull(2,3),pfull(2,4)
!      write(6,*) pfull(3,1),pfull(3,2),pfull(3,3),pfull(3,4)      
!      write(6,*) pfull(4,1),pfull(4,2),pfull(4,3),pfull(4,4)
!      write(6,*) pfull(5,1),pfull(5,2),pfull(5,3),pfull(5,4)      
!      write(6,*) pfull(6,1),pfull(6,2),pfull(6,3),pfull(6,4)
!      pause
   

      return 
      end
