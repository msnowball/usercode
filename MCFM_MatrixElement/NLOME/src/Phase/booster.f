      subroutine booster(lflag,q,ph,p)
      implicit none
      integer lflag,i
      double precision rsq,q2,dp,c1
      double precision q(4),ph(4),p(4)
!---                                        _                       
!---   Boost of a 4-vector ( relative speed q/q(0) ):               
!---                                                                
!---   ph is the 4-vector in the rest frame of q                    
!---   p is the corresponding 4-vector in the lab frame             
!---                                                                  
!---              INPUT       OUTPUT         
!---                                                            
!---  lflag= 0:   q, ph       p             
!---                                                                
!---  lflag= 1:   q, p        ph            
!---
      q2=q(4)*q(4)-q(1)*q(1)-q(2)*q(2)-q(3)*q(3)
      rsq=sqrt(q2)
      if (lflag==0) then
         dp = (q(1)*ph(1)+q(2)*ph(2)+q(3)*ph(3)+q(4)*ph(4))/rsq
         c1=(dp+ph(4))/(rsq+q(4))
         do i=1,3
            p(i)=ph(i)+c1*q(i)
         enddo
         p(4)=dp
      else
         dp=(q(4)*p(4)-q(1)*p(1)-q(2)*p(2)-q(3)*p(3))/rsq
         c1=(dp+p(4))/(rsq+q(4))
         do i=1,3
            ph(i)=p(i)-c1*q(i)
         enddo
         ph(4)=dp
      endif
      return
      end
