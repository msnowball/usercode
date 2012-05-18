
      subroutine transform_fb_if(p,q,x,ip,jp,kp)
************************************************************************
*     Author: R.K. Ellis                                               *
*     September, 1999.                                                 *
*     Given p (-p1 + -p2 --> p3 ... px .. p_(npart+2))                 *
*     produce q (-q1 + -q2 --> q3 ... qx .. q_(npart+1))               *
*     by Lorentz transformation with jp denoting the vector            *
*     which is removed (ie all components if q(jp) set to zero)        *
*     ip is the emitter, kp is the specatator                          *
*     Correct branch chosen automatically                              *
*     x is x for ii,if,fi and y for ff                                 *
************************************************************************
      implicit none 
      include 'constants.f'
      include 'npart.f'
      include 'interface_settings.f' 
      double precision p(mxpart,4),q(mxpart,4),x,omx,y,omy,
     . k(4),kt(4),ks(4),kDk,ksDks,kDp(3:mxpart),ksDp(3:mxpart)
      integer ip,kp,j,nu,jp,ipart,inp,fip

      npart=n_final+1
      do j=1,npart+2
      do nu=1,4
        q(j,nu)=0d0
      enddo
      enddo

!======= Modified for NLOME, need a different transformation for 
!======= initial-initial  

      if ((ip .le. 2) .and. (kp .le. 2)) then
c---initial-initial
         ipart=1
         omx=one-x
         do j=1,npart+2
            do nu=1,4
               if (j.eq.ip) then
                  q(ipart,nu)=x*p(ip,nu)
               elseif (j.eq.jp) then
                  goto 19
               elseif (j.eq.kp) then
                  q(ipart,nu)=p(kp,nu)+p(jp,nu)+omx*p(ip,nu)
               else
                  q(ipart,nu)=p(j,nu)
               endif
            enddo
            ipart=ipart+1
 19         continue
          enddo
        
          return
          
      elseif (((ip .le. 2) .and. (kp .gt. 2)) .or.
     .        ((ip .gt. 2) .and. (kp .le. 2))) then
c---  initial-final or final-initial
         if (ip .le. 2) then
            inp=ip
            fip=kp  
         else
            inp=kp
            fip=ip
         endif
         ipart=1
         omx=one-x
        do j=1,npart+2
           do nu=1,4
              if (j.eq.inp) then
                 q(ipart,nu)=x*p(inp,nu)
              elseif (j.eq.jp) then
                 goto 20
              elseif (j.eq.fip) then
                 q(ipart,nu)=p(jp,nu)+p(fip,nu)+omx*p(inp,nu)
              else
                   q(ipart,nu)=p(j,nu)
                endif
             enddo
             ipart=ipart+1
 20          continue
          enddo
          return
          
       elseif ((ip .gt. 2) .and. (kp .gt. 2)) then
          
! 888     continue
c---final-final
        ipart=1
        y=x
        omy=one-y
           do j=1,npart+2
                do nu=1,4
                   if (j.eq.ip) then
                   
                         q(ipart,nu)=p(jp,nu)+p(ip,nu)-y/omy*p(kp,nu)
                     
                   elseif (j.eq.jp) then
                   goto 21
                   elseif (j.eq.kp) then
                         q(ipart,nu)=p(kp,nu)/omy
                    
                   else
                   q(ipart,nu)=p(j,nu)
                   endif
                enddo
                ipart=ipart+1
 21             continue
           enddo
        return
      endif
      end
