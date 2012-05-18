      subroutine NLO_LIK(p,mur,muf,ansNLO,errNLO) 
      implicit none 
!======= LIKELIHOOD DRIVER, will take a point p, and calculate the likelihood 
!======= n_its times, where n_its is read in from the params file 
!======= then provided that ansNLO is accurate to within errNLO 
!======= the code returns, otherwise it will have another go with increased vegas params

      include 'constants.f' 
      include 'interface_settings.f'
      include 'calc_settings.f' 
      include 'scale.f'
      include 'facscale.f' 
      
      double precision p(n_final,4),mur,muf,ansNLO,errNLO
      double precision NLOME,LOME
      double precision ansit(n_its) 
      
      integer itmxb,itmxv,itmxr
      integer ncallb,ncallv,ncallr 
      integer int_flag,int_fr
      integer idum 
      integer it,i,i_o

      common/ranno/idum 
      common/my_itsNLOME/itmxb,itmxv,itmxr
      common/my_ncsNLOME/ncallb,ncallv,ncallr 
      common/int_setup/int_flag,int_fr

      
      
!----- integration params 
      i_o=idum 
      int_flag=1
      int_fr=2

!----- set up scale
      scale=mur
      facscale=muf

!---- Vegas parameters 
      it=0
      call set_vegas_params(it)
           
     
    
      if(int_fr.eq.1) then 
         ncallr=ncallr/10
      endif
     
 111  continue 
      ansNLO=0d0 
      errNLO=0d0 
!      write(6,*) 'it',it,ncallr
      if(it.gt.n_itsMX) return 
      do i=1,n_its 
         call set_vegas_params(it)
         ansit(i)=0d0 
!         write(6,*) 'Using seed = ',idum 
!===== Calcualte Likelihoods          
         ansit(i)=NLOME(p) 
         ansNLO=ansNLO+ansit(i)
!        write(6,*) 'Found ',i,ansit(i)
!====== Change seed 
         call new_seed
     
         if(idum.eq.0) then 
            idum=i_o*2 
         endif
      enddo
!====== Calculate average and err  

      ansNLO=ansNLO/dfloat(n_its) 
      do i=1,n_its 
         errNLO=errNLO+(ansNLO-ansit(i))**2
      enddo

      errNLO=dsqrt(errNLO/dfloat(n_its))
!====== Check we were good enough 
!      write(6,*) ansNLO,errNLO

      if(ansNLO.lt.1d-30) return 
      if((dabs(ansNLO)+errNLO)/dabs(ansNLO).gt.(one+err_max)) then 
       
         it=it+1 
         if(it.lt.n_itsMX) then 
            call set_vegas_params(it) 
            goto 111 
         endif
      endif

      return 
      end 


      subroutine set_vegas_params(n) 
      implicit none 
      include 'alfacut.f' 
      integer itmxb,itmxv,itmxr
      integer ncallb,ncallv,ncallr 
      common/my_itsNLOME/itmxb,itmxv,itmxr
      common/my_ncsNLOME/ncallb,ncallv,ncallr 
      integer n,idum
      logical reset_ai 
      data reset_ai /.true./ 
      save reset_ai 

      itmxb=2
      ncallb=200
      itmxv=2
      ncallv=400
      itmxr=5
      ncallr=20000

      if(n.eq.0) then 

         return

      else
         itmxv=(itmxv+(n+1)**3) 
         ncallv=ncallv*(n+1)**2
         itmxr=(itmxr+(n+1)**3) 
         ncallr=ncallr*(n+1)**2
     
      endif

      return 
      end 

      subroutine new_seed
      implicit none 
      integer idum,iloop 
      common/ranno/idum 
      integer now(3) 
      double precision ran2
      
      iloop=0
      call itime(now) 
      if(now(3).gt.0) then 
         idum=(now(1)+now(2))*now(3)**2 
      else
         idum=(now(1)+now(2))**2
      endif

 !     write(6,*) now(1),now(2),now(3),idum
 11   continue 
      if(idum.eq.0) then 
         call itime(now)  
         if(now(3).gt.0) then 
            idum=(now(1)+now(2))*now(3)**2 
         else
            idum=(now(1)+now(2))**2
         endif
      endif

      if(idum.lt.int(5d3)) then 
         idum=idum*(now(2)+now(1)+1)
         iloop=iloop+1
!         write(6,*) 'small',idum,iloop
         iloop=iloop+1 
!         write(6,*) iloop
         goto 11
      elseif(idum.gt.int(1d8)) then
         if(now(1)+now(2).gt.0) then 
            idum=idum/(now(1)+now(2)+1)
         else
            idum=idum/1000
         endif

         iloop=iloop+1
!         write(6,*) 'big',idum,iloop
         goto 11
      endif

      idum=-abs(idum)
      return 
      end
