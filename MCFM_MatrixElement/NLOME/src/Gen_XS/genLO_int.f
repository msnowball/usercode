      double precision function genLO_int(vector,wgt) 
      implicit none 
      include 'constants.f'
      include 'n_born.f' 
      include 'vegas_common_cr.f'
      include 'npart.f' 
      include 'interface_settings.f' 
      include 'facscale.f' 
      include 'limits.f' 
      include 'process.f'
      include 'scale.f' 
      include 'outputstring.f' 
      
      double precision vector(mxdim_cr),wgt 
      integer j,nu,ih1,ih2,k,i
      double precision pborn(mxpart,4)
      double precision pin(n_final,4),ptemp(n_final,4),pswt
      double precision val,LOME
      integer n_events,n_max,n_start,n_wrote,n_temp 
      common/n_events/n_events,n_start,n_wrote
      double precision wt_store,ran2
      common/wt_store/wt_store
      logical failed,passed_cuts,cuts
      logical flat
      logical bin_cr,write_out
      common/write_out/write_out
      character*72 file_name
      character*4 part,my_part
      logical dynamicscale
      double precision ptgamsq,mZsq
      double precision LO_pin
      logical passed 
      common/cuts/cuts
      common/part/part
      common/bin_cr/bin_cr
      common/pborn/pborn 
      double precision LOME_MET 
 
      integer nmom,zero
      character*200 strformat
      double precision my_simp_born
      logical first
      common/flat/flat
      data first/.true./
      save first
      integer n_it,n_call 

      genLO_int=0d0
!----- Vegas parameters  

!-------- Phase Space 
      if(case.eq.'Z_only') then 
         if(flat) then 
            call gen2flat(pin,vector,pswt,*999)
         else
            call gen2l(pin,vector,pswt,*999) 
         endif
      elseif(case.eq.'W_only') then 
         call genW(pin,vector,pswt,*999) 
         write(6,*) pin 
         pause
      elseif((case.eq.'Zgamma').or.(case.eq.'Wgamma')) then
         if (flat) then
            call gen3flat(pin,vector,pswt,*999)
         else
            call gen2l_gam(pin,vector,pswt,*999)
!            call gen_V_phots_jets(vector,1,0,pin,pswt,*999)
!            call gen2l_gam_new(pin,vector,pswt,*999)
         endif
      elseif((case.eq.'ZZlept').or.(case.eq.'HZZ_4l')) then 
         if (flat) then
            call gen4flat(pin,vector,pswt,*999)
         else
            call gen4l(pin,vector,pswt,*999)
         endif
      endif
     
      if(inc_MET) then 
 !        genLO_int=LOME_MET(pin)*pswt 
         goto 8181 
      endif
!---------- mass cuts 
      call boost_Q(pin,ptemp) 
      call recastE(ptemp) 
      call fixx1x2(ptemp,pborn)
      call masscuts(pborn,*999)  
!      call writeout(pborn) 
!      pause

      if(cuts) call docuts(0,pborn,passed) 
      if(passed.eqv..false.) goto 999

      genLO_int=LOME(pin)*pswt
 8181 continue 
      val=genLO_int*wgt

!----------- writeout 
      if(write_out.eqv..false.) goto 999
      if(bin_cr) then 
!------Metropolis 
!------Initialize variables 
         if(first) then 
            file_name='/scratch/ciaran/'//outputstring//'.DAT'          
            first=.false.  
            open(unit=333,file=file_name,status='unknown')     
         endif
               
         nmom=n_final+2
         zero=ichar('0')
         strformat='('//
     &        char(nmom+zero)//
     &        '(E20.12,1X,E20.12,1X,E20.12,1X,E20.12,1X),(E20.12))'

         write(333,strformat) pborn(3:nmom,1),pborn(3:nmom,2),
     &        pborn(3:nmom,3),pborn(3:nmom,4),val/dfloat(itmx_cr)

      endif


 999  continue 
      return 

      end 
