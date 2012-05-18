      subroutine write_i2mcfm_info(unit,lstring,rstring,tag)
      implicit none 
      include 'constants.f'
      include 'interface_settings.f' 
      character*72 f96,f97,f98,f99
      character*(*) tag,lstring,rstring 
      integer unit 
      logical writeall


      
c---  f96 character format            
      f96='('''//lstring//''',a20,12x,''['',a,'']'','''//rstring//''')' 
c---  f97 integer format      
      f97='('''//lstring//''',i20,12x,''['',a,'']'','''//rstring//''')' 
c---  f98 logical format      
      f98='('''//lstring//''',L20,12x,''['',a,'']'','''//rstring//''')' 
c---  f99 floating point format
      f99='('''//lstring//''',f20.4,12x,''['',a,'']'','''
     &     //rstring//''')' 
    
      writeall=.false.
      if (tag .eq. 'WRITEALL') writeall=.true.
      

      if((tag.eq.'ret eps poles').or.(tag.eq.'WRITEALL')) then 
         write(unit,fmt=f98) ret_poles, 'ret eps poles'
      endif

      if ((tag.eq.'efficient').or.(tag.eq.'WRITEALL')) then 
         write(unit,fmt=f98) efficient,'efficient' 
      endif
      if ((tag.eq.'params').or.(tag.eq.'WRITEALL')) then 
         write(unit,fmt=f98) params,'params' 
      endif
       if ((tag.eq.'read_as').or.(tag.eq.'WRITEALL')) then 
         write(unit,fmt=f98) read_as,'read_as' 
      endif
      
       if ((tag.eq.'read_mv').or.(tag.eq.'WRITEALL')) then 
         write(unit,fmt=f98) read_mv,'read_mv' 
      endif
      
       if ((tag.eq.'read_ew').or.(tag.eq.'WRITEALL')) then 
         write(unit,fmt=f98) read_ew,'read_ew' 
      endif
      
      if ((tag.eq.'read_mq').or.(tag.eq.'WRITEALL')) then 
         write(unit,fmt=f98) read_mq,'read_mq' 
      endif
      
       if ((tag.eq.'read_ml').or.(tag.eq.'WRITEALL')) then 
         write(unit,fmt=f98) read_ml,'read_ml' 
      endif
      
       if ((tag.eq.'read_ckm').or.(tag.eq.'WRITEALL')) then 
         write(unit,fmt=f98) read_ckm,'read_ckm' 
      endif
      if ((tag.eq.'read_scheme').or.(tag.eq.'WRITEALL')) then 
         write(unit,fmt=f98) read_scheme,'read_scheme' 
      endif
      

       if ((tag.eq.'alphaew_eq1').or.(tag.eq.'WRITEALL')) then 
         write(unit,fmt=f98) alphaew_eq1,'alpha_ew = 1' 
      endif
       if ((tag.eq.'alphas_eq1').or.(tag.eq.'WRITEALL')) then 
         write(unit,fmt=f98) alphas_eq1,'alpha_s = 1' 
      endif

      if((tag.eq.'inherit alphas').or.(tag.eq.'WRITEALL')) then 
         write(unit,fmt=f98) inherit_as, 'Inherit alpha_S'
      endif

      return 
      end subroutine 
