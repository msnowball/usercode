
      subroutine check_scheme 
      implicit none 
      include 'scheme.f'
      include 'interface_settings.f'

      if(inscheme.ne.scheme) then 
         write(6,21) 'Warning setting scheme from ',inscheme,
     &        ' to ',scheme 
         inscheme=scheme
      endif
      
 21   format(1x,a28,a4,a4,a4)
      return 
      end subroutine 
