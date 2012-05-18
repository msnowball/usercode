
      subroutine read_interface_settings 
      implicit none 
      include 'interface_settings.f' 
      character*200 instring
      character*90 line
      instring=
     &'HERE'
     % //'/interface_settings.DAT'

      open(unit=2,file=instring,status='old',err=30) 
      
!      write(6,*) '************  Reading Interface Settings ************' 
!      write(6,*) '*                                                   *'
!      read(2,*) line 
!      read(2,*) efficient
!      call write_i2mcfm_info(6,' *   ','     *','efficient')
!      read(2,*) ret_poles 
!      call write_i2mcfm_info(6,' *   ','     *','eps poles')

      read(2,*) line 
!      write(6,*) '*     ',line,'                 *'
      read(2,*) params
!      call write_i2mcfm_info(6,' *   ','        *','params')
      read(2,*) read_as
!      call write_i2mcfm_info(6,' *   ','       *','read_as')
      read(2,*) read_ew
!      call write_i2mcfm_info(6,' *   ','       *','read_ew')
      read(2,*) read_mv
!      call write_i2mcfm_info(6,' *   ','       *','read_mv')
      read(2,*) read_mq
!      call write_i2mcfm_info(6,' *   ','       *','read_mq')
      read(2,*) read_ml
!      call write_i2mcfm_info(6,' *   ','       *','read_ml')
      read(2,*) read_ckm
!      call write_i2mcfm_info(6,' *   ','      *','read_ckm')
      read(2,*) read_scheme
!      call write_i2mcfm_info(6,' *   ','   *','read_scheme')
      
      read(2,*) line  
      read(2,*) alphas_eq1 
!      call write_i2mcfm_info(6,' *   ','   *','alphas_eq1')
      read(2,*) alphaew_eq1
!      call write_i2mcfm_info(6,' *   ','  *','alphaew_eq1')
      read(2,*) inherit_as
!      call write_i2mcfm_info(6,'*  ','   *','inherit alphas')
      write(6,*) '*****************************************************'
      close(2) 
    
 9    format(a90)

      return 
 30   write(6,*) 'Error opening interface_settings.DAT, check manual' 
      stop 
      end subroutine

