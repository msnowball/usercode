
      subroutine read_cuts
      implicit none 
      include 'user_cuts.f' 
      include 'limits.f'
      include 'frag.f' 
      include 'process.f'
     
      character*72 my_string,line

      
      
      my_string=
     &'/scratch/osghpc/snowball'
     &//'/UF_HZZ4L_Analysis/FullAnalysis/NLOME/hpc/NLOME'
     &//'/Particle_Cuts.DAT' 

!---- Open file
      open(unit=2,file=my_string,status='old',err=40) 

!---- proc will choose which cuts to read in, 
!---- proc = 1 => Z=> ll 
!---- proc = 2 => Z+gamma => ll gamma
!---- proc = 3 => ZZ => 4l 


!---- always read in jet parameters
      
      read(2,*) line 
      read(2,*) jet_pt_max 
      read(2,*) jet_veto 
      
      if((case.eq.'Z_only').or.(case.eq.'W_only')) then 
         read(2,*) line
         read(2,*) lept_pt_1 
         read(2,*) lept_pt_2
         read(2,*) eta_lept
         read(2,*) mll_s 
         read(2,*) mll_h 
         wsqmin=mll_s**2 
         wsqmax=mll_h**2
         goto 111
      else
         read(2,*) line
         read(2,*) line 
         read(2,*) line 
         read(2,*) line
         read(2,*) line
         read(2,*) line
      endif

      if(case.eq.'Zgamma') then 
         read(2,*) line 
         read(2,*) lept_pt_1
!         write(6,*) 'Read 1',lept_pt_1
         read(2,*) lept_pt_2
!         write(6,*) 'Read 2',lept_pt_2
         read(2,*) gammpt 
!         write(6,*) 'Read 3',gammpt
         read(2,*) eta_lept
         read(2,*) eta_phot
         read(2,*) mll_s 
!         write(6,*) 'Read 4',mll_s
         read(2,*) mll_h
!         write(6,*) 'Read 5',mll_h
         read(2,*) mllgam_min 
!         write(6,*) 'Read 6',mllgam_min
         read(2,*) mllgam_max
!         write(6,*) 'Read 7',mllgam_max
         read(2,*) mlgam_s
!         write(6,*) 'Read 8',mlgam_s
         read(2,*) mlgam_h
!          write(6,*) 'Read 9',mlgam_h
         read(2,*) R_lgam
!         write(6,*) 'Read 10',R_lgam
         read(2,*) frag 
!          write(6,*) 'Read 11',frag          
         read(2,*) fragset
!         write(6,*) 'Read 12',fragset
         read(2,*) cone_ang
!          write(6,*) 'Read 13',cone_ang
         read(2,*) epsilon_h 
!          write(6,*) 'Read 14',epsilon_h
         wsqmin=mll_s**2 
!         pause
         wsqmax=mll_h**2
         goto 111
      else
         read(2,*) line 
         read(2,*) line
         read(2,*) line 
         read(2,*) line
         read(2,*) line
         read(2,*) line 
         read(2,*) line
         read(2,*) line 
         read(2,*) line
         read(2,*) line
         read(2,*) line 
         read(2,*) line
         read(2,*) line
         read(2,*) line 
         read(2,*) line
         read(2,*) line
         read(2,*) line
      endif

      if((case.eq.'ZZlept').or.(case.eq.'HZZ_4l')
     &     .or.(case.eq.'ZZ4lsb')) then 
         read(2,*) line
         read(2,*) lept_pt_1
         read(2,*) lept_pt_2
         read(2,*) lept_pt_3
         read(2,*) lept_pt_4
         read(2,*) eta_lept
         read(2,*) mll_s
         read(2,*) mll_h
         read(2,*) mll2_s
         read(2,*) mll2_h 
         wsqmin=mll_s**2
         wsqmax=mll_h**2
         bbsqmin=mll2_s**2
         bbsqmax=mll2_h**2
         goto 111
      else
         read(2,*) line
         read(2,*) line
         read(2,*) line 
         read(2,*) line
         read(2,*) line
         read(2,*) line 
         read(2,*) line          
         read(2,*) line
         read(2,*) line
         read(2,*) line
      endif

 111  continue
      return 
 40   write(6,*) 'Couldnt open Cuts' 
      

      end
