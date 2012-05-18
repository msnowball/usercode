!----------------------- Interface to MCFM ------------------------
!-                                                                -
!- This routine reads in and intializes the physical parameters   - 
!-  used by MCFM for masses, widths ew couplings etc. 
!-  
!-  Individual sections can be turned off using variables in 
!-  Interface settings if own subroutines are being used 
!-                                                                -
!-  C. Williams 2011
!------------------------------------------------------------------

      subroutine read_params 
      implicit none 
      include 'masses.f'
      include 'nflav.f' 
      include 'nores.f'
      include 'ewinput.f' 
      include 'useet.f'  
      include 'qcdcouple.f'
      include 'interface_settings.f' 
      include 'pdlabel.f' 
      include 'outputstring.f'
      include 'anomcoup.f'
      include 'calc_settings.f' 

      integer j,ih1,ih2
      character*72 string
      character*4 part
      character*90 line
      character*7 pdlabel_2
      double precision Vtd,Vts,Vtb
      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      double precision cmass,bmass,sqrts
      double precision amz
      logical spira
      logical theory_mode,total_cr,write_out 
      logical flat,cuts
      integer idum,origij,idum_cr
      common/ranno_cr/idum_cr
      common/ranno/idum 
      save /ranno/
      save /ranno_cr/

      common/origij/origij
      common/flat/flat
      common/cuts/cuts
      common/write_out/write_out
      common/techP/theory_mode,total_cr
      common/spira/spira
      common/part/part
      common/energy/sqrts
      common/pdlabel_2/pdlabel_2
      common/couple/amz
      common/density/ih1,ih2
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/qmass/cmass,bmass
      data useEt/.false./
      data mtau,mtausq/1.777d0,3.157729d0/

      string=
     &'HERE'
     &//'/NLOME_params.DAT' 

!---- Open file
      open(unit=2,file=string,status='old',err=40) 

     

      read(2,*) line 
      read(2,*) nproc_1
      
      read(2,*) line 
      read(2,*) ih1
      read(2,*) ih2 
      read(2,*) sqrts

      read(2,*) line 
      read(2,*) pdlabel
      
!------ Vegas settings 
      read(2,*) line 
      read(2,*) origij 

      idum=-abs(origij)
      idum_cr=-abs(origij+1)


!-------------------------------------

!        HIGGS MASS 
      read(2,*) line 
      read(2,*) hmass 
      read(2,*) spira

!-------------------------------------------------------------------------
!
!             ALPHA_S(M_Z) 
!
!--------------------------------------------------------------------------

      read(2,*) line 
      if(read_as) then 
         read(2,*) amz
      else
         read(2,*) line 
      endif
!-------------------------------------------------------------------------
!
!             EW PARAMETERS
!
!--------------------------------------------------------------------------

!--- First common block read in is ewinput.f' 
!--- Controlled by read_ew which is set in interface_settings.DAT
      if(read_ew) then 
!-----first line is [electroweak parameters] 
         read(2,*) line 
!-----Read ewscheme (internal parameter) 
         read(2,*) ewscheme 
!-----Read G_F 
         read(2,*) Gf_inp 
!------Read alpha_em(mz) 
         read(2,*) aemmz_inp 
!------Read sin^2(theta_W) 
         read(2,*) xw_inp
!------ Read n flavours 
         read(2,*) nflav
      else 
!----- Dont read these parameters from the file 
         do j=1,6
            read(2,*) line 
         enddo
      endif

   
!-------------------------------------------------------------------------
!
!             Masses and widths of Vector Bosons 
!
!--------------------------------------------------------------------------

!---- Second common block is masses and widths of vector bosons 
      if(read_mv) then 
!----- first line is 
         read(2,*) line 
!------ wmas 
         read(2,*) wmass_inp
         wmass=wmass_inp 
!-------w_width 
         read(2,*) wwidth 
!------ zmass 
         read(2,*) zmass_inp
         zmass=zmass_inp 
!------ z_width 
         read(2,*) zwidth 
      else
         do j=1,5 
            read(2,*) line 
         enddo
      endif


!-------------------------------------------------------------------------
!
!             Masses of Quarks
!
!--------------------------------------------------------------------------

      if(read_mq) then 
         read(2,*) line 
!----- Masses of quarks in order d,u,s,c,b,t 
         read(2,*) md
         read(2,*) mu 
         read(2,*) ms
         read(2,*) mc
         cmass=mc
         mcsq=mc**2 
         read(2,*) mb
         bmass=mb
         mbsq=mb**2
         read(2,*) mt 
        
      else
         do j=1,7
            read(2,*) line 
         enddo
      endif


!-------------------------------------------------------------------------
!
!             Masses of Leptons
!
!--------------------------------------------------------------------------

!      if(read_ml) then 
!         read(2,*) line 
!----- Masses of quarks in order e,mu,tau
!         read(2,*) mel
!         read(2,*) mmu 
!         read(2,*) mtau
!         mtausq=mtau**2
!------- Tau width 
!         read(2,*) tauwidth
!      else
!         do j=1,5
!            read(2,*) line 
!         enddo
!      endif

!-------------------------------------------------------------------------
!
!             CKM MATRIX 
!  Read in ckm 
!--------------------------------------------------------------------------

      if(read_ckm) then 
         read(2,*) line 
!----- CKM MATRIX ELEMENTS        
          read(2,*) Vud
         read(2,*) Vus 
         read(2,*) Vub
         read(2,*) Vcd
         read(2,*) Vcs
         read(2,*) Vcb
         read(2,*) Vtd 
         read(2,*) Vts 
         read(2,*) Vtb 
      else
         do j=1,10
            read(2,*) line 
         enddo
      endif
      
      
      
      read(2,*) line
      read(2,*) inscheme         
     
      read(2,*) line
      read(2,*) delg1_z    
      read(2,*) delk_z
      read(2,*) delk_g
      read(2,*) lambda_z
      read(2,*) lambda_g
      read(2,*) h1Z
      read(2,*) h1gam
      read(2,*) h2Z
      read(2,*) h2gam
      read(2,*) h3Z
      read(2,*) h3gam
      read(2,*) h4Z
      read(2,*) h4gam
      read(2,*) tevscale
     

      read(2,*) 
      read(2,*) outputstring
      read(2,*) write_out

      read(2,*) line 
      read(2,*) cuts 
      read(2,*) flat
      read(2,*) n_its 
      read(2,*) err_max

      if((Vtd.ne.0d0).or.(Vts.ne.0d0)) then 
         write(6,*) 'WARNING Vtd Vts are = 0 in MCFM' 
      endif

      
   

      close(2) 

      return 
 40   write(6,*) 'Problems openining i2mcfm_params.DAT ' 
      stop 
      return 
      end subroutine 

      
