

!---------------- NLOME FORTRAN EXAMPLE -----------------------------------------
!
!   THIS IS AN EXAMPLE ROUTINE FOR NLOME : PLEASE FOLLOWING THE FOLLOWING SETPS: 
!   
!   a) Ensure you have run the script Fortran_example.sh, which will Install and make NLOME 
!   b) copy the files NLO_ZZ_params.DAT and Cuts_ZZ.DAT to NLOME/NLOME_params.DAT and NLOME/Particle_Cuts.DAT 
!      to ensure correct cuts and settings 
!   c) Run on the Terminal through ./NLO_ZZ_EX ZZ_example_events 

!---- This example will take in 10 4 lepton events and provide the NLO weights for signal and background hypothesis 
!---- We will demonstrate two signal hypothesis, 1) m_H=200 GeV and 2) m_H=400 GeV. 

!---- We will store the events, weights and erros in an output file NLO_ZZ_EX_OUT 

!---- In addition we will calculate the Log-Liklihood difference between signal only and signal + background 

!---- Note that the code does not rescale by the cross section, i.e. the code returns V(x)+R(x) in the notation of 1204.4424 


!---- This code is in current need of documentation and is still a beta release 
!---- so if you have any questions/suggestions please get in touch (ciaran@fnal.gov) 


      program NLO_ZZ_example 
      implicit none 
!------- Need this in this version to change m_H, will be upgraded to subrotine call in released version 
      include 'NLOME/src/Inc/masses.f'

!------ in this example we know the number of events 
      integer n_events,i,j
      parameter(n_events=1) 
      
!----- array for signal weights 
      double precision P_sig_200(n_events,1:2),P_sig_400(n_events,1:2)
!------ array for backgroud weights 
      double precision P_back(n_events,1:2) 
!------ cross sections for likelihoods (calcualted using MCFM using Cuts in Cuts_ZZ.DAT) 
      double precision sig_xs_200,sig_xs_400,back_xs
      parameter(sig_xs_200=1.408d0,sig_xs_400=0.213d0,back_xs=8.937d0) 
!------ no of final state particles 
      integer n_p,errno
      parameter(n_p=4)
!------- phase space point 
      double precision p(n_p,4),mu_R,mu_F
      double precision Lik_back,Lik_200,Lik_400
      double precision xs_tot_200,xs_tot_400
      double precision m4l

!-----Initalize arrays
      do i=1,n_events 
         do j=1,2
            P_back(i,j)=0d0 
            P_sig_200(i,j)=0d0 
            P_sig_400(i,j)=0d0 
         enddo
      enddo

!------- intilize NLOME (to scale 2mz) NOTE MH=200 has been set in NLOME_params.DAT so initilzed Higgs is fine 

      

      open(unit=22,file='ZZ_example_events',status='old',err=44) 


!------ SIGNAL = 125 GeV
      mu_R=125d0 
      mu_F=125d0
      call init_NLOME(mu_R,mu_F) 
      hmass=125d0
      call chooser 

!------Read in from file (note formating here)

      do i=1,n_events 
         read(22,*,iostat=errno) p(1,1),p(1,2),p(1,3),p(1,4),
     &        p(2,1),p(2,2),p(2,3),p(2,4),
     &        p(3,1),p(3,2),p(3,3),p(3,4), 
     &        p(4,1),p(4,2),p(4,3),p(4,4),m4l

!----- calculate weight mh=200
         if(i.eq.1) write(6,*) 'Calculating For Higgs singal = 125 GeV' 
!------ first entrant in P_sig_200 is the weight, second is the MC error
         call LO_LIK(p,mu_R,mu_F,P_sig_200(i,1),P_sig_200(i,2)) 
         write(6,*) 'Have done ',i,' events so far'
       enddo
         
      close(22) 
      open(unit=22,file='ZZ_example_events',status='old',err=44) 
      
      mu_R=m4l 
      mu_F=m4l 
!------ SIGNAL = 400 GeV
      call init_NLOME(mu_R,mu_F) 
!======= CHANGE HIGGS MASS (WILL BE UPGRADED TO BE EASIER IN RELEASED CODE) 
         hmass=m4l
         call chooser

!------Read in from file (note formating here)
      do i=1,n_events 
         read(22,*,iostat=errno) p(1,1),p(1,2),p(1,3),p(1,4),
     &        p(2,1),p(2,2),p(2,3),p(2,4),
     &        p(3,1),p(3,2),p(3,3),p(3,4), 
     &        p(4,1),p(4,2),p(4,3),p(4,4),m4l
         if(i.eq.1) write(6,*) 'Calculating For Higgs singal = ',m4l
          
!------ first entrant in P_sig_400 is the weight, second is the MC error
         call LO_LIK(p,mu_R,mu_F,P_sig_400(i,1),P_sig_400(i,2)) 
         write(6,*) 'Have done ',i,' events so far'

      enddo
         
      close(22) 
      open(unit=22,file='ZZ_example_events',status='old',err=44) 

!------- BACKGROUND 
       mu_R=2*91.1876d0 
       mu_F=mu_R
       call init_NLOME(mu_R,mu_F) 
!----- change process 
       call swap_process(81,114) 
       call chooser

!------Read in from file (note formating here)

      do i=1,n_events 
         read(22,*,iostat=errno)  p(1,1),p(1,2),p(1,3),p(1,4),
     &        p(2,1),p(2,2),p(2,3),p(2,4),
     &        p(3,1),p(3,2),p(3,3),p(3,4), 
     &        p(4,1),p(4,2),p(4,3),p(4,4),m4l
                 
         if(i.eq.1) write(6,*) 'Calculating For Background Model ' 
           
!------ first entrant in P_back is the weight, second is the MC error
         call LO_LIK(p,mu_R,mu_F,P_back(i,1),P_back(i,2))
           write(6,*) 'Have done ',i,' events so far'

      enddo
        

!------- Calculate Likelihoods 
      
      xs_tot_200=sig_xs_200+back_xs
      xs_tot_400=sig_xs_400+back_xs
       
      Lik_200=0d0 
      Lik_400=0d0 
      Lik_back=0d0
      do i=1,n_events 
         Lik_200=Lik_200+dlog(P_back(i,1)+P_sig_200(i,1))
     &        -dlog(xs_tot_200)
         Lik_400=Lik_400+dlog(P_back(i,1)+P_sig_400(i,1))
     &        -dlog(xs_tot_400)
         Lik_back=Lik_back+dlog(P_back(i,1))
     &        -dlog(back_xs)
      enddo
     
!------- display information 
      
      write(6,*) '****************** OUTPUT *********************' 
      write(6,*) '* Calcualted ',n_events,'                     *'
      write(6,*) '*  P_S (125) = ',P_sig_200(1,1),P_sig_200(1,2) 
      write(6,*) '*  P_S (m4l) = ',P_sig_400(1,1),P_sig_400(1,2) 
      write(6,*) '*  P_B  = ',P_back(1,1),P_back(1,2)
      write(6,*) '* log(Lik(back)/Lik(sig+back)):                *' 
      write(6,*) '*   mh=125 GeV ',-Lik_200+Lik_back,'           *' 
      write(6,*) '*   mh=m4l GeV ',-Lik_400+Lik_back,'           *' 
      write(6,*) '***********************************************'






      return 
!----- errors
 44   write(6,*) 'Couldnt open file'
      stop 

      end
      









      

      
