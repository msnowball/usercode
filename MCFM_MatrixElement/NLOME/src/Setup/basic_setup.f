!----------------------- Interface to MCFM ------------------------
!-                                                                -
!- Basic setup Of NLOME  - Prints Header and Citation file       - 
!- Sets up process specific quantities                            -
!-                                                                -
!-                                                                -
!-                                                                -
!-  C. Williams and J. M. Campbell 2011                           -
!------------------------------------------------------------------

      subroutine basic_setup(n1) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f' 
      integer nproc,n1
      integer number_of_particles
      double precision i2mcfm_version 
      logical first 
      
      common/version/i2mcfm_version
      common/nproc/nproc 
!------ Print Header 
 !     call print_header(i2mcfm_version)
      
!------- Determine n_particles 
      n_particles=number_of_particles(n1) 
      n_final=n_particles-2 
!     ------- Setup of process specific variables this will need to happen after 
!------- After things like masses have been defined 
      nproc=n1

      call print_header(1)
      call chooser
      

      return 
      end subroutine 




      subroutine print_header(n2) 
      implicit none 
      integer n1
      double precision n2
      

      write(6,91) '************** NLOME  - version ',n2,'**************'
      write(6,92) '*                                                  *'
      write(6,93) '* NLOME ',n2,'                  May 2012           *'
      write(6,92) '*                                                  *'
      write(6,92) '* Authors: John Campbell, Walter Giele,            *'
      write(6,92) '*          Ciaran Williams                         *'
      write(6,92) '*                                                  *'
      write(6,92) '*  See: 1204.4424 for details on method            *'
      write(6,92) '*                                                  *'
      write(6,92) '* Based on MCFM (Campbell, Ellis, Williams)        *'
      write(6,92) '* please cite following:                           *'
      write(6,92) '****************************************************'
     


 91   format(a32,f3.1,a17)
 93   format(a9,f3.1,a40)
 92   format(a52)

      return 
      end subroutine 



      integer function number_of_particles(n1) 
      implicit none 
      include 'interface_settings.f' 
      integer n1 

      inc_MET=.false.
      number_of_particles=0 
  
      if((n1.eq.1).or.(n1.eq.6)) then 
!     1  '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))' 'N'
!     6  '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))' 'N'
         number_of_particles=4 
         inc_MET=.true. 
         no_neutrinos=1
      elseif((n1.ge.11).and.(n1.le.19)) then 
!     11 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+f(p5)'  'N'
!     13 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+c~(p5)' 'N'
!     14 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+c~(p5) [massless]' 'L'
!     16 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+f(p5)' 'N'
!     18 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+c(p5)' 'N'
!     19 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+c(p5) [massless]' 'L'
         number_of_particles=5 
      elseif(((n1.ge.20).and.(n1.le.22)).or.((n1.eq.26).or.(n1.eq.27)))
     &        then
!     20 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +b(p5)+b~(p6) [massive]' 'N'
!     21 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +b(p5)+b~(p6)' 'N'
!     22   '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +f(p5)+f(p6)' 'N'
!     25 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) +b(p5)+b~(p6) [massive]' 'N'
!     26 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) +b(p5)+b~(p6)' 'N'
!     27 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) +f(p5)+f(p6)' 'N'
         number_of_particles=6
      elseif((n1.ge.31).and.(n1.le.35)) then 
         number_of_particles=4 
!     31 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))' 'N'
!     32 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4)))' 'N'
!     33   '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4))' 'N'
!     34 '  f(p1)+f(p2) --> Z^0(-->3*(d(p5)+d~(p6)))' 'N'
!     35 '  f(p1)+f(p2) --> Z^0(-->2*(u(p5)+u~(p6)))' 'N'
      elseif((n1.ge.41).and.(n1.le.43)) then 
         number_of_particles=5
!     41 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+f(p5)' 'N'
!     42 '  f(p1)+f(p2) --> Z_0(-->3*(nu(p3)+nu~(p4)))+f(p5)' 'N'
!     43 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4))+f(p5)' 'N'
      elseif((n1.ge.61).and.(n1.le.89)) then 
!     61 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->e^-(p5)+nu~(p6))' 'N'
!     62 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->q(p5)+q~(p6))' 'N'
!     63 '  f(p1)+f(p2) --> W^+(--> q(p3)+ q~(p4)) +W^-(-->e^-(p5)+nu~(p6))' 'N'
!     64 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->e^-(p5)+nu~(p6)) [no pol]' 'L'
!66   '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->e^-(p5)+nu~(p6))+f(p7)' 'L'
!
!     71 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->e^-(p5)+e^+(p6))' 'N'
!     72 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->3*(nu_e(p5)+nu~_e(p6)))' 'N'
!     73 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->b(p5)+b~(p6))' 'N'
!     74 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->3*(d(p5)+d~(p6)))' 'N'
!     75 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->2*(u(p5)+u~(p6)))' 'N'
!
!     76 '  f(p1)+f(p2) --> W^-(-->mu^-(p3)+nu~(p4))+Z^0(-->e^-(p5)+e^+(p6))' 'N'
!77   '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+Z^0(-->3*(nu(p5)+nu~(p6)))' 'N'
!     78 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+Z^0(-->b(p5)+b~(p6))' 'N'
!     79 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+Z^0(-->3*(d(p5)+d~(p6)))' 'N'
!     80 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+Z^0(-->2*(u(p5)+u~(p6)))' 'N'
!
!     81 '  f(p1)+f(p2) --> Z^0(-->mu^-(p3)+mu^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))' 'N'
!     82 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->3*(nu(p5)+nu~(p6)))' 'N'
!     83 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->b(p5)+b~(p6))' 'N'
!     84 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4)) + Z^0(-->3*(nu(p5)+nu~(p6)))' 'N'
!
!     86 '  f(p1)+f(p2) --> Z^0(-->mu^-(p3)+mu^+(p4))+Z^0(-->e^-(p5)+e^+(p6))[no gamma^*]' 'N'
!     87 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->3*(nu(p5)+nu~(p6))) [no gamma^*]' 'N'
!88   '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+Z^0(-->b(p5)+b~(p6)) [no gamma^*]' 'N'
!     89 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4)) + Z^0(-->3*(nu(p5)+nu~(p6))) [no gamma^*]' 'N'
         number_of_particles=6 
      elseif((n1.eq.91).or.(n1.eq.96)) then 
!     91 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + H(-->b(p5)+b~(p6))' 'N'
!     96 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + H(-->b(p5)+b~(p6))'  'N'
         number_of_particles=6
      elseif((n1.eq.92).or.(n1.eq.97)) then
!     92 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + H(-->W^+(nu(p5),e^+(p6))W^-(e^-(p7),nu~(p8)))' 'N'
!     97 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + H(-->W^+(nu(p5),e^+(p6))W^-(e^-(p7),nu~(p8)))' 'N'
         number_of_particles=8 
      elseif((n1.ge.101).and.(n1.le.103)) then 
!     101 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + H(-->b(p5)+b~(p6))' 'N'
!     102 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))) + H(-->b(p5)+b~(p6))' 'N'
!     103 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4)) + H(-->b(p5)+b~(p6))'      'N'
         number_of_particles=6 
      elseif((n1.ge.106).and.(n1.le.108)) then 
!       106 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + H(-->W^+(nu(p5),e^+(p6))W^-(e^-(p7),nu~(p8)))' 'N'
!     107 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))) + H(-->W^+(nu(p5),e^+(p6))W^-(e^-(p7),nu~(p8)))' 'N'
!     108 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4)) + H(-->W^+(nu(p5),e^+(p6))W^-(e^-(p7),nu~(p8)))'      'N'
         number_of_particles=8 
      elseif((n1.eq.111).or.(n1.eq.112)) then 
         number_of_particles=4 
!     111  '  f(p1)+f(p2) --> H(-->b(p3)+b~(p4))' 'N'
!     112 '  f(p1)+f(p2) --> H(-->tau^-(p3)+tau^+(p4))' 'N'
      elseif((n1.ge.113).and.(n1.le.116).or.(n1.eq.118)) then 
!     113 '  f(p1)+f(p2) --> H(--> W^+(nu(p3)+e^+(p4)) + W^-(e^-(p5)+nu~(p6)))' 'N'
!     114 '  f(p1)+f(p2) --> H(-->Z^0(mu^-(p3)+mu^+(p4)) + Z^0(e^-(p5)+e^+(p6))' 'N'
!     115 '  f(p1)+f(p2) --> H(-->Z^0(3*(nu(p3)+nu~(p4)))+ Z^0(e^-(p5)+e^+(p6))' 'N'
!     116 '  f(p1)+f(p2) --> H(-->Z^0(mu^-(p3)+mu^+(p4)) + Z^0(b(p5)+b~(p6))' 'N'
         number_of_particles=6
      elseif((n1.ge.151).and.(n1.le.153)) then 
         number_of_particles=8
!     151 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->b~(p6)+e^-(p7)+nu~(p8))' 'N'
!     152 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->b~(p6)+q(p7)+q~(p8))' 'N'
!     153 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->b~(p6)+e^-(p7)+nu~(p8)) (uncorr)' 'N'
      elseif((n1.ge.161).and.(n1.le.168)) then 
         number_of_particles=6
!     161 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+q(p6) [t-channel]' 'N'
!     162 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+q(p6) [decay]' 'N'
!     163 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+q(p6) [t-channel] mb>0' 'N'
!     166 '  f(p1)+f(p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+q(p6) [t-channel]' 'N'
!     167 '  f(p1)+f(p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+q(p6) [decay]' 'N'
!     168 '  f(p1)+f(p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+q(p6) [t-channel] mb>0' 'N'

      elseif((n1.eq.203).or.(n1.eq.204)) then 
!     203 '  f(p1)+f(p2) --> H(-->b(p3)+b~(p4)) + f(p5)' 'N'
!     204 '  f(p1)+f(p2) --> H(-->tau^-(p3)+tau^+(p4)) + f(p5)' 'N'
         number_of_particles=5
      elseif((n1.eq.208).or.(n1.eq.209)) then 
!     208 '  f(p1)+f(p2) --> H(-->W^+(nu(p3),e^+(p4))W^-(e^-(p5),nu~(p6)))+f(p7)' 'N'
!     209 '  f(p1)+f(p2) --> H(-->Z^+(e^-(p3),e^+(p4))Z(mu^-(p5),mu^+(p6)))+f(p7)' 'N'
         number_of_particles=7 
      elseif((n1.eq.271).or.(n1.eq.272)) then 
!     271 '  f(p1)+f(p2) --> H(b(p3)+b~(p4))+f(p5)+f(p6) [in heavy top limit]' 'N'
!     272 '  f(p1)+f(p2) --> H(tau^-(p3)+tau^+(p4))+f(p5)+f(p6) [in heavy top limit]' 'N'
         number_of_particles=6
      elseif((n1.eq.273).or.(n1.eq.274)) then 
         number_of_particles=8 
!     273 '  f(p1)+f(p2) --> H(-->W^+(nu(p3),e^+(p4))W^-(e^-(p5),nu~(p6)))+f(p7)+f(p8)' 'N'
!     274 '  f(p1)+f(p2) --> H(-->Z^+(e^-(p3),e^+(p4))Z(mu^-(p5),mu^+(p6)))+f(p7)+f(p8)' 'N'
         elseif((n1.eq.300)) then 
            number_of_particles=5 
            

      else
         write(6,*) 'Unrecognised number of particles' 
         stop
      endif

      if(number_of_particles.lt.1) then 
         write(6,*) 'Couldnt read the number of particles. ' 
         write(6,*) 'Please convert process number to MCFM notation' 
         stop 
      endif
         

      return 
      end function
