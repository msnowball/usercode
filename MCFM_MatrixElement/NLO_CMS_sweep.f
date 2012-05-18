!------ program for CMS 
!------ read in events 
!------ for each event calculate Higgs at 125, Higgs at m4l Background 
!------ write outputs 

      program NLO_CMS_ex
      implicit none 
      include 'NLOME/src/Inc/masses.f'
   
      integer i,nargs,istem
      character*200 inputfile,outfile_sig1,outfile_sig2
      character*200 outfile_back,datafile
      integer nu
      double precision mu_r,mu_f
      character*200 path_base,dir
      integer ipath,idir
      integer n_calc 
      double precision m4l,m4l_calc 
      double precision p(4,4),nlo_wt,nlo_err
      logical nomore 
      character*1 ic

      nomore=.false. 
      nargs=iargc() 

      if(nargs.eq.1) then 
         call getarg(1,datafile) 
      else
         stop 
      endif

      path_base='/scratch/ciaran/CMS/DATA/' 
      dir='/scratch/ciaran/CMS/Out/'
      
      idir=LEN_TRIM(dir)
      ipath=LEN_TRIM(path_base)
      istem=LEN_TRIM(datafile) 

      inputfile=path_base(1:ipath)//datafile
 
      
      do i=1,9
         n_calc=0
         write(ic,'(I1)') i
         outfile_sig1=dir(1:idir)//datafile(1:istem)//'_sig_12'//ic 
         
!-------- open files      
         open(unit=22,file=inputfile,status='old',err=93) 
         open(unit=23,file=outfile_sig1,status='unknown') 
!------setup Higgs mass at mh
         m4l=dfloat(120+i) 
         mu_r=m4l 
         mu_f=m4l
         call init_NLOME(mu_r,mu_f) 
         hmass=m4l 
         call chooser 
         
!--------------- Signal generated at each m_4l 

 11      continue 
      
!---- Read event 
         call read_cmsstyle(22,p,m4l,nomore)
         if(nomore) goto 22 
!---- check m4l 
         m4l_calc=0d0
         do nu=1,3 
            m4l_calc=m4l_calc-(p(1,nu)+p(2,nu)+p(3,nu)+p(4,nu))**2 
         enddo
         nu=4 
         m4l_calc=m4l_calc+(p(1,nu)+p(2,nu)+p(3,nu)+p(4,nu))**2 
         if((dsqrt(m4l_calc)-m4l).gt.1d-1) then 
            write(6,*) 'Bugger, looks like m4l is wrong',m4l,m4l_calc
            stop
         endif
      
      
      call NLO_LIK(p,mu_r,mu_f,nlo_wt,nlo_err) 
      call write_cms(23,p,nlo_wt,nlo_err) 
      write(6,*) nlo_wt,nlo_err
  
      n_calc=n_calc+1 
      write(6,*) 
      write(6,*) 'Calculated ',n_calc,' so far ' 
      write(6,*) 
      goto 11
 22   continue 

      close(22)
      close(23) 
  
      enddo

      return 
 93   write(6,*) 'Couldnt open inputfile' 

      end 


      subroutine read_cmsstyle(unit,p,m4l,nomore) 
      implicit none 
      integer unit,errno 
      double precision m4l,p(4,4) 
      logical nomore 

      nomore=.false. 
      read(unit,*,iostat=errno) p(1,1),p(1,2),p(1,3),p(1,4),
     &  p(2,1),p(2,2),p(2,3),p(2,4),
     &  p(3,1),p(3,2),p(3,3),p(3,4),
     &  p(4,1),p(4,2),p(4,3),p(4,4),m4l 

      if(errno.ne.0) then 
         nomore=.true. 
      endif
      return 
      end 

      

      subroutine write_cms(unit,p,nlo_lik,nlo_err) 
      implicit none 
      integer unit,errno
      double precision nlo_lik,nlo_err,p(4,4)  

       
      write(unit,*,iostat=errno) p(1,1),p(1,2),p(1,3),p(1,4),
     &  p(2,1),p(2,2),p(2,3),p(2,4),
     &  p(3,1),p(3,2),p(3,3),p(3,4),
     &  p(4,1),p(4,2),p(4,3),p(4,4),nlo_lik,nlo_err 

     
      return 
      end 
 




      

      


