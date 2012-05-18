

      subroutine swap_process(n_swap) 
      implicit none 
      include 'process.f'
      integer nproc 
      common/nproc/nproc 
      integer n_old,n_swap 
      logical first 
      data first /.true./ 
      save first 


      n_old = nproc 
      nproc = n_swap
      if(nproc.eq.n_old) return 

      if(first) then 
         write(6,*) '************************************************' 
         write(6,*) '* Swaping nrpoc from ',n_old,' to ',n_swap,' *' 
         write(6,*) '* (Information only displayed on first call)    *' 
         write(6,*) '************************************************'
         first=.false. 
      endif 

!---- currently only implemented for Higgs (ZZ 4l) 
   
      if((n_old.eq.81).and.(nproc.eq.114)) then 
 !        write(6,*) 'in case ',case 
         case='HZZ_4l'
 !        write(6,*) 'out case ',case
      elseif((n_old.eq.114).and.(nproc.eq.81)) then 
         case='ZZlept' 
      else
         write(6,*) 'Sorry couldnt swap from ',n_old, ' to ',nproc
         write(6,*) 'Check doccumentation ' 
         stop
      endif

      
      return 
      end 
      
      
      
