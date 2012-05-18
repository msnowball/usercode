      subroutine writeinfo(unitno,commchars,xsec,xsec_err,itno)
************************************************************************
*   Routine to write out run information to a desired unit             *
************************************************************************
      implicit none
      include 'PDFerrors.f'
      include 'process.f'
      integer unitno,j,k,itno
      double precision xsec,xsec_err
      double precision lord_bypart(-1:1,-1:1),lordnorm,rescale
      double precision ggpart,gqpart,qgpart,qqpart,qqbpart,
     . gqbpart,qbgpart,qbqbpart,qbqpart
      
      character*2 commchars
      character*4 part
      character*30 runstring
      logical creatent,dswhisto,dryrun,makecuts
      integer nproc,ih1,ih2,itmx1,itmx2,ncall1,ncall2,origij
      integer NPTYPE,NGROUP,NSET
      double precision sqrts
      double precision Rcut
 
      common/outputflags/creatent,dswhisto      

      common/nproc/nproc
      common/part/part
      common/runstring/runstring
      common/energy/sqrts
      common/density/ih1,ih2
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/dryrun/dryrun
      
      common/pdflib/NPTYPE,NGROUP,NSET
      
      common/Rcut/Rcut
      common/makecuts/makecuts

      common/origij/origij

      common/bypart/lord_bypart
      common/finalpart/ggpart,gqpart,qgpart,qqpart,qqbpart

      return 
c$$$      if (itno .gt. 0) then
c$$$c--- write warning that result is only intermediate; populate the
c$$$c--- variables in finalpart (normally done in mcfm_exit)
c$$$      write(unitno,*) commchars//
c$$$     & ' Intermediate result for iteration',itno,')'
c$$$      lordnorm=0d0
c$$$      do j=-1,1
c$$$      do k=-1,1
c$$$        lordnorm=lordnorm+lord_bypart(j,k)
c$$$      enddo
c$$$      enddo
c$$$      ggpart=lord_bypart( 0, 0)/lordnorm
c$$$      gqpart=lord_bypart( 0,+1)/lordnorm
c$$$      gqbpart=lord_bypart( 0,-1)/lordnorm
c$$$      qgpart=lord_bypart(+1, 0)/lordnorm
c$$$      qbgpart=lord_bypart(-1, 0)/lordnorm
c$$$      qqpart=lord_bypart(+1,+1)/lordnorm
c$$$      qbqbpart=lord_bypart(-1,-1)/lordnorm
c$$$      qqbpart=lord_bypart(+1,-1)/lordnorm
c$$$      qbqpart=lord_bypart(-1,+1)/lordnorm
c$$$      endif
c$$$      write(unitno,55) commchars//
c$$$     & ' Cross-section is: ',xsec,' +/-',xsec_err,')'
c$$$      write(unitno,*)
c$$$
c$$$c--- for gg->H+X processes, also write out the cross section
c$$$c---  normalized by sigma(gg->H, finite mt)/sigma(gg->H, mt-> infinity)
c$$$      if (((case(1:5) .eq. 'ggfus') .or. (case(1:3) .eq. 'HWW')
c$$$     & .or.(case(1:3) .eq. 'HZZ')) .and. (case .ne. 'HWWint')
c$$$     &  .and. (case .ne. 'HWW_tb') .and. (case .ne. 'HZZint')
c$$$     &  .and. (case .ne. 'HZZ_tb') ) then
c$$$        call finitemtcorr(rescale)
c$$$        write(unitno,55) commchars//'Rescaled x-sec is:',
c$$$     &     xsec*rescale,' +/-',xsec_err*rescale,')'
c$$$        write(unitno,*)
c$$$      endif
c$$$     
c$$$      write(unitno,*) commchars,
c$$$     &                ' Contribution from parton sub-processes:'
c$$$      write(unitno,95)commchars,'   GG    ',ggpart*xsec,ggpart*100d0
c$$$      write(unitno,95)commchars,'   GQ    ',gqpart*xsec,gqpart*100d0
c$$$      write(unitno,95)commchars,'   GQB   ',gqbpart*xsec,gqbpart*100d0
c$$$      write(unitno,95)commchars,'   QG    ',qgpart*xsec,qgpart*100d0
c$$$      write(unitno,95)commchars,'   QBG   ',qbgpart*xsec,qbgpart*100d0
c$$$      write(unitno,95)commchars,'   QQ    ',qqpart*xsec,qqpart*100d0
c$$$      write(unitno,95)commchars,'   QBQB  ',qbqbpart*xsec,qbqbpart*100d0
c$$$      write(unitno,95)commchars,'   QQB   ',qqbpart*xsec,qqbpart*100d0
c$$$      write(unitno,95)commchars,'   QBQ   ',qbqpart*xsec,qbqpart*100d0
c$$$      write(unitno,*)
c$$$
c$$$      if (PDFerrors) then
c$$$        do j=0,maxPDFsets
c$$$          write(unitno,56) j,PDFxsec(j)
c$$$        enddo
c$$$        write(unitno,*)
c$$$      endif
c$$$
c$$$      if (commchars .eq. ' (') then
c$$$c--- new routine for writing out contents of input file
c$$$        call writeinput(unitno,' (',' )','WRITEALL')
c$$$      else
c$$$        call writeinput(unitno,commchars,'  ','WRITEALL')
c$$$      endif
c$$$
c$$$c--- old lines for writing out inputs
c$$$
c$$$c      write(unitno,*) '( Run corresponds to this input file)'
c$$$c      write(unitno,*)
c$$$c      write(unitno,*)
c$$$c     . '( [Flags to specify the mode in which MCFM is run] )'
c$$$c      write(unitno,98) evtgen,'evtgen'
c$$$c      write(unitno,98) creatent,'creatent'
c$$$c      write(unitno,98) skipnt,'skipnt'
c$$$c      write(unitno,98) dswhisto,'dswhisto'
c$$$c
c$$$c      write(unitno,*)
c$$$c      write(unitno,*)
c$$$c     . '( [General options to specify the process and execution] )'
c$$$c      write(unitno,97) nproc,'nproc'
c$$$c      write(unitno,96) part,'part'
c$$$c      write(unitno,96) runstring,'runstring'
c$$$c      write(unitno,99) sqrts,'sqrts'
c$$$c      write(unitno,97) ih1,'ih1'
c$$$c      write(unitno,97) ih2,'ih2'
c$$$c      write(unitno,99) hmass,'hmass'
c$$$cc--- catch special scale choices for stop+b process
c$$$c      if ( (nproc .eq. 231) .or. (nproc .eq. 236)      
c$$$c     . .or.(nproc .eq. 241) .or. (nproc .eq. 246)      
c$$$c     . .or.(nproc .eq. 242) .or. (nproc .eq. 247) ) then
c$$$c         write(unitno,99) renscale_L,'renscale_L'
c$$$c         write(unitno,99) facscale_L,'facscale_L'
c$$$c         write(unitno,99) renscale_H,'renscale_H'
c$$$c         write(unitno,99) facscale_H,'facscale_H'
c$$$c      else
c$$$c         write(unitno,99) scale,'scale'
c$$$c         write(unitno,99) facscale,'facscale'
c$$$c      endif
c$$$c
c$$$c      write(unitno,98) dynamicscale,'dynamicscale'
c$$$c      write(unitno,98) zerowidth,'zerowidth'
c$$$c      write(unitno,98) removebr,'removebr'
c$$$c      write(unitno,97) itmx1,'itmx1'
c$$$c      write(unitno,97) ncall1,'ncall1'
c$$$c      write(unitno,97) itmx2,'itmx2'
c$$$c      write(unitno,97) ncall2,'ncall2'
c$$$c      write(unitno,97) origij,'ij'
c$$$c      write(unitno,98) dryrun,'dryrun'
c$$$c      write(unitno,98) Qflag,'Qflag'
c$$$c      write(unitno,98) Gflag,'Gflag'
c$$$c      
c$$$c      write(unitno,*)
c$$$c      write(unitno,*) 
c$$$c     . '( [Heavy quark masses] )'
c$$$c      write(unitno,99) mt,'top mass'
c$$$c      write(unitno,99) mb,'bottom mass'
c$$$c      write(unitno,99) mc,'charm mass'
c$$$c
c$$$c      write(unitno,*)
c$$$c      write(unitno,*) 
c$$$c     . '( [Pdf selection] )'
c$$$c      write(unitno,96) pdlabel,'pdlabel '
c$$$c      write(unitno,97) NGROUP,'NGROUP'
c$$$c      write(unitno,97) NSET,'NSET'
c$$$c      write(unitno,96) PDFname,'LHAPDF group'
c$$$c      write(unitno,97) PDFmember,'LHAPDF set'
c$$$c
c$$$c      write(unitno,*)
c$$$c      write(unitno,*)
c$$$c     . '( [Jet definition and event cuts] )'
c$$$c      write(unitno,99) dsqrt(wsqmin),'m34min'
c$$$c      write(unitno,99) dsqrt(wsqmax),'m34max'
c$$$c      write(unitno,99) dsqrt(bbsqmin),'m56min'
c$$$c      write(unitno,99) dsqrt(bbsqmax),'m56max'
c$$$c      write(unitno,98) inclusive,'inclusive'
c$$$c      write(unitno,96) algorithm,'algorithm'
c$$$c      write(unitno,99) ptjetmin,'ptjetmin'
c$$$c      write(unitno,99) etajetmax,'etajetmax'
c$$$c      write(unitno,99) Rcut,'Rcut'
c$$$c      write(unitno,98) makecuts,'makecuts'
c$$$c      write(unitno,99) leptpt,'leptpt'
c$$$c      write(unitno,99) leptrap,'leptrap'
c$$$c      write(unitno,99) misspt,'misspt'
c$$$c      write(unitno,99) leptpt2,'leptpt2'
c$$$c      write(unitno,99) leptrap2,'leptrap2'
c$$$c      write(unitno,99) Rjlmin,'Rjlmin'
c$$$c      write(unitno,99) Rllmin,'Rllmin'
c$$$c      write(unitno,99) delyjjmin,'delyjjmin'
c$$$c      write(unitno,98) jetsopphem,'jetsopphem'
c$$$c      write(unitno,97) lbjscheme,'lbjscheme'
c$$$c      write(unitno,99) gammpt,'gammpt'
c$$$c      write(unitno,99) gammrap,'gammrap'
c$$$c      write(unitno,99) gammcone,'gammcone'
c$$$c      write(unitno,99) gammcut,'gammcut'
c$$$c
c$$$c      write(unitno,*)
c$$$c      write(unitno,*)
c$$$c     . '( [Anomalous couplings of the W and Z] )'
c$$$c      write(unitno,99) delg1_z,'delg1_z'
c$$$c      write(unitno,99) delk_z,'delk_z'
c$$$c      write(unitno,99) delk_g,'delk_g'
c$$$c      write(unitno,99) lambda_z,'lambda_z'
c$$$c      write(unitno,99) lambda_g,'lambda_g'
c$$$c      write(unitno,99) tevscale,'tevscale'
c$$$c
c$$$c      write(unitno,*)
c$$$c      write(unitno,*) 
c$$$c     . '( [How to resume/save a run] )'
c$$$c      write(unitno,98) readin,'readin'
c$$$c      write(unitno,98) writeout,'writeout'
c$$$c      write(unitno,96) ingridfile,'ingridfile'
c$$$c      write(unitno,96) outgridfile,'outgridfile'
c$$$c
c$$$c      write(unitno,99)
c$$$
c$$$      return
c$$$
c$$$c--- 55 format
c$$$   55 format(a20,f24.10,a4,f24.10,a1)
c$$$c--- 56 character format
c$$$   56 format('( PDF error set ',i3,'  --->',f13.3,' fb  )')
c$$$c--- 95 character format
c$$$   95 format(a2,5x,a9,' |',f18.5,f8.2,'%')
c$$$c--- 96 character format      
c$$$   96 format(' (',a20,12x,'[',a,']',' )')  
c$$$c--- 97 integer format      
c$$$   97 format(' (',i20,12x,'[',a,']',' )')  
c$$$c--- 98 logical format      
c$$$   98 format(' (',L20,12x,'[',a,']',' )')  
c$$$c--- 99 floating point format
c$$$   99 format(' (',f20.4,12x,'[',a,']',' )')  
c$$$      
      end
      
      
