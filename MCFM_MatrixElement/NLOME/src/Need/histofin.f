      block data linlog_data
      implicit none
      include 'nplot.f'
      data linlog/nplot*'lin'/
      end

      subroutine histofin(xsec,xsec_err,itno,itmx)
c--- This outputs the final histograms for itno=0
c--- For itno>0, this is an intermediate result only
      implicit none
      include 'verbose.f'
      include 'PDFerrors.f'
      include 'histo.f'
      include 'outputoptions.f'
      include 'vanillafiles.f'
      integer j,nlength,itno,itmx,nplotmax,nempty
      character*255 runname,outfiledat,outfiletop,outfileerr
c--F  Add gnuplot output and root output
      character*255, outfilegnuplot, outfileps
      character*255, outfileroot, outfilerootC
c--F
      character*3 oldbook
      character mop
      double precision xsec,xsec_err,scalefac,itscale
      double precision EHIST(4,1000,100)   
      integer IHISTOMATCH(100),ICOUNTHISTO  
      logical scaleplots                  
      common/runname/runname
      common/nlength/nlength
      common/nplotmax/nplotmax
      COMMON/EHISTO/EHIST,IHISTOMATCH,ICOUNTHISTO
      common/scaleplots/scalefac,scaleplots
      
      if (itno .eq. 0) then
      write(6,*)
      write(6,*) '****************************************************'
      if (vanillafiles) then
        write(6,*) 'output file name normally  ',runname(1:nlength)
        runname='mcfm-output'
	nlength=11
	write(6,*)
        write(6,*) '  but renamed to: ',runname
      else
        write(6,*) 'output files  ',runname(1:nlength)
      endif
      write(6,*) '****************************************************'
      call flush(6)
      scaleplots=.false.
      else
      scaleplots=.true.
      scalefac=1d0/dfloat(itno)
      endif

      outfiledat=runname
      outfiletop=runname
      outfilegnuplot=runname
      outfileps=runname
      outfileroot=runname
      outfilerootC=runname
      outfileerr=runname
      outfiledat(nlength+1:nlength+4)='.dat'
      outfiletop(nlength+1:nlength+4)='.top'
      outfilegnuplot(nlength+1:nlength+4)='.gnu'
      outfileps(nlength+1:nlength+3)='.ps'
      outfileroot(nlength+1:nlength+5)='.root'
      outfilerootC(nlength+1:nlength+2)='.C'
      outfileerr(nlength+1:nlength+10)='_error.top'
      

      if ((PDFerrors) .and. (ICOUNTHISTO .gt. 0)) then
        open(unit=95,file=outfileerr,status='unknown')
      endif

      if (writedat) then
        open(unit=98,file=outfiledat,status='unknown')
      endif
      if (writetop) then
      open(unit=99,file=outfiletop,status='unknown')
      endif
      if (writegnu) then
      open(unit=97, file=outfilegnuplot,status='unknown')
      endif
      if (writeroot) then
      open(unit=96, file=outfilerootC, status='unknown')
      endif
      
c--- write out run info to top of files
      if (writedat) then
      call writeinfo(98,' (',xsec,xsec_err,itno)      
      endif
      if (writetop) then
      call writeinfo(99,' (',xsec,xsec_err,itno)      
      endif
      if (writegnu) then
      call writeinfo(97,'# ',xsec,xsec_err,itno) 
      write(97,120) outfileps(1:nlength+3)
      endif
      if (writeroot) then
      call writeinfo(96,'//',xsec,xsec_err,itno) 
      write(96,121) outfileroot(1:nlength+5)
      endif

  120 FORMAT (/1x,
     & ' set terminal postscript col enhanced', /1x,
     & ' set output "', A, '"', /1X,
     & ' set style data points', /1X,
     & ' set key off')
  121 FORMAT (/1x,
     & ' {', /1x,
     & ' mcfmhisto = new TFile("', A, '", "recreate");',/1x,
     & ' mcfmhisto -> cd();',/1x,
     & ' histos = new TObjArray(0);',/1x)


c--- make sure to scale results by the maximum number of iterations
      if (itno .eq. 0) then
        itscale=1d0/dfloat(itmx)
        mop='V'
      else
        itscale=1d0/dfloat(itno)
        mop='U'
      endif
      
c--- calculate the errors in each plot (and store in 2*maxhisto+j)    
      do j=1,nplotmax
      if (verbose) then
c        write(6,*) 'Calculating errors for plot ',j
        call flush(6)
      endif
      call mopera(j,mop,maxhisto+j,2*maxhisto+j,itscale,1d0)
      enddo

      do j=1,nplotmax
      if (verbose) then
c        write(6,*) 'Finalizing plot ',j
        call flush(6)
      endif
c--- ensure that MFINAL doesn't turn off booking for intermediate results
      oldbook=book(j)
      call mfinal(j)
      if (itno .gt. 0) then
      book(j)=oldbook
      endif
      enddo

      nempty=0
      do j=1,nplotmax
      if (verbose) then
c        write(6,*) 'Writing .dat for plot ',j
        call flush(6)
      endif
      if (writedat) then
        call mprint(j,2*maxhisto+j)
      endif
c      if (book(j) .ne. 'YES') nempty=nempty+1
      enddo
      
      if (writedat) close(unit=98)

c---generate topdrawer file - only for non-empty plots
c--F also gnuplot files
      do j=1,nplotmax-nempty
      if (verbose) then
c        write(6,*) 'Writing .top for plot ',j
        call flush(6)
      endif
      if (writetop) then
        call mtop(j,2*maxhisto+j,'x','y',linlog(j))
      endif
      if (writegnu) then
        call mgnuplot(j,2*maxhisto+j,'x','y',linlog(j))
      endif
      if (writeroot) then
        call mrootplot(j,2*maxhisto+j,'x','y')
      endif
      if ((PDFerrors) .and. (IHISTOMATCH(j) .ne. 0)) then
        call emtop(j,2*maxhisto,'x','y',linlog(j))
      endif
      enddo
      if (writetop) close(unit=99)
      if (writegnu) close(unit=97)

      if (writeroot) then
c--F  closing statements for root file
        write(96,*) ' mcfmhisto -> cd();'
        write(96,*) ' if (histos -> GetEntries() > 0 ) then  {'
        write(96,*) '  histos->Write();'
        write(96,*) '  mcfmhisto -> Close();'
        write(96,*) ' }'
        write(96,*) '}'
        close(unit=96)
      endif
      
c---generate error file
      if ((PDFerrors) .and. (ICOUNTHISTO .gt. 0)) then
        do j=1,nplotmax
          if (IHISTOMATCH(j) .ne. 0) then
            if (verbose) then
c              write(6,*) 'Writing .top for plot ',j
              call flush(6)
            endif
            call etop(j,2*maxhisto,'x','y',linlog(j))
          endif
        enddo
        close (unit=95)
      endif
      
      return
      end

