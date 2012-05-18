c --- Common block for keeping track of weights, used
c --- if unweighting is selected :
      double precision wtmax_cr,newwt_cr
      logical evtgen_cr
      logical unweight_cr
      logical skipnt_cr
      common/maxwt_cr/wtmax_cr,newwt_cr,evtgen_cr,skipnt_cr,unweight_cr

c --- Useful local variables where weights are being checked :
      double precision wtabs_cr
