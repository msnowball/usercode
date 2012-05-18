VPATH= $(DIRS)
FC = gfortran
CXX = g++
#FFLAGS   = -g  -ffpe-trap=invalid,zero,overflow
FFLAGS  += -fno-automatic -fno-f2c -O2 -I$(INCPATH) -Iobj

HERE   = $(PWD)
SOURCEDIR= $(HERE)/src
INCPATH=$(SOURCEDIR)/Inc
LIBDIR=$(HERE)/Lib 
LIBRARY       = libNLOME.a
OUTPUT_OPTION = -o $(HERE)/obj/$@ 

DIRS = 	$(HERE):\
	$(HERE)/obj:\
        $(SOURCEDIR)/Need:$(SOURCEDIR)/W\
        $(SOURCEDIR)/Setup:$(SOURCEDIR)/Gen_NLO_Weight:\
	$(SOURCEDIR)/ZZ:$(SOURCEDIR)/H:$(SOURCEDIR)/Phase\
	$(SOURCEDIR)/Parton:$(SOURCEDIR)/Integrate:\
	$(SOURCEDIR)/Main:$(SOURCEDIR)/Z:$(SOURCEDIR)/User\
	$(SOURCEDIR)/Zgamma:$(SOURCEDIR)/Wgamma:$(SOURCEDIR)/Red\
	$(SOURCEDIR)/FragDips\
	$(SOURCEDIR)/ggHgaga:$(SOURCEDIR)/BoostFunc\
	$(SOURCEDIR)/Main:$(SOURCEDIR)/Gen_XS\
	$(SOURCEDIR)/Gen_Lik:

GENLIKFILES=\
LO_Lik.o\
NLO_Lik.o


GENXSFILES=\
genLO_int.o\
gen_LO_xs.o\
genNLO_int.o\
gen_NLO_xs.o

USERFILES=\
Cuts_Drell_Yan.o\
Cuts_Vgamma.o\
Cuts_ZZ_4l.o\
photo_iso.o\
frix.o\
iso.o\
my_ptyrap.o

MAINFILES=\
LOME.o\
NLOME.o\
lord_ME.o\
init_NLOME.o\
LO_int_vegas.o\
lord_vegas.o\
virt_vegas.o\
real_vegas.o\
virt_int_vegas.o\
vdips_int.o\
lv_int.o\
vdips_xint.o\
VirtME.o\
setzdips.o\
real_int_3dv.o\
real_int_2dv.o\
RealME.o\
real_int.o\


NEEDFILES=\
docuts.o\
dips_fb_fi.o\
dips_fb_if.o\
transform_fb_if.o\
transform_fb_fi.o\
set_dipsfbps.o\
determine_xminmax.o\
determine_xminmax_r.o\
swap_process.o\
includedipole.o\
readwriteps.o\
masscuts.o\
my_simp_born.o\
donothing_gvec.o\
storeptilde.o\
zeromsq.o\
transform.o\
scaleset.o\
smalls.o\
Gen_mcfm_p.o\
writeinfo.o\
my_simp.o\
alfamz.o\
branch.o\
dipoles.o\
dipolesub.o\
histofin.o\
inv_cuts.o\
boost_Q.o\
check_scheme.o\
couplz.o\
ckmfill.o\
coupling.o\
coupling2.o\
dot.o\
dotem.o\
ddilog.o\
dclaus.o\
getet.o\
higgsp.o\
ptyrap.o\
int_dips.o\
lnrat.o\
i3m.o\
pgen.o\
recastE.o\
read_params.o\
read_cuts.o\
set_epinv.o\
sethparams.o\
spinoru.o\
writeout.o


INTEGRATEFILES=\
dpolint.o\
dromo.o\
dromb.o\
trapzd.o\
midpnt.o\
ebook.o\
lenocc.o\
mbook.o\
ran0.o\
ran1.o\
rn.o\
vegas_cr.o\
vegas.o\

HFILES=\
qqb_hzz.o\
qqb_hzz_v.o\
qqb_hzz_gs.o\
qqb_hzz_g.o\
qqb_hzz_gvec.o\
qqb_hzz_z.o


PARTONFILES = \
checkpath.o \
newton1.o


PARTONFILES +=\
   alfamz.o \
   Ctq4Fn.o \
   Ctq5Par.o \
   Ctq5Pdf.o \
   Cteq6Pdf-2008.o \
   cteq3.o \
   mrs96.o \
   mrs98.o \
   mrs98lo.o \
   mrs98ht.o \
   mrs99.o \
   mrsebh.o \
   mrsg.o \
   mrst2001lo.o \
   jeppelo.o \
   mrst2001.o \
   mrst2002.o \
   mrst2004.o \
   mrst2004f3.o \
   mrst2004f4.o \
   mstwpdf.o \
   mt.o \
   eks98r.o \
   fdist_linux.o \
   pdfwrap_linux.o


PHASEFILES=\
boost.o\
breitw.o\
phi3m0.o\
phase3.o\
phase4.o\
phi3m.o\
phi1_2m.o\
booster.o\
new_gen2l_gam.o\
gen2l_gam.o\
gen2l.o\
gen4l.o\
genP4.o\
genP3.o\
gen4h.o\
gen2flat.o\
gen3flat.o\
gen4flat.o\
genW.o\
phi1_2.o\
fixx1x2.o\
FBPS.o\
phase_weight.o\
my_phase6.o\

GGHGAGAFILES=\
Ftriangle.o\
gg_hgamgam_gs.o\
gg_hgamgam_v.o\
gg_hgamgamg.o\
gg_hgamgam.o\
gg_hgamgam_gvec.o\
gg_hgamgam_z.o\
msqgamgam.o

SETUPFILES=\
basic_setup.o\
chooser.o\
i2mcfm_fill.o\
read_interface_settings.o\
write_i2mcfm_info.o


WGAMMAFILES=\
fagamma.o\
qqb_wfrag.o\
qqb_wgam_fragdips.o\
qqb_wgam_gs.o\
qqb_w_g.o\
qqb_wgam_z.o\
fbgamma.o\
qqb_wgam.o\
qqb_wgam_g.o\
qqb_wgam_v.o

WFILES=\
qqb_w.o\
qqb_w_g.o\
qqb_w_gs.o\
qqb_w_v.o\
qqb_w_z.o


ZFILES=\
qqb_z.o\
qqb_z1jet.o\
qqb_z_gs.o\
qqb_z_gsd.o\
qqb_z_v.o\
qqb_z_z.o

ZZFILES=\
a6g.o\
a61g.o\
a6treeg.o\
a6tree.o\
a6loop.o\
fsc.o\
fcc.o\
fsc1.o\
fsc2.o\
fsc3.o\
fsc4.o\
fsc5.o\
fsc6.o\
fsc7.o\
fsc8.o\
fcc_qpgmgpqm.o\
fcc_qpgpgmqm.o\
fcc_qpgpgpqm.o\
fcc_qpgpqmgm.o\
fcc_qpgpqmgp.o\
fcc_qpqmgmgp.o\
fcc_qpqmgpgm.o\
fcc_qpqmgpgp.o\
fvf.o\
fa.o\
fax.o\
faxsl.o\
L34_12.o\
BigT.o\
xzqqgg.o\
xzqqgg_v.o\
fvs.o\
t.o\
lfunctions.o\
gg_ZZ.o\
qqb_zz.o\
qqb_zz_g.o\
qqb_zz_gs.o\
qqb_zz_v.o\
qqb_zz_z.o\
vvg.o\
vpole.o\
zzamp.o\
zzamps.o\
zzgamp.o

REDFILES=\
branch.o\
ff_alpha.o\
pdfwrap_lhapdf.o\
readcoup.o\
topwidth.o\
checkorder.o\
higgsp.o\
pdfwrap_linux.o\
sethparams.o\
checkversion.o\
newton1.o\
pdfwrap_pdflib.o\
setmb_msbar.o

WGAMMAFILES=\
fagamma.o\
qqb_wfrag.o\
qqb_wgam_fragdips.o\
qqb_wgam_gs.o\
qqb_w_g.o\
qqb_wgam_z.o\
fbgamma.o\
qqb_wgam.o\
qqb_wgam_g.o\
qqb_wgam_v.o

FRAGDIPFILES=\
GGdR_frag.o\
NP_FragSetII.o\
dipolesfrag.o\
fragdriver.o\
locate.o\
store_zdip.o\
NP_FragSetI.o\
Pert_Frag.o\
dopols.o\
int_dips_ga.o\
rescale_pjet.o\
transformfrag.o



ZGAMMAFILES=\
ggZgam_vec.o\
gg_zgam.o\
qqb_zfrag.o\
qqb_zgam.o\
qqb_zgam_fragdips.o\
qqb_zgam_g.o\
qqb_zgam_gs.o\
qqb_zgam_v.o\
qqb_zgam_z.o


OURCODE=$(NEEDFILES) $(SETUPFILES) $(ZZFILES) $(WFILES) $(HFILES)\
	$(PHASEFILES) $(PARTONFILES) $(INTEGRATEFILES) \
	$(ZFILES) $(USERFILES)  $(ZGAMMAFILES)  $(REDFILES)\
	$(WGAMMAFILES) $(FRAGDIPFILES) $(MAINFILES)\
	$(GENXSFILES) $(GGHGAGAFILES) $(GENLIKFILES) \


$(LIBRARY):     $(OURCODE)
		ar cru $(LIBRARY) $(OURCODE) 
		ranlib $(LIBRARY)

%.o: %.f        
		$(FC) $(FFLAGS) -c $< 
# -----------------------------------------------------------------------------




# -----------------------------------------------------------------------------

clean:
	- rm -f *.o *.s *.prj *~ core 

# -----------------------------------------------------------------------------

# DO NOT DELETE
 
