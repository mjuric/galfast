bin_PROGRAMS = simulate.x
EXTRA_DIST = ri_prior.txt

##--compiler-options
## SUFFIXES = .cu
if HAVE_CUDA
.cu.o:
	$(NVCC) -o $@ -c $(NVCCFLAGS) -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $< -v
	$(NVCC) -M -arch sm_13 -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $< > $(DEPDIR)/$*.Tpo
	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po
	##$(NVCC) -o $<.ptx -ptx --ptxas-options=-v -arch sm_13 -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $<

.cu.lo:
	$(top_srcdir)/cudalt.py $@ $(NVCC) -c $(NVCCFLAGS) -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $<
	$(NVCC) -M -arch sm_13 -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $< > $(DEPDIR)/$*.Tpo
	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po
	##$(NVCC) -o $<.ptx -ptx --ptxas-options=-v -arch sm_13 -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $<

else
.cu.o:
	touch $@
endif

# set the include path found by configure
INCLUDES = -I$(top_srcdir)/src/common $(all_includes)

# the library search path.
SUBDIRS = common gpc

# explicit dependency for CUDA sources
include ./$(DEPDIR)/cuda.Po

noinst_HEADERS = model_J08.h model_densityCube.h model_expDisk.h \
	model_powerLawEllipsoid.h simulate.h simulate_base.h simulate_gpu.cu.h skyconfig_impl.h \
	skygen.cu.h skygen.h
simulate_x_SOURCES = cuda.cu model_J08.cpp model_densityCube.cpp \
	model_expDisk.cpp model_powerLawEllipsoid.cpp observe.cpp observe_FeH.cpp \
	observe_kinTMIII.cpp observe_vel2pm.cpp os_skygen.cpp simulate.cpp simulate.h simulate2.cpp \
	simulate_cpu.cpp
simulate_x_LDADD = $(top_builddir)/src/common/libcommon.a \
	$(top_builddir)/src/gpc/libgpc.a
