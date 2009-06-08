bin_PROGRAMS = simulate.x
EXTRA_DIST = ri_prior.txt

##--compiler-options
## SUFFIXES = .cu
if HAVE_CUDA
.cu.o:
	$(NVCC) -o $<.ptx -ptx --ptxas-options=-v -arch sm_13 -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $<
	$(NVCC) -o $@ -c $(NVCCFLAGS) -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $<
	$(NVCC) -M -arch sm_13 -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $< > $(DEPDIR)/$*.Tpo
	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po

.cu.lo:
	$(NVCC) -o $<.ptx -ptx --ptxas-options=-v -arch sm_13 -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $<
	$(top_srcdir)/cudalt.py $@ $(NVCC) -c $(NVCCFLAGS) -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $<
	$(NVCC) -M -arch sm_13 -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $< > $(DEPDIR)/$*.Tpo
	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po

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

noinst_HEADERS = simulate.h simulate_base.h skygen.h simulate_gpu.cu.h skygen.cu.h 
simulate_x_SOURCES = observe.cpp observe_FeH.cpp observe_kinTMIII.cpp \
	observe_vel2pm.cpp os_skygen.cpp simulate.cpp simulate.h simulate2.cpp simulate_cpu.cpp \
	cuda.cu
simulate_x_LDADD = $(top_builddir)/src/gpc/libgpc.a \
	$(top_builddir)/src/common/libcommon.a
