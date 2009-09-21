bin_PROGRAMS = galfast.x
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
INCLUDES = -I$(top_srcdir)/src/common -I$(top_srcdir)/src/modules \
	-I$(top_srcdir)/src/skygen $(all_includes)

# the library search path.
SUBDIRS = common gpc modules skygen

# explicit dependency for CUDA sources
include ./$(DEPDIR)/kernels_gpu.Po

noinst_HEADERS = simulate_gpu.cu.h sph_polygon.h pipeline.h
galfast_x_SOURCES = footprint.cpp galfast.cpp kernels_cpu.cpp kernels_gpu.cu \
	pipeline.cpp tests.cpp
galfast_x_LDADD = $(top_builddir)/src/skygen/libskygen.la \
	$(top_builddir)/src/common/libcommon.a $(top_builddir)/src/gpc/libgpc.a \
	$(top_builddir)/src/modules/libmodules.la
