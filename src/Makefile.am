bin_PROGRAMS = galfast.x
EXTRA_DIST = ri_prior.txt

BUILT_SOURCES = version.stamp
DISTCLEANFILES = version.h version.h.tmp

## Generate a version string from git hashes
version.stamp:
	@ echo '#ifndef version_h__' >  version.h.tmp
	@ echo '#define version_h__' >>  version.h.tmp
	@ echo '#define VERSION_BASE "'`git describe`'"' >> version.h.tmp
	@ echo '#define LINES_MODIFIED_SINCE_CHECKOUT '`git diff | wc --lines` >> version.h.tmp
	@ echo '#define DIFF_SHA1 "'`git diff | sha1sum | cut -b 1-8`'"' >> version.h.tmp
	@ echo 'static const char *VERSION_STRING = LINES_MODIFIED_SINCE_CHECKOUT ? VERSION_BASE "-diffsha=" DIFF_SHA1 : VERSION_BASE;' >> version.h.tmp
	@ echo '#endif' >>  version.h.tmp
	@ test -f version.h || touch version.h
	@ cmp -s version.h.tmp version.h || mv version.h.tmp version.h

##--compiler-options
## SUFFIXES = .cu
if HAVE_CUDA
.cu.o:
	$(NVCC) -o $@ -c $(NVCCFLAGS) -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $< --ptxas-options=-v
	$(NVCC) -M -arch sm_13 -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $< > $(DEPDIR)/$*.Tpo
	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po
	##$(NVCC) -o $<.ptx -ptx --ptxas-options=-v -arch sm_13 -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $<

.cu.lo:
	$(top_srcdir)/cudalt.py $@ $(NVCC) -c $(NVCCFLAGS) -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $< --ptxas-options=-v
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

noinst_HEADERS = pipeline.h simulate_gpu.cu.h sph_polygon.h
galfast_x_SOURCES = footprint.cpp galfast.cpp kernels_cpu.cpp kernels_gpu.cu \
	pipeline.cpp tests.cpp
galfast_x_LDADD = $(top_builddir)/src/skygen/libskygen.la \
	$(top_builddir)/src/common/libcommon.a $(top_builddir)/src/gpc/libgpc.a \
	$(top_builddir)/src/modules/libmodules.la
