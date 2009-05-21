bin_PROGRAMS = simulate.x
EXTRA_DIST = ri_prior.txt

## SUFFIXES = .cu
if HAVE_CUDA
.cu.o:
	$(NVCC) --host-compilation 'C++' -o $<.ptx -ptx --ptxas-options=-v -arch sm_13 -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $<
	$(NVCC) -o $@ -c $(NVCCFLAGS) -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $<

.cu.lo:
	$(NVCC) --host-compilation 'C++' -o $<.ptx -ptx --ptxas-options=-v -arch sm_13 -I$(top_builddir) $(INCLUDES) $(CPPFLAGS) $<
	$(top_srcdir)/cudalt.py $@ $(NVCC) -c $(NVCCFLAGS) -I$(top_builddir) $(INCLUDES) $(CPPFLAGS)  $<
else
.cu.o:
	touch $@
endif

# set the include path found by configure
INCLUDES = -I$(top_srcdir)/src/common -I/usr/include/sqlplus/ -I/usr/include/mysql/ $(all_includes)

# the library search path.
SUBDIRS = common gpc



noinst_HEADERS = simulate.h simulate_base.h
simulate_x_SOURCES = observe.cpp observe_FeH.cpp observe_kinTMIII.cpp \
	observe_vel2pm.cpp simulate.cpp simulate.h simulate2.cpp simulate_cpu.cpp simulate_gpu.cu
simulate_x_LDADD = $(top_builddir)/src/gpc/libgpc.a \
	$(top_builddir)/src/common/libcommon.a
