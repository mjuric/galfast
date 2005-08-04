bin_PROGRAMS = galaxy unique.x binning.x tobinary.x schlegel.x pm.x dm raytrace.x bin_run.x bin_volume.x merge_maps.x mlcolor.x selector.x model.x explore_maps.x merge_volumes.x fitscatread.x simulate.x locusfit.x locus_demo.x fitter.x
galaxy_SOURCES = galaxy.cpp

# set the include path found by configure
INCLUDES = -I$(top_srcdir)/src/common -I/usr/include/sqlplus/ -I/usr/include/mysql/ $(all_includes)

# the library search path.
galaxy_LDFLAGS = $(all_libraries) 
SUBDIRS = common
unique_x_SOURCES =  unique.cpp
unique_x_LDADD = $(top_builddir)/src/common/libcommon.a
binning_x_SOURCES =  binning.cpp
binning_x_LDADD = $(top_builddir)/src/common/libcommon.a
tobinary_x_SOURCES =  tobinary.cpp
tobinary_x_LDADD = $(top_builddir)/src/common/libcommon.a
schlegel_x_SOURCES =  schlegel.cpp
schlegel_x_LDADD = $(top_builddir)/src/common/libcommon.a
pm_x_SOURCES = pm.cpp
pm_x_LDADD = $(top_builddir)/src/common/libcommon.a
dm_SOURCES = dm.cpp
dm_LDADD = $(top_builddir)/src/common/libcommon.a


raytrace_x_LDADD = $(top_builddir)/src/common/libcommon.a

bin_run_x_SOURCES =  bin_run.cpp
bin_run_x_LDADD = $(top_builddir)/src/common/libcommon.a
bin_volume_x_SOURCES =  bin_volume.cpp
bin_volume_x_LDADD = $(top_builddir)/src/common/libcommon.a
raytrace_x_SOURCES =  raytrace_main.cpp
merge_maps_x_LDADD = $(top_builddir)/src/common/libcommon.a
merge_maps_x_SOURCES =  merge_maps.cpp
mlcolor_x_SOURCES =  mlcolor.cpp
selector_x_SOURCES =  selector.cpp
selector_x_LDADD = $(top_builddir)/src/common/libcommon.a
model_x_SOURCES =  model.cpp
explore_maps_x_SOURCES =  explore_maps.cpp
explore_maps_x_LDADD = $(top_builddir)/src/common/libcommon.a
noinst_HEADERS =  phObjc.h remap.h vrml.h tests.cpp vrml.cpp
merge_volumes_x_SOURCES =  merge_volumes.cpp
merge_volumes_x_LDADD = $(top_builddir)/src/common/libcommon.a
fitscatread_x_SOURCES =  fitscatread.cpp
fitscatread_x_LDADD = $(top_builddir)/src/common/libcommon.a
mlcolor_x_LDADD = $(top_builddir)/src/common/libcommon.a
simulate_x_SOURCES = simulate.cpp
simulate_x_LDADD = $(top_builddir)/src/common/libcommon.a
locusfit_x_SOURCES =  locusfit.cpp
locusfit_x_LDADD = $(top_builddir)/src/common/libcommon.a
locus_demo_x_SOURCES =  locus_demo.cpp
fitter_x_SOURCES =  fitter.cpp
fitter_x_LDADD = $(top_builddir)/src/common/libcommon.a
