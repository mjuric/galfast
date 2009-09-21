INCLUDES = $(all_includes) -I$(top_srcdir)/src/common
METASOURCES = AUTO
lib_LTLIBRARIES = libmodules.la
libmodules_la_SOURCES = FeH_gpu.cu.h FeH_host.cpp fixedFeH_gpu.cu.h \
	fixedFeH_host.cpp gal2other_gpu.cu.h gal2other_host.cpp kinTMIII_gpu.cu.h \
	kinTMIII_host.cpp photometricErrors_gpu.cu.h photometricErrors_host.cpp \
	photometry_gpu.cu.h photometry_host.cpp unresolvedMultiples_gpu.cu.h \
	unresolvedMultiples_host.cpp vel2pm_gpu.cu.h vel2pm_host.cpp
libmodules_la_LDFLAGS = -module
noinst_HEADERS = module_lib.h
