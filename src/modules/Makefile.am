INCLUDES = $(all_includes) -I$(top_srcdir)/src/common
METASOURCES = AUTO
lib_LTLIBRARIES = libmodules.la
libmodules_la_SOURCES = FeH_gpu.cu.h FeH_host.cpp fixedFeH_gpu.cu.h \
	fixedFeH_host.cpp unresolvedMultiples_gpu.cu.h unresolvedMultiples_host.cpp
libmodules_la_LDFLAGS = -module
