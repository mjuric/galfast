INCLUDES = -I$(top_srcdir)/src/common -I$(top_srcdir)/src/gpc \
	-I$(top_srcdir)/src/modules $(all_includes)
METASOURCES = AUTO
lib_LTLIBRARIES = libskygen.la
libskygen_la_LDFLAGS = -module
libskygen_la_SOURCES = model_J08.cpp model_densityCube.cpp \
	model_expDisk.cpp model_powerLawEllipsoid.cpp \
	os_skygen.cpp skyconfig_impl.h skygen.cu.h
noinst_HEADERS = model_J08.h model_densityCube.h model_expDisk.h \
	model_powerLawEllipsoid.h skygen.h model_lib.h
