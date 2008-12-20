dnl RS_BOOST([MINIMUM-VERSION], [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl Test for the Boost C++ libraries of a particular version (or newer)
dnl Defines:
dnl   BOOST_CPPFLAGS to the set of flags required to compiled Boost
AC_DEFUN([RS_BOOST], 
[
  AC_SUBST(BOOST_CPPFLAGS)
  AC_SUBST(BOOST_LIBS)
  BOOST_CPPFLAGS=""
  path_given="no"

dnl Extract the path name from a --with-boost=PATH argument
  AC_ARG_WITH(boost,
    AC_HELP_STRING([--with-boost=PATH],[absolute path name where the Boost C++ libraries reside, or `int', for internal library. Alternatively, the BOOST_ROOT environment variable will be used]),
    [
    if test "$withval" != yes ; then
        BOOST_ROOT="$withval"
    fi
    ]
  )

  if test "x$BOOST_ROOT" = x ; then
    BOOST_ROOT="/usr/local"
  fi

  boost_min_version=ifelse([$1], ,1.20.0,[$1])
  WANT_BOOST_MAJOR=`expr $boost_min_version : '\([[0-9]]\+\)'`
  WANT_BOOST_MINOR=`expr $boost_min_version : '[[0-9]]\+\.\([[0-9]]\+\)'`
  WANT_BOOST_SUB_MINOR=`expr $boost_min_version : '[[0-9]]\+\.[[0-9]]\+\.\([[0-9]]\+\)'`

  BOOST_CPPFLAGS="-I$BOOST_ROOT/include/boost-${WANT_BOOST_MAJOR}_$WANT_BOOST_MINOR"

  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  OLD_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
  AC_MSG_CHECKING([for the Boost C++ libraries, version $boost_min_version or newer])
  AC_TRY_COMPILE(
    [
#include <boost/version.hpp>
],
    [],
    [
      have_boost="yes"
    ],
    [
      AC_MSG_RESULT(no)
      have_boost="no"
      ifelse([$3], , :, [$3])
    ])

  if test "$have_boost" = "yes"; then
    WANT_BOOST_VERSION=`expr $WANT_BOOST_MAJOR \* 100000 \+ $WANT_BOOST_MINOR \* 100 \+ $WANT_BOOST_SUB_MINOR`

    AC_TRY_COMPILE(
      [
#include <boost/version.hpp>
],
      [
#if BOOST_VERSION >= $WANT_BOOST_VERSION
// Everything is okay
#else
#  error Boost version is too old
#endif

],
      [
	AC_MSG_RESULT(yes)
dnl	if test "$target_os" = "mingw32"; then
dnl	   boost_libsuff=mgw
dnl	else
dnl	   boost_libsuff=gcc
	   boost_libsuff=
dnl	fi
	boost_libsuff_r=$boost_libsuff-mt;
dnl	if test "x$enable_debug" = xyes ; then
dnl	    boost_libsuff=$boost_libsuff-d;
dnl	    boost_libsuff_r=$boost_libsuff_r-d;
dnl	fi
dnl	boost_libsuff=$boost_libsuff-${WANT_BOOST_MAJOR}_$WANT_BOOST_MINOR
dnl	boost_libsuff_r=$boost_libsuff_r-${WANT_BOOST_MAJOR}_$WANT_BOOST_MINOR
dnl        AC_MSG_RESULT([Suffix: $boost_libsuff])
        ifelse([$2], , :, [$2])
      ],
      [
        AC_MSG_RESULT([no, version of installed Boost libraries is too old])
        ifelse([$3], , :, [$3])
      ])
  fi
  CPPFLAGS="$OLD_CPPFLAGS"
  AC_LANG_RESTORE
])


dnl RS_BOOST_THREAD([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl Test for the Boost thread library
dnl Defines
dnl   BOOST_LDFLAGS to the set of flags required to compile boost_thread
AC_DEFUN([RS_BOOST_THREAD], 
[
    AC_REQUIRE([RS_BOOST])
  AC_MSG_CHECKING([whether we can use boost_thread library])
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  OLD_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$BOOST_CPPFLAGS -D_REENTRANT"
  OLD_LIBS="$LIBS"
  LIBS="-lboost_thread$boost_libsuff_r"
    AC_TRY_COMPILE(
	    [ 
		#include <boost/thread.hpp> 
		bool bRet = 0;
		void thdfunc() { bRet = 1; }
	    ],
	    [
		boost::thread thrd(&thdfunc);
		thrd.join();
		return bRet == 1;
	    ], 
	    [
		AC_MSG_RESULT([yes])
		ifelse([$1], , :, [$1])
	    ],
	    [ 
		AC_MSG_RESULT([no])
		ifelse([$2], , :, [$2])
	    ])

    AC_SUBST(BOOST_CPPFLAGS)
    AC_SUBST(BOOST_LIBS_R)
    BOOST_CPPFLAGS="$CPPFLAGS"
    BOOST_LIBS_R="$LIBS $BOOST_LIBS_R"
    CPPFLAGS="$OLD_CPPFLAGS"
    LIBS="$OLD_LIBS"
    AC_LANG_RESTORE
])
	    
dnl RS_BOOST_DATETIME([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl Test for the Boost datetime library
dnl Defines
dnl   BOOST_LDFLAGS to the set of flags required to compile boost_datetime
AC_DEFUN([RS_BOOST_DATETIME], 
[
    AC_REQUIRE([RS_BOOST])
  AC_MSG_CHECKING([whether we can use boost_datetime library])
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  OLD_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$BOOST_CPPFLAGS"
  OLD_LIBS="$LIBS"
  LIBS="-lboost_date_time$boost_libsuff"
    AC_TRY_COMPILE(
	    [ 
		#include <boost/date_time/gregorian/gregorian.hpp> 
	    ],
	    [
		using namespace boost::gregorian;
		date d = from_string("1978-01-27");
		return d == boost::gregorian::date(1978, Jan, 27);
	    ], 
	    [
		AC_MSG_RESULT([yes])
		ifelse([$1], , :, [$1])
	    ],
	    [ 
		AC_MSG_RESULT([no])
		ifelse([$2], , :, [$2])
	    ])

    AC_SUBST(BOOST_CPPFLAGS)
    AC_SUBST(BOOST_LIBS)
    AC_SUBST(BOOST_LIBS_R)
    BOOST_CPPFLAGS="$CPPFLAGS"
    BOOST_LIBS="$LIBS $BOOST_LIBS"
    BOOST_LIBS_R="$LIBS $BOOST_LIBS_R"
    CPPFLAGS="$OLD_CPPFLAGS"
    LIBS=$OLD_LIBS
    AC_LANG_RESTORE
])
	    
dnl RS_BOOST_REGEX([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl Test for the Boost regex library
dnl Defines
dnl   BOOST_LDFLAGS to the set of flags required to compile boost_datetime
AC_DEFUN([RS_BOOST_REGEX], 
[
    AC_REQUIRE([RS_BOOST])
  AC_MSG_CHECKING([whether we can use boost_regex library])
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  OLD_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$BOOST_CPPFLAGS"
  OLD_LIBS="$LIBS"
  LIBS="-lboost_regex$boost_libsuff"
    AC_TRY_COMPILE(
	    [ 
		#include <boost/regex.hpp> 
	    ],
	    [
		using namespace boost;
		cmatch what;
		if(!regex_match("27/01/78",what,regex("(\\\d+)/(\\\d+)/(\\\d+)")))
		    return 0;

		return what[1]=="27" && what[2]=="01" && what[3]=="78";
	    ], 
	    [
		AC_MSG_RESULT([yes])
		ifelse([$1], , :, [$1])
	    ],
	    [ 
		AC_MSG_RESULT([no])
		ifelse([$2], , :, [$2])
	    ])

    AC_SUBST(BOOST_CPPFLAGS)
    AC_SUBST(BOOST_LIBS)
    AC_SUBST(BOOST_LIBS_R)
    BOOST_CPPFLAGS="$CPPFLAGS"
    BOOST_LIBS="$LIBS $BOOST_LIBS"
    BOOST_LIBS_R="$LIBS $BOOST_LIBS_R"
    CPPFLAGS="$OLD_CPPFLAGS"
    LIBS="$OLD_LIBS"
    AC_LANG_RESTORE
])

dnl RS_BOOST_PROGRAM_OPTIONS([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl Test for the Boost regex library
dnl Defines
dnl   BOOST_LDFLAGS to the set of flags required to compile boost_datetime
AC_DEFUN([RS_BOOST_PROGRAM_OPTIONS], 
[
    AC_REQUIRE([RS_BOOST])
  AC_MSG_CHECKING([whether we can use boost_program_options library])
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  OLD_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$BOOST_CPPFLAGS"
  OLD_LIBS="$LIBS"
  LIBS="-lboost_program_options$boost_libsuff"
    AC_TRY_COMPILE(
	    [ 
		#include <boost/program_options.hpp> 
	    ],
	    [
		using namespace boost::program_options;
		return 0;
	    ], 
	    [
		AC_MSG_RESULT([yes])
		ifelse([$1], , :, [$1])
	    ],
	    [ 
		AC_MSG_RESULT([no])
		ifelse([$2], , :, [$2])
	    ])

    AC_SUBST(BOOST_CPPFLAGS)
    AC_SUBST(BOOST_LIBS)
    AC_SUBST(BOOST_LIBS_R)
    BOOST_CPPFLAGS="$CPPFLAGS"
    BOOST_LIBS="$LIBS $BOOST_LIBS"
    BOOST_LIBS_R="$LIBS $BOOST_LIBS_R"
    CPPFLAGS="$OLD_CPPFLAGS"
    LIBS="$OLD_LIBS"
    AC_LANG_RESTORE
])

dnl RS_BOOST_IOSTREAMS([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl Test for the Boost iostreams library
dnl Defines
dnl   BOOST_LDFLAGS to the set of flags required to compile boost_iostreams
AC_DEFUN([RS_BOOST_IOSTREAMS], 
[
    AC_REQUIRE([RS_BOOST])
  AC_MSG_CHECKING([whether we can use boost_iostreams library])
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  OLD_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$BOOST_CPPFLAGS"
  OLD_LIBS="$LIBS"
  LIBS="-lboost_iostreams$boost_libsuff"
    AC_TRY_COMPILE(
	    [ 
		#include <boost/iostreams/device/file.hpp> 
	    ],
	    [
		using namespace boost::iostreams;
		file_sink fs("bla.bin");
		return 0;
	    ], 
	    [
		AC_MSG_RESULT([yes])
		ifelse([$1], , :, [$1])
	    ],
	    [ 
		AC_MSG_RESULT([no])
		ifelse([$2], , :, [$2])
	    ])

    AC_SUBST(BOOST_CPPFLAGS)
    AC_SUBST(BOOST_LIBS)
    AC_SUBST(BOOST_LIBS_R)
    BOOST_CPPFLAGS="$CPPFLAGS"
    BOOST_LIBS="$LIBS $BOOST_LIBS"
    BOOST_LIBS_R="$LIBS $BOOST_LIBS_R"
    CPPFLAGS="$OLD_CPPFLAGS"
    LIBS="$OLD_LIBS"
    AC_LANG_RESTORE
])

dnl @synopsis AC_caolan_CHECK_PACKAGE(PACKAGE, FUNCTION, LIBRARY , HEADERFILE [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Provides --with-PACKAGE, --with-PACKAGE-include and --with-PACKAGE-libdir
dnl options to configure. Supports the now standard --with-PACKAGE=DIR 
dnl approach where the package's include dir and lib dir are underneath DIR,
dnl but also allows the include and lib directories to be specified seperately
dnl
dnl adds the extra -Ipath to CFLAGS if needed 
dnl adds extra -Lpath to LD_FLAGS if needed
dnl searches for the FUNCTION in the LIBRARY with 
dnl AC_CHECK_LIBRARY and thus adds the lib to LIBS
dnl
dnl defines HAVE_PKG_PACKAGE if it is found, (where PACKAGE in the 
dnl HAVE_PKG_PACKAGE is replaced with the actual first parameter passed)
dnl note that autoheader will complain of not having the HAVE_PKG_PACKAGE and you 
dnl will have to add it to acconfig.h manually
dnl
dnl @version $Id: acinclude.m4,v 1.1 2005/12/15 09:06:00 mjuric Exp $
dnl @author Caolan McNamara <caolan@skynet.ie>
dnl
dnl with fixes from...
dnl Alexandre Duret-Lutz <duret_g@lrde.epita.fr>

AC_DEFUN(AC_caolan_CHECK_PACKAGE,
[

AC_ARG_WITH($1,
[  --with-$1[=DIR]	root directory of $1 installation],
with_$1=$withval 
if test "${with_$1}" != yes; then
	$1_include="$withval/include" 
	$1_libdir="$withval/lib"
fi
)

AC_ARG_WITH($1-include,
[  --with-$1-include=DIR        specify exact include dir for $1 headers],
$1_include="$withval")

AC_ARG_WITH($1-libdir,
[  --with-$1-libdir=DIR        specify exact library dir for $1 library
  --without-$1        disables $1 usage completely], 
$1_libdir="$withval")

if test "${with_$1}" != no ; then
	OLD_LIBS=$LIBS
	OLD_LDFLAGS=$LDFLAGS
	OLD_CFLAGS=$CFLAGS
	OLD_CPPFLAGS=$CPPFLAGS

	if test "${$1_libdir}" ; then
		LDFLAGS="$LDFLAGS -L${$1_libdir}"
	fi
	if test "${$1_include}" ; then
		CPPFLAGS="$CPPFLAGS -I${$1_include}"
		CFLAGS="$CFLAGS -I${$1_include}"
	fi

	AC_CHECK_LIB($3,$2,,no_good=yes)
	AC_CHECK_HEADER($4,,no_good=yes)
	if test "$no_good" = yes; then
dnl	broken
		ifelse([$6], , , [$6])
		
		LIBS=$OLD_LIBS
		LDFLAGS=$OLD_LDFLAGS
		CPPFLAGS=$OLD_CPPFLAGS
		CFLAGS=$OLD_CFLAGS
	else
dnl	fixed
		ifelse([$5], , , [$5])

		AC_DEFINE(HAVE_PKG_$1, 1, Define if you have the package )
	fi

fi

])


dnl @synopsis ACX_PTHREAD([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl @summary figure out how to build C programs using POSIX threads
dnl
dnl This macro figures out how to build C programs using POSIX threads.
dnl It sets the PTHREAD_LIBS output variable to the threads library and
dnl linker flags, and the PTHREAD_CFLAGS output variable to any special
dnl C compiler flags that are needed. (The user can also force certain
dnl compiler flags/libs to be tested by setting these environment
dnl variables.)
dnl
dnl Also sets PTHREAD_CC to any special C compiler that is needed for
dnl multi-threaded programs (defaults to the value of CC otherwise).
dnl (This is necessary on AIX to use the special cc_r compiler alias.)
dnl
dnl NOTE: You are assumed to not only compile your program with these
dnl flags, but also link it with them as well. e.g. you should link
dnl with $PTHREAD_CC $CFLAGS $PTHREAD_CFLAGS $LDFLAGS ... $PTHREAD_LIBS
dnl $LIBS
dnl
dnl If you are only building threads programs, you may wish to use
dnl these variables in your default LIBS, CFLAGS, and CC:
dnl
dnl        LIBS="$PTHREAD_LIBS $LIBS"
dnl        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
dnl        CC="$PTHREAD_CC"
dnl
dnl In addition, if the PTHREAD_CREATE_JOINABLE thread-attribute
dnl constant has a nonstandard name, defines PTHREAD_CREATE_JOINABLE to
dnl that name (e.g. PTHREAD_CREATE_UNDETACHED on AIX).
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a threads
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands to
dnl run it if it is not found. If ACTION-IF-FOUND is not specified, the
dnl default action will define HAVE_PTHREAD.
dnl
dnl Please let the authors know if this macro fails on any platform, or
dnl if you have any other suggestions or comments. This macro was based
dnl on work by SGJ on autoconf scripts for FFTW (www.fftw.org) (with
dnl help from M. Frigo), as well as ac_pthread and hb_pthread macros
dnl posted by Alejandro Forero Cuervo to the autoconf macro repository.
dnl We are also grateful for the helpful feedback of numerous users.
dnl
dnl @category InstalledPackages
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl @version 2005-06-15
dnl @license GPLWithACException

AC_DEFUN([ACX_PTHREAD], [
AC_REQUIRE([AC_CANONICAL_HOST])
AC_LANG_SAVE
AC_LANG_C
acx_pthread_ok=no

# We used to check for pthread.h first, but this fails if pthread.h
# requires special compiler flags (e.g. on True64 or Sequent).
# It gets checked for in the link test anyway.

# First of all, check if the user has set any of the PTHREAD_LIBS,
# etcetera environment variables, and if threads linking works using
# them:
if test x"$PTHREAD_LIBS$PTHREAD_CFLAGS" != x; then
        save_CFLAGS="$CFLAGS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
        save_LIBS="$LIBS"
        LIBS="$PTHREAD_LIBS $LIBS"
        AC_MSG_CHECKING([for pthread_join in LIBS=$PTHREAD_LIBS with CFLAGS=$PTHREAD_CFLAGS])
        AC_TRY_LINK_FUNC(pthread_join, acx_pthread_ok=yes)
        AC_MSG_RESULT($acx_pthread_ok)
        if test x"$acx_pthread_ok" = xno; then
                PTHREAD_LIBS=""
                PTHREAD_CFLAGS=""
        fi
        LIBS="$save_LIBS"
        CFLAGS="$save_CFLAGS"
fi

# We must check for the threads library under a number of different
# names; the ordering is very important because some systems
# (e.g. DEC) have both -lpthread and -lpthreads, where one of the
# libraries is broken (non-POSIX).

# Create a list of thread flags to try.  Items starting with a "-" are
# C compiler flags, and other items are library names, except for "none"
# which indicates that we try without any flags at all, and "pthread-config"
# which is a program returning the flags for the Pth emulation library.

acx_pthread_flags="pthreads none -Kthread -kthread lthread -pthread -pthreads -mthreads pthread --thread-safe -mt pthread-config"

# The ordering *is* (sometimes) important.  Some notes on the
# individual items follow:

# pthreads: AIX (must check this before -lpthread)
# none: in case threads are in libc; should be tried before -Kthread and
#       other compiler flags to prevent continual compiler warnings
# -Kthread: Sequent (threads in libc, but -Kthread needed for pthread.h)
# -kthread: FreeBSD kernel threads (preferred to -pthread since SMP-able)
# lthread: LinuxThreads port on FreeBSD (also preferred to -pthread)
# -pthread: Linux/gcc (kernel threads), BSD/gcc (userland threads)
# -pthreads: Solaris/gcc
# -mthreads: Mingw32/gcc, Lynx/gcc
# -mt: Sun Workshop C (may only link SunOS threads [-lthread], but it
#      doesn't hurt to check since this sometimes defines pthreads too;
#      also defines -D_REENTRANT)
#      ... -mt is also the pthreads flag for HP/aCC
# pthread: Linux, etcetera
# --thread-safe: KAI C++
# pthread-config: use pthread-config program (for GNU Pth library)

case "${host_cpu}-${host_os}" in
        *solaris*)

        # On Solaris (at least, for some versions), libc contains stubbed
        # (non-functional) versions of the pthreads routines, so link-based
        # tests will erroneously succeed.  (We need to link with -pthreads/-mt/
        # -lpthread.)  (The stubs are missing pthread_cleanup_push, or rather
        # a function called by this macro, so we could check for that, but
        # who knows whether they'll stub that too in a future libc.)  So,
        # we'll just look for -pthreads and -lpthread first:

        acx_pthread_flags="-pthreads pthread -mt -pthread $acx_pthread_flags"
        ;;
esac

if test x"$acx_pthread_ok" = xno; then
for flag in $acx_pthread_flags; do

        case $flag in
                none)
                AC_MSG_CHECKING([whether pthreads work without any flags])
                ;;

                -*)
                AC_MSG_CHECKING([whether pthreads work with $flag])
                PTHREAD_CFLAGS="$flag"
                ;;

		pthread-config)
		AC_CHECK_PROG(acx_pthread_config, pthread-config, yes, no)
		if test x"$acx_pthread_config" = xno; then continue; fi
		PTHREAD_CFLAGS="`pthread-config --cflags`"
		PTHREAD_LIBS="`pthread-config --ldflags` `pthread-config --libs`"
		;;

                *)
                AC_MSG_CHECKING([for the pthreads library -l$flag])
                PTHREAD_LIBS="-l$flag"
                ;;
        esac

        save_LIBS="$LIBS"
        save_CFLAGS="$CFLAGS"
        LIBS="$PTHREAD_LIBS $LIBS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"

        # Check for various functions.  We must include pthread.h,
        # since some functions may be macros.  (On the Sequent, we
        # need a special flag -Kthread to make this header compile.)
        # We check for pthread_join because it is in -lpthread on IRIX
        # while pthread_create is in libc.  We check for pthread_attr_init
        # due to DEC craziness with -lpthreads.  We check for
        # pthread_cleanup_push because it is one of the few pthread
        # functions on Solaris that doesn't have a non-functional libc stub.
        # We try pthread_create on general principles.
        AC_TRY_LINK([#include <pthread.h>],
                    [pthread_t th; pthread_join(th, 0);
                     pthread_attr_init(0); pthread_cleanup_push(0, 0);
                     pthread_create(0,0,0,0); pthread_cleanup_pop(0); ],
                    [acx_pthread_ok=yes])

        LIBS="$save_LIBS"
        CFLAGS="$save_CFLAGS"

        AC_MSG_RESULT($acx_pthread_ok)
        if test "x$acx_pthread_ok" = xyes; then
                break;
        fi

        PTHREAD_LIBS=""
        PTHREAD_CFLAGS=""
done
fi

# Various other checks:
if test "x$acx_pthread_ok" = xyes; then
        save_LIBS="$LIBS"
        LIBS="$PTHREAD_LIBS $LIBS"
        save_CFLAGS="$CFLAGS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"

        # Detect AIX lossage: JOINABLE attribute is called UNDETACHED.
	AC_MSG_CHECKING([for joinable pthread attribute])
	attr_name=unknown
	for attr in PTHREAD_CREATE_JOINABLE PTHREAD_CREATE_UNDETACHED; do
	    AC_TRY_LINK([#include <pthread.h>], [int attr=$attr; return attr;],
                        [attr_name=$attr; break])
	done
        AC_MSG_RESULT($attr_name)
        if test "$attr_name" != PTHREAD_CREATE_JOINABLE; then
            AC_DEFINE_UNQUOTED(PTHREAD_CREATE_JOINABLE, $attr_name,
                               [Define to necessary symbol if this constant
                                uses a non-standard name on your system.])
        fi

        AC_MSG_CHECKING([if more special flags are required for pthreads])
        flag=no
        case "${host_cpu}-${host_os}" in
            *-aix* | *-freebsd* | *-darwin*) flag="-D_THREAD_SAFE";;
            *solaris* | *-osf* | *-hpux*) flag="-D_REENTRANT";;
        esac
        AC_MSG_RESULT(${flag})
        if test "x$flag" != xno; then
            PTHREAD_CFLAGS="$flag $PTHREAD_CFLAGS"
        fi

        LIBS="$save_LIBS"
        CFLAGS="$save_CFLAGS"

        # More AIX lossage: must compile with cc_r
        AC_CHECK_PROG(PTHREAD_CC, cc_r, cc_r, ${CC})
else
        PTHREAD_CC="$CC"
fi

AC_SUBST(PTHREAD_LIBS)
AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_CC)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_pthread_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_PTHREAD,1,[Define if you have POSIX threads libraries and header files.]),[$1])
        :
else
        acx_pthread_ok=no
        $2
fi
AC_LANG_RESTORE
])dnl ACX_PTHREAD
