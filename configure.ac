dnl Copyright (C) 2009, 2010, 2011, 2012, 2013
dnl Computer Graphics Group, University of Siegen, Germany.
dnl Written by Martin Lambers <martin.lambers@uni-siegen.de>.
dnl See http://www.cg.informatik.uni-siegen.de/ for contact information.
dnl
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program.  If not, see <http://www.gnu.org/licenses/>.

dnl Autotools init stuff
AC_PREREQ([2.65])
AC_INIT([ecmdb], [0.7.0], [martin.lambers@uni-siegen.de], [ecmdb], [http://ecm-planet.sourceforge.net/])
AC_CONFIG_SRCDIR([src/ecmdb/main.cpp])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AC_CANONICAL_TARGET
AC_USE_SYSTEM_EXTENSIONS
AC_SYS_LARGEFILE
AM_INIT_AUTOMAKE([1.11.1 color-tests silent-rules tar-pax no-texinfo.tex -Wall -Wno-portability])
AM_SILENT_RULES([yes])
LT_PREREQ([2.2.6])
LT_INIT([win32-dll])
AC_PROG_CXX
AC_PROG_INSTALL
AC_LANG([C++])

dnl We want C++0x, and blindly add -std=gnu++0x. This only works for gcc (and possibly clang)!
CXXFLAGS="-std=gnu++0x $CXXFLAGS"

dnl Useful definitions
AC_DEFINE_UNQUOTED([PLATFORM], ["${target}"], [Platform triplet])

dnl Special handling of W32
case "${target}" in *-*-mingw32*) w32="yes" ;; *) w32="no " ;; esac
if test "$w32" = "yes"; then
    LDFLAGS="-mwindows $LDFLAGS"
    AC_CHECK_TOOL([WINDRES], [windres], [:])
    if test "$WINDRES" = ":"; then
        AC_MSG_ERROR([windres not found])
    fi
    AC_CHECK_TOOL([MAKENSIS], [makensis], [:])
    if test "$MAKENSIS" = ":"; then
        AC_MSG_WARN([makensis not found: you will not be able to create a W32 installer])
    fi
else
    MAKENSIS=":"
fi
AM_CONDITIONAL([HAVE_MAKENSIS], [test "$MAKENSIS" != ":"])
AM_CONDITIONAL([W32], [test "$w32" = "yes"])
if test "$w32" = "yes"; then W32=1; else W32=0; fi
AC_DEFINE_UNQUOTED([W32], [$W32], [Use W32 API?])

dnl Configure switches
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug], [Enable debugging checks. Disabled by default.])],
    [if test "$enableval" = "yes"; then ndebug=0; else ndebug=1; fi], [ndebug=1])
if test "$ndebug" = "1"; then
    AC_DEFINE_UNQUOTED([NDEBUG], [1], [Define to 1 to disable assertions and other debugging checks.])
fi

dnl Feature checks needed by the base modules
AC_LANG_PUSH([C])
dnl - dbg
AC_CHECK_FUNCS([backtrace sigaction])
dnl - fio
case "${target}" in *-*-mingw*) LIBS="$LIBS -lshlwapi" ;; esac
AC_CHECK_FUNCS([fdatasync fnmatch fseeko ftello getpwuid link mmap posix_fadvise readdir_r symlink])
dnl - opt
case "${target}" in *-*-mingw*) CPPFLAGS="$CPPFLAGS -D_BSD_SOURCE" ;; esac
AC_CHECK_DECLS([optreset], [], [], [#include <getopt.h>])
dnl - str
AM_ICONV([])
AC_CHECK_FUNCS([nl_langinfo vasprintf wcswidth])
dnl - sys
AC_CHECK_FUNCS([nanosleep sysconf sched_yield])
dnl - thread
AC_MSG_CHECKING([for GCC atomic builtins])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]],
    [[int b, x = 0; b = __sync_bool_compare_and_swap(&x, 0, 1);]])],
    HAVE_GCC_ATOMIC_BUILTINS=1, HAVE_GCC_ATOMIC_BUILTINS=0)
if test "$HAVE_GCC_ATOMIC_BUILTINS" = "1"; then
    AC_MSG_RESULT([yes])
else
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([GCC atomic builtins not found. Please use a newer compiler.])
fi
AC_SEARCH_LIBS([pthread_create], [pthread])
CPPFLAGS_bak="$CPPFLAGS"
LDFLAGS_bak="$LDFLAGS"
CPPFLAGS="-pthread $CPPFLAGS"
LDFLAGS="-pthread $LDFLAGS"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])], [], [CPPFLAGS="$CPPFLAGS_bak"; LDFLAGS="$LDFLAGS_bak"])
dnl - timer
AC_SEARCH_LIBS([clock_gettime], [rt])
AC_CHECK_FUNCS([clock_gettime])
AC_LANG_POP([C])

dnl pkg-config (required to detect libraries)
PKG_PROG_PKG_CONFIG([])
if test -z "$PKG_CONFIG"; then
    AC_MSG_WARN([pkg-config not found])
    AC_MSG_WARN([required libraries cannot be detected])
    AC_MSG_WARN([please install pkg-config; Debian package: pkg-config])
fi

dnl libecmdb
PKG_CHECK_MODULES([libecmdb], [ecmdb >= 0.0], [HAVE_LIBECMDB=1], [HAVE_LIBECMDB=0])
if test "$HAVE_LIBECMDB" != "1"; then
    AC_MSG_WARN([library libecmdb not found:])
    AC_MSG_WARN([$libecmdb_PKG_ERRORS])
fi

dnl UUID
if test "$w32" = "yes"; then
    libuuid_CFLAGS=""
    libuuid_LIBS="-lrpcrt4"
else
    PKG_CHECK_MODULES([libuuid], [uuid >= 0.0], [HAVE_LIBUUID=1], [HAVE_LIBUUID=0])
    if test "$HAVE_LIBUUID" != "1"; then
        AC_MSG_WARN([library libuuid not found:])
        AC_MSG_WARN([$libuuid_PKG_ERRORS])
        AC_MSG_WARN([libuuid is provided by util-linux-ng; Debian package: uuid-dev])
    fi
    AC_DEFINE_UNQUOTED([HAVE_LIBUUID], [$HAVE_LIBUUID], [Have libuuid?])
fi

dnl GDAL
AC_ARG_VAR([GDALCONFIG], [gdal-config command])
if test -z "$GDALCONFIG"; then
    AC_CHECK_TOOLS([GDALCONFIG], [gdal-config], [])
fi
if test -n "$GDALCONFIG"; then
    libgdal_CFLAGS="`$GDALCONFIG --cflags`"
    libgdal_LIBS="`$GDALCONFIG --libs`"
    AC_SUBST([libgdal_CFLAGS])
    AC_SUBST([libgdal_LIBS])
    HAVE_LIBGDAL="1"
else
    AC_MSG_WARN([GDAL library not found.])
    HAVE_LIBGDAL="0"
fi

dnl Check if all libraries were found
if test -z "$PKG_CONFIG" \
    -o "$HAVE_LIBECMDB" != "1" \
    -o \( "$w32" != "yes" -a "$HAVE_LIBUUID" != "1" \) \
    -o "$HAVE_LIBGDAL" != "1" ; then
    AC_MSG_ERROR([One or more libraries were not found. See messages above.])
fi

dnl Output
AC_CONFIG_FILES([Makefile \
	doc/Makefile \
	src/base/Makefile \
	src/uuid/Makefile \
	src/ecmdb/Makefile \
	src/Makefile
	])
AC_OUTPUT
