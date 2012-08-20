# Copyright (C) 2011, 2012
# Martin Lambers <martin.lambers@uni-siegen.de>
#
# Copying and distribution of this file, with or without modification, are
# permitted in any medium without royalty provided the copyright notice and this
# notice are preserved. This file is offered as-is, without any warranty.

bin_PROGRAMS = ecmdb

ecmdb_SOURCES = \
	quadlist.h quadlist.cpp \
	srgb.h srgb.cpp \
	metadata.h metadata.cpp \
        compression-info.h compression-info.cpp \
	cmd_add.cpp \
	cmd_commit.cpp \
	cmd_create.cpp \
	cmd_filename.cpp \
	cmd_finalize.cpp \
	cmd_resolutions.cpp \
	cmds.h cmds.cpp \
	main.cpp

ecmdb_CPPFLAGS = \
	-I$(top_srcdir)/src/base \
	-I$(top_srcdir)/src/uuid \
	$(libecm_CFLAGS) \
	$(libgdal_CFLAGS) \
	-fopenmp

ecmdb_LDADD = \
	../uuid/libuuid.la \
	$(libuuid_LIBS) \
	../base/libbase.la \
	$(libecm_LIBS) \
	$(libgdal_LIBS) \
	$(LIBICONV)

ecmdb_LDFLAGS = -fopenmp
