# Copyright (C) 2010, 2011, 2012
# Martin Lambers <martin.lambers@uni-siegen.de>
#
# Copying and distribution of this file, with or without modification, are
# permitted in any medium without royalty provided the copyright notice and this
# notice are preserved. This file is offered as-is, without any warranty.

man_MANS = ecmdb.1

info_TEXINFOS = ecmdb.texi

IMAGES = \
	earth-bmng.jpg \
	earth-natearth.jpg \
	earth-srtm4.jpg \
	earth-terrasarx.jpg \
	earth-utah.jpg \
	moon-wacdtm.jpg \
	mars-moc.jpg \
	mars-mola.jpg \
	mars-hirise.jpg

ecmdb_TEXINFOS = ecmdb.texi $(IMAGES)

AM_MAKEINFOHTMLFLAGS = --no-split

dist_html_DATA = ecmdb.html $(IMAGES)

EXTRA_DIST = $(man_MANS)
