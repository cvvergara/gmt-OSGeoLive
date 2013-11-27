#!/bin/bash
#
# $Id: bswap32.sh 11490 2013-05-16 06:26:21Z pwessel $

# test byteswapping 32 bit integers

# generate random binary file
BYTES=$((1024*1024)) # 1 MiB
head -c $BYTES /dev/urandom > data.b

# swap
gmt gmtconvert -bi1I -bo1Iw data.b | gmt gmtconvert -bi1Iw -bo1I > out.u
gmt gmtconvert -bi1i -bo1iw data.b | gmt gmtconvert -bi1iw -bo1i > out.d

# compare result
diff -q data.b out.u
diff -q data.b out.d

# swap using gmt xyz2grd
gmt xyz2grd -Sswapped_tmp.b -ZI data.b
gmt xyz2grd -Sswapped.b -ZI swapped_tmp.b

# compare result
diff -q data.b swapped.b

# files must differ
if diff -q swapped_tmp.b swapped.b; then
  false # return with error
fi
