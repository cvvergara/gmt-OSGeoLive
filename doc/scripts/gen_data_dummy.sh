#!/bin/bash
#
# $Id: gen_data_dummy.sh 11490 2013-05-16 06:26:21Z pwessel $
# This script makes the dummy data sets needed in GMT_linear.sh GMT_log.sh GMT_pow.sh
#
gmt gmtmath -T0/100/1  T SQRT = sqrt.d
gmt gmtmath -T0/100/10 T SQRT = sqrt.d10
