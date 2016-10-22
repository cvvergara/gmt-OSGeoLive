#!/bin/bash
#	$Id: GMT_App_F_stand+.sh 16386 2016-05-07 08:53:41Z pwessel $
#
#	Makes the octal code charts in Appendix F
#
# Use the row, col values to generate the octal code needed and
# plot it with gmt pstext, including the header row and left column

bash "${src:-.}"/func_F_stand+.sh > GMT_App_F_stand+.ps
