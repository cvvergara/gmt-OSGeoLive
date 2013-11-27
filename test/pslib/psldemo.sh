#!/bin/bash
#
# $Id: psldemo.sh 12114 2013-09-03 19:19:00Z fwobbe $
#
# Purpose:      Test all PSL functions at least once
# GMT progs:    libpslib, psldemo
# Unix progs:   -
#
ps=psldemo.ps
export PSL_SHAREDIR="$GMT_SHAREDIR"
psldemo > $ps
