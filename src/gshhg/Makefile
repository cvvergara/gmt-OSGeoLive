# $Id: Makefile 9817 2012-04-29 22:15:19Z pwessel $ # # makefile for gshhg directory
# Does not depend on GMT src or libraries
#
#-------------------------------------------------------------------------------
# !! STOP EDITING HERE, THE REST IS FIXED !!
#-------------------------------------------------------------------------------

GMTSRCDIR = ../
include $(GMTSRCDIR)config.mk
include $(GMTSRCDIR)common.mk

PROGS_O		= gshhg.o gshhg_dp.o gshhgtograss.o
MAN1		= $(PROGS_O:.o=.1)

#-------------------------------------------------------------------------------
#	software targets
#-------------------------------------------------------------------------------

all:		$(PROGS)

install:	all
		$(INSTALL) $(PROGS) $(bindir)

uninstall:
		cd $(bindir); \rm -f $(PROGS)

install-man:	$(MAN1)
		$(INSTALL) -m 644 $(MAN1) $(mandir)/man1

uninstall-man:
		cd $(mandir)/man1; \rm -f $(MAN1)

spotless::	clean

clean:	
		\rm -f *.o *% $(PROGS)

#-------------------------------------------------------------------------------
#	object file dependencies
#-------------------------------------------------------------------------------

$(PROGS_O):	gshhg.h

#-------------------------------------------------------------------------------
#	program rules
#-------------------------------------------------------------------------------

gshhg$(EXE):		gshhg.o
gshhg_dp$(EXE):		gshhg_dp.o
gshhgtograss$(EXE):	gshhgtograss.o

$(PROGS):
		$(CC) $(LDFLAGS) $(@:.exe=).o $(NETCDF_LIB) $(GDAL_LIB) $(LIBS) -o $@
		$(COMPRESS) $@
