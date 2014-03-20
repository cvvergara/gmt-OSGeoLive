/*--------------------------------------------------------------------
 *	$Id: gmtlegs.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *    Copyright (c) 1991-2014 by P. Wessel and W. H. F. Smith
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * gmtlegs - A program that locates all legs within a given
 * area, or only those legs that have certain types of data.
 * The program differs from LOC on VAX in that it checks each
 * 1 by 1 degree block to see if the leg actually crossed it.
 *
 * Author:	Paul Wessel
 * Date:	25-MAY-1987
 * Revised:	13-FEB-1991
 * Version 2.0: 21-Jun-1991
 * Version 3.2: 10-MAR-1999
 */

#include "gmt.h"

#define MSIZE 64000

struct INFO {
	char legname[16];
	int leg_id;
	int gmt;
	GMT_LONG mark;
	struct INFO *next_info;
};

int main (int argc, char **argv)
{
	int i, ninfo, n_found = 0, bin, leg, no, flagin, nlegs, gmt, gmtflag = 0;
	size_t not_used = 0;
	char lname[16], line[BUFSIZ], file[BUFSIZ];
	double west, east, south, north;
	FILE *fleg = NULL, *fbin = NULL;
	GMT_LONG error = FALSE, full = FALSE, verbose = FALSE, inside(double w, double e, double s, double n, int binnr);
	struct INFO *leginfo[MSIZE], *ptr[MSIZE], *make_info(char *name, int id_no, int flag);
	static char *data[8] = { "-", "G", "M", "GM", "T", "GT", "MT", "GMT" };

	argc = (int)GMT_begin (argc, argv);

	west = east = south = north = 0.0;
  	gmtflag = 0;

	for (i = 1; i < argc; i++) {
  		if (argv[i][0] == '-') {
  			switch(argv[i][1]) {
  				case 'G':
  					gmtflag += 1;
  					break;
  				case 'L':
  					full = TRUE;
  					break;
  				case 'M':
  					gmtflag += 2;
  					break;
				case 'R':
					sscanf (&argv[i][2], "%lf/%lf/%lf/%lf", &west, &east, &south, &north);
					break;
  				case 'T':
  					gmtflag += 4;
  					break;
				case 'V':
					verbose = TRUE;
					break;
  				default:
  					error = TRUE;
  					break;
  			}
  		}
  		else
  			error = TRUE;
  	}
  	while (west > east) west -= 360.0;
  	if (south >= north) error = TRUE;
  	if (west >= east) error = TRUE;
  	if (argc == 1 || error) {
  		fprintf(stderr,"usage: gmtlegs -R<west>/<east>/<south>/<north> [-G] [-L] [-M] [-T] [-V]\n");
  		fprintf(stderr,"	-R <west>, <east>, <south>, <north> is region in degrees.\n");
  		fprintf(stderr,"	-G (gravity) -M (magnetics) -T (topography). Choose any combination.\n");
  		fprintf(stderr,"	   Legs having at least the specified data are returned. [Default is any data]\n");
  		fprintf(stderr, "	-L gives a long listing (legname and datatypes). [Default is legnames only]\n");
  		fprintf(stderr, "	-V (verbose) prints number of legs found\n");
  		exit (EXIT_FAILURE);
  	}

  	GMT_getsharepath ("mgg", "gmt_legs", ".d", file);
  	if ((fleg = fopen (file, "r")) == NULL) {
		fprintf (stderr,"gmtlegs: Could not open %s\n", file);
		exit (EXIT_FAILURE);
	}
  	GMT_getsharepath ("mgg", "gmt_index", ".b", file);
	if ((fbin = fopen (file, "rb")) == NULL) {
		fprintf(stderr,"gmtlegs: Could not open %s\n", file);
		exit (EXIT_FAILURE);
	}

	/* Make sure that WESN is in whole degrees, truncate if necessary */

	west = floor (west);
	east = ceil (east);
	south = floor (south);
	north = ceil (north);

	/* Read info about each leg */

	ninfo = 0;
	while (fgets (line, BUFSIZ, fleg)) {
		sscanf (line, "%s %d %d", lname, &no, &flagin);
		ptr[no] = leginfo[ninfo++] = make_info (lname, no, flagin);
	}
	fclose (fleg);

	/* Start reading the index-file */

	while (fread ((void *)&bin, (size_t)4, (size_t)1, fbin) != 0) {
		not_used = fread ((void *)&nlegs, (size_t)4, (size_t)1, fbin);
		if (inside (west, east, south, north, bin)) {
			for (i = 0; i < nlegs; i++) {
				not_used = fread ((void *)&leg, (size_t)4, (size_t)1, fbin);
				gmt = leg & 15;
				leg >>= 4;
				if ((gmt & gmtflag) == gmtflag) ptr[leg]->mark = TRUE;
			}
		}
		else
			fseek(fbin, (long)(nlegs*4), 1);
	}
	fclose(fbin);

	for (i = 0; i < ninfo; i++) {
		if (leginfo[i]->mark) {
			(full) ? (printf("%s\t%s\n",leginfo[i]->legname,
					  data[leginfo[i]->gmt])) :
				 (printf("%s\n", leginfo[i]->legname));
			n_found++;
		}
	}
	if (verbose) fprintf (stderr, "gmtlegs: found %d legs\n", n_found);
	exit (EXIT_SUCCESS);
}

struct INFO *make_info(char *name, int id_no, int flag)
{
	struct INFO *new;
	new = (struct INFO *) GMT_memory (VNULL, (size_t)1, sizeof (struct INFO), "gmtlegs");
	strcpy (new->legname,name);
	new->leg_id = id_no;
	new->gmt = flag;
	new->mark = FALSE;
	new->next_info = 0;
	return (new);
}

GMT_LONG inside (double w, double e, double s, double n, int binnr)
{
	double lat, lon;
	lon = (double) (binnr % 360);
	while (lon > e) lon -= 360.0;
	if (lon < w || lon >= e) return (FALSE);
	lat = (double)(binnr / 360) - 90.0;
	if (lat < s || lat >= n) return (FALSE);
	return (TRUE);
}
