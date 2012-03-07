/*--------------------------------------------------------------------
 *	$Id: dat2gmt.c,v 1.22 2011/07/11 19:22:04 guru Exp $
 *
 *    Copyright (c) 1991-2011 by P. Wessel and W. H. F. Smith
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * dat2gmt will read an *.dat file (created by gmt2dat), decode it
 * and create a new *.gmt file. To edit a *.gmt file, run gmt2dat
 * to make an ASCII listing, use the editor to change the contents,
 * and then run dat2gmt to create a new binary *.gmtfile.
 *
 * Author:	Paul Wessel
 * Date:	8-FEB-1988
 * Version:	2.0
 * Version:	3.2 12-APR-1999 (DOS compatible)
 *
 */

#include "gmt.h"
#include "gmt_mgg.h"

#define MEM_SIZE 65536

int main (int argc, char **argv)
{
	int leg_year, n_records, rec, time, mag, top, yr, mo, dy, hr, mi, sc, n_alloc = MEM_SIZE, conv;
	double lat, lon, grv;
	char line[BUFSIZ], agency[10], stime[20], sgrv[10], smag[10], stop[10], *not_used = NULL;
	FILE *fpi = NULL, *fpo = NULL;
	struct GMTMGG_TIME *gmt;
	struct GMTMGG_REC *record = NULL;

	argc = (int)GMT_begin (argc, argv);
	
	if (argc != 3 || ((fpi = fopen(argv[1], "r")) == NULL) || ((fpo = fopen(argv[2], "wb")) == NULL)) {
		fprintf(stderr, "usage: dat2gmt asciifile file.gmt\n");
		exit (EXIT_FAILURE);
	}
	
	record = (struct GMTMGG_REC *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (struct GMTMGG_REC), "dat2gmt");
	
	/* Read first record of file containing start-year, n_records and info */
	
	not_used = fgets (line, BUFSIZ, fpi);
	sscanf (line, "%d %d %s", &leg_year, &n_records, agency);
	gmt = gmtmgg_init (leg_year);
	
	rec = 0;
	while (fgets (line, BUFSIZ, fpi)) {
		sscanf (line, "%s %lf %lf %s %s %s", stime, &lat, &lon, sgrv, smag, stop);
		sscanf (stime, "%d/%d/%d/%d:%d:%d", &yr, &mo, &dy, &hr, &mi, &sc);
		record[rec].time = gmtmgg_time (&time, yr, mo, dy, hr, mi, sc, gmt);
		record[rec].lat = (int) (lat * 1000000.0);
		record[rec].lon = (int) (lon * 1000000.0);
		conv = sscanf(sgrv, "%lf", &grv);
		record[rec].gmt[0] = (conv > 0 && !GMT_is_dnan (grv)) ? (int) rint(grv * 10.0) : GMTMGG_NODATA;
		conv = sscanf(smag, "%d", &mag);
		record[rec].gmt[1] = (conv > 0) ? mag : GMTMGG_NODATA;
		conv = sscanf(stop, "%d", &top);
		record[rec].gmt[2] = (conv > 0) ? top : GMTMGG_NODATA;
		rec++;
		if (rec == n_alloc) {
			n_alloc <<= 1;
			record = (struct GMTMGG_REC *) GMT_memory((void *)record, (size_t)n_alloc, sizeof (struct GMTMGG_REC), "dat2gmt");
		}
	}
	fclose (fpi);
	
	GMT_free ((void *)gmt);
	n_records = rec;
	if (fwrite((void *)(&leg_year), (size_t)4, (size_t)1, fpo) != 1) {
		fprintf(stderr,"dat2gmt: Error while writing first year\n");
		exit (EXIT_FAILURE);
	}
	if (fwrite((void *)(&n_records), (size_t)4, (size_t)1, fpo) != 1) {
		fprintf(stderr,"dat2gmt: Error while writing no of records\n");
		exit (EXIT_FAILURE);
	}
	if (fwrite((void *)agency, (size_t)10, (size_t)1, fpo) != 1) {
		fprintf(stderr,"dat2gmt: Error while writing info-header\n");
		exit (EXIT_FAILURE);
	}
	for (rec = 0; rec < n_records; rec++) {
		if (fwrite ((void *)(&record[rec]), (size_t)18, (size_t)1, fpo) != 1) {
			fprintf (stderr,"dat2gmt: Error writing data record no %d\n",rec);
			exit (EXIT_FAILURE);
		}
	}
	fclose (fpo);
	
	GMT_free ((void *)record);

	GMT_end (argc, argv);
	
	exit (EXIT_SUCCESS);
}
