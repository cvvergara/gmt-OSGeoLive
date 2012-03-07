/*--------------------------------------------------------------------
 *	$Id: gmt2dat.c,v 1.21 2011/07/11 19:22:05 guru Exp $
 *
 *    Copyright (c) 1991-2011 by P. Wessel and W. H. F. Smith
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * gmt2dat will read a gmt file and create an ascii listing
 * of all the information stored in the *.gmt file. The ASCII
 * file can then be edited manually. Run dat2gmt to put the
 * edited file into *.gmt binary format.
 *
 * Author:	Paul Wessel
 * Date:	8-FEB-1988
 * Version:	3.2, 12-APR-1999 (DOS compatible)
 *
 */

#include "gmt.h"
#include "gmt_mgg.h"

int main (int argc, char **argv)
{
	int n_records, leg_year;
	int rec, yr, mo, dy, hr, mi, sc;
	double lat, lon;
	char agency[10];
	FILE *fp = NULL;
	struct GMTMGG_TIME *gmt = NULL;
	struct GMTMGG_REC *record = NULL;

	argc = (int)GMT_begin (argc, argv);

	if (argc != 2 || (fp = fopen(argv[1], "rb")) == NULL) {
		fprintf(stderr, "usage: gmt2dat file.gmt > asciifile\n");
		exit (EXIT_FAILURE);
	}
	
	/* Read first record of file containing start-year, n_records and info */
	
	if (fread((void *)(&leg_year), (size_t)4, (size_t)1, fp) != 1) {
		fprintf(stderr,"gmt2dat: Error while reading first year\n");
		exit (EXIT_FAILURE);
	}
	if (fread((void *)(&n_records), (size_t)4, (size_t)1, fp) != 1) {
		fprintf(stderr,"gmt2dat: Error while reading no of records\n");
		exit (EXIT_FAILURE);
	}
	if (fread((void *)agency, (size_t)10, (size_t)1, fp) != 1) {
		fprintf(stderr,"gmt2dat: Error while reading info-header\n");
		exit (EXIT_FAILURE);
	}
	
	record = (struct GMTMGG_REC *) GMT_memory (VNULL, (size_t)n_records, sizeof (struct GMTMGG_REC), "gmt2dat");
	
	gmt = gmtmgg_init (leg_year);
	
	printf ("%d %d %s\n", leg_year, n_records, agency);
	
	for (rec = 0; rec < n_records; rec++) {
		if (fread ((void *)(&record[rec]), (size_t)18, (size_t)1, fp) != 1) {
			fprintf (stderr,"gmt2dat: Error reading data record no %d\n",rec);
			exit (EXIT_FAILURE);
		}
	}
	fclose(fp);

	for (rec = 0; rec < n_records; rec++) {
		gmtmgg_date (record[rec].time, &yr, &mo, &dy, &hr, &mi, &sc, gmt);
		lat = (double) record[rec].lat * MDEG2DEG;
		lon = (double) record[rec].lon * MDEG2DEG;
		printf("%d/%.2d/%.2d/%.2d:%.2d:%.2d%10.5f%10.5f", yr, mo, dy, hr, mi, sc, lat, lon);
		(record[rec].gmt[0] == GMTMGG_NODATA) ? printf ("      NaN") : printf ("%9.2f", (float) record[rec].gmt[0] * 0.1);
		(record[rec].gmt[1] == GMTMGG_NODATA) ? printf ("    NaN") : printf ("%7d", record[rec].gmt[1]);
		(record[rec].gmt[2] == GMTMGG_NODATA) ? printf ("    NaN\n") : printf ("%7d\n", record[rec].gmt[2]);
	}
	
	GMT_free ((void *)gmt);
	GMT_free ((void *)record);

	GMT_end (argc, argv);
	
	exit (EXIT_SUCCESS);
}
