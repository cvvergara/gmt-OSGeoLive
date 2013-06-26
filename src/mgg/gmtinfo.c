/*--------------------------------------------------------------------
 *	$Id: gmtinfo.c 9923 2012-12-18 20:45:53Z pwessel $
 *
 *    Copyright (c) 1991-2013 by P. Wessel and W. H. F. Smith
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * gmtinfo prints out a string of information pertaining to each specified leg.
 * The parameters are: legname, year, no of records in file, no of gravity,
 * magnetics, and bathymetry points, area covered in lat/lon, and the d/m/yr
 * of the first and last records in the file(s).
 *
 * Author:	Paul Wessel
 * Date:	3-FEB-1988
 * Version:	3.2, 12-APR-1999 (DOS compatible)
 *
 */
#include "gmt.h"
#include "gmt_mgg.h"

int main (int argc, char **argv)
{
	int leg_year, n_records;
	int rec, n_grv, n_mag, n_top;
	int begin_day, begin_mo, begin_yr, end_day, end_mo, end_yr;
	int hour, minute, second, nfiles=0, i, verbose = FALSE;
	double xmin1, xmax1, xmin2, xmax2, ymin, ymax, xmin, xmax, lat, lon, old_x = 0.0;
	char agency[10], greenwich, gmtfile[100];
	FILE *fp = NULL;
	struct GMTMGG_REC record;
	struct GMTMGG_TIME *gmt = NULL;
	
	argc = (int)GMT_begin (argc, argv);

	gmtmggpath_init();
	
	for (i = 1; i < argc; i++) if (!strcmp(argv[i], "-v")) verbose = TRUE;
	
	for (i = 1; i < argc; i++) {	/* Loop over gmt files */
		if (argv[i][0] == '-') continue;
  		if (gmtmggpath_func (gmtfile, argv[i])) {
   			fprintf(stderr,"gmtinfo : Cannot find leg %s\n", argv[i]);
     			continue;
  		}
		if ((fp = fopen (gmtfile, "rb")) == NULL) {
			fprintf(stderr,"gmtinfo: Could not open file %s\n", gmtfile);
			continue;
		}
		nfiles++;
	
		if (verbose) fprintf (stderr, "Now reading leg %s\n", argv[i]);
		
		/* Read first record of file containing start-year, n_records and info */
	
		if (fread((void *)(&leg_year), (size_t)4, (size_t)1, fp) != 1) {
			fprintf(stderr,"gmtinfo: Error while reading first year\n");
			exit (EXIT_FAILURE);
		}
		if (fread((void *)(&n_records), (size_t)4, (size_t)1, fp) != 1) {
			fprintf(stderr,"gmtinfo: Error while reading no of records\n");
			exit (EXIT_FAILURE);
		}
		if (fread((void *)agency, (size_t)10, (size_t)1, fp) != 1) {
			fprintf(stderr,"gmtinfo: Error while reading info-header\n");
			exit (EXIT_FAILURE);
		}
	
		gmt = gmtmgg_init (leg_year);	/* Initialize gmt_structure */
	
		/* Start reading data from file */
	
		xmin1 = xmin2 = 360.0;
		xmax1 = xmax2 = -360.0;
		ymin = 180.0;
		ymax = -180.0;
		n_grv = n_mag = n_top = 0;
		greenwich = FALSE;
		for (rec = 0; rec < n_records; rec++) {
			if (fread((void *)(&record), (size_t)18, (size_t)1, fp) != 1) {
				fprintf(stderr,"gmtinfo: Error reading data record no %d\n",rec);
				exit (EXIT_FAILURE);
			}
			if (rec == 0)
				gmtmgg_date(record.time,&begin_yr,&begin_mo,&begin_day,&hour,&minute,&second,gmt);

			if (record.gmt[0] != GMTMGG_NODATA) n_grv++;
			if (record.gmt[1] != GMTMGG_NODATA) n_mag++;
			if (record.gmt[2] != GMTMGG_NODATA) n_top++;
			lon = record.lon * 0.000001;
			lat = record.lat * 0.000001;
			if (lon < 180.0) {
				xmin1 = MIN (lon, xmin1);
				xmax1 = MAX (lon, xmax1);
			}
			else {
				xmin2 = MIN (lon, xmin2);
				xmax2 = MAX (lon, xmax2);
			}
			ymin = MIN (lat, ymin);
			ymax = MAX (lat, ymax);
			if (rec > 0 && (fabs(old_x-lon) > 180.0))
				greenwich = TRUE;
			old_x = lon;
		}
		fclose (fp);
		
		gmtmgg_date (record.time,&end_yr,&end_mo,&end_day,&hour,&minute,&second,gmt);
		GMT_free ((void *)gmt);

		if (greenwich) {
			xmin = MAX(xmin1,xmin2);
			xmax = MIN(xmax1,xmax2);
		}
		else {
			xmin = MIN(xmin1,xmin2);
			xmax = MAX(xmax1,xmax2);
		}
		if (xmin > xmax) xmin -= 360.0;
		printf("Info: %s Year: %d N_recs: %d n_g: %d n_m: %d n_t: %d W: %.5f E: %.5f S: %.5f N: %.5f Begins: %d %d %d Ends: %d %d %d\n",
			agency, leg_year, n_records, n_grv, n_mag, n_top, xmin, xmax, ymin, ymax,
			begin_day, begin_mo, begin_yr, end_day, end_mo, end_yr);
	}
	
	if (nfiles == 0) {
		fprintf(stderr, "usage: gmtinfo leg(s)\n");
		exit (EXIT_FAILURE);
	}

	gmtmgg_end ();

	GMT_end (argc, argv);
	
	exit (EXIT_SUCCESS);
}
