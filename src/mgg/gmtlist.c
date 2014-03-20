/*--------------------------------------------------------------------
 *	$Id: gmtlist.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *    Copyright (c) 1991-2014 by P. Wessel and W. H. F. Smith
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * gmtlist produces ASCII listings of <legid>.gmt files. The *.gmt files
 * contains time(s), latitude(y), longitude(x), gravity(g), magnetics(m),
 * and bathymetry(t), and the user may extract any combination of these 6
 * parameters + distance (in km), heading, velocity, and weight by using
 * the option -Fsxygmtdhvw.  The sequence in which the flag characters
 * appear determines the sequence in which the parameters will be printed
 * out. If no options is specified, the default is -Fsxygmtdhvw.  E.g. to
 * create an input file for surface, use -Fxyg (for gravity).
 * If upper case letters are used for gmt (GMT), then only records that have
 * that particular data are written out.  E.g -Fxyg gives
 * lon/lat/grav, whereas -FxyG gives lon/lat/grav where there is gravity data.
 * To select a section of the track, specify the start/endpoints by:
 *	1) Start-time (mm/dd/yyyy/hh:mm) OR start-distance (km)
 *	2) Stop-time (mm/dd/yyyy/hh:mm) OR stop-distance (km)
 * To select data inside an area, use the -R option.
 * To start output with a header string, use -H.
 * Several formats for time is available. The default when using -Fs is
 * seconds from Jan 1 the year the cruise started. -Fsc (calender) gives
 * yyyymmddhhmms output, and -Fsj (julian) gives yyyyjjhhmmss output.
 *
 * Author:	Paul Wessel
 * Date:	19-APR-1988
 * Version:	2.1 1-JUL-1992
 *		3.2 10-MAR-1999
 *		3.3.1 25-JUN-1999
 *
 *
 */

#include "gmt.h"
#include "gmt_mgg.h"
#include "x_system.h"

#define KMPRDEG 111.1949e-6
#define MAXLEGS 5000
#define S_PR_DAY 86400

int main (int argc, char **argv)
{
	int leg_year, n_records, rec, i, j, argno, id, no[10], nval = 0;
	int start_time = 0, stop_time = 2000000000, dlon, last_lon = 0;
	int mon1, day1 = 0, year1, hour1, min1, time, last_lat = 0, mon2, day2 = 0, year2, hour2, min2;
	int yy, mm, dd, hh, mi, ss, jd, dt, last_time = 0, n_cruises = 0;
	size_t not_used = 0;

	GMT_LONG error = FALSE, wantgmt, want_all = FALSE, geodetic = TRUE;
	GMT_LONG correct = FALSE, tsec = FALSE, calender = FALSE, do_heading = FALSE;
	GMT_LONG no_g, no_m, no_t, greenwich = FALSE, do_speed = FALSE, binary = FALSE;

	char gmtfile[BUFSIZ], corrfile[BUFSIZ], agency[10], ndata, g, m, t;

	double lat, lon, grv, mag, top, dist, val, start_dist = 0., stop_dist = 1.0E100;
	double ds, dx, dy, west = 0.0, east = 360.0, south = -90.0, north = 90.0;
	double heading, speed, weight = 1.0;

	struct CORR **bin = NULL;

	int binsize = sizeof (struct CORR);
	int nlegs = 0;

	int get_id (char *name, int nlegs, struct CORR **bin);

	struct GMTMGG_TIME *gmt = NULL;

	FILE *fp = NULL, *fpc = NULL;

	struct GMTMGG_REC record;

	g = m = t = FALSE;
	no_g = no_m = no_t = FALSE;

	/* Check and interpret the command line arguments */

	argc = (int)GMT_begin (argc, argv);

	for (i =1; !error && i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {

				case 'H':
				case 'R':
				case 'V':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				case 'b':	/* Binary output */
					binary = TRUE;
					break;

				case 'C':
					if (argv[i][2] == 0)	/* Use default corrfile */
						GMT_getsharepath ("mgg", "cxx_corrections", ".b", corrfile);
					else
						strcpy (corrfile, &argv[i][2]);
					correct = TRUE;
					break;

				case 'D':		/* Assign start/stop times for sub-section */
					if (argv[i][2] == 'a') {	/* Start date */
						sscanf(&argv[i][3], "%d/%d/%d/%d:%d",
							&mon1, &day1, &year1, &hour1, &min1);
					}
					else if (argv[i][2] == 'b')	 {	/* Stop date */
						sscanf(&argv[i][3], "%d/%d/%d/%d:%d",
							&mon2, &day2, &year2, &hour2, &min2);
					}
					else
						error = TRUE;
					break;

				case 'F':	/* Selected output fields */
					for (j = 2; argv[i][j]; j++) {
						switch (argv[i][j]) {
							case 'G':	/* Records with gravity != GMTMGG_NODATA requested */
								no_g = TRUE;
	       						case 'g':		/* Gravity is requested */
								no[nval++] = 3;
								g = 1;
								break;
							case 'M':	/* Records with magnetics != GMTMGG_NODATA requested */
								no_m = TRUE;
							case 'm':		/* Magnetics is requested */
								no[nval++] = 4;
								m = 1;
								break;
							case 'T':	/* Records with topo != GMTMGG_NODATA requested */
								no_t = TRUE;
							case 't':		/* Topography is requested */
								no[nval++] = 5;
								t = 1;
								break;
							case 'x':		/* Longitude is requested */
								no[nval++] = 1;
								break;
							case 'y':		/* Latitude is requested */
								no[nval++] = 2;
								break;
							case 's':		/* Time (in sec) is requested */
								no[nval++] = 0;
								if (argv[i][j+1] == 'c') {
									calender = TRUE;
									j++;
								}
								else if (argv[i][j+1] == 'j')
									j++;
								else
									tsec = TRUE;
								break;
							case 'd':		/* Distance (in km) is requested */
								no[nval++] = 6;
								break;
							case 'h':		/* Heading is requested */
								no[nval++] = 7;
								do_heading = TRUE;
								break;
							case 'v':		/* velocity (in m/s) is requested */
								no[nval++] = 8;
								do_speed = TRUE;
								break;
							case 'w':		/* weights (Set with -W) is requested */
								no[nval++] = 9;
								break;
							default:
								error = TRUE;
								break;
						}
					}
					break;

				case 'S':		/* Assign start/stop position for sub-section */
					if (argv[i][2] == 'a')	/* Start position */
						start_dist = atof(&argv[i][3]);
					else if (argv[i][2] == 'b')	/* Stop position */
						stop_dist = atof(&argv[i][3]);
					else
						error = TRUE;
					break;

				case 'G':
					geodetic = FALSE;
					break;
				case 'W':		/* Assign a weight to these data */
					weight = atof (&argv[i][2]);
					break;

				default:		/* Options not recognized */
					error = TRUE;
					break;
			}
		}
		else
			n_cruises++;
	}

	/* Check that the options selected are mutually consistent */

	if (nval > 10) error = TRUE;
	if (start_dist > stop_dist || start_time > stop_time) error = TRUE;
	if ((day1 > 0 && start_dist > 0.) || (day2 > 0 && stop_dist < 1.0e100)) error = TRUE;
	if (east < west || south > north) error = TRUE;
	if (n_cruises == 0) error = TRUE;
	if (weight <= 0.0) error = TRUE;

	if (error || argc == 1) {	/* Display usage */
		fprintf(stderr,"usage: gmtlist <cruise(s)> [-C<corrfile>] [-Da<startdate>] [-Db<stopdate>] [-F<dataflags>]\n");
		fprintf(stderr,"	[-G] [-H] [%s] [-Sa<startdist>] [-Sb<stopdist>] [-V] [-W<Weight>] [-b]\n\n", GMT_Rgeo_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf(stderr,"	<cruises> is one or more legnames, e.g., c2104 v3206 etc.\n");
		fprintf(stderr,"	OPTIONS:\n\n");
		fprintf(stderr,"	-C<file> applies crossover corrections to the data. If no file name\n");
		fprintf(stderr,"	   is given, the default correction file is assumed\n");
		fprintf(stderr,"	-Da<date> lists from date (given as mm/dd/yr/hh:mm)\n");
		fprintf(stderr,"	-Db<date> lists up to date (given as mm/dd/yr/hh:mm)\n");
		fprintf(stderr,"	-F Dataflags is a string made up of 1 or more of these characters:\n");
		fprintf(stderr,"	  s means list time in seconds, sc gives dates, sj gives Julian day\n");
		fprintf(stderr,"	  x means list longitude (degrees)\n");
		fprintf(stderr,"	  y means list latitude (degrees)\n");
		fprintf(stderr,"	  g means list gravity (mGal)\n");
		fprintf(stderr,"	  m means list magnetics (nTesla)\n");
		fprintf(stderr,"	  t means list topography (m)\n");
		fprintf(stderr,"	  d means list distance (km)\n");
		fprintf(stderr,"	  h means list heading (Degrees east from north)\n");
		fprintf(stderr,"	  v means list velocity (m/s)\n");
		fprintf(stderr,"	  w means list weight (see -W)\n");
		fprintf(stderr,"	  If G, M, or T is used instead of g, m, or t, then only the records\n");
		fprintf(stderr,"	  that have that combination of data will be listed\n");
		fprintf(stderr,"	  The data is written out in the order specified in <dataflags>\n");
		fprintf(stderr,"	  [Default is -Fsxygmtdhvw and all records]\n");
		fprintf(stderr,"	-G force geographical longitudes (-180/+180) [Default is 0-360]\n");
		fprintf(stderr,"	-H write header record\n");
		fprintf(stderr,"	-R only return data inside the specified region\n");
		fprintf(stderr,"	-Sa<dist> lists from dist (in km)\n");
		fprintf(stderr,"	-Sb<dist> lists up to dist (in km)\n");
		fprintf(stderr,"	-V verbose, report progress\n");
		fprintf(stderr,"	-W sets weight for these data\n");
		fprintf(stderr,"	-b means binary (double) output [ascii]\n");
		exit (EXIT_FAILURE);
	}

	gmtmggpath_init();

	if ((west < 0.0 && east > 0.0) || (west < 360.0 && east > 360.0)) greenwich = TRUE;
	if (!geodetic) greenwich = TRUE;

	if (correct) {	/* Read correction table */
		if ((fpc = fopen(corrfile, "rb")) == NULL) {
			fprintf(stderr, "Could not read correction file %s\n", corrfile);
			exit (EXIT_FAILURE);
		}
		bin = (struct CORR **) GMT_memory (VNULL, (size_t)MAXLEGS, sizeof (struct CORR *), "gmtlist");
		i = 0;
		bin[i] = (struct CORR *) GMT_memory (VNULL, (size_t)1, sizeof (struct CORR), "gmtlist");
		while (fread((void *)bin[i], (size_t)binsize, (size_t)1, fpc) == 1) {
			i++;
			bin[i] = (struct CORR *) GMT_memory (VNULL, (size_t)1, sizeof (struct CORR), "gmtlist");
		}
		fclose(fpc);
		nlegs = i;
	}


	/* Sort the  order in which the parameters appear */

	if (nval == 0) {		/* Nothing selected, default used */
		g = m = t = TRUE;	/* No data was specified so all data [default] is output */
		for (i = 0; i < 10; i++) no[i] = i;
		nval = 10;
		want_all = tsec = TRUE;
		do_heading = do_speed = TRUE;
	}

	if (!binary && GMT_io.io_header[GMT_OUT]) {
		/* Write out header record */
		for (i = 0; i < nval; i++) {
			switch(no[i]) {
				case 0:	/* Print out time header */
					if (tsec)
						printf("time(s)");
					else if (calender)
						printf("yyyy\tmm\tdd\thh\tmm\tss");
					else
						printf("yyyy\tjd\thh\tmm\tss");
					break;
				case 1:	/* Print out longitude */
					printf("lon");
					break;
				case 2:	/* Print out latitude */
					printf("lat");
					break;
				case 3:	/* Print out gravity */
		 			printf("faa(mGal)");
		 			break;
				case 4:	/* Print out magnetics */
		 			printf("mag(nTesla)");
					break;
				case 5:	/* Print out bathymetry */
					printf("topo(m)");
					break;
				case 6:	/* Print out distance */
					printf("dist(km)");
					break;
				case 7:	/* Print out heading */
					printf("heading");
					break;
				case 8:	/* Print out velocity */
					printf("speed(m/s)");
					break;
				case 9:	/* Print out weights */
					printf("weight");
					break;
			}
			((i+1) < nval) ? printf("\t") : printf("\n");
		}
	}

	for (argno = 1; argno < argc; argno++) {	/* Loop over all the files */

		if (argv[argno][0] == '-') continue;
  		if (gmtmggpath_func (gmtfile, argv[argno])) {
   			fprintf (stderr, "gmtlist : Cannot find leg %s\n", argv[argno]);
     			continue;
  		}
		if ((fp = fopen (gmtfile, "rb")) == NULL) {
			fprintf (stderr,"gmtlist: Could not open %s\n", gmtfile);
			continue;
		}

		if (gmtdefs.verbose) fprintf (stderr, "gmtlist: Now processing cruise %s\n",  argv[argno]);

		/* Read first record of file containing start-year, n_records and info */

		if (fread ((void *)(&leg_year), (size_t)4, (size_t)1, fp) != 1) {
			fprintf (stderr,"gmtlist: Error while reading first year\n");
			exit (EXIT_FAILURE);
		}
		if (fread ((void *)(&n_records), (size_t)4, (size_t)1, fp) != 1) {
			fprintf (stderr,"gmtlist: Error while reading no of records\n");
			exit (EXIT_FAILURE);
		}
		if (fread ((void *)agency, (size_t)10, (size_t)1, fp) != 1) {
			fprintf (stderr,"gmtlist: Error while reading info-header\n");
			exit (EXIT_FAILURE);
		}

		gmt = gmtmgg_init (leg_year);	/* Initialize gmt_structure */
		if (correct)
			id = get_id(argv[argno], nlegs, bin);
		else
			id = -1;

		/* Decode date to time in sec if needed */

		if (day1 > 0) gmtmgg_time (&start_time, year1, mon1, day1, hour1, min1, 0, gmt);
		if (day2 > 0) gmtmgg_time (&stop_time, year2, mon2, day2, hour2, min2, 0, gmt);

		dist = 0.0;
		time = 0;

		wantgmt = (g || m || t) ? TRUE : FALSE;
		if (want_all) wantgmt = FALSE;

		/* Start reading data from file */

		for (rec = 0; rec < n_records && dist < stop_dist && time < stop_time; rec++) {
			if (fread ((void *)(&record), (size_t)18, (size_t)1, fp) != 1) {
				fprintf (stderr,"gmtlist: Error reading data record no %d\n",rec);
				exit (EXIT_FAILURE);
			}

			/* Compute accumulated distance along track (Flat Earth) */

			if (rec == 0) {
				last_lon = record.lon;
				last_lat = record.lat;
				last_time = record.time;
				ds = 0.0;
				dt = 0;
				heading = speed = GMTMGG_NODATA;
			}
			else {
				dlon = record.lon - last_lon;
				if (abs (dlon) > 180000000) dlon = irint(copysign ((double) (360000000 - abs (dlon)), (double)dlon));
				dx = (double) dlon * cosd (0.5e-06*(double)(record.lat+last_lat));
				dy = (double) (record.lat - last_lat);
				ds = KMPRDEG * hypot (dx, dy);
				dt = record.time - last_time;
				if (do_heading) {
					heading = (dx == 0.0 && dy == 0.0) ? GMTMGG_NODATA : 90.0 - R2D * atan2 (dy, dx);
					if (heading < 0.0) heading += 360.0;
				}
				if (do_speed) speed = (dt == 0) ? GMTMGG_NODATA : 1000.0 * ds / dt;
				last_lon = record.lon;
				last_lat = record.lat;
				last_time = record.time;
			}
			dist += ds;

			/* Check if record has the required fields */

			if (no_g && record.gmt[0] == GMTMGG_NODATA) continue;
			if (no_m && record.gmt[1] == GMTMGG_NODATA) continue;
			if (no_t && record.gmt[2] == GMTMGG_NODATA) continue;

			/* Check if time or dist falls outside specified range */

			if (dist < start_dist) continue;
			if (dist > stop_dist) continue;
			if (record.time < start_time) continue;
			if (record.time > stop_time) continue;

			time = record.time;
			lat = (double) record.lat*0.000001;
			lon = (double) record.lon*0.000001;

			/* Check is lat/lon is outside specified area */

			if (lat < south || lat > north) continue;
			while (lon > east) lon -= 360.0;
			while (lon < west) lon += 360.0;
			if (lon > east) continue;
			while (lon > 360.0) lon -= 360.0;

			grv = (record.gmt[0] != GMTMGG_NODATA) ? (double) record.gmt[0]*0.1 : GMTMGG_NODATA;
			mag = record.gmt[1];
			top = record.gmt[2];

			ndata = (g && record.gmt[0] != GMTMGG_NODATA) ? 1 : 0;
			ndata += (m && record.gmt[1] != GMTMGG_NODATA) ? 1 : 0;
			ndata += (t && record.gmt[2] != GMTMGG_NODATA) ? 1 : 0;
			if (ndata == 0 && wantgmt) continue;

			/* This record will now be printed out */

			for (i = 0; i < nval; i++) {
				switch(no[i]) {
					case 0:	/* Print out time */
						if (binary) {
							val = record.time;
							not_used = fwrite ((void *)&val, sizeof (double), (size_t)1, stdout);
						}
						else if (tsec)
							printf ("%d", record.time);
						else if (calender) {
							gmtmgg_date (record.time,&yy,&mm,&dd,&hh,&mi,&ss,gmt);
							printf ("%d\t%2d\t%2d\t%2d\t%2d\t%2d",
								yy, mm, dd, hh, mi, ss);
						}
						else {	/* julian day etc. */
							jd = gmtmgg_date (record.time,&yy,&mm,&dd,&hh,&mi,&ss,gmt);
							printf ("%d\t%3d\t%2d\t%2d\t%2d",
								yy, jd, hh, mi, ss);
						}
						break;
					case 1:	/* Print out longitude */
						if (lon > 180.0 && greenwich) lon -= 360.0;
						if (binary)
							not_used = fwrite ((void *)&lon, sizeof (double), (size_t)1, stdout);
						else
							(GMT_io.io_header[GMT_OUT]) ? printf("%.5f", lon) : printf("%9.5f", lon);
						break;
					case 2:	/* Print out latitude */
						if (binary)
							not_used = fwrite ((void *)&lat, sizeof (double), (size_t)1, stdout);
						else
							(GMT_io.io_header[GMT_OUT]) ? printf("%.5f", lat) : printf("%8.5f", lat);
						break;
					case 3:	/* Print out gravity */
						if (record.gmt[0] == GMTMGG_NODATA)
							(binary) ? not_used = fwrite ((void *)&GMT_f_NaN, sizeof (double), (size_t)1, stdout) : printf ("NaN");
						else {
							if (id >= 0) grv -= bin[id]->dc_shift_gmt[0] + bin[id]->drift_rate_gmt[0] * record.time;
							if (binary)
								not_used = fwrite ((void *)&grv, sizeof (double), (size_t)1, stdout);
							else
		 						(GMT_io.io_header[GMT_OUT]) ? printf("%.2f", grv) : printf("%9.2f", grv);
						}
		 				break;
					case 4:	/* Print out magnetics */
						if (record.gmt[1] == GMTMGG_NODATA)
							(binary) ? not_used = fwrite ((void *)&GMT_f_NaN, sizeof (double), (size_t)1, stdout) : printf ("NaN");
						else {
							if (id >= 0) mag -= bin[id]->dc_shift_gmt[1] + bin[id]->drift_rate_gmt[1] * record.time;
							if (binary)
								not_used = fwrite ((void *)&mag, sizeof (double), (size_t)1, stdout);
							else
		 						(GMT_io.io_header[GMT_OUT]) ? printf("%.1f", mag) : printf("%8.1f", mag);
						}
						break;
					case 5:	/* Print out bathymetry */
						if (record.gmt[2] == GMTMGG_NODATA)
							(binary) ? not_used = fwrite ((void *)&GMT_f_NaN, sizeof (double), (size_t)1, stdout) : printf ("NaN");
						else {
							if (id >= 0) top -= bin[id]->dc_shift_gmt[2] + bin[id]->drift_rate_gmt[2] * record.time;
							if (binary)
								not_used = fwrite ((void *)&top, sizeof (double), (size_t)1, stdout);
							else
								(GMT_io.io_header[GMT_OUT]) ? printf("%.1f", top) : printf("%8.1f", top);
						}
						break;
					case 6:	/* Print out distance */
						if (binary)
							not_used = fwrite ((void *)&dist, sizeof (double), (size_t)1, stdout);
						else
							(GMT_io.io_header[GMT_OUT]) ? printf("%.3f", dist) : printf("%9.3f", dist);
						break;
					case 7:	/* Print out heading */
						if (heading == GMTMGG_NODATA)
							(binary) ? not_used = fwrite ((void *)&GMT_f_NaN, sizeof (double), (size_t)1, stdout) : printf ("NaN");
						else {
							if (binary)
								not_used = fwrite ((void *)&heading, sizeof (double), (size_t)1, stdout);
							else
								(GMT_io.io_header[GMT_OUT]) ? printf("%.1f", heading) : printf("%6.1f", heading);
						}
						break;
					case 8:	/* Print out velocity */
						if (speed == GMTMGG_NODATA)
							(binary) ? not_used = fwrite ((void *)&GMT_f_NaN, sizeof (double), (size_t)1, stdout) : printf ("NaN");
						else {
							if (binary)
								not_used = fwrite ((void *)&speed, sizeof (double), (size_t)1, stdout);
							else
								(GMT_io.io_header[GMT_OUT]) ? printf("%.2f", speed) : printf("%6.2f", speed);
						}
						break;
					case 9:	/* Print out weight */
						if (binary)
							not_used = fwrite ((void *)&weight, sizeof (double), (size_t)1, stdout);

						else
							(GMT_io.io_header[GMT_OUT]) ? printf("%g", weight) : printf("%6.2f", weight);
						break;
				}
				if (!binary) ((i+1) < nval) ? printf("\t") : printf("\n");
			}
		}
		fclose (fp);
		GMT_free ((void *)gmt);
	}

	gmtmgg_end ();
	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

int get_id (char *name, int nlegs, struct CORR **bin)
{
	int left, right, mid, cmp;

	left = 0;
	right = nlegs-1;
	while (left <= right) {
		mid = (left + right)/2;
		cmp = strcmp(name, bin[mid]->name);
		if (cmp < 0)
			right = mid-1;
		else if (cmp > 0)
			left = mid+1;
		else
			return (mid);
	}
	return (-1);
}
