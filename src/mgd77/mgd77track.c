/*--------------------------------------------------------------------
 *	$Id: mgd77track.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *    Copyright (c) 2004-2014 by P. Wessel
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * mgd77track will read *.mgd77-files and and write PostScript code
 * that will create a navigation plot on a basemap using the
 * projection specified by the user. WESN must be specified on the 
 * command line along with other options for scaling, annotation etc.
*
 * To select a sub-section of the track, specify the start/endpoints by:
 *	1) Start-time (yyyy-mm-ddT[hh:mm:ss]) OR start-distance (km)
 *	2) Stop-time  (yyyy-mm-ddT[hh:mm:ss]) OR stop-distance (km)
 *
 * Author:	Paul Wessel
 * Date:	19-AUG-2004
 * Version:	1.0 Based alot on the old gmttrack.c
 *
 *
 */
 
#include "mgd77.h"
#include "pslib.h"

#define MGD77TRACK_ANSIZE 0.125
#define MGD77TRACK_MARK_NEWDAY	0
#define MGD77TRACK_MARK_SAMEDAY	1
#define MGD77TRACK_MARK_DIST	2

#define ANNOT	0
#define LABEL	1

struct MGD77TRACK_ANNOT {
	double annot_int_dist;
	double tick_int_dist;
	double annot_int_time;
	double tick_int_time;
};

struct MGD77TRACK_LEG_ANNOT {	/* Structure used to annotate legs after clipping is terminated */
	double x, y;
	double lon, lat;
	double angle;
	char text[16];
};

struct MGD77TRACK_MARKER {
	double marker_size, font_size;
	struct GMT_FILL f, s;	/* Font and symbol colors */
	GMT_LONG font_no;
};

int main (int argc, char **argv)
{
	int i, j, n_id = 0, mrk = 0, annot_name = 0, dist_flag = 2, use, n_paths;
	
	GMT_LONG rec, first_rec, last_rec, n_alloc_c = GMT_SMALL_CHUNK;
	GMT_LONG this_julian = 0, last_julian, argno, n_cruises = 0;
	GMT_LONG error = FALSE, first, both = FALSE, no_clip = FALSE;
	GMT_LONG  annot_tick[2] = {0, 0}, draw_tick[2] = {0, 0}, dist_gap = FALSE, time_gap = FALSE;
	
	char *start_date = NULL, *stop_date = NULL, ms[GMT_TEXT_LEN], mc[GMT_TEXT_LEN];
       	char label[GMT_LONG_TEXT], date[GMT_TEXT_LEN], clock[GMT_TEXT_LEN], comment[BUFSIZ];
	char mfs[GMT_TEXT_LEN], mf[GMT_TEXT_LEN], mfc[GMT_TEXT_LEN], name[GMT_TEXT_LEN], **list = NULL, *p = NULL;

	double west, east, south, north, start_time, stop_time, start_dist, gap_d = 0.0, gap_t = 0.0;
	double stop_dist, x, y, annot_dist[2] = {0, 0}, tick_dist[2] = {0, 0}, annot_time[2] = {0, 0}, tick_time[2] = {0, 0};
	double *track_dist = NULL, angle, plot_x, plot_y, *lon = NULL, *lat = NULL, *track_time = NULL;
	double annotsize = MGD77TRACK_ANSIZE, factor = 0.001, scale = 1.0, c_angle;
	
	struct MGD77_CONTROL M;
	struct MGD77_DATASET D;
	struct MGD77TRACK_LEG_ANNOT *cruise_id = NULL;
	struct GMT_PEN pen;
	struct GMT_gcal calendar;
	struct MGD77TRACK_MARKER marker[3];	/* Time (2) and distance marker information */
	struct MGD77TRACK_ANNOT info[2];

	int get_annotinfo (char *args, struct MGD77TRACK_ANNOT *info);
	void annot_legname (double x, double y, double lon, double lat, double angle, char *text, double size);
	double heading (GMT_LONG rec, double *lon, double *lat, GMT_LONG n_records);
	GMT_LONG bad_coordinates (double lon, double lat);

	argc = (int)GMT_begin (argc, argv);		/* Initialize GMT Machinery */
	
	GMT_init_pen (&pen, GMT_PENWIDTH);
	GMT_init_fill (&marker[MGD77TRACK_MARK_NEWDAY].s, 0, 0, 0);	/* Black color for new Julian day marker */
	GMT_init_fill (&marker[MGD77TRACK_MARK_SAMEDAY].s, 255, 255, 255);	/* White color for other time markers in same day */
	GMT_init_fill (&marker[MGD77TRACK_MARK_DIST].s, 0, 0, 0);		/* Black color for distance markers */
	GMT_init_fill (&marker[MGD77TRACK_MARK_NEWDAY].f, 0, 0, 0);	/* Black color for all text */
	GMT_init_fill (&marker[MGD77TRACK_MARK_SAMEDAY].f, 0, 0, 0);	/* Black color for all text */
	GMT_init_fill (&marker[MGD77TRACK_MARK_DIST].f, 0, 0, 0);		/* Black color for all text */
	if (gmtdefs.measure_unit == GMT_CM) {
		marker[MGD77TRACK_MARK_NEWDAY].marker_size = marker[MGD77TRACK_MARK_SAMEDAY].marker_size = 0.1 / 2.54;	/* 1 mm */
		marker[MGD77TRACK_MARK_DIST].marker_size = 0.15 / 2.54;	/* 1.5 mm */
	}
	else {	/* Assume we think in inches */
		marker[MGD77TRACK_MARK_NEWDAY].marker_size = marker[MGD77TRACK_MARK_SAMEDAY].marker_size = 0.04;
		marker[MGD77TRACK_MARK_DIST].marker_size = 0.06;
	}
	marker[MGD77TRACK_MARK_NEWDAY].font_size = marker[MGD77TRACK_MARK_SAMEDAY].font_size = marker[MGD77TRACK_MARK_DIST].font_size = gmtdefs.annot_font_size[0] * GMT_u2u[GMT_PT][GMT_INCH];
	marker[MGD77TRACK_MARK_NEWDAY].font_no = GMT_font_lookup ("Times-BoldItalic", GMT_font, GMT_N_FONTS);
	marker[MGD77TRACK_MARK_SAMEDAY].font_no = GMT_font_lookup ("Times-Italic", GMT_font, GMT_N_FONTS);
	marker[MGD77TRACK_MARK_DIST].font_no = GMT_font_lookup ("Times-Roman", GMT_font, GMT_N_FONTS);

	/* Initialize MGD77 output order and other parameters*/
	
	MGD77_Init (&M);			/* Initialize MGD77 Machinery */
	west = 0.0;	east = 360.0;	south = -90.0;	north = 90.0;
	start_time = start_dist = 0.0;
	stop_time = stop_dist = DBL_MAX;
	memset ((void *)info, 0, 2*sizeof (struct MGD77TRACK_ANNOT));

	for (i = 1; !error && i < argc; i++) {	/* Process input options */
		if (argv[i][0] != '-') continue;

		switch(argv[i][1]) {

			case 'B':
			case 'P':
			case 'O':
			case 'K':
			case 'R':
			case 'J':
			case 'U':
			case 'V':
			case 'X':
			case 'x':
			case 'Y':
			case 'y':
			case 'c':
			case '\0':
				error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
				break;
					
			case 'A':
				annot_name = 1;
				j = 2;
				if (argv[i][2] == 'c') j = 3, annot_name = 2;
				if (argv[i][j] && argv[i][j] != ',') annotsize = atof (&argv[i][j]) * GMT_u2u[GMT_PT][GMT_INCH];
				if ((p = strchr (&argv[i][2], ','))) {	/* Want label at regular spacing */
					p++;	/* Skip the comma */
					error += get_annotinfo (p, &info[LABEL]);
					annot_name = -annot_name;	/* Flag to tell machinery not to annot at entry */
				}
				break;

			case 'C':	/* Distance calculation flag */
				if (argv[i][2] == 'f') dist_flag = 1;
				if (argv[i][2] == 'g') dist_flag = 2;
				if (argv[i][2] == 'e') dist_flag = 3;
				if (dist_flag < 1 || dist_flag > 3) {
					fprintf(stderr, "%s: ERROR -C: Flag must be f, g, or e\n", GMT_program);
					error++;
				}
				break;
			case 'D':		/* Assign start/stop times for sub-section */
				if (argv[i][2] == 'a') {		/* Start date */
					start_date = &argv[i][3];
				}
				else if (argv[i][2] == 'b') {		/* Stop date */
					stop_date = &argv[i][3];
				}
				else
					error = TRUE;
				break;

			case 'F':	/* Do NOT apply bitflags */
				M.use_flags[MGD77_M77_SET] = M.use_flags[MGD77_CDF_SET] = FALSE;
				break;
				
			case 'G':
				switch (argv[i][2]) {
					case 'd':	/* Distance gap in km */
						gap_d = atof (&argv[i][3]) * 1000.0;	/* Gap converted to m from km */
						dist_gap = TRUE;
						break;
					case 't':	/* Distance gap in minutes */
						gap_t = atof (&argv[i][3]) * 60.0;	/* Gap converted to seconds from minutes */
						time_gap = TRUE;
						break;
					default:
						fprintf(stderr, "%s: ERROR -G: Requires t|d and a positive value in km (d) or minutes (t)\n", GMT_program);
						error++;
						break;
				}
				break;
				
			case 'I':
				MGD77_Process_Ignore (argv[i][1], &argv[i][2]);
				break;
			case 'L':
				error += get_annotinfo (&argv[i][2], &info[ANNOT]);
				break;
	
			case 'N':
				no_clip = TRUE;
				break;
	
			case 'S':		/* Assign start/stop position for sub-section (in meters) */
				if (argv[i][2] == 'a') {		/* Start position */
					MGD77_Set_Unit (&argv[i][3], &scale, 1);
					start_dist = atof (&argv[i][3]) * scale;
				}
				else if (argv[i][2] == 'b') {		/* Stop position */
					MGD77_Set_Unit (&argv[i][3], &scale, 1);
					stop_dist = atof (&argv[i][3]) * scale;
				}
				else
					error = TRUE;
				break;
	
			case 'T':	/* Marker attributes */
				switch (argv[i][2]) {
					case 'T':	/* New day marker */
						mrk = MGD77TRACK_MARK_NEWDAY;
						break;
					case 't':	/* Same day marker */
						mrk = MGD77TRACK_MARK_SAMEDAY;
						break;
					case 'd':	/* Distance marker */
						mrk = MGD77TRACK_MARK_DIST;
						break;
					default:
						error = TRUE;
						break;
				}
				if (error) {
					fprintf(stderr, "%s: ERROR: Unrecognized modifier %c given to -T\n", GMT_program, argv[i][2]);
					exit (EXIT_FAILURE);
				}
				strcpy (comment, &argv[i][3]);
				for (j = 0; j < (int)strlen (comment); j++) if (comment[j] == ',') comment[j] = ' ';	/* Replace commas with spaces */
				j = sscanf (comment, "%s %s %s %s %s", ms, mc, mfs, mf, mfc);
				if (j != 5) {
					fprintf(stderr, "%s: ERROR: -TT|t|d takes 5 arguments\n", GMT_program);
					exit (EXIT_FAILURE);
				}
				marker[mrk].marker_size = GMT_convert_units (ms, GMT_INCH);
				GMT_getfill (mc, &marker[mrk].s);
				marker[mrk].font_size = atof (mfs) * GMT_u2u[GMT_PT][GMT_INCH];
				marker[mrk].font_no = GMT_font_lookup (mf, GMT_font, GMT_N_FONTS);
				GMT_getfill (mfc, &marker[mrk].f);
				break;
			
			case 'W':
				GMT_getpen (&argv[i][2], &pen);
				break;
	
			default:		/* Options not recognized */
				error = TRUE;
				break;
		}
	}
	
	/* Check that the options selected are mutually consistent */
	
	if (start_date && start_dist > 0.0) {
		fprintf(stderr, "%s: ERROR: Cannot specify both start time AND start distance\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (stop_date && stop_dist < DBL_MAX) {
		fprintf(stderr, "%s: ERROR: Cannot specify both stop time AND stop distance\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (east < west || south > north) {
		fprintf(stderr, "%s: ERROR: Region set incorrectly\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (start_date && GMT_verify_expectations (GMT_IS_ABSTIME, GMT_scanf (start_date, GMT_IS_ABSTIME, &start_time), start_date)) {
		fprintf (stderr, "%s: ERROR: Start time (%s) in wrong format\n", GMT_program, start_date);
		exit (EXIT_FAILURE);
	}
	if (stop_date && GMT_verify_expectations (GMT_IS_ABSTIME, GMT_scanf (stop_date, GMT_IS_ABSTIME, &stop_time), stop_date)) {
		fprintf (stderr, "%s: ERROR: Stop time (%s) in wrong format\n", GMT_program, stop_date);
		exit (EXIT_FAILURE);
	}
	if (start_time > stop_time) {
		fprintf(stderr, "%s: ERROR: Start time exceeds stop time!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (start_dist > stop_dist) {
		fprintf(stderr, "%s: ERROR: Start distance exceeds stop distance!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (dist_gap && gap_d <= 0.0) {
		fprintf(stderr, "%s: ERROR: -Gd: Must specify a positive gap distance in km\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (time_gap && gap_t <= 0.0) {
		fprintf(stderr, "%s: ERROR: -Gt: Must specify a positive gap distance in minutes\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	
	if (GMT_give_synopsis_and_exit || argc == 1) {	/* Display usage */
		fprintf (stderr, "mgd77track %s - Plot track-line map of MGD77 cruises\n\n", MGD77_VERSION);
		fprintf (stderr, "usage: mgd77track cruise(s) %s %s [-A[c][size]][,<inc><unit>] [%s]\n", GMT_Rgeo_OPT, GMT_J_OPT, GMT_B_OPT);
		fprintf (stderr, "\t[-Cf|g|e] [-Da<startdate>] [-Db<stopdate>] [-F] [-Gt|d<gap>] [-I<code>] [-K] [-L<trackticks>] [-N] [-O] [-P] [-Sa<startdist>[unit]]\n");
		fprintf (stderr, "\t[-Sb<stopdist>[unit]] [-TT|t|d<ms,mc,mfs,mf,mfc>] [%s] [-V] [-W<pen>] [%s]\n", GMT_U_OPT, GMT_X_OPT);
		fprintf (stderr, "\t[%s] [%s]\n\n", GMT_Y_OPT, GMT_c_OPT);
      
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
              
		MGD77_Cruise_Explain ();
		GMT_explain_option ('J');
		GMT_explain_option ('R');
		fprintf (stderr, "	OPTIONS:\n\n");
		fprintf(stderr, "	-A will Annotate legs when they enter the grid. Append c for cruise ID [Default is file prefix]\n");
		fprintf(stderr, "	   <size> is optional text size in points [9].  The font used is controlled by LABEL_FONT\n");
		fprintf(stderr, "	   Optionally, append ,<inc>[unit] to place label every <inc> units apart. <unit> may be\n");
		fprintf(stderr, "	   k (km), n (nautical miles) or d (days), h (hours).\n");
		GMT_explain_option ('B');
		fprintf (stderr,"	-C Select procedure for along-track distance calculations:\n");
		fprintf (stderr,"	   f Flat Earth\n");
		fprintf (stderr,"	   g Great circle [Default]\n");
		fprintf (stderr,"	   e Ellipsoidal (geodesic) using current ellipsoid\n");
		fprintf (stderr, "	-Da plots from <startdate> (given as yyyy-mm-ddT[hh:mm:ss]) [Start of cruise]\n");
		fprintf (stderr, "	-Db plots up to <stopdate> (given as yyyy-mm-ddT[hh:mm:ss]) [End of cruise]\n");
		fprintf(stderr,"	-F Do NOT apply bitflags to MGD77+ cruises [Default applies error flags stored in the file]\n");
		fprintf(stderr,"	-G Consider point separations exceeding d<gap> (km) or t<gap> (minutes) to indicate a gap (do not draw) [0]\n");
		fprintf(stderr,"	-I Ignore certain data file formats from consideration. Append combination of act to ignore\n");
		fprintf (stderr, "	   (a) MGD77 ASCII, (c) MGD77+ netCDF, (m) MGD77T ASCII, or (t) plain table files. [Default ignores none]\n");
		GMT_explain_option ('K');
		fprintf (stderr, "	-L puts time/distance log marks on the track. E.g. a500ka24ht6h means (a)nnotate\n");
		fprintf (stderr, "	   every 500 km (k) and 24 h(ours), with (t)ickmarks every 500 km and 6 (h)ours\n");
		fprintf (stderr, "	   Units of n(autical miles) and d(ays) are also recognized.\n");
		fprintf (stderr, "	-N Do Not clip leg name annotation that fall outside map border [Default will clip]\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "	-Sa plots from <startdist> (in m; append k, m, or n) [Start of the cruise]\n");
		fprintf (stderr, "	-Sb plots up to <stopdist> (in m; append k, m, or n) [End of the cruise]\n");
		fprintf (stderr, "	-T sets attributes of marker items. Append T for new day marker, t for same\n");
		fprintf (stderr, "	   day marker, and d for distance marker.  Then, append 5 comma-separated items:\n");
		fprintf (stderr, "	   <markersize>[unit],<markercolor>,<markerfontsize,<markerfont>,<markerfontcolor>\n");
		fprintf (stderr, "	   Default settings for the three marker types are:\n");
		fprintf (stderr, "	     -TT%g,black,%g,%ld,black\n", marker[MGD77TRACK_MARK_NEWDAY].marker_size, marker[MGD77TRACK_MARK_NEWDAY].font_size, marker[MGD77TRACK_MARK_NEWDAY].font_no);
		fprintf (stderr, "	     -Tt%g,white,%g,%ld,black\n", marker[MGD77TRACK_MARK_SAMEDAY].marker_size, marker[MGD77TRACK_MARK_SAMEDAY].font_size, marker[MGD77TRACK_MARK_SAMEDAY].font_no);
		fprintf (stderr, "	     -Td%g,black,%g,%ld,black\n", marker[MGD77TRACK_MARK_DIST].marker_size, marker[MGD77TRACK_MARK_DIST].font_size, marker[MGD77TRACK_MARK_DIST].font_no);
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		fprintf (stderr, "	-W sets track pen attributes [width = %gp, color = (%d/%d/%d), texture = solid line].\n", 
			pen.width, pen.rgb[0], pen.rgb[1], pen.rgb[2]);
		GMT_explain_option ('X');
		GMT_explain_option ('c');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	n_paths = MGD77_Path_Expand (&M, argv, argc, &list);	/* Get list of requested IDs */

	if (n_paths == 0) {
		fprintf(stderr, "%s: ERROR: No cruises given\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	use = (M.original) ? MGD77_ORIG : MGD77_REVISED;

	if (east < west) west -= 360.0;
		
	GMT_err_fail (GMT_map_setup (west, east, south, north), "");
	
	GMT_plotinit (argc, argv);
	
	GMT_map_clip_on (GMT_no_rgb, 3);
	GMT_setpen (&pen);
	both = (info[ANNOT].annot_int_time && info[ANNOT].annot_int_dist);
	
	if (no_clip) {
		cruise_id = (struct MGD77TRACK_LEG_ANNOT *) GMT_memory (VNULL, (size_t)n_alloc_c, sizeof (struct MGD77TRACK_LEG_ANNOT), "mgd77track");
	}
	MGD77_Select_Columns ("time,lon,lat", &M, MGD77_SET_ALLEXACT);	/* This sets up which columns to read */

	for (argno = 0; argno < n_paths; argno++) {		/* Process each ID */
	
		if (MGD77_Open_File (list[argno], &M, MGD77_READ_MODE)) continue;

		if (gmtdefs.verbose) fprintf (stderr, "%s: Now processing cruise %s\n", GMT_program, list[argno]);
		
		if (MGD77_Read_Header_Record (list[argno], &M, &D.H)) {
			fprintf (stderr, "%s: Error reading header sequence for cruise %s\n", GMT_program, list[argno]);
			exit (EXIT_FAILURE);
		}
		rec = 0;
		last_julian = -1;
		
		if (abs(annot_name) == 2) {	/* Use MGD77 cruise ID */
			strcpy (name, D.H.mgd77[use]->Survey_Identifier);
		}
		else {			/* Use file name prefix */
			strcpy (name, list[argno]);
			for (i = 0; i < (int)strlen (name); i++) if (name[i] == '.') name[i] = '\0';
		}
	
		/* Start reading data from file */
	
		if (MGD77_Read_Data (list[argno], &M, &D)) {
			fprintf (stderr, "%s: Error reading data set for cruise %s\n", GMT_program, list[argno]);
			exit (EXIT_FAILURE);
		}
		MGD77_Close_File (&M);
		track_time = (double*)D.values[0];
		lon = (double*)D.values[1];
		lat = (double*)D.values[2];
		GMT_err_fail (GMT_distances (lon, lat, D.H.n_records, 1.0, dist_flag, &track_dist), "");	/* Work internally in meters */
		for (rec = 0; rec < D.H.n_records && bad_coordinates(lon[rec], lat[rec]) && track_time[rec] < start_time && track_dist[rec] < start_dist; rec++);	/* Find first record of interest */
		first_rec = rec;
		for (rec = D.H.n_records - 1; rec && track_time[rec] > stop_time && bad_coordinates(lon[rec], lat[rec]) && track_dist[rec] > stop_dist; rec--);	/* Find last record of interest */
		last_rec = rec;
		if (gmtdefs.verbose) fprintf (stderr, "mgd77track: Plotting %s [%s]\n", list[argno], D.H.mgd77[use]->Survey_Identifier);
		sprintf (comment, "Tracking %s", list[argno]);
		ps_comment (comment);
		
		/* First draw the track line, clip segments outside the area */
		
		if (dist_gap || time_gap) {
			GMT_LONG start, stop;
			start = first_rec;
			while (start < last_rec && ((dist_gap && (track_dist[start+1] - track_dist[start]) > gap_d) || (time_gap && (track_time[start+1] - track_time[start]) > gap_t))) {	/* First start of first segment */
				lon[start] = GMT_d_NaN;	/* Flag to make sure we do not plot this gap later */
				start++;
			}
			while (start <= last_rec) {
				stop = start;
				while (stop < last_rec && ((dist_gap && (track_dist[stop+1] - track_dist[stop]) < gap_d) || (time_gap && (track_time[stop+1] - track_time[stop]) < gap_t))) stop++;	/* stop will be last point in segment */
				GMT_n_plot = GMT_geo_to_xy_line (&lon[start], &lat[start], stop-start+1);
				GMT_plot_line (GMT_x_plot, GMT_y_plot, GMT_pen, GMT_n_plot);
				start = stop + 1;
				while (start < last_rec && ((dist_gap && (track_dist[start+1] - track_dist[start]) > gap_d) || (time_gap && (track_time[start+1] - track_time[start]) > gap_t))) {	/* First start of first segment */
					lon[start] = GMT_d_NaN;	/* Flag to make sure we do not plot this gap later */
					start++;
				}
			}
		}
		else {	/* Plot the whole shabang */
			GMT_n_plot = GMT_geo_to_xy_line (lon, lat, D.H.n_records);
			GMT_plot_line (GMT_x_plot, GMT_y_plot, GMT_pen, GMT_n_plot);
		}

		first = TRUE;
		for (rec = first_rec; rec <= last_rec; rec++) {
			if (bad_coordinates(lon[rec], lat[rec]) || GMT_outside (lon[rec], lat[rec])) {
				first = TRUE;
				continue;
			}
			GMT_geo_to_xy (lon[rec], lat[rec], &x, &y);
			if (first) {
				if (annot_name > 0) {
					c_angle = heading (rec, lon, lat, D.H.n_records);
					if (no_clip) {	/* Keep these in a list to plot after clipping is turned off */
						cruise_id[n_id].x = x;
						cruise_id[n_id].y = y;
						cruise_id[n_id].lon = lon[rec];
						cruise_id[n_id].lat = lat[rec];
						cruise_id[n_id].angle = c_angle;

						strcpy (cruise_id[n_id].text, name);
						n_id++;
						if (n_id == n_alloc_c) {
							n_alloc_c <<= 1;
							cruise_id = (struct MGD77TRACK_LEG_ANNOT *) GMT_memory ((void *)cruise_id, (size_t)n_alloc_c, sizeof (struct MGD77TRACK_LEG_ANNOT), "mgd77track");
						}
					}
					else
						annot_legname (x, y, lon[rec], lat[rec], c_angle, name, GMT_u2u[GMT_INCH][GMT_PT] * 1.25 * annotsize);
				}
				first = FALSE;
				for (i = 0; i < 2; i++) {
					if (info[i].annot_int_dist > 0) annot_dist[i] = (track_dist[rec] / info[i].annot_int_dist + 1) * info[i].annot_int_dist;
					if (info[i].tick_int_dist > 0) tick_dist[i] = (track_dist[rec] / info[i].tick_int_dist + 1) * info[i].tick_int_dist;
					if (info[i].annot_int_time > 0) annot_time[i] = ceil (track_time[rec] / info[i].annot_int_time) * info[i].annot_int_time;
					if (info[i].tick_int_time > 0) tick_time[i] = ceil (track_time[rec] / info[i].tick_int_time) * info[i].tick_int_time;
				}
			}
			
			/* See if we need to annotate/tick the trackline for time/km and/or ID marks */
			
			for (i = 0; i < 2; i++) {
				if (info[i].annot_int_time && (track_time[rec] >= annot_time[i])) {
					annot_time[i] += info[i].annot_int_time;
					annot_tick[i] = 1;
				}
				if (info[i].annot_int_dist && (track_dist[rec] >= annot_dist[i])) {
					annot_dist[i] += info[i].annot_int_dist;
					annot_tick[i] += 2;
				}
				if (info[i].tick_int_time && (track_time[rec] >= tick_time[i])) {
					tick_time[i] += info[i].tick_int_time;
					draw_tick[i] = 1;
				}
				if (info[i].tick_int_dist && (track_dist[rec] >= tick_dist[i])) {
					tick_dist[i] += info[i].tick_int_dist;
					draw_tick[i] += 2;
				}
			}
			if (annot_tick[ANNOT]) {
				angle = heading (rec, lon, lat, D.H.n_records);
				if (angle < 0.0)
					angle += 90.0;
				else
					angle -= 90.0;
				if (annot_tick[ANNOT] & 1) {	/* Time mark */
					GMT_gcal_from_dt (annot_time[ANNOT], &calendar);			/* Convert t to a complete calendar structure */
					GMT_format_calendar (date, clock, &GMT_plot_calclock.date, &GMT_plot_calclock.clock, FALSE, 1, annot_time[ANNOT]);
					this_julian = calendar.day_y;
					if (this_julian != last_julian) {
						mrk = MGD77TRACK_MARK_NEWDAY;
						sprintf (label, "%s+%s", date, clock);
						ps_circle (x, y, marker[mrk].marker_size, marker[mrk].s.rgb, TRUE);
					}
					else {
						mrk = MGD77TRACK_MARK_SAMEDAY;
						sprintf (label, "+%s", clock);
						ps_circle (x, y, marker[mrk].marker_size, marker[mrk].s.rgb, TRUE);
					}
					ps_setfont (marker[mrk].font_no);
					ps_setpaint (marker[mrk].f.rgb);
					plot_x = x;	plot_y = y;
					GMT_smart_justify (5, angle, 0.5 * marker[mrk].font_size, 0.5 * marker[mrk].font_size, &plot_x, &plot_y);
					ps_text (plot_x, plot_y, GMT_u2u[GMT_INCH][GMT_PT] * marker[mrk].font_size, label, angle, 5, 0);
					last_julian = calendar.day_y;
				}
				if (annot_tick[ANNOT] & 2) {	/* Distance mark */
					mrk = MGD77TRACK_MARK_DIST;
					sprintf (label, "%d km  ", (int)((annot_dist[ANNOT] - info[ANNOT].annot_int_dist) * factor));
					ps_square (x, y, marker[mrk].marker_size, marker[mrk].s.rgb, TRUE);
					ps_setpaint (marker[mrk].f.rgb);
					plot_x = x;	plot_y = y;
					GMT_smart_justify (7, angle, 0.5 * marker[mrk].font_size, 0.5 * marker[mrk].font_size, &plot_x, &plot_y);
					ps_text (plot_x, plot_y, GMT_u2u[GMT_INCH][GMT_PT] * marker[mrk].font_size, label, angle, 7, 0);
				}
			}
			if (both && !(annot_tick[ANNOT] & 1) && (draw_tick[ANNOT] & 1)) {
				mrk = (this_julian != last_julian) ? MGD77TRACK_MARK_NEWDAY : MGD77TRACK_MARK_SAMEDAY;
				ps_circle (x, y, marker[mrk].marker_size, marker[mrk].s.rgb, TRUE);
			}
			if (both && !(annot_tick[ANNOT] & 2) && (draw_tick[ANNOT] & 2)) {
				mrk = (this_julian != last_julian) ? MGD77TRACK_MARK_NEWDAY : MGD77TRACK_MARK_SAMEDAY;
				ps_square (x, y, marker[mrk].marker_size, marker[mrk].s.rgb, TRUE);
			}
			if (draw_tick[ANNOT]) {
				mrk = MGD77TRACK_MARK_DIST;
				ps_setpaint (marker[mrk].s.rgb);
				ps_cross (x, y, marker[mrk].marker_size);
			}
			if (annot_tick[ANNOT] || draw_tick[ANNOT]) annot_tick[ANNOT] = draw_tick[ANNOT] = FALSE;
			if (annot_tick[LABEL]) {
				angle = heading (rec, lon, lat, D.H.n_records);
				if (angle < 0.0)
					angle += 90.0;
				else
					angle -= 90.0;
				annot_legname (x, y, lon[rec], lat[rec], angle, name, GMT_u2u[GMT_INCH][GMT_PT] * 1.25 * annotsize);
				annot_tick[LABEL] = FALSE;
			}
		}
		GMT_free ((void *)track_dist);
		n_cruises++;
	}
		
	GMT_map_clip_off();

	GMT_map_basemap ();
	
	if (annot_name > 0 && no_clip) {	/* Plot leg names after clipping is terminated ( see -N) */
		int id;
		double size;
		size = GMT_u2u[GMT_INCH][GMT_PT] * 1.25 * annotsize;
		for (id = 0; id < n_id; id++) annot_legname (cruise_id[id].x, cruise_id[id].y, cruise_id[id].lon, cruise_id[id].lat, cruise_id[id].angle, cruise_id[id].text, size);
		GMT_free ((void *)cruise_id);
	}

	GMT_plotend ();
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Plotted %ld cruises\n", GMT_program, n_cruises);

	MGD77_Path_Free (n_paths, list);
	MGD77_end (&M);
	
	exit (EXIT_SUCCESS);
}

double heading (GMT_LONG rec, double *lon, double *lat, GMT_LONG n_records)
{
	GMT_LONG i1, i2, j;
	double angle, x1, x0, y1, y0;
	double sum_x2 = 0.0, sum_xy = 0.0, sum_y2 = 0.0, dx, dy;
	
	i1 = rec - 10;
	if (i1 < 0) i1 = 0;
	i2 = i1 + 10;
	if (i2 > (n_records-1)) i2 = n_records - 1;
	GMT_geo_to_xy (lon[rec], lat[rec], &x0, &y0);
	for (j = i1; j <= i2; j++) {	/* L2 fit for slope over this range of points */
		GMT_geo_to_xy (lon[j], lat[j], &x1, &y1);
		dx = x1 - x0;
		dy = y1 - y0;
		sum_x2 += dx * dx;
		sum_y2 += dy * dy;
		sum_xy += dx * dy;
	}
	if (sum_y2 < GMT_CONV_LIMIT)	/* Line is horizontal */
		angle = 0.0;
	else if (sum_x2 < GMT_CONV_LIMIT)	/* Line is vertical */
		angle = 90.0;
	else
		angle = (GMT_IS_ZERO (sum_xy)) ? 90.0 : d_atan2d (sum_xy, sum_x2);
	if (angle > 90.0)
		angle -= 180;
	else if (angle < -90.0)
		angle += 180.0;
	return (angle);
}

int get_annotinfo (char *args, struct MGD77TRACK_ANNOT *info)
{
	int i1, i2, error = FALSE;
	int flag1, flag2, type;
	double value;
	
	info->annot_int_dist = info->tick_int_dist = 0;
	info->annot_int_time = info->tick_int_time = 0;

	i1 = 0;
	while (args[i1]) {
		flag1 = 'a';
		if (isalpha ((int)args[i1])) {
			flag1 = args[i1];
			if (isupper(flag1)) flag1 = tolower (flag1);
			i1++;
		}
		i2 = i1;
		while (args[i2] && strchr ("athkmnd", (int)args[i2]) == NULL) i2++;
		value = atof (&args[i1]);
		flag2 = args[i2];
		if (isupper(flag2)) flag2 = tolower (flag2);
		if (flag2 == 'd') {		/* Days */
			value *= GMT_DAY2SEC_F;
			type = 't';
		}
		else if (flag2 == 'h') {		/* Hours */
			value *= GMT_HR2SEC_F;
			type = 't';
		}
		else if (flag2 == 'k') {	/* kilometers */
			value *= 1000;
			type = 'd';
		}
		else if (flag2 == 'n') {	/* Nautical miles */
			value *= MGD77_METERS_PER_NM;
			type = 'd';
		}
		else if (flag2 == 'm') {	/* Minutes */
			value *= GMT_MIN2SEC_F;
			type = 't';
		}
		else				/* Default is seconds */
			type = 't';
		i2++;
		if (flag1 == 'a') {	/* Annotation interval */
			if (type == 'd')	/* Distance */
				info->annot_int_dist = (int)value;
			else
				info->annot_int_time = (int)value;
		}
		else {			/* Tickmark interval */
			if (type == 'd')	/* Distance */
				info->tick_int_dist = (int)value;
			else
				info->tick_int_time = (int)value;
		}
		i1 = i2;
	}
	if (info->annot_int_dist <= 0 && info->tick_int_dist <= 0 && info->annot_int_time <= 0 && info->tick_int_time <= 0)
		error = TRUE;
	if (info->annot_int_dist <= 0)
		info->annot_int_dist = info->tick_int_dist;
	else if (info->tick_int_dist <= 0)
		info->tick_int_dist = info->annot_int_dist;
	if (info->annot_int_time <= 0)
		info->annot_int_time = info->tick_int_time;
	else if (info->tick_int_time <= 0)
		info->tick_int_time = info->annot_int_time;
	return (error);
}

void annot_legname (double x, double y, double lon, double lat, double angle, char *text, double size)
{
	int just;
	
	if (lat < project_info.s)
		just = (angle >= 0.0) ? 1 : 3;
	else if (lat > project_info.n)
		just = (angle >= 0.0) ? 11 : 9;
	else if (lon < project_info.w)
		just = (angle >= 0.0) ? 9 : 1;
	else
		just = (angle >= 0.0) ? 3 : 11;
	ps_setfont (gmtdefs.label_font);
	ps_setpaint (gmtdefs.basemap_frame_rgb);
	GMT_smart_justify (just, angle, GMT_u2u[GMT_PT][GMT_INCH] * 0.15 * size, GMT_u2u[GMT_PT][GMT_INCH] * 0.15 * size, &x, &y);
	ps_text (x, y, size, text, angle, just, 0);
}

GMT_LONG bad_coordinates (double lon, double lat) {
	return (GMT_is_dnan (lon) || GMT_is_dnan (lat));
}
