/*--------------------------------------------------------------------
 *	$Id: backtracker.c 9923 2012-12-18 20:45:53Z pwessel $
 *
 *   Copyright (c) 1999-2013 by P. Wessel
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; version 2 or any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   Contact info: www.soest.hawaii.edu/pwessel
 *--------------------------------------------------------------------*/
/*
 * Program for moving points along small circles on a sphere given a
 * set of plate motion stage (Euler) poles.
 * backtracker can move a point forward or backward in time.
 * It can do so either along flowlines or hotspot tracks.
 * It can move point to final position or generate a track between
 * starting point and final position.  The output track is a GMT
 * multisegment file and can be plotted with psxy -M.
 *
 * Author:	Paul Wessel, SOEST, Univ. of Hawaii, Honolulu, HI, USA
 * Date:	29-DEC-1999
 * Version:	1.1
 *
 *-------------------------------------------------------------------------
 * The ASCII Euler file must have following format:
 *
 * 1. Any number of comment lines starting with # in first column
 * 2. Any number of blank lines (just carriage return, no spaces)
 * 2. Any number of stage pole records which each have the format:
 *    lon(deg)  lat(deg)  tstart(Ma)  tstop(Ma)  ccw-angle(deg)
 * 3. stage records must go from oldest to youngest rotation
 * 4. Note tstart is larger (older) that tstop for each record
 * 5. No gaps allowed: tstart must equal the previous records tstop
 *
 * Example: Duncan & Clague [1985] Pacific-Hotspot rotations:
 *
 * # Time in Ma, angles in degrees
 * # lon  lat	tstart	tend	ccw-angle
 * 165     85	150	100	24.0
 * 284     36	100	74	15.0
 * 265     22	74	65	7.5
 * 253     17	65	42	14.0
 * 285     68	42	0	34.0
 *
 * ASCII point location file(s) must have the following format:
 *
 * 1. Any number of comment lines starting with # in first column
 * 2. Any number of blank lines (just carriage return, no spaces)
 * 3. For special header records, see -H
 * 4. Any number of data recordswhich each have the format:
 *    lon lat age    (or lat lon age, see -: option). age in Ma.
 *
 * Binary files cannot have header records, and data fields must all be
 * either single or double precision (see -bi option)
 */

#include "spotter.h"

int main (int argc, char **argv)
{
	struct EULER *p = NULL;			/* Pointer to array of stage poles */

	GMT_LONG n_points;			/* Number of data points read */
	GMT_LONG n_chunk;			/* Total length or array returned by libeuler functions */
	GMT_LONG n_track;			/* Number of points in a track segment */
	GMT_LONG n_stages = 0;			/* Number of stage poles */
	GMT_LONG n_segments;			/* Number of path segments written out */
	GMT_LONG n_skipped = 0;			/* Number of points skipped because t < 0 */
	GMT_LONG track_limit = 0;		/* Flag for limiting the output track extent [no limit] */
	GMT_LONG n_m = 0;			/* Number of points in hotspot drift table, if used */
	GMT_LONG n_args, nm_alloc;
	GMT_LONG n_files = 0, fno;
	GMT_LONG n_fields, n_expected_fields;
	GMT_LONG n_read = 0;
	GMT_LONG n_out;
	GMT_LONG i, j, k;			/* Misc. counters */

	GMT_LONG error = FALSE;		/* Set to TRUE if arguments are inconsistent */
	GMT_LONG flowline = FALSE;	/* TRUE means we want flowlines, FALSE we want hotspot tracks */
	GMT_LONG forward = FALSE;	/* TRUE we go FROM hotspot to seamount, FALSE is reverse */
	GMT_LONG make_path = FALSE;	/* TRUE means create continuous path, FALSE works on discrete points */
	GMT_LONG seg_files = FALSE;	/* TRUE will write individual files for each segment */
	GMT_LONG finite = FALSE;	/* TRUE if stage pole file contains finite rotation poles instead */
	GMT_LONG fix_all = FALSE;	/* TRUE is all data shall be of a fixed age */
	GMT_LONG stage_id = FALSE;
	GMT_LONG first = TRUE;
	GMT_LONG done = FALSE;
	GMT_LONG nofile = TRUE;
	GMT_LONG confidence = FALSE;
	GMT_LONG finite_out = FALSE;	/* FALSE unless -W is set */
	GMT_LONG single_rotation = FALSE;	/* FALSE unless -e is set */
	GMT_LONG do_track = FALSE;	/* TRUE if -L is used */

	double d_km = 0.0;		/* Step interval along calculated tracks */
	double t_zero = 0.0;		/* Current age in Ma */
	double t_fix = 0.0;		/* Fixed age in Ma */
	double upper_age = 0.0;		/* Extend oldest age back to this time, in Ma */
	double *c = NULL;		/* Array of track chunks returned by libeuler routines */
	double lon, lat;		/* Seamounts location in decimal degrees */
	double *dlon = NULL, *dlat = NULL, *dt = NULL;	/* Drifting hotspot locations (t) in decimal degrees */
	double age, t;			/* Age of seamount, in Ma */
	double *in = NULL, *out = NULL;	/* i/o arrays used by GMT */
	double t_low, t_high;		/* upper/lower age/stage to ouput for track segments */
	double t_end;
	double rot_lon, rot_lat, rot_w;	/* A finite rotation */
	double R[3][3];			/* Rotation matrix */
	double x[3], y[3];		/* Two 3-D unit vectors */

	char *euler_file = CNULL;	/* Name pointer for file with stage poles */
	char *file_stem = CNULL;	/* Name pointer for file stem for individual track files */
	char *drift_file = CNULL;	/* Name pointer for file with hotspot motion */
	char buffer[BUFSIZ];		/* Input buffer for reading data */
	char type[50];			/* What kind of line (flowline or hotspot track) */
	char dir[8];			/* From or To */
	char conf_flag = 0;		/* t for time, a for angle, blank for nothing */
	char txt_a[GMT_TEXT_LEN], txt_b[GMT_TEXT_LEN];
	char *not_used = NULL;

	FILE *fp = NULL;		/* File pointer for input data */
	FILE *fpo = NULL;		/* File pointer for output data */

	PFI spot_func = NULL;			/* Pointer to the required forth/back track function */

	out = (double *)NULL;

	argc = (int)GMT_begin (argc, argv);

#ifdef DEBUG
	if (gmtdefs.verbose) fprintf (stderr, "DEBUG mode\n");
#endif
	GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
	GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;


	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'H':
				case 'M':
				case 'V':
				case ':':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], NULL, NULL, NULL, NULL);
					break;

				/* Supplemental parameters */

				case 'b':
					error += GMT_parse_b_option (&argv[i][2]);
					break;

				case 'A':	/* Output only an age-limited segment of the track */
					if (argv[i][2]) {	/* Gave specific limits for all input points */
						sscanf (&argv[i][2], "%lf/%lf", &t_low, &t_high);
						track_limit = 1;
					}
					else {	/* Limits for each input point given in columns 4 and 5 */
						track_limit = 2;
					}
					break;

				case 'C':	/* Use total reconstruction poles */
					finite = TRUE;
					break;

				case 'D':	/* Specify in which direction we should project */
					switch (argv[i][2]) {
						case 'B':	/* Go FROM hotspot TO seamount */
						case 'b':
							forward = FALSE;
							break;
						case 'F':	/* Go FROM seamount TO hotspot */
						case 'f':
							forward = TRUE;
							break;
						default:
							error++;
							fprintf (stderr, "%s ERROR Option -D: Append b or f\n", GMT_program);
							break;
					}
					break;

				case 'L':	/* Specify what kind of track to project */
					do_track = TRUE;
					switch (argv[i][2]) {
						case 'F':	/* Calculate flowlines */
							stage_id = TRUE;
						case 'f':
							flowline = TRUE;
							break;
						case 'B':	/* Calculate hotspot tracks */
							stage_id = TRUE;
						case 'b':
							flowline = FALSE;
							break;
						default:
							error++;
							fprintf (stderr, "%s ERROR Option -L: Append f or b\n", GMT_program);
							break;
					}
					d_km = (argv[i][3]) ? atof (&argv[i][3]) : -1.0;
					break;

				case 'E':	/* File with stage poles */
					euler_file  = &argv[i][2];
					break;

				case 'e':	/* Apply a fixed total reconstruction rotation to all input points  */
					single_rotation  = TRUE;
					sscanf (&argv[i][2], "%[^/]/%[^/]/%lg", txt_a, txt_b, &rot_w);
					error += GMT_verify_expectations (GMT_io.in_col_type[0], GMT_scanf_arg (txt_a, GMT_io.in_col_type[0], &rot_lon), txt_a);
					error += GMT_verify_expectations (GMT_io.in_col_type[1], GMT_scanf_arg (txt_b, GMT_io.in_col_type[1], &rot_lat), txt_b);
					break;

				case 'F':	/* File with hotspot motions */
					drift_file  = &argv[i][2];
					break;

				case 'Q':	/* Fixed age for all points */
					t_fix = atof (&argv[i][2]);
					fix_all = TRUE;
					break;

				case 'S':	/* Set file stem for individual output files */
					if (argv[i][2]) {
						file_stem = &argv[i][2];
						seg_files = TRUE;
					}
					else {
						fprintf (stderr, "%s ERROR Option -S: Append a file stem\n", GMT_program);
						error++;
					}
					break;

				case 'T':	/* Current age [0 Ma] */
					t_zero = atof (&argv[i][2]);
					break;

				case 'W':	/* Report confidence ellipses */
					confidence = finite_out = TRUE;
					conf_flag = (argv[i][2]) ? argv[i][2] : 0;
					break;

				case 'N':	/* Extend oldest stage back to this time [no extension] */
					upper_age = atof (&argv[i][2]);
					break;

				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			n_files++;
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "%s %s - Forward and backward flowlines and hotspot tracks\n\n", GMT_program, SPOTTER_VERSION);
		fprintf (stderr, "usage: %s [infile(s)] -E<euler.d> OR -eplon/plat/prot [-A[young/old]] [-C] [-Df|b] [-F<driftfile] [%s] [-Lf|b<d_km>]\n", GMT_program, GMT_H_OPT);
		fprintf (stderr, "\t[-N<upper_age>] [-Q<t_fix>] [-S<stem>] [-T<t_zero>] [-V] [-W] [%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_b_OPT, GMT_m_OPT);
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\tinfiles (in ASCII or binary) has 3 or more columns.  If no file(s) is given, standard input is read.\n");
		fprintf (stderr, "\tFirst 3 columns must have lon, lat (or lat, lon, see -:) and age (Ma).\n");
		fprintf (stderr, "\t-E specifies the rotations to be used (see man page for format).\n");
		fprintf (stderr, "\t-e Alternatively, specify a single finite rotation (in degrees) to be applied to all input points.\n");
		fprintf (stderr, "\tOPTIONS:\n\n");
		fprintf (stderr, "\t-A Output tracks for ages (or stages, see -L) between young and old [Default is entire track].\n");
		fprintf (stderr, "\t   If no limit is given, then each seamount should have their limits in columns 4 and 5 instead.\n");
		fprintf (stderr, "\t   Only applicable in conjunction with the -L option.\n");
		fprintf (stderr, "\t-C The file given with -E contains total reconstruction poles [Default is stage poles].\n");
		fprintf (stderr, "\t-Db Backtrack mode: move forward in time (from older to younger positions) [Default].\n");
		fprintf (stderr, "\t-Df Flowline mode: move backward in time (from younger to older positions).\n");
		fprintf (stderr, "\t-F File with lon, lat, time records describing motion of hotspot responsible for\n");
		fprintf (stderr, "\t   the seamount/path we are concerned with [fixed hotspots].  If given, then the\n");
		fprintf (stderr, "\t   input lon, lat is replaced by the position of the drifting hotspot at the given age.\n");
		fprintf (stderr, "\t   If -F is used the <d_km> in -L is assumed to be point spacing in Ma.\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-Lb Compute hotspot tracks sampled every <d_km> interval [Default projects single points].\n");
		fprintf (stderr, "\t-Lf Compute flowline for seamounts of unknown but maximum age [Default projects single points].\n");
		fprintf (stderr, "\t    If no <d_km> is given, the start/stop points for each stage are returned.\n");
		fprintf (stderr, "\t    If B and F is used instead, stage id is returned as z-value [Default is predicted ages].\n");
		fprintf (stderr, "\t-N extends earliest stage pole back to <upper_age> [no extension].\n");
		fprintf (stderr, "\t-Q Assigned a fixed age to all input points.\n");
		fprintf (stderr, "\t-S writes tracks to individual files <stem>.# instead of to stdout (requires -L).\n");
		fprintf (stderr, "\t-T sets the current age in Ma [0].\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t-W Return projected point and confidence ellipse for the finite rotation.\n");
		fprintf (stderr, "\t   The input time must exactly match the age of a finite rotation or else we skip the point.\n");
		fprintf (stderr, "\t   Output record will be lon,lat,az,major,minor.\n");
		fprintf (stderr, "\t   -Wt will output lon,lat,time,az,major,minor.\n");
		fprintf (stderr, "\t   -Wa will output lon,lat,angle,az,major,minor.\n");
		fprintf (stderr, "\t   Use -D to specify which direction to rotate [forward in time].\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 3 input columns\n");
		GMT_explain_option ('o');
		GMT_explain_option ('.');
		fprintf (stderr, "\t   Output produced by -Lf|t can be plotted with psxy using the -m option.\n");
		GMT_explain_option ('m');
		exit (EXIT_FAILURE);
	}

	if (confidence && do_track) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  -W cannot be set if -Lf or -Lb are set\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 3 + ((track_limit == 2) ? 2 : 0);
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 3) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 3 columns\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

	GMT_lat_swap_init ();	/* Initialize auxiliary latitude machinery */
	
	if (drift_file && (fp = fopen (drift_file, "r")) != NULL) {
		nm_alloc = GMT_SMALL_CHUNK;
		dlon = (double *) GMT_memory (VNULL, (size_t)nm_alloc, sizeof (double), GMT_program);
		dlat = (double *) GMT_memory (VNULL, (size_t)nm_alloc, sizeof (double), GMT_program);
		dt   = (double *) GMT_memory (VNULL, (size_t)nm_alloc, sizeof (double), GMT_program);
		k = 0;
		while (fgets (buffer, 512, fp) != NULL) {
			if (buffer[0] == '#' || buffer[0] == '\n') continue;
			sscanf (buffer, "%lf %lf %lf", &dlon[k], &dlat[k], &dt[k]);
			dlat[k] = GMT_lat_swap (dlat[k], GMT_LATSWAP_G2O);	/* Convert to geocentric */
			k++;
			if (k == nm_alloc) {
				nm_alloc <<= 1;
				dlon = (double *) GMT_memory ((void *)dlon, (size_t)nm_alloc, sizeof (double), GMT_program);
				dlat = (double *) GMT_memory ((void *)dlat, (size_t)nm_alloc, sizeof (double), GMT_program);
				dt   = (double *) GMT_memory ((void *)dt,   (size_t)nm_alloc, sizeof (double), GMT_program);
			}
		}
		GMT_fclose (fp);
		n_m = k;
		dlon = (double *) GMT_memory ((void *)dlon, (size_t)k, sizeof (double), GMT_program);
		dlat = (double *) GMT_memory ((void *)dlat, (size_t)k, sizeof (double), GMT_program);
		dt   = (double *) GMT_memory ((void *)dt,   (size_t)k, sizeof (double), GMT_program);
	}
		
	if (single_rotation) {	/* Get rotation matrix R */
		spotter_make_rot_matrix (rot_lon, rot_lat, rot_w, R);
	}
	else {	/* Load in the stage poles */
		n_stages = spotter_init (euler_file, &p, (int)flowline, finite, finite_out, &upper_age, gmtdefs.verbose);
		spot_func = ((flowline + forward) == 1) ? spotter_forthtrack : spotter_backtrack;

		if (fabs (d_km) > GMT_SMALL) {		/* User wants to interpolate tracks rather than project individual points */
			make_path = TRUE;
			(flowline) ? sprintf (type, "Flowline") : sprintf (type, "Hotspot track");
			(forward) ? sprintf (dir, "from") : sprintf (dir, "to");
		}
		else if (track_limit > 0) {	/* Limits on track set but track option not selected */
	                fprintf (stderr, "%s: GMT SYNTAX ERROR.  -A requires -L\n", GMT_program);
			exit (EXIT_FAILURE);
		}
	}

	fpo = GMT_stdout;
	n_out = (seg_files) ? 4 : 3;	/* Append smt id number as 4th column when individual files are requested */

	/* Read the seamount data from file or stdin */

	n_points = n_segments = 0;

	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;
	n_args = (argc > 1) ? argc : 2;

	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 3 + ((track_limit == 2) ? 2 : 0);
	if (fix_all && !GMT_io.ncol[GMT_IN]) n_expected_fields = 2;	/* Lon, lat only; use fixed t */
	if (confidence) n_out = 5 + !(conf_flag == 0);
	if (single_rotation && !GMT_io.binary[GMT_IN]) n_expected_fields = GMT_MAX_COLUMNS;	/* Allow any input for -F mode */
	if (do_track) GMT_io.multi_segments[GMT_OUT] = TRUE;	/* Since we are producing segments */

	for (fno = 1; !done && fno < n_args; fno++) {	/* Loop over input files, if any */
		if (!nofile && argv[fno][0] == '-') continue;

		if (nofile) {	/* Just read standard input */
			fp = GMT_stdin;
			done = TRUE;
			if (gmtdefs.verbose) fprintf (stderr, "%s: Reading from standard input\n", GMT_program);
		}
		else if ((fp = GMT_fopen (argv[fno], GMT_io.r_mode)) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[fno]);
			continue;
		}

		if (!nofile && gmtdefs.verbose) fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[fno]);

		if (GMT_io.io_header[GMT_IN]) {
			for (i = 0; i < GMT_io.n_header_recs; i++) {
				not_used = GMT_fgets (buffer, BUFSIZ, fp);
				if (first && !GMT_io.binary[GMT_OUT] && GMT_io.io_header[GMT_OUT]) GMT_fputs (buffer, GMT_stdout);
			}
			first = FALSE;
		}

		while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */
			n_read++;
			while ((GMT_io.status & GMT_IO_SEGMENT_HEADER) && !(GMT_io.status & GMT_IO_EOF)) {
				if (!make_path) GMT_write_segmentheader (GMT_stdout, n_expected_fields);
				n_fields = GMT_input (fp, &n_expected_fields, &in);
				n_read++;
			}
			if (GMT_io.status & GMT_IO_EOF) continue;

			if (GMT_io.status & GMT_IO_MISMATCH) {
				fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
				continue;
			}
			
			if (!out) out = (double *) GMT_memory (VNULL, (size_t)MAX(n_out, n_expected_fields), sizeof (double), GMT_program);

			in[GMT_Y] = GMT_lat_swap (in[GMT_Y], GMT_LATSWAP_G2O);	/* Convert to geocentric */

			if (single_rotation) {	/* Simple reconstruction, then exit */
				GMT_geo_to_cart (in[GMT_Y], in[GMT_X], x, TRUE);	/* Get x-vector */
				spotter_matrix_vect_mult (R, x, y);			/* Rotate the x-vector */
				GMT_cart_to_geo (&out[GMT_Y], &out[GMT_X], y, TRUE);	/* Recover lon lat representation; TRUE to get degrees */
				out[GMT_Y] = GMT_lat_swap (out[GMT_Y], GMT_LATSWAP_O2G);	/* Convert back to geodetic */
				memcpy ((void *)&out[GMT_Z], (void *)&in[GMT_Z], (n_fields - 2) * sizeof (double));
				GMT_output (fpo, n_fields, out);
				continue;
			}
			
			if (fix_all) in[2] = t_fix;
			if (in[2] < 0.0) {	/* Negative ages are flags for points to be skipped */
				n_skipped++;
				continue;
			}

			if (n_m) {	/* Must account for hotspot drift */
				GMT_intpol (dt, dlon, n_m, 1, &age, &lon, gmtdefs.interpolant);
				GMT_intpol (dt, dlat, n_m, 1, &age, &lat, gmtdefs.interpolant);
			}
			else {	/* Use intput location */
				lon = in[0];
				lat = in[1];
			}
			lon *= D2R;	lat *= D2R;
			if (track_limit) {
				if (track_limit == 2) t_low = in[3], t_high = in[4];
				age = t_high;	/* No point working more than necessary */
			}
			else
				age = in[2];

			if (age > upper_age) {	/* Points older than oldest stage cannot be used */
				fprintf (stderr, "%s: Seamount at line %ld has age (%g) > oldest stage (%g) (skipped)\n", GMT_program, n_read, in[2], upper_age);
				n_skipped++;
				continue;
			}

			if (make_path) {	/* Asked for paths, now write out several multiple segment tracks */
				if (seg_files) {
					sprintf (buffer, "%s.%ld", file_stem, n_points);
					if ((fpo = GMT_fopen (buffer, GMT_io.w_mode)) == NULL) {
						fprintf (stderr, "%s: Error - cannot create file %s\n", GMT_program, buffer);
						exit (EXIT_FAILURE);
					}
					out[3] = (double)n_points;	/* Put the seamount id number in 4th column */
				}
				sprintf (GMT_io.segment_header, "> %s %s %g %g\n", type, dir, in[0], in[1]);
				GMT_write_segmentheader (fpo, n_out);
				if (n_m) {	/* Must generate intermediate points in time */
					fprintf (stderr, "%s: Using drift file %s\n", GMT_program, drift_file);
					t = (track_limit) ? t_low : 0.0;
					t_end = (track_limit) ? t_high : age;
					GMT_intpol (dt, dlon, n_m, 1, &t, &lon, gmtdefs.interpolant);
					GMT_intpol (dt, dlat, n_m, 1, &t, &lat, gmtdefs.interpolant);
					lon *= D2R;	lat *= D2R;
					n_chunk = (*spot_func) (&lon, &lat, &t, 1, p, n_stages, 0.0, t_zero, TRUE + stage_id, NULL, &c);
					lat = GMT_lat_swap (lat * R2D, GMT_LATSWAP_O2G);	/* Convert back to geodetic */
					out[GMT_X] = lon * R2D;
					out[GMT_Y] = lat;
					out[GMT_Z] = t;
					GMT_output (fpo, n_out, out);
					t += d_km;	/* dt, actually */
					while (t < t_end) {
						GMT_intpol (dt, dlon, n_m, 1, &t, &lon, gmtdefs.interpolant);
						GMT_intpol (dt, dlat, n_m, 1, &t, &lat, gmtdefs.interpolant);
						lon *= D2R;	lat *= D2R;
						n_chunk = (*spot_func) (&lon, &lat, &t, 1, p, n_stages, 0.0, t_zero, TRUE + stage_id, NULL, &c);
						lat = GMT_lat_swap (lat * R2D, GMT_LATSWAP_O2G);	/* Convert back to geodetic */
						out[GMT_X] = lon * R2D;
						out[GMT_Y] = lat;
						out[GMT_Z] = t;
						GMT_output (fpo, n_out, out);
						t += d_km;	/* dt, actually */
					}
					GMT_intpol (dt, dlon, n_m, 1, &t_end, &lon, gmtdefs.interpolant);
					GMT_intpol (dt, dlat, n_m, 1, &t_end, &lat, gmtdefs.interpolant);
					lon *= D2R;	lat *= D2R;
					n_chunk = (*spot_func) (&lon, &lat, &t_end, 1, p, n_stages, 0.0, t_zero, TRUE + stage_id, NULL, &c);
					lat = GMT_lat_swap (lat * R2D, GMT_LATSWAP_O2G);	/* Convert back to geodetic */
					out[GMT_X] = lon * R2D;
					out[GMT_Y] = lat;
					out[GMT_Z] = t_end;
					GMT_output (fpo, n_out, out);
				}
				else {
					if (!confidence) n_chunk = (*spot_func) (&lon, &lat, &age, 1, p, n_stages, d_km, t_zero, TRUE + stage_id, NULL, &c);

					i = 0;
					n_track = irint (c[i++]);
					for (j = 0; j < n_track; j++, i += 3) {
						out[GMT_Z] = c[i+2];
						if (track_limit && (out[GMT_Z] < t_low || out[GMT_Z] > t_high)) continue;
						out[GMT_X] = c[i] * R2D;
						out[GMT_Y] = GMT_lat_swap (c[i+1] * R2D, GMT_LATSWAP_O2G);	/* Convert back to geodetic */
						GMT_output (fpo, n_out, out);
					}
				}
				if (seg_files) GMT_fclose (fpo);

				GMT_free ((void *)c);
			}
			else {	/* Just return the projected locations */
				if (confidence) {	/* Asked for confidence ellipses on reconstructed points */
					if (spotter_conf_ellipse (in[GMT_X], in[GMT_Y], age, p, n_stages, conf_flag, forward, out)) {
						fprintf (stderr, "%s: Confidence ellipses only for the age of rotations.  Point with age %g skipped\n", GMT_program, age);
						continue;
					}
				}
				else {
					n_chunk = (*spot_func) (&lon, &lat, &age, 1, p, n_stages, 0.0, t_zero, TRUE + stage_id, NULL, &c);
					out[GMT_X] = lon * R2D;	out[GMT_Y] = lat * R2D;
					for (k = 2; k < n_expected_fields; k++) out[k] = in[k];
				}
				out[GMT_Y] = GMT_lat_swap (out[GMT_Y], GMT_LATSWAP_O2G);	/* Convert back to geodetic */
				GMT_output (fpo, n_out, out);
			}

			n_points++;

		}
		if (fp != stdin) GMT_fclose (fp);
	}

	if (gmtdefs.verbose) {
		if (make_path)
			fprintf (stderr, "%s: %ld segments written\n", GMT_program, n_points);
		else
			fprintf (stderr, "%s: %ld points projected\n", GMT_program, n_points);
	}

	if (gmtdefs.verbose && n_skipped) fprintf (stderr, "%s: %ld points skipped because age < 0\n", GMT_program, n_skipped);

	/* Clean up and exit */

	if (!single_rotation) GMT_free ((void *)p);
	if (n_m) {
		GMT_free ((void *)dlon);
		GMT_free ((void *)dlat);
		GMT_free ((void *)dt);
	}
	if (out) GMT_free ((void *)out);
	
	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}
