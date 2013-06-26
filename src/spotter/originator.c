/*--------------------------------------------------------------------
 *	$Id: originator.c 9923 2012-12-18 20:45:53Z pwessel $
 *
 *   Copyright (c) 2000-2013 by P. Wessel
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
 * originator reads file of seamount locations and tries to match each
 * seamount with a probable hotspot by drawing flowines back in time and
 * keeping track of which hotspot is closest to each flowline.  It then
 * reports the closest hotspot, the stage of the flowline involved, the
 * implied pseudo-age of the seamount, and the minimum distance between
 * the flowline and hotspot (in km).
 *
 * Author:	Paul Wessel, SOEST, Univ. of Hawaii, Honolulu, HI, USA
 * Date:	29-DEC-1999
 * Version:	1.0
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
 * ASCII seamount location file(s) must have the following format:
 *
 * 1. Any number of comment lines starting with # in first column
 * 2. Any number of blank lines (just carriage return, no spaces)
 * 3. For special header records, see -H
 * 4. Any number of data records which each have the format:
 *    lon lat height radius crustal_age    (or lat lon ..., see -: option).
 *    crustal_age in Ma, height and radius are not used by originator but
 *    are used by hotspotter.
 *
 * Binary files cannot have header records, and data fields must all be
 * either single or double precision (see -bi option).  Output file will
 * be ASCII since it contains a text string (hotspot ID).
 *
 * The file with a list of hotspots must have the following format:
 *
 * 1. Any number of comment lines starting with # in first column
 * 2. Any number of blank lines (just carriage return, no spaces)
 * 2. Any number of hotspot records which each have the format:
 *    lon(deg)  lat(deg)  id  name
 *    the id is a 3-character tag (e.g., HWI), the name is the
 *    full name of the hotspot (e.g., Hawaii).
 *
 * Example: 
 *
 * # Partial list (Pacific) of HotSpots from Table 1 of Yamaji, 1992
 * #Lon		Lat	Abbreviation	Hotspot_name
 * 167		3	CRL	Caroline
 * 230		46	COB	Cobb
 * 205		20	HWI	Hawaii
 * 221.9	-50.9	LSV	Louisville
 * 220		-29	MDN	MacDonald
 * 221		-11	MRQ	Marquesas
 * 231		-27	PTC	Pitcairn
 * 254		-27	SLG	Sala y Gomez
 * 192		-15	SAM	Samoa
 * 212		-18	SOC	Society
 */
 
#include "spotter.h"

#define KM_PR_RAD (R2D * project_info.DIST_KM_PR_DEG)

struct HOTSPOT_ORIGINATOR {
	struct HOTSPOT *h;	/* Pointer to regular HOTSPOT structure */
	/* Extra variables needed for this program */
	double np_dist;		/* Distance to nearest point on the current flowline */
	double np_sign;		/* "Sign" of this distance (see code) */
	double np_time;		/* Predicted time at nearest point */
	double np_lon;		/* Longitude of nearest point on the current flowline */
	double np_lat;		/* Latitude  of nearest point on the current flowline */
	GMT_LONG nearest;	/* Point id of current flowline node points closest to hotspot */
	GMT_LONG stage;		/* Stage to which seamount belongs */
};

int main (int argc, char **argv)
{
	GMT_LONG x, y, n_max_spots, fno, n_files = 0, n_args, n_input = 5;
	GMT_LONG n_fields, n_expected_fields, n_best_hs = 1, tdz_output = 0, n_out = 0;
	GMT_LONG  i, j, k, n, kk, ns, nh, nc, np, n_read, n_skipped = 0;

	GMT_LONG	error = FALSE, done, nofile = TRUE, truncate_ages = FALSE, first = TRUE, better, use_ID = FALSE;
	GMT_LONG finite = FALSE, degree = FALSE;		/* fitite is TRUE if stage pole file contains finite rotation poles instead */

	double x_smt, y_smt, z_smt, r_smt, t_smt, *c = NULL, *in = NULL, d_km = 5.0, upper_age = 180.0, dist, out[5];
	double max_dist = 1.0e100, hx_dist, hx_dist_km, dist_NA, dist_NX, del_dist, dt = 0.0, A[3], H[3], N[3], X[3];
	double dlon, fix_r, fix_t, yg, yc;

	char *hsfile = CNULL, *efile = CNULL, age[GMT_TEXT_LEN], buffer[BUFSIZ], *not_used = NULL;
	char fmt1[BUFSIZ], fmt2[BUFSIZ];

	FILE *fp = NULL;

	struct EULER *p = NULL;
	struct HOTSPOT *orig_hotspot = NULL;
	struct HOTSPOT_ORIGINATOR *hotspot = NULL, *hot = NULL;

	int comp_hs (const void *p1, const void *p2);

	argc = (int)GMT_begin (argc, argv);

	GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
	GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;

	/* Check command line arguments */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'H':
				case 'V':
				case ':':
				case '\0':
					error += GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'b':
					error += GMT_parse_b_option (&argv[i][2]);
					break;

				case 'C':	/* Use finite rotation poles */
					finite = TRUE;
					break;

				case 'D':
					d_km = atof (&argv[i][2]);
					break;
				case 'E':
					efile = &argv[i][2];
					break;
				case 'F':
					hsfile = &argv[i][2];
					break;
				case 'L':
					n_out = 3;
					switch (argv[i][2]) {
						case 'L':
							degree = TRUE;
						case 'l':
							tdz_output = 3;
							n_out = 5;
							break;
						case 'w':
						case 'W':
							tdz_output = 2;
							break;
						case 'T':
							degree = TRUE;
						case 't':
							tdz_output = 1;
							break;
						default:
							tdz_output = 1;
							break;
					}
					break;
				case 'N':
					upper_age = atof (&argv[i][2]);
					break;
				case 'Q':
					sscanf (&argv[i][2], "%lg/%lg", &fix_r, &fix_t);
					n_input = 3;
					break;
				case 'S':
					n_best_hs = atoi (&argv[i][2]);
					break;
				case 'T':
					truncate_ages = TRUE;
					break;
				case 'W':
					max_dist = atof (&argv[i][2]);
					break;
				case 'Z':
					use_ID = TRUE;
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

	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = n_input;
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < n_input) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least %ld columns\n", GMT_program, n_input);
		error++;
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "%s %s - Associate seamounts with hotspot point sources\n\n", GMT_program, GMT_VERSION);
		fprintf (stderr, "usage: %s [<xyfiles>] -E<euler_file> -F<hotspot_file> [-C] [-D<d_km>]\n", GMT_program);
		fprintf (stderr, "	[-H] [-L[flag]] [-N<upper_age>] [-Qr/t] [-S<n_hs>] [-T] [-V] [-W<maxdist>] [-Z] [%s]\n\n", GMT_t_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\txyfiles is one or more seamount (x,y,z,r,t) files\n");
		fprintf (stderr, "\t-E specifies the rotations to be used (see man page for format)\n\n");
		fprintf (stderr, "\t-F Specify file name for hotspot locations.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-C The file given with -E contains finite rotation poles [Default is stage poles]\n");
		fprintf (stderr, "\t-D set sampling interval in km along tracks [5].\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-L Output information for closest approach for nearest hotspot only (ignores -S).\n");
		fprintf (stderr, "\t   -Lt gives (time, dist, z) [Default].\n");
		fprintf (stderr, "\t   -Lw gives (omega, dist, z).\n");
		fprintf (stderr, "\t   -Ll gives (lon, lat, time, dist, z).\n");
		fprintf (stderr, "\t   dist is in km; use upper case T,W,L to get dist in spherical degrees.\n");
		fprintf (stderr, "\t-N set age (in m.y.) for seafloor where age == NaN [180].\n");
		fprintf (stderr, "\t-Q input files has (x,y,z) only. Append constant r/t to use.\n");
		fprintf (stderr, "\t-S Report the <n_hs> closest hotSpots [1].\n");
		fprintf (stderr, "\t-T Truncate seamount ages exceeding the upper age set with -N [no truncation] \n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t-W Only report seamounts whose closest encounter to a hotspot is less than <maxdist> km\n");
		fprintf (stderr, "\t   [Default reports for all seamounts] \n");
		fprintf (stderr, "\t-Z Write hotspot ID number rather than hotspot TAG\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "	   Default is 5 input columns\n");
		exit (EXIT_FAILURE);
	}

	if (!hsfile) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -F:  Must specify hotspot file\n", GMT_program);
		error = TRUE;
	}
	if (!efile) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E:  Must specify Euler pole file\n", GMT_program);
		error = TRUE;
	}
	if (d_km <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -D:  Must specify a positive interval\n", GMT_program);
		error = TRUE;
	}
	if (max_dist <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -W:  Must specify a positive distance in km\n", GMT_program);
		error = TRUE;
	}

	if (error) exit (EXIT_FAILURE);

	GMT_lat_swap_init ();	/* Initialize auxiliary latitude machinery */

	nh = spotter_hotspot_init (hsfile, &orig_hotspot);
	if (n_best_hs <= 0 || n_best_hs > nh) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S:  Give value between 1 and %ld\n", GMT_program, nh);
		exit (EXIT_FAILURE);
	}
	n_max_spots = MIN (n_best_hs, nh);

	hotspot = (struct HOTSPOT_ORIGINATOR *) GMT_memory (VNULL, (size_t)nh, sizeof (struct HOTSPOT_ORIGINATOR), GMT_program);
	for (i = 0; i < nh; i++) {
		hotspot[i].h = &orig_hotspot[i];	/* Point to the original hotspot structures */
		hotspot[i].np_dist = 1.0e100;
	}
	
	ns = spotter_init (efile, &p, TRUE, finite, FALSE, &upper_age, gmtdefs.verbose);

	hot = (struct HOTSPOT_ORIGINATOR *) GMT_memory (VNULL, (size_t)nh, sizeof (struct HOTSPOT_ORIGINATOR), GMT_program);

	x = (gmtdefs.xy_toggle[1]) ? 1 : 0;	y = 1 - x;		/* Set up which columns have x and y on output */

	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;
	n_args = (argc > 1) ? argc : 2;

	sprintf (fmt1, "%s%s%s%s%s%s%s%s%%s%s", gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.field_delimiter);
	if (use_ID)
		sprintf (fmt2, "%s%%d%s%%ld%s%s%s%s", gmtdefs.field_delimiter, gmtdefs.field_delimiter, gmtdefs.field_delimiter, gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.d_format);
	else
		sprintf (fmt2, "%s%%s%s%%ld%s%s%s%s", gmtdefs.field_delimiter, gmtdefs.field_delimiter, gmtdefs.field_delimiter, gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.d_format);
	
	done = FALSE;
	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : n_input;
	n = 0;
	if (!tdz_output) n_out = n_expected_fields;
	for (fno = 1; !done && fno < n_args; fno++) {	/* Loop over all input files */

		if (!nofile && argv[fno][0] == '-') continue;
		if (nofile) {
			fp = GMT_stdin;
			done = TRUE;
		}
		else if ((fp = GMT_fopen (argv[fno], GMT_io.r_mode)) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[fno]);
			continue;
		}

		if (!nofile && gmtdefs.verbose) fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[fno]);

		n_read = 0;
		if (GMT_io.io_header[GMT_IN]) {
			for (i = 0; i < GMT_io.n_header_recs; i++) {
				not_used = GMT_fgets (GMT_io.segment_header, BUFSIZ, fp);
				if (first && !GMT_io.binary[GMT_OUT] && GMT_io.io_header[GMT_OUT]) GMT_fputs (GMT_io.segment_header, GMT_stdout);
				n_read++;
			}
			first = FALSE;
		}

		while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */
			n_read++;
			while ((GMT_io.status & GMT_IO_SEGMENT_HEADER) && !(GMT_io.status & GMT_IO_EOF)) {
				GMT_write_segmentheader (GMT_stdout, n_out);
				n_fields = GMT_input (fp, &n_expected_fields, &in);
				n_read++;
			}
			if (GMT_io.status & GMT_IO_EOF) continue;

			if (GMT_io.status & GMT_IO_MISMATCH) {
				fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
				continue;
			}
			
			if (n_input == 3) {	/* set constant r,t values */
				in[3] = fix_r;
				in[4] = fix_t;
			}
			if (GMT_is_dnan (in[4])) {	/* Age is NaN, assign value */
				t_smt = upper_age;
			}
			else {			/* Assign given value, truncate if necessary */
				t_smt = in[4];
				if (t_smt > upper_age) {
					if (truncate_ages) {
						t_smt = upper_age;
					}
					else {
						fprintf (stderr, "%s: Seamounts near line %ld has age (%g) > oldest stage (%g) (skipped)\n", GMT_program, n_read, t_smt, upper_age);
						continue;
					}
				}
			}
			if (t_smt < 0.0) {	/* Negative ages are flags for points to be skipped */
				n_skipped++;
				continue;
			}
						
			x_smt = in[0] * D2R;
			y_smt = GMT_lat_swap (in[GMT_Y], GMT_LATSWAP_G2O) * D2R;	/* Convert to geocentric */
			z_smt = in[2];
			r_smt = in[3];

			if (gmtdefs.verbose && !(n % 10)) fprintf (stderr, "%s: Working on seamount # %5ld\r", GMT_program, n);

			nc = spotter_forthtrack (&x_smt, &y_smt, &t_smt, (GMT_LONG)1, p, ns, d_km, 0.0, TRUE, NULL, &c);

			np = (GMT_LONG) c[0];

			memcpy ((void *)hot, (void *)hotspot, nh * sizeof (struct HOTSPOT_ORIGINATOR));

			for (kk = 0, k = 1; kk < np; kk++, k += 3) {	/* For this seamounts track */
				yg = GMT_lat_swap (R2D*c[k+1], GMT_LATSWAP_O2G);	/* Convert back to geodetic */
				for (j = 0; j < nh; j++) {	/* For all hotspots */
					dist = GMT_great_circle_dist (hot[j].h->lon, hot[j].h->lat, R2D * c[k], yg);
					if (!degree) dist *= project_info.DIST_KM_PR_DEG;
					if (dist < hot[j].np_dist) {
						hot[j].np_dist = dist;
						hot[j].nearest = kk;	/* Index of nearest point on the flowline */
					}
				}
			}
			for (j = 0; j < nh; j++) {

				yc = GMT_lat_swap (hot[j].h->lat, GMT_LATSWAP_G2O);	/* Convert to geocentric */
				GMT_geo_to_cart (yc, hot[j].h->lon, H, TRUE);	/* 3-D Cartesian vector of this hotspot */

				/* Fine-tune the nearest point by considering intermediate points along greatcircle between knot points */

				k = 3 * hot[j].nearest + 1;			/* Corresponding index for x into the (x,y,t) array c */
				GMT_geo_to_cart (c[k+1], c[k], N, FALSE);	/* 3-D vector of nearest node to this hotspot */
				better = FALSE;
				if (hot[j].nearest > 0) {	/* There is a point along the flowline before the nearest node */
					GMT_geo_to_cart (c[k-2], c[k-3], A, FALSE);	/* 3-D vector of end of this segment */
					if (GMT_great_circle_intersection (A, N, H, X, &hx_dist) == 0) {	/* X is between A and N */
						hx_dist_km = d_acos (hx_dist) * KM_PR_RAD;
						if (hx_dist_km < hot[j].np_dist) {	/* This intermediate point is even closer */
							GMT_cart_to_geo (&hot[j].np_lat, &hot[j].np_lon, X, TRUE);
							hot[j].np_dist = hx_dist_km;
							dist_NA = d_acos (fabs (GMT_dot3v (A, N))) * KM_PR_RAD;
							dist_NX = d_acos (fabs (GMT_dot3v (X, N))) * KM_PR_RAD;
							del_dist = dist_NA - dist_NX;
							dt = (del_dist > 0.0) ? (c[k+2] - c[k-1]) * dist_NX / del_dist : 0.0;
							better = TRUE;
						}
					}
				}
				if (hot[j].nearest < (np-1) ) {	/* There is a point along the flowline after the nearest node */
					GMT_geo_to_cart (c[k+4], c[k+3], A, FALSE);	/* 3-D vector of end of this segment */
					if (GMT_great_circle_intersection (A, N, H, X, &hx_dist) == 0) {	/* X is between A and N */
						hx_dist_km = d_acos (hx_dist) * KM_PR_RAD;
						if (hx_dist_km < hot[j].np_dist) {	/* This intermediate point is even closer */
							GMT_cart_to_geo (&hot[j].np_lat, &hot[j].np_lon, X, TRUE);
							hot[j].np_dist = hx_dist_km;
							dist_NA = d_acos (fabs (GMT_dot3v (A, N))) * KM_PR_RAD;
							dist_NX = d_acos (fabs (GMT_dot3v (X, N))) * KM_PR_RAD;
							del_dist = dist_NA - dist_NX;
							dt = (del_dist > 0.0) ? (c[k+5] - c[k+2]) * dist_NX / del_dist : 0.0;
							better = TRUE;
						}
					}
				}
				if (better) {	/* Point closer to hotspot was found between nodes */
					hot[j].np_time = c[k+2] + dt;	/* Add time adjustment */
				}
				else {	/* Just use node coordinates */
					hot[j].np_lon  = c[k] * R2D;	/* Longitude of the flowline's closest approach to hotspot */
					hot[j].np_lat  = GMT_lat_swap (c[k+1] * R2D, GMT_LATSWAP_O2G);	/* Geodetic latitude  of the flowline's closest approach to hotspot */
					hot[j].np_time = c[k+2];	/* Predicted time at the flowline's closest approach to hotspot */
				}

				/* Assign sign to distance: If the vector from the hotspot pointing up along the trail is positive
				 * x-axis and y-axis is normal to that, flowlines whose closest approach point's longitude is
				 * further east are said to have negative distance. */
				 
				dlon = fmod (hot[j].h->lon - hot[j].np_lon, 360.0);
				if (fabs (dlon) > 180.0) dlon = copysign (360.0 - fabs (dlon), -dlon);
				hot[j].np_sign = copysign (1.0, dlon);
				 
				/* Assign stage id for this point on the flowline */

				k = 0;
				while (k < ns && hot[j].np_time <= p[k].t_stop) k++;
				hot[j].stage = ns - k;
				if (hot[j].stage == 0) hot[j].stage++;
			}

			if (nh > 1) qsort ((void *)hot, (size_t)nh, sizeof(struct HOTSPOT_ORIGINATOR), comp_hs);

			if (hot[0].np_dist < max_dist) {
				if (tdz_output == 1) {	/* Want time, dist, z output */
					out[0] = hot[0].np_time;
					out[1] = hot[0].np_dist * hot[0].np_sign;
					out[2] = z_smt;
					GMT_output (GMT_stdout, n_out, out);
				}
				else if (tdz_output == 2) {	/* Want omega, dist, z output */
					out[0] = spotter_t2w (p, ns, hot[0].np_time);
					out[1] = hot[0].np_dist * hot[0].np_sign;
					out[2] = z_smt;
					GMT_output (GMT_stdout, n_out, out);
				}
				else if (tdz_output == 3) {	/* Want x, y, time, dist, z output */
					out[0] = hot[0].np_lon;
					out[1] = hot[0].np_lat;
					out[2] = hot[0].np_time;
					out[3] = hot[0].np_dist * hot[0].np_sign;
					out[4] = z_smt;
					GMT_output (GMT_stdout, n_out, out);
				}
				else {	/* Conventional originator output */
					if (t_smt == 180.0)
						strcpy (age, "NaN");
					else
						sprintf (age, "%g", t_smt);
					sprintf (buffer, fmt1, in[x], in[y], z_smt, r_smt, age);
					GMT_fputs (buffer, GMT_stdout);
					if (use_ID)
						for (j = 0; j < n_max_spots; j++) {
							sprintf (buffer, fmt2, hot[j].h->id, hot[j].stage, hot[j].np_time, hot[j].np_dist);
							GMT_fputs (buffer, GMT_stdout);
						}
					else
						for (j = 0; j < n_max_spots; j++) {
							sprintf (buffer, fmt2, hot[j].h->abbrev, hot[j].stage, hot[j].np_time, hot[j].np_dist);
							GMT_fputs (buffer, GMT_stdout);
						}
					GMT_fputs ("\n", GMT_stdout);
				}
			}

			GMT_free ((void *)c);
			n++;
		}

		if (fp != GMT_stdin) GMT_fclose (fp);
	}

	if (gmtdefs.verbose) fprintf (stderr, "%s: Working on seamount # %5ld\n", GMT_program, n);

	GMT_free ((void *)hotspot);
	GMT_free ((void *)orig_hotspot);
	GMT_free ((void *)hot);
	GMT_free ((void *)p);

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

int comp_hs (const void *p1, const void *p2)
{
	struct HOTSPOT_ORIGINATOR *a, *b;

	a = (struct HOTSPOT_ORIGINATOR *) p1;
	b = (struct HOTSPOT_ORIGINATOR *) p2;
	if (a->np_dist < b->np_dist) return (-1);
	if (a->np_dist > b->np_dist) return (1);
	return (0);
}
