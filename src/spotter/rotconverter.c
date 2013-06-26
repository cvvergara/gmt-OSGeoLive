/*--------------------------------------------------------------------
 *	$Id: rotconverter.c 9923 2012-12-18 20:45:53Z pwessel $
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
 * Program for converting between finite and stage poles or to add rotations.
 *
 * Author:	Paul Wessel, SOEST, Univ. of Hawaii, Honolulu, HI, USA
 * Date:	24-OCT-2001
 * Version:	1.0
 *		31-MAR-2006: Changed -H to -C to avoid clash with GMT -H
 *
 *-------------------------------------------------------------------------
 * An ASCII stage pole (Euler) file must have following format:
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
 * AN ASCII finite pole must have the following format:
*
 * 1. Any number of comment lines starting with # in first column
 * 2. Any number of blank lines (just carriage return, no spaces)
 * 2. Any number of finite pole records which each have the format:
 *    lon(deg)  lat(deg)  tstop(Ma)  ccw-angle(deg)
 * 3. finite rotations must go from youngest to oldest rotation
 *
 * Example: Duncan & Clague [1985] Pacific-Hotspot rotations:
 *
 * # Time in Ma, angles in degrees
 * #longitude	latitude	time(My)	angle(deg)
 * 285.00000	 68.00000	 42.0000	 34.0000
 * 275.66205	 53.05082	 65.0000	 43.5361
 * 276.02501	 48.34232	 74.0000	 50.0405
 * 279.86436	 46.30610	100.0000	 64.7066
 * 265.37800	 55.69932	150.0000	 82.9957
  */

#include "spotter.h"

int main (int argc, char **argv)
{
	struct EULER *p = NULL;		/* Pointer to array of stage poles */
	struct EULER *a = NULL, *b = NULL;		/* Pointer to arrays of stage poles */

	GMT_LONG i, j, k, n_p, n_a = 1, n_b;	/* Misc. counters */
	GMT_LONG last_sign, n_slash;
	GMT_LONG rate_or_rad = 0;		/* Default unit is degrees/My */

	GMT_LONG error = FALSE;		/* Set to TRUE if arguments are inconsistent */
	GMT_LONG geodetic = TRUE;	/* Want 0-360 range */
	GMT_LONG finite_in = TRUE;	/* Default format is finite rotation poles for both input and output */
	GMT_LONG finite_out = TRUE;	/* FALSE gives stage poles */
	GMT_LONG first = TRUE;		/* TRUE for first input file */
	GMT_LONG south = FALSE;		/* TRUE if we want to report poles in southern hemisphere */
	GMT_LONG north = FALSE;		/* TRUE if we want to report poles in northern hemisphere */
	GMT_LONG online_rot = FALSE;	/* TRUE if we gave a rotation on the commandline rather than file name */
	GMT_LONG no_time = FALSE;	/* TRUE if we gave a rotation on the commandline as lon/lat/angle only */
	GMT_LONG header = FALSE;		/* TRUE to write out a header record */
	GMT_LONG transpose = FALSE;	/* TRUE if we want to change the sign of the final rotation */
	GMT_LONG reduce_angle = FALSE;	/* TRUE to scale stage pole angles by reduce_fact */

	double zero = 0.0;		/* Needed to pass to spotter_init */
	double lon, lat;		/* Pole location for online rotations */
	double t0, t1;			/* Start, stop times for online rotations */
	double angle;			/* Rotation angle for online rotations */
	double reduce_fact = 0.5;	/* To get half-angles */

	char *start_text[2] = {"tstart(My)", "astart(deg)"};	/* Misc. column titles for rates or angles */
	char *end_text[2] = {"tend(My)", "aend(deg)"};
	char *time_text[2] = {"ttime(My)", "tangle(deg)"};
	char fmt[BUFSIZ];

	argc = (int)GMT_begin (argc, argv);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-' && argv[i][1]) {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'V':
				case ':':
					error += GMT_parse_common_options (argv[i], NULL, NULL, NULL, NULL);
					break;

				/* Supplemental parameters */

				case 'F':
					if (strlen (argv[i]) != 4) {
						fprintf (stderr, "%s: ERROR: Must specify -F<in><out>\n", GMT_program);
						error++;
						continue;
					}
					switch (argv[i][2]) {	/* Input/output format */
						case 'F':	/* Finite rotations */
						case 'f':
							finite_in = TRUE;
							break;
						case 'S':	/* Stage rotations */
						case 's':
							finite_in = FALSE;
							break;
						default:
							fprintf (stderr, "%s: ERROR: Must specify f|s\n", GMT_program);
							error++;
							break;
					}
					switch (argv[i][3]) {	/* Output format */
						case 'F':	/* Finite rotations */
						case 'f':
							finite_out = TRUE;
							break;
						case 'S':	/* Stage rotations */
						case 's':
							finite_out = FALSE;
							break;
						default:
							fprintf (stderr, "%s: ERROR: Must specify f|s\n", GMT_program);
							error++;
							break;
					}
					break;

				case 'C':	/* Write column header record */
					header = TRUE;
					if (argv[i][2] == 'a') rate_or_rad = 1;
					break;

				case 'D':	/* Convert to finite rotation poles instead */
					geodetic = FALSE;
					break;

				case 'E':	/* Convert to finite rotation poles instead */
					reduce_angle = TRUE;
					if (argv[i][2]) reduce_fact = atof (&argv[i][2]);
					break;

				case 'N':	/* Ensure all poles reported are in northern hemisphere */
					north = TRUE;
					break;

				case 'S':	/* Ensure all poles reported are in southern hemisphere */
					south = TRUE;
					break;

				case 'T':	/* Transpose the final result (i.e., change sign of rotation) */
					transpose = TRUE;
					break;
					
				case '0': case '1': case '2': case '3': case '4': case '5': case '6':
					case '7': case '8': case '9': case '.':
					break;	/* Probably a rotation lon/lat/angle with negative longitude */
				default:					
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "%s %s - Manipulate finite and stage rotations\n\n", GMT_program, SPOTTER_VERSION);
		fprintf (stderr, "usage: %s [+][-] rotA [[+][-] rotB] [[+][-] rotC] ... [-C[a|t]] [-D] [-E[<factor>]] [-F<in><out>]\n", GMT_program);
		fprintf (stderr, "\t[-N] [-S] [-T] [-V] > outfile\n\n");
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\trotA, rotB, etc. are finite or stage rotation pole files (see -F to set which kind).\n");
		fprintf (stderr, "\t   Alternatively, they can be a single rotation in lon/lat[/tstart[/tstop]]/angle format.\n");
		fprintf (stderr, "\t   For finite rotations with no time info, give lon/lat/angle only.\n");
		fprintf (stderr, "\t   The rotations will be added/subtracted in the order given.\n");
		fprintf (stderr, "\tOPTIONS:\n\n");
		fprintf (stderr, "\t-C Write a column header record with names of each column [No header].\n");
		fprintf (stderr, "\t   Append a if opening angles and t if opening rates [Default]\n");
		fprintf (stderr, "\t-D Report longitudes in -180/+180 range [ Default is 0-360]\n");
		fprintf (stderr, "\t-E Reduce opening angles for stage rotations by factor [0.5]\n");
		fprintf (stderr, "\t   Typically used to get half-rates needed for flowlines\n");
		fprintf (stderr, "\t-F Set input and output file type: f for finite and s for stage rotations [Default is -Fff]\n");
		fprintf (stderr, "\t-N Ensure all poles are in northern hemisphere [ Default ensures positive opening angles/rates]\n");
		fprintf (stderr, "\t-S Ensure all poles are in southern hemisphere [ Default ensures positive opening angles/rates]\n");
		fprintf (stderr, "\t-T Transpose the result (i.e., change sign of rotation angle)\n");
		GMT_explain_option ('V');
		exit (EXIT_FAILURE);
	}

	if (south && north) {
		fprintf (stderr, "%s: ERROR: Cannot specify both -N and -S!\n", GMT_program);
		error++;
	}
	if (reduce_angle && finite_out) {
		fprintf (stderr, "%s: ERROR: -E requires stage rotations on output\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	if (finite_out && no_time)
		sprintf (fmt, "%s%s%s%s%s", gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.d_format);
	else if (finite_out)
		sprintf (fmt, "%s%s%s%s%s%s%s", gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.d_format);
	else
		sprintf (fmt, "%s%s%s%s%s%s%s%s%s", gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.d_format);
	last_sign = +1;
	for (i = 1; i < argc; i++) {
		if (!strncmp (argv[i], "-D", (size_t)2)) continue;
		if (!strncmp (argv[i], "-E", (size_t)2)) continue;
		if (!strncmp (argv[i], "-V", (size_t)2)) continue;
		if (!strncmp (argv[i], "-F", (size_t)2)) continue;
		if (!strncmp (argv[i], "-S", (size_t)2)) continue;
		if (!strncmp (argv[i], "-N", (size_t)2)) continue;
		if (!strncmp (argv[i], "-T", (size_t)2)) continue;

		if (argv[i][0] == '-' && argv[i][1] == '\0') {
			last_sign = -1;
			continue;
		}
		if (argv[i][0] == '+' && argv[i][1] == '\0') {
			last_sign = +1;
			continue;
		}

		if (GMT_access (argv[i], R_OK)) {	/* Not a readable file, is it a lon/lat/t0[/t1]/omega specification? */
			for (j = n_slash = 0; argv[i][j]; j++) if (argv[i][j] == '/') n_slash++;
			if ((n_slash + finite_in) < 3 || (n_slash + finite_in) > 4) {	/* No way it can be a online rotation, cry foul */
				fprintf (stderr, "%s: ERROR: Cannot read file %s\n", GMT_program, argv[i]);
				exit (EXIT_FAILURE);
			}
			else {	/* Tru to decode as a single rotation */

				j = sscanf (argv[i], "%lf/%lf/%lf/%lf/%lf", &lon, &lat, &t0, &t1, &angle);
				if (finite_in) angle = t1, t1 = 0.0;			/* Only 4 input values */
				if (n_slash == 2) angle = t0, t0 = 1.0, t1 = 0.0, no_time = TRUE;	/* Quick lon/lat/angle finite rotation, no time */
				if (t0 < t1) {
					fprintf (stderr, "%s: ERROR: Online rotation has t_start (%g) younger than t_stop (%g)\n", GMT_program, t0, t1);
					exit (EXIT_FAILURE);
				}
				if (angle == 0.0) {
					fprintf (stderr, "%s: ERROR: Online rotation has zero opening angle\n", GMT_program);
					exit (EXIT_FAILURE);
				}
				online_rot = TRUE;
			}
		}
		else
			online_rot = FALSE;

		if (first) {	/* First time loading a rotation model */
			if (online_rot) {
				n_a = 1;
				a = (struct EULER *) GMT_memory (VNULL, (size_t)1, sizeof (struct EULER), GMT_program);
				a[0].lon = lon;	a[0].lat = lat;
				a[0].t_start = t0;	a[0].t_stop = t1;
				a[0].duration = t0 - t1;
				a[0].omega = angle / a[0].duration;
				if (!finite_in) spotter_stages_to_finite (a, n_a, TRUE, TRUE);
			}
			else
				n_a = spotter_init (argv[i], &a, FALSE, finite_in, TRUE, &zero, gmtdefs.verbose);	/* Return finite rotations */
			zero = 0.0;
			if (last_sign == -1) {	/* Leading - sign, simply reverse the rotation angles */
				for (j = 0; j < n_a; j++) {
					a[j].omega = -a[j].omega;
					spotter_cov_of_inverse (&a[j], a[j].C);
				}
				last_sign = 1;
			}
			first = FALSE;
		}
		else {			/* For additional times, load a second model and add/subtract them */
			if (online_rot) {
				n_b = 1;
				b = (struct EULER *) GMT_memory (VNULL, (size_t)1, sizeof (struct EULER), GMT_program);
				b[0].lon = lon;	b[0].lat = lat;
				b[0].t_start = t0;	b[0].t_stop = t1;
				b[0].duration = t0 - t1;
				b[0].omega = angle / b[0].duration;
				if (!finite_in) spotter_stages_to_finite (b, n_b, TRUE, TRUE);
			}
			else
				n_b = spotter_init (argv[i], &b, FALSE, finite_in, TRUE, &zero, gmtdefs.verbose);	/* Return finite rotations */
			zero = 0.0;
			spotter_add_rotations (a, n_a, b, last_sign * n_b, &p, &n_p);			/* Add the two finite rotations sets, returns finite rotations in p */
			GMT_free ((void *)a);
			GMT_free ((void *)b);
			a = p;
			n_a = n_p;
		}
	}
	if (finite_out && no_time && header)
		printf ("#longitude%slatitude%sangle(deg)\n", gmtdefs.field_delimiter, gmtdefs.field_delimiter);
	else if (finite_out && header)	/* Easy, simply output what we've got following a header*/
		printf ("#longitude%slatitude%s%s%sangle(deg)\n", gmtdefs.field_delimiter, gmtdefs.field_delimiter, time_text[rate_or_rad], gmtdefs.field_delimiter);
	else if (finite_out)		/* Easy, simply output what we've got without a header */
		i = 0;	/* Do nothing here really */
	else {	/* Convert finite to stages before output */
		spotter_finite_to_stages (a, n_a, TRUE, TRUE);				/* To ensure we have the right kind of poles for output */
		if (header) printf ("#longitude%slatitude%s%s%s%s%sangle(deg)\n", gmtdefs.field_delimiter, gmtdefs.field_delimiter, start_text[rate_or_rad], gmtdefs.field_delimiter, end_text[rate_or_rad], gmtdefs.field_delimiter);
	}

	for (i = 0; i < n_a; i++) {
		if (transpose) a[i].omega = -a[i].omega;
		if ((south && a[i].lat > 0.0) || (north && a[i].lat < 0.0) || (!(south || north) && a[i].omega < 0.0))	/* flip to antipole */
			a[i].lat = -a[i].lat, a[i].lon += 180.0, a[i].omega = -a[i].omega;
		while (a[i].lon >= 360.0) a[i].lon -= 360.0;				/* Force geodetic longitude range */
		while (a[i].lon <= -180.0) a[i].lon += 360.0;				/* Force geodetic longitude range */
		while (geodetic && a[i].lon < 0.0) a[i].lon += 360.0;			/* Force geodetic longitude range */
		while (!geodetic && a[i].lon > 180.0) a[i].lon -= 360.0;		/* Force a geographic longitude range */
		if (reduce_angle) a[i].omega *= reduce_fact;
		if (finite_out && no_time)
			printf (fmt, a[i].lon, a[i].lat, a[i].omega * a[i].duration);
		else if (finite_out)
			printf (fmt, a[i].lon, a[i].lat, a[i].t_start, a[i].omega * a[i].duration);
		else
			printf (fmt, a[i].lon, a[i].lat, a[i].t_start, a[i].t_stop, a[i].omega * a[i].duration);
		if (a[i].has_cov) {
			double K[9];
			spotter_covar_to_record (&a[i], K);
			for (k = 0; k < 9; k++) {
				printf ("%s", gmtdefs.field_delimiter);
				printf (gmtdefs.d_format, K[k]);
			}
		}
		printf ("\n");
	}

	GMT_free ((void *)a);
	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}
