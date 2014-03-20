/*--------------------------------------------------------------------
*	$Id: mapproject.c 10173 2014-01-01 09:52:34Z pwessel $
*
*	Copyright (c) 1991-2014 by P. Wessel and W. H. F. Smith
*	See LICENSE.TXT file for copying and redistribution conditions.
*
*	This program is free software; you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation; version 2 or any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	Contact info: gmt.soest.hawaii.edu
*--------------------------------------------------------------------*/
/*
* mapproject reads a pair of coordinates [+ optional data fields] from
* standard input or file(s) and transforms the coordinates according to the
* map projection selected. See the man page for projections currently supported.
*
* The default is to expect longitude, latitude, [and optional datavalues],
* and return x, y, [ and optional datavalues], but if the -I option is used,
* the reverse is true.  Specifying -C means that the origin of projected coordinates
* should be set to origin of projection  [Default origin is lower left corner of "map"].
* If your data is lat lon instead of lon lat [Default], use -: to toggle x/y -> y/x.
* Note that only unprojected values are affected by the -: switch.  True x,y values are
* always printed out as x,y.  Option -G allows calculation of distances along track or
* to a fixed point, while -L calculates shortest distances to a line.
* Finally, datum conversions can also be done, alone or in series with a
* map projection.
*
*
* Author:	Paul Wessel
* Date:	1-MAR-1990
* Version:	4
*
*/

#include "gmt.h"
#include "gmt_proj.h"   /* For auxiliary latitude stuff */

struct MAPPROJECT_CTRL {	/* All control options for this program (except common args) */
	/* active is TRUE if the option has been activated */
	struct A {	/* -Ab|B|f|F<lon0>/<lat0> */
		GMT_LONG active;
		GMT_LONG azims;
		GMT_LONG reverse;	/* TRUE if we want back-azimuths instead of regular azimuths */
		GMT_LONG geodesic;	/* TRUE if we want geodesic azimuths [Default is great circle azimuths] */
		double lon, lat;	/* Fixed point of reference */
	} A;
	struct C {	/* -C[<false_easting>/<false_northing>] */
		GMT_LONG active;
		double easting, northing;	/* Shifts */
	} C;
	struct D {	/* -D<c|i|m|p> */
		GMT_LONG active;
		char unit;
	} D;
	struct E {	/* -E[<datum>] */
		GMT_LONG active;
		struct GMT_DATUM datum;	/* Contains a, f, xyz[3] */
	} E;
	struct F {	/* -F[k|m|n|i|c|p] */
		GMT_LONG active;
		char unit;
	} F;
	struct G {	/* -G<lon0>/<lat0>[/e|E|k|K|m|M|n|N|c|C|d|D] */
		GMT_LONG active;
		GMT_LONG mode;		/* 1 = distance to fixed point, 2 = cumulative distances, 3 = incremental distances, 4 = 2nd point in cols 3/4 */
		double lon, lat;	/* Fixed point of reference */
		char unit;
	} G;
	struct I {	/* -I */
		GMT_LONG active;
	} I;
	struct L {	/* -L<line.xy>[/<unit>] */
		GMT_LONG active;
		GMT_LONG mode;	/* 0 = dist to nearest point, 1 = also get the point, 2 = instead get seg#, pt# */
		char *file;	/* Name of file with lines */
		char unit;
	} L;
	struct Q {	/* -Q[e|d] */
		GMT_LONG active;
		GMT_LONG mode;	/* 1 = print =Qe, 2 print -Qd, 3 print both */
	} Q;
	struct S {	/* -S */
		GMT_LONG active;
	} S;
	struct T {	/* -T[h]<from>[/<to>] */
		GMT_LONG active;
		GMT_LONG heights;	/* True if we have heights */
		struct GMT_DATUM from;	/* Contains a, f, xyz[3] */
		struct GMT_DATUM to;	/* Contains a, f, xyz[3] */
	} T;
};

int main (int argc, char **argv)
{
	GMT_LONG k, fno, n_files = 0, x, y, n_args, unit = 0, n_slash;
	GMT_LONG n_fields, n_expected_fields, proj_type = 0, *n_output = NULL, save[2] = {0,0}, two, n_out;
	GMT_LONG fmt[2], pos, i, j, n = 0, n_read = 0, n_read_in_seg;

	GMT_LONG error = FALSE, do_geo_conv = FALSE;
	GMT_LONG nofile = TRUE, done = FALSE, first = TRUE;
	GMT_LONG pure_ascii = FALSE, line_start = TRUE, long_verbose = FALSE;
	GMT_LONG geodetic_calc = FALSE, greenwich = FALSE;
	GMT_LONG datum_conv_only = FALSE, double_whammy = FALSE, shift_xy = FALSE;

	double west = 0.0, east = 0.0, south = 0.0, north = 0.0, x_in = 0.0, y_in = 0.0, d = 0.0, s = 0.0;
	double xmin, xmax, ymin, ymax, *in = NULL, *out = NULL, fwd_scale, inv_scale, xtmp, ytmp;
	double x_in_min, x_in_max, y_in_min, y_in_max, inch_to_unit, unit_to_inch;
	double x_out_min, x_out_max, y_out_min, y_out_max, u_scale, d_scale;
	double xnear, ynear, lon_prev = 0, lat_prev = 0;

	char line[BUFSIZ], format[BUFSIZ], unit_name[GMT_TEXT_LEN], scale_unit_name[GMT_TEXT_LEN];
	char txt_a[GMT_LONG_TEXT], txt_b[GMT_LONG_TEXT], p[BUFSIZ], c;
	char from[GMT_LONG_TEXT], to[GMT_LONG_TEXT];

	FILE *fp = NULL;

	struct GMT_TABLE *xyline = NULL;
	struct MAPPROJECT_CTRL *Ctrl = NULL;

	PFD distance_func;
	PFI near_a_line;
	PFD azimuth_func;

	void *New_mapproject_Ctrl (), Free_mapproject_Ctrl (struct MAPPROJECT_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct MAPPROJECT_CTRL *) New_mapproject_Ctrl ();		/* Allocate and initialize defaults in a new control structure */
	
	out = (double *)NULL;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'V':
					if (argv[i][2] == 'l') long_verbose = TRUE;
				case 'R':
				case 'H':
				case 'J':
				case 'M':
				case ':':
				case 'b':
				case 'f':
				case 'g':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */

				case 'A':
					Ctrl->A.active = TRUE;
					n = sscanf (&argv[i][2], "%c%[^/]/%s", &c, txt_a, txt_b);
					if (n < 1) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR.  Expected -Ab|B|f|F[<lon0>/<lat0]>\n", GMT_program);
						error++;
					}
					else {
						switch (c) {
							case 'B':
								Ctrl->A.geodesic = TRUE;
							case 'b':
								Ctrl->A.reverse = TRUE;
								break;
							case 'F':
								Ctrl->A.geodesic = TRUE;
							case 'f':
								break;
							default:
								fprintf (stderr, "%s: GMT SYNTAX ERROR.  Expected -Ab|B|f|F[<lon0>/<lat0>]\n", GMT_program);
								error++;
								break;
						}
						if (n == 3) {
							error += GMT_verify_expectations (GMT_io.in_col_type[0], GMT_scanf_arg (txt_a, GMT_io.in_col_type[0], &Ctrl->A.lon), txt_a);
							error += GMT_verify_expectations (GMT_io.in_col_type[1], GMT_scanf_arg (txt_b, GMT_io.in_col_type[1], &Ctrl->A.lat), txt_b);
						}
						else
							Ctrl->A.azims = TRUE;
					}
					break;
				case 'C':
					Ctrl->C.active = TRUE;
					if (argv[i][2]) {	/* Also gave shifts */
						n = sscanf (&argv[i][2], "%lf/%lf", &Ctrl->C.easting, &Ctrl->C.northing);
						if (n != 2) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR.  Expected -C[<false_easting>/<false_northing>]\n", GMT_program);
							error++;
						}
						shift_xy = TRUE;
					}
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					Ctrl->D.unit = argv[i][2];
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					if (GMT_set_datum (&argv[i][2], &Ctrl->E.datum) == -1) error++;
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					Ctrl->F.unit = argv[i][2];
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					for (n_slash = 0, k = 2; argv[i][k]; k++) if (argv[i][k] == '/') n_slash++;
					if (n_slash == 2 || n_slash == 1) {	/* Got -Glon0/lat0[/units] */
						Ctrl->G.mode = 1;
						n = sscanf (&argv[i][2], "%[^/]/%[^/]/%c", txt_a, txt_b, &c);
						if (n < 2) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR.  Expected -G<lon0>/<lat0>[/e|E|k|K|m|M|n|N|c|C|d|D]\n", GMT_program);
							error++;
						}
						if (n_slash == 2) Ctrl->G.unit = c;
						if (Ctrl->G.unit == 'c') GMT_io.in_col_type[0] = GMT_io.in_col_type[1] = GMT_IS_FLOAT;	/* Cartesian data only */
						error += GMT_verify_expectations (GMT_io.in_col_type[0], GMT_scanf_arg (txt_a, GMT_io.in_col_type[0], &Ctrl->G.lon), txt_a);
						error += GMT_verify_expectations (GMT_io.in_col_type[1], GMT_scanf_arg (txt_b, GMT_io.in_col_type[1], &Ctrl->G.lat), txt_b);
					}
					else if (argv[i][2] == '+') {				/* Got -G+[units] */
						Ctrl->G.mode = 4;
						Ctrl->G.unit = argv[i][3];
					}
					else if (argv[i][2] == '-') {				/* Got -G-[units] */
						Ctrl->G.mode = 3;
						Ctrl->G.unit = argv[i][3];
					}
					else {				/* Got -G[units] */
						Ctrl->G.mode = 2;
						Ctrl->G.unit = argv[i][2];
					}
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					break;
				case 'L':
					Ctrl->L.active = TRUE;
					Ctrl->L.file = strdup (&argv[i][2]);
					k = (int)strlen (Ctrl->L.file) - 1;
					if (Ctrl->L.file[k] == '+') {	/* Flag to get point number instead of coordinates at nearest point on line */
						Ctrl->L.mode = 3;
						Ctrl->L.file[k] = '\0';
						k--;
					}
					k--;
					if (k >= 0 && Ctrl->L.file[k] == '/' && strchr ("ekmndcC", Ctrl->L.file[k+1])) {
						Ctrl->L.unit = Ctrl->L.file[k+1];
						Ctrl->L.file[k] = 0;
					}
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
					if (argv[i][2] == 'e') Ctrl->Q.mode |= 1;
					if (argv[i][2] == 'd') Ctrl->Q.mode |= 2;
					if (argv[i][2] == '\0') Ctrl->Q.mode = 3;
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					k = 2;
					if (argv[i][k] == 'h') {	/* We will process lon, lat, height data */
						k = 3;
						Ctrl->T.heights = TRUE;	/* If FALSE we set height = 0 */
					}

					if (strchr (&argv[i][k], '/')) {	/* Gave from/to */
						sscanf (&argv[i][k], "%[^/]/%s", from, to);
					}
					else {	/* to not given, set to - which means WGS-84 */
						strcpy (to, "-");
						strcpy (from, &argv[i][k]);
					}
					if (GMT_set_datum (to, &Ctrl->T.to) == -1 || GMT_set_datum (from, &Ctrl->T.from) == -1) {
						error++;
						fprintf (stderr, "%s: GMT SYNTAX ERROR -T: Usage -T[h]<from>[/<to>]\n", GMT_program);
					}
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
		fprintf (stderr, "mapproject %s - Forward and Inverse map transformations and geodesy\n\n", GMT_VERSION);
		fprintf (stderr, "usage: mapproject <infiles> %s %s [-C[<dx/dy>]]\n", GMT_J_OPT, GMT_Rgeo_OPT);
		fprintf (stderr, "\t[-Ab|B|f|F[<lon0/lat0>]] [-Dc|i|m|p] [-E[<datum>]] [-F[k|m|n|i|c|p]] [-G[<lon0/lat0>/][+|-][<unit>]\n");
		fprintf (stderr, "\t[%s] [-I] [-L<line.xy>[/<unit>]][+] [-Q[e|d]] [-S] [-T[h]<from>[/<to>]\n", GMT_H_OPT);
		fprintf (stderr, "\t[-V[l]] [%s] [%s] [%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_b_OPT, GMT_f_OPT, GMT_g_OPT, GMT_m_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\tinfiles (in ASCII or binary) has 2 or more columns.  If no file(s) is given, standard input is read.\n");
		GMT_explain_option ('J');
		GMT_explain_option ('R');
		fprintf (stderr, "\t   If UTM and -C are used then -R is optional (automatically set to match UTM zone).\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A Calculate azimuths from previous point in the input data with -Af. If a specified\n");
		fprintf (stderr, "\t   point is provided, all azimuths are computed with respect to that point.\n");
		fprintf (stderr, "\t   Use -Ab to calculate backazimuths from data to previous or the specified point.\n");
		fprintf (stderr, "\t   Upper case B or F gives azimuths of geodesics using current ellipsoid.\n");
		fprintf (stderr, "\t-C returns x/y relative to projection center [Default is relative to lower left corner].\n");
		fprintf (stderr, "\t   Optionally append dx/dy to add (or subtract if -I) (i.e., false easting & northing) [0/0].\n");
		fprintf (stderr, "\t   Units are plot units unless -F[unit] is set in which case the unit is always meters.\n");
		fprintf (stderr, "\t-D Temporarily reset MEASURE_UNIT to be c (cm), i (inch), m (meter), or p (point).\n");
		fprintf (stderr, "\t   Cannot be used if -F is set.\n");
		fprintf (stderr, "\t-E Convert (lon, lat, h) to Earth Centered Earth Fixed (ECEF) coordinates [-I for inverse].\n");
		fprintf (stderr, "\t   Specify <datum> using datum ID (see -Qd or man page) or as <ellipsoid>:<dx,dy,dz>\n");
		fprintf (stderr, "\t   where <ellipsoid> may be ellipsoid ID (see -Qe or man page) or <semimajor>[,<inv_flattening>].\n");
		fprintf (stderr, "\t   If <datum> = - or not given we assume WGS-84.\n");
		fprintf (stderr, "\t-F force projected values to be in actual meters [Default uses the given plot scale].\n");
		fprintf (stderr, "\t   Specify unit by appending k (km), m (miles), n (nautical miles), i (inch), c (cm), or p (points).\n");
		fprintf (stderr, "\t-G Calculate distances to specified point OR cumulative distances along track (if point not given).\n");
		fprintf (stderr, "\t   Use -G+[<unit>] to get provide <lon0> <lat0> from two extra input columns.\n");
		fprintf (stderr, "\t   Use -G-[<unit>] to get distance increments rather than cumulate distances along track.\n");
		fprintf (stderr, "\t   Specify unit as m(e)ter, (k)m, (m)ile, (n)autical mile, (d)egree, or (c)artesian in user units.\n");
		fprintf (stderr, "\t   Unit C means Cartesian distances after first projecting the input coordinates (-R, -J).\n");
		fprintf (stderr, "\t   Units E, K, M, N, D mean geodesic distance using current ellipsoid [lower case is spherical].\n");
		fprintf (stderr, "\t   Default is meters on spherical earth.\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I means Inverse, i.e., get lon/lat from x/y input [Default is lon/lat -> x/y].\n");
		fprintf (stderr, "\t-L Calculate minimum distances to specified line(s) in the file <line.xy>.\n");
		fprintf (stderr, "\t   Specify unit as m(e)ter, (k)m, (m)ile, (n)autical mile, (d)egree, or (c)artesian in user units.\n");
		fprintf (stderr, "\t   Unit C means Cartesian distances after first projecting the input coordinates (-R, -J).\n");
		fprintf (stderr, "\t   Calculations uses spherical approximations.  Default unit is meters.\n");
		fprintf (stderr, "\t   Three columns are added on output: min dist and lon, lat of the closest point on the line.\n");
		fprintf (stderr, "\t   Append + to get line segment id and fractional point number instead of lon/lat.\n");
		fprintf (stderr, "\t-Q list projection parameters and exit.  For subsets [Default is all] use\n");
		fprintf (stderr, "\t   -Qe shows ellipsoid parameters,\n");
		fprintf (stderr, "\t   -Qd shows datum parameters.\n");
		fprintf (stderr, "\t-S means Suppress points outside region.\n");
		fprintf (stderr, "\t-T means coordinate transformation from datum <from> to datum <to>.\n");
		fprintf (stderr, "\t   Prepend h if input data are lon, lat, height [Default sets height = 0].\n");
		fprintf (stderr, "\t   Specify datums using datum ID (see -Qd or man page) or as <ellipsoid>:<dx,dy,dz>\n");
		fprintf (stderr, "\t   where <ellipsoid> may be ellipsoid ID (see -Qe or man page) or <semimajor>[,<inv_flattening>].\n");
		fprintf (stderr, "\t   <from> = - means WGS-84.  If /<to> is not given we assume WGS-84.\n");
		fprintf (stderr, "\t   -T can be used as pre- or post- (-I) processing for -J -R.\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t   Append l for long verbose, reporting every 1000 points.\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2 input columns.\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('g');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (Ctrl->Q.mode & 1) {	/* List ellipsoid parameters */
		fprintf (stderr, "GMT supports %d ellipsoids, given below (-> indicates default setting)\n", GMT_N_ELLIPSOIDS);
		fprintf (stderr, "  ID                      Date        a           1/f\n");
		fprintf (stderr, "-----------------------------------------------------------\n");
		for (i = 0; i < GMT_N_ELLIPSOIDS; i++) {
			(i == gmtdefs.ellipsoid) ? fprintf (stderr, "->") : fprintf (stderr, "  ");
			fprintf (stderr, "%-23s %4ld %13.3f %14.9f\n", gmtdefs.ref_ellipsoid[i].name, gmtdefs.ref_ellipsoid[i].date, gmtdefs.ref_ellipsoid[i].eq_radius, 1.0/gmtdefs.ref_ellipsoid[i].flattening);
		}
		fprintf (stderr, "-----------------------------------------------------------\n");
	}
	if (Ctrl->Q.mode & 2) {	/* List datum parameters */
		fprintf (stderr, "GMT supports %d datums, given below (-> indicates default setting)\n", GMT_N_DATUMS);
		fprintf (stderr, "  ID  Name                               Ellipsoid                 x     y     z   Region\n");
		fprintf (stderr, "-----------------------------------------------------------------------------------------\n");
		for (i = 0; i < GMT_N_DATUMS; i++) {
			(!strcmp (gmtdefs.datum[i].name, "WGS 1984")) ? fprintf (stderr, "->") : fprintf (stderr, "  ");
			fprintf (stderr, "%3ld %-34s %-23s %5.0f %5.0f %5.0f %s\n", i, gmtdefs.datum[i].name, gmtdefs.datum[i].ellipsoid, gmtdefs.datum[i].xyz[0], gmtdefs.datum[i].xyz[1], gmtdefs.datum[i].xyz[2], gmtdefs.datum[i].region);
		}
		fprintf (stderr, "-----------------------------------------------------------------------------------------\n");
	}
	if (Ctrl->Q.mode) exit (EXIT_FAILURE);

	geodetic_calc = (Ctrl->G.mode || Ctrl->A.active || Ctrl->L.active);

	if (Ctrl->T.active && (Ctrl->G.mode + Ctrl->E.active + Ctrl->L.active) > 0) {	/* No good... */
		fprintf (stderr, "%s: GMT SYNTAX ERROR: -T cannot work with -E, -G or -L\n", GMT_program);
		error++;
	}
	if (geodetic_calc && Ctrl->I.active) {	/* No good... */
		fprintf (stderr, "%s: GMT SYNTAX ERROR: -A, -G, and -L cannot work with -I\n", GMT_program);
		error++;
	}
	if (project_info.projection == GMT_NO_PROJ && (Ctrl->G.mode || Ctrl->L.active) && proj_type == 2) {	/* Must have -J */
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -J option with selected form of -G or -L\n", GMT_program);
		error++;
	}
	if (!project_info.region_supplied && project_info.projection == GMT_UTM && Ctrl->C.active) {	/* Set default UTM region from zone info */
		if (!project_info.utm_zoney && project_info.utm_hemisphere == 0) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  -Ju need zone specification with latitude or hemisphere selection\n", GMT_program);
			error++;
		}
		if (GMT_UTMzone_to_wesn (project_info.utm_zonex, project_info.utm_zoney, project_info.utm_hemisphere, &west, &east, &south, &north)) {
			error++;
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Bad UTM zone\n", GMT_program);
		}
		else if (gmtdefs.verbose)
			fprintf (stderr, "%s: UTM zone used to generate region %g/%g/%g/%g\n", GMT_program, west, east, south, north);
		
		project_info.region_supplied = TRUE;
	}
	if (!project_info.region_supplied && !(geodetic_calc || Ctrl->T.active || Ctrl->E.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	if ((Ctrl->D.active + Ctrl->F.active) == 2) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Can specify only one of -D and -F\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 2) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 2 columns\n", GMT_program);
		error++;
	}
	if (((Ctrl->T.active && GMT_datum.h_given) || Ctrl->E.active) && GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 3) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  For -E or -T, binary input data (-bi) must have at least 3 columns\n", GMT_program);
		error++;
	}

	if (Ctrl->D.active) GMT_err_fail (GMT_set_measure_unit (Ctrl->D.unit), "-D");
	if (Ctrl->G.active) error += GMT_get_dist_scale (Ctrl->G.unit, &d_scale, &proj_type, &distance_func);
	if (Ctrl->L.active) error += GMT_get_dist_scale (Ctrl->L.unit, &d_scale, &proj_type, &distance_func);
	if (Ctrl->T.active) GMT_datum_init (&Ctrl->T.from, &Ctrl->T.to, Ctrl->T.heights);

	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

#ifdef SET_IO_MODE
	GMT_setmode (GMT_OUT);
#endif

	if (Ctrl->E.active) GMT_ECEF_init (&Ctrl->E.datum);
	if (Ctrl->F.active) unit = GMT_check_scalingopt ('F', Ctrl->F.unit, scale_unit_name);

	if (Ctrl->T.active && project_info.projection != GMT_LINEAR && project_info.region_supplied) {	/* Do datum shift & project coordinates */
		double_whammy = TRUE;
		if (Ctrl->I.active) {	/* Need to set the ellipsoid to that of the old datum */
			if (GMT_datum.from.ellipsoid_id < 0) {
				gmtdefs.ellipsoid = GMT_N_ELLIPSOIDS - 1;
				gmtdefs.ref_ellipsoid[i].eq_radius = GMT_datum.from.a;
				gmtdefs.ref_ellipsoid[i].flattening = GMT_datum.from.f;
			}
			else
				gmtdefs.ellipsoid = GMT_datum.from.ellipsoid_id;
		}
		else {	/* Need to set the ellipsoid to that of the new datum */
			if (GMT_datum.to.ellipsoid_id < 0) {
				gmtdefs.ellipsoid = GMT_N_ELLIPSOIDS - 1;
				gmtdefs.ref_ellipsoid[i].eq_radius = GMT_datum.to.a;
				gmtdefs.ref_ellipsoid[i].flattening = GMT_datum.to.f;
			}
			else
				gmtdefs.ellipsoid = GMT_datum.to.ellipsoid_id;
		}
	}
	else
		datum_conv_only = Ctrl->T.active;

	GMT_init_scales (unit, &fwd_scale, &inv_scale, &inch_to_unit, &unit_to_inch, unit_name);

	if (Ctrl->G.mode) {	/* save output format in case -J changes it */
		save[0] = GMT_io.out_col_type[0];
		save[1] = GMT_io.out_col_type[1];
	}
	u_scale = (Ctrl->I.active) ? inv_scale : fwd_scale;

	pure_ascii = !(GMT_io.binary[GMT_IN] || GMT_io.binary[GMT_OUT]);

	if (project_info.projection == GMT_NO_PROJ) {	/* Supply dummy linear proj */
		GMT_parse_J_option ("x1d");
		if (!project_info.region_supplied) {
			west = 0.0;	east = 360.0;
			south = -90.0;	north = 90.0;
		}
	}
	GMT_err_fail (GMT_map_setup (west, east, south, north), "");

	if (Ctrl->G.active && Ctrl->G.unit == 'c') GMT_io.in_col_type[0] = GMT_io.in_col_type[1] = GMT_IS_FLOAT;	/* Cartesian data only */
	if (Ctrl->L.active && Ctrl->L.unit == 'c') GMT_io.in_col_type[0] = GMT_io.in_col_type[1] = GMT_IS_FLOAT;	/* Cartesian data only */

	azimuth_func = (GMT_IS_SPHERICAL || !Ctrl->A.geodesic) ? GMT_az_backaz_sphere : GMT_az_backaz_geodesic;
	
	if (Ctrl->G.mode && proj_type < 2) {	/* Ensure we use the selected output coordinates */
		GMT_io.out_col_type[0] = save[0];
		GMT_io.out_col_type[1] = save[1];
	}

	if (gmtdefs.verbose && !(geodetic_calc || Ctrl->T.active)) {
		sprintf (format, "%s/%s/%s/%s", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		xmin = (Ctrl->C.active) ? project_info.xmin - project_info.x0 : project_info.xmin;
		xmax = (Ctrl->C.active) ? project_info.xmax - project_info.x0 : project_info.xmax;
		ymin = (Ctrl->C.active) ? project_info.ymin - project_info.y0 : project_info.ymin;
		ymax = (Ctrl->C.active) ? project_info.ymax - project_info.y0 : project_info.ymax;
		if (Ctrl->F.active) {	/* Convert to GMT inches */
			strncpy (unit_name, scale_unit_name, (size_t)GMT_TEXT_LEN);
			xmin /= project_info.x_scale;
			xmax /= project_info.x_scale;
			ymin /= project_info.y_scale;
			ymax /= project_info.y_scale;
		}

		/* Convert inches to chosen MEASURE */
		xmin *= inch_to_unit;
		xmax *= inch_to_unit;
		ymin *= inch_to_unit;
		ymax *= inch_to_unit;

		if (shift_xy) {
			xmin += Ctrl->C.easting;
			xmax += Ctrl->C.easting;
			ymin += Ctrl->C.northing;
			ymax += Ctrl->C.northing;
		}

		fprintf (stderr, "%s:  Transform ", GMT_program);
		fprintf (stderr, format, project_info.w, project_info.e, project_info.s, project_info.n);
		(Ctrl->I.active) ? fprintf (stderr, " <- ") : fprintf (stderr, " -> ");
		fprintf (stderr, format, xmin, xmax, ymin, ymax);
		fprintf (stderr, " [%s]\n", unit_name);
	}

	if (GMT_io.in_col_type[0] & GMT_IS_GEO && proj_type == 0) {	/* Geographic data */
		GMT_distance_func = (PFD) GMT_great_circle_dist;
		near_a_line = (PFI) GMT_near_a_line_spherical;
		greenwich = (west < 0.0 && east > 0.0);
	}
	else {
		GMT_distance_func = (PFD) GMT_cartesian_dist;
		near_a_line = (PFI) GMT_near_a_line_cartesian;
	}
	if (Ctrl->L.active) {
		GMT_import_table ((void *)Ctrl->L.file, GMT_IS_FILE, &xyline, 0.0, greenwich, FALSE, FALSE);
		if (proj_type == 2) {	/* Must convert the line points first */
			for (i = 0; i < xyline->n_segments; i++) {
				for (j = 0; j < xyline->segment[i]->n_rows; j++) {
					GMT_geo_to_xy (xyline->segment[i]->coord[GMT_X][j], xyline->segment[i]->coord[GMT_Y][j], &xtmp, &ytmp);
					xyline->segment[i]->coord[GMT_X][j] = xtmp;
					xyline->segment[i]->coord[GMT_Y][j] = ytmp;
				}
			}
		}
		else if (GMT_io.in_col_type[0] & GMT_IS_GEO && proj_type == 0 && !GMT_IS_SPHERICAL) {
			do_geo_conv = TRUE;
			/* Will need spherical trig so convert to geocentric latitudes */
			for (i = 0; i < xyline->n_segments; i++) {
				for (j = 0; j < xyline->segment[i]->n_rows; j++) {
					xyline->segment[i]->coord[GMT_Y][j] = GMT_lat_swap (xyline->segment[i]->coord[GMT_Y][j], GMT_LATSWAP_G2O);	/* Convert to geocentric */
				}
			}
		}
	}

	/* Now we are ready to take on some input values */

	x = (gmtdefs.xy_toggle[GMT_OUT]) ? 1 : 0;	y = 1 - x;		/* Set up which columns have x and y for output only*/
	if ((GMT_IS_MAPPING || Ctrl->E.active) && Ctrl->I.active) {
		GMT_io.out_col_type[0] = GMT_IS_LON;	GMT_io.out_col_type[1] = GMT_IS_LAT;	/* Inverse projection expects x,y and gives lon, lat */
		GMT_io.in_col_type[0] = GMT_io.in_col_type[1] = GMT_IS_FLOAT;
	}
	if (datum_conv_only) {	/* Both in and out are geographic */
		GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
		GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
		GMT_io.in_col_type[2] = GMT_io.out_col_type[2] = GMT_IS_FLOAT;
	}

	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : GMT_MAX_COLUMNS;
	if (GMT_io.ncol[GMT_OUT])	/* Want binary output to be limitied to first ncol[1] columns */
		n_output = &GMT_io.ncol[GMT_OUT];
	else	/* Default to the number of input fields */
		n_output = &n_expected_fields;

	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;
	n_args = (argc > 1) ? argc : 2;

	x_in_min = y_in_min = x_out_min = y_out_min = DBL_MAX;
	x_in_max = y_in_max = x_out_max = y_out_max = -DBL_MAX;

	two = (Ctrl->E.active || (Ctrl->T.active && GMT_datum.h_given)) ? 3 : 2;	/* # of output points from conversion */

	if (shift_xy && Ctrl->F.active) {	/* Use same units in -C and -F */
		if (Ctrl->I.active) {
			Ctrl->C.easting /= u_scale;
			Ctrl->C.northing /= u_scale;
		}
		else {
			Ctrl->C.easting *= u_scale;
			Ctrl->C.northing *= u_scale;
		}
	}
	if (Ctrl->G.mode >= 2 && proj_type == 2) {	/* Must project the fixed point here */
		GMT_geo_to_xy (Ctrl->G.lon, Ctrl->G.lat, &xtmp, &ytmp);
		if (Ctrl->C.active) {	/* Change origin from lower left to projection center */
			xtmp -= project_info.x0;
			ytmp -= project_info.y0;
		}
		if (Ctrl->F.active) {	/* Convert to 1:1 scale */
			xtmp /= project_info.x_scale;
			ytmp /= project_info.y_scale;
			if (unit) {
				xtmp *= u_scale;
				ytmp *= u_scale;
			}
		}
		else if (gmtdefs.measure_unit != GMT_INCH) {	/* Convert from inch to whatever */
			xtmp *= inch_to_unit;
			ytmp *= inch_to_unit;
		}
		if (shift_xy) {
			xtmp += Ctrl->C.easting;
			ytmp += Ctrl->C.northing;
		}
		Ctrl->G.lon = xtmp;
		Ctrl->G.lat = ytmp;
	}

	if (Ctrl->L.mode == 3)
		fmt[0] = fmt[1] = 2;
	else {
		fmt[0] = 0;
		fmt[1] = 1;
	}

	n = n_read_in_seg = 0;
	for (fno = 1; !done && fno < n_args; fno++) {	/* Loop over input files, if any */
		if (!nofile && argv[fno][0] == '-') continue;

		if (nofile) {	/* Just read standard input */
			fp = GMT_stdin;
			done = TRUE;
			if (gmtdefs.verbose) fprintf (stderr, "%s: Reading from standard input\n", GMT_program);
#ifdef SET_IO_MODE
			GMT_setmode (GMT_IN);
#endif
		}
		else if ((fp = GMT_fopen (argv[fno], GMT_io.r_mode)) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[fno]);
			continue;
		}

		if (!nofile && gmtdefs.verbose) fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[fno]);

		if (GMT_io.io_header[GMT_IN]) {
			for (i = 0; i < GMT_io.n_header_recs; i++) {
				GMT_fgets (line, BUFSIZ, fp);
				if (first && !GMT_io.binary[GMT_OUT] && GMT_io.io_header[GMT_OUT]) GMT_fputs (line, GMT_stdout);
			}
			first = FALSE;
		}
		if (GMT_IO_GAP_CHECKING) GMT_write_segmentheader (GMT_stdout, *n_output);

		if (Ctrl->I.active) {		/* Do inverse transformation */

			n_fields = GMT_input (fp, &n_expected_fields, &in);
			while (!GMT_REC_IS_EOF) {	/* Not yet EOF */

				while (GMT_REC_IS_NEW_SEGMENT && !GMT_REC_IS_EOF) {
					GMT_write_segmentheader (GMT_stdout, n_expected_fields);
					n_fields = GMT_input (fp, &n_expected_fields, &in);
					n_read_in_seg = 0;
				}
				if (GMT_REC_IS_GAP) {
					GMT_write_segmentheader (GMT_stdout, *n_output);
					GMT_io.status = 0;	/* Done with gap */
					line_start = TRUE;
					n_read_in_seg = 0;
				}
				if (GMT_REC_IS_EOF) continue;

				while (!GMT_REC_IS_LINE_BREAK && !GMT_REC_IS_EOF) {	/* Keep going until FALSE or = 2 segment header */

					n_read++;
					n_read_in_seg++;
					if (GMT_REC_IS_ERROR && n_fields < 2) {
						fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
						continue;
					}
					if (!out) out = (double *) GMT_memory (VNULL, (size_t)n_expected_fields, sizeof (double), GMT_program);

					if (gmtdefs.verbose) {
						x_in = in[GMT_X];
						y_in = in[GMT_Y];
					}
					if (shift_xy) {
						in[GMT_X] -= Ctrl->C.easting;
						in[GMT_Y] -= Ctrl->C.northing;
					}
					if (Ctrl->E.active) {
						GMT_ECEF_inverse (in, out);
					}
					else {
						if (Ctrl->F.active) {	/* Convert from 1:1 scale */
							if (unit) {
								in[GMT_X] *= u_scale;
								in[GMT_Y] *= u_scale;
							}
							in[GMT_X] *= project_info.x_scale;
							in[GMT_Y] *= project_info.y_scale;
						}
						else if (gmtdefs.measure_unit != GMT_INCH) {	/* Convert from whatever to inch */
							in[GMT_X] *= unit_to_inch;
							in[GMT_Y] *= unit_to_inch;
						}
						if (Ctrl->C.active) {	/* Then correct so lower left corner is (0,0) */
							in[GMT_X] += project_info.x0;
							in[GMT_Y] += project_info.y0;
						}
						GMT_xy_to_geo (&out[GMT_X], &out[GMT_Y], in[GMT_X], in[GMT_Y]);
					}
					if (double_whammy) {	/* Now apply datum shift */
						in[GMT_X] = out[GMT_X];
						in[GMT_Y] = out[GMT_Y];
						GMT_conv_datum (in, out);
					}
					if (Ctrl->S.active && GMT_map_outside (out[GMT_X], out[GMT_Y])) {
						n_fields = GMT_input (fp, &n_expected_fields, &in);
						continue;
					}
					if (gmtdefs.verbose) {
						x_in_min = MIN (x_in_min, x_in);
						x_in_max = MAX (x_in_max, x_in);
						y_in_min = MIN (y_in_min, y_in);
						y_in_max = MAX (y_in_max, y_in);
						x_out_min = MIN (x_out_min, out[GMT_X]);
						x_out_max = MAX (x_out_max, out[GMT_X]);
						y_out_min = MIN (y_out_min, out[GMT_Y]);
						y_out_max = MAX (y_out_max, out[GMT_Y]);
					}

					if (pure_ascii && n_expected_fields > 2) {
						/* Special case: ASCII i/o and at least 3 columns:
						  Columns beyond first two could be text strings */

						/* We will use GMT_strtok to step past the first 2 [or 3] columns.  The remainder
						* will then be the user text that we want to preserve.  Since strtok places
						* 0 to indicate start of next token we count our way to the start of the text. */

						strcpy (line, GMT_io.current_record);
						GMT_chop (line);
						pos = 0;
						GMT_strtok (line, " \t,", &pos, p);	/* Returns xstring and update pos */
						GMT_strtok (line, " \t,", &pos, p);	/* Returns ystring and update pos */
						GMT_ascii_output_one (GMT_stdout, out[x], 0);
						GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
						GMT_ascii_output_one (GMT_stdout, out[y], 1);
						GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
						if (Ctrl->E.active) {
							GMT_strtok (line, " \t,", &pos, p);	/* Returns zstring and update pos */
							GMT_ascii_output_one (GMT_stdout, out[2], 2);
							GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
						}
						GMT_fputs (&line[pos], GMT_stdout);	GMT_fputs("\n", GMT_stdout);	/* Start of user text */
					}
					else {	/* Simply copy other columns and output */
						for (k = two; k < *n_output; k++) out[k] = in[k];
						GMT_output (GMT_stdout, *n_output, out);
					}
					n++;
					if (long_verbose && (n%1000) == 0) fprintf (stderr, "%s: Projected %ld points\r", GMT_program, n);
					n_fields = GMT_input (fp, &n_expected_fields, &in);
				}
			}
		}
		else {		/* Do forward transformation */

			n_fields = GMT_input (fp, &n_expected_fields, &in);
			while (!GMT_REC_IS_EOF) {	/* Not yet EOF */

				while (GMT_REC_IS_NEW_SEGMENT && !GMT_REC_IS_EOF) {
					GMT_write_segmentheader (GMT_stdout, n_expected_fields);
					n_fields = GMT_input (fp, &n_expected_fields, &in);
					line_start = TRUE;
					n_read_in_seg = 0;
				}
				if (GMT_REC_IS_GAP) {
					GMT_write_segmentheader (GMT_stdout, *n_output);
					GMT_io.status = 0;	/* Done with gap */
					line_start = TRUE;
					n_read_in_seg = 0;
				}
				if (GMT_REC_IS_EOF) continue;

				while (!GMT_REC_IS_LINE_BREAK && !GMT_REC_IS_EOF) {	/* Keep going until FALSE or = 2 segment header */
					if (GMT_REC_IS_ERROR && n_fields < 2) {
						fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
						continue;
					}
					if (!out) out = (double *) GMT_memory (VNULL, (size_t)n_expected_fields+3, sizeof (double), GMT_program);
					n_read++;
					n_read_in_seg++;

					/* Because -: is processed in GMT_input we use [0] for lon and [1] for y always */

					if (Ctrl->S.active && GMT_map_outside (in[GMT_X], in[GMT_Y])) {
						n_fields = GMT_input (fp, &n_expected_fields, &in);
						continue;
					}

					if (datum_conv_only) {
						GMT_conv_datum (in, out);
					}
					else if (Ctrl->E.active) {
						GMT_ECEF_forward (in, out);
					}
					else {
						if (double_whammy) {	/* Apply datum shift first */
							GMT_conv_datum (in, out);
							in[GMT_X] = out[GMT_X];
							in[GMT_Y] = out[GMT_Y];
						}
						GMT_geo_to_xy (in[GMT_X], in[GMT_Y], &out[GMT_X], &out[GMT_Y]);
						if (Ctrl->C.active) {	/* Change origin from lower left to projection center */
							out[GMT_X] -= project_info.x0;
							out[GMT_Y] -= project_info.y0;
						}
						if (Ctrl->F.active) {	/* Convert to 1:1 scale */
							out[GMT_X] /= project_info.x_scale;
							out[GMT_Y] /= project_info.y_scale;
							if (unit) {
								out[GMT_X] *= u_scale;
								out[GMT_Y] *= u_scale;
							}
						}
						else if (gmtdefs.measure_unit != GMT_INCH) {	/* Convert from inch to whatever */
							out[GMT_X] *= inch_to_unit;
							out[GMT_Y] *= inch_to_unit;
						}
						if (shift_xy) {
							out[GMT_X] += Ctrl->C.easting;
							out[GMT_Y] += Ctrl->C.northing;
						}
					}
					if (gmtdefs.verbose) {
						x_in_min = MIN (x_in_min, in[GMT_X]);
						x_in_max = MAX (x_in_max, in[GMT_X]);
						y_in_min = MIN (y_in_min, in[GMT_Y]);
						y_in_max = MAX (y_in_max, in[GMT_Y]);
					}
					if (geodetic_calc) {	/* Get either distances or azimuths */
						if (Ctrl->G.mode) {	/* Cumulative distances along track */
							if (Ctrl->G.mode == 4) {
								if (proj_type == 2) {	/* Calculate Cartesian distances using projected units */
									GMT_geo_to_xy (in[2], in[3], &xtmp, &ytmp);
									s = hypot (xtmp - out[GMT_X], ytmp - out[GMT_Y]);
								}
								else
									s = d_scale * (*distance_func) (in[GMT_X], in[GMT_Y], in[2], in[3]);
							}
							else if (Ctrl->G.mode >= 2 && line_start)
								s = d = 0.0;
							else if (proj_type == 2)	/* Calculate Cartesian distances using projected units */
								s = hypot (Ctrl->G.lon - out[GMT_X], Ctrl->G.lat - out[GMT_Y]);
							else if (proj_type == 1)	/* Plain Cartesian distances using input points */
								s = hypot (Ctrl->G.lon - in[GMT_X], Ctrl->G.lat - in[GMT_Y]);
							else				/* Great circle distances */
								s = d_scale * (*distance_func) (Ctrl->G.lon, Ctrl->G.lat, in[GMT_X], in[GMT_Y]);
							if (Ctrl->G.mode >= 2) {
								line_start = FALSE;
								if (Ctrl->G.mode >= 3)	/* Increments */
									d = s;
								else
									d += s;		/* Cumulative */
								if (proj_type == 2) {	/* Calculate distances using projected units */
									Ctrl->G.lon = out[GMT_X];
									Ctrl->G.lat = out[GMT_Y];
								}
								else {
									Ctrl->G.lon = in[GMT_X];
									Ctrl->G.lat = in[GMT_Y];
								}
							}
							else
								d = s;
						}
						else if (Ctrl->L.active) {	/* Compute closest distance to line */
							if (proj_type == 2)	/* Using projected coordinates */
								(void) near_a_line (out[GMT_X], out[GMT_Y], xyline, Ctrl->L.mode, &d, &xnear, &ynear);
							else {			/* Using input coordinates */
								y_in = (do_geo_conv) ? GMT_lat_swap (in[GMT_Y], GMT_LATSWAP_G2O) : in[GMT_Y];				/* Convert to geocentric */
								(void) near_a_line (in[GMT_X], y_in, xyline, Ctrl->L.mode, &d, &xnear, &ynear);
								if (do_geo_conv && Ctrl->L.mode != 3) ynear = GMT_lat_swap (ynear, GMT_LATSWAP_O2G);		/* Convert back to geodetic */
							}
							d *= d_scale;
						}
						else if (Ctrl->A.azims) {	/* Azimuth from previous point */
							if (n_read_in_seg == 1) {	/* First point has undefined azimuth since there is no previous point */
								d = GMT_d_NaN;
							}
							else {
								d = (*azimuth_func) (lon_prev, lat_prev, in[GMT_X], in[GMT_Y], Ctrl->A.reverse);
							}
							lon_prev = in[GMT_X];
							lat_prev = in[GMT_Y];
						}
						else {	/* Azimuths with respect to a fixed point */
							d = (*azimuth_func) (Ctrl->A.lon, Ctrl->A.lat, in[GMT_X], in[GMT_Y], Ctrl->A.reverse);
						}
						if (pure_ascii && n_expected_fields > 2) {
							/* Special case: Ascii input and at least 3 columns:
							 Columns beyond first two could be text strings */

							/* We will use GMT_strtok to step past the first 2 [or 3] columns.  The remainder
							* will then be the user text that we want to preserve.  Since strtok places
							* 0 to indicate start of next token we count our way to the start of the text. */

							strcpy (line, GMT_io.current_record);
							GMT_chop (line);
							pos = 0;
							GMT_strtok (line, " \t,", &pos, p);	/* Returns xstring */
							GMT_strtok (line, " \t,", &pos, p);	/* Returns ystring */

							GMT_ascii_output_one (GMT_stdout, in[x], 0);
							GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
							GMT_ascii_output_one (GMT_stdout, in[y], 1);
							GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
							GMT_fputs (&line[pos], GMT_stdout);
							GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
							GMT_ascii_output_one (GMT_stdout, d, 2);
							if (Ctrl->L.active) {
								GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
								GMT_ascii_output_one (GMT_stdout, xnear, fmt[0]);
								GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
								GMT_ascii_output_one (GMT_stdout, ynear, fmt[1]);
							}
							GMT_fputs ("\n", GMT_stdout);
						}
						else {	/* Simply copy other columns and output */
							for (k = 0; k < *n_output; k++) out[k] = in[k];
							n_out = *n_output;
							out[n_out++] = d;
							if (Ctrl->L.active) {
								out[n_out++] = xnear;
								out[n_out++] = ynear;
							}
							GMT_output (GMT_stdout, n_out, out);
						}
					}
					else {
						if (gmtdefs.verbose) {
							x_out_min = MIN (x_out_min, out[GMT_X]);
							x_out_max = MAX (x_out_max, out[GMT_X]);
							y_out_min = MIN (y_out_min, out[GMT_Y]);
							y_out_max = MAX (y_out_max, out[GMT_Y]);
						}
						if (pure_ascii && n_expected_fields > 2) {
							/* Special case: Ascii input and at least 3 columns:
							  Columns beyond first two could be text strings */

							/* We will use GMT_strtok to step past the first 2 [or 3] columns.  The remainder
							* will then be the user text that we want to preserve.  Since strtok places
							* 0 to indicate start of next token we count our way to the start of the text. */

							strcpy (line, GMT_io.current_record);
							GMT_chop (line);
							pos = 0;
							GMT_strtok (line, " \t,", &pos, p);	/* Returns xstring and update pos */
							GMT_strtok (line, " \t,", &pos, p);	/* Returns ystring and update pos */
							GMT_ascii_output_one (GMT_stdout, out[x], 0);
							GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
							GMT_ascii_output_one (GMT_stdout, out[y], 1);
							GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
							if (Ctrl->E.active || (Ctrl->T.active && GMT_datum.h_given)) {
								GMT_strtok (line, " \t,", &pos, p);	/* Returns zstring and update pos */
								GMT_ascii_output_one (GMT_stdout, out[2], 2);
								GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
							}
							GMT_fputs (&line[pos], GMT_stdout);	GMT_fputs("\n", GMT_stdout);
						}
						else {	/* Simply copy other columns and output */
							for (k = two; k < *n_output; k++) out[k] = in[k];
							GMT_output (GMT_stdout, *n_output, out);
						}
					}
					n++;
					if (long_verbose && (n%1000) == 0) fprintf (stderr, "%s: Projected %ld points\r", GMT_program, n);
					n_fields = GMT_input (fp, &n_expected_fields, &in);
				}
			}
		}
		if (fp != GMT_stdin) GMT_fclose (fp);
	}

	if (gmtdefs.verbose && n_read > 0) {
		fprintf (stderr, "%s: Projected %ld points\n", GMT_program, n);
		sprintf (format, "%%s: Input extreme values:  Xmin: %s Xmax: %s Ymin: %s Ymax %s\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		fprintf (stderr, format, GMT_program, x_in_min, x_in_max, y_in_min, y_in_max);
		if (!geodetic_calc) {
			sprintf (format, "%%s: Output extreme values:  Xmin: %s Xmax: %s Ymin: %s Ymax %s\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
			fprintf (stderr, format, GMT_program, x_out_min, x_out_max, y_out_min, y_out_max);
			if (Ctrl->I.active) {
				if (Ctrl->E.active)
					fprintf (stderr, "%s: Mapped %ld ECEF coordinates [m] to (lon,lat,h)\n", GMT_program, n);
				else
					fprintf (stderr, "%s: Mapped %ld x-y pairs [%s] to lon-lat\n", GMT_program, n, unit_name);
			}
			else if (Ctrl->T.active && GMT_datum.h_given)
				fprintf (stderr, "%s: Datum-converted %ld (lon,lat,h) triplets\n", GMT_program, n);
			else if (Ctrl->T.active)
				fprintf (stderr, "%s: Datum-converted %ld (lon,lat) pairs\n", GMT_program, n);
			else if (Ctrl->E.active)
				fprintf (stderr, "%s: Mapped %ld (lon,lat,h) triplets to ECEF coordinates [m]\n", GMT_program, n);
			else if (GMT_IS_MAPPING)
				fprintf (stderr, "%s: Mapped %ld lon-lat pairs to x-y [%s]\n", GMT_program, n, unit_name);
			else
				fprintf (stderr, "%s: Mapped %ld data pairs to x-y [%s]\n", GMT_program, n, unit_name);
		}
		if (Ctrl->S.active && n != n_read) fprintf (stderr, "%s: %ld fell outside region\n", GMT_program, n_read - n);
	}

	if (out) GMT_free ((void *)out);

	if (Ctrl->L.active) GMT_free_table (xyline);

	Free_mapproject_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_mapproject_Ctrl () {	/* Allocate and initialize a new control structure */
	struct MAPPROJECT_CTRL *C;
	
	C = (struct MAPPROJECT_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct MAPPROJECT_CTRL), "New_mapproject_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	
	C->L.mode = 2;
	
	return ((void *)C);
}

void Free_mapproject_Ctrl (struct MAPPROJECT_CTRL *C) {	/* Deallocate control structure */
	if (C->L.file) free ((void *)C->L.file);	
	GMT_free ((void *)C);	
}
