/*--------------------------------------------------------------------
 *	$Id: project.c 9923 2012-12-18 20:45:53Z pwessel $
 *
 *	Copyright (c) 1991-2013 by P. Wessel and W. H. F. Smith
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
 * project.c
 * reads (x,y,[z]) data and writes some combination of (x,y,z,p,q,u,v),
 * where p,q is the distance along,across the track of the projection of (x,y),
 * and u,v are the un-transformed (x,y) coordinates of the projected position.
 * Can also create (x,y) along track.
  
   Author: 	Walter H. F. Smith
   Date:	19 April, 1988.
   Modified:	4 December 1988, to be more flexible.
   Complete rebuild 22 June, 1989 to use vector products and do more things.
   version 2.0
   		23-FEB-1998	PW: Added support for multiple files, multi-segment formats
				and binary i/o.  Old -M renamed -Q.
		03-NOV-1998	PW: Can read any number of data columns; z in -Fz refers to
				all these columns in the output.
   Version:	3.4		PW: Fixed problem with small circle distances
   Version:	4.1.x
*/

#include "gmt.h"

#define PROJECT_N_FARGS	7

struct PROJECT_CTRL {	/* All control options for this program (except common args) */
	/* active is TRUE if the option has been activated */
	struct A {	/* -A<azimuth> */
		GMT_LONG active;
		double azimuth;
	} A;
	struct C {	/* -C<ox>/<oy> */
		GMT_LONG active;
		double x, y;
	} C;
	struct D {	/* -D<d_or_g> */
		GMT_LONG active;
		char set;
	} D;
	struct E {	/* -E<bx>/<by> */
		GMT_LONG active;
		double x, y;
	} E;
	struct F {	/* -F<flags> */
		GMT_LONG active;
		char col[PROJECT_N_FARGS];	/* Character codes for desired output in the right order */
	} F;
	struct G {	/* -G<inc> */
		GMT_LONG active;
		double inc;
	} G;
	struct L {	/* -L[w][<l_min>/<l_max>] */
		GMT_LONG active;
		GMT_LONG constrain;
		double min, max;
	} L;
	struct N {	/* -N */
		GMT_LONG active;
	} N;
	struct Q {	/* -Q */
		GMT_LONG active;
	} Q;
	struct S {	/* -S */
		GMT_LONG active;
	} S;
	struct T {	/* -T<px>/<py> */
		GMT_LONG active;
		double x, y;
	} T;
	struct W {	/* -W<w_min>/<w_max> */
		GMT_LONG active;
		double min, max;
	} W;
};

struct PROJECT_DATA {
        double  a[6];
	double *z;
	char *t;
};

int main (int argc, char **argv)
{
	GMT_LONG	i, n_used, n_total_read, n_total_used = 0, n_alloc = GMT_CHUNK;
	GMT_LONG	j, k, n_outputs;
	GMT_LONG	n_files = 0, fno, n_args, n_fields, n_expected_fields;
	GMT_LONG	n_z = 0, n_items, output_choice[PROJECT_N_FARGS];

	double	xx, yy, cos_theta, sin_theta, sin_lat_to_pole = 1.0;
	double	theta = 0.0, d_along, *in = NULL, *out = (double *)NULL;

	double	a[3], b[3], x[3], xt[3], pole[3], center[3], e[9];

	GMT_LONG	dateline = FALSE, error = FALSE, find_new_point = FALSE, first = TRUE, z_first = TRUE, pure_ascii, skip;
	GMT_LONG greenwich = FALSE, nofile = TRUE, done = FALSE, want_z_output = FALSE;

	FILE *fp = NULL;

	char	modifier, record_str[BUFSIZ], heading[PROJECT_N_FARGS][GMT_TEXT_LEN], buffer[GMT_TEXT_LEN];
	char	txt_a[GMT_LONG_TEXT], txt_b[GMT_LONG_TEXT];

	struct PROJECT_DATA *p_data = NULL;
	struct PROJECT_CTRL *Ctrl = NULL;

	int	compare_distances (const void *point_1, const void *point_2);
	double	oblique_setup (double plat, double plon, double *p, double *clat, double *clon, double *c, GMT_LONG c_given, GMT_LONG generate);
	void	oblique_transform (double xlat, double xlon, double *x_t_lat, double *x_t_lon, double *p, double *c);
	void	make_euler_matrix (double *p, double *e, double theta);
	void	matrix_3v (double *a, double *x, double *b);
	void	matrix_2v (double *a, double *x, double *b);
	void	sphere_project_setup (double alat, double alon, double *a, double blat, double blon, double *b, double azim, double *p, double *c, GMT_LONG two_pts);
	void	flat_project_setup (double alat, double alon, double blat, double blon, double plat, double plon, double *azim, double *e, GMT_LONG two_pts, GMT_LONG pole_set);
	void	copy_text_from_col3 (char *line, char *z_cols);
	void	*New_project_Ctrl (), Free_project_Ctrl (struct PROJECT_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct PROJECT_CTRL *) New_project_Ctrl ();		/* Allocate and initialize defaults in a new control structure */

	memset ((void *)output_choice, 0, PROJECT_N_FARGS * sizeof (GMT_LONG));

	/* New in GMT 4: Must process -N before the rest to ensure col_types are set properly (there is no -H here)  */
	for (i = 1; i < argc; i++) if (!strcmp (argv[i], "-N")) Ctrl->N.active = TRUE;
	if (!Ctrl->N.active) {
		GMT_io.in_col_type[0] = GMT_IS_LON;
		GMT_io.in_col_type[1] = GMT_IS_LAT;
	}

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {


				/* Common parameters */

				case 'H':
				case 'M':
				case 'V':
				case ':':
				case 'b':
				case 'f':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'F':
					Ctrl->F.active = TRUE;
					for (j = 2, k = 0; argv[i][j]; j++, k++) {
						if (k < PROJECT_N_FARGS)
							Ctrl->F.col[k] = argv[i][j];
						else {
							error++;
							fprintf (stderr, "%s: GMT SYNTAX ERROR -F option: Too many output columns selected\n", GMT_program);
							fprintf (stderr, "%s: GMT SYNTAX ERROR -F option: Choose from -Fxyzpqrs\n", GMT_program);
						}
					}
					break;
				case 'A':
					Ctrl->A.active = TRUE;
					Ctrl->A.azimuth = atof(&argv[i][2]);
					break;
				case 'C':
					Ctrl->C.active = TRUE;
					if (sscanf(&argv[i][2], "%[^/]/%s", txt_a, txt_b) != 2) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -C option.  Correct syntax: -C<lon0>/<lat0>\n", GMT_program);
						error++;
					}
					else {
						error += GMT_verify_expectations (GMT_io.in_col_type[0], GMT_scanf_arg (txt_a, GMT_io.in_col_type[0], &Ctrl->C.x), txt_a);
						error += GMT_verify_expectations (GMT_io.in_col_type[1], GMT_scanf_arg (txt_b, GMT_io.in_col_type[1], &Ctrl->C.y), txt_b);
						if (error) fprintf (stderr, "%s: GMT SYNTAX ERROR -C option:  Undecipherable argument %s\n", GMT_program, &argv[i][2]);
					}
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					Ctrl->D.set = argv[i][2];
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					if (sscanf(&argv[i][2], "%[^/]/%s", txt_a, txt_b) != 2) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -E option.  Correct syntax: -E<lon1>/<lat1>\n", GMT_program);
						error++;
					}
					else {
						error += GMT_verify_expectations (GMT_io.in_col_type[0], GMT_scanf_arg (txt_a, GMT_io.in_col_type[0], &Ctrl->E.x), txt_a);
						error += GMT_verify_expectations (GMT_io.in_col_type[1], GMT_scanf_arg (txt_b, GMT_io.in_col_type[1], &Ctrl->E.y), txt_b);
						if (error) fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Undecipherable argument %s\n", GMT_program, &argv[i][2]);
					}
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.inc = atof(&argv[i][2]);
					break;
				case 'L':
					Ctrl->L.active = TRUE;
					modifier = argv[i][2]; 
					if (modifier == 'W' || modifier == 'w') {
						Ctrl->L.constrain = TRUE;
					}
					else {
						if (sscanf(&argv[i][2], "%lf/%lf", &Ctrl->L.min, &Ctrl->L.max) != 2) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR -L option.  Correct syntax: -L[w | <min>/<max>]\n", GMT_program);
							error++;
						}
					}
					break;
				case 'N': /* Handled above but still in argv */
					Ctrl->N.active = TRUE;
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					if (sscanf(&argv[i][2], "%[^/]/%s", txt_a, txt_b) != 2) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -T option.  Correct syntax: -T<lonp>/<latp>\n", GMT_program);
						error++;
					}
					else {
						error += GMT_verify_expectations (GMT_io.in_col_type[0], GMT_scanf_arg (txt_a, GMT_io.in_col_type[0], &Ctrl->T.x), txt_a);
						error += GMT_verify_expectations (GMT_io.in_col_type[1], GMT_scanf_arg (txt_b, GMT_io.in_col_type[1], &Ctrl->T.y), txt_b);
						if (error) fprintf (stderr, "%s: GMT SYNTAX ERROR -T option:  Undecipherable argument %s\n", GMT_program, &argv[i][2]);
					}
					break;
				case 'W':
					Ctrl->W.active = TRUE;
					if (sscanf(&argv[i][2], "%lf/%lf", &Ctrl->W.min, &Ctrl->W.max) != 2) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -W option.  Correct syntax: -W<min>/<max>\n", GMT_program);
						error++;
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
		fprintf (stderr, "project %s - Project data onto line or great circle, generate track, or translate coordinates\n\n", GMT_VERSION);
		fprintf(stderr,"usage:	project [files] -C<ox>/<oy> [-A<azimuth>] [-D<d_or_g>] [-E<bx>/<by>]\n");
		fprintf(stderr,"\t[-F<flags>] [-G<dist>] [%s] [-L[w][<l_min>/<l_max>]]\n", GMT_H_OPT);
		fprintf(stderr,"\t[-N] [-Q] [-S] [-T<px>/<py>] [-V] [-W<w_min>/<w_max>]\n");
		fprintf(stderr,"\t[%s] [%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_b_OPT, GMT_f_OPT, GMT_m_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf(stderr,"\tproject will read stdin or file, and does not want input if -G option.\n");
		fprintf(stderr,"\tThe projection may be defined in (only) one of three ways:\n");
		fprintf(stderr,"\t  (1) by a center -C and an azimuth -A,\n");
		fprintf(stderr,"\t  (2) by a center -C and end point of the path -E,\n");
		fprintf(stderr,"\t  (3) by a center -C and a roTation pole position -T.\n");
		fprintf(stderr,"\t  In a spherical projection [default], all cases place the central meridian\n");
		fprintf(stderr,"\t  of the transformed coordinates (p,q) through -C (p = 0 at -C).  The equator\n");
		fprintf(stderr,"\t  of the (p,q) system (line q = 0) passes through -C and makes an angle\n");
		fprintf(stderr,"\t  <azimuth> with North (case 1), or passes through -E (case 2), or is\n");
		fprintf(stderr,"\t  determined by the pole -T (case 3).  In (3), point -C need not be on equator.\n");
		fprintf(stderr,"\t  In a cartesian [-N option] projection, p = q = 0 at -O in all cases;\n");
		fprintf(stderr,"\t  (1) and (2) orient the p axis, while (3) orients the q axis.\n\n");
		fprintf(stderr,"\t-C<ox>/<oy> sets the location of the center.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf(stderr,"\t-A<azimuth> sets the option (1) Azimuth, (degrees CW from North).\n");
		fprintf(stderr,"\t-D will force the location of the Discontinuity in the r coordinate;\n");
		fprintf(stderr,"\t  -Dd (dateline) means [-180 < r < 180], -Dg (greenwich) means [0 < r < 360].\n");
		fprintf(stderr,"\t  The default does not check; in spherical case this usually results in [-180,180].\n");
		fprintf(stderr,"\t-E<bx>/<by> sets the option (2) location of end point E.\n");
		fprintf(stderr,"\t-Fflags: Indicate what output you want as one or more of xyzpqrs in any order;\n");
		fprintf(stderr,"\t  where x,y,[z] refer to input data locations and optional values,\n");
		fprintf(stderr,"\t  p,q are the coordinates of x,y in the projection's coordinate system,\n");
		fprintf(stderr,"\t  r,s is the projected position of x,y (taking q = 0) in the (x,y) coordinate system.\n");
		fprintf(stderr,"\t  p,q may be scaled from degrees into kilometers by the -Q option.  See -L, -Q, -W.\n");
		fprintf(stderr,"\t  Note z refers to all input data columns beyond the required x,y\n");
		fprintf(stderr,"\t  [Default is all fields, i.e., -Fxyzpqrs].\n");
		fprintf(stderr,"\t  If -G is set, -F is not available and output defaults to rsp.\n");
		fprintf(stderr,"\t-G means Generate (r,s,p) points along profile every <dist> units (No input data used).\n");
		fprintf(stderr,"\t   If E given, will generate from C to E; else must give -L<l_min>/<l_max> for length.\n");
		GMT_explain_option ('H');
		fprintf(stderr,"\t-L Check the Length along the projected track and use only certain points.\n");
		fprintf(stderr,"\t  -Lw will use only those points Within the span from C to E (Must have set -E).\n");
		fprintf(stderr,"\t  -L<l_min>/<l_max> will only use points whose p is [l_min <= p <= l_max].\n");
		fprintf(stderr,"\t  Default uses all points.  Note p = 0 at C and increases toward E in azim direction.\n");
		fprintf(stderr,"\t-N means Flat_earth; a cartesian projection is made.  Default is spherical.\n");
		fprintf(stderr,"\t-Q means convert to Map units, so x,y,r,s are degrees,\n");
		fprintf(stderr,"\t  while p,q,dist,l_min,l_max,w_min,w_max are km.\n");
		fprintf(stderr,"\t  If not set, then p,q,dist,l_min,l_max,w_min,w_max are assumed to be in same units as x,y,r,s.\n");
		fprintf(stderr,"\t-S means the output should be Sorted into increasing p value.\n");
		fprintf(stderr,"\t-T<px>/<py> sets the option (3) location of the roTation pole to the projection.\n");
		GMT_explain_option ('V');
		fprintf(stderr,"\t-W Check the width across the projected track and use only certain points.\n");
		fprintf(stderr,"\t  This will use only those points whose q is [w_min <= q <= w_max].\n");
		fprintf(stderr,"\t  Note that q is positive to your LEFT as you walk from C toward E in azim direction.\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf(stderr,"\t  Default is 2 input columns (x, y).\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (Ctrl->L.active && !Ctrl->L.constrain && Ctrl->L.min >= Ctrl->L.max) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -L option.  w_min must be < w_max\n", GMT_program);
		error++;
	}
	if (Ctrl->W.active && Ctrl->W.min >= Ctrl->W.max) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -W option.  w_min must be < w_max\n", GMT_program);
		error++;
	}
	if ((Ctrl->A.active + Ctrl->E.active + Ctrl->T.active) > 1) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR: Specify only one of -A, -E, and -T\n", GMT_program);
		error++;
	}
	if (Ctrl->E.active && (Ctrl->C.x == Ctrl->E.x) && (Ctrl->C.y == Ctrl->E.y)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option: Second point must differ from origin!\n", GMT_program);
		error++;
	}
	if (Ctrl->G.active && Ctrl->L.min == Ctrl->L.max && !Ctrl->E.active) {	/* We don't know how long to generate  */
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option: Must also specify -Lmin/max or use -E instead\n", GMT_program);
		error++;
	}
	if (Ctrl->G.active && Ctrl->F.active) {	/* -F not allowed with -G  */
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option: -F not allowed [Defaults to rsp]\n", GMT_program);
		error++;
	}
	if (Ctrl->G.active && Ctrl->G.inc <= 0.0) {	/* No increment given  */
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option: Must specify a positive increment\n", GMT_program);
		error++;
	}
	if (Ctrl->L.constrain && !Ctrl->E.active) {	/* Same problem.  */
		fprintf (stderr, "%s: GMT SYNTAX ERROR -L option: Must specify -Lmin/max or use -E instead\n", GMT_program);
		error++;
	}
	/* Convert user's -F choices to internal parameters */
	for (k = n_outputs = 0; k < PROJECT_N_FARGS && Ctrl->F.col[k]; k++) {
		switch (Ctrl->F.col[k]) {
			case 'z':	/* Special flag, can mean any number of z columns */
				output_choice[k] = -1;
				want_z_output = TRUE;
				break;
			case 'x':
				output_choice[k] = 0;
				break;
			case 'y':
				output_choice[k] = 1;
				break;
			case 'p':
				output_choice[k] = 2;
				break;
			case 'q':
				output_choice[k] = 3;
				break;
			case 'r':
				output_choice[k] = 4;
				find_new_point = TRUE;
				break;
			case 's':
				output_choice[k] = 5;
				find_new_point = TRUE;
				break;
			default:
				fprintf (stderr, "%s: GMT SYNTAX ERROR -F option:  Unrecognized choice %c\n", GMT_program, Ctrl->F.col[k]);
				error = TRUE;
		}
		n_outputs++;
	}

	if (n_outputs > PROJECT_N_FARGS) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -F option: Too many output columns selected (%ld)\n", GMT_program, n_outputs);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 2) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 2 columns\n", GMT_program);
		error++;
	}
	if (Ctrl->D.active && Ctrl->D.set) {
		if (Ctrl->D.set == 'D' || Ctrl->D.set == 'd')
			dateline = TRUE;
		else if (Ctrl->D.set == 'g' || Ctrl->D.set == 'G')
			greenwich = TRUE;
		else {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -D option:  Unrecognized modifier %c\n", GMT_program, Ctrl->D.set);
			error = TRUE;
		}
	}
	
	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

#ifdef SET_IO_MODE
	GMT_setmode (GMT_OUT);
#endif

	pure_ascii = !(GMT_io.binary[GMT_IN] || GMT_io.binary[GMT_OUT]);

	if (n_outputs == 0 && !(Ctrl->G.active) ) {	/* Generate default -F setting (all) */
		n_outputs = PROJECT_N_FARGS;
		for (i = 0; i < 2; i++) output_choice[i] = (int)i;
		output_choice[2] = -1;
		for (i = 3; i < n_outputs; i++) output_choice[i] = (int)i - 1;
		find_new_point = TRUE;
	}

	p_data = (struct PROJECT_DATA *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (struct PROJECT_DATA), GMT_program);

	if (Ctrl->G.active && Ctrl->E.active && (Ctrl->L.min == Ctrl->L.max) ) Ctrl->L.constrain = TRUE;	/* Default generate from A to B  */

	/* Set up rotation matrix e for flat earth, or pole and center for spherical; get Ctrl->L.min, Ctrl->L.max if stay_within  */

	if (Ctrl->N.active) {
		theta = Ctrl->A.azimuth;
		flat_project_setup(Ctrl->C.y, Ctrl->C.x, Ctrl->E.y, Ctrl->E.x, Ctrl->T.y, Ctrl->T.x, &theta, e, Ctrl->E.active, Ctrl->T.active);
		/* Azimuth (theta) is now cartesian in degrees */
		if (Ctrl->L.constrain) {
			Ctrl->L.min = 0.0;
			xx = Ctrl->E.x - Ctrl->C.x;
			yy = Ctrl->E.y - Ctrl->C.y;
			Ctrl->L.max = d_sqrt(xx*xx + yy*yy);
			if (Ctrl->Q.active) Ctrl->L.max *= project_info.DIST_KM_PR_DEG;
		}
	}
	else {
		if (Ctrl->T.active) {
			sin_lat_to_pole = oblique_setup(Ctrl->T.y, Ctrl->T.x, pole, &Ctrl->C.y, &Ctrl->C.x, center, Ctrl->C.active, Ctrl->G.active);
		}
		else {
			sphere_project_setup(Ctrl->C.y, Ctrl->C.x, a, Ctrl->E.y, Ctrl->E.x, b, Ctrl->A.azimuth, pole, center, Ctrl->E.active);
		}
		if (Ctrl->L.constrain) {
			Ctrl->L.min = 0.0;
			Ctrl->L.max = GMT_dot3v(a,b);
			Ctrl->L.max = d_acosd(Ctrl->L.max);
			if (Ctrl->Q.active) Ctrl->L.max *= project_info.DIST_KM_PR_DEG;
		}
	}

	/* Now things are initialized.  We will work in degrees for awhile, so we convert things:  */

	if (Ctrl->Q.active) {
		Ctrl->G.inc /= project_info.DIST_KM_PR_DEG;
		Ctrl->L.min /= project_info.DIST_KM_PR_DEG;
		Ctrl->L.max /= project_info.DIST_KM_PR_DEG;
		Ctrl->W.min /= project_info.DIST_KM_PR_DEG;
		Ctrl->W.max /= project_info.DIST_KM_PR_DEG;
	}

	/*  Now we are ready to work  */

	n_used = n_total_read = 0;

	if (Ctrl->G.active) {	/* Not input data expected, just generate track from arguments given */
		n_outputs = 3;
		output_choice[0] = 4;
		output_choice[1] = 5;
		output_choice[2] = 2;
		out = (double *) GMT_memory (VNULL, (size_t)n_outputs, sizeof (double), GMT_program);

		d_along = Ctrl->L.min;
		while ((Ctrl->L.max - d_along) > (GMT_CONV_LIMIT*Ctrl->G.inc)) {
			p_data[n_used].a[2] = d_along;
			p_data[n_used].t = NULL;	/* Initialize since that is not done by realloc */
			p_data[n_used].z = NULL;	/* Initialize since that is not done by realloc */
			n_used++;
			d_along = Ctrl->L.min + n_used * Ctrl->G.inc;
			if (n_used == (n_alloc-1)) {
				n_alloc <<= 1;
				p_data = (struct PROJECT_DATA *) GMT_memory ((void *)p_data, (size_t)n_alloc, sizeof (struct PROJECT_DATA), GMT_program);
			}
		}
		p_data[n_used].a[2] = Ctrl->L.max;
		p_data[n_used].t = NULL;	/* Initialize since that is not done by realloc */
		p_data[n_used].z = NULL;	/* Initialize since that is not done by realloc */
		n_used ++;

		/* We need to find r,s  */

		if (Ctrl->N.active) {
			sincosd (90.0 + theta, &sin_theta, &cos_theta);
			for (i = 0; i < n_used; i++) {
				p_data[i].a[4] = Ctrl->C.x + p_data[i].a[2] * cos_theta;
				p_data[i].a[5] = Ctrl->C.y + p_data[i].a[2] * sin_theta;
				while (greenwich && p_data[i].a[4] < 0.0) p_data[i].a[4] += 360.0;
				while (dateline && p_data[i].a[4] > 180.0) p_data[i].a[4] -= 360.0;
			}
		}
		else {
			double C[3], N[3];
			GMT_geo_to_cart(Ctrl->C.y, Ctrl->C.x, C, TRUE);
			GMT_cross3v (pole, C, N);		/* This is vector normal to meridian plan */
			GMT_normalize3v (N);			/* Make it a unit vector */
			make_euler_matrix (N, e, 90.0);		/* Rotation matrix e about N */
			matrix_3v (e, pole, x);			/* x is the generating vector for our circle */
			for (i = 0; i < n_used; i++) {
				make_euler_matrix (pole, e, p_data[i].a[2]);
				matrix_3v(e,x,xt);
				GMT_cart_to_geo(&(p_data[i].a[5]), &(p_data[i].a[4]), xt, TRUE);
				while (greenwich && p_data[i].a[4] < 0.0) p_data[i].a[4] += 360.0;
				while (dateline && p_data[i].a[4] > 180.0) p_data[i].a[4] -= 360.0;
			}
		}

		/* At this stage, all values are still in degrees.  */

		if (Ctrl->Q.active) {
			for (i = 0; i < n_used; i++) {
				p_data[i].a[2] *= project_info.DIST_KM_PR_DEG;
				p_data[i].a[3] *= project_info.DIST_KM_PR_DEG;
			}
		}

		/* Now output generated track */

		if (!GMT_io.binary[GMT_OUT]) {
			if (GMT_io.io_header[GMT_OUT]) {sprintf (buffer, "lon%slat%sdist\n", gmtdefs.field_delimiter, gmtdefs.field_delimiter);	GMT_fputs(buffer, GMT_stdout);}
		}

		for (i = 0; i < n_used; i++) {
			for (j = 0; j < n_outputs; j++) out[j] = p_data[i].a[output_choice[j]];
			GMT_output (GMT_stdout, n_outputs, out);
		}
	}

	else {	/* Must read input file */

		n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : GMT_MAX_COLUMNS;

		if (n_files > 0)
			nofile = FALSE;
		else
			n_files = 1;

		n_args = (argc > 1) ? argc : 2;

		for (fno = 1; !done && fno < n_args; fno++) {	/* Loop over input files, if any */
			if (!nofile && argv[fno][0] == '-') continue;

			if (nofile) {	/* Just read standard input */
				fp = GMT_stdin;
				done = TRUE;
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
				GMT_fgets (record_str, BUFSIZ, fp);
				sscanf (record_str, "%s %s %s", heading[0], heading[1], heading[6]);
				if (! (heading[6]) ) strcpy (heading[6],"Z");
				strcpy (heading[2],"p");
				strcpy (heading[3],"q");
				strcpy (heading[4],"r");
				strcpy (heading[5],"s");
				for (i = 1; i < GMT_io.n_header_recs; i++) GMT_fgets (record_str, BUFSIZ, fp);
			}

			n_fields = GMT_input (fp, &n_expected_fields, &in);
			n_used = 0;

			while (! (GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

				while (GMT_io.status & GMT_IO_SEGMENT_HEADER) {
					GMT_write_segmentheader (GMT_stdout, n_expected_fields);
					n_fields = GMT_input (fp, &n_expected_fields, &in);
				}
				if (z_first) {
					n_z = n_expected_fields - 2;
					if (n_z == 0 && want_z_output) {
						fprintf (stderr, "%s: No data columns, cannot use z flag in -F\n", GMT_program);
						exit (EXIT_FAILURE);
					}
					z_first = FALSE;
				}

				while (! (GMT_io.status & (GMT_IO_SEGMENT_HEADER | GMT_IO_EOF))) {	/* Keep going until FALSE or = 2 segment header */

					n_total_read ++;
					if (GMT_io.status & GMT_IO_MISMATCH) {
						fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_total_read);
						continue;
					}

					xx = in[GMT_X];
					yy = in[GMT_Y];

					if (Ctrl->N.active) {
						x[0] = xx - Ctrl->C.x;
						x[1] = yy - Ctrl->C.y;
						matrix_2v (e,x,xt);
					}
					else {
						oblique_transform(yy, xx, &xt[1], &xt[0], pole, center);
					}

					skip = ((Ctrl->L.active && (xt[0] < Ctrl->L.min || xt[0] > Ctrl->L.max)) || (Ctrl->W.active && (xt[1] < Ctrl->W.min || xt[1] > Ctrl->W.max)));

					if (skip) {
						n_fields = GMT_input (fp, &n_expected_fields, &in);
						continue;
					}

					p_data[n_used].a[0] = xx;
					p_data[n_used].a[1] = yy;
					p_data[n_used].a[2] = xt[0];
					p_data[n_used].a[3] = xt[1];
					p_data[n_used].t = NULL;	/* Initialize since that is not done by realloc */
					p_data[n_used].z = NULL;	/* Initialize since that is not done by realloc */
					if (n_z) {	/* Copy over z column(s) */
#ifdef DEBUG
						GMT_memtrack_off (GMT_mem_keeper);	/* Because it gives way too many pointers to search through */
#endif
						if (pure_ascii) {	/* Must store all text beyond x,y columns */
							p_data[n_used].t = (char *) GMT_memory (VNULL, strlen (GMT_io.current_record), sizeof (char), GMT_program);
							copy_text_from_col3 (GMT_io.current_record, p_data[n_used].t);
						}
						else {
							p_data[n_used].z = (double *) GMT_memory (VNULL, (size_t)n_z, sizeof (double), GMT_program);
							memcpy ((void *)p_data[n_used].z, (void *)&in[GMT_Z], (size_t)(n_z * sizeof(double)));
						}
#ifdef DEBUG
						GMT_memtrack_on (GMT_mem_keeper);
#endif
					}
					n_used++;
					if (n_used == n_alloc) {
						n_alloc <<= 1;
						p_data = (struct PROJECT_DATA *) GMT_memory ((void *)p_data, (size_t)n_alloc, sizeof (struct PROJECT_DATA), GMT_program);
					}

					n_fields = GMT_input (fp, &n_expected_fields, &in);
				}

				if (Ctrl->S.active) qsort ((void *)p_data, (size_t)n_used, sizeof (struct PROJECT_DATA), compare_distances);

/*				Get here when all data are loaded with p,q and p is in increasing order if desired.  */

				if (find_new_point) {	/* We need to find r,s  */

					if (Ctrl->N.active) {
						sincosd (theta, &sin_theta, &cos_theta);
						for (i = 0; i < n_used; i++) {
							p_data[i].a[4] = Ctrl->C.x + p_data[i].a[2] * cos_theta;
							p_data[i].a[5] = Ctrl->C.y + p_data[i].a[2] * sin_theta;
							while (greenwich && p_data[i].a[4] < 0.0) p_data[i].a[4] += 360.0;
							while (dateline && p_data[i].a[4] > 180.0) p_data[i].a[4] -= 360.0;
						}
					}
					else {
						GMT_geo_to_cart(Ctrl->C.y, Ctrl->C.x, x, TRUE);
						for (i = 0; i < n_used; i++) {
							make_euler_matrix(pole, e, p_data[i].a[2]);
							matrix_3v(e,x,xt);
							GMT_cart_to_geo(&(p_data[i].a[5]), &(p_data[i].a[4]), xt, TRUE);
							while (greenwich && p_data[i].a[4] < 0.0) p_data[i].a[4] += 360.0;
							while (dateline && p_data[i].a[4] > 180.0) p_data[i].a[4] -= 360.0;
						}
					}
				}

				/* At this stage, all values are still in degrees.  */

				if (Ctrl->Q.active) {
					for (i = 0; i < n_used; i++) {
						p_data[i].a[2] *= project_info.DIST_KM_PR_DEG;
						p_data[i].a[3] *= project_info.DIST_KM_PR_DEG;
					}
				}

				/* Now output  */

				if (!GMT_io.binary[GMT_OUT]) {	/* First do header */
					if (first && GMT_io.io_header[GMT_OUT]) {
						for (j = 0; j < n_outputs; j++) {
							if (output_choice[j] == -1)
								GMT_fputs (heading[6], GMT_stdout);
							else
								GMT_fputs (heading[output_choice[j]], GMT_stdout);
							(j == (n_outputs - 1)) ? GMT_fputs ("\n", GMT_stdout) : GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
						}
						first = FALSE;
					}
				}

				n_items = n_outputs + ((want_z_output && n_z) ? n_z - 1 : 0);
				if (!out) out = (double *) GMT_memory (VNULL, (size_t)n_items, sizeof (double), GMT_program);

				/* Special case for pure ascii since we may pass text */

				if (n_z && pure_ascii) {
					for (i = 0; i < n_used; i++) {
						for (j = 0; j < n_outputs; j++) {
							if (output_choice[j] == -1)	/* Output all z columns as one string */
								GMT_fputs (p_data[i].t, GMT_stdout);
							else {
								sprintf (buffer, gmtdefs.d_format, p_data[i].a[output_choice[j]]);
								GMT_fputs (buffer, GMT_stdout);
							}
							(j == (n_outputs - 1)) ? GMT_fputs ("\n", GMT_stdout) : GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
						}
					}
				}
				else {	/* Any other i/o combination */
					for (i = 0; i < n_used; i++) {
						for (j = k = 0; j < n_outputs; j++) {
							if (output_choice[j] == -1) {	/* Copy over all z columns */
								memcpy ((void *)&out[k], (void *)p_data[i].z, (size_t)(n_z * sizeof (double)));
								k += n_z;
							}
							else
								out[k++] = p_data[i].a[output_choice[j]];
						}
						GMT_output (GMT_stdout, n_items, out);
					}
				}

				n_total_used += n_used;
				n_used = 0;

			}

			if (fp != GMT_stdin) GMT_fclose(fp);
		}
	}

	if (gmtdefs.verbose) fprintf(stderr, "%s: %ld read, %ld used\n", GMT_program, n_total_read, n_total_used);

#ifdef DEBUG
	if (n_z) GMT_memtrack_off (GMT_mem_keeper);
#endif
	for (i = 0; i < n_alloc; i++) {
		if (p_data[i].t) GMT_free ((void *)p_data[i].t);
		if (p_data[i].z) GMT_free ((void *)p_data[i].z);
	}
#ifdef DEBUG
	if (n_z) GMT_memtrack_on (GMT_mem_keeper);
#endif
	GMT_free ((void *)p_data);
	GMT_free ((void *)out);

	Free_project_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

int compare_distances (const void *point_1, const void *point_2)
{
	double	d_1, d_2;

	d_1 = ((struct PROJECT_DATA *)point_1)->a[2];
	d_2 = ((struct PROJECT_DATA *)point_2)->a[2];

	if (d_1 < d_2)
		return (-1);
	if (d_1 > d_2)
		return (1);
	else
		return (0);
}

double oblique_setup (double plat, double plon, double *p, double *clat, double *clon, double *c, GMT_LONG c_given, GMT_LONG generate)
{
	/* routine sets up a unit 3-vector p, the pole of an 
	   oblique projection, given plat, plon, the position 
	   of this pole in the usual coordinate frame.
	   c_given = TRUE means that clat, clon are to be used
	   as the usual coordinates of a point through which the
	   user wants the central meridian of the oblique
	   projection to go.  If such a point is not given, then
	   the central meridian will go through p and the usual
	   N pole.  In either case, a unit 3-vector c is created
	   which is the directed normal to the plane of the central
	   meridian, pointing in the positive normal (east) sense.
	   Latitudes and longitudes are in degrees. */

	double	s[3];  /* s points to the south pole  */
	double	x[3];  /* tmp vector  */
	double cp, sin_lat_to_pole;

	s[0] = s[1] = 0.0;
	s[2] = -1.0;

	GMT_geo_to_cart(plat, plon, p, TRUE);

	if (c_given) {	/* s points to user's clat, clon  */
		GMT_geo_to_cart(*clat, *clon, s, TRUE);
	}
	GMT_cross3v(p, s, x);
	GMT_normalize3v(x);
	GMT_cross3v(x, p, c);
	GMT_normalize3v(c);
	cp = GMT_dot3v (p, c);
	if (!generate) memcpy ((void *)c, (void *)x, 3*sizeof(double));
	if (!c_given) GMT_cart_to_geo(clat, clon, c, TRUE);	/* return the possibly adjusted center  */
	sin_lat_to_pole = d_sqrt (1.0 - cp * cp);
	return (sin_lat_to_pole);
}

#if 0
double oblique_setup (double plat, double plon, double *p, double *clat, double *clon, double *c, GMT_LONG c_given)
{
	/* routine sets up a unit 3-vector p, the pole of an 
	   oblique projection, given plat, plon, the position 
	   of this pole in the usual coordinate frame.
	   c_given = TRUE means that clat, clon are to be used
	   as the usual coordinates of a point through which the
	   user wants the central meridian of the oblique
	   projection to go.  If such a point is not given, then
	   the central meridian will go through p and the usual
	   N pole.  In either case, a unit 3-vector c is created
	   which is the directed normal to the plane of the central
	   meridian, pointing in the positive normal (east) sense.
	   Latitudes and longitudes are in degrees. */

	double	s[3];  /* s points to the south pole  */
	double cp, sin_lat_to_pole;

	s[0] = s[1] = 0.0;
	s[2] = -1.0;

	GMT_geo_to_cart(plat, plon, p, TRUE);

	if (c_given) {	/* s points to user's clat, clon  */
		GMT_geo_to_cart(*clat, *clon, s, TRUE);
	}
	GMT_cross3v(p, s, c);
	GMT_normalize3v(c);
	cp = GMT_dot3v (p, s);
	sin_lat_to_pole = d_sqrt (1.0 - cp * cp);
	return (sin_lat_to_pole);
}
#endif

void oblique_transform (double xlat, double xlon, double *x_t_lat, double *x_t_lon, double *p, double *c)
{
	/* routine takes the point x at conventional (xlat, xlon) and
	   computes the transformed coordinates (x_t_lat, x_t_lon) in
	   an oblique reference frame specified by the unit 3-vectors
	   p (the pole) and c (the directed normal to the oblique
	   central meridian).  p and c have been computed earlier by
	   the routine oblique_setup().
	   Latitudes and longitudes are in degrees. */

	double	x[3], p_cross_x[3], temp1, temp2;

	GMT_geo_to_cart(xlat, xlon, x, TRUE);

	temp1 = GMT_dot3v(x,p);
	*x_t_lat = d_asind(temp1);

	GMT_cross3v(p,x,p_cross_x);
	GMT_normalize3v(p_cross_x);

	temp1 = GMT_dot3v(p_cross_x, c);
	temp2 = GMT_dot3v(x, c);
	*x_t_lon = copysign(d_acosd(temp1), temp2);
}

void make_euler_matrix (double *p, double *e, double theta)
{
	/* Routine to fill an euler matrix e with the elements
	   needed to rotate a 3-vector about the pole p through
	   an angle theta (in degrees).  p is a unit 3-vector.
	   Latitudes and longitudes are in degrees. */

	double	cos_theta, sin_theta, one_minus_cos_theta;
	double	pxsin, pysin, pzsin, temp;

	sincosd (theta, &sin_theta, &cos_theta);
	one_minus_cos_theta = 1.0 - cos_theta;

	pxsin = p[0] * sin_theta;
	pysin = p[1] * sin_theta;
	pzsin = p[2] * sin_theta;

	temp = p[0] * one_minus_cos_theta;
	e[0] = temp * p[0] + cos_theta;
	e[1] = temp * p[1] - pzsin;
	e[2] = temp * p[2] + pysin;

	temp = p[1] * one_minus_cos_theta;
	e[3] = temp * p[0] + pzsin;
	e[4] = temp * p[1] + cos_theta;
	e[5] = temp * p[2] - pxsin;

	temp = p[2] * one_minus_cos_theta;
	e[6] = temp * p[0] - pysin;
	e[7] = temp * p[1] + pxsin;
	e[8] = temp * p[2] + cos_theta;
}

void matrix_3v (double *a, double *x, double *b)
{
	/* routine to find b, where Ax = b, A is a 3 by 3 square matrix,
	   and x and b are 3-vectors.  A is stored row wise, that is:
	   
	   A = { a11, a12, a13, a21, a22, a23, a31, a32, a33 }  */

	b[0] = x[0]*a[0] + x[1]*a[1] + x[2]*a[2];
	b[1] = x[0]*a[3] + x[1]*a[4] + x[2]*a[5];
	b[2] = x[0]*a[6] + x[1]*a[7] + x[2]*a[8];
}

void matrix_2v (double *a, double *x, double *b)
{
	/* routine to find b, where Ax = b, A is a 2 by 2 square matrix,
	   and x and b are 2-vectors.  A is stored row wise, that is:
	   
	   A = { a11, a12, a21, a22 }  */

	b[0] = x[0]*a[0] + x[1]*a[1];
	b[1] = x[0]*a[2] + x[1]*a[3];
}

void sphere_project_setup (double alat, double alon, double *a, double blat, double blon, double *b, double azim, double *p, double *c, GMT_LONG two_pts)
{
	/* routine to initialize a pole vector, p, and a central meridian 
	   normal vector, c, for use in projecting points onto a great circle.
	   
	   The great circle is specified in either one of two ways:
	   if (two_pts), then the user has given two points, a and b,
	   which specify the great circle (directed from a to b);
	   if !(two_pts), then the user has given one point, a, and an azimuth,
	   azim, clockwise from north, which defines the projection.

	   The strategy is to use the oblique_transform operations above,
	   in such a way that the great circle of the projection is the
	   equator of an oblique transform, and the central meridian goes
	   through a.  Then the transformed longitude gives the distance
	   along the projection circle, and the transformed latitude gives
	   the distance normal to the projection circle.

	   If (two_pts), then p = normalized(a X b).  If not, we temporarily
	   create p_temp = normalized(a X n), where n is the north pole.
	   p_temp is then rotated about a through the angle azim to give p.
	   After p is found, then c = normalized(p X a).

	   Latitudes and longitudes are in degrees.
	*/

	double	e[9];	/* Euler rotation matrix, if needed  */

	/* First find p vector  */

	if (two_pts) {
		GMT_geo_to_cart(alat, alon, a, TRUE);
		GMT_geo_to_cart(blat, blon, b, TRUE);
		GMT_cross3v(a, b, p);
		GMT_normalize3v(p);
	}
	else {
		GMT_geo_to_cart(alat, alon, a, TRUE);
		b[0] = b[1] = 0.0;	/* set b to north pole  */
		b[2] = 1.0;
		GMT_cross3v(a, b, c);	/* use c for p_temp  */
		GMT_normalize3v(c);
		make_euler_matrix(a, e, -azim);
		matrix_3v(e, c, p);	/* c (p_temp) rotates to p  */
	}

	/* Now set c vector  */

	GMT_cross3v(p, a, c);
	GMT_normalize3v(c);
}

void flat_project_setup (double alat, double alon, double blat, double blon, double plat, double plon, double *azim, double *e, GMT_LONG two_pts, GMT_LONG pole_set)
{
	/* Sets up stuff for rotation of cartesian 2-vectors, analogous
	   to the spherical three vector stuff above.
	   Output is the negative cartesian azimuth in degrees.
	   Latitudes and longitudes are in degrees. */

	if (two_pts)
		*azim = 90.0 - d_atan2d(blat - alat, blon - alon);
	else if (pole_set)
		*azim = 180.0 - d_atan2d(plat - alat, plon - alon);

	*azim = -(*azim);
	e[0] = e[3] = cosd(*azim);
	e[1] = sind(*azim);
	e[2] = -e[1];
}

void copy_text_from_col3 (char *line, char *z_cols)
{	/* returns the input line starting at the 3rd column */

	GMT_LONG i;

	/* First replace any commas with spaces */

	for (i = 0; line[i]; i++) if (line[i] == ',') line[i] = ' ';

	sscanf (line, "%*s %*s %[^\n]", z_cols);
}

void *New_project_Ctrl () {	/* Allocate and initialize a new control structure */
	struct PROJECT_CTRL *C;
	
	C = (struct PROJECT_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PROJECT_CTRL), "New_project_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	
	return ((void *)C);
}

void Free_project_Ctrl (struct PROJECT_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);	
}
