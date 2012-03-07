/*--------------------------------------------------------------------
 *	$Id: grdproject.c,v 1.62 2011/07/08 21:27:06 guru Exp $
 *
 *	Copyright (c) 1991-2011 by P. Wessel and W. H. F. Smith
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
 * grdproject reads a geographical grid file and evaluates the grid at new grid positions
 * specified by the map projection and new dx/dy values using a weighted average of all
 * points within the search radius. Optionally, grdproject may perform the inverse
 * transformation, going from rectangular coordinates to geographical.
 *
 * Author:	Paul Wessel
 * Date:	15-JUL-2000
 * Ver:		4
 *
 */

#include "gmt.h"

struct GRDPROJECT_CTRL {
	struct A {	/* -A[k|m|n|i|c|p] */
		GMT_LONG active;
		char unit;
	} A;
	struct C {	/* -C[<dx/dy>] */
		GMT_LONG active;
		double easting, northing;
	} C;
	struct D {	/* -Ddx[/dy] */
		GMT_LONG active;
		double xinc, yinc;
	} D;
	struct E {	/* -E<dpi> */
		GMT_LONG active;
		GMT_LONG dpi;
	} E;
	struct F {	/* -F */
		GMT_LONG active;
	} F;
	struct G {	/* -G */
		GMT_LONG active;
		char *file;
	} G;
	struct I {	/* -I */
		GMT_LONG active;
	} I;
	struct M {	/* -Mc|i|m */
		GMT_LONG active;
		char unit;
	} M;
	struct N {	/* -N<nx/ny> */
		GMT_LONG active;
		GMT_LONG nx, ny;
	} N;
	struct S {	/* -S[-]b|c|l|n[/<threshold>] */
		GMT_LONG active;
		GMT_LONG antialias;
		GMT_LONG interpolant;
		double threshold;
	} S;
};

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE, set_n = FALSE, shift_xy = FALSE, offset;

	char *infile = NULL, format[BUFSIZ], unit_name[GRD_UNIT_LEN], scale_unit_name[GRD_UNIT_LEN];

	GMT_LONG i, j, unit = 0;
	
	GMT_LONG nm;
	
	float *geo = NULL, *rect = NULL;
	
	double w, e, s, n;
	double xmin, xmax, ymin, ymax, inch_to_unit, unit_to_inch, fwd_scale, inv_scale;

	struct GRD_HEADER g_head, r_head;
	struct GMT_EDGEINFO edgeinfo;
	struct GRDPROJECT_CTRL *Ctrl = NULL;

	void *New_grdproject_Ctrl (), Free_grdproject_Ctrl (struct GRDPROJECT_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDPROJECT_CTRL *)New_grdproject_Ctrl ();	/* Allocate and initialize a new control structure */
	
	infile = CNULL;
	w = e = s = n = 0.0;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'J':
				case 'R':
				case 'V':
				case '\0':
					error += GMT_parse_common_options (argv[i], &w, &e, &s, &n);
					break;

				/* Supplemental parameters */

				case 'C':
					Ctrl->C.active = TRUE;
					if (argv[i][2]) {	/* Also gave shifts */
	 					n = sscanf (&argv[i][2], "%lf/%lf", &Ctrl->C.easting, &Ctrl->C.northing);
						if (n != 2) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR.  Expected -C[<false_easting>/<false_northing>]\n", GMT_program);
							error++;
						}
					}
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					if (GMT_getinc (&argv[i][2], &Ctrl->D.xinc, &Ctrl->D.yinc)) {
						GMT_inc_syntax ('D', 1);
						error = TRUE;
					}
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					Ctrl->E.dpi = atoi (&argv[i][2]);
					break;
				case 'A':
					Ctrl->A.active = TRUE;
					Ctrl->A.unit = argv[i][2];
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					break;
				case 'G':
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					break;
				case 'M':	/* Directly specify units */
					Ctrl->M.active = TRUE;
					Ctrl->M.unit = argv[i][2];
					break;
				case 'N':
					sscanf (&argv[i][2], "%" GMT_LL "d/%" GMT_LL "d", &Ctrl->N.nx, &Ctrl->N.ny);
					if (Ctrl->N.ny == 0) Ctrl->N.ny = Ctrl->N.nx;
					Ctrl->N.active = TRUE;
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					for (j = 2; j < 5 && argv[i][j]; j++) {
						switch (argv[i][j]) {
							case '-':
								Ctrl->S.antialias = FALSE; break;
							case 'n':
								Ctrl->S.interpolant = BCR_NEARNEIGHBOR; break;
							case 'l':
								Ctrl->S.interpolant = BCR_BILINEAR; break;
							case 'b':
								Ctrl->S.interpolant = BCR_BSPLINE; break;
							case 'c':
								Ctrl->S.interpolant = BCR_BICUBIC; break;
							case '/':
								Ctrl->S.threshold = atof (&argv[i][j+1]);
								j = 5; break;
							default:
								fprintf (stderr, "%s: Warning: The -S option has changed meaning. Use -S[-]b|c|l|n[/threshold] to specify interpolation mode.\n", GMT_program);
								j = 5; break;
						}
					}
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else 
			infile = argv[i];
	}

	if ((Ctrl->D.active + Ctrl->E.active + Ctrl->N.active) == 0) Ctrl->N.active = set_n = TRUE;

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "grdproject %s - Project geographical grid to/from rectangular grid\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdproject <in_grdfile> -G<out_grdfile> %s\n", GMT_J_OPT);
		fprintf (stderr, "\t[-A[k|m|n|i|c|p]] [-C[<dx/dy>]] [-D%s] [-E<dpi>] [-F]\n", GMT_inc_OPT);
		fprintf (stderr, "\t[-I] [-Mc|i|m] [-N<nx/ny>] [%s]\n", GMT_Rgeo_OPT);
		fprintf (stderr, "\t[-S[-]b|c|l|n[/<threshold>]] [-V]\n\n");

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<in_grdfile> is data set to be transformed\n");
		fprintf (stderr, "\t-G name of output grid\n");
		GMT_explain_option ('J');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A force projected values to be in actual meters [Default uses the given map scale]\n");
		fprintf (stderr, "\t   Specify another unit by appending k (km), m (miles), n (nautical miles), i (inch), c (cm), or p (points)\n");
		fprintf (stderr, "\t-C coordinates relative to projection center [Default is relative to lower left corner]\n");
		fprintf (stderr, "\t   Optionally append dx/dy to add (or subtract if -I) (i.e., false easting & northing) [0/0]\n");
		GMT_inc_syntax ('D', 0);
		fprintf (stderr, "\t-E sets dpi for output grid\n");
		fprintf (stderr, "\t-F toggle between pixel and grid registration  [Default is same as input]\n");
		fprintf (stderr, "\t-I Inverse transformation from rectangular to geographical\n");
		fprintf (stderr, "\t-M Temporarily reset MEASURE_UNIT to be c (cm), i (inch), m (meter), or p (point)\n");
		fprintf (stderr, "\t   Cannot be used if -A is set.\n");
		fprintf (stderr, "\t-N sets the number of nodes for the new grid\n");
		fprintf (stderr, "\t   Only one of -D, -E, and -N can be specified!\n");
		fprintf (stderr, "\t   If none are specified, nx,ny of the input file is used\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t-S Determines the interpolation mode (b = B-spline, c = bicubic, l = bilinear,\n");
		fprintf (stderr, "\t   n = nearest-neighbor) [Default: bicubic]\n");
		fprintf (stderr, "\t   Optionally, prepend - to switch off antialiasing [Default: on]\n");
		fprintf (stderr, "\t   Append /<threshold> to change the minimum weight in vicinity of NaNs. A threshold of\n");
		fprintf (stderr, "\t   1.0 requires all nodes involved in interpolation to be non-NaN; 0.5 will interpolate\n");
		fprintf (stderr, "\t   about half way from a non-NaN to a NaN node [Default: 0.5]\n");
		GMT_explain_option ('V');

		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&Ctrl->D.xinc, &Ctrl->D.yinc, &Ctrl->F.active, &Ctrl->D.active);

	if (!infile) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input file\n", GMT_program);
		error++;
	}
	if (!Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option:  Must specify output file\n", GMT_program);
		error++;
	}
	/*if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}*/
	if ((Ctrl->M.active + Ctrl->A.active) == 2) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Can specify only one of -A and -M\n", GMT_program);
		error++;
	}
	if ((Ctrl->D.active + Ctrl->E.active + Ctrl->N.active) != 1) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify only one of -D, -E, or -N\n", GMT_program);
		error++;
	}
	if (Ctrl->D.active && (Ctrl->D.xinc <= 0.0 || Ctrl->D.yinc < 0.0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -D option.  Must specify positive increment(s)\n", GMT_program);
		error++;
	}
	if (Ctrl->N.active && !set_n && (Ctrl->N.nx <= 0 || Ctrl->N.ny <= 0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -N option.  Must specify positive integers\n", GMT_program);
		error++;
	}
	if (Ctrl->E.active && Ctrl->E.dpi <= 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option.  Must specify positive dpi\n", GMT_program);
		error++;
	}
	if (Ctrl->S.active && (Ctrl->S.threshold < 0.0 || Ctrl->S.threshold > 1.0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  threshold must be in [0,1] range\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (Ctrl->M.active) GMT_err_fail (GMT_set_measure_unit (Ctrl->M.unit), "-M");
	shift_xy = !(Ctrl->C.easting == 0.0 && Ctrl->C.northing == 0.0);
	
	unit = GMT_check_scalingopt ('A', Ctrl->A.unit, scale_unit_name);
	GMT_init_scales (unit, &fwd_scale, &inv_scale, &inch_to_unit, &unit_to_inch, unit_name);

	if (Ctrl->I.active) {	/* Must flip the column types since in is Cartesian and out is geographic */
		GMT_io.out_col_type[0] = GMT_IS_LON;	GMT_io.out_col_type[1] = GMT_IS_LAT;	/* Inverse projection expects x,y and gives lon, lat */
		GMT_io.in_col_type[0] = GMT_io.in_col_type[1] = GMT_IS_FLOAT;
		project_info.degree[0] = project_info.degree[1] = FALSE;
	}

	if (!project_info.region_supplied) {
		double ww, ee, ss, nn;
		char opt_R[BUFSIZ];
		struct GRD_HEADER head;
		if (project_info.utm_hemisphere == 0) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -J option: When -R is not provided you have to specify the hemisphere\n", GMT_program);
			exit (EXIT_FAILURE);
		}
		GMT_err_fail (GMT_read_grd_info (infile, &head), infile);
		w = head.x_min;		e = head.x_max;
		s = head.y_min;		n = head.y_max;
		if (!Ctrl->I.active) {
			sprintf (opt_R, "-R%.12f/%.12f/%.12f/%.12f", w, e, s, n);
			GMT_parse_R_option (opt_R, &ww, &ee, &ss, &nn);
			GMT_err_fail (GMT_map_setup (w, e, s, n), "");
		}
		else {			/* Do inverse transformation */
			double x_c, y_c, lon_t, lat_t;
			/* Obtain a first crude estimation of the good -R */
			x_c = (w + e) / 2.; 		/* mid point of projected coords */
			y_c = (s + n) / 2.;
			if (project_info.projection == GMT_UTM && !project_info.north_pole && y_c > 0) y_c *= -1;
			if (y_c > 0)
				GMT_parse_R_option ("-R-180/180/0/80", &ww, &ee, &ss, &nn);
			else
				GMT_parse_R_option ("-R-180/180/-80/0", &ww, &ee, &ss, &nn);
			if (project_info.projection == GMT_UTM && !project_info.north_pole && y_c < 0) y_c *= -1;	/* Undo the *-1 (only for the UTM case) */
			if (shift_xy) {
				x_c -= Ctrl->C.easting;
				y_c -= Ctrl->C.northing;
			}
			/* Convert from 1:1 scale */ 
			if (unit) {
				x_c *= fwd_scale;
				y_c *= fwd_scale;
			}

			GMT_err_fail (GMT_map_setup (ww, ee, ss, nn), "");

			x_c *= project_info.x_scale;
			y_c *= project_info.y_scale;

			if (Ctrl->C.active) {	/* Then correct so lower left corner is (0,0) */
				x_c += project_info.x0;
				y_c += project_info.y0;
			}
			GMT_xy_to_geo (&lon_t, &lat_t, x_c, y_c);
			sprintf (opt_R, "-R%.12f/%.12f/%.12f/%.12f", lon_t-1, lon_t+1, lat_t-1, lat_t+1);
			if (gmtdefs.verbose) fprintf (stderr, "First opt_R\t %s\t%g\t%g\n", opt_R, x_c, y_c);
			GMT_parse_R_option (opt_R, &ww, &ee, &ss, &nn);
			project_info.region = 0;	/* We need to reset this to not fall into non-wanted branch deeper down */
			GMT_err_fail (GMT_map_setup (ww, ee, ss, nn), "");

			/* Finally obtain the good limits */
			if (shift_xy) {
				w -= Ctrl->C.easting;	e -= Ctrl->C.easting;
				s -= Ctrl->C.northing;	n -= Ctrl->C.northing;
			}
			if (unit) {
				w *= fwd_scale;		e *= fwd_scale;
				s *= fwd_scale;		n *= fwd_scale;
			}
			w *= project_info.x_scale;	e *= project_info.x_scale;
			s *= project_info.y_scale;	n *= project_info.y_scale;

			if (Ctrl->C.active) {
				w += project_info.x0;	e += project_info.x0;
				s += project_info.y0;	n += project_info.y0;
			}

			GMT_xy_to_geo (&ww, &ss, w, s);		/* SW corner */
			GMT_xy_to_geo (&ee, &nn, e, n);		/* NE corner */
			sprintf (opt_R, "-R%.12f/%.12f/%.12f/%.12fr", ww, ss, ee, nn);
			if (gmtdefs.verbose) fprintf (stderr, "Second opt_R\t %s\n", opt_R);
			GMT_parse_common_options (opt_R, &ww, &ee, &ss, &nn);
			w = ww;		e = ee;
			s = ss;		n = nn;
		}
	}

	GMT_err_fail (GMT_map_setup (w, e, s, n), "");

	xmin = (Ctrl->C.active) ? project_info.xmin - project_info.x0 : project_info.xmin;
	xmax = (Ctrl->C.active) ? project_info.xmax - project_info.x0 : project_info.xmax;
	ymin = (Ctrl->C.active) ? project_info.ymin - project_info.y0 : project_info.ymin;
	ymax = (Ctrl->C.active) ? project_info.ymax - project_info.y0 : project_info.ymax;
	if (Ctrl->A.active) {	/* Convert to chosen units */
		strncpy (unit_name, scale_unit_name, (size_t)GRD_UNIT_LEN);
		xmin /= project_info.x_scale;
		xmax /= project_info.x_scale;
		ymin /= project_info.y_scale;
		ymax /= project_info.y_scale;
		if (unit) {	/* Change the 1:1 unit used */
			xmin *= fwd_scale;
			xmax *= fwd_scale;
			ymin *= fwd_scale;
			ymax *= fwd_scale;
		}
	}
	else {	/* Convert inches to chosen MEASURE */
		xmin *= inch_to_unit;
		xmax *= inch_to_unit;
		ymin *= inch_to_unit;
		ymax *= inch_to_unit;
	}
	if (shift_xy) {
		xmin += Ctrl->C.easting;
		xmax += Ctrl->C.easting;
		ymin += Ctrl->C.northing;
		ymax += Ctrl->C.northing;
	}

	GMT_grd_init (&r_head, argc, argv, FALSE);
	GMT_grd_init (&g_head, argc, argv, FALSE);

	sprintf (format, "(%s/%s/%s/%s)", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);

	if (Ctrl->I.active) {	/* Transforming from rectangular projection to geographical */

		/* if (!project_info.region) d_swap (s, e); */  /* Got w/s/e/n, make into w/e/s/n */

		g_head.x_min = w;	g_head.x_max = e;	g_head.y_min = s;	g_head.y_max = n;

		GMT_err_fail (GMT_read_grd_info (infile, &r_head), infile);

		GMT_boundcond_init (&edgeinfo);
		nm = GMT_get_nm (4 + r_head.nx, 4 + r_head.ny);
		GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 2;

		rect = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
		GMT_err_fail (GMT_read_grd (infile, &r_head, rect, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), infile);
		offset = r_head.node_offset;		/* Same as input */
		if (Ctrl->F.active) offset = !offset;	/* Toggle */
		if (set_n) {
			Ctrl->N.nx = r_head.nx;
			Ctrl->N.ny = r_head.ny;
		}
		GMT_err_fail (GMT_grdproject_init (&g_head, Ctrl->D.xinc, Ctrl->D.yinc, Ctrl->N.nx, Ctrl->N.ny, Ctrl->E.dpi, offset), Ctrl->G.file);
		nm = GMT_get_nm (g_head.nx, g_head.ny);
		geo = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
		if (gmtdefs.verbose) {
			fprintf (stderr, "%s: Transform ", GMT_program);
			fprintf (stderr, format, g_head.x_min, g_head.x_max, g_head.y_min, g_head.y_max);
			fprintf (stderr, " <-- ");
			fprintf (stderr, format, xmin, xmax, ymin, ymax);
			fprintf (stderr, " [%s]\n", unit_name);
		}

		/* Modify input rect header if -A, -C, -M have been set */

		if (shift_xy) {
			r_head.x_min -= Ctrl->C.easting;
			r_head.x_max -= Ctrl->C.easting;
			r_head.y_min -= Ctrl->C.northing;
			r_head.y_max -= Ctrl->C.northing;

		}
		if (Ctrl->A.active) {	/* Convert from 1:1 scale */
			if (unit) {	/* Undo the 1:1 unit used */
				r_head.x_min *= inv_scale;
				r_head.x_max *= inv_scale;
				r_head.y_min *= inv_scale;
				r_head.y_max *= inv_scale;
			}
			r_head.x_min *= project_info.x_scale;
			r_head.x_max *= project_info.x_scale;
			r_head.y_min *= project_info.y_scale;
			r_head.y_max *= project_info.y_scale;
		}
		else if (gmtdefs.measure_unit != GMT_INCH) {	/* Convert from inch to whatever */
			r_head.x_min *= unit_to_inch;
			r_head.x_max *= unit_to_inch;
			r_head.y_min *= unit_to_inch;
			r_head.y_max *= unit_to_inch;
		}
		if (Ctrl->C.active) {	/* Then correct so lower left corner is (0,0) */
			r_head.x_min += project_info.x0;
			r_head.x_max += project_info.x0;
			r_head.y_min += project_info.y0;
			r_head.y_max += project_info.y0;
		}
		r_head.x_inc = GMT_get_inc (r_head.x_min, r_head.x_max, r_head.nx, r_head.node_offset);
		r_head.y_inc = GMT_get_inc (r_head.y_min, r_head.y_max, r_head.ny, r_head.node_offset);

		GMT_grd_project (rect, &r_head, geo, &g_head, &edgeinfo, Ctrl->S.antialias, Ctrl->S.interpolant, Ctrl->S.threshold, TRUE);

		GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 0;
		sprintf (g_head.x_units, "longitude [degrees_east]");
		sprintf (g_head.y_units, "latitude [degrees_north]");
		GMT_err_fail (GMT_write_grd (Ctrl->G.file, &g_head, geo, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);
	}
	else {	/* Forward projection from geographical to rectangular grid */

		GMT_err_fail (GMT_read_grd_info (infile, &g_head), infile);

		GMT_boundcond_init (&edgeinfo);
		nm = GMT_get_nm (4 + g_head.nx, 4 + g_head.ny);
		GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 2;

		geo = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
		GMT_err_fail (GMT_read_grd (infile, &g_head, geo, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), infile);

		r_head.x_min = project_info.xmin;	r_head.x_max = project_info.xmax;
		r_head.y_min = project_info.ymin;	r_head.y_max = project_info.ymax;
		if (Ctrl->A.active) {	/* Convert from 1:1 scale */
			if (unit) {	/* Undo the 1:1 unit used */
				Ctrl->D.xinc *= inv_scale;
				Ctrl->D.yinc *= inv_scale;
			}
			Ctrl->D.xinc *= project_info.x_scale;
			Ctrl->D.yinc *= project_info.y_scale;
		}
		else if (gmtdefs.measure_unit != GMT_INCH) {	/* Convert from inch to whatever */
			Ctrl->D.xinc *= unit_to_inch;
			Ctrl->D.yinc *= unit_to_inch;
		}
		if (set_n) {
			Ctrl->N.nx = g_head.nx;
			Ctrl->N.ny = g_head.ny;
		}

		if (gmtdefs.verbose) {
			fprintf (stderr, "%s: Transform ", GMT_program);
			fprintf (stderr, format, g_head.x_min, g_head.x_max, g_head.y_min, g_head.y_max);
			fprintf (stderr, " --> ");
			fprintf (stderr, format, xmin, xmax, ymin, ymax);
			fprintf (stderr, " [%s]\n", unit_name);
		}

		offset = g_head.node_offset;		/* Same as input */
		if (Ctrl->F.active) offset = !offset;	/* Toggle */

		GMT_err_fail (GMT_grdproject_init (&r_head, Ctrl->D.xinc, Ctrl->D.yinc, Ctrl->N.nx, Ctrl->N.ny, Ctrl->E.dpi, offset), Ctrl->G.file);
		nm = GMT_get_nm (r_head.nx, r_head.ny);
		rect = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
		GMT_grd_project (geo, &g_head, rect, &r_head, &edgeinfo, Ctrl->S.antialias, Ctrl->S.interpolant, Ctrl->S.threshold, FALSE);

		/* Modify output rect header if -A, -C, -M have been set */

		if (Ctrl->C.active) {	/* Change origin from lower left to projection center */
			r_head.x_min -= project_info.x0;
			r_head.x_max -= project_info.x0;
			r_head.y_min -= project_info.y0;
			r_head.y_max -= project_info.y0;
		}
		if (Ctrl->A.active) {	/* Convert to 1:1 scale */
			r_head.x_min /= project_info.x_scale;
			r_head.x_max /= project_info.x_scale;
			r_head.y_min /= project_info.y_scale;
			r_head.y_max /= project_info.y_scale;
			if (unit) {	/* Change the 1:1 unit used */
				r_head.x_min *= fwd_scale;
				r_head.x_max *= fwd_scale;
				r_head.y_min *= fwd_scale;
				r_head.y_max *= fwd_scale;
			}
		}
		else if (gmtdefs.measure_unit != GMT_INCH) {	/* Convert from inch to whatever */
			r_head.x_min /= unit_to_inch;
			r_head.x_max /= unit_to_inch;
			r_head.y_min /= unit_to_inch;
			r_head.y_max /= unit_to_inch;
		}
		if (shift_xy) {
			r_head.x_min += Ctrl->C.easting;
			r_head.x_max += Ctrl->C.easting;
			r_head.y_min += Ctrl->C.northing;
			r_head.y_max += Ctrl->C.northing;

		}
		r_head.x_inc = GMT_get_inc (r_head.x_min, r_head.x_max, r_head.nx, r_head.node_offset);
		r_head.y_inc = GMT_get_inc (r_head.y_min, r_head.y_max, r_head.ny, r_head.node_offset);

		/* rect xy values are here in GMT projected units chosen by user */

		GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 0;
		strcpy (r_head.x_units, unit_name);
		strcpy (r_head.y_units, unit_name);
		GMT_err_fail (GMT_write_grd (Ctrl->G.file, &r_head, rect, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);
	}

	GMT_free ((void *)geo);
	GMT_free ((void *)rect);

	Free_grdproject_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_grdproject_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDPROJECT_CTRL *C;
	
	C = (struct GRDPROJECT_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDPROJECT_CTRL), "New_grdproject_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->S.antialias = TRUE; C->S.interpolant = BCR_BICUBIC; C->S.threshold = 0.5;
		
	return ((void *)C);
}

void Free_grdproject_Ctrl (struct GRDPROJECT_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}
