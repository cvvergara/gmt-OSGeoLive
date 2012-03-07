/*--------------------------------------------------------------------
 *	$Id: grdfilter.c,v 1.91 2011/07/08 22:39:06 guru Exp $
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
 * grdfilter.c  reads a grid file and creates filtered grid file
 *
 * user selects boxcar, gaussian, cosine arch, median, or mode filter
 * user selects distance scaling as appropriate for map work, etc.
 *
 * Author:  W.H.F. Smith
 * Date: 	16 Aug 88
 *
 * Modified:	27 Sep 88 by W.H.F. Smith, to use new median routine.
 *
 * Updated:	PW: 13-Jun-1991 to v2.0
 *		PW: 13-Jul-1992 to actually do v2 i/o.
 *		PW: 15-Jul-1992 to handle arbitrary new -R -I
 *		PW: 03-Jan-1995 to offer mode-filter (from filter1d)
 *		WS: 21-May-1998 to make warnings silent unless -V on.
 *		PW: 03-Jun-1998 upgrade to GMT 3.1
 *		PW: 02-Jun-1999 upgrade to GMT 3.3 + added SINCOS option
 *		PW: 18-Oct-1999 Use sincos directly
 *		PW: 18-JUN-2000 3.3.5
 * Version:	4
 *		PW: 27-MAY-2004 Added extreme values filter options l, L, u, U
 *		PW: 16-NOV-2005 Added spherical filtering support & highpass option
 *		PW: 27-APR-2007 Added -D5 to handle Mercator (img) grids
*/

#define GMT_WITH_NO_PS
#include "gmt.h"

struct GRDFILTER_CTRL {
	struct D {	/* -D<distflag> */
		GMT_LONG active;
		GMT_LONG mode;
	} D;
	struct F {	/* <type>[-]<filter_width>[<mode>] */
		GMT_LONG active;
		GMT_LONG highpass;
		char filter;	/* Character codes for the filter */
		double width;
		GMT_LONG mode;
	} F;
	struct G {	/* -G<file> */
		GMT_LONG active;
		char *file;
	} G;
	struct I {	/* -Idx[/dy] */
		GMT_LONG active;
		double xinc, yinc;
	} I;
	struct N {	/* -Np|i|r */
		GMT_LONG active;
		GMT_LONG mode;	/* 0 is default (i), 1 is replace (r), 2 is preserve (p) */
	} N;
	struct T {	/* -T */
		GMT_LONG active;
	} T;
};

#define IMG2LAT(img) (2.0*atand(exp((img)*D2R))-90.0)
#define LAT2IMG(lat) (R2D*log(tand(0.5*((lat)+90.0))))

#define GRDFILTER_WIDTH		0
#define GRDFILTER_HALF_WIDTH	1
#define GRDFILTER_X_SCALE	2
#define GRDFILTER_Y_SCALE	3
#define GRDFILTER_INV_R_SCALE	4

#define GRDFILTER_N_FILTERS	9

#define NAN_IGNORE	0
#define NAN_REPLACE	1
#define NAN_PRESERVE	2

struct FILTER_INFO {
	GMT_LONG nx;			/* The max number of filter weights in x-direction */
	GMT_LONG ny;			/* The max number of filter weights in y-direction */
	GMT_LONG x_half_width;	/* Number of filter nodes to either side needed at this latitude */
	GMT_LONG y_half_width;	/* Number of filter nodes above/below this point (ny_f/2) */
	GMT_LONG d_flag;
	double dx, dy;		/* Grid spacing in original units */
	double y_min, y_max;	/* Grid limits in y(lat) */
	double *x, *y;		/* Distances in original units along x and y to distance nodes */
	double par[5];		/* [0] is filter width, [1] is 0.5*filter_width, [2] is xscale, [3] is yscale, [4] is 1/r_half for filter */
	double x_off, y_off;	/* Offsets relative to original grid */
	char *visit;		/* BOOLEAN array to keep track of which longitude nodes we have already visited once */
	PFD weight_func, radius_func;
};

int main (int argc, char **argv)
{
	GMT_LONG	nx_out, ny_out, n_in_median, n_nan = 0, visit_check = FALSE;
	GMT_LONG	j_origin, i_out, j_out, half_nx, i_orig = 0, go_on, nm;
	GMT_LONG	i_in, j_in, ii, jj, i, j, ij_in, ij_out, ij_wt, effort_level;
	GMT_LONG	filter_type, one_or_zero = 1, GMT_n_multiples = 0, nx_wrap = 0;
	
	GMT_LONG *i_origin = NULL;


	GMT_LONG error, new_range, fast_way, slow = FALSE, same_grid = FALSE;
	GMT_LONG wrap_case_x = FALSE, wrap_case_y = FALSE;
	GMT_LONG full_360, full_180;

	float	*input = NULL, *output = NULL, *A = NULL;

	double	west_new, east_new, south_new, north_new, merc_range, lat_out, w;
	double	x_scale = 1.0, y_scale = 1.0, x_width, y_width, y, par[5];
	double	x_out, y_out, wt_sum, value, last_median, this_median, xincnew2, yincnew2;
	double	xincold2, yincold2, y_shift = 0.0, x_fix = 0.0, y_fix = 0.0, max_lat;
	double	*weight = NULL, *work_array = VNULL, *x_shift = VNULL;

	char	*fin = CNULL, c, filter_code[GRDFILTER_N_FILTERS] = {'b', 'c', 'g', 'm', 'p', 'l', 'L', 'u', 'U'};
	char	*filter_name[GRDFILTER_N_FILTERS] = {"Boxcar", "Cosine Arch", "Gaussian", "Median", "Mode", "Lower", "Lower+", "Upper", "Upper-"};

	struct	GRD_HEADER h, test_h;
	struct	FILTER_INFO F;
	struct	GRDFILTER_CTRL *Ctrl = NULL;

	void   set_weight_matrix (struct FILTER_INFO *F, double *weight, double output_lat, double par[], double x_off, double y_off);
	void init_area_weights (struct GRD_HEADER *G, GMT_LONG mode, float *A);
	double UnitWeight (double r, double par[]);
	double CosBellWeight (double r, double par[]);
	double GaussianWeight (double r, double par[]);
	double CartRadius (double x0, double y0, double x1, double y1, double par[]);
	double CartScaledRadius (double x0, double y0, double x1, double y1, double par[]);
	double FlatEarthRadius (double x0, double y0, double x1, double y1, double par[]);
	double SphericalRadius (double x0, double y0, double x1, double y1, double par[]);
	void *New_grdfilter_Ctrl (), Free_grdfilter_Ctrl (struct GRDFILTER_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDFILTER_CTRL *)New_grdfilter_Ctrl ();	/* Allocate and initialize a new control structure */
	
	error = new_range = FALSE;
	fin = NULL;
	west_new = east_new = 0.0;
	filter_type = -1;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'R':
				case 'V':
				case 'f':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west_new, &east_new, &south_new, &north_new);
					break;

				case 'D':
					Ctrl->D.active = TRUE;
					Ctrl->D.mode = atoi(&argv[i][2]);
					break;
				case 'F':
					if (strchr ("bcgmpLlUu", argv[i][2])) {	/* OK filter code */
						Ctrl->F.active = TRUE;
						Ctrl->F.filter = argv[i][2];
						Ctrl->F.width = atof (&argv[i][3]);
						if (Ctrl->F.width < 0.0) Ctrl->F.highpass = TRUE;
						Ctrl->F.width = fabs (Ctrl->F.width);
						if (Ctrl->F.filter == 'p') {	/* Check for some futher info in case of mode filtering */
							c = argv[i][strlen(argv[i])-1];
							if (c == '-') Ctrl->F.mode = -1;
							if (c == '+') Ctrl->F.mode = +1;
						}
					}
					else {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -F option.  Correct syntax: -FX<width>, X one of bcgmplLuU\n", GMT_program);
						error++;
					}
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					if (GMT_getinc (&argv[i][2], &Ctrl->I.xinc, &Ctrl->I.yinc)) {
						GMT_inc_syntax ('I', 1);
						error = TRUE;
					}
					break;
				case 'N':
					if (!argv[i][2]) {	/* Pixel registration OBSOLETE but BACKWARD COMPATIBLE for now */
						one_or_zero = 0;
					}
					else {	/* Treatment of NaNs */
						Ctrl->N.active = TRUE;
						switch (argv[i][2]) {
							case 'i':
								Ctrl->N.mode = NAN_IGNORE;	/* Default */
								break;
							case 'r':
								Ctrl->N.mode = NAN_REPLACE;	/* Replace */
								break;
							case 'p':
								Ctrl->N.mode = NAN_PRESERVE;	/* Preserve */
								break;
							default:
								fprintf (stderr, "%s: GMT SYNTAX ERROR -N option.  Correct syntax: -Ni|p|r\n", GMT_program);
								break;
						}
					}
					break;
				case 'T':	/* Toggle registration */
					Ctrl->T.active = TRUE;
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			fin = argv[i];
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "grdfilter %s - Filter a 2-D grid file in the space (or time) domain\n\n", GMT_VERSION);
		fprintf(stderr,"usage: grdfilter input_file -D<distance_flag> -F<type>[-]<filter_width>[<mode>]\n");
		fprintf(stderr,"\t-G<output_file> [%s] [-Ni|p|r] [%s] [-T] [-V] [%s]\n", GMT_I_OPT, GMT_Rgeo_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf(stderr,"\tDistance flag determines how grid (x,y) maps into distance units of filter width as follows:\n");
		fprintf(stderr,"\t   -D0 grid x,y same units as <filter_width>, cartesian Distances.\n");
		fprintf(stderr,"\t   -D1 grid x,y in degrees, <filter_width> in km, cartesian Distances.\n");
		fprintf(stderr,"\t   -D2 grid x,y in degrees, <filter_width> in km, x_scaled by cos(middle y), cartesian Distances.\n");
		fprintf(stderr,"\t   These first three options are faster; they allow weight matrix to be computed only once.\n");
		fprintf(stderr,"\t   Next three options are slower; weights must be recomputed for each scan line.\n");
		fprintf(stderr,"\t   -D3 grid x,y in degrees, <filter_width> in km, x_scale varies as cos(y), cartesian Distances.\n");
		fprintf(stderr,"\t   -D4 grid x,y in degrees, <filter_width> in km, spherical Distances.\n");
		fprintf(stderr,"\t   -D5 grid x,y in Mercator units (-Jm1), <filter_width> in km, spherical Distances.\n");
		fprintf(stderr,"\t-F sets the low-pass filter type and full diameter (6 sigma) filter-width.  Choose between\n");
		fprintf(stderr,"\t   convolution-type filters which differ in how weights are assigned and geospatial\n");
		fprintf(stderr,"\t   filters that seek to return a representative value.\n");
		fprintf(stderr,"\t   Give negative filter width to select highpass filtering [lowpass].\n");
		fprintf(stderr,"\t   Convolution filters:\n");
		fprintf(stderr,"\t     b: Boxcar : a simple averaging of all points inside filter radius.\n");
		fprintf(stderr,"\t     c: Cosine arch : a weighted averaging with cosine arc weights\n");
		fprintf(stderr,"\t     g: Gaussian : weighted averaging with Gaussian weights.\n");
		fprintf(stderr,"\t   Geospatial filters:\n");
		fprintf(stderr,"\t     l: Lower : return minimum of all points.\n");
		fprintf(stderr,"\t     L: Lower+ : return minimum of all +ve points.\n");
		fprintf(stderr,"\t     m: Median : return the median value of all points.\n");
		fprintf(stderr,"\t     p: Maximum likelihood probability estimator : return mode of all points.\n");
		fprintf(stderr,"\t        By default, we return the average if more than one mode is found.\n");
		fprintf(stderr,"\t        Append - or + to the width to instead return the smallest or largest mode.\n");
		fprintf(stderr,"\t     u: Upper : return maximum of all points.\n");
		fprintf(stderr,"\t     U: Upper- : return maximum of all -ve points.\n");
		fprintf(stderr,"\t-G sets output filename for filtered grid\n");
		fprintf(stderr, "\n\tOPTIONS:\n");
		GMT_inc_syntax ('I', 0);
		fprintf(stderr,"\t   The new xinc and yinc should be divisible by the old ones (new lattice is subset of old).\n");
		fprintf(stderr, "\t-N specifies how NaNs in the input grid should be treated.  There are three options:\n");
		fprintf(stderr, "\t   -Ni skips all NaN values and returns a filtered value unless all are NaN [Default]\n");
		fprintf(stderr, "\t   -Np sets filtered output to NaN is any NaNs are found inside filter circle.\n");
		fprintf(stderr, "\t   -Nr sets filtered output to NaN if the corresponding input node was NaN.\n");
		fprintf(stderr, "\t      (only possible if the input and output grids are coregistered).\n");
		fprintf(stderr, "\t-T Toggles between grid and pixel registration for output grid [Default is same as input registration]\n");
		fprintf(stderr, "\t-R for new Range of output grid; enter <WESN> (xmin, xmax, ymin, ymax) separated by slashes.\n");
		GMT_explain_option ('V');
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&Ctrl->I.xinc, &Ctrl->I.yinc, NULL, &Ctrl->I.active);

	if (!Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option:  Must specify output file\n", GMT_program);
		error++;
	}
	if (!fin) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input file\n", GMT_program);
		error++;
	}
	if (Ctrl->D.mode < 0 || Ctrl->D.mode > 5) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -D option:  Choose from the range 0-5\n", GMT_program);
		error++;
	}
	if (!Ctrl->F.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR: -F option is required:\n", GMT_program);
		error++;
	}
	if (Ctrl->F.active && Ctrl->F.width == 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -F option:  filter fullwidth must be nonzero:\n", GMT_program);
		error++;
	}
	if (Ctrl->I.active && (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
		error++;
	}
	if (Ctrl->T.active && one_or_zero == 0) {	/* Both -N and -T set, not good */
		fprintf (stderr, "%s: GMT SYNTAX ERROR -T option:  Not allowed with obsolete -N option\n", GMT_program);
		error++;
	}
	if (project_info.region_supplied && Ctrl->I.active && Ctrl->F.highpass) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -F option:  Highpass filtering requires original -R -I\n", GMT_program);
		error++;
	}
	
	if (error) exit (EXIT_FAILURE);

	/* Assign filter_type number */
	
	for (filter_type = 0; filter_type < GRDFILTER_N_FILTERS && filter_code[filter_type] != Ctrl->F.filter; filter_type++);
	
	if (project_info.region_supplied) new_range = TRUE;

	GMT_err_fail (GMT_read_grd_info (fin, &h), fin);
	GMT_grd_init (&h, argc, argv, TRUE);	/* Update command history only */

	if (Ctrl->T.active) {	/* Make output grid of the opposite registration */
		one_or_zero = (h.node_offset) ? 1 : 0;
	}
	else
		one_or_zero = (h.node_offset) ? 0 : 1;

	/* Read the input grid file and close it  */

	nm = GMT_get_nm (h.nx, h.ny);
	input = (float *) GMT_memory (VNULL, (size_t)nm, sizeof(float), GMT_program);

	GMT_err_fail (GMT_read_grd (fin, &h, input, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), fin);
	
	full_360 = (Ctrl->D.mode && GMT_360_RANGE (h.x_max, h.x_min));	/* Periodic geographic grid */
	full_180 = (Ctrl->D.mode && GMT_180_RANGE (h.y_min, h.y_max));	/* Full latitude range for geographic grid */

	last_median = 0.5 * (h.z_min + h.z_max);

	A = (float *) GMT_memory (VNULL, (size_t)nm, sizeof(float), GMT_program);

	/* Check range of output area and set i,j offsets, etc.  */

	if (!new_range) {
		west_new = h.x_min;
		east_new = h.x_max;
		south_new = h.y_min;
		north_new = h.y_max;
	}
	if (!Ctrl->I.active) {
		Ctrl->I.xinc = h.x_inc;
		Ctrl->I.yinc = h.y_inc;
	}

	if (!full_360) {
		if (west_new < h.x_min) error = TRUE;
		if (east_new > h.x_max) error = TRUE;
	}
	if (south_new < h.y_min) error = TRUE;
	if (north_new > h.y_max) error = TRUE;
	if (Ctrl->I.xinc <= 0.0) error = TRUE;
	if (Ctrl->I.yinc <= 0.0) error = TRUE;

	if (error) {
		fprintf(stderr,"%s: New WESN incompatible with old.\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	/* Make sure output grid is kosher */

	test_h.x_min = west_new;	test_h.x_max = east_new;	test_h.x_inc = Ctrl->I.xinc;
	test_h.y_min = south_new;	test_h.y_max = north_new;	test_h.y_inc = Ctrl->I.yinc;
	GMT_RI_prepare (&test_h);	/* Ensure -R -I consistency and set nx, ny */
	GMT_err_fail (GMT_grd_RI_verify (&test_h, 1), Ctrl->G.file);
	/* Copy back in case test_h were changed */
	west_new  = test_h.x_min;	east_new  = test_h.x_max;	Ctrl->I.xinc = test_h.x_inc;
	south_new = test_h.y_min;	north_new = test_h.y_max;	Ctrl->I.yinc = test_h.y_inc;

	/* We can save time by computing a weight matrix once [or once pr scanline] only
	   if new grid spacing is a multiple of old spacing */

	fast_way = (fabs (fmod (Ctrl->I.xinc / h.x_inc, 1.0)) < GMT_SMALL && fabs (fmod (Ctrl->I.yinc / h.y_inc, 1.0)) < GMT_SMALL);
	same_grid = !(new_range  || Ctrl->I.active || h.node_offset == one_or_zero);
	if (!fast_way && gmtdefs.verbose) {
		fprintf (stderr, "%s: Warning - Your output grid spacing is such that filter-weights must\n", GMT_program);
		fprintf (stderr, "be recomputed for every output node, so expect this run to be slow.  Calculations\n");
		fprintf (stderr, "can be speeded up significantly if output grid spacing is chosen to be a multiple\n");
		fprintf (stderr, "of the input grid spacing.  If the odd output grid is necessary, consider using\n");
		fprintf (stderr, "a \'fast\' grid for filtering and then resample onto your desired grid with grdsample.\n");
	}
	if (Ctrl->N.mode == NAN_REPLACE && !same_grid) {
		fprintf (stderr, "%s: Warning: -Nr requires co-registered input/output grids, option is ignored\n", GMT_program);
		Ctrl->N.mode = NAN_IGNORE;
	}
	nx_out = one_or_zero + irint ( (east_new - west_new) / Ctrl->I.xinc);
	ny_out = one_or_zero + irint ( (north_new - south_new) / Ctrl->I.yinc);

	nm = GMT_get_nm (nx_out, ny_out);
	output = (float *) GMT_memory (VNULL, (size_t)nm, sizeof(float), GMT_program);
	i_origin = (GMT_LONG *) GMT_memory (VNULL, (size_t)nx_out, sizeof(GMT_LONG), GMT_program);
	if (!fast_way) x_shift = (double *) GMT_memory (VNULL, (size_t)nx_out, sizeof(double), GMT_program);

	xincnew2 = (one_or_zero) ? 0.0 : 0.5 * Ctrl->I.xinc;
	yincnew2 = (one_or_zero) ? 0.0 : 0.5 * Ctrl->I.yinc;
	xincold2 = (h.node_offset) ? 0.5 * h.x_inc : 0.0;
	yincold2 = (h.node_offset) ? 0.5 * h.y_inc : 0.0;

	if (fast_way && h.node_offset == one_or_zero) {	/* multiple grid but one is pix, other is grid */
		x_fix = 0.5 * h.x_inc;
		y_fix = 0.5 * h.y_inc;
	}
	init_area_weights (&h, Ctrl->D.mode, A);

	/* Set up the distance scalings for lon and lat, and assign pointer to distance function  */

	switch (Ctrl->D.mode) {
		case 0:	/* Plain, unscaled isotropic Cartesian distances */
			x_scale = y_scale = 1.0;
			F.radius_func = CartRadius;
			break;
		case 1:	/* Plain, scaled (degree to km) isotropic Cartesian distances */
			x_scale = y_scale = project_info.DIST_KM_PR_DEG;
			F.radius_func = CartScaledRadius;
			break;
		case 2:	/* Flat Earth Cartesian distances, xscale fixed at mid latitude */
			x_scale = project_info.DIST_KM_PR_DEG * cosd (0.5 * (north_new + south_new));
			y_scale = project_info.DIST_KM_PR_DEG;
			F.radius_func = FlatEarthRadius;
			if (full_360) wrap_case_x = TRUE;		/* For periodic boundaries */
			if (full_180) wrap_case_y = wrap_case_x;	/* For spherical caps */
			break;
		case 3:	/* Flat Earth Cartesian distances, xscale reset for each latitude */
			x_scale = project_info.DIST_KM_PR_DEG * ((fabs (south_new) > north_new) ? cosd (south_new) : cosd (north_new));
			y_scale = project_info.DIST_KM_PR_DEG;
			F.radius_func = FlatEarthRadius;
			if (full_360) wrap_case_x = TRUE;		/* For periodic boundaries */
			if (full_180) wrap_case_y = wrap_case_x;	/* For spherical caps */
			break;
		case 4:	/* Great circle distances */
			x_scale = 0.0;
			y_scale = project_info.DIST_KM_PR_DEG;
			F.radius_func = SphericalRadius;
			if (full_360) wrap_case_x = TRUE;		/* For spherical filtering */
			if (full_180) wrap_case_y = wrap_case_x;	/* For spherical filtering */
			break;
		case 5:	/* Great circle distances with Mercator coordinates */
			/* Get the max |lat| extent of the grid */
			max_lat = IMG2LAT (MAX (fabs (h.y_min), fabs (h.y_max)));
			merc_range = LAT2IMG (max_lat + (0.5 * Ctrl->F.width / project_info.DIST_KM_PR_DEG)) - LAT2IMG (max_lat);
			x_scale = y_scale = 0.5 * Ctrl->F.width / merc_range;
			F.radius_func = SphericalRadius;
			if (full_360) wrap_case_x = TRUE;		/* For spherical filtering; there is no polar caps in img files */
			break;
	}
			
	switch (filter_type) {
		case 1:	/*  Cosine-bell filter weights */
			par[GRDFILTER_INV_R_SCALE] = 2.0 / Ctrl->F.width;
			F.weight_func = CosBellWeight;
			break;
		case 2:	/*  Gaussian filter weights */
			par[GRDFILTER_INV_R_SCALE] = -18.0 / (Ctrl->F.width * Ctrl->F.width);
			F.weight_func = GaussianWeight;
			break;
		default:	/* Everything else uses unit weights */
			F.weight_func = UnitWeight;
			break;
	}
	
	/* Set up miscellaneous filter parameters needed when computing the weights */
	
	par[GRDFILTER_WIDTH] = Ctrl->F.width;
	par[GRDFILTER_HALF_WIDTH] = 0.5 * Ctrl->F.width;
	par[GRDFILTER_X_SCALE] = x_scale;
	par[GRDFILTER_Y_SCALE] = (Ctrl->D.mode == 5) ? project_info.DIST_KM_PR_DEG : y_scale;
	F.d_flag = Ctrl->D.mode;
	F.dx = h.x_inc;
	F.dy = h.y_inc;
	F.y_min = h.y_min;
	F.y_max = h.y_max;
	x_width = Ctrl->F.width / (h.x_inc * x_scale);
	y_width = Ctrl->F.width / (h.y_inc * y_scale);
	F.y_half_width = (GMT_LONG) (ceil(y_width) / 2.0);
	F.x_half_width = (GMT_LONG) (ceil(x_width) / 2.0);
	F.nx = 2 * F.x_half_width + 1;
	F.ny = 2 * F.y_half_width + 1;
	if (x_scale == 0.0 || F.nx < 0 || F.nx > h.nx) {	/* Safety valve when x_scale -> 0.0 */
		F.nx = h.nx;
		F.x_half_width = (F.nx - 1) / 2;
		if ((F.nx - 2 * F.x_half_width - 1) > 0) F.x_half_width++;	/* When nx is even we may come up short by 1 */
		visit_check = ((2 * F.x_half_width + 1) >= h.nx);		/* Must make sure we only visit each node once along a row */
	}
	if (F.ny < 0 || F.ny > h.ny) {
		F.ny = h.ny;
		F.y_half_width = h.ny / 2;
	}
	F.x = (double *) GMT_memory (VNULL, (size_t)(F.x_half_width+1), sizeof (double), GMT_program);
	F.y = (double *) GMT_memory (VNULL, (size_t)(F.y_half_width+1), sizeof (double), GMT_program);
	F.visit = (char *)GMT_memory (VNULL, h.nx, sizeof (char), GMT_program);
	for (i = 0; i <= F.x_half_width; i++) F.x[i] = i * F.dx;
	for (j = 0; j <= F.y_half_width; j++) F.y[j] = j * F.dy;
	
	weight = (double *) GMT_memory (VNULL, (size_t)(F.nx*F.ny), sizeof(double), GMT_program);

	if (filter_type >= 3) {	/* These filters are not convolutions; they require sorting or comparisons */
		slow = TRUE;
		work_array = (double *) GMT_memory (VNULL, (size_t)(F.nx*F.ny), sizeof(double), GMT_program);
	}

	if (wrap_case_x) {	/* Data on a sphere so must check for both periodic and polar wrap-arounds */
		nx_wrap = h.nx - !h.node_offset;	/* So we basically bypass the duplicate point at east for gridline-registered grids */
	}	
	if (gmtdefs.verbose) {
		fprintf(stderr,"%s: Input nx,ny = (%d %d), output nx,ny = (%ld %ld), filter nx,ny = (%ld %ld)\n", GMT_program, h.nx, h.ny, nx_out, ny_out, F.nx, F.ny);
		fprintf(stderr,"%s: Filter type is %s.\n", GMT_program, filter_name[filter_type]);
	}

	/* Compute nearest xoutput i-indices and shifts once */

	for (i_out = 0; i_out < nx_out; i_out++) {
		x_out = west_new + i_out * Ctrl->I.xinc + xincnew2;
		i_origin[i_out] = GMT_x_to_i (x_out, h.x_min, h.x_inc, h.xy_off, h.nx);
		if (!fast_way) x_shift[i_out] = x_out - (h.x_min + i_origin[i_out] * h.x_inc + xincold2);
	}

	/* Determine how much effort to compute weights:
		1 = Compute weights once for entire grid
		2 = Compute weights once per scanline
		3 = Compute weights for every output point [slow]
	*/

	if (fast_way && Ctrl->D.mode <= 2)
		effort_level = 1;
	else if (fast_way && Ctrl->D.mode > 2)
		effort_level = 2;
	else 
		effort_level = 3;
	
	if (effort_level == 1) set_weight_matrix (&F, weight, 0.0, par, x_fix, y_fix);
	half_nx = (h.node_offset) ? h.nx / 2 : (h.nx - 1) / 2;
	
	for (j_out = 0; j_out < ny_out; j_out++) {

		if (gmtdefs.verbose) fprintf (stderr, "%s: Processing output line %ld\r", GMT_program, j_out);
		y_out = north_new - j_out * Ctrl->I.yinc - yincnew2;
		lat_out = (Ctrl->D.mode == 5) ? IMG2LAT (y_out) : y_out;
		j_origin = GMT_y_to_j (y_out, h.y_min, h.y_inc, h.xy_off, h.ny);
		if (Ctrl->D.mode == 3) par[GRDFILTER_X_SCALE] = project_info.DIST_KM_PR_DEG * cosd (lat_out);	/* Update flat-earth longitude scale */

		if (Ctrl->D.mode > 2) {	/* Update max filterweight nodes to deal with at this latitude */
			y = fabs (lat_out);
			if (Ctrl->D.mode == 4) y += (par[GRDFILTER_HALF_WIDTH] / par[GRDFILTER_Y_SCALE]);
			F.x_half_width = (y < 90.0) ? MIN ((F.nx - 1) / 2, irint (par[GRDFILTER_HALF_WIDTH] / (F.dx * par[GRDFILTER_Y_SCALE] * cosd (y)))) : (F.nx - 1) / 2;
			if (y > 90.0 && (F.nx - 2 * F.x_half_width - 1) > 0) F.x_half_width++;	/* When nx is even we may come up short by 1 */
			visit_check = ((2 * F.x_half_width + 1) >= h.nx);	/* Must make sure we only visit each node once along a row */
		}
			
		if (effort_level == 2) set_weight_matrix (&F, weight, y_out, par, x_fix, y_fix);
		if (!fast_way) y_shift = y_out - (h.y_max - j_origin * h.y_inc - yincold2);

		for (i_out = 0; i_out < nx_out; i_out++) {

			if (effort_level == 3) set_weight_matrix (&F, weight, y_out, par, x_shift[i_out], y_shift);
			wt_sum = value = 0.0;
			n_in_median = 0;	go_on = TRUE;
			ij_out = j_out * nx_out + i_out;
			if (Ctrl->N.mode == NAN_REPLACE && GMT_is_fnan (input[ij_out])) {	/* [Here we know ij_out == ij_in]. Since output will be NaN we bypass the filter loop */
				output[ij_out] = GMT_f_NaN;
				n_nan++;
				continue;	/* Done with this node */
			}

			/* Now loop over the filter domain and collect those points that should be considered by the filter operation */
			
			for (jj = -F.y_half_width; go_on && jj <= F.y_half_width; jj++) {
				j_in = j_origin + jj;
				if ( (j_in < 0) || (j_in >= h.ny) ) continue;
				if (visit_check) memset ((void *)F.visit, 0, h.nx);	/* Reset our longitude visit counter */

				for (ii = -F.x_half_width; go_on && ii <= F.x_half_width; ii++) {
					i_in = i_origin[i_out] + ii;
					if (wrap_case_x) {	/* Just wrap around the globe */
						if (i_in < 0) i_in += nx_wrap;
						else if (i_in >= nx_wrap) i_in -= nx_wrap;
						i_orig = i_in;
					}
					if ( (i_in < 0) || (i_in >= h.nx)) continue;
				
					if (visit_check) {	/* Make sure we never include the same node twice along a given row */
						if (F.visit[i_in]) continue;		/* Already been used */
						F.visit[i_in] = 1;			/* Now marked as visited */
					}
					ij_wt = (jj + F.y_half_width) * F.nx + ii + F.x_half_width;
					if (weight[ij_wt] <= 0.0) continue;

					ij_in = j_in*h.nx + i_in;
					if (GMT_is_fnan (input[ij_in])) {
						if (Ctrl->N.mode == NAN_PRESERVE) go_on = FALSE;	/* -Np means no NaNs are allowed */
						continue;
					}

					/* Get here when point is usable  */
					if (slow) {
						work_array[n_in_median] = input[ij_in];
						n_in_median++;
					}
					else {
						w = weight[ij_wt] * A[ij_in];
						value += input[ij_in] * w;
						wt_sum += w;
					}
				}
			}

			/* Now we have done the convolution and we can get the value  */

			if (!go_on) {	/* -Np in effect and there were NaNs inside circle */
				output[ij_out] = GMT_f_NaN;
				n_nan++;
			}
			else if (slow) {
				if (n_in_median) {
					switch (filter_type) {
						case 3:	/* Median */
							GMT_median (work_array, n_in_median, h.z_min, h.z_max, last_median, &this_median);
							last_median = this_median;
							break;
						case 4:	/* Mode */
							GMT_mode (work_array, n_in_median, n_in_median/2, TRUE, Ctrl->F.mode, &GMT_n_multiples, &this_median);
							break;
						case 5:	/* Lowest of all */
							this_median = GMT_extreme (work_array, n_in_median, DBL_MAX, 0, -1);
							break;
						case 6:	/* Lowest of positive values */
							this_median = GMT_extreme (work_array, n_in_median, 0.0, +1, -1);
							break;
						case 7:	/* Upper of all values */
							this_median = GMT_extreme (work_array, n_in_median, -DBL_MAX, 0, +1);
							break;
						case 8:	/* Upper of negative values */
							this_median = GMT_extreme (work_array, n_in_median, 0.0, -1, +1);
							break;
					}
					output[ij_out] = (float)this_median;
				}
				else {
					output[ij_out] = GMT_f_NaN;
					n_nan++;
				}
			}
			else {
				if (wt_sum == 0.0) {	/* Assign value = GMT_f_NaN */
					output[ij_out] = GMT_f_NaN;
					n_nan++;
				}
				else
					output[ij_out] = (float)(value / wt_sum);
			}
		}
	}
	if (gmtdefs.verbose) fprintf (stderr, "\n");

	if (Ctrl->F.highpass) {
		if (gmtdefs.verbose) fprintf (stderr, "%s: Subtracting lowpass-filtered data from grid to obtain high-pass filtered data\n", GMT_program);
		for (ij_out = 0; ij_out < nx_out * ny_out; ij_out++) output[ij_out] = input[ij_out] - output[ij_out];
	}
	
	/* At last, that's it!  Output: */

	if (n_nan && gmtdefs.verbose) fprintf (stderr, "%s: Unable to estimate value at %ld nodes, set to NaN\n", GMT_program, n_nan);
	if (GMT_n_multiples > 0 && gmtdefs.verbose) fprintf (stderr, "%s: WARNING: %ld multiple modes found\n", GMT_program, GMT_n_multiples);

	h.nx = (int)nx_out;
	h.ny = (int)ny_out;
	h.x_min = west_new;
	h.x_max = east_new;
	h.y_min = south_new;
	h.y_max = north_new;
	h.x_inc = Ctrl->I.xinc;
	h.y_inc = Ctrl->I.yinc;
	h.node_offset = !one_or_zero;

	GMT_err_fail (GMT_write_grd (Ctrl->G.file, &h, output, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);

	GMT_free ((void *) A);
	GMT_free ((void *) input);
	GMT_free ((void *) output);
	GMT_free ((void *) weight);
	GMT_free ((void *) i_origin);
	GMT_free ((void *) F.x);
	GMT_free ((void *) F.y);
	GMT_free ((void *) F.visit);
	if (slow) GMT_free ((void *) work_array);
	if (!fast_way) GMT_free ((void *) x_shift);

	Free_grdfilter_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void set_weight_matrix (struct FILTER_INFO *F, double *weight, double output_lat, double par[], double x_off, double y_off)
{
	/* x_off and y_off give offset between output node and 'origin' input node for this window (0,0 for integral grids).
	 * fast is TRUE when input/output grids are offset by integer values in dx/dy.
	 * Here, par[0] = filter_width, par[1] = filter_width / 2, par[2] = x_scale, part[3] = y_scale, and
	 * par[4] is the normalization distance needed for the Cosine-bell or Gaussian weight function.
	 */

	GMT_LONG i, j, ij;
	double x, y, yc, y0, r;

	yc = y0 = output_lat - y_off;		/* Input latitude of central point (i,j) = (0,0) */
	if (F->d_flag == 5) yc = IMG2LAT (yc);	/* Recover actual latitude in IMG grid at this center point */
	for (j = -F->y_half_width; j <= F->y_half_width; j++) {
		y = y0 + ((j < 0) ? F->y[-j] : -F->y[j]);	/* y or latitude at this row */
		if (F->d_flag > 2 && (y < F->y_min || y > F->y_max)) {		/* This filter row is outside input grid domain */
			for (i = -F->x_half_width, ij = (j + F->y_half_width) * F->nx; i <= F->x_half_width; i++, ij++) weight[ij] = -1.0;
			continue;	/* Done with this row */
		}
		if (F->d_flag == 5) y = IMG2LAT (y);	/* Recover actual latitudes */
		for (i = -F->x_half_width; i <= F->x_half_width; i++) {
			x = (i < 0) ? -F->x[-i] : F->x[i];
			ij = (j + F->y_half_width) * F->nx + i + F->x_half_width;
			r = F->radius_func (x_off, yc, x, y, par);
			weight[ij] = (r > par[GRDFILTER_HALF_WIDTH]) ? -1.0 : F->weight_func (r, par);
		}
	}
}

double CartRadius (double x0, double y0, double x1, double y1, double par[])
{	/* Plain Cartesian distance */
	return (hypot (x0 - x1, y0 - y1));
}

double CartScaledRadius (double x0, double y0, double x1, double y1, double par[])
{	/* Plain scaled Cartesian distance (xscale = yscale) */
	return (par[GRDFILTER_X_SCALE] * hypot (x0 - x1, y0 - y1));
}

double FlatEarthRadius (double x0, double y0, double x1, double y1, double par[])
{	/* Cartesian radius with different scales */
	return (hypot (par[GRDFILTER_X_SCALE] * (x0 - x1), par[GRDFILTER_Y_SCALE] * (y0 - y1)));
}

double SphericalRadius (double x0, double y0, double x1, double y1, double par[])
{	/* Great circle distance with polar wrap-around test on 2nd point */
	if (fabs (y1) > 90.0) {	/* Must find the point across the pole */
		y1 = copysign (180.0 - fabs (y1), y1);
		x1 += 180.0;
	}
	return (GMT_great_circle_dist_km (x0, y0, x1, y1));
}

double UnitWeight (double r, double par[])
{
	/* Return unit weight since we know r is inside radius */

	return (1.0);
}

double CosBellWeight (double r, double par[])
{
	/* Return the cosine-bell filter weight for given r.
	 * The parameter r_f_half is the 5th parameter passed.
	 */

	return (1.0 + cos (M_PI * r * par[GRDFILTER_INV_R_SCALE]));
}

double GaussianWeight (double r, double par[])
{
	/* Return the Gaussian filter weight for given r.
	 * The parameter sig_2 is the 5th parameter passed.
	 */

	return (exp (r * r * par[GRDFILTER_INV_R_SCALE]));
}

void init_area_weights (struct GRD_HEADER *G, GMT_LONG mode, float *A)
{
	/* Precalculate the area weight of each node.  There are several considerations:
	 * 1. Mercator img grids (d_flag == 5) has irregular latitude spacing so we
	 *    must compute the n and s latitudes for each cell to get area.
	 * 2. Poles for grid-registered grids require a separate weight formula
	 * 3. Grid-registered grids have boundary nodes that only apply to 1/2 the area
	 *    (and the four corners (unless poles) only 1/4 the area of other cells).
	 */
	GMT_LONG row, col, ij;
	double row_weight, col_weight, dy_half = 0.0, dx, y, lat, lat_s, lat_n, s2 = 0.0;
	
	if (mode) {	/* Geographic data */
		if (mode == 5) dy_half = 0.5 * G->y_inc;	/* Half img y-spacing */
		dx = G->x_inc * R2D;				/* Longitude increment in radians */
		s2 = sind (0.5 * G->y_inc);			/* Holds sin (del_y/2) */
	}
	else
		dx = G->x_inc;
	for (row = 0; row < G->ny; row++) {
		if (mode == 5) {		/* Adjust lat if IMG grid.  Note: these grids do not reach a pole. */
			y = GMT_j_to_y (row, G->y_min, G->y_max, G->y_inc, G->xy_off, G->ny);	/* Current input Merc y */
			lat = IMG2LAT (y);			/* Get actual latitude */
			lat_s = IMG2LAT (y - dy_half);		/* Bottom grid cell latitude */
			lat_n = IMG2LAT (y + dy_half);		/* Top grid cell latitude */
			row_weight = sind (lat_n) - sind (lat_s);
		}
		else if (mode) {	/* Geographic data, watch for poles */
			lat = GMT_j_to_y (row, G->y_min, G->y_max, G->y_inc, G->xy_off, G->ny);	/* Current input latitude */
			if (fabs (lat) == 90.0)	/* Poles are different */
				row_weight = 1.0 - cosd (0.5 * G->y_inc);
			else {	/* All other points */
				row_weight = 2.0 * cosd (lat) * s2;
				/* Note: No need for special weight-sharing between w/e gridline-reg grids since we explicitly only use the west node */
			}
		}
		else {	/* If not geographic then row_weight is a constant 1 except for gridline-registered grids at the ends */
			row_weight = (G->node_offset == 0 && (row == 0 || row == (G->ny-1))) ? 0.5 : 1.0;	/* Share weight with repeat point */
			row_weight *= G->y_inc;
		}

		for (col = 0, ij = row * G->nx; col < G->nx; col++, ij++) {
			col_weight = dx * ((G->node_offset == 0 && (col == 0 || col == (G->nx-1))) ? 0.5  : 1.0);
			A[ij] = (float)(row_weight * col_weight);
		}
	}
}

void *New_grdfilter_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDFILTER_CTRL *C;
	
	C = (struct GRDFILTER_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDFILTER_CTRL), "New_grdfilter_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->D.mode = -1;	
	return ((void *)C);
}

void Free_grdfilter_Ctrl (struct GRDFILTER_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}
