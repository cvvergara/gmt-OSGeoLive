/*--------------------------------------------------------------------
 *	$Id: greenspline.c 9923 2012-12-18 20:45:53Z pwessel $
 *
 *	Copyright (c) 2008-2013 by P. Wessel
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
 * greenspline grids data using Green's functions for a selected spline.
 * The data may be Cartesian or geographical, and gridding can be done
 * in a Cartesian 1-D, 2-D, 3-D space or on a sphere.  The spline may be evaluated on
 * a grid, on parts of a grid, or at specified arbitrary locations.
 * Five classes of splines (for a total of 10 splines) are implemented:
 * 1. Minimum curvature Cartesian spline [Sandwell, Geophys. Res. Lett, 1987]
 * 2. Minimum curvature Cartesian spline in tension [Wessel & Bercovici, Math., Geol., 1998]
 * 3. Regularized Cartesian spline in tension [Mitasova & Mitas, Math. Geol., 1993]
 * 4. Minimum curvature spherical spline [Parker, "Geophysical Inverse Theory", 1994]
 * 5. Minimum curvature spherical spline in tension [Wessel & Becker, Geophys. J. Int, 2008]
 *
 * Originally published as:
 *   "Wessel, P., 2009. A general-purpose Green's function-based interpolator,
 *	Computers & Geosciences, 35: 1247–1254".
 *
 * Author:	Paul Wessel
 * Date:	05-DEC-2008
 * Version:	1.0
 *
 */

#include "gmt.h"

#define SANDWELL_1987_1D		0
#define SANDWELL_1987_2D		1
#define SANDWELL_1987_3D		2
#define WESSEL_BERCOVICI_1998_1D	3
#define WESSEL_BERCOVICI_1998_2D	4
#define WESSEL_BERCOVICI_1998_3D	5
#define MITASOVA_MITAS_1993_2D		6
#define MITASOVA_MITAS_1993_3D		7
#define PARKER_1994			8
#define WESSEL_BECKER_2008		9

#define N_METHODS			10

#ifndef M_LOG_2
#define M_LOG_2 0.69314718055994530942
#endif
#ifndef M_GAMMA
#define M_GAMMA 0.577215664901532860606512
#endif
#ifndef M_SQRT_PI
#define M_SQRT_PI 1.772453850905516027298167483341
#endif
#ifndef M_INV_SQRT_PI
#define M_INV_SQRT_PI (1.0 / M_SQRT_PI)
#endif

#define N_X 	100001

typedef double REAL;	/* Change this to float and recompile if you only want single precision calculations */
GMT_LONG TEST = FALSE;	/* Global variable for TESTING only, used to reproduce all the Green's functions and their derivatives (Fig. 1, in Wessel, 2009) */

struct GREENSPLINE_CTRL {
	struct A {	/* -A<gradientfile> */
		GMT_LONG active;
		GMT_LONG mode;	/* 0 = azimuths, 1 = directions, 2 = dx,dy components, 3 = dx, dy, dz components */
		char *file;
	} A	;
	struct C {	/* -C<cutoff> */
		GMT_LONG active;
		double value;
		char *file;
	} C;
	struct D {	/* -D<distflag> */
		GMT_LONG active;
		GMT_LONG mode;
	} D;
	struct F {	/* -F */
		GMT_LONG active;
	} F;
	struct G {	/* -G<output_grdfile> */
		GMT_LONG active;
		char *file;
	} G;
	struct I {	/* -Idx[/dy[/dz]] */
		GMT_LONG active;
		double inc[3];
	} I;
	struct L {	/* -L */
		GMT_LONG active;
	} L;
	struct N {	/* -N<outputnode_file> */
		GMT_LONG active;
		char *file;
	} N;
	struct Q {	/* -Qdaz */
		GMT_LONG active;
		double az;
		double dir[3];
	} Q;
	struct S {	/* -S<mode>[/args] */
		GMT_LONG active, fast;
		GMT_LONG mode;
		double value[2];
		char *arg;
	} S;
	struct T {	/* -T<mask_grdfile> */
		GMT_LONG active;
		char *file;
	} T;
};

struct ZGRID {
	GMT_LONG nz;
	double z_min, z_max, z_inc;
};

#ifdef DEBUG
void dump_green (PFD G, PFD D, double par[], double x0, double x1, GMT_LONG N, double *zz, double *gg);
#endif

/* Functions for complex math */

static void Cdiv (double A[], double B[], double C[])
{	/* Complex division */
	double i_denom;
	i_denom = 1.0 / (B[0]*B[0] + B[1]*B[1]);
	C[0] = (A[0]*B[0] + A[1]*B[1]) * i_denom;
	C[1] = (A[1]*B[0] - A[0]*B[1]) * i_denom;
}

static void Cmul (double A[], double B[], double C[])
{	/* Complex multiplication */
	C[0] = A[0]*B[0] - A[1]*B[1];
	C[1] = A[0]*B[1] + A[1]*B[0];
}

static void Ccot (double Z[], double cotZ[])
{	/* Complex cot(z) */
	double sx, cx, e, A[2], B[2];
	
	sincos (2.0*Z[0], &sx, &cx);
	e = exp (-2.0*Z[1]);
	A[0] = -e * sx;		A[1] = B[0] = e * cx;
	A[1] += 1.0;	B[0] -= 1.0;	B[1] = -A[0];
	Cdiv (A, B, cotZ);
}
void GMT_check_lattice_cpy (struct GMT_GRD_INFO *info, double *x_inc, double *y_inc, GMT_LONG *pixel, GMT_LONG *active);

PFD GMT_azimuth_func;

int main (int argc, char **argv)
{
	GMT_LONG i, j, k, p, ii, n, m, nm, fno, n_files = 0, n_fields, n_read, n_args, dimension = 0;
	GMT_LONG error = 0, n_expected_fields, normalize = 1, unit = 0, n_items;
	GMT_LONG old_n_alloc, n_alloc, ij, ji, nxy, n_ok = 0;
	GMT_LONG do_grid, done, nofile = TRUE;
	char line[BUFSIZ];
	char *method[N_METHODS] = {"minimum curvature Cartesian spline [1-D]",
		"minimum curvature Cartesian spline [2-D]",
		"minimum curvature Cartesian spline [3-D]",
		"continuous curvature Cartesian spline in tension [1-D]",
		"continuous curvature Cartesian spline in tension [2-D]",
		"continuous curvature Cartesian spline in tension [3-D]",
		"regularized Cartesian spline in tension [2-D]",
		"regularized Cartesian spline in tension [3-D]",
		"minimum curvature spherical spline",
		"continuous curvature spherical spline in tension"};
	char *mem_unit[3] = {"kb", "Mb", "Gb"}, txt[6][GMT_TEXT_LEN];
	float *w = NULL;
	double *obs, **D = NULL, **X, *alpha, *in = NULL, r, par[11], norm[7], az, grad;
	double *WB_z = NULL, *WB_g = NULL, mem, r_min = -1.0, r_max = 1.0, part, C, p_val;
	REAL *A = NULL;
	double x0 = 0.0, x1 = 5.0;
	struct GRD_HEADER h;
	struct ZGRID Z;
	struct GMT_TABLE *T = NULL, *S = NULL;
	struct GMT_GRD_INFO info;
	FILE *fp = NULL;

	struct GREENSPLINE_CTRL *Ctrl = NULL;

	void *New_greenspline_Ctrl (), Free_greenspline_Ctrl (struct GREENSPLINE_CTRL *C);

	double spline1d_sandwell (double r, double par[], double *z);
	double spline1d_Wessel_Bercovici (double r, double par[], double *z);
	double gradspline1d_sandwell (double r, double par[], double *z);
	double gradspline1d_Wessel_Bercovici (double r, double par[], double *z);
	double spline2d_sandwell (double r, double par[], double *z);
	double spline2d_Wessel_Bercovici (double r, double par[], double *z);
	double spline2d_Mitasova_Mitas (double r, double par[], double *z);
	double spline2d_Parker (double x, double par[], double *z);
	double spline2d_Wessel_Becker (double x, double par[], double *z);
	double gradspline2d_sandwell (double r, double par[], double *z);
	double gradspline2d_Wessel_Bercovici (double r, double par[], double *z);
	double gradspline2d_Mitasova_Mitas (double r, double par[], double *z);
	double gradspline2d_Parker (double x, double par[], double *z);
	double gradspline2d_Wessel_Becker (double x, double par[], double *z);
	double spline2d_lookup (double x, double par[], double *y);
	double spline2d_Wessel_Becker_lookup (double x, double par[], double *z);
	double gradspline2d_Wessel_Becker_lookup (double x, double par[], double *z);
	void spline2d_Wessel_Becker_init (double par[], double *z, double *g, GMT_LONG grad);
	double spline3d_sandwell (double r, double par[], double *G);
	double spline3d_Wessel_Bercovici (double r, double par[], double *G);
	double spline3d_Mitasova_Mitas (double r, double par[], double *G);
	double gradspline3d_sandwell (double r, double par[], double *G);
	double gradspline3d_Wessel_Bercovici (double r, double par[], double *G);
	double gradspline3d_Mitasova_Mitas (double r, double par[], double *G);
	double undo_normalization (double *X, double q, GMT_LONG mode, double *coeff);
	void do_normalization (double **X, double *obs, GMT_LONG n, GMT_LONG mode, double *coeff);
	double get_radius (double *X0, double *X1, GMT_LONG dim);
	double get_dircosine (double *D, double *X0, double *X1, GMT_LONG dim, GMT_LONG baz);
	
	PFD G = NULL, dGdr = NULL;
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GREENSPLINE_CTRL *)New_greenspline_Ctrl ();	/* Allocate and initialize a new control structure */

	GMT_grd_init (&h, argc, argv, FALSE);
	memset ((void *)par,  0, (size_t)(7*sizeof(double)));
	memset ((void *)norm, 0, (size_t)(7*sizeof(double)));
	memset ((void *)&info, 0, sizeof(struct GMT_GRD_INFO));
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */
                      
				case 'H':
				case 'V':
				case ':':
				case 'b':
				case '\0':
					error += GMT_parse_common_options (argv[i], NULL, NULL, NULL, NULL);
					break;
                              
				case 'R':	/* Must be handled separately since it can take 1,2,3 dimensions */
					project_info.region_supplied = TRUE;
				
					if (argv[i][2] == 'g' && argv[i][3] == '\0') {	/* Got -Rg */
						h.x_min = 0.0;	h.x_max = 360.0;	h.y_min = -90.0;	h.y_max = 90.0;
						break;
					}
					if (argv[i][2] == 'd' && argv[i][3] == '\0') {	/* Got -Rd */
						h.x_min = -180.0;	h.x_max = 180.0;	h.y_min = -90.0;	h.y_max = 90.0;
						break;
					}
					if (!GMT_access (&argv[i][2], R_OK)) {	/* Gave a readable file, presumably a grid */		
						GMT_err_fail (GMT_read_grd_info (&argv[i][2], &info.grd), &argv[i][2]);
						h.x_min = info.grd.x_min; h.x_max = info.grd.x_max;
						h.y_min = info.grd.y_min; h.y_max = info.grd.y_max;
						info.active = TRUE;
						break;
					}
					
					n_items = sscanf (&argv[i][2], "%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%s", txt[0], txt[1], txt[2], txt[3], txt[4], txt[5]);
					if (!(n_items == 2 || n_items == 4 || n_items == 6)) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -R option:  Give 2, 4, or 6 coordinates\n", GMT_program);
						exit (EXIT_FAILURE);
					}
					error += GMT_verify_expectations (GMT_io.in_col_type[GMT_X], GMT_scanf_arg (txt[0], GMT_io.in_col_type[GMT_X], &h.x_min), txt[0]);
					error += GMT_verify_expectations (GMT_io.in_col_type[GMT_X], GMT_scanf_arg (txt[1], GMT_io.in_col_type[GMT_X], &h.x_max), txt[1]);
					if (n_items > 2) {
						error += GMT_verify_expectations (GMT_io.in_col_type[GMT_Y], GMT_scanf_arg (txt[2], GMT_io.in_col_type[GMT_Y], &h.y_min), txt[2]);
						error += GMT_verify_expectations (GMT_io.in_col_type[GMT_Y], GMT_scanf_arg (txt[3], GMT_io.in_col_type[GMT_Y], &h.y_max), txt[3]);
					}
					if (n_items == 6) {
						error += GMT_verify_expectations (GMT_io.in_col_type[GMT_Z], GMT_scanf_arg (txt[4], GMT_io.in_col_type[GMT_Z], &Z.z_min), txt[4]);
						error += GMT_verify_expectations (GMT_io.in_col_type[GMT_Z], GMT_scanf_arg (txt[5], GMT_io.in_col_type[GMT_Z], &Z.z_max), txt[5]);
					}
					break;
								
				/* Supplemental parameters */
                              
				case 'A':
					Ctrl->A.active = TRUE;
					j = 2;
					if (strchr (argv[i], ',')) {	/* Specified a particular format with -A<mode>,<file> */
						Ctrl->A.mode = (int)(argv[i][2] - '0');
						j = 4;
					}
					Ctrl->A.file = strdup (&argv[i][j]);
					break;
				case 'C':
					Ctrl->C.active = TRUE;
					if (strchr (argv[i], '/')) {
						char tmp[BUFSIZ];
						sscanf (&argv[i][2], "%lf/%s", &Ctrl->C.value, tmp);
						Ctrl->C.file = strdup (tmp);
					}
					else
						Ctrl->C.value = atof (&argv[i][2]);
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					Ctrl->D.mode = atoi(&argv[i][2]);	/* Since I added 0 to be 1-D later so no it is -1 */
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					k = GMT_getincn (&argv[i][2], Ctrl->I.inc, 3);
					if (k < 1) {
						GMT_inc_syntax ('I', 1);
						error = TRUE;
					}
					if (Ctrl->I.inc[GMT_Y] == 0.0) Ctrl->I.inc[GMT_Y] = Ctrl->I.inc[GMT_X];
					if (Ctrl->I.inc[GMT_Z] == 0.0) Ctrl->I.inc[GMT_Z] = Ctrl->I.inc[GMT_X];
					break;
				case 'L':
					Ctrl->L.active = TRUE;
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					Ctrl->N.file = strdup (&argv[i][2]);
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
					if (strchr (argv[i], '/')) {	/* Got 3-D vector components */
						k = sscanf (&argv[i][2], "%lf/%lf/%lf", &Ctrl->Q.dir[0], &Ctrl->Q.dir[1], &Ctrl->Q.dir[2]);
						if (k != 3) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR -Q option:  Append azimuth (2-D) or x/y/z components (3-D)\n", GMT_program);
							exit (EXIT_FAILURE);
						}
						GMT_normalize3v (Ctrl->Q.dir);	/* Normalize to unit vector */
					}
					else if (!argv[i][2]) {	/* No argument given*/
						fprintf (stderr, "%s: GMT SYNTAX ERROR -Q option:  Append azimuth (2-D) or x/y/z components (3-D)\n", GMT_program);
						exit (EXIT_FAILURE);
					}
					else	/* 2-D azimuth */
						Ctrl->Q.az = atof(&argv[i][2]);
					break;
				case 'S':
					Ctrl->S.arg = strdup (argv[i]);
					switch (argv[i][2]) {
						case 'c':
							Ctrl->S.mode = SANDWELL_1987_1D;
							break;
						case 't':
							Ctrl->S.mode = WESSEL_BERCOVICI_1998_1D;
							if (strchr (argv[i], '/'))
								sscanf (&argv[i][3], "%lf/%lf", &Ctrl->S.value[0], &Ctrl->S.value[1]);
							else
								Ctrl->S.value[0] = atof (&argv[i][3]);
							break;
						case 'r':
							Ctrl->S.mode = MITASOVA_MITAS_1993_2D;
							if (strchr (argv[i], '/'))
								sscanf (&argv[i][3], "%lf/%lf", &Ctrl->S.value[0], &Ctrl->S.value[1]);
							else
								Ctrl->S.value[0] = atof (&argv[i][3]);
							break;
						case 'p':
							Ctrl->S.mode = PARKER_1994;
							break;
						case 'Q':
							Ctrl->S.mode = WESSEL_BECKER_2008;
							Ctrl->S.fast = TRUE; 
							if (strchr (argv[i], '/')) {
								k = sscanf (&argv[i][3], "%lf/%lf/%lf/%lf", &Ctrl->S.value[0], &Ctrl->S.value[1], &r_min, &r_max);
								if (k == 2) {
									r_min = -1.0;
									r_max = +1.0;
								}
							}
							else {
								Ctrl->S.value[0] = atof (&argv[i][3]);
								Ctrl->S.value[1] = (double)N_X;
							}
							break;
						case 'q':
							Ctrl->S.mode = WESSEL_BECKER_2008;
							Ctrl->S.value[0] = atof (&argv[i][3]);
							if (Ctrl->S.value[0] == 0.0) {
								Ctrl->S.mode = PARKER_1994;
							}
							break;
						default:
							fprintf (stderr, "%s: GMT SYNTAX ERROR -D option:  Append c|t|g|p|q\n", GMT_program);
							exit (EXIT_FAILURE);
						break;
					}
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					Ctrl->T.file = strdup (&argv[i][2]);
					break;
				case '+':	/* Undocumented option for testing only */
					TEST = TRUE;
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
		fprintf (stderr, "greenspline %s - Interpolation using a Green's function for splines in 1-3 dimensions\n\n", GMT_VERSION);
		fprintf (stderr, "usage: greenspline [xyzfile] [-A[<format>,]<gradientfile>] [-R<xmin>/<xmax[/<ymin>/<ymax>[/<zmin>/<zmax>]]] [-I<dx>[/<dy>[/<dz>]] -Goutfile\n");
		fprintf (stderr, "\t[-C<cut>[/<file>]] [-D<mode>] [-F] [%s] [%s] [-L] [-N<nodes>]\n", GMT_H_OPT, GMT_I_OPT);
		fprintf (stderr, "\t[-Q<az>] [-Sc|t|r|p|q[<pars>]] [-T<maskfile>] [-V] [%s] [%s]\n\n", GMT_t_OPT, GMT_bi_OPT);
		
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
		
		fprintf (stderr, "\tgreenspline will read from standard input or specified ASCII xyz-file(s) (see -bi if binary).\n\n");
		fprintf (stderr, "\tChoose one of three ways to specify where to evaluate the spline:\n");
		fprintf (stderr, "\t1. Specify a rectangular grid domain with options -R, -I and optionally -F.\n");
		fprintf (stderr, "\t2. Supply a mask file via -T whose values are NaN or 0.  The spline will then\n");
		fprintf (stderr, "\t   only be evaluated at the nodes originally set to zero.\n");
		fprintf (stderr, "\t3. Specify a set of output locations via the -N option.\n\n");
		fprintf (stderr, "\t-G Output data. Give name of output file.\n");

		fprintf (stderr, "\n\tOPTIONS:\n");

		fprintf (stderr, "\t-A ASCII file with surface gradients V to use in the modeling.  Specify format:\n");
		fprintf (stderr, "\t   (0) For 1-D: x, slope, (1) X, Vmagnitude, Vazimuth(s), (2) X, Vazimuth(s), Vmagnitude,\n");
		fprintf (stderr, "\t   (3) X, Vmagnitude, Vangle(s), (4) X, Vcomponents, or (5) X, Vunit-vector, Vmagnitude.\n");
		fprintf (stderr, "\t   Here, X = (x, y[, z]) is the position vector, V is the gradient vector.\n");
		fprintf (stderr, "\t-C Solve by SVD and eliminate eigenvalues whose ratio to largest eigenvalue is less than <cut>.\n");
		fprintf (stderr, "\t   Optionally append /<filename> to save the eigenvalues to this file.\n");
		fprintf (stderr, "\t   A negative cutoff will stop execution after saving the eigenvalues\n");
		fprintf (stderr, "\t   [Default uses Gauss-Jordan elimination to solve the linear system].\n");
		fprintf (stderr,"\t-D Distance flag determines how we calculate distances between (x,y) points:\n");
		fprintf (stderr,"\t   Option 0 applies to Cartesian 1-D spline interpolation:\n");
		fprintf (stderr,"\t     -D0 x in user units, Cartesian distances.\n");
		fprintf (stderr,"\t   Options 1-3 apply to Cartesian 2-D surface spline interpolation:\n");
		fprintf (stderr,"\t     -D1 x,y in user units, Cartesian distances.\n");
		fprintf (stderr,"\t     -D2 x,y in degrees, flat Earth distances.\n");
		fprintf (stderr,"\t     -D3 x,y in degrees, spherical distances in km.\n");
		fprintf (stderr,"\t   Option 4 applies to 2-D spherical surface spline interpolation:\n");
		fprintf (stderr,"\t     -D4 x,y in degrees, use cosine of spherical distances in degrees.\n");
		fprintf (stderr,"\t   Option 5 applies to Cartesian 3-D volume interpolation:\n");
		fprintf (stderr,"\t     -D5 x,y,z in user units, Cartesian distances.\n");
		fprintf (stderr,"\t   For option 3-4, use ELLIPSOID to select geodesic or great cicle arcs.\n");
		fprintf (stderr, "\t-F will force pixel registration [Default is grid registration] (requires -R, -I).\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I Specifies a regular set of output locations.  Give equidistant increment for each dimension.\n");
		fprintf (stderr, "\t   Requires -R for specifying the output domain.\n");
		fprintf (stderr, "\t-L Leave trend alone.  Do not remove least squares plane from data before spline fit.\n");
		fprintf (stderr, "\t   Only applies to -D0-2 [Default removes linear trend, fits residuals, and restores trend].\n");
		fprintf (stderr, "\t-N ASCII file with desired output locations.\n");
		fprintf (stderr, "\t   The resulting ASCII coordinates and interpolation are written to file given in -G\n");
		fprintf (stderr, "\t   or stdout if no file specified (see -bo for binary output).\n");
		fprintf (stderr, "\t-Q Calculate the directional derivative in the <az> direction and return it instead of surface elevation.\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t-R Specify a regular set of output locations.  Give min and max coordinates for each dimension.\n");
		fprintf (stderr, "\t   Requires -I for specifying equidistant increments.  For 2D-gridding a gridfile may be given;\n");
		fprintf (stderr, "\t   this then also sets -I (and perhaps -F); use those options to override the grid settings.\n");
		fprintf (stderr, "\t-S specifies which spline to use (if needed, normalized <tension> must be between 0 and 1):\n");
		fprintf (stderr, "\t   -Sc is minimum curvature spline (Sandwell, 1987) [Default].\n");
		fprintf (stderr, "\t   -St<tension>[/<scale>] is spline in tension (Wessel & Bercovici, 1998).\n");
		fprintf (stderr, "\t      Optionally, specify a length-scale [Default is grid spacing].\n");
		fprintf (stderr, "\t   -Sr<tension> is a regularized spline in tension (Mitasova & Mitas, 1993).\n");
		fprintf (stderr, "\t      Optionally, specify a length-scale [Default is grid spacing].\n");
		fprintf (stderr, "\t   -Sp is a spherical surface spline (Parker, 1994); automatically sets -D4.\n");
		fprintf (stderr, "\t   -Sq is a spherical surface spline in tension (Wessel & Becker, 2008); automatically sets -D4.\n");
		fprintf (stderr, "\t      Use -SQ to speed up calculations by using precalculated lookup tables.\n");
		fprintf (stderr, "\t      Append /n to set the (odd) number of points in the spline [%d].\n", N_X);
		fprintf (stderr, "\t-T Mask grid file whose values are NaN or 0; its header implicitly sets -R, -I (and -F).\n");
		GMT_explain_option ('V');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2-4 input columns depending on dimensionality (see -D).\n");
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (! (project_info.region_supplied || Ctrl->N.active || Ctrl->T.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  No output locations specified (use either [-R -I], -N, or -T)\n", GMT_program);
		error++;
	}
	dimension = (Ctrl->D.mode == 0) ? 1 : ((Ctrl->D.mode == 5) ? 3 : 2);
	if (info.active && dimension != 2) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  The -R<gridfile> option only applies to 2-D gridding\n", GMT_program);
		error++;
	}
	if (dimension == 2) GMT_check_lattice_cpy (&info, &Ctrl->I.inc[GMT_X], &Ctrl->I.inc[GMT_Y], &Ctrl->F.active, &Ctrl->I.active);
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = dimension + 1;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] <= dimension) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least %ld columns\n", GMT_program, dimension);
		error++;
	}
	if (Ctrl->S.value[0] < 0.0 || Ctrl->S.value[0] >= 1.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option.  Tension must be in range 0 <= t < 1\n", GMT_program);
		error = TRUE;
	}
	if (!(Ctrl->S.mode == PARKER_1994 || Ctrl->S.mode == WESSEL_BECKER_2008) && Ctrl->D.mode == 3) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Sc|t|r option.  Cannot select -D3\n", GMT_program);
		error = TRUE;
	}
	if (Ctrl->I.active && (Ctrl->I.inc[GMT_X] <= 0.0 || (dimension > 1 && Ctrl->I.inc[GMT_Y] <= 0.0) || (dimension == 3 && Ctrl->I.inc[GMT_Z] <= 0.0))) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
		error = TRUE;
	}
	if (dimension == 2 && !Ctrl->N.active && !(Ctrl->G.active  || Ctrl->G.file)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option.  Must specify output grid file name\n", GMT_program);
		error = TRUE;
	}
	if (Ctrl->C.active && Ctrl->C.value < 0.0 && !Ctrl->C.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -C option.  Must specify file name for eigenvalues if cut < 0\n", GMT_program);
		error = TRUE;
	}
	if (Ctrl->T.active && !Ctrl->T.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -T option.  Must specify mask grid file name\n", GMT_program);
		error = TRUE;
	}
	if (Ctrl->N.active && !Ctrl->N.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -N option.  Must specify node file name\n", GMT_program);
		error = TRUE;
	}
	if ((do_grid = (Ctrl->I.active + project_info.region_supplied)) == 1 && dimension == 2) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Must specify -R, -I, [-F], -G for gridding\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

#ifdef SET_IO_MODE
	GMT_setmode (GMT_OUT);
#endif

	if (Ctrl->S.mode == SANDWELL_1987_1D || Ctrl->S.mode == WESSEL_BERCOVICI_1998_1D) Ctrl->S.mode += (dimension - 1);	
	if (Ctrl->S.mode == MITASOVA_MITAS_1993_2D ) Ctrl->S.mode += (dimension - 2);	
	if (Ctrl->S.mode == PARKER_1994 || Ctrl->S.mode == WESSEL_BECKER_2008) Ctrl->D.mode = 4;	/* Automatically set */

	GMT_io.in_col_type[GMT_X] = GMT_IS_LON;
	GMT_io.in_col_type[GMT_Y] = GMT_IS_LAT;

	Ctrl->D.mode--;	/* Since I added 0 to be 1-D later so no it is -1 */
	switch (Ctrl->D.mode) {	/* Set pointers to distance functions */
		case -1:	/* Cartesian 1-D x data */
			GMT_distance_func = NULL;
			GMT_azimuth_func = NULL;
			GMT_io.in_col_type[GMT_X] = GMT_IS_FLOAT;
			GMT_io.in_col_type[GMT_Y] = GMT_IS_FLOAT;
			normalize = 2;
			break;
		case 0:	/* Cartesian 2-D x,y data */
			GMT_distance_func = GMT_cartesian_dist;
			GMT_azimuth_func = GMT_az_backaz_cartesian;
			GMT_io.in_col_type[GMT_X] = GMT_IS_FLOAT;
			GMT_io.in_col_type[GMT_Y] = GMT_IS_FLOAT;
			normalize = 2;
			break;
		case 1:	/* 2-D lon, lat data, but scale to Cartesian flat earth */
			GMT_distance_func = GMT_flatearth_dist_km;
			GMT_azimuth_func = GMT_az_backaz_flatearth;
			normalize = 2;
			break;
		case 2:	/* 2-D lon, lat data, use spherical distances (geodesic if ELLIPSOID is nor sphere) */
			GMT_distance_func = (PFD)((GMT_IS_SPHERICAL) ? GMT_geodesic_dist_km : GMT_great_circle_dist_km);
			GMT_azimuth_func = (GMT_IS_SPHERICAL) ? GMT_az_backaz_sphere : GMT_az_backaz_geodesic;
			break;
		case 3:	/* 2-D lon, lat data, and Green's function needs cosine of spherical or geodesic distance */
			GMT_distance_func = (PFD)((GMT_IS_SPHERICAL) ? GMT_geodesic_dist_cos : GMT_great_circle_dist_cos);
			GMT_azimuth_func = (GMT_IS_SPHERICAL) ? GMT_az_backaz_sphere : GMT_az_backaz_geodesic;
			break;
		case 4:	/* 3-D Cartesian x,y,z data */
			GMT_distance_func = NULL;
			GMT_azimuth_func = NULL;
			break;
		default:	/* Cannot happen unless we make a bug */
			fprintf (stderr, "%s: BUG since D (=%ld) cannot be outside 0-5 range\n", GMT_program, Ctrl->D.mode+1);
			break;
	}

	if (Ctrl->D.mode <= 1 && Ctrl->L.active)
		normalize = 1;	/* Do not de-plane, just remove mean and normalize */
	else if (Ctrl->D.mode > 1 && Ctrl->L.active)
		fprintf (stderr, "%s: WARNING: -L ignored for -D modes 2 and 3\n", GMT_program);
	
	if (Ctrl->Q.active && dimension == 2) sincosd (Ctrl->Q.az, &Ctrl->Q.dir[GMT_X], &Ctrl->Q.dir[GMT_Y]);
	
	/* Now we are ready to take on some input values */

	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;
	done = FALSE;
	n_args = (argc > 1) ? argc : 2;
	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : dimension + 1;

	n_alloc = GMT_CHUNK;
	X = (double **) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double *), GMT_program);
	for (k = 0; k < n_alloc; k++) X[k] = (double *) GMT_memory (VNULL, dimension, sizeof (double), GMT_program);
	obs = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);

	n = m = n_read = 0;
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

		if (!nofile && gmtdefs.verbose) fprintf (stderr, "%s: Reading file %s\n", GMT_program, argv[fno]);

		if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) GMT_fgets (line, BUFSIZ, fp);

		while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

			n_read++;

			if (GMT_io.status & GMT_IO_MISMATCH) {
				fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
				continue;
			}

			for (k = 0; k < dimension; k++) X[n][k] = in[k];
			obs[n] = in[dimension];
			n++;

			if (n == n_alloc) {	/* Get more memory */
				old_n_alloc = n_alloc;
				n_alloc <<= 1;
				X = (double **) GMT_memory ((void *)X, n_alloc, sizeof (double *), GMT_program);
				for (k = (int)old_n_alloc; k < (int)n_alloc; k++) X[k] = (double *) GMT_memory (VNULL, dimension, sizeof (double), GMT_program);
				obs = (double *) GMT_memory ((void *)obs, n_alloc, sizeof (double), GMT_program);
			}
		}

		if (fp != GMT_stdin) GMT_fclose (fp);
	}

	for (k = n; k < (int)n_alloc; k++) GMT_free ((void *)X[k]);	/* Remove what was not used */
	X = (double **) GMT_memory ((void *)X, n, sizeof (double *), GMT_program);
	obs = (double *) GMT_memory ((void *)obs, n, sizeof (double), GMT_program);
	nm = n;

	if (Ctrl->A.active) {	/* Read gradient constraints from file */
		if (GMT_io.binary[GMT_IN]) GMT_io.ncol[GMT_IN]++;	/* Must assume it is just one extra column */
		GMT_import_table ((void *)Ctrl->A.file, GMT_IS_FILE, &S, 0.0, FALSE, FALSE, FALSE);
		m = S->n_records;	/* Total number of gradient constraints */
		nm += m;	/* New total of linear equations to solve */
		X = (double **) GMT_memory ((void *)X, (size_t)nm, sizeof (double *), GMT_program);
		for (k = n; k < nm; k++) X[k] = (double *) GMT_memory (VNULL, (size_t)dimension, sizeof (double), GMT_program);
		obs = (double *) GMT_memory ((void *)obs, (size_t)nm, sizeof (double), GMT_program);
		D = (double **) GMT_memory (VNULL, (size_t)m, sizeof (double *), GMT_program);
		for (k = 0; k < m; k++) D[k] = (double *) GMT_memory (VNULL, (size_t)dimension, sizeof (double), GMT_program);
		for (i = k = 0, p = n; i < S->n_segments; i++) {
			for (j = 0; j < S->segment[i]->n_rows; j++, k++, p++) {
				for (ii = 0; ii < dimension; ii++) X[p][ii] = S->segment[i]->coord[ii][j];
				switch (dimension) {
					case 1:	/* 1-D */
						switch (Ctrl->A.mode) {
							case 0:	/* x, slope */
								D[k][0] = 1.0;	/* Dummy since there is no direction for 1-D spline (the gradient is in the x-y plane) */
								obs[p] = S->segment[i]->coord[dimension][j];
								break;
							default:
								fprintf (stderr, "%s: Bad gradient mode selected for 1-D data (%ld) - aborting!\n", GMT_program, Ctrl->A.mode);
								exit (EXIT_FAILURE);
								break;
						}
						break;
					case 2:	/* 2-D */
						switch (Ctrl->A.mode) {
							case 1:	/* (x, y, az, gradient) */
								az = D2R * S->segment[i]->coord[2][j];
								obs[p] = S->segment[i]->coord[3][j];
								break;
							case 2:	/* (x, y, gradient, azimuth) */
								az = D2R * S->segment[i]->coord[3][j];
								obs[p] = S->segment[i]->coord[2][j];
								break;
							case 3:	/* (x, y, direction, gradient) */
								az = M_PI_2 - D2R * S->segment[i]->coord[2][j];
								obs[p] = S->segment[i]->coord[3][j];
								break;
							case 4:	/* (x, y, gx, gy) */
								az = atan2 (S->segment[i]->coord[2][j], S->segment[i]->coord[3][j]);		/* Get azimuth of gradient */
								obs[p] = hypot (S->segment[i]->coord[3][j], S->segment[i]->coord[3][j]);	/* Get magnitude of gradient */
								break;
							case 5:	/* (x, y, nx, ny, gradient) */
								az = atan2 (S->segment[i]->coord[2][j], S->segment[i]->coord[3][j]);		/* Get azimuth of gradient */
								obs[p] = S->segment[i]->coord[4][j];	/* Magnitude of gradient */
								break;
							default:
								fprintf (stderr, "%s: Bad gradient mode selected for 2-D data (%ld) - aborting!\n", GMT_program, Ctrl->A.mode);
								exit (EXIT_FAILURE);
								break;
						}
						sincos (az, &D[k][GMT_X], &D[k][GMT_Y]);
						break;
					case 3:	/* 3-D */
						switch (Ctrl->A.mode) {
							case 4:	/* (x, y, z, gx, gy, gz) */
								for (ii = 0; ii < 3; ii++) D[k][ii] = S->segment[i]->coord[3+ii][j];	/* Get the gradient vector */
								obs[p] = GMT_mag3v (D[k]);	/* This is the gradient magnitude */
								GMT_normalize3v (D[k]);		/* These are the direction cosines of the gradient */
								break;
							case 5: /* (x, y, z, nx, ny, nz, gradient) */
								for (ii = 0; ii < 3; ii++) D[k][ii] = S->segment[i]->coord[3+ii][j];	/* Get the unit vector */
								obs[p] = S->segment[i]->coord[6][j];	/* This is the gradient magnitude */
								break;
							default:
								fprintf (stderr, "%s: Bad gradient mode selected for 3-D data (%ld) - aborting!\n", GMT_program, Ctrl->A.mode);
								exit (EXIT_FAILURE);
								break;
						}
						break;
					default:
						fprintf (stderr, "%s: Bad dimension selected (%ld) - aborting!\n", GMT_program, dimension);
						exit (EXIT_FAILURE);
						break;
				}
			}
		}
		GMT_free_table (S);
	}
	
	if (m > 0 && normalize > 1) {
		normalize &= 1;	/* Only allow taking out data mean for mixed z/slope data */
		fprintf (stderr, "%s: Only remove/restore mean z in mixed {z, grad(z)} data sets\n", GMT_program);
	}
	
	if (gmtdefs.verbose) {
		if (m == 0)
			fprintf (stderr, "%s: Found %ld data points, yielding a %ld by %ld set of linear equations\n", GMT_program, n, nm, nm);
		else
			fprintf (stderr, "%s: Found %ld data points and %ld gradients, yielding a %ld by %ld set of linear equations\n",
				GMT_program, n, m, nm, nm);
	}
		
	if (Ctrl->T.file) {	/* Will have zeros and NaNs, only */
		if (GMT_read_grd_info (Ctrl->T.file, &h)) {
			fprintf (stderr, "%s: Error opening file %s\n", GMT_program, Ctrl->T.file);
			exit (EXIT_FAILURE);
		}
		nxy = n_ok = GMT_get_nm (h.nx, h.ny);
		w = (float *) GMT_memory (NULL, nxy, sizeof (float), GMT_program);
		if (GMT_read_grd (Ctrl->T.file, &h, w, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE)) {;
			fprintf (stderr, "%s: Error reading %s\n", GMT_program, Ctrl->T.file);
			exit (EXIT_FAILURE);
		}
		GMT_grd_init (&h, argc, argv, TRUE);
		for (ij = 0; ij < nxy; ij++) if (GMT_is_fnan (w[ij])) n_ok--;
	}
	else if (Ctrl->N.active) {	/* Read output locations from file */
		GMT_import_table ((void *)Ctrl->N.file, GMT_IS_FILE, &T, 0.0, FALSE, FALSE, FALSE);
	}
	else {	/* Fill in an equidistant output table or grid */
		h.node_offset = (int)Ctrl->F.active;
		h.x_inc = Ctrl->I.inc[GMT_X];
		Z.nz = h.ny = 1;	/* So that output logic will work for lower dimensions */
		if (dimension > 1) {
			h.y_inc = Ctrl->I.inc[GMT_Y];
			GMT_RI_prepare (&h);	/* Ensure -R -I consistency and set nx, ny */
			GMT_err_fail (GMT_grd_RI_verify (&h, 1), Ctrl->G.file);
			if (dimension == 3) {	/* Also set nz */
				Z.z_inc = Ctrl->I.inc[GMT_Z];
				Z.nz = GMT_get_n (Z.z_min, Z.z_max, Z.z_inc, h.node_offset);
			}
		}
		else
			h.nx = irint ((h.x_max - h.x_min) / h.x_inc) + h.node_offset;
		nxy = n_ok = h.nx * h.ny * Z.nz;
		if (dimension == 2) w = (float *) GMT_memory (NULL, nxy, sizeof (float), GMT_program);
	}

	switch (Ctrl->S.mode) {	/* Assing pointers to Green's functions and the gradient and set up required parameters */
		case SANDWELL_1987_1D:
			G = spline1d_sandwell;
			dGdr = gradspline1d_sandwell;
			break;
		case SANDWELL_1987_2D:
			G = spline2d_sandwell;
			dGdr = gradspline2d_sandwell;
			break;
		case SANDWELL_1987_3D:
			G = spline3d_sandwell;
			dGdr = gradspline3d_sandwell;
			break;
		case WESSEL_BERCOVICI_1998_1D:
			if (Ctrl->S.value[1] == 0.0 && h.x_inc > 0.0) Ctrl->S.value[1] = h.x_inc;
			if (Ctrl->S.value[1] == 0.0) Ctrl->S.value[1] = 1.0;
			par[0] = sqrt (Ctrl->S.value[0] / (1.0 - Ctrl->S.value[0])) / Ctrl->S.value[1];
			par[1] = 2.0 / par[0];
			G = spline1d_Wessel_Bercovici;
			dGdr = gradspline1d_Wessel_Bercovici;
			break;
		case WESSEL_BERCOVICI_1998_2D:
			if (Ctrl->S.value[1] == 0.0 && h.x_inc > 0.0) Ctrl->S.value[1] = 0.5 * (h.x_inc + h.y_inc);
			if (Ctrl->S.value[1] == 0.0) Ctrl->S.value[1] = 1.0;
			par[0] = sqrt (Ctrl->S.value[0] / (1.0 - Ctrl->S.value[0])) / Ctrl->S.value[1];
			par[1] = 2.0 / par[0];
			G = spline2d_Wessel_Bercovici;
			dGdr = gradspline2d_Wessel_Bercovici;
			break;
		case WESSEL_BERCOVICI_1998_3D:
			if (Ctrl->S.value[1] == 0.0 && h.x_inc > 0.0) Ctrl->S.value[1] = (h.x_inc + h.y_inc + Z.z_inc) / 3.0;
			if (Ctrl->S.value[1] == 0.0) Ctrl->S.value[1] = 1.0;
			par[0] = sqrt (Ctrl->S.value[0] / (1.0 - Ctrl->S.value[0])) / Ctrl->S.value[1];
			par[1] = 2.0 / par[0];
			G = spline3d_Wessel_Bercovici;
			dGdr = gradspline3d_Wessel_Bercovici;
			break;
		case MITASOVA_MITAS_1993_2D:
			/* par[0] = Ctrl->S.value[0]; */
			if (Ctrl->S.value[1] == 0.0) Ctrl->S.value[1] = 1.0;
			p_val = sqrt (Ctrl->S.value[0] / (1.0 - Ctrl->S.value[0])) / Ctrl->S.value[1];
			if (gmtdefs.verbose) fprintf (stderr, "p_val = %g\n", p_val);
			par[0] = p_val;
			par[1] = 0.25 * par[0] * par[0];
			G = spline2d_Mitasova_Mitas;
			dGdr = gradspline2d_Mitasova_Mitas;
			break;
		case MITASOVA_MITAS_1993_3D:
			/* par[0] = Ctrl->S.value[0]; */
			if (Ctrl->S.value[1] == 0.0) Ctrl->S.value[1] = 1.0;
			p_val = sqrt (Ctrl->S.value[0] / (1.0 - Ctrl->S.value[0])) / Ctrl->S.value[1];
			if (gmtdefs.verbose) fprintf (stderr, "p_val = %g\n", p_val);
			par[0] = p_val;
			par[1] = 0.25 * par[0] * par[0];
			G = spline3d_Mitasova_Mitas;
			dGdr = gradspline3d_Mitasova_Mitas;
			break;
		case PARKER_1994:
			par[0] = 6.0 / (M_PI*M_PI);
			G = spline2d_Parker;
			dGdr = gradspline2d_Parker;
			if (TEST) x0 = -1.0, x1 = 1.0;
			break;
		case WESSEL_BECKER_2008:
			par[0] = -0.5;
			p_val = sqrt (Ctrl->S.value[0] / (1.0 - Ctrl->S.value[0]));
			if (p_val <= 0.5) {	/* nu is real */
				double z[2];
				par[0] += sqrt (0.25 - p_val * p_val);
				par[1] = 0.0;
				par[4] = sin (M_PI * par[0]);
				z[0] = par[0] + 1.0;	z[1] = 0.0;
				par[3] = (M_PI / tan (M_PI * par[0])) - M_LOG_2 + 2.0 * (M_GAMMA + GMT_psi (z, NULL));
			}
			else {	/* nu is complex */
				double z[2], cot_piv[2], psi[2];
				par[1] = sqrt (p_val * p_val - 0.25);
				par[4] = -cosh (M_PI * par[1]);
				z[0] = par[0] * M_PI;	z[1] = par[1] * M_PI;
				Ccot (z, cot_piv);
				cot_piv[0] *= M_PI;	cot_piv[1] *= M_PI;
				z[0] = par[0] + 1.0;	z[1] = par[1];
				(void) GMT_psi (z, psi);
				psi[0] += M_GAMMA;
				psi[0] *= 2.0;	psi[1] *= 2.0;
				z[0] = cot_piv[0] + psi[0] - M_LOG_2;	z[1] = cot_piv[1] + psi[1];
				par[3] = z[0];	/* Ignore complex parts which cancels out */
			}
			par[2] = (M_PI / par[4]) - M_LOG_2;
			par[6] = 1.0 / (par[3] - par[2]);
			if (Ctrl->S.fast) {
				GMT_LONG nx;
				par[7] = Ctrl->S.value[1];
				par[8] = (r_max - r_min) / (par[7] - 1.0);
				par[9] = 1.0 / par[8];
				par[10] = r_min;
				nx = (GMT_LONG)irint (par[7]);
				if (gmtdefs.verbose) fprintf (stderr, "Precalculate -SQ lookup table with %ld items from %g to %g...", nx, r_min, r_max);
				WB_z = (double *) GMT_memory (VNULL, nx, sizeof (double), GMT_program);
				if (Ctrl->A.active) WB_g = (double *) GMT_memory (VNULL, nx, sizeof (double), GMT_program);
				spline2d_Wessel_Becker_init (par, WB_z, WB_g, Ctrl->A.active);
				if (gmtdefs.verbose) fprintf (stderr, "done\n");
				G = spline2d_Wessel_Becker_lookup;
				dGdr = gradspline2d_Wessel_Becker_lookup;
			}
			else {
				G = spline2d_Wessel_Becker;
				dGdr = gradspline2d_Wessel_Becker;
			}
			if (TEST) x0 = -1.0, x1 = 1.0;
			break;
	}

#ifdef DEBUG
	if (TEST) {	/* Just dump green's functions and exit */
		fprintf (stderr, "Warning: greenspline running in TEST mode for %s\n", method[Ctrl->S.mode]);
		printf ("# %s\n", method[Ctrl->S.mode]);
		dump_green (G, dGdr, par, x0, x1, 10001, WB_z, WB_g);
		GMT_free ((void *)Ctrl);
		if (dimension == 2) GMT_free ((void *)w);
		GMT_end (argc, argv);

		exit (EXIT_SUCCESS);
	}
#endif

	/* Remove mean (or LS plane) from data (we will add it back later) */

	do_normalization (X, obs, n, normalize, norm);
		
	/* Set up linear system Ax = z */
	
	mem = ((double)nm * (double)nm * (double)sizeof (REAL)) / 1024.0;	/* In kb */
	unit = 0;
	if (mem > 1024.0) {	/* Convert to Mb */
		mem /= 1024.0;
		unit = 1;
	}
	if (mem > 1024.0) {	/* Convert to Gb */
		mem /= 1024.0;
		unit = 2;
	}
	if (gmtdefs.verbose) fprintf (stderr, "%s: Square matrix requires %.1f %s\n", GMT_program, mem, mem_unit[unit]);
	A = (REAL *) GMT_memory (VNULL, nm * nm, sizeof (REAL), GMT_program);
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Build linear system using %s\n", GMT_program, method[Ctrl->S.mode]);

	r_min = DBL_MAX;	r_max = -DBL_MAX;
	for (j = 0; j < nm; j++) {	/* For each value or slope constraint */
		for (i = j; i < nm; i++) {
			ij = GMT_IJ (j, i, nm);
			ji = GMT_IJ (i, j, nm);
			r = get_radius (X[i], X[j], dimension);
			if (r < r_min) r_min = r;
			if (r > r_max) r_max = r;
			if (j < n) {	/* Value constraint */
				A[ij] = G (r, par, WB_z);
				if (ij == ji) continue;	/* Do the diagonal terms only once */
				if (i < n)
					A[ji] = A[ij];
				else {
					/* Get D, the directional cosine between the two points */
					/* Then get C = GMT_dot3v (D, dataD); */
					/* A[ji] = dGdr (r, par, WB_g) * C; */
					C = get_dircosine (D[i-n], X[i], X[j], dimension, FALSE);
					grad = dGdr (r, par, WB_g);
					A[ji] = grad * C;
				}
			}
			else if (i > n) {	/* Remaining gradient constraints */
				if (ij == ji) continue;	/* Diagonal gradient terms are zero */
				C = get_dircosine (D[j-n], X[i], X[j], dimension, TRUE);
				grad = dGdr (r, par, WB_g);
				A[ij] = grad * C;
				C = get_dircosine (D[i-n], X[i], X[j], dimension, FALSE);
				A[ji] = grad * C;
			}
		}
	}

	if (Ctrl->C.active) {		/* Solve using svd decomposition */
		GMT_LONG n_use;
		REAL *v, *s, *b, eig_max = 0.0;
		char format[GMT_LONG_TEXT];
		void svdcmp (REAL *a, GMT_LONG m, GMT_LONG n, REAL *w, REAL *v);
		GMT_LONG solve_svd (REAL *u, GMT_LONG m, GMT_LONG n, REAL *v, REAL *w, REAL *b, GMT_LONG k, REAL *x, double cutoff);
		
		if (gmtdefs.verbose) fprintf (stderr, "%s: Solve linear equations by SVD ", GMT_program);
	
		v = (REAL *) GMT_memory (VNULL, nm * nm, sizeof (REAL), GMT_program);
		s = (REAL *) GMT_memory (VNULL, nm, sizeof (REAL), GMT_program);
		svdcmp (A, nm, nm, s, v);
		if (Ctrl->C.file) {	/* Save the eigen-values for study */
			REAL *eig;
		 	eig = (REAL *) GMT_memory (VNULL, nm, sizeof (REAL), GMT_program);
			memcpy ((void *)eig, (void *)s, (size_t)(nm * sizeof (REAL)));
			if (gmtdefs.verbose) fprintf (stderr, "Eigen-value rations s(i)/s(0) saved to %s\n", Ctrl->C.file);
			if ((fp = GMT_fopen (Ctrl->C.file, "w")) == NULL) {
				fprintf (stderr, "%s: Error creating file %s\n", GMT_program, Ctrl->C.file);
				exit (EXIT_FAILURE);
			}
			sprintf (format, "%%d%s%s\n", gmtdefs.field_delimiter, gmtdefs.d_format);
			/* Sort eigenvalues into ascending order */
			if (sizeof (REAL) == sizeof (double))
				qsort ((void *)eig, nm, sizeof (double), GMT_comp_double_asc);
			else
				qsort ((void *)eig, nm, sizeof (float), GMT_comp_float_asc);
			eig_max = eig[nm-1];
			for (i = 0, j = nm-1; i < nm; i++, j--) fprintf (fp, format, i, eig[j] / eig_max);
			GMT_fclose (fp);
			GMT_free ((void *)eig);
			if (Ctrl->C.value < 0.0) exit (EXIT_SUCCESS);
		}
		b = (REAL *) GMT_memory (VNULL, nm, sizeof (REAL), GMT_program);
		if (sizeof (REAL) == sizeof (double)) {
			memcpy ((void *)b, (void *)obs, (size_t)(nm * sizeof (REAL)));
			n_use = solve_svd (A, nm, nm, v, s, b, 1, obs, Ctrl->C.value);
		}
		else {	/* Must use temporary float array to capture result */
			REAL *z4;
			for (i = 0; i < nm; i++) b[i] = (REAL)obs[i];
			z4 = (REAL *) GMT_memory (VNULL, nm, sizeof (REAL), GMT_program);
			n_use = solve_svd (A, nm, nm, v, s, b, 1, z4, Ctrl->C.value);
			for (i = 0; i < nm; i++) obs[i] = (double)z4[i];
			GMT_free ((void *)z4);
		}
		if (gmtdefs.verbose) fprintf (stderr, "[%ld of %ld eigen-values used]\n", n_use, nm);
			
		GMT_free ((void *)s);
		GMT_free ((void *)v);
		GMT_free ((void *)b);
	}
	else {				/* Gauss-Jordan elimination */
		GMT_LONG gaussj (REAL *a, GMT_LONG n, GMT_LONG ndim, REAL *b, GMT_LONG m, GMT_LONG mdim);
		if (gmtdefs.verbose) fprintf (stderr, "%s: Solve linear equations by Gauss-Jordan elimination\n", GMT_program);
		if (sizeof (REAL) == sizeof (double))
			gaussj (A, nm, nm, obs, 1, 1);
		else {
			REAL *z4;
			z4 = (REAL *) GMT_memory (VNULL, nm, sizeof (REAL), GMT_program);
			for (i = 0; i < nm; i++) z4[i] = (REAL)obs[i];
			gaussj (A, nm, nm, z4, 1, 1);
			for (i = 0; i < nm; i++) obs[i] = (double)z4[i];
			GMT_free ((void *)z4);
		}
	}
	alpha = obs;	/* Just a different name since the obs vector now holds the alpha factors */
		
	GMT_free ((void *)A);

	if (Ctrl->N.file) {	/* Specified nodes only */
		GMT_LONG n_output;
		double out[4];
		FILE *fp;
	
		if (!Ctrl->G.active) {
			fp = GMT_stdout;
#ifdef SET_IO_MODE
			GMT_setmode (GMT_OUT);
#endif
		}
		else if ((fp = GMT_fopen (Ctrl->G.file, GMT_io.w_mode)) == NULL) {
			fprintf (stderr, "%s: Error creating file %s\n", GMT_program, Ctrl->G.file);
			exit (EXIT_FAILURE);
		}
	
		n_output = dimension + 1;
		if (gmtdefs.verbose) fprintf (stderr, "%s: Evaluate spline at %ld given locations\n", GMT_program, T->n_records);
		for (i = 0; i < T->n_segments; i++) {
			for (j = 0; j < T->segment[i]->n_rows; j++) {
				for (ii = 0; ii < dimension; ii++) out[ii] = T->segment[i]->coord[ii][j];
				out[dimension] = 0.0;
				for (p = 0; p < nm; p++) {
					r = get_radius (out, X[p], dimension);
					if (r < r_min) r_min = r;
					if (r > r_max) r_max = r;
					if (Ctrl->Q.active) {
						C = get_dircosine (Ctrl->Q.dir, out, X[p], dimension, FALSE);
						part = dGdr (r, par, WB_z) * C;
					}
					else
						part = G (r, par, WB_z);
					out[dimension] += alpha[p] * part;
				}
				out[dimension] = undo_normalization (out, out[dimension], normalize, norm);
				GMT_output (fp, n_output, out);
			}
		}
		GMT_free_table (T);
		if (fp != GMT_stdout) GMT_fclose (fp);
	}
	else {
		GMT_LONG nz_off, nxy, n_output;
		double *xp, *yp = NULL, wp, V[4];
		FILE *fp = NULL;
		if (gmtdefs.verbose) fprintf (stderr, "%s: Evaluate spline at %ld equidistant output locations\n", GMT_program, n_ok);
		/* Precalculate coordinates */
		xp = (double *) GMT_memory (NULL, h.nx, sizeof (double), GMT_program);
		if (dimension > 1) yp = (double *) GMT_memory (NULL, h.ny, sizeof (double), GMT_program);
		for (i = 0; i < h.nx; i++) xp[i] = GMT_i_to_x (i, h.x_min, h.x_max, h.x_inc, h.xy_off, h.nx);
		if (dimension > 1) for (j = 0; j < h.ny; j++) yp[j] = GMT_j_to_y (j, h.y_min, h.y_max, h.y_inc, h.xy_off, h.ny);
		nxy = h.nx * h.ny;
		n_output = dimension + 1;
		if (dimension != 2) {	/* Write ascii table to named file or stdout */
			fp = (Ctrl->G.active) ? GMT_fopen (Ctrl->G.file, "w") : GMT_stdout;
		}
		for (k = nz_off = 0; k < Z.nz; k++, nz_off += nxy) {
			if (dimension == 3) V[GMT_Z] = GMT_i_to_x (k, Z.z_min, Z.z_max, Z.z_inc, h.xy_off, Z.nz);
			for (j = 0; j < h.ny; j++) {
				if (dimension > 1) V[GMT_Y] = yp[j];
				for (i = 0; i < h.nx; i++) {
					ij = GMT_IJ (j, i, h.nx) + nz_off;
					if (dimension == 2 && GMT_is_fnan (w[ij])) continue;	/* Only do solution where mask is not NaN */
					V[GMT_X] = xp[i];
					/* Here, V holds the current output coordinates */
					for (p = 0, wp = 0.0; p < nm; p++) {
						r = get_radius (V, X[p], dimension);
						if (r < r_min) r_min = r;
						if (r > r_max) r_max = r;
						if (Ctrl->Q.active) {
							C = get_dircosine (Ctrl->Q.dir, V, X[p], dimension, FALSE);
							part = dGdr (r, par, WB_z) * C;
						}
						else
							part = G (r, par, WB_z);
						wp += alpha[p] * part;
					}
					V[dimension] = (float)undo_normalization (V, wp, normalize, norm);
					if (dimension == 2)	/* Special 2-D grid output */
						w[ij] = (float)V[dimension];
					else	/* Crude dump for now for both 1-D and 3-D */
						GMT_output (fp, n_output, V);
				}
			}
		}
		if (dimension == 2) {
			sprintf (h.remark, "Method: %s (%s)", method[Ctrl->S.mode], Ctrl->S.arg);
			GMT_write_grd (Ctrl->G.file, &h, w, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE);
			GMT_free ((void *)w);
		}
		else if (Ctrl->G.active)
			GMT_fclose (fp);
			
		GMT_free ((void *)xp);
		if (dimension > 1) GMT_free ((void *)yp);
	}
	
	for (p = 0; p < nm; p++) GMT_free ((void *)X[p]);
	GMT_free ((void *)X);
	GMT_free ((void *)obs);
	if (m) {
		for (p = 0; p < m; p++) GMT_free ((void *)D[p]);
		GMT_free ((void *)D);
	}
	if (Ctrl->S.fast) {
		GMT_free ((void *)WB_z);
		if (Ctrl->A.active) GMT_free ((void *)WB_g);
	}
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Done [ R min/max = %g/%g]\n", GMT_program, r_min, r_max);

	GMT_free ((void *)Ctrl);
	GMT_end (argc, argv);
	
	exit (EXIT_SUCCESS);
}

#ifdef DEBUG
/* Dump a table of x, G, dGdx for test purposes */
void dump_green (PFD G, PFD D, double par[], double x0, double x1, GMT_LONG N, double *zz, double *gg)
{
	GMT_LONG i;
	double x, dx, dy, y, t, ry, rdy;
	double min_y, max_y, min_dy, max_dy;
	
	min_y = min_dy = DBL_MAX;
	max_y = max_dy = -DBL_MAX;
	
	dx = (x1 - x0) / (N - 1);
	for (i = 0; i < N; i++) {
		x = x0 + i * dx;
		t = (x0 < 0.0) ? acosd (x) : x;
		y = G(x, par, zz);
		dy = D(x, par, gg);
		if (y < min_y) min_y = y;
		if (y > max_y) max_y = y;
		if (dy < min_dy) min_dy = dy;
		if (dy > max_dy) max_dy = dy;
	}
	ry = max_y - min_y;
	rdy = max_dy - min_dy;
	for (i = 0; i < N; i++) {
		x = x0 + i * dx;
		t = (x0 < 0.0) ? acosd (x) : x;
		y = G(x, par, zz);
		dy = D(x, par, gg);
		dy = (rdy > 0.0) ? (dy - min_dy)/rdy : 1.0;
		printf ("%g\t%g\t%g\t%g\n", x, (y - min_y) / ry, dy, t);
	}
}
#endif

/*----------------------  ONE DIMENSION ---------------------- */
/* spline1d_sandwell computes the Green function for a 1-d spline
 * as per Sandwell [1987], G(r) = r^3.  All r must be >= 0.
 */

double spline1d_sandwell (double r, double par[], double *G)
{
	if (r == 0.0) return (0.0);

	return (pow (r, 3.0));	/* Just regular spline; par not used */
}

double gradspline1d_sandwell (double r, double par[], double *G)
{
	return (r);	/* Just regular spline; par not used */
}

double spline1d_Wessel_Bercovici (double r, double par[], double *G)
/* spline1d_Wessel_Bercovici computes the Green function for a 1-d spline
 * in tension as per Wessel and Bercovici [1988], G(u) = exp(-u) + u - 1,
 * where u = par[0] * r and par[0] = sqrt (t/(1-t)).
 * All r must be >= 0. par[0] = c
 */
{
	double cx;

	if (r == 0.0) return (0.0);

	cx = par[0] * r;
	return (exp (-cx) + cx - 1.0);
}

double gradspline1d_Wessel_Bercovici (double r, double par[], double *G)
{
	double cx;

	if (r == 0.0) return (0.0);

	cx = par[0] * r;
	return (1.0 - exp (-cx));
}

/*----------------------  TWO DIMENSIONS ---------------------- */

double spline2d_sandwell (double r, double par[], double *G)
{
	if (r == 0.0) return (0.0);

	return (r * r * (log (r) - 1.0));	/* Just regular spline; par not used */
}

double gradspline2d_sandwell (double r, double par[], double *G)
{
	if (r == 0.0) return (0.0);

	return (r * (2.0 * log (r) - 1.0));	/* Just regular spline; par not used */
}

/* spline2d_Wessel_Bercovici computes the Green function for a 2-d spline
 * in tension as per Wessel and Bercovici [1988], G(u) = K(u) - log(u),
 * where u = par[0] * r and par[0] = sqrt (t/(1-t)).
 * K is the modified Bessel function K of order zero.
 * All r must be >= 0.
 * par[0] = c
 * par[1] = 2/c
 */

double spline2d_Wessel_Bercovici (double r, double par[], double *G)
{
	double y, z, cx, t, g;

	if (r == 0.0) return (0.0);

	cx = par[0] * r;
	if (r <= par[1]) {
		y = 0.25 * (t = cx * cx);
		z = t / 14.0625;
		g = (-log(0.5*cx) * (z * (3.5156229 + z * (3.0899424 + z * (1.2067492 + z * (0.2659732 + 
			z * (0.360768e-1 + z * 0.45813e-2))))))) + (y * (0.42278420 + y * (0.23069756 + 
			y * (0.3488590e-1 + y * (0.262698e-2 + y * (0.10750e-3 + y * 0.74e-5))))));
	}
	else {
		y = par[1] / r;
		g = (exp (-cx) / sqrt (cx)) * (1.25331414 + y * (-0.7832358e-1 + y * (0.2189568e-1 + 
			y * (-0.1062446e-1 + y * (0.587872e-2 + y * (-0.251540e-2 + y * 0.53208e-3))))))
			+ log (cx) - M_LOG_2 + M_GAMMA;
	}
	return (g);
}

double gradspline2d_Wessel_Bercovici (double r, double par[], double *G)
{
	double y, z, cx, t, dgdr;

	if (r == 0.0) return (0.0);

	cx = par[0] * r;
	if (r <= par[1]) {
		y = 0.25 * (t = cx * cx);
		z = t / 14.0625;
		dgdr = -((log(0.5*cx) * (cx * (0.5 + z * (0.87890594 + z * (0.51498869 + z * (0.15084934 + 
			z * (0.2658733e-1 + z * (0.301532e-2 + z * 0.32411e-3)))))))) + (1.0/cx) * (y * (0.15443144 +
			y * (-0.67278579 + y * (-0.18156897 + y * (-0.1919402e-1 + y * (-0.110404e-2 + y * (-0.4686e-4))))))));
	}
	else {
		y = par[1] / r;
		dgdr = 0.5*y - ((exp (-cx) / sqrt (cx)) * (1.25331414 + y * (0.23498619 + y * (-0.3655620e-1 +
			y * (0.1504268e-1 + y * (-0.780353e-2 + y * (0.325614e-2 + y * (-0.68245e-3))))))));
	}
	return (dgdr*par[0]);
}

/* spline2d_Mitasova_Mitas computes the regularized Green function for a 2-d spline
 * in tension as per Mitasova and Mitas [1993], G(u) = log (u) + Ei(u),
 * where u = par[1] * r^2. All r must be >= 0.
 * par[0] = phi (par[0] unused by function)
 * par[1] = phi^2/4
 */
 
double spline2d_Mitasova_Mitas (double r, double par[], double *G)
{
	double x, z, g, En, Ed;
	
	if (r == 0.0) return (0.0);
	
	x = z = par[1] * r * r;
	if (z <= 1.0) {
		g = 0.99999193 * x;
		x *= x;
		g -= 0.24991055 * x;
		x *= x;
		g += 0.05519968 * x;
		x *= x;
		g -= 0.00976004 * x;
		x *= x;
		g += 0.00107857 * x;
	}
	else {
		g = log(x) + M_GAMMA;
		En = 0.2677737343 +  8.6347608925 * x;
		Ed = 3.9584869228 + 21.0996530827 * x;
		x *= x;
		En += 18.0590169730 * x;
		Ed += 25.6329561486 * x;
		x *= x;
		En += 8.5733287401 * x;
		Ed += 9.5733223454 * x;
		x *= x;
		En += x;
		Ed += x;
		g += (En / Ed) / (z * exp(z));
	}
	return (g);
}

double gradspline2d_Mitasova_Mitas (double r, double par[], double *G)
{
	double u, dgdr;
	
	if (r == 0.0) return (0.0);

	u = par[1] * r * r;
	dgdr = 2.0 * (1.0 - exp(-u))/r;
	return (dgdr);
}

/*----------------------  TWO DIMENSIONS (SPHERE) ---------------------- */

/* spline2d_Parker computes the  Green function for a 2-d surface spline
 * as per Parker [1994], G(x) = dilog(),
 * where x is cosine of distances. All x must be -1 <= x <= +1.
 * Parameters passed are:
 * par[0] = 6/M_PI^2 (to normalize results)
 */
 
double spline2d_Parker (double x, double par[], double *G)
{	/* Normalized to 0-1 */
	if (x == +1.0) return (1.0);
	if (x == -1.0) return (0.0);
	return (GMT_dilog (0.5 - 0.5 * x) * par[0]);
}

double gradspline2d_Parker (double x, double par[], double *G)
{	/* Normalized to 0-1 */
	if (x == +1.0 || x == -1.0) return (0.0);
	return (log(0.5 - 0.5 * x) * sqrt ((1.0 - x) / (1.0 + x)));
}

/* spline2d_Wessel_Becker computes the Green function for a 2-d surface spline
 * in tension as per Wessel and Becker [2008], G(x) = M_PI * Pv(-x)/sin (v*x) - log (1-x),
 * where x is cosine of distances.  All x must be -1 <= x <= +1.
 * Parameters passed are:
 * par[0] = Real(nu)
 * par[1] = Imag(nu)
 * par[2] = G(-1)
 * par[3] = G(+1)
 * par[4] = Real {sin (nu * M_PI}
 * par[5] = Imag {sin (nu * M_PI)} == 0
 * par[6] = 1 / (par[3] - par[2])
 * par[7-9] is used by the lookup macinery
 */
 
double spline2d_Wessel_Becker (double x, double par[], double *G)
{	/* g = M_PI * Pv(-x)/sin (v*x) - log (1-x) normalized to 0-1 */
	GMT_LONG n;
	double z[2], pq[4];
	
	if (x >= (1.0-GMT_CONV_LIMIT)) return (1.0);
	if (x <= (GMT_CONV_LIMIT-1.0)) return (0.0);

	GMT_PvQv (-x, par, pq, &n);	/* Get P_nu(-x) */
	Cdiv (pq, &par[4], z);		/* Get P_nu(-x) / sin (nu*M_PI) */
	pq[0] = M_PI * z[0] - log (1.0 - x);
#ifdef DEBUG
	if (TEST) {
		pq[1] = M_PI * z[1];
		if (gmtdefs.verbose && fabs (pq[1]) > 1.0e-6) fprintf (stderr, "Im{G(%g)} = %g\n", x, pq[1]);
	}
#endif
	
	return ((pq[0] - par[2]) * par[6]);	/* Normalizes to yield 0-1 */
}

double gradspline2d_Wessel_Becker (double x, double par[], double *G)
{	/* g = -M_PI * (v+1)*[x*Pv(-x)+Pv+1(-x)]/(sin (v*x)*sin(theta)) - 1/(1-x), normalized to 0-1 */
	GMT_LONG n;
	double z[2], v1[2], pq[4], s;
	
	if (fabs (x) >= (1.0 - GMT_CONV_LIMIT)) return (0.0);

	GMT_PvQv (-x, par, pq, &n);			/* Get P_nu(-x) */
	z[0] = pq[0] * x;	z[1] = pq[1] * x;	/* Get x*P_nu(-x) */
	v1[0] = par[0] + 1.0;	v1[1] = par[1];		/* Get nu+1 */
	GMT_PvQv (-x, v1, pq, &n);			/* Get P_(nu+1)(-x) */
	z[0] += pq[0];	z[1] += pq[1];			/* Get x*P_nu(-x) + P_(nu+1)(-x) */
	Cdiv (z, &par[4], pq);				/* Get ---"--- / sin (nu*M_PI) */
	Cmul (pq, v1, z);				/* Mul by nu + 1 */
	s = M_PI / sqrt (1.0 - x*x);			/* Mul by pi/sin(theta) */
	z[0] *= s;
#ifdef DEBUG
	if (TEST) {
		z[1] *= s;
		if (gmtdefs.verbose && fabs (z[1]) > 1.0e-6) fprintf (stderr, "Im{G(%g)} = %g\n", x, z[1]);
	}
#endif
	z[0] += sqrt ((1.0 + x)/(1.0 - x));		/* Add in last term */
	
	return (-z[0]);
}

/* Given the lookup tables, this is how we use these functions
 * Here, par[7] is number of points in spline
 *	 par[8] is spline spacing dx
 *	 par[9] is inverse, 1/dx
 *	par[10] is min x
 */

double spline2d_lookup (double x, double par[], double *y)
{
	GMT_LONG k;
	double f, f0, df;
	
	f = (x - par[10]) * par[9];	/* Floating point index */
	f0 = floor (f);
	df = f - f0;
	k = (GMT_LONG)f0;
	if (df == 0.0) return (y[k]);
	return (y[k]*(1.0 - df) + y[k+1] * df);
}

double spline2d_Wessel_Becker_lookup (double x, double par[], double *y)
{
	return (spline2d_lookup (x, par, y));
}

double gradspline2d_Wessel_Becker_lookup (double x, double par[], double *y)
{
	return (spline2d_lookup (x, par, y));
}

void spline2d_Wessel_Becker_init (double par[], double *z, double *g, GMT_LONG grad)
{
	int i, nx;
	double x;
#ifdef DUMP
	FILE *fp;
	size_t n_out;
	double out[3];
	fp = fopen ("greenspline.b", "wb");
	n_out = (grad) ? 3 : 2;
#endif
	nx = (GMT_LONG)irint (par[7]);
	for (i = 0; i < nx; i++) {
		x = par[10] + i * par[8];
		z[i] = spline2d_Wessel_Becker (x, par, NULL);
		if (grad) g[i] = gradspline2d_Wessel_Becker (x, par, NULL);
#ifdef DUMP
		out[0] = x;	out[1] = z[i];	if (grad) out[2] = g[i];
		fwrite ((void *)out, sizeof (double), n_out, fp);
#endif
	}
#ifdef DUMP
	fclose (fp);
#endif
}

/*----------------------  THREE DIMENSIONS ---------------------- */
		
/* spline3d_sandwell computes the Green function for a 3-d spline
 * as per Sandwell [1987], G(r) = r.  All r must be >= 0.
 */

double spline3d_sandwell (double r, double par[], double *G)
{
	if (r == 0.0) return (0.0);

	return (r);	/* Just regular spline; par not used */
}

double gradspline3d_sandwell (double r, double par[], double *G)
{
	return (1.0);	/* Just regular spline; par not used */
}

double spline3d_Wessel_Bercovici (double r, double par[], double *G)
/* spline1d_Wessel_Bercovici computes the Green function for a 3-d spline
 * in tension as per Wessel and Bercovici [1988], G(u) = [exp(-u) -1]/u,
 * where u = par[0] * r and par[0] = sqrt (t/(1-t)).
 * All r must be >= 0. par[0] = c
 */
{
	double cx;

	if (r == 0.0) return (0.0);

	cx = par[0] * r;
	return (((exp (-cx) - 1.0) / cx) + 1.0);
}

double gradspline3d_Wessel_Bercovici (double r, double par[], double *G)
{
	double cx;

	if (r == 0.0) return (0.0);

	cx = par[0] * r;
	return ((1.0 - exp (-cx) * (cx + 1.0))/ (cx * r));
}

/* spline3d_Mitasova_Mitas computes the regularized Green function for a 3-d spline
 * in tension as per Mitasova and Mitas [1993], G(u) = erf (u/2)/u - 1/sqrt(pi),
 * where u = par[0] * r. All r must be >= 0. par[0] = phi
 */
 
double spline3d_Mitasova_Mitas (double r, double par[], double *G)
{
	double x;
	
	if (r == 0.0) return (0.0);
	
	x = par[0] * r;
	return ((erf (0.5 * x) / x) - M_INV_SQRT_PI);
}

double gradspline3d_Mitasova_Mitas (double r, double par[], double *G)
{
	double u, dgdr;
	
	if (r == 0.0) return (0.0);

	u = par[0] * r;
	dgdr = ((u/M_SQRT_PI)*exp (-u*u) - erf(0.5*u)) / (u * r);
	return (dgdr);
}

/* GENERAL NUMERICAL FUNCTIONS */

/* Modified from similar function in Numerical Recipes */

GMT_LONG gaussj (REAL *a, GMT_LONG n, GMT_LONG ndim, REAL *b, GMT_LONG m, GMT_LONG mdim)
{
	GMT_LONG i, j, k, l, ll, *ipiv, *indxc, *indxr, irow = 0, icol = 0;
	REAL big, dum, pivinv;
	
	ipiv  = (GMT_LONG *) GMT_memory (VNULL, (size_t)n, sizeof (GMT_LONG), "gaussj");
	indxc = (GMT_LONG *) GMT_memory (VNULL, (size_t)n, sizeof (GMT_LONG), "gaussj");
	indxr = (GMT_LONG *) GMT_memory (VNULL, (size_t)n, sizeof (GMT_LONG), "gaussj");
	
	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++) {
			if (ipiv[j] != 1) {
				for (k = 0; k < n; k++) {
					if (ipiv[k] == 0) {
						if ((dum = fabs(a[j*ndim+k])) >= big) {
							big = dum;
							irow = j;
							icol = k;
						}
					}
					else if (ipiv[k] > 1) {
						fprintf (stderr, "gaussj: Singular matrix!\n");
						GMT_free ((void *) ipiv);
						GMT_free ((void *) indxc);
						GMT_free ((void *) indxr);
						return (1);
					}
				}
			}
		}
		ipiv[icol]++;
		
		if (irow != icol) {
			for (l = 0; l < n; l++) {
				dum = a[irow*ndim+l];
				a[irow*ndim+l] = a[icol*ndim+l];
				a[icol*ndim+l] = dum;
			}
			for (l = 0; l < m; l++) {
				dum = b[irow*mdim+l];
				b[irow*mdim+l] = b[icol*mdim+l];
				b[icol*mdim+l] = dum;
			}
		}
		
		indxr[i] = irow;
		indxc[i] = icol;
		if (a[icol*ndim+icol] == 0.0) {
			fprintf (stderr, "gaussj: Singular matrix!\n");
			GMT_free ((void *) ipiv);
			GMT_free ((void *) indxc);
			GMT_free ((void *) indxr);
			return (1);
		}
		pivinv = 1.0 / a[icol*ndim+icol];
		a[icol*ndim+icol] = 1.0;
		for (l = 0; l < n; l++) a[icol*ndim+l] *= pivinv;
		for (l = 0; l < m; l++)  b[icol*mdim+l] *= pivinv;
		for (ll = 0; ll < n; ll++) {
			if (ll != icol) {
				dum = a[ll*ndim+icol];
				a[ll*ndim+icol] = 0.0;
				for (l = 0; l < n; l++) a[ll*ndim+l] -= a[icol*ndim+l] * dum;
				for (l = 0; l < m; l++) b[ll*mdim+l] -= b[icol*mdim+l] * dum;
			}
		}
	}
	for (l = n-1; l >= 0; l--) {
		if (indxr[l] != indxc[l]) {
			for (k = 0; k < n; k++) {
				dum = a[k*ndim+indxr[l]];
				a[k*ndim+indxr[l]] = a[k*ndim+indxc[l]];
				a[k*ndim+indxc[l]] = dum;
			}
		}
	}
	
	GMT_free ((void *) ipiv);
	GMT_free ((void *) indxc);
	GMT_free ((void *) indxr);
	
	return (0);
}

/* Given a matrix a[0..m-1][0...n-1], this routine computes its singular
	value decomposition, A=UWVt.  The matrix U replaces a on output.
	The diagonal matrix of singular values W is output as a vector
	w[0...n-1].  The matrix V (Not V transpose) is output as
	v[0...n-1][0....n-1].  m must be greater than or equal to n; if it is
	smaller, then a should be filled up to square with zero rows.
	
	Modified from Numerical Recipes -> page 68.
	
*/


#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void svdcmp (REAL *a, GMT_LONG m, GMT_LONG n, REAL *w, REAL *v)
{
	/* void svdcmp(REAL *a,int m,int n,REAL *w,REAL *v) */
	
	GMT_LONG flag,i,its,j,jj,k,l=0,nm=0;
	REAL c,f,h,s,x,y,z;
	REAL anorm=0.0,tnorm, g=0.0,scale=0.0;
	REAL *rv1;
	
	if(m < n)
	{
		fprintf(stderr,"Error in SVCMP: m < n augment A with additional rows\n");
		exit(-1);
	}
	
		/* allocate work space */
		
	rv1=(REAL *)calloc(n,sizeof(REAL));
	if(rv1 == NULL)
	{
		fprintf(stderr,"Error in SVCMP: Can't allocate work space\n");
		exit(-1);
	}
	
		/* do householder reduction to bidiagonal form */
		
	for(i=0;i<n;i++)
	{
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if(i < m)
		{
			for(k=i;k<m;k++)
				scale += fabs(a[k*n+i]);		/* a[k][i] */
			if(scale)
			{
				for(k=i;k<m;k++)
				{
					a[k*n+i] /= scale;	/* a[k][i] */
					s += a[k*n+i]*a[k*n+i];	/* a[k][i] */
				}
				f=a[i*n+i];	/* a[i][i] */
				g= -1.0*SIGN(sqrt(s),f);
				h=f*g-s;
				a[i*n+i]=f-g;	/* a[i][i] */
				if(i != n-1)
				{
					for(j=l;j<n;j++)
					{
						for(s=0.0,k=i;k<m;k++)
							s += a[k*n+i]*a[k*n+j];	/* a[k][i] a[k][j] */
						f=s/h;
						for(k=i;k<m;k++)
							a[k*n+j] += f*a[k*n+i];	/* a[k][j] a[k][i] */
					}
				}
				for(k=i;k<m;k++)
					a[k*n+i] *= scale;	/* a[k][i] */
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if(i <= m-1 && i != n-1)
		{
			for(k=l;k<n;k++)
				scale += fabs(a[i*n+k]);	/* a[i][k] */
			if(scale)
			{
				for(k=l;k<n;k++)
				{
					a[i*n+k] /= scale;	/* a[i][k] */
					s += a[i*n+k]*a[i*n+k];	/* a[i][k] */
				}
				f=a[i*n+l];	/* a[i][l] */
				g = -1.0*SIGN(sqrt(s),f);
				h=f*g-s;
				a[i*n+l]=f-g;	/* a[i][l] */
				for(k=l;k<n;k++)
					rv1[k]=a[i*n+k]/h;	/* a[i][k] */
				if(i != m-1)
				{
					for(j=l;j<m;j++)
					{
						for(s=0.0,k=l;k<n;k++)
							s += a[j*n+k]*a[i*n+k];	/*a[j][k] a[i][k] */
						for(k=l;k<n;k++)
							a[j*n+k] += s*rv1[k];	/* a[j][k] */
					}
				}
				for(k=l;k<n;k++)
					a[i*n+k] *= scale;	/* a[i][k] */
			}
		}
		tnorm=fabs(w[i])+fabs(rv1[i]);
		anorm=MAX(anorm,tnorm);
	}
						
	/* accumulation of right-hand transforms */
		
	for(i=n-1;i>=0;i--)
	{
		if(i < n-1)
		{
			if(g)
			{
				for(j=l;j<n;j++)
					v[j*n+i]=(a[i*n+j]/a[i*n+l])/g;	/* v[j][i] a[i][j] a[i][l] */
				for(j=l;j<n;j++)
				{
					for(s=0.0,k=l;k<n;k++)
						s += a[i*n+k]*v[k*n+j];	/* a[i][k] v[k][j] */
					for(k=l;k<n;k++)
						v[k*n+j] += s*v[k*n+i];	/* v[k][j] v[k][i] */
				}
			}
			for(j=l;j<n;j++)
				v[i*n+j]=v[j*n+i]=0.0;	/* v[i][j] v[j][i] */
		}
		v[i*n+i]=1.0;	/* v[i][i] */
		g=rv1[i];
		l=i;
	}
	
	/* accumulation of left-hand transforms */
		
	for(i=n-1;i>=0;i--)
	{
		l=i+1;
		g=w[i];
		if(i < n-1)
			for(j=l;j<n;j++)
				a[i*n+j]=0.0;	/* a[i][j] */
		if(g)
		{
			g=1.0/g;
			if(i != n-1)
			{
				for(j=l;j<n;j++)
				{
					for(s=0.0,k=l;k<m;k++)
						s += a[k*n+i]*a[k*n+j];	/* a[k][i] a[k][j] */
					f=(s/a[i*n+i])*g;	/* a[i][i] */
					for(k=i;k<m;k++)
						a[k*n+j] += f*a[k*n+i];	/* a[k][j] a[k][i] */
				}
			}
			for(j=i;j<m;j++)
				a[j*n+i] *= g;	/* a[j][i] */
		}
		else
		{
			for(j=i;j<m;j++)
				a[j*n+i]=0.0;	/* a[j][i] */
		}
		++a[i*n+i];	/* a[i][i] */
	}
	
	/* diagonalization of the bidiagonal form */
		
	for(k=n-1;k>=0;k--)			/* loop over singular values */
	{
		for(its=1;its<=30;its++)	/* loop over allowed iterations */
		{
			flag=1;
			for(l=k;l>=0;l--)		/* test for splitting */
			{
				nm=l-1;
				if(fabs(rv1[l])+anorm == anorm)
				{
					flag=0;
					break;
				}
				if(fabs(w[nm])+anorm == anorm)
					break;
			}
			if(flag)
			{
				c=0.0;			/* cancellation of rv1[l] if l > 1 */
				s=1.0;
				for(i=l;i<=k;i++)
				{
					f=s*rv1[i];
					if(fabs(f)+anorm != anorm)
					{
						g=w[i];
						h=hypot(f,g);
						w[i]=h;
						h=1.0/h;
						c=g*h;
						s=(-1.0*f*h);
						for(j=0;j<m;j++)
						{
							y=a[j*n+nm];	/* a[j][nm] */
							z=a[j*n+i];	/* a[j][i] */
							a[j*n+nm]=(y*c)+(z*s);	/* a[j][nm] */
							a[j*n+i]=(z*c)-(y*s);	/* a[j][i] */
						}
					}
				}
			}
			z=w[k];
			if(l == k)		/* convergence */
			{
				if(z < 0.0)	/* singular value is made positive */
				{
					w[k]= -1.0*z;
					for(j=0;j<n;j++)
						v[j*n+k] *= (-1.0);	/* v[j][k] */
				}
				break;
			}
			if (its == 30)
			{
				fprintf(stderr,"Error in SVDCMP: No convergence in 30 iterations\n");
				exit(-1);
			}
			x=w[l];		/* shift from bottom 2-by-2 minor */
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=hypot(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			
				/* next QR transformation */
			
			c=s=1.0;
			for(j=l;j<=nm;j++)
			{
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=hypot(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=(x*c)+(g*s);
				g=(g*c)-(x*s);
				h=y*s;
				y=y*c;
				for(jj=0;jj<n;jj++)
				{
					x=v[jj*n+j];	/* v[jj][j] */
					z=v[jj*n+i];	/* v[jj][i] */
					v[jj*n+j]=(x*c)+(z*s);	/* v[jj][j] */
					v[jj*n+i]=(z*c)-(x*s);	/* v[jj][i] */
				}
				z=hypot(f,h);
				w[j]=z;		/* rotation can be arbitrary if z=0 */
				if(z)
				{
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for(jj=0;jj<m;jj++)
				{
					y=a[jj*n+j];	/* a[jj][j] */
					z=a[jj*n+i];	/* a[jj][i] */
					a[jj*n+j]=(y*c)+(z*s);	/* a[jj][j] */
					a[jj*n+i]=(z*c)-(y*s);	/* a[jj][i] */
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free((void *)rv1);
}

/* Given the singular value decomposition of a matrix a[0...m-1][0...n-1]
	solve the system of equations ax=b for x.  Input the matrices 
	U[0....m-1][0...n-1],w[0...n-1], and V[0...n-1][0...n-1] determined from
	svdcmp.  Also input the matrix b[0...m-1][0....k-1] and the solution vector
	x[0....k-1][0....n-1] is output. Singular values whose ratio to the maximum
	singular value are smaller than cutoff are zeroed out. The matrix U is
	overwritten.
	
*/
	

GMT_LONG solve_svd (REAL *u, GMT_LONG m, GMT_LONG n, REAL *v, REAL *w, REAL *b, GMT_LONG k, REAL *x, double cutoff)
{
	REAL *ut,sing_max;
	GMT_LONG i,j, n_use = 0;

	void mat_trans (REAL a[], GMT_LONG mrow, GMT_LONG ncol, REAL at[]);
	void mat_mult (REAL a[], GMT_LONG mrow, GMT_LONG ncol, REAL b[], GMT_LONG kcol, REAL c[]);
	   
		/* allocate work space */
		
	ut=(REAL *)calloc(n*m,sizeof(REAL));	/* space for the transpose */
	if(ut == NULL)
	{
		fprintf(stderr,"Error in solve_svd: Can't allocate work space\n");
		exit(-1);
	}
	
		/* find maximum singular value */
	
	sing_max=w[0];
	for(i=1;i<n;i++) sing_max=MAX(sing_max,w[i]);
		
	/* loop through singular values removing small ones */
		
	for(i=0;i<n;i++)
	{
		if ((w[i]/sing_max) > cutoff) {
			w[i]=1.0/w[i];
			n_use++;
		}
		else
			w[i]=0.0;
	}
	
	/* multiply V by 1/w */
	
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			v[j*n+i] *= w[i];
			
	/* get transpose of U */
		
	mat_trans(u,m,n,ut);
	
	/* multiply v(1/w)ut  -> this overwrites the matrix U */
		
	mat_mult(v,n,n,ut,m,u);
	
	/* multiply this result by b to get x */
		
	mat_mult(u,n,m,b,k,x);

	/* free work space */

	free((void *)ut);
	
	return (n_use);
}

/* Normalization parameters are stored in the coeff array which holds up to 7 terms
 * coeff[0]:	The mean x coordinate
 * coeff[1]:	The mean y coordinate
 * coeff[2]:	The mean w coordinate
 * coeff[3]:	The linear x slope
 * coeff[4]:	The linear y slope
 * coeff[5]:	The offset for the detrended data (combined into coeff[2])
 * coeff[6]:	The normalizing scale for the detrended data
 */

void do_normalization_1d (double **X, double *obs, GMT_LONG n, GMT_LONG mode, double *coeff)
{	/* mode == 1 norm w; mode == 2: also remove linear trend, norm & 4: also normalize result by range */

	GMT_LONG	i;
	double	d;
	
	memset ((void *)coeff, 0, 5*sizeof(double));
	for (i = 0; i < n; i++) {	/* Find mean w-value */
		coeff[GMT_Z] += obs[i];
		if (mode == 1) continue;
		coeff[GMT_X] += X[i][GMT_X];
	}
	coeff[GMT_Z] /= n;

	if (mode & 2) {	/* Solve for LS linera trend using deviations from (0, 0, 0) */
		double	xx, zz, sxx, sxz;
		sxx = sxz = 0.0;
		coeff[GMT_X] /= n;
		for (i = 0; i < n; i++) {
			xx = X[i][GMT_X] - coeff[GMT_X];
			zz = obs[i] - coeff[GMT_Z];
			sxx += (xx * xx);
			sxz += (xx * zz);
		}
		if (sxx != 0.0) coeff[3] = sxz/ sxx;
	}
	
	/* Remove linear trend (or mean) */
	
	coeff[5] = DBL_MAX;	coeff[6] = -DBL_MAX;
	for (i = 0; i < n; i++) {	/* Also find min/max */
		obs[i] -= coeff[GMT_Z];
		if (mode == 2) obs[i] -= (coeff[3] * (X[i][GMT_X] - coeff[GMT_X]));
		if (obs[i] < coeff[5]) coeff[5] = obs[i];
		if (obs[i] > coeff[6]) coeff[6] = obs[i];
	}
	if (mode & 4) {	/* Normalize by range */
		coeff[6] -= coeff[5];	/* Range */
		d = (coeff[6] == 0.0) ? 1.0 : 1.0 / coeff[6];
		for (i = 0; i < n; i++) obs[i] = (obs[i] - coeff[5]) * d;	/* Normalize 0-1 */
		coeff[GMT_Z] += coeff[5];	/* Combine the two constants in one place */
	}
	
	/* Recover obs(x) = w_norm(x) * coeff[6] + coeff[GMT_Z] + coeff[3]*(x-coeff[GMT_X]) */

}

void do_normalization (double **X, double *obs, GMT_LONG n, GMT_LONG mode, double *coeff)
{	/* mode == 1 norm z; mode == 2: also remove plane, norm & 4: also normalize result by range */

	GMT_LONG	i;
	double	d;
	
	memset ((void *)coeff, 0, 5*sizeof(double));
	for (i = 0; i < n; i++) {	/* Find mean z-value */
		coeff[GMT_Z] += obs[i];
		if (mode == 1) continue;
		coeff[GMT_X] += X[i][GMT_X];
		coeff[GMT_Y] += X[i][GMT_Y];
	}
	coeff[GMT_Z] /= n;

	if (mode & 2) {	/* Solve for LS plane using deviations from (0, 0, 0) */
		double	xx, yy, zz, sxx, sxy, sxz, syy, syz;
		sxx = sxy = sxz = syy = syz = 0.0;
		coeff[GMT_X] /= n;
		coeff[GMT_Y] /= n;
		for (i = 0; i < n; i++) {

			xx = X[i][GMT_X] - coeff[GMT_X];
			yy = X[i][GMT_Y] - coeff[GMT_Y];
			zz = obs[i] - coeff[GMT_Z];

			sxx += (xx * xx);
			sxz += (xx * zz);
			sxy += (xx * yy);
			syy += (yy * yy);
			syz += (yy * zz);
		}

		d = sxx*syy - sxy*sxy;
		if (d != 0.0) {
			coeff[3] = (sxz*syy - sxy*syz)/d;
			coeff[4] = (sxx*syz - sxy*sxz)/d;
		}
	}
	
	/* Remove plane (or mean) */
	
	coeff[5] = DBL_MAX;	coeff[6] = -DBL_MAX;
	for (i = 0; i < n; i++) {	/* Also find min/max */
		obs[i] -= coeff[GMT_Z];
		if (mode == 2) obs[i] -= (coeff[3] * (X[i][GMT_X] - coeff[GMT_X]) + coeff[4] * (X[i][GMT_Y] - coeff[GMT_Y]));
		if (obs[i] < coeff[5]) coeff[5] = obs[i];
		if (obs[i] > coeff[6]) coeff[6] = obs[i];
	}
	if (mode & 4) {	/* Normalize by range */
		coeff[6] -= coeff[5];	/* Range */
		d = (coeff[6] == 0.0) ? 1.0 : 1.0 / coeff[6];
		for (i = 0; i < n; i++) obs[i] = (obs[i] - coeff[5]) * d;	/* Normalize 0-1 */
		coeff[GMT_Z] += coeff[5];	/* Combine the two constants in one place */
	}
	
	/* Recover obs(x,y) = w_norm(x,y) * coeff[6] + coeff[GMT_Z] + coeff[3]*(x-coeff[GMT_X]) + coeff[4]*(y-coeff[GMT_Y]) */

}

double undo_normalization (double *X, double w_norm, GMT_LONG mode, double *coeff)
{
	double w;
	w = w_norm;
	if (mode & 4) w_norm *= coeff[6];
	w += coeff[GMT_Z];
	if (mode & 2) w += coeff[3] * (X[GMT_X] - coeff[GMT_X]) + coeff[4] * (X[GMT_Y] - coeff[GMT_Y]);
	return (w);
}

void mat_trans (REAL a[], GMT_LONG mrow, GMT_LONG ncol, REAL at[])
{
	/* Return the transpose of a */
	GMT_LONG i, j;
	for (i = 0; i < ncol; i++) for (j = 0; j < mrow; j++) at[mrow*i+j] = a[ncol*j+i];
}

void mat_mult (REAL a[], GMT_LONG mrow, GMT_LONG ncol, REAL b[], GMT_LONG kcol, REAL c[])
{
	/* Matrix multiplication a * b = c */
	
	GMT_LONG i, j, k, ij;
	
	for (i = 0; i < kcol; i++) {
		for (j = 0; j < mrow; j++) {
			ij = j * kcol + i;
			c[ij] = 0.0;
			for (k = 0; k < ncol; k++) c[ij] += a[j * ncol + k] * b[k * kcol + i];
		}
	}
}


double get_radius (double *X0, double *X1, GMT_LONG dim)
{
	double r = 0.0;
	/* Get distance between the two points */
	switch (dim) {
		case 1:	/* 1-D, just get x difference */
			r = fabs (X0[GMT_X] - X1[GMT_X]);
			break;
		case 2:	/* 2-D or spherical surface */
			r = GMT_distance_func (X0[GMT_X], X0[GMT_Y], X1[GMT_X], X1[GMT_Y]);
			break;
		case 3:	/* 3-D Cartesian */
			r = hypot (X0[GMT_X] - X1[GMT_X], X0[GMT_Y] - X1[GMT_Y]);
			r = hypot (r, X0[GMT_Z] - X1[GMT_Z]);
			break;
	}
	return (r);
}

double get_dircosine (double *D, double *X0, double *X1, GMT_LONG dim, GMT_LONG baz)
{
	/* D, the directional cosines of the observed gradient:
	 * X0: Observation point.
	 * X1: Prediction point.
	 * Compute N, the direction cosine of X1-X2, then C = D dot N.
	 */
	GMT_LONG ii;
	double az, C = 0.0, N[3];
	
	switch (dim) {
		case 1:	/* 1-D, always 1 */
			C = 1.0;
			break;
		case 2:	/* 2-D */
			az = GMT_azimuth_func (X0[GMT_X], X0[GMT_Y], X1[GMT_X], X1[GMT_Y], baz);
			sincosd (az, &N[GMT_X], &N[GMT_Y]);
			for (ii = 0; ii < 2; ii++) C += D[ii] * N[ii];	/* Dot product of 2-D unit vectors */
			break;
		case 3:	/* 3-D */
			for (ii = 0; ii < 3; ii++) N[ii] = X1[ii] - X0[ii];	/* Difference vector */
			GMT_normalize3v (N);	/* Normalize to unit vector */
			C = GMT_dot3v (D, N);	/* Dot product of 3-D unit vectors */
			if (baz) C = -C;	/* The opposite direction for X0-X1 */
			break;
	}
	return (C);
}

void GMT_check_lattice_cpy (struct GMT_GRD_INFO *info, double *x_inc, double *y_inc, GMT_LONG *pixel, GMT_LONG *active)
{	/* Uses provided settings to initialize the lattice settings from
	 * the -R<grdfile> if it was given; else it does nothing.
	 */
	if (!info->active) return;	/* -R<grdfile> was not used; use existing settings */
	
	/* Here, -R<grdfile> was used and we will use the settings supplied by the grid file (unless overridden) */
	if (!active || *active == FALSE) {	/* -I not set separately */
		*x_inc = info->grd.x_inc;
		*y_inc = info->grd.y_inc;
	}
	if (pixel) {	/* An pointer not NULL was passed that indicates grid registration */
		/* If a -F like option was set then toggle grid setting, else use grid setting */
		*pixel = (*pixel) ? !info->grd.node_offset : info->grd.node_offset;
	}
	if (active) *active = TRUE;	/* When 4th arg is not NULL it is set to TRUE (for Ctrl->active args) */
}

void *New_greenspline_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GREENSPLINE_CTRL *C;
	
	C = (struct GREENSPLINE_CTRL *) GMT_memory (VNULL, 1, sizeof (struct GREENSPLINE_CTRL), "New_greenspline_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->S.mode = SANDWELL_1987_2D;
	return ((void *)C);
}

void Free_greenspline_Ctrl (struct GREENSPLINE_CTRL *C) {	/* Deallocate control structure */
	if (C->C.file) GMT_free ((void *)C->C.file);	
	if (C->G.file) GMT_free ((void *)C->G.file);	
	if (C->N.file) GMT_free ((void *)C->N.file);	
	if (C->S.arg) GMT_free ((void *)C->S.arg);	
	if (C->T.file) GMT_free ((void *)C->T.file);	
	GMT_free ((void *)C);	
}
