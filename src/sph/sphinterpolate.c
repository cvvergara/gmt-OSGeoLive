/*--------------------------------------------------------------------
 *	$Id: sphinterpolate.c 9923 2012-12-18 20:45:53Z pwessel $
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
 * Spherical gridding in tension.  We read input data and want to create
 * a grid using various interpolants on a sphere.  This program relies
 * on two Fortan F77 libraries by Renka:
 * Renka, R, J,, 1997, Algorithm 772: STRIPACK: Delaunay Triangulation
 *  and Voronoi Diagram on the Surface of a Sphere, AMC Trans. Math.
 *  Software, 23 (3), 416-434.
 * Renka, R, J,, 1997, Algorithm 773: SSRFPACK: Interpolation of scattered
 *  data on the Surface of a Sphere with a surface under tension, AMC
 *  Trans. Math. Software, 23 (3), 435-442.
 * We translate to C using f2c -r8 and link with -lf2c
 *
 * Author:      Paul Wessel
 * Date:        9-APR-2008
 *
 */
 
#include "gmt.h"
#include "sph.h"

struct SPHINTERPOLATE_CTRL {
	struct F {	/* -F */
		GMT_LONG active;
	} F;
	struct G {	/* -G<grdfile> */
		GMT_LONG active;
		char *file;
	} G;
	struct I {	/* -Idx[/dy] */
		GMT_LONG active;
		double xinc, yinc;
	} I;
	struct Q {	/* -Q<interpolation> */
		GMT_LONG active;
		int mode;
		double value[2];
	} Q;
	struct T {	/* -T for variable tension */
		GMT_LONG active;
	} T;
	struct Z {	/* -Z to scale data */
		GMT_LONG active;
	} Z;
};

int main (int argc, char **argv)
{
	GMT_LONG i, j, ij, ij_f, n = 0;
	GMT_LONG n_expected_fields, n_args, m, n_fields, fno, n_files = 0;
	size_t n_alloc = GMT_CHUNK, nm;

	GMT_LONG error = FALSE, nofile = TRUE, done = FALSE;

	double sx, sy, cx, cy, w_min, w_max, sf = 1.0;
	double *xx = NULL, *yy = NULL, *zz = NULL, *ww = NULL, *surfd = NULL, *in = NULL;
	float *surf = NULL;

	char line[BUFSIZ], *not_used = NULL;
	FILE *fp = NULL;

	struct SPHINTERPOLATE_CTRL *Ctrl = NULL;
	struct GRD_HEADER h;
	void *New_sphinterpolate_Ctrl (), Free_sphinterpolate_Ctrl (struct SPHINTERPOLATE_CTRL *C);
	GMT_LONG get_args (char *arg, double pars[], char *msg);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct SPHINTERPOLATE_CTRL *)New_sphinterpolate_Ctrl ();	/* Allocate and initialize a new control structure */
	GMT_grd_init (&h, argc, argv, FALSE);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
        
				/* Common parameters */
        
				case 'H':
				case 'M':
				case 'R':
				case 'V':
				case ':':
				case 'b':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], &h.x_min, &h.x_max, &h.y_min, &h.y_max);
					break;
        
				/* Supplemental parameters */
        
				case 'F':
					Ctrl->F.active = TRUE;
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
				case 'Q':
					Ctrl->Q.active = TRUE;
					switch (argv[i][2]) {
						case '0':	/* Linear */
							Ctrl->Q.mode = 0;
							break;
						case '1':
							Ctrl->Q.mode = 1;
							break;
						case '2':
							Ctrl->Q.mode = 2;
							if (argv[i][3] == '/') {	/* Gave optional n/m/dgmx */
								if ((m = get_args (&argv[i][4], Ctrl->Q.value, "-Q3/N[/M[/U]]")) < 0) error = TRUE;
							}
							break;
						case '3':
							Ctrl->Q.mode = 3;
							if (argv[i][3] == '/') {	/* Gave optional e/sm/niter */
								if ((m = get_args (&argv[i][4], Ctrl->Q.value, "-Q3/E[/U[/niter]]")) < 0) error = TRUE;
							}
							break;
						default:
							error = TRUE;
							fprintf (stderr, "GMT ERROR %s: -%c Mode must be in 0-3 range\n", GMT_program, argv[i][1]);
							break;
					}
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					break;
				case 'Z':
					Ctrl->Z.active = TRUE;
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
		fprintf (stderr, "sphinterpolate %s - Spherical gridding in tension of data on a sphere\n", GMT_VERSION);
		fprintf (stderr, "==> The hard work is done by algorithms 772 (STRIPACK) & 773 (SSRFPACK) by R. J. Renka [1997] <==\n\n");
		fprintf (stderr, "usage: sphinterpolate [<infiles>] -G<grdfile> %s [-F] [%s]\n", GMT_I_OPT, GMT_H_OPT);
		fprintf (stderr, "\t[-Q<mode>][/args] [-T] [-V] [-Z] [%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_b_OPT, GMT_m_OPT);
                
		fprintf (stderr, "\t-G Specify file name for the final gridded solution.\n");
		GMT_inc_syntax ('I', 0);
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
                
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\tinfiles (in ASCII) has 3 or more columns.  If no file(s) is given, standard input is read.\n");
		fprintf (stderr, "\t-F Force pixel registration for output grid [Default is gridline orientation]\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-Q Select tension factors to achive the following [Default is no tension]:\n");
		fprintf (stderr, "\t   0: Piecewise linear interpolation ; no tension [Default]\n");
		fprintf (stderr, "\t   1: Smooth interpolation with local gradient estimates.\n");
		fprintf (stderr, "\t   2: Smooth interpolation with global gradient estimates and tension.  Optionally append /N/M/U:\n");
		fprintf (stderr, "\t      N = Number of iterations to converge solutions for gradients and variable tensions (-T only) [3]\n");
		fprintf (stderr, "\t      M = Number of Gauss-Seidel iterations when determining gradients [10]\n");
		fprintf (stderr, "\t      U = Maximum change in a gradient at the last iteration [0.01]\n");
		fprintf (stderr, "\t   3: Smoothing.  Optionally append /E/U/N, where\n");
		fprintf (stderr, "\t      E = Expected squared error in a typical (scaled) data value [0.01]\n");
		fprintf (stderr, "\t      U = Upper bound on  weighted sum of squares of deviations from data [npoints]\n");
		fprintf (stderr, "\t      N = Number of iterations to converge solutions for gradients and variable tensions (-T only) [3]\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t   If no region is specified we default to the entire world [-Rg]\n");
		fprintf (stderr, "\t-T Use variable tension (ignored for -Q0) [constant]\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t-Z Scale data by 1/(max-min) prior to gridding [no scaling]\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 3 input columns\n");
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&Ctrl->I.xinc, &Ctrl->I.yinc, &Ctrl->F.active, &Ctrl->I.active);

	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 3;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 3) {
            fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 3 columns\n", GMT_program);
		error++;
	}
	if (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
		error = TRUE;
	}
	if (!Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G:  Must specify output file\n", GMT_program);
		error = TRUE;
	}
	if (Ctrl->Q.mode < 0 || Ctrl->Q.mode > 3) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -T:  Must specify a mode in the 0-3 range\n", GMT_program);
		error = TRUE;
	}
	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}
	if (!project_info.region_supplied) {	/* Default is global region */
		h.x_min = 0.0;	h.x_max = 360.0;	h.y_min = -90.0;	h.y_max = 90.0;
	}

#ifdef SET_IO_MODE
	GMT_setmode (GMT_OUT);
#endif

	/* Now we are ready to take on some input values */

	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;
	n_args = (argc > 1) ? argc : 3;
	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 3;

	n_alloc = GMT_CHUNK;
	xx = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	yy = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	zz = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	ww = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	n = 0;
	w_min = DBL_MAX;	w_max = -DBL_MAX;
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

		if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) not_used = GMT_fgets (line, BUFSIZ, fp);

		while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

			if (GMT_io.status & GMT_IO_MISMATCH) {
				fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld\n", GMT_program, n_fields, n_expected_fields, n);
				exit (EXIT_FAILURE);
			}
			while (GMT_io.status & GMT_IO_SEGMENT_HEADER) {	/* Segment header, get next record */
				n_fields = GMT_input (fp, &n_expected_fields, &in);
			}
			sincosd (in[GMT_Y], &sy, &cy);
			sincosd (in[GMT_X], &sx, &cx);
			xx[n] = cy * cx;
			yy[n] = cy * sx;
			zz[n] = sy;
			ww[n] = in[GMT_Z];
			if (Ctrl->Z.active) {
				if (ww[n] < w_min) w_min = ww[n];
				if (ww[n] > w_max) w_max = ww[n];
			}
			n++;

			if (n == (int)n_alloc) {	/* Get more memory */
				n_alloc <<= 1;
				xx = (double *) GMT_memory ((void *)xx, (size_t)n_alloc, sizeof (double), GMT_program);
				yy = (double *) GMT_memory ((void *)yy, (size_t)n_alloc, sizeof (double), GMT_program);
				zz = (double *) GMT_memory ((void *)zz, (size_t)n_alloc, sizeof (double), GMT_program);
				ww = (double *) GMT_memory ((void *)ww, (size_t)n_alloc, sizeof (double), GMT_program);
			}
		}

		if (fp != GMT_stdin) GMT_fclose (fp);
	}

	xx = (double *) GMT_memory ((void *)xx, (size_t)n, sizeof (double), GMT_program);
	yy = (double *) GMT_memory ((void *)yy, (size_t)n, sizeof (double), GMT_program);
	zz = (double *) GMT_memory ((void *)zz, (size_t)n, sizeof (double), GMT_program);
	ww = (double *) GMT_memory ((void *)ww, (size_t)n, sizeof (double), GMT_program);

	if (gmtdefs.verbose) fprintf (stderr, "%s: Do Delaunay triangulation using %ld points\n", GMT_program, n);

	if (Ctrl->Z.active && w_max > w_min) {	/* Scale the data */
		sf = 1.0 / (w_max - w_min);
		for (i = 0; i < n; i++) ww[i] *= sf;
	}
	
	/* Set up output grid */
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Evaluate output grid\n", GMT_program);
	h.node_offset = (int)Ctrl->F.active;
	h.x_inc = Ctrl->I.xinc;
	h.y_inc = Ctrl->I.yinc;
	GMT_RI_prepare (&h);	/* Ensure -R -I consistency and set nx, ny */
	GMT_err_fail (GMT_grd_RI_verify (&h, 1), Ctrl->G.file);

	h.xy_off = 0.5 * h.node_offset;
	nm = GMT_get_nm (h.nx, h.ny);
	surfd = (double *) GMT_memory (VNULL, (size_t)nm, sizeof(double), GMT_program);
	
	/* Do the interpolation */
	
	ssrfpack_grid (xx, yy, zz, ww, n, Ctrl->Q.mode, Ctrl->Q.value, Ctrl->T.active, gmtdefs.verbose, &h, surfd);
	GMT_free ((void *)xx);
	GMT_free ((void *)yy);
	GMT_free ((void *)zz);
	GMT_free ((void *)ww);
	
	/* Convert the doubles to float and unto the Fortran transpose order */
	
	sf = (w_max - w_min);
	surf = (float *) GMT_memory (VNULL, (size_t)nm, sizeof(float), GMT_program);
	for (j = ij = 0; j < h.ny; j++) {
		for (i = 0; i < h.nx; i++, ij++) {
			ij_f = i * h.ny + j;
			surf[ij] = (float)surfd[ij_f];
			if (Ctrl->Z.active) surf[ij] *= (float)sf;
		}
	}
	GMT_free ((void *)surfd);
	
	/* Write solution */
	
	GMT_err_fail (GMT_write_grd (Ctrl->G.file, &h, surf, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);
	
	/* Free variables */
	
	GMT_free ((void *)surf);

	if (gmtdefs.verbose) fprintf (stderr, "%s: Done!\n", GMT_program);

	Free_sphinterpolate_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_sphinterpolate_Ctrl () {	/* Allocate and initialize a new control structure */
	struct SPHINTERPOLATE_CTRL *C;
	
	C = (struct SPHINTERPOLATE_CTRL *) GMT_memory (VNULL, 1, sizeof (struct SPHINTERPOLATE_CTRL), "New_sphinterpolate_Ctrl");
	return ((void *)C);
}

void Free_sphinterpolate_Ctrl (struct SPHINTERPOLATE_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}

GMT_LONG get_args (char *arg, double par[], char *msg)
{
	GMT_LONG m;
	char txt_a[32], txt_b[32], txt_c[32];
	m = sscanf (arg, "%[^/]/%[^/]/%s", txt_a, txt_b, txt_c);
	if (m < 1) {
		fprintf (stderr, "GMT ERROR %s: %s\n", GMT_program, msg);
		m = -1;
	}
	par[0] = atof (txt_a);
	if (m >= 2) par[1] = atof (txt_b);
	if (m == 3) par[2] = atof (txt_c);
	return (m);
}
