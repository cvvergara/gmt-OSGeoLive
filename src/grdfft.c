/*--------------------------------------------------------------------
 *	$Id: grdfft.c 9923 2012-12-18 20:45:53Z pwessel $
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
 *  grdfft.c
 *
 * A program to do various operations on grid files in the frequency domain.
 *
 * W.H.F. Smith, 7 March 1990.
 *
 * Version 2.0.2 has option to do power spectral estimates.
 * added by WHFSmith, 4 Feb 1992.
 *
 * 1 Sept 1998:  Fixed bug in -E argv[i][2] should be [i][j].
 *   Note to myself:  users want azimuthal power spectra, and
 *   cross-spectra between two grids.  --whfs
 * 5/17/2006 PW: Added enhanced -F option for Gaussian band-pass.
 * 6/17/2009 PW: Added enhanced -F option for Butterworth band-pass.
 *
 * Version:	4
 */

#include "gmt.h"

#ifndef M_LN2
#define M_LN2          0.69314718055994530942  /* log_e 2 */
#endif

struct GRDFFT_CTRL {
	GMT_LONG n_op_count, n_par;
	GMT_LONG *operation;
	double	*par;

	struct A {	/* -A<azimuth> */
		GMT_LONG active;
		double value;
	} A;
	struct C {	/* -C<zlevel> */
		GMT_LONG active;
		double value;
	} C;
	struct D {	/* -D[<scale>|g] */
		GMT_LONG active;
		double value;
	} D;
	struct E {	/* -E[x_or_y][w] */
		GMT_LONG active;
		GMT_LONG mode;
	} E;
	struct F {	/* -F[x_or_y]<lc>/<lp>/<hp>/<hc> or -F[x_or_y]<lo>/<hi> */
		GMT_LONG active;
		GMT_LONG mode;
		double lc, lp, hp, hc;
	} F;
	struct G {	/* -G<outfile> */
		GMT_LONG active;
		char *file;
	} G;
	struct I {	/* -I[<scale>|g] */
		GMT_LONG active;
		double value;
	} I;
	struct L {	/* -L */
		GMT_LONG active;
	} L;
	struct M {	/* -M */
		GMT_LONG active;
	} M;
	struct N {	/* -N<stuff> */
		GMT_LONG active;
		double value;
	} N;
	struct S {	/* -S<scale> */
		GMT_LONG active;
		double value;
	} S;
	struct T {	/* -T<te/rl/rm/rw/ri> */
		GMT_LONG active;
		double te, rhol, rhom, rhow, rhoi;
	} T;
};

#ifndef FSIGNIF
#define FSIGNIF 24
#endif
#define UP_DOWN_CONTINUE 0
#define AZIMUTHAL_DERIVATIVE 1
#define DIFFERENTIATE	2
#define INTEGRATE	3
#define ISOSTASY	4
#define FILTER_EXP	5
#define FILTER_BW	6
#define FILTER_COS	7
#define SPECTRUM	8

/* Macro definition ij_data(i,j) finds the array index to an element
	containing the real data(i,j) in the padded complex array:  */

#define	ij_data(i,j) (2*(nx2*((j)+j_data_start)+(i)+i_data_start))


GMT_LONG	narray[2];
GMT_LONG	i, j, k, n, nx2 = 0, ny2 = 0, ndatac, i_data_start, j_data_start;

GMT_LONG	map_units = FALSE, force_narray = FALSE, suggest_narray = FALSE, leave_trend_alone = FALSE, n_user_set = FALSE;

double	data_var, data2_var, delta_kx, delta_ky;
double	a[3];	/* Plane fitting coefficients  */
double	mGal_at_45 = 980619.9203; /* Moritz's 1980 IGF value for gravity in mGal at 45 degrees latitude */

char	*infile = NULL;

struct GRD_HEADER h;
struct F_INFO {
	double	lc[3];		/* Low-cut frequency for r, x, and y	*/
	double	lp[3];		/* Low-pass frequency for r, x, and y	*/
	double	hp[3];		/* High-pass frequency for r, x, and y	*/
	double	hc[3];		/* High-cut frequency for r, x, and y	*/
	double	ltaper[3];	/* Low taper width for r, x, and y	*/
	double	htaper[3];	/* High taper width for r, x, and y	*/
	double	llambda[3];	/* Low full-wavelength where Gauss amp = 0.5 for r, x, and y	*/
	double	hlambda[3];	/* High full-wavelength where Gauss amp = 0.5  for r, x, and y	*/
	double	bw_order;	/* Order, N, of Butterworth filter	*/
	PFD filter;		/* Points to the correct filter function */
	GMT_LONG	do_this[3];	/* T/F this filter wanted for r, x, and y	*/
	GMT_LONG	set_already;
	GMT_LONG	kind;		/* FILTER_EXP, FILTER_BW, FILTER_COS  */
	GMT_LONG	arg;		/* 0 = Gaussian, 1 = Butterworth, 2 = cosine taper,  */
} f_info;

struct FFT_SUGGESTION {
	GMT_LONG	nx;
	GMT_LONG	ny;
	GMT_LONG	worksize;	/* # single-complex elements needed in work array  */
	GMT_LONG	totalbytes;	/* (8*(nx*ny + worksize))  */
	double	run_time;
	double	rms_rel_err;
}	fft_sug[3];	/* [0] holds fastest, [1] most accurate, [2] least storage  */


float	*datac, *workc;
double	scale_out = 1.0;

double cosine_weight (double freq, int j), gauss_weight (double freq, int j), bw_weight (double freq, int j);

int main (int argc, char **argv)
{
	GMT_LONG op_count = 0, par_count = 0, filter_type = 0;

	GMT_LONG error = FALSE, stop, give_wavelength = FALSE;
	struct GRDFFT_CTRL *Ctrl = NULL;
	GMT_LONG parse_f_string (char *c), count_slashes (char *txt);
	void read_data(struct GRDFFT_CTRL *Ctrl, int argc, char **argv);
	void do_continuation(double *zlevel), do_azimuthal_derivative(double *azim), do_differentiate(double *par), do_integrate(double *par);
	void do_filter(void), do_spectrum(double *par, GMT_LONG give_wavelength), do_isostasy(double *par);
	double get_filter_weight (GMT_LONG k);

	void *New_grdfft_Ctrl (), Free_grdfft_Ctrl (struct GRDFFT_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDFFT_CTRL *) New_grdfft_Ctrl ();	/* Allocate and initialize a new control structure */

	for (i = 0; i < 3; i++) {
		f_info.lc[i] = f_info.lp[i] = -1.0;	/* Set negative, below valid frequency range  */
		f_info.hp[i] = f_info.hc[i] = DBL_MAX;	/* Set huge positive, above valid frequency range  */
		f_info.ltaper[i] = f_info.htaper[i] = 0.0;	/* 1/width of taper, zero when not used  */
		f_info.do_this[i] = FALSE;
	}
	f_info.set_already = FALSE;
	
        for (i = 1; i < argc; i++) {
                if (argv[i][0] == '-') {
                        switch (argv[i][1]) {
                        
				/* Common parameters */
			
				case 'V':
				case '\0':
					error += GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;
				
				/* Supplemental parameters */
			
                                case 'A':
					Ctrl->n_op_count++;
					Ctrl->operation = (GMT_LONG *) GMT_memory ((void *)Ctrl->operation, (size_t)Ctrl->n_op_count, sizeof(GMT_LONG), GMT_program);
                                        Ctrl->operation[op_count] = AZIMUTHAL_DERIVATIVE;
					op_count++;
					Ctrl->n_par++;
					Ctrl->par = (double *) GMT_memory ((void *)Ctrl->par, (size_t)Ctrl->n_par, sizeof(double), GMT_program);
					if ((sscanf(&argv[i][2], "%lf", &Ctrl->par[par_count])) != 1) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -A option:  Cannot read azimuth\n", GMT_program);
						error++;
					}
					par_count++;
                                        break;
                                case 'C':
					Ctrl->n_op_count++;
					Ctrl->operation = (GMT_LONG *) GMT_memory ((void *)Ctrl->operation, (size_t)Ctrl->n_op_count, sizeof(GMT_LONG), GMT_program);
                                        Ctrl->operation[op_count] = UP_DOWN_CONTINUE;
					op_count++;
					Ctrl->n_par++;
					Ctrl->par = (double *) GMT_memory ((void *)Ctrl->par, (size_t)Ctrl->n_par, sizeof(double), GMT_program);
					if ((sscanf(&argv[i][2], "%lf", &Ctrl->par[par_count])) != 1) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -C option:  Cannot read zlevel\n", GMT_program);
						error++;
					}
					par_count++;
                                        break;
                                case 'D':
					Ctrl->n_op_count++;
					Ctrl->operation = (GMT_LONG *) GMT_memory ((void *)Ctrl->operation, (size_t)Ctrl->n_op_count, sizeof(GMT_LONG), GMT_program);
                                        Ctrl->operation[op_count] = DIFFERENTIATE;
					op_count++;
					Ctrl->n_par++;
					Ctrl->par = (double *) GMT_memory ((void *)Ctrl->par, (size_t)Ctrl->n_par, sizeof(double), GMT_program);
					Ctrl->par[par_count] = 0.0;
                                        if (argv[i][2]) {
						Ctrl->par[par_count] = (argv[i][2] == 'g' || argv[i][2] == 'G') ? mGal_at_45 : atof (&argv[i][2]);
						if (Ctrl->par[par_count] == 0.0) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR -D option:  scale must be nonzero\n", GMT_program);
							error++;
						}
					}
					par_count++;
                                        break;
                                case 'E':
					Ctrl->n_op_count++;
					Ctrl->operation = (GMT_LONG *) GMT_memory ((void *)Ctrl->operation, (size_t)Ctrl->n_op_count, sizeof(GMT_LONG), GMT_program);
                                        Ctrl->operation[op_count] = SPECTRUM;
					op_count++;
					Ctrl->n_par++;
					Ctrl->par = (double *) GMT_memory ((void *)Ctrl->par, (size_t)Ctrl->n_par, sizeof(double), GMT_program);
					Ctrl->par[par_count] = 0.0;
					j = 2;
					while(argv[i][j]) {
 						if (argv[i][j] == 'x' || argv[i][j] == 'X')
							Ctrl->par[par_count] = 1.0;
						else if (argv[i][j] == 'y' || argv[i][j] == 'Y')
							Ctrl->par[par_count] = -1.0;
						else if (argv[i][j] == 'w' || argv[i][j] == 'W')
							give_wavelength = TRUE;
						j++;
					}
					par_count++;
					Ctrl->E.active = TRUE;
                                        break;
                                case 'F':
                                	if (!(f_info.set_already)) {
						filter_type = count_slashes (&argv[i][2]);
						f_info.kind = FILTER_EXP + (filter_type - 1);
						Ctrl->n_op_count++;
						Ctrl->operation = (GMT_LONG *) GMT_memory ((void *)Ctrl->operation, (size_t)Ctrl->n_op_count, sizeof(GMT_LONG), GMT_program);
                                        	Ctrl->operation[op_count] = f_info.kind;
                                        	op_count++;
                                        	f_info.set_already = TRUE;
                                        }
                                        if (parse_f_string(&argv[i][2])) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -F option: ", GMT_program);
                                        	error++;
                                        }
                                        break;
                                case 'G':
                                        Ctrl->G.file = strdup (&argv[i][2]);
                                        break;
                                case 'I':
					Ctrl->n_op_count++;
					Ctrl->operation = (GMT_LONG *) GMT_memory ((void *)Ctrl->operation, (size_t)Ctrl->n_op_count, sizeof(GMT_LONG), GMT_program);
                                        Ctrl->operation[op_count] = INTEGRATE;
					op_count++;
					Ctrl->n_par++;
					Ctrl->par = (double *) GMT_memory ((void *)Ctrl->par, (size_t)Ctrl->n_par, sizeof(double), GMT_program);
   					Ctrl->par[par_count] = 0.0;
					if (argv[i][2]) {
						Ctrl->par[par_count] = (argv[i][2] == 'g' || argv[i][2] == 'G') ? mGal_at_45 : atof (&argv[i][2]);
						if (Ctrl->par[par_count] == 0.0) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR -I option:  scale must be nonzero\n", GMT_program);
							error++;
						}
					}
					par_count++;
                                        break;
#ifdef DEBUG
				case 'Q':
	 				Ctrl->n_op_count++;
					Ctrl->operation = (GMT_LONG *) GMT_memory ((void *)Ctrl->operation, (size_t)Ctrl->n_op_count, sizeof(GMT_LONG), GMT_program);
	                                Ctrl->operation[op_count] = -1;
					op_count++;
                                        break;
#endif
                                case 'L':
                                        Ctrl->L.active = TRUE;
                                        break;
                                case 'M':
                                        Ctrl->M.active = TRUE;
                                        break;
                                case 'N':
					if (argv[i][2] == 'f' || argv[i][2] == 'F')
						force_narray = TRUE;
					else if (argv[i][2] == 'q' || argv[i][2] == 'Q')
						suggest_narray = TRUE;
					else {
	                                        sscanf(&argv[i][2], "%" GMT_LL "d/%" GMT_LL "d", &nx2, &ny2);
						n_user_set = TRUE;
					}
                                        break;
                                case 'S':
                                        Ctrl->S.active = TRUE;
                                        scale_out = (argv[i][2] == 'd' || argv[i][2] == 'D') ? 1.0e6: atof (&argv[i][2]);
                                        break;
                                case 'T':
                                        Ctrl->T.active = TRUE;
					Ctrl->n_op_count++;
					Ctrl->operation = (GMT_LONG *) GMT_memory ((void *)Ctrl->operation, (size_t)Ctrl->n_op_count, sizeof(GMT_LONG), GMT_program);
                                        Ctrl->operation[op_count] = ISOSTASY;
					op_count++;
 					Ctrl->n_par += 5;
					Ctrl->par = (double *) GMT_memory ((void *)Ctrl->par, (size_t)Ctrl->n_par, sizeof(double), GMT_program);
                                        n = sscanf (&argv[i][2], "%lf/%lf/%lf/%lf/%lf", &Ctrl->par[par_count],
						&Ctrl->par[par_count+1], &Ctrl->par[par_count+2], &Ctrl->par[par_count+3], &Ctrl->par[par_count+4]);
					for (j = 1, k = 0; j < 5; j++) if (Ctrl->par[par_count+j] < 0.0) k++;
					if (n != 5 || k > 0) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -T option.  Correct syntax:\n", GMT_program);
						fprintf (stderr, "\t-T<te>/<rhol>/<rhom>/<rhow>/<rhoi>, all densities >= 0\n");
						error++;
					}
					par_count += 5;
                                        Ctrl->L.active = TRUE;
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

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "grdfft %s - Perform mathematical operations on grid files in the wavenumber (or frequency) domain\n\n", GMT_VERSION);
		fprintf (stderr,"usage: grdfft <in_grdfile> [-G<out_grdfile>] [-A<azimuth>] [-C<zlevel>]\n");
		fprintf (stderr,"\t[-D[<scale>|g]] [-E[x_or_y][w]] [-F[x_or_y]<parameters>] [-I[<scale>|g]] [-L] [-M]\n");
		fprintf (stderr,"\t[-N<stuff>] [-S<scale>] [-T<te/rl/rm/rw/ri>] [-V]\n\n");
		
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
		
		fprintf (stderr,"\tin_grdfile is the input netCDF grid file.\n");
		fprintf (stderr, "\tOPTIONS:\n");
		fprintf (stderr,"\t-G filename for output netCDF grid file.\n");
		fprintf (stderr,"\t-A<azimuth> Take azimuthal derivative along line <azimuth> degrees CW from North.\n");
		fprintf (stderr,"\t-C<zlevel>  Continue field upward (+) or downward (-) to zlevel (meters).\n");
		fprintf (stderr,"\t-D Differentiate, i.e., multiply by kr [ * scale].  Use -Dg to get mGal from m.\n");
		fprintf (stderr,"\t-E Estimate spEctrum of r [x] [y].  Write f, power[f], 1 std dev(power[f]) to stdout.\n");
		fprintf (stderr,"\t\tAppend w to write wavelength instead of frequency.\n");
		fprintf (stderr,"\t-F Filter r [x] [y] freq according to one of three kinds of filter specifications:\n");
		fprintf (stderr,"\t   a) Cosine band-pass: Append four wavelengths <lc>/<lp>/<hp>/<hc>.\n");
		fprintf (stderr,"\t      freq outside <lc>/<hc> are cut; inside <lp>/<hp> are passed, rest are tapered.\n");
		fprintf (stderr,"\t      Replace wavelength by - to skip, e.g., -F-/-/500/100 is a low-pass filter.\n");
		fprintf (stderr,"\t   b) Gaussian band-pass: Append two wavelengths <lo>/<hi> where filter amplitudes = 0.5.\n");
		fprintf (stderr,"\t      Replace wavelength by - to skip, e.g., -F300/- is a high-pass Gaussian filter.\n");
		fprintf (stderr,"\t   c) Butterworth band-pass: Append two wavelengths and order <lo>/<hi>/<order> where filter amplitudes = 0.5.\n");
		fprintf (stderr,"\t      Replace wavelength by - to skip, e.g., -F300/-/2 is a high-pass 2nd-order Butterworth filter.\n");
		fprintf (stderr,"\t-I Integrate, i.e., divide by kr [ * scale].  Use -Ig to get m from mGal.\n");
		fprintf (stderr,"\t-L Leave trend alone.  Do not remove least squares plane from data [Default removes plane].\n");
		fprintf (stderr,"\t-M Map units used.  Convert grid dimensions from degrees to meters.\n");
		fprintf (stderr,"\t-N<stuff>  Choose or inquire about suitable grid dimensions for FFT.\n");
		fprintf (stderr,"\t\t-Nf will force the FFT to use the dimensions of the data.\n");
		fprintf (stderr,"\t\t-Nq will inQuire about more suitable dimensions.\n");
		fprintf (stderr,"\t\t-N<nx>/<ny> will do FFT on array size <nx>/<ny> (Must be >= grid size).\n");
		fprintf (stderr,"\t\tDefault chooses dimensions >= data which optimize speed, accuracy of FFT.\n");
		fprintf (stderr,"\t\tIf FFT dimensions > grid dimensions, data are extended and tapered to zero.\n");
		fprintf (stderr,"\t-S multiply field by scale after inverse FFT [1.0].\n");
		fprintf (stderr,"\t   Give -Sd to convert deflection of vertical to micro-radians.\n");
		fprintf (stderr,"\t-T Compute isostatic response.  Input file is topo load. Append elastic thickness,\n");
		fprintf (stderr,"\t   and densities of load, mantle, water, and infill, all in SI units.\n");
		fprintf (stderr,"\t   It also implicitly sets -L.\n");
		GMT_explain_option ('V');
		fprintf (stderr,"\tList operations in the order desired for execution.\n");
		exit (EXIT_FAILURE);
	}

	if (!(Ctrl->n_op_count)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify at least one operation\n", GMT_program);
		error++;
	}
	if (n_user_set && (nx2 <= 0 || ny2 <= 0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  nx2 and/or ny2 <= 0\n", GMT_program);
		error++;
	}
	if (scale_out == 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  scale must be nonzero\n", GMT_program);
		error++;
	}
	if (!infile) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input file\n", GMT_program);
		error++;
	}
	if (Ctrl->E.active && Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  -G ignored (stdout used for ASCII output)\n", GMT_program);
		error++;
	}
	if (!Ctrl->E.active && !Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option:  Must specify output file\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	read_data (Ctrl, argc, argv);

	/* Check that no NaNs are present */
	
	stop = FALSE;
	for (j = 0; !stop && j < h.ny; j++) for (i = 0; !stop && i < h.nx; i++) stop = GMT_is_fnan (datac[ij_data(i,j)]);
	if (stop) {
		fprintf (stderr, "%s: Input grid cannot have NaNs!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

#ifdef FTEST
	{
		double f = 0.0;
		while (f < 3.0) {
			printf ("%g\t%g\n", f, f_info.filter (f, 0));	/* Radial filter */
			f += 0.01;
		}
		exit (-1);
	}
#endif
	if (gmtdefs.verbose) fprintf (stderr, "%s: forward FFT...", GMT_program);
	
	GMT_fourt (datac, narray, 2, -1, 1, workc);

	for (op_count = par_count = 0; op_count < Ctrl->n_op_count; op_count++) {
		switch (Ctrl->operation[op_count]) {
			case UP_DOWN_CONTINUE:
				if (gmtdefs.verbose) ((Ctrl->par[par_count] < 0.0) ? fprintf (stderr, "downward continuation...") : fprintf (stderr, "upward continuation..."));
				do_continuation (&Ctrl->par[par_count]);
				par_count++;
				break;
			case AZIMUTHAL_DERIVATIVE:
				if (gmtdefs.verbose) fprintf (stderr, "azimuthal derivative...");
				do_azimuthal_derivative (&Ctrl->par[par_count]);
				par_count++;
				break;
			case DIFFERENTIATE:
				if (gmtdefs.verbose) fprintf (stderr, "differentiate...");
				do_differentiate (&Ctrl->par[par_count]);
				par_count++;
				break;
			case INTEGRATE:
				if (gmtdefs.verbose) fprintf (stderr, "integrate...");
				do_integrate (&Ctrl->par[par_count]);
				par_count++;
				break;
			case ISOSTASY:
				if (gmtdefs.verbose) fprintf (stderr, "isostasy...");
				do_isostasy (&Ctrl->par[par_count]);
				par_count += 5;
				break;
			case FILTER_COS:
				if (gmtdefs.verbose) fprintf (stderr, "cosine filter...");
				do_filter ();
				break;
			case FILTER_EXP:
				if (gmtdefs.verbose) fprintf (stderr, "Gaussian filter...");
				do_filter ();
				break;
			case FILTER_BW:
				if (gmtdefs.verbose) fprintf (stderr, "Butterworth filter...");
				do_filter ();
				break;
			case SPECTRUM:
				if (gmtdefs.verbose) fprintf (stderr, "spectrum...");
				do_spectrum (&Ctrl->par[par_count], give_wavelength);
				par_count += 2;
				break;
		}
	}

	if (Ctrl->G.file) {

		if (gmtdefs.verbose) fprintf (stderr, "inverse FFT...");

		GMT_fourt (datac, narray, 2, 1, 1, workc);

		scale_out *= (2.0 / ndatac);
		for (i = 0; i < ndatac; i+= 2) datac[i] *= (float)scale_out;

		/* The data are in the middle of the padded array:  */

		GMT_err_fail (GMT_write_grd (Ctrl->G.file, &h, datac, h.x_min, h.x_max, h.y_min, h.y_max, GMT_pad, TRUE), Ctrl->G.file);
	}
	
	GMT_free ((void *)workc);
	GMT_free ((void *)datac);

	if (gmtdefs.verbose) fprintf (stderr, "Done\n");

	Free_grdfft_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void read_data (struct GRDFFT_CTRL *Ctrl, int argc, char **argv)
{
	GMT_LONG	worksize, factors[32];
	double	tdummy, edummy;
	void	suggest_fft(GMT_LONG nx, GMT_LONG ny, struct FFT_SUGGESTION *fft_sug, GMT_LONG do_print);
	void	fourt_stats(GMT_LONG nx, GMT_LONG ny, GMT_LONG *f, double *r, GMT_LONG *s, double *t);
	void remove_plane(void), taper_edges(void);
	
	GMT_err_fail (GMT_read_grd_info (infile, &h), infile);
	
	GMT_grd_init (&h, argc, argv, TRUE);

	/* Get dimensions as may be appropriate:  */
	if (n_user_set) {
		if (nx2 < h.nx || ny2 < h.ny) {
			fprintf(stderr,"%s: Error: You specified -Nnx/ny smaller than input grid.  Ignored.\n", GMT_program);
			n_user_set = FALSE;
		}
	}
	if (!(n_user_set) ) {
		if (force_narray) {
			nx2 = h.nx;
			ny2 = h.ny;
		}
		else {
			suggest_fft(h.nx, h.ny, fft_sug, (gmtdefs.verbose || suggest_narray));
			if (fft_sug[1].totalbytes < fft_sug[0].totalbytes) {
				/* The most accurate solution needs same or less storage
				 * as the fastest solution; use the most accurate's dimensions:  */
				nx2 = fft_sug[1].nx;
				ny2 = fft_sug[1].ny;
			}
			else {
				/* Use the sizes of the fastest solution  */
				nx2 = fft_sug[0].nx;
				ny2 = fft_sug[0].ny;
			}
		}
	}

	/* Get here when nx2 and ny2 are set to the vals we will use.  */
	narray[0] = nx2;
	narray[1] = ny2;
	fourt_stats(nx2, ny2, factors, &edummy, &worksize, &tdummy);
	if (gmtdefs.verbose) fprintf(stderr,"%s:  Data dimension %d %d\tFFT dimension %ld %ld\n", GMT_program,
		h.nx, h.ny, nx2, ny2);

	/* Make an array of floats 2 * nx2 * ny2 for complex data:  */
	ndatac = 2 * nx2 * ny2;
        datac = (float *) GMT_memory (VNULL, (size_t)ndatac, sizeof(float), GMT_program);
	memset ((void *)datac, 0, (size_t)(ndatac*sizeof(float))); 
	if (worksize) {
                if (worksize < nx2) worksize = nx2;
                if (worksize < ny2) worksize = ny2;
		worksize *= 2;
        	workc = (float *) GMT_memory (VNULL, (size_t)worksize, sizeof(float), GMT_program);
		memset ((void *)workc, 0, (size_t)(worksize*sizeof(float))); 
	}
	else {
		workc = (float *) GMT_memory (VNULL, (size_t)4, sizeof(float), GMT_program);
		memset ((void *)workc, 0, (size_t)(4*sizeof(float))); 
	}

	/* Put the data in the middle of the padded array:  */

	i_data_start = GMT_pad[0] = (nx2 - h.nx)/2;	/* zero if nx2 < h.nx+1  */
	j_data_start = GMT_pad[3] = (ny2 - h.ny)/2;
	GMT_pad[1] = nx2 - h.nx - GMT_pad[0];
	GMT_pad[2] = ny2 - h.ny - GMT_pad[3];

	GMT_err_fail (GMT_read_grd (infile, &h, datac, h.x_min, h.x_max, h.y_min, h.y_max, GMT_pad, TRUE), infile);

	if (!(Ctrl->L.active) )remove_plane();
	if (!(force_narray) )taper_edges();

	delta_kx = 2 * M_PI / (nx2 * h.x_inc);
	delta_ky = 2 * M_PI / (ny2 * h.y_inc);
	if (Ctrl->M.active) {
		/* Give delta_kx, delta_ky units of 2pi/meters  */
		delta_kx /= (project_info.DIST_M_PR_DEG * cosd (0.5 * (h.y_min + h.y_max)) );
		delta_ky /= project_info.DIST_M_PR_DEG;
	}
}

void remove_plane (void)
{
	/* Remove the best-fitting plane by least squares.

	Let plane be z = a0 + a1 * x + a2 * y.  Choose the
	center of x,y coordinate system at the center of 
	the array.  This will make the Normal equations 
	matrix G'G diagonal, so solution is trivial.  Also,
	spend some multiplications on normalizing the 
	range of x,y into [-1,1], to avoid roundoff error.  */

	GMT_LONG	one_or_zero;
	double	x_half_length, one_on_xhl, y_half_length, one_on_yhl;
	double	sumx2, sumy2, x, y, z;

	one_or_zero = (h.node_offset) ? 0 : 1;
	x_half_length = 0.5 * (h.nx - one_or_zero);
	one_on_xhl = 1.0 / x_half_length;
	y_half_length = 0.5 * (h.ny - one_or_zero);
	one_on_yhl = 1.0 / y_half_length;

	sumx2 = sumy2 = data_var = 0.0;
	a[2] = a[1] = a[0] = 0.0;

	for (j = 0; j < h.ny; j++) {
		y = one_on_yhl * (j - y_half_length);
		for (i = 0; i < h.nx; i++) {
			x = one_on_xhl * (i - x_half_length);
			z = datac[ij_data(i,j)];
			a[0] += z;
			a[1] += z*x;
			a[2] += z*y;
			sumx2 += x*x;
			sumy2 += y*y;
		}
	}
	a[0] /= (h.nx*h.ny);
	a[1] /= sumx2;
	a[2] /= sumy2;
	for (j = 0; j < h.ny; j++) {
		y = one_on_yhl * (j - y_half_length);
		for (i = 0; i < h.nx; i++) {
			x = one_on_xhl * (i - x_half_length);
			datac[ij_data(i,j)] -= (float)(a[0] + a[1]*x + a[2]*y);
			data_var += (datac[ij_data(i,j)] * datac[ij_data(i,j)]);
		}
	}
	data_var = sqrt(data_var / (h.nx*h.ny - 1));
	/* Rescale a1,a2 into user's units, in case useful later:  */
	a[1] *= (2.0/(h.x_max - h.x_min));
	a[2] *= (2.0/(h.y_max - h.y_min));
	if(gmtdefs.verbose)fprintf (stderr,"%s: Plane removed.  Mean, S.D., Dx, Dy:  %.8g\t%.8g\t%.8g\t%.8g\n", GMT_program, a[0],data_var,a[1],a[2]);
}

void taper_edges (void)
{
	GMT_LONG	im, jm, il1, ir1, il2, ir2, jb1, jb2, jt1, jt2;
	double	scale, cos_wt;

	/* Note that if nx2 = h.nx+1 and ny2 = h.ny + 1, then this routine
		will do nothing; thus a single row/column of zeros may be
		added to the bottom/right of the input array and it cannot
		be tapered.  But when (nx2 - h.nx)%2 == 1 or ditto for y,
		this is zero anyway.  */


	/* First reflect about xmin and xmax, point symmetric about edge point:  */

	for (im = 1; im <= i_data_start; im++) {
		il1 = -im;	/* Outside xmin; left of edge 1  */
		ir1 = im;	/* Inside xmin; right of edge 1  */
		il2 = il1 + h.nx - 1;	/* Inside xmax; left of edge 2  */
		ir2 = ir1 + h.nx - 1;	/* Outside xmax; right of edge 2  */
		for (j = 0; j < h.ny; j++) {
			datac[ij_data(il1,j)] = (float)2.0*datac[ij_data(0,j)] - datac[ij_data(ir1,j)];
			datac[ij_data(ir2,j)] = (float)2.0*datac[ij_data((h.nx-1),j)] - datac[ij_data(il2,j)];
		}
	}

	/* Next, reflect about ymin and ymax.
		At the same time, since x has been reflected,
		we can use these vals and taper on y edges:  */

	scale = M_PI / (j_data_start + 1);

	for (jm = 1; jm <= j_data_start; jm++) {
		jb1 = -jm;	/* Outside ymin; bottom side of edge 1  */
		jt1 = jm;	/* Inside ymin; top side of edge 1  */
		jb2 = jb1 + h.ny - 1;	/* Inside ymax; bottom side of edge 2  */
		jt2 = jt1 + h.ny - 1;	/* Outside ymax; bottom side of edge 2  */
		cos_wt = 0.5 * (1.0 + cos(jm * scale) );
		for (i = -i_data_start; i < nx2 - i_data_start; i++) {
			datac[ij_data(i,jb1)] = (float)(cos_wt * (2.0*datac[ij_data(i,0)] - datac[ij_data(i,jt1)]));
			datac[ij_data(i,jt2)] = (float)(cos_wt * (2.0*datac[ij_data(i,(h.ny-1))] - datac[ij_data(i,jb2)]));
		}
	}

	/* Now, cos taper the x edges:  */

	scale = M_PI / (i_data_start + 1);
	for (im = 1; im <= i_data_start; im++) {
		il1 = -im;
		ir1 = im;
		il2 = il1 + h.nx - 1;
		ir2 = ir1 + h.nx - 1;
		cos_wt = 0.5 * (1.0 + cos(im * scale) );
		for (j = -j_data_start; j < ny2 - j_data_start; j++) {
			datac[ij_data(il1,j)] *= (float)cos_wt;
			datac[ij_data(ir2,j)] *= (float)cos_wt;
		}
	}
}

double kx (GMT_LONG k)
{
	/* Return the value of kx given k,
		where kx = 2 pi / lambda x,
		and k refers to the position
		in the datac array, datac[k].  */

	GMT_LONG	ii;

	ii = (k/2)%nx2;
	if (ii > nx2/2) ii -= nx2;
	return(ii * delta_kx);
}

double ky (GMT_LONG k)
{
	/* Return the value of ky given k,
		where ky = 2 pi / lambda y,
		and k refers to the position
		in the datac array, datac[k].  */

	GMT_LONG	jj;

	jj = (k/2)/nx2;
	if (jj > ny2/2) jj -= ny2;
	return(jj * delta_ky);
}

double modk (GMT_LONG k)
{
	/* Return the value of sqrt(kx*kx + ky*ky),
		given k, where k is array position.  */

	return (hypot (kx(k), ky(k)));
}

void do_differentiate (double *par)
{
	GMT_LONG	k;
	double scale, fact;

	/* Differentiate in frequency domain by multiplying by kr [scale optional] */
	
	scale = (*par != 0.0) ? *par : 1.0;
	datac[0] = datac[1] = 0.0;
	for (k = 2; k < ndatac; k+= 2) {
		fact = scale * modk(k);
		datac[k] *= (float)fact;
		datac[k+1] *= (float)fact;
	}
}

void do_integrate (double *par)
{
	/* Integrate in frequency domain by dividing by kr [scale optional] */
	GMT_LONG	k;
	double fact, scale;

	scale = (*par != 0.0) ? *par : 1.0;
	datac[0] = datac[1] = 0.0;
	for (k = 2; k < ndatac; k+= 2) {
		fact = 1.0 / (scale * modk(k) );
		datac[k] *= (float)fact;
		datac[k+1] *= (float)fact;
	}
}

void do_continuation (double *zlevel)
{
	GMT_LONG	k;
	double tmp;

	/* If z is positive, the field will be upward continued using exp[- k z].  */

	for (k = 2; k < ndatac; k+= 2) {
		tmp = exp (-(*zlevel) * modk(k) );
		datac[k] *= (float)tmp;
		datac[k+1] *= (float)tmp;
	}
}

void do_azimuthal_derivative (double *azim)
{
	GMT_LONG	k;
	float	tempr, tempi, fact;
	double cos_azim, sin_azim;

	cos_azim = cosd (*azim);
	sin_azim = sind (*azim);

	datac[0] = datac[1] = 0.0;
	for (k = 2; k < ndatac; k+= 2) {
		fact = (float)(sin_azim * kx(k) + cos_azim * ky(k));
		tempr = -(datac[k+1] * fact);
		tempi = (datac[k] * fact);
		datac[k] = tempr;
		datac[k+1] = tempi;
	}
}

void do_isostasy (double *par)
{

	/* Do the isostatic response function convolution in the Freq domain.
	All units assumed to be in SI (that is kx, ky, modk wavenumbers in m**-1,
	densities in kg/m**3, Te in m, etc.
	rw, the water density, is used to set the Airy ratio and the restoring
	force on the plate (rm - ri)*gravity if ri = rw; so use zero for topo in air.  */
	GMT_LONG	k;
	double	airy_ratio, rigidity_d, d_over_restoring_force, mk, k2, k4, transfer_fn;

	double	te;	/* Elastic thickness, SI units (m)  */
	double	rl;	/* Load density, SI units  */
	double	rm;	/* Mantle density, SI units  */
	double	rw;	/* Water density, SI units  */
	double	ri;	/* Infill density, SI units  */
	
	double	youngs_modulus = 1.0e11;	/* Pascal = Nt/m**2  */
	double	normal_gravity = 9.80619203;	/* m/s**2  */
	double	poissons_ratio = 0.25;

	te = par[0];	rl = par[1];	rm = par[2];	rw = par[3];	ri = par[4];
	rigidity_d = (youngs_modulus * te * te * te) / (12.0 * (1.0 - poissons_ratio * poissons_ratio));
	d_over_restoring_force = rigidity_d / ( (rm - ri) * normal_gravity);
	airy_ratio = -(rl - rw)/(rm - ri);

	if (te == 0.0) {	/* Airy isostasy; scale global variable scale_out and return */
		scale_out *= airy_ratio;
		return;
	}
	
	for (k = 0; k < ndatac; k+= 2) {
		mk = modk(k);
		k2 = mk * mk;
		k4 = k2 * k2;
		transfer_fn = airy_ratio / ( (d_over_restoring_force * k4) + 1.0);
		datac[k] *= (float)transfer_fn;
		datac[k+1] *= (float)transfer_fn;
	}
}

void do_filter (void)
{
	GMT_LONG	k;
	double	weight;
	double get_filter_weight (GMT_LONG k);
	for (k = 0; k < ndatac; k += 2) {
		weight = get_filter_weight(k);
		datac[k] *= (float)weight;
		datac[k+1] *= (float)weight;
	}
}

double get_filter_weight (GMT_LONG k)
{
	double	freq, return_value = 1.0;
	GMT_LONG	j;
		
	for (j = 0; j < 3; j++) {
		if (!(f_info.do_this[j])) continue;	/* Only do one of x, y, or r filtering */
		switch (j) {
			case 0:
				freq = modk(k);
				break;
			case 1:
				freq = kx(k);
				break;
			default:
				freq = ky(k);
				break;
		}
		return_value = f_info.filter (freq, j);
	}

	return (return_value);
}

double gauss_weight (double freq, int j) {
	double hi, lo;
	lo = (f_info.llambda[j] == -1.0) ? 0.0 : exp (-M_LN2 * pow (freq * f_info.llambda[j], 2.0));	/* Low-pass part */
	hi = (f_info.hlambda[j] == -1.0) ? 1.0 : exp (-M_LN2 * pow (freq * f_info.hlambda[j], 2.0));	/* Hi-pass given by its complementary low-pass */
	return (hi - lo);
}

double bw_weight (double freq, int j) {
	double hi, lo;
	lo = (f_info.llambda[j] == -1.0) ? 0.0 : sqrt (1.0 / (1.0 + pow (freq * f_info.llambda[j], f_info.bw_order)));	/* Low-pass part */
	hi = (f_info.hlambda[j] == -1.0) ? 1.0 : sqrt (1.0 / (1.0 + pow (freq * f_info.hlambda[j], f_info.bw_order)));	/* Hi-pass given by its complementary low-pass */
	return (hi - lo);
}

double cosine_weight (double freq, int j) {
	if (freq <= f_info.lc[j] || freq >= f_info.hc[j]) return(0.0);	/* In fully cut range.  Weight is zero.  */
	if (freq > f_info.lc[j] && freq < f_info.lp[j]) return (0.5 * (1.0 + cos (M_PI * (freq - f_info.lp[j]) * f_info.ltaper[j])));
	if (freq > f_info.hp[j] && freq < f_info.hc[j]) return (0.5 * (1.0 + cos (M_PI * (freq - f_info.hp[j]) * f_info.htaper[j])));
	return (1.0);	/* Freq is in the fully passed range, so weight is multiplied by 1.0  */
}

GMT_LONG count_slashes (char *txt)
{
	GMT_LONG i, n;
	for (i = n = 0; txt[i]; i++) if (txt[i] == '/') n++;
	return (n);
}

GMT_LONG parse_f_string (char *c)
{
	GMT_LONG	i, j, n_tokens, pos, descending;
	double	fourvals[4];
	char	line[GMT_LONG_TEXT], p[GMT_LONG_TEXT];
	
	/* Syntax is either -F[x|y]lc/hc/lp/hp (Cosine taper), -F[x|y]lo/hi (Gaussian), or 0F[x|y]lo/hi/order (Butterworth) */
	
	strcpy(line, c);
	i = j = 0;	/* j is Filter type:  r=0, x=1, y=2  */
	
	if (line[i] == 'x') {
		j = 1;
		i++;
	}
	else if (line[i] == 'y') {
		j = 2;
		i++;
	}
	
	f_info.do_this[j] = TRUE;
	fourvals[0] = fourvals[1] = fourvals[2] = fourvals[3] = -1.0;
	
	n_tokens = pos = 0;
	while ((GMT_strtok (&line[i], "/", &pos, p))) {
		if (n_tokens > 3) {
			fprintf(stderr,"%s: Too many slashes in -F.\n", GMT_program);
			return(TRUE);
		}
		if(p[0] == '-') {
			fourvals[n_tokens] = -1.0;
		}
		else {
			if ((sscanf(p, "%lf", &fourvals[n_tokens])) != 1) {
				fprintf(stderr,"%s:  Cannot read token %ld.\n", GMT_program, n_tokens);
				return(TRUE);
			}
		}
		n_tokens++;
	}
	
	if (!(n_tokens == 2 || n_tokens == 3 || n_tokens == 4)) {
		fprintf(stderr,"%s: -F Cannot find 2-4 tokens separated by slashes.\n", GMT_program);
		return(TRUE);
	}
	descending = TRUE;
	if (f_info.kind == FILTER_BW && n_tokens == 3) n_tokens = 2;	/* So we dont check the order as a wavelength */

	for (i = 1; i < n_tokens; i++) {
		if (fourvals[i] == -1.0 || fourvals[i-1] == -1.0) continue;
		if (fourvals[i] > fourvals[i-1]) descending = FALSE;
	}
	if (!(descending)) {
		fprintf(stderr,"%s: -F Wavelengths are not in descending order.\n", GMT_program);
		return(TRUE);
	}
	if (f_info.kind == FILTER_COS) {	/* Cosine band-pass specification */
		if ( (fourvals[0] * fourvals[1]) < 0.0 || (fourvals[2] * fourvals[3]) < 0.0) {
			fprintf(stderr,"%s: -F Pass/Cut specification error.\n", GMT_program);
			return(TRUE);
		}
	
		/* Now everything is OK  */
	
		if (fourvals[0] >= 0.0 || fourvals[1] >= 0.0) {	/* Lower end values are set  */
			f_info.lc[j] = (2.0 * M_PI)/fourvals[0];
			f_info.lp[j] = (2.0 * M_PI)/fourvals[1];
			if (fourvals[0] != fourvals[1]) f_info.ltaper[j] = 1.0/(f_info.lc[j] - f_info.lp[j]);
		}
	
		if (fourvals[2] >= 0.0 || fourvals[3] >= 0.0) {	/* Higher end values are set  */
			f_info.hp[j] = (2.0 * M_PI)/fourvals[2];
			f_info.hc[j] = (2.0 * M_PI)/fourvals[3];
			if (fourvals[2] != fourvals[3]) f_info.htaper[j] = 1.0/(f_info.hc[j] - f_info.hp[j]);
		}
		f_info.filter = (PFD) cosine_weight;
	}
	else if (f_info.kind == FILTER_BW) {	/* Butterworth specification */
		f_info.llambda[j] = (fourvals[0] == -1.0) ? -1.0 : fourvals[0] / TWO_PI;	/* TWO_PI is used to counteract the 2*pi in the wavenumber */
		f_info.hlambda[j] = (fourvals[1] == -1.0) ? -1.0 : fourvals[1] / TWO_PI;
		f_info.bw_order = 2.0 * fourvals[2];
		f_info.filter = (PFD) bw_weight;
	}
	else {	/* Gaussian half-amp specifications */
		f_info.llambda[j] = (fourvals[0] == -1.0) ? -1.0 : fourvals[0] / TWO_PI;	/* TWO_PI is used to counteract the 2*pi in the wavenumber */
		f_info.hlambda[j] = (fourvals[1] == -1.0) ? -1.0 : fourvals[1] / TWO_PI;
		f_info.filter = (PFD) gauss_weight;
	}
	f_info.arg = f_info.kind - FILTER_EXP;
	return(FALSE);
}

void do_spectrum (double *par, GMT_LONG give_wavelength)
{
	/* This is modeled on the 1-D case, using the following ideas:
	 *	In 1-D, we ensemble average over samples of length L = 
	 *	n * dt.  This gives n/2 power spectral estimates, at
	 *	frequencies i/L, where i = 1, n/2.  If we have a total
	 *	data set of ndata, we can make nest=ndata/n independent
	 *	estimates of the power.  Our standard error is then
	 *	1/sqrt(nest).
	 *	In making 1-D estimates from 2-D data, we may choose
	 *	n and L from nx2 or ny2 and delta_kx, delta_ky as 
	 *	appropriate.  In this routine, we are giving the sum over
	 * 	all frequencies in the other dimension; that is, an
	 *	approximation of the integral.
	 */

	GMT_LONG	k, nk, nused, ifreq;
	double	delta_k, r_delta_k, freq, *power, eps_pow;
	PFD	get_k;
	char	format[GMT_TEXT_LEN], buffer[GMT_LONG_TEXT];
	/*  Added by WHFS in version 3.001 to normalize power  :  */
	double	powfactor;

	if (*par > 0.0) {
		/* X spectrum desired  */
		delta_k = delta_kx;
		nk = nx2/2;
		get_k = (PFD)kx;
	}
	else if (*par < 0.0) {
		/* Y spectrum desired  */
		delta_k = delta_ky;
		nk = ny2/2;
		get_k = (PFD)ky;
	}
	else {
		/* R spectrum desired  */
		if (delta_kx < delta_ky) {
			delta_k = delta_kx;
			nk = nx2/2;
		}
		else {
			delta_k = delta_ky;
			nk = ny2/2;
		}
		get_k = (PFD)modk;
	}

	/* Get an array for summing stuff:  */
	power = (double *) GMT_memory (VNULL, (size_t)nk, sizeof(double), GMT_program);
	for (k = 0; k < nk; k++) power[k] = 0.0;

	/* Loop over it all, summing and storing, checking range for r:  */

	r_delta_k = 1.0 / delta_k;
	
	for (nused = 0, k = 2; k < ndatac; k+= 2) {
		freq = (*get_k)(k);
		ifreq = irint(fabs(freq)*r_delta_k) - 1;
		if (ifreq < 0) ifreq = 0;	/* Might happen when doing r spectrum  */
		if (ifreq >= nk) continue;	/* Might happen when doing r spectrum  */
		power[ifreq] += (datac[k]*datac[k] + datac[k+1]*datac[k+1]);
		nused++;
	}

	/* Now get here when array is summed.  */
	eps_pow = 1.0/sqrt((double)nused/(double)nk);
	delta_k /= (2.0*M_PI);	/* Write out frequency, not wavenumber  */
	sprintf (format, "%s\t%s\t%s\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
	powfactor = 4.0 / pow ((double)ndatac, 2.0);
	for (k = 0; k < nk; k++) {
		freq = (k + 1) * delta_k;
		if (give_wavelength) freq = 1.0/freq;
		power[k] *= powfactor;
		sprintf (buffer, format, freq, power[k], eps_pow * power[k]);
		GMT_fputs (buffer, GMT_stdout);
	}
}

void fourt_stats (GMT_LONG nx, GMT_LONG ny, GMT_LONG *f, double *r, GMT_LONG *s, double *t)
{
	/* Find the proportional run time, t, and rms relative error, r,
	 * of a Fourier transform of size nx,ny.  Also gives s, the size
	 * of the workspace that will be needed by the transform.
	 * To use this routine for a 1-D transform, set ny = 1.
	 * 
	 * This is all based on the comments in Norman Brenner's code
	 * FOURT, from which our C codes are translated.
	 * Brenner says:
	 * r = 3 * pow(2, -FSIGNIF) * sum{ pow(prime_factors, 1.5) }
	 * where FSIGNIF is the smallest bit in the floating point fraction.
	 * 
	 * Let m = largest prime factor in the list of factors.
	 * Let p = product of all primes which appear an odd number of
	 * times in the list of prime factors.  Then the worksize needed
	 * s = max(m,p).  However, we know that if n is radix 2, then no
	 * work is required; yet this formula would say we need a worksize
	 * of at least 2.  So I will return s = 0 when max(m,p) = 2.
	 *
	 * I have two different versions of the comments in FOURT, with
	 * different formulae for t.  The simple formula says 
	 * 	t = n * (sum of prime factors of n).
	 * The more complicated formula gives coefficients in microsecs
	 * on a cdc3300 (ancient history, but perhaps proportional):
	 *	t = 3000 + n*(500 + 43*s2 + 68*sf + 320*nf),
	 * where s2 is the sum of all factors of 2, sf is the sum of all
	 * factors greater than 2, and nf is the number of factors != 2.
	 * We know that factors of 2 are very good to have, and indeed,
	 * Brenner's code calls different routines depending on whether
	 * the transform is of size 2 or not, so I think that the second
	 * formula is more correct, using proportions of 43:68 for 2 and
	 * non-2 factors.  So I will use the more complicated formula.
	 * However, I realize that the actual numbers are wrong for today's
	 * architectures, and the relative proportions may be wrong as well.
	 * 
	 * W. H. F. Smith, 26 February 1992.
	 *  */

	GMT_LONG	n_factors, i, sum2, sumnot2, nnot2;
	GMT_LONG	nonsymx, nonsymy, nonsym, storage;
	GMT_LONG ntotal;
	double	err_scale, sig_bits = FSIGNIF;
	GMT_LONG get_prime_factors(GMT_LONG n, GMT_LONG *f), get_non_symmetric_f(GMT_LONG *f, GMT_LONG n);

	/* Find workspace needed.  First find non_symmetric factors in nx, ny  */
	n_factors = get_prime_factors((GMT_LONG)nx, f);
	nonsymx = get_non_symmetric_f(f, n_factors);
	n_factors = get_prime_factors((GMT_LONG)ny, f);
	nonsymy = get_non_symmetric_f(f, n_factors);
	nonsym = MAX(nonsymx, nonsymy);

	/* Now get factors of ntotal  */
	ntotal = ((size_t)nx) * ((size_t)ny);
        n_factors = get_prime_factors (ntotal, f);
	storage = MAX(nonsym, f[n_factors - 1]);
	if (storage == 2)
		*s = 0;
	else
		*s = storage;

	/* Now find time and error estimates:  */

	err_scale = 0.0;
	sum2 = 0;
	sumnot2 = 0;
	nnot2 = 0;
	for(i = 0; i < n_factors; i++) {
		if (f[i] == 2)
			sum2 += f[i];
		else {
			sumnot2 += f[i];
			nnot2++;
		}
		err_scale += pow((double)f[i], 1.5);
	}
	*t = 1.0e-06*(3000.0 + ntotal * (500.0 + 43.0*sum2 + 68.0*sumnot2 + 320.0*nnot2));
	*r = err_scale * 3.0 * pow(2.0, -sig_bits);
	return;
} 

GMT_LONG get_non_symmetric_f (GMT_LONG *f, GMT_LONG n)
{
	/* Return the product of the non-symmetric factors in f[]  */
	GMT_LONG	i = 0, j = 1, retval = 1;

	if (n == 1) return(f[0]);

	while(i < n) {
		while(j < n && f[j] == f[i]) j++;
		if ((j-i)%2) retval *= f[i];
		i = j;
		j = i + 1;
	}
	if (retval == 1) retval = 0;	/* There are no non-sym factors  */
	return(retval);
}

void suggest_fft (GMT_LONG nx, GMT_LONG ny, struct FFT_SUGGESTION *fft_sug, GMT_LONG do_print)
{
	GMT_LONG	f[64], xstop, ystop;
	GMT_LONG	nx_best_t, ny_best_t;
	GMT_LONG	nx_best_e, ny_best_e;
	GMT_LONG	nx_best_s, ny_best_s;
        GMT_LONG     nxg, nyg;       /* Guessed by this routine  */
        GMT_LONG     nx2, ny2, nx3, ny3, nx5, ny5;   /* For powers  */
	double	current_time, best_time, given_time, s_time, e_time;
	GMT_LONG	current_space, best_space, given_space, e_space, t_space;
	double	current_err, best_err, given_err, s_err, t_err;
	void	fourt_stats(GMT_LONG nx, GMT_LONG ny, GMT_LONG *f, double *r, GMT_LONG *s, double *t);


	fourt_stats(nx, ny, f, &given_err, &given_space, &given_time);
	given_space += nx*ny;
	given_space *= 8;
	if (do_print) fprintf(stderr,"%s:  Data dimension\t%ld %ld\ttime factor %.8g\trms error %.8e\tbytes %ld\n", GMT_program,
		nx, ny, given_time, given_err, given_space);

	best_err = s_err = t_err = given_err;
	best_time = s_time = e_time = given_time;
	best_space = t_space = e_space = given_space;
	nx_best_e = nx_best_t = nx_best_s = nx;
	ny_best_e = ny_best_t = ny_best_s = ny;

	xstop = 2 * nx;
	ystop = 2 * ny;

        for (nx2 = 2; nx2 <= xstop; nx2 *= 2) {
          for (nx3 = 1; nx3 <= xstop; nx3 *= 3) {
            for (nx5 = 1; nx5 <= xstop; nx5 *= 5) {
                nxg = nx2 * nx3 * nx5;
                if (nxg < nx || nxg > xstop) continue;

                for (ny2 = 2; ny2 <= ystop; ny2 *= 2) {
                  for (ny3 = 1; ny3 <= ystop; ny3 *= 3) {
                    for (ny5 = 1; ny5 <= ystop; ny5 *= 5) {
                        nyg = ny2 * ny3 * ny5;
                        if (nyg < ny || nyg > ystop) continue;

			fourt_stats(nxg, nyg, f, &current_err, &current_space, &current_time);
			current_space += nxg*nyg;
			current_space *= 8;
			if (current_err < best_err) {
				best_err = current_err;
				nx_best_e = nxg;
				ny_best_e = nyg;
				e_time = current_time;
				e_space = current_space;
			}
			if (current_time < best_time) {
				best_time = current_time;
				nx_best_t = nxg;
				ny_best_t = nyg;
				t_err = current_err;
				t_space = current_space;
			}
			if (current_space < best_space) {
				best_space = current_space;
				nx_best_s = nxg;
				ny_best_s = nyg;
				s_time = current_time;
				s_err = current_err;
			}

		    }
		  }
		}

	    }
	  }
	}

	if (do_print) {
		fprintf(stderr,"%s:  Highest speed\t%ld %ld\ttime factor %.8g\trms error %.8e\tbytes %ld\n", GMT_program,
			nx_best_t, ny_best_t, best_time, t_err, t_space);
		fprintf(stderr,"%s:  Most accurate\t%ld %ld\ttime factor %.8g\trms error %.8e\tbytes %ld\n", GMT_program,
			nx_best_e, ny_best_e, e_time, best_err, e_space);
		fprintf(stderr,"%s:  Least storage\t%ld %ld\ttime factor %.8g\trms error %.8e\tbytes %lds\n", GMT_program,
			nx_best_s, ny_best_s, s_time, s_err, best_space);
	}
	/* Fastest solution:  */
	fft_sug[0].nx = nx_best_t;
	fft_sug[0].ny = ny_best_t;
	fft_sug[0].worksize = (t_space/8) - (nx_best_t * ny_best_t);
	fft_sug[0].totalbytes = t_space;
	fft_sug[0].run_time = best_time;
	fft_sug[0].rms_rel_err = t_err;
	/* Most accurate solution:  */
	fft_sug[1].nx = nx_best_e;
	fft_sug[1].ny = ny_best_e;
	fft_sug[1].worksize = (e_space/8) - (nx_best_e * ny_best_e);
	fft_sug[1].totalbytes = e_space;
	fft_sug[1].run_time = e_time;
	fft_sug[1].rms_rel_err = best_err;
	/* Least storage solution:  */
	fft_sug[2].nx = nx_best_s;
	fft_sug[2].ny = ny_best_s;
	fft_sug[2].worksize = (best_space/8) - (nx_best_s * ny_best_s);
	fft_sug[2].totalbytes = best_space;
	fft_sug[2].run_time = s_time;
	fft_sug[2].rms_rel_err = s_err;

	return;
}

GMT_LONG get_prime_factors (GMT_LONG n, GMT_LONG *f)
{
	/* Fills the integer array f with the prime factors of n.
	 * Returns the number of locations filled in f, which is
	 * one if n is prime.
	 *
	 * f[] should have been malloc'ed to enough space before
	 * calling prime_factors().  We can be certain that f[32]
	 * is enough space, for if n fits in a long, then n < 2**32,
	 * and so it must have fewer than 32 prime factors.  I think
	 * that in general, ceil(log2((double)n)) is enough storage
	 * space for f[].
	 *
	 * Tries 2,3,5 explicitly; then alternately adds 2 or 4
	 * to the previously tried factor to obtain the next trial
	 * factor.  This is done with the variable two_four_toggle.
	 * With this method we try 7,11,13,17,19,23,25,29,31,35,...
	 * up to a maximum of sqrt(n).  This shortened list results
	 * in 1/3 fewer divisions than if we simply tried all integers
	 * between 5 and sqrt(n).  We can reduce the size of the list
	 * of trials by an additional 20% by removing the multiples
	 * of 5, which are equal to 30m +/- 5, where m >= 1.  Starting
	 * from 25, these are found by alternately adding 10 or 20.
	 * To do this, we use the variable ten_twenty_toggle.
	 *
	 * W. H. F. Smith, 26 Feb 1992, after D.E. Knuth, vol. II  */

	GMT_LONG	current_factor;	/* The factor currently being tried  */
	GMT_LONG	max_factor;	/* Don't try any factors bigger than this  */
	GMT_LONG	n_factors = 0;	/* Returned; one if n is prime  */
	GMT_LONG	two_four_toggle = 0;	/* Used to add 2 or 4 to get next trial factor  */
	GMT_LONG	ten_twenty_toggle = 0;	/* Used to add 10 or 20 to skip_five  */
	GMT_LONG	skip_five = 25;	/* Used to skip multiples of 5 in the list  */
	GMT_LONG	m;	/* Used to keep a working copy of n  */


	/* Initialize m and max_factor  */
	
	m = GMT_abs(n);
	if (m < 2) return(0);
	max_factor = (GMT_LONG)floor(sqrt((double)m));

	/* First find the 2s  */
	current_factor = 2;
	while(!(m%current_factor)) {
		m /= current_factor;
		f[n_factors] = (GMT_LONG)current_factor;
		n_factors++;
	}
	if (m == 1) return(n_factors);

	/* Next find the 3s  */
	current_factor = 3;
	while(!(m%current_factor)) {
		m /= current_factor;
		f[n_factors] = (GMT_LONG)current_factor;
		n_factors++;
	}
	if (m == 1) return(n_factors);

	/* Next find the 5s  */
	current_factor = 5;
	while(!(m%current_factor)) {
		m /= current_factor;
		f[n_factors] = (GMT_LONG)current_factor;
		n_factors++;
	}
	if (m == 1) return(n_factors);

	/* Now try all the rest  */

	while (m > 1 && current_factor <= max_factor) {

		/* Current factor is either 2 or 4 more than previous value  */

		if (two_four_toggle) {
			current_factor += 4;
			two_four_toggle = 0;
		}
		else {
			current_factor += 2;
			two_four_toggle = 1;
		}

		/* If current factor is a multiple of 5, skip it.  But first,
			set next value of skip_five according to 10/20 toggle:  */

		if (current_factor == skip_five) {
			if (ten_twenty_toggle) {
				skip_five += 20;
				ten_twenty_toggle = 0;
			}
			else {
				skip_five += 10;
				ten_twenty_toggle = 1;
			}
			continue;
		}

		/* Get here when current_factor is not a multiple of 2,3 or 5:  */

		while(!(m%current_factor)) {
			m /= current_factor;
			f[n_factors] = (GMT_LONG)current_factor;
			n_factors++;
		}
	}

	/* Get here when all factors up to floor(sqrt(n)) have been tried.  */

	if (m > 1) {
		/* m is an additional prime factor of n  */
		f[n_factors] = (GMT_LONG)m;
		n_factors++;
	}
	return (n_factors);
}

void *New_grdfft_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDFFT_CTRL *C;
	
	C = (struct GRDFFT_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDFFT_CTRL), "New_grdfft_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */

	return ((void *)C);
}

void Free_grdfft_Ctrl (struct GRDFFT_CTRL *C) {	/* Deallocate control structure */
	if (C->operation) GMT_free ((void *)C->operation);	
	if (C->par) GMT_free ((void *)C->par);	
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}
