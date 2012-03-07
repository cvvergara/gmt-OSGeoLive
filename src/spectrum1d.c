/*--------------------------------------------------------------------
 *	$Id: spectrum1d.c,v 1.53 2011/07/08 21:27:06 guru Exp $
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
 *
 *	Compute auto and cross spectra using Welch's method of
 *	multiple overlapped windows.  Find 1 standard error bars
 *	following expressions in Bendat & Piersol.
 *
 *	-D<dt>	set delta_time, the sampling of the timeseries.
 *		[Default <dt> = 1.0]
 *
 *	-N<name_stem>	name stem for filenames.  Files will be
 *		created called <name_stem>.xpower, etc.
 *		[Default <name_stem> = "spectrum"]
 *
 *	-S<segment_size>	give an integer radix 2 window width.
 *		Total timeseries will be split into pieces of
 *		length <segment_size>.  Std errors in spectra are
 *		approximately 1/sqrt(n_data/segment_size).
 *
 *	-V	Verbose operation; write info to stderr.
 *		[Default is silent]
 *
 *	-W	write Wavelength in col 1 instead of frequency.
 *		[Default writes frequency in cycles/dt_units]
 *
 *	-C<output flags>	input has 2 cols, X(t) Y(t); do cross-spectra.
 *		[Default is one column; do X power spectrum only]
 *		Optional string of 0 to 8 output flags { x y c n p a g o }
 *		0 or 8 produces all cross-spectra in output files.  x = x-power
 *		spectrum, y = y-power spectrum, c = coherent power spectrum,
 *		n = noise power spectrum, p = phase spectrum, a = admittance
 *		power spectrum, g = gain spectrum, o = squared coherency.
 *
 *	-b	input and/or output data are in binary form
 *
 *	Author:		W. H. F. Smith
 *	Date:		11 April 1991-2000
 *	Revised:	5 June 1991-2000 to add W and N options in prep for
 *			GMT v2.0 release.
 *			PW: Upgrade to GMT 3.1 w/ -b
 *			BCH-J: Upgraded -C: optional output flags
 *	References:	Julius S. Bendat & Allan G. Piersol,
 *			"Random Data", 2nd revised edition, 566pp.,
 *			1986, John Wiley & Sons, New York. [ B&P below]
 *
 *			Peter D. Welch, "The use of Fast Fourier
 *			Transform for the estimation of power spectra:
 *			a method based on time averaging over short,
 *			modified periodograms", IEEE Transactions on
 *			Audio and Electroacoustics, Vol AU-15, No 2,
 *			June, 1967.
 *  Version:	4.1.2	No longer uses global variables.
 */

#include "gmt.h"

#define SPECTRUM1D_N_OUTPUT_CHOICES 8

struct SPECTRUM1D_CTRL {
	struct C {	/* -C[<xycnpago>] */
		GMT_LONG active;
		char col[SPECTRUM1D_N_OUTPUT_CHOICES];	/* Character codes for desired output in the right order */
	} C;
	struct D {	/* -D<inc> */
		GMT_LONG active;
		double inc;
	} D;
	struct N {	/* -N<namestem> */
		GMT_LONG active;
		char *name;
	} N;
	struct S {	/* -S<segment_size> */
		GMT_LONG active;
		GMT_LONG size;
	} S;
	struct W {	/* -W */
		GMT_LONG active;
	} W;
};

struct SPECTRUM1D_INFO {	/* Control structure for spectrum1d */
	GMT_LONG	y_given;
	GMT_LONG	n_spec, window, window_2;
	float	*datac;
	double	dt, x_variance, y_variance, d_n_windows, y_pow;
	struct	SPEC {
		double	xpow;	/* PSD in X(t)  */
		double	ypow;	/* PSD in Y(t)  */
		double	gain;	/* Amplitude increase X->Y in Optimal Response Function  */
		double	phase;	/* Phase of Y w.r.t. X in Optimal Response Function  */
		double	coh;	/* (squared) Coherence of Y and X; SNR of Y = coh/(1-coh)  */
		double	radmit;	/* Real part of Admittance; used e.g., for gravity/topography  */
	} *spec;
};

int main (int argc, char **argv)
{
	GMT_LONG	error = FALSE;
	
	char	buffer[BUFSIZ], *not_used = NULL;

	GMT_LONG	i, k, j, n_expected_fields, n_fields, n_read, n_outputs, window_test = 2;
	GMT_LONG	n_alloc, n_data;

	float	*x = NULL, *y = NULL;

	double	*in = NULL;
	
	FILE	*fp = NULL;
	
	struct SPECTRUM1D_INFO C;
	struct SPECTRUM1D_CTRL *Ctrl = NULL;

	void alloc_arrays (struct SPECTRUM1D_INFO *C);
	void compute_spectra (struct SPECTRUM1D_INFO *C, float *x, float *y, GMT_LONG n_data);
	void write_output (struct SPECTRUM1D_INFO *C, char *cols, GMT_LONG n_outputs, GMT_LONG write_wavelength, char *namestem);
	void free_space (struct SPECTRUM1D_INFO *C);
	void *New_spectrum1d_Ctrl (), Free_spectrum1d_Ctrl (struct SPECTRUM1D_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct SPECTRUM1D_CTRL *)New_spectrum1d_Ctrl ();	/* Allocate and initialize a new control structure */

	memset ((void *)&C, 0, sizeof (struct SPECTRUM1D_INFO));
	
	n_outputs = 0;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'H':
				case 'V':
				case 'b':
				case 'f':
				case '\0':
					error += (GMT_LONG)GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'C':
					Ctrl->C.active = TRUE;
					if (argv[i][2]) memset ((void *)Ctrl->C.col, 0, SPECTRUM1D_N_OUTPUT_CHOICES * sizeof (char));	/* Reset and read options */
					for (j = 2, k = 0; argv[i][j]; j++, k++) {
						if (k < SPECTRUM1D_N_OUTPUT_CHOICES)
							Ctrl->C.col[k] = argv[i][j];
						else {
							error++;
							fprintf (stderr, "%s: GMT SYNTAX ERROR -C option: Too many output columns selected\n", GMT_program);
							fprintf (stderr, "%s: GMT SYNTAX ERROR -C option: Choose from -Cxycnpago\n", GMT_program);
						}
					}
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					Ctrl->D.inc = atof(&argv[i][2]);
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					if (argv[i][2]) {
						free ((void *)Ctrl->N.name);
						Ctrl->N.name = strdup (&argv[i][2]);
					}
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					Ctrl->S.size = atoi(&argv[i][2]);
					while (window_test < Ctrl->S.size) {
						window_test += window_test;
					}
					break;
				case 'W':
					Ctrl->W.active = TRUE;
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else {
			if ((fp = GMT_fopen(argv[i], GMT_io.r_mode)) == NULL) {
				fprintf(stderr,"%s:  Cannot open r %s\n", GMT_program, argv[i]);
				error = TRUE;
			}
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {	/* Display usage */
		fprintf (stderr, "spectrum1d %s - Compute auto- [and cross- ] spectra from one [or two] timeseries\n\n", GMT_VERSION);
		fprintf(stderr,"usage:  spectrum1d -S<segment_size> [-C[<xycnpago>]] [-D<dt>] [%s] [-N<name_stem>]\n", GMT_H_OPT);
		fprintf(stderr,"\t[-V] [-W] [%s] [%s]\n\n", GMT_b_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf(stderr,"\t-S Use data subsets of <segment_size> elements.\n");
		fprintf(stderr,"\t   <segment_size> must be radix 2;\n");
		fprintf(stderr,"\t   std. err. = 1/sqrt(n_data/segment_size).\n");
		fprintf(stderr,"\tOptions:\n");
		fprintf(stderr,"\t-C[<xycnpago>] 2 column X(t),Y(t) input; estimate Cross-spectra\n\t   [Default 1 col, X power only].\n");
		fprintf(stderr,"\t   Optionally specify cross-spectra output(s)  [Default is all].\n");
		fprintf(stderr,"\t   x = xpower, y = ypower, c = coherent power, n = noise power\n");
		fprintf(stderr,"\t   p = phase, a = admittance, g = gain, o = squared coherency.\n\n");
		fprintf(stderr,"\t-D set delta_time sampling interval of data [Default = 1.0].\n");
		GMT_explain_option ('H');
		fprintf(stderr,"\t-N supply name stem for files [Default = 'spectrum'].\n");
		fprintf(stderr,"\t   Output files will be named <name_stem>.xpower, etc.\n");
		GMT_explain_option ('V');
		fprintf(stderr,"\t-W write Wavelength of spectral estimate in col 1 [Default = frequency].\n");
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf(stderr,"\t   Default is 2 input columns.\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (Ctrl->S.size <= 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option: segment size must be positive\n", GMT_program);
		error++;
	}
	if (window_test != Ctrl->S.size) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option: Segment size not radix 2.  Try %ld or %ld\n", GMT_program,
			(window_test/2), window_test);
		error++;
	}
	if (Ctrl->D.inc <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -D option: Sampling interval must be positive\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = (Ctrl->C.active + 1);
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < (Ctrl->C.active + 1)) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data must have at least %ld columns\n", GMT_program, Ctrl->C.active + 1);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

	C.dt = Ctrl->D.inc;
	C.y_given = Ctrl->C.active;
	C.window = Ctrl->S.size;
	for (k = n_outputs = 0; k < SPECTRUM1D_N_OUTPUT_CHOICES && Ctrl->C.col[k]; k++) {
		if (!strchr ("xycnpago", Ctrl->C.col[k])) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -C option.  Unrecognized output choice %c\n", GMT_program, Ctrl->C.col[k]);
			error++;
		}
		n_outputs++;
	}

	if (!Ctrl->C.active) {		/* ensure x-power output */
		Ctrl->C.col[0] = 'x';
		n_outputs = 1;
		Ctrl->C.active = TRUE;
	}

#ifdef SET_IO_MODE
	GMT_setmode (GMT_OUT);
#endif

	alloc_arrays(&C);

	if (fp == NULL) {
		fp = GMT_stdin;
#ifdef SET_IO_MODE
		GMT_setmode (GMT_IN);
#endif
	}
	n_alloc = GMT_CHUNK;
	n_data = n_read = 0;
	x = (float *) GMT_memory (VNULL, (size_t)n_alloc, sizeof(float), GMT_program);
	if (C.y_given) y = (float *)GMT_memory (VNULL, (size_t)n_alloc, sizeof(float), GMT_program);
	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : C.y_given + 1;

	if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) not_used = GMT_fgets (buffer, BUFSIZ, fp);

	while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

		n_read++;

		if (GMT_io.status & GMT_IO_MISMATCH) {
			fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
			continue;
		}

                x[n_data] = (float)in[GMT_X];
		if (C.y_given) y[n_data] = (float)in[GMT_Y];
                n_data++;

                if (n_data == n_alloc) {
                        n_alloc <<= 1;
			x = (float *) GMT_memory ((void *)x, (size_t)n_alloc, sizeof(float), GMT_program);
			if (C.y_given) y = (float *)GMT_memory ((void *)y, (size_t)n_alloc, sizeof(float), GMT_program);
                }
        }
        if (fp != GMT_stdin) GMT_fclose(fp);
	x = (float *) GMT_memory ((void *)x, (size_t)n_data, sizeof(float), GMT_program);
	if (C.y_given) y = (float *)GMT_memory ((void *)y, (size_t)n_data, sizeof(float), GMT_program);
	if (gmtdefs.verbose) fprintf(stderr,"Read %ld data points.\n", n_data);

	compute_spectra (&C, x, y, n_data);

	GMT_free ((void *)x);
	if (C.y_given) GMT_free ((void *)y);

	write_output (&C, Ctrl->C.col, n_outputs, Ctrl->W.active, Ctrl->N.name);

	free_space (&C);

	Free_spectrum1d_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void alloc_arrays (struct SPECTRUM1D_INFO *C)
{
	C->n_spec = C->window/2;	/* This means we skip zero frequency; data are detrended  */
	C->window_2 = 2 * C->window;		/* This is for complex array stuff  */

	C->spec = (struct SPEC *) GMT_memory (VNULL, (size_t)C->n_spec, sizeof(struct SPEC), GMT_program);
	C->datac = (float *) GMT_memory (VNULL, (size_t)C->window_2, sizeof(float), GMT_program);
}

void compute_spectra (struct SPECTRUM1D_INFO *C, float *x, float *y, GMT_LONG n_data)
{
	GMT_LONG	n_windows, w, i, t_start, t_stop, t, f;
	GMT_LONG	narray;
	float	work = 0.0;
	double	dw, spec_scale, x_varp, y_varp = 1.0, one_on_nw, co_quad;
	double	xreal, ximag, yreal, yimag, xpower, ypower, co_spec, quad_spec;
	char	format[BUFSIZ];

	void detrend_and_hanning (struct SPECTRUM1D_INFO *C);
	
	narray = C->window;

	/* Scale factor for spectral estimates should be 1/4 of amount given in
		Bendat & Piersol eqn 11-102 because I compute 2 * fft in my
		one-sided code below.  However, tests show that I need 1/8 of
		their equation to match variances approximately:  */

	/* This used to read:  spec_scale = 0.5 / (C->window_2 * dt);  */
	spec_scale = C->dt / (C->window_2);

	C->d_n_windows = (double)n_data / (double)C->window;

	n_windows = irint (2.0 * C->d_n_windows) - 1;
	one_on_nw = 1.0 / (double)n_windows;
	dw = (n_windows > 1) ? (double)(n_data - C->window) / (double)(n_windows - 1) : 1.0;

	for (w = 0; w < n_windows; w++) {
		t_start = (GMT_LONG)floor (0.5 + w * dw);
		t_stop = t_start + C->window;
		if (C->y_given) {
			for (t = t_start, i = 0; t < t_stop; t++, i+=2) {
				C->datac[i] = x[t];
				C->datac[i+1] = y[t];
			}
		}
		else {
			for (t = t_start, i = 0; t < t_stop; t++, i+=2) {
				C->datac[i] = x[t];
				C->datac[i+1] = 0.0;
			}
		}

		detrend_and_hanning (C);

		GMT_fourt (C->datac, &narray, 1, -1, 1, &work);

		/* Get one-sided estimates:  */

		x_varp = spec_scale * (C->datac[0] * C->datac[0]);
		if (C->y_given) {
			y_varp = spec_scale * (C->datac[1] * C->datac[1]);
			for (i = 0, f = 2; i < C->n_spec; i++, f+=2) {
				xreal = (i == C->n_spec - 1) ? C->datac[f] : C->datac[f] + C->datac[C->window_2 - f];
				ximag = (i == C->n_spec - 1) ? 0.0 : C->datac[f+1] - C->datac[C->window_2 - f + 1];
				yreal = (i == C->n_spec - 1) ? C->datac[f+1] : C->datac[f+1] + C->datac[C->window_2 - f + 1];
				yimag = (i == C->n_spec - 1) ? 0.0 : C->datac[C->window_2 - f] - C->datac[f];
				xpower = spec_scale * (xreal * xreal + ximag * ximag);
				ypower = spec_scale * (yreal * yreal + yimag * yimag);
				co_spec = spec_scale * (xreal * yreal + ximag * yimag);
				quad_spec = spec_scale * (ximag * yreal - yimag * xreal);

				x_varp += xpower;
				y_varp += ypower;
				C->spec[i].xpow += xpower;
				C->spec[i].ypow += ypower;
				/* Temporarily store co-spec in gain:  */
				C->spec[i].gain += co_spec;
				/* Temporarily store quad-spec in phase:  */
				C->spec[i].phase += quad_spec;
			}
			x_varp *= (C->dt/C->n_spec);
			y_varp *= (C->dt/C->n_spec);
		}
		else {
			for (i = 0, f = 2; i < C->n_spec; i++, f+=2) {
				xreal = C->datac[f] + C->datac[C->window_2 - f];
				ximag = C->datac[f+1] - C->datac[C->window_2 - f + 1];
				xpower = spec_scale * (xreal * xreal + ximag * ximag);
				x_varp += xpower;
				C->spec[i].xpow += xpower;
			}
			x_varp *= (C->dt/C->n_spec);
		}

		if (gmtdefs.verbose) {
			C->y_pow = (C->y_given) ? C->y_variance/y_varp : 0.0;
			fprintf(stderr,"Window %ld from %ld to %ld\n", w, t_start, t_stop);
			sprintf(format, "X var:  %s  X pow:  %s  ratio:  %s  Y var:  %s  Y pow:  %s  ratio:  %s\n",
				gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
			fprintf(stderr, format, C->x_variance, x_varp, (C->x_variance/x_varp), C->y_variance, y_varp, C->y_pow);
		}
	}
	/* Now we can divide by n_windows for the ensemble average.
		The cross spectral stuff needs to be computed:  */

	if (C->y_given ) {
		for (i = 0; i < C->n_spec; i++) {
			C->spec[i].xpow *= one_on_nw;
			C->spec[i].ypow *= one_on_nw;
			co_spec = C->spec[i].gain * one_on_nw;
			quad_spec = C->spec[i].phase * one_on_nw;
			C->spec[i].phase = d_atan2(quad_spec, co_spec);
			co_quad = co_spec * co_spec + quad_spec * quad_spec;
			C->spec[i].coh = co_quad / (C->spec[i].xpow * C->spec[i].ypow);
			C->spec[i].gain = sqrt(co_quad) / C->spec[i].xpow;
			C->spec[i].radmit = co_spec / C->spec[i].xpow;
		}
	}
	else {
		for (i = 0; i < C->n_spec; i++) {
			C->spec[i].xpow *= one_on_nw;
		}
	}
}

void write_output (struct SPECTRUM1D_INFO *C, char *col, GMT_LONG n_outputs, GMT_LONG write_wavelength, char *namestem)
{
	GMT_LONG	i, j;
	double	delta_f, eps_pow, out[3], *f_or_w;
	char	fname[GMT_LONG_TEXT];
	FILE	*fpout;

	delta_f = 1.0 / (C->window * C->dt);
	eps_pow = 1.0 / sqrt(C->d_n_windows);	/* Multiplicative error bars for power spectra  */

	f_or_w = (double *) GMT_memory (VNULL, (size_t)C->n_spec, sizeof (double), GMT_program);
	for (i = 0; i < C->n_spec; i++) f_or_w[i] = (write_wavelength) ? 1.0 / ((i + 1) * delta_f) : (i + 1) * delta_f;

	/* loop through output choices */
	for (j = 0; j < n_outputs; j++) {
		switch (col[j]) {
		case 'x':		/* write x power [ B&P 2nd Ed. eqn. 9.32 ] */
			sprintf(fname, "%s.xpower", namestem);
			if ( (fpout = GMT_fopen (fname, GMT_io.w_mode)) == NULL) {
				fprintf(stderr,"%s:  Cannot open w %s\n", GMT_program, fname);
				exit (EXIT_FAILURE);
			}
			if (gmtdefs.verbose) fprintf(stderr,"%s:  Writing %s\n", GMT_program, fname);
			for (i = 0; i < C->n_spec; i++) {
				out[GMT_X] = f_or_w[i];
				out[GMT_Y] = C->spec[i].xpow;
				out[2] = eps_pow * C->spec[i].xpow;
				GMT_output (fpout, 3, out);
			}
			GMT_fclose(fpout);
			break;

		case 'y':		/* Write y power [ B&P 2nd Ed. eqn. 9.32 ] */
			sprintf(fname, "%s.ypower", namestem);
			if ( (fpout = GMT_fopen(fname, GMT_io.w_mode)) == NULL) {
				fprintf(stderr,"%s:  Cannot open w %s\n", GMT_program, fname);
				exit (EXIT_FAILURE);
			}
			if (gmtdefs.verbose) fprintf(stderr,"%s:  Writing %s\n", GMT_program, fname);
			for (i = 0; i < C->n_spec; i++) {
				out[GMT_X] = f_or_w[i];
				out[GMT_Y] = C->spec[i].ypow;
				out[2] = eps_pow * C->spec[i].ypow;
				GMT_output (fpout, 3, out);
			}
			GMT_fclose(fpout);
			break;
		case 'c':		/* Write Coherent Output power [ B&P 2nd Ed. eqn. 9.71 ] */
			sprintf(fname, "%s.cpower", namestem);
			if ( (fpout = GMT_fopen(fname, GMT_io.w_mode)) == NULL) {
				fprintf(stderr,"%s:  Cannot open w %s\n", GMT_program, fname);
				exit (EXIT_FAILURE);
			}
			if (gmtdefs.verbose) fprintf(stderr,"%s:  Writing %s\n", GMT_program, fname);
			for (i = 0; i < C->n_spec; i++) {
				out[GMT_X] = f_or_w[i];
				out[GMT_Y] = C->spec[i].ypow * C->spec[i].coh;
				out[2] = out[GMT_Y] * eps_pow * sqrt( (2.0 - C->spec[i].coh) / C->spec[i].coh);
				GMT_output (fpout, 3, out);
			}
			GMT_fclose(fpout);
			break;
		case 'n':		/* Write Noise Output power [ B&P 2nd Ed. eqn. 9.73 & Table 9.6 ] */
			sprintf(fname, "%s.npower", namestem);
			if ( (fpout = GMT_fopen(fname, GMT_io.w_mode)) == NULL) {
				fprintf(stderr,"%s:  Cannot open w %s\n", GMT_program, fname);
				exit (EXIT_FAILURE);
			}
			if (gmtdefs.verbose) fprintf(stderr,"%s:  Writing %s\n", GMT_program, fname);
			for (i = 0; i < C->n_spec; i++) {
				out[GMT_X] = f_or_w[i];
				out[GMT_Y] = C->spec[i].ypow * (1.0 - C->spec[i].coh);
				out[2] = out[GMT_Y] * eps_pow;
				GMT_output (fpout, 3, out);
			}
			GMT_fclose(fpout);
			break;
		case 'g':		/* Write Gain spectrum [ B&P 2nd Ed. eqn. 9.90 & Table 9.6 ] */
			sprintf(fname, "%s.gain", namestem);
			if ( (fpout = GMT_fopen(fname, GMT_io.w_mode)) == NULL) {
				fprintf(stderr,"%s:  Cannot open w %s\n", GMT_program, fname);
				exit (EXIT_FAILURE);
			}
			if (gmtdefs.verbose) fprintf(stderr,"%s:  Writing %s\n", GMT_program, fname);
			for (i = 0; i < C->n_spec; i++) {
				out[GMT_X] = f_or_w[i];
				out[GMT_Y] = C->spec[i].gain;
				out[2] = out[GMT_Y] * eps_pow * sqrt( (1.0 - C->spec[i].coh) / (2.0 * C->spec[i].coh) );
				GMT_output (fpout, 3, out);
			}
			GMT_fclose(fpout);
			break;
		case 'a':		/* Write Real Admittance spectrum
			We don't know the correct error estimate and it is not
			in Bendat and Piersol, or Priestly, or any standard text.
			Smith needs to derive this, and should make a note to
			check the expression given by Marcia Maia et al in Geophys. 
			J. Int, 100, 337-348, 1990, equation 10, page 341.
			Meanwhile we will default to use the expression related to
			that for the gain spectrum:
		*/
			sprintf(fname, "%s.admit", namestem);
			if ( (fpout = GMT_fopen(fname, "w")) == NULL) {
				fprintf(stderr,"%s:  Cannot open w %s\n", GMT_program, fname);
				exit (EXIT_FAILURE);
			}
			if (gmtdefs.verbose) fprintf(stderr,"%s:  Writing %s\n", GMT_program, fname);
			for (i = 0; i < C->n_spec; i++) {
				out[GMT_X] = f_or_w[i];
				out[GMT_Y] = C->spec[i].radmit;
				out[2] = fabs (eps_pow * sqrt( (1.0 - C->spec[i].coh) / (2.0 * C->spec[i].coh) ) * out[GMT_Y]);
				GMT_output (fpout, 3, out);
			}
			GMT_fclose(fpout);
			break;
		case 'p':		/* Write Phase spectrum [ B&P 2nd Ed. eqn. 9.91 & Table 9.6 ] */
			sprintf(fname, "%s.phase", namestem);
			if ( (fpout = GMT_fopen(fname, GMT_io.w_mode)) == NULL) {
				fprintf(stderr,"%s:  Cannot open w %s\n", GMT_program, fname);
				exit (EXIT_FAILURE);
			}
			if (gmtdefs.verbose) fprintf(stderr,"%s:  Writing %s\n", GMT_program, fname);
			for (i = 0; i < C->n_spec; i++) {
				out[GMT_X] = f_or_w[i];
				out[GMT_Y] = C->spec[i].phase;
				out[2] = eps_pow * sqrt( (1.0 - C->spec[i].coh) / (2.0 * C->spec[i].coh) );
				GMT_output (fpout, 3, out);
			}
			GMT_fclose(fpout);
		case 'o':		/* Write Coherency spectrum [ B&P 2nd Ed. eqn. 9.82 ] */
			sprintf(fname, "%s.coh", namestem);
			if ( (fpout = GMT_fopen(fname, GMT_io.w_mode)) == NULL) {
				fprintf(stderr,"%s:  Cannot open w %s\n", GMT_program, fname);
				exit (EXIT_FAILURE);
			}
			if (gmtdefs.verbose) fprintf(stderr,"%s:  Writing %s\n", GMT_program, fname);
			for (i = 0; i < C->n_spec; i++) {
				out[GMT_X] = f_or_w[i];
				out[GMT_Y] = C->spec[i].coh;
				out[2] = out[GMT_Y] * eps_pow * (1.0 - C->spec[i].coh) * sqrt(2.0 / C->spec[i].coh);
				GMT_output (fpout, 3, out);
			}
			GMT_fclose(fpout);
			break;
		}
	}

	GMT_free ((void *)f_or_w);
}

void free_space (struct SPECTRUM1D_INFO *C)
{
	GMT_free ((void *)C->spec);
	GMT_free ((void *)C->datac);
}

void detrend_and_hanning (struct SPECTRUM1D_INFO *C)
{
	GMT_LONG	i, t;
	double	sumx, sumtx, sumy, sumty, sumt2, x_slope, x_mean, y_slope, y_mean;
	double	t_factor, h_period, h_scale, hc, hw, tt;
	sumx = 0.0;
	sumtx = 0.0;
	sumy = 0.0;
	sumty = 0.0;
	sumt2 = 0.0;
	C->x_variance = 0.0;
	C->y_variance = 0.0;
	t_factor = 2.0 / (C->window - 1);
	h_period = M_PI / (double)C->window;	/* For Hanning window  */
	h_scale = sqrt(8.0/3.0);		/* For Hanning window  */

	if (C->y_given) {
		for (i = 0, t = 0; i < C->window_2; i+=2, t++) {
			tt = t * t_factor - 1.0;
			sumt2 += (tt * tt);
			sumx += C->datac[i];
			sumtx += (tt * C->datac[i]);
			sumy += C->datac[i+1];
			sumty += (tt * C->datac[i+1]);
		}
	}
	else {
		for (i = 0, t = 0; i < C->window_2; i+=2, t++) {
			tt = t * t_factor - 1.0;
			sumt2 += (tt * tt);
			sumx += C->datac[i];
			sumtx += (tt * C->datac[i]);
		}
	}
	x_slope = sumtx / sumt2;
	x_mean = sumx / C->window;
	if (C->y_given) {
		y_slope = sumty / sumt2;
		y_mean = sumy / C->window;
		for (i = 0, t = 0; i < C->window_2; i+=2, t++) {
			hc = cos(t * h_period);
			hw = h_scale * (1.0 - hc * hc);
			tt = t * t_factor - 1.0;
			C->datac[i] -= (float)(x_mean + tt * x_slope);
			C->datac[i] *= (float)hw;
			C->x_variance += (C->datac[i] * C->datac[i]);
			C->datac[i+1] -= (float)(y_mean + tt * y_slope);
			C->datac[i+1] *= (float)hw;
			C->y_variance += (C->datac[i+1] * C->datac[i+1]);
		}
		C->x_variance /= C->window;
		C->y_variance /= C->window;
	}
	else {
		for (i = 0, t = 0; i < C->window_2; i+=2, t++) {
			hc = cos(t * h_period);
			hw = h_scale * (1.0 - hc * hc);
			tt = t * t_factor - 1.0;
			C->datac[i] -= (float)(x_mean + tt * x_slope);
			C->datac[i] *= (float)hw;
			C->x_variance += (C->datac[i] * C->datac[i]);
		}
		C->x_variance /= C->window;
	}
}
void *New_spectrum1d_Ctrl () {	/* Allocate and initialize a new control structure */
	struct SPECTRUM1D_CTRL *C;
	
	C = (struct SPECTRUM1D_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct SPECTRUM1D_CTRL), "New_spectrum1d_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->D.inc = 1.0;
	C->C.col[0] = 'x';
	C->C.col[1] = 'y';
	C->C.col[2] = 'c';
	C->C.col[3] = 'n';
	C->C.col[4] = 'p';
	C->C.col[5] = 'a';
	C->C.col[6] = 'g';
	C->C.col[7] = 'o';
	C->N.name = strdup ("spectrum");
	return ((void *)C);
}

void Free_spectrum1d_Ctrl (struct SPECTRUM1D_CTRL *C) {	/* Deallocate control structure */
	if (C->N.name) free ((void *)C->N.name);	
	GMT_free ((void *)C);	
}
