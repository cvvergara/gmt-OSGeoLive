/*--------------------------------------------------------------------
 *	$Id: splitxyz.c 9923 2012-12-18 20:45:53Z pwessel $
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
 * 
 * read a file of lon, lat, zvalue[, distance, azimuth]
 * and split it into profile segments.
 * 
 * Author:	W. H. F. Smith
 * Date:	24 July, 1991-2000
 * Version:	4
 * 4.1.2:	Added -Q option to select outut columns [PW. 08-MAR-2006]
 */

#include "gmt.h"

#define SPLITXYZ_F_RES 1000	/* Number of points in filter halfwidth  */
#define SPLITXYZ_N_OUTPUT_CHOICES 5

struct SPLITXYZ_CTRL {
	struct A {	/* -A<azimuth>/<tolerance> */
		GMT_LONG active;
		double azimuth, tolerance;
	} A;
	struct C {	/* -C<course_change> */
		GMT_LONG active;
		double value;
	} C;
	struct D {	/* -D<mindist> */
		GMT_LONG active;
		double value;
	} D;
	struct F {	/* -F<xy_filter>/<z_filter> */
		GMT_LONG active;
		double xy_filter, z_filter;
	} F;
	struct G {	/* -D<gap> */
		GMT_LONG active;
		double value;
	} G;
	struct M {	/* -M */
		GMT_LONG active;
	} M;
	struct N {	/* -N<namestem> */
		GMT_LONG active;
		char *name;
	} N;
	struct Q {	/* -Q[<xyzdg>] */
		GMT_LONG active;
		char col[SPLITXYZ_N_OUTPUT_CHOICES];	/* Character codes for desired output in the right order */
	} Q;
	struct S {	/* -S */
		GMT_LONG active;
	} S;
	struct Z{	/* -Z */
		GMT_LONG active;
	} Z;
};

struct  DATA    {
	double val[SPLITXYZ_N_OUTPUT_CHOICES];	/* 0 = x, 1 = y, 2 = z, 3 = d, 4 = a */
        double  w;
};

int main(int argc, char **argv)
{
	GMT_LONG     i, ndata, n_alloc, nprofiles, begin, end;
	GMT_LONG	j, k, n_required_fields, output_choice[SPLITXYZ_N_OUTPUT_CHOICES], n_expected_fields, n_outputs = 0, n_read, n_fields;

	GMT_LONG	error = FALSE, ok, hilow = FALSE, z_selected = FALSE;

	double	dy, dx, last_d, last_c, last_s, csum, ssum, this_c, this_s, dotprod;
	double	mean_azim, fwork[SPLITXYZ_F_RES], *in = NULL, out[SPLITXYZ_N_OUTPUT_CHOICES];

 	char    buffer[BUFSIZ], filename[BUFSIZ];

	struct  DATA *data = NULL;
	struct SPLITXYZ_CTRL *Ctrl = NULL;

	FILE    *fp = NULL, *fpout = NULL;

	void	filterz(struct DATA *data, GMT_LONG ndata, double z_filter, double *fwork, GMT_LONG hilow), filterxy_setup(double *fwork), filterxy(struct DATA *data, GMT_LONG begin, GMT_LONG end, double xy_filter, double *fwork, GMT_LONG hilow);
	void *New_splitxyz_Ctrl (), Free_splitxyz_Ctrl (struct SPLITXYZ_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);
	
	Ctrl = (struct SPLITXYZ_CTRL *)New_splitxyz_Ctrl ();	/* Allocate and initialize a new control structure */

	fpout = GMT_stdout;
	memset ((void *)output_choice, 0, SPLITXYZ_N_OUTPUT_CHOICES * sizeof (GMT_LONG));

        for (i = 1; i < argc; i++) {
                if (argv[i][0] == '-') {
                        switch (argv[i][1]) {
 
				/* Common parameters */

				case 'H':
				case 'V':
				case ':':
				case 'b':
				case 'f':
				case '\0':
					error += (GMT_LONG)GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'A':
					Ctrl->A.active = TRUE;
                                        if ( (sscanf(&argv[i][2], "%lf/%lf", &Ctrl->A.azimuth, &Ctrl->A.tolerance)) != 2) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -A option: Can't decipher values\n", GMT_program);
						error = TRUE;
					}
                                        break;
				case 'C':
					Ctrl->C.active = TRUE;
                                        if ( (sscanf(&argv[i][2], "%lf", &Ctrl->C.value)) != 1) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -C option: Can't decipher value\n", GMT_program);
						error = TRUE;
					}
                                        break;
                                case 'D':
					Ctrl->D.active = TRUE;
                                        if ( (sscanf(&argv[i][2], "%lf", &Ctrl->D.value)) != 1) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -D option: Can't decipher value\n", GMT_program);
						error = TRUE;
					}
                                        break;
				case 'F':
					Ctrl->F.active = TRUE;
                                        if ( (sscanf(&argv[i][2], "%lf/%lf", &Ctrl->F.xy_filter, &Ctrl->F.z_filter)) != 2) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -F option: Can't decipher values\n", GMT_program);
						error = TRUE;
					}
                                        break;
                                case 'G':
					Ctrl->G.active = TRUE;
                                        if ( (sscanf(&argv[i][2], "%lf", &Ctrl->G.value)) != 1) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -G option: Can't decipher value\n", GMT_program);
						error = TRUE;
					}
                                        break;
				case 'M':
                                        Ctrl->M.active = TRUE;
                                        break;
				case 'N':
					Ctrl->N.active = TRUE;
 					if (argv[i][2])
						Ctrl->N.name = strdup (&argv[i][2]);
					else {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -N option: Append a name stem\n", GMT_program);
						error = TRUE;
					}
                                        break;
				case 'Q':
					Ctrl->Q.active = TRUE;
					for (j = 2, k = 0; argv[i][j]; j++, k++) {
						if (k < SPLITXYZ_N_OUTPUT_CHOICES)
							Ctrl->Q.col[k] = argv[i][j];
						else {
							error++;
							fprintf (stderr, "%s: GMT SYNTAX ERROR -Q option: Too many output columns selected\n", GMT_program);
							fprintf (stderr, "%s: GMT SYNTAX ERROR -Q option: Choose from -Qxyzdg\n", GMT_program);
						}
					}
					break;
				case 'S':
                                        Ctrl->S.active = TRUE;
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
                else {
                        if ( (fp = GMT_fopen(argv[i], GMT_io.r_mode)) == NULL) {
                                fprintf(stderr,"%s:  Cannot open %s\n", GMT_program, argv[i]);
                                error = TRUE;
                        }
                }
        }

	if (argc == 1 || GMT_give_synopsis_and_exit) {	/* Display usage */
 		fprintf (stderr, "splitxyz %s - Split xyz[dh] files into segments\n\n", GMT_VERSION);
		fprintf(stderr,"usage:  splitxyz [<xyz[dh]file>] -C<course_change> [-A<azimuth>/<tolerance>]\n");
		fprintf(stderr,"\t[-D<minimum_distance>] [-F<xy_filter>/<z_filter>] [-G<gap>] [%s] [-M]\n", GMT_H_OPT);
		fprintf(stderr,"\t[-N<namestem>] [-Q<flags>] [-S] [-V] [-Z] [-:] [%s] [%s]\n\n", GMT_b_OPT, GMT_f_OPT);
 
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf(stderr,"\tGive xyz[dh]file name or read stdin.\n");
		fprintf(stderr,"\t-C  Profile ends when change of heading exceeds <course_change>.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf(stderr,"\t-A Only write profile if mean direction is w/in +/- <tolerance>\n");
		fprintf(stderr,"\t     of <azimuth>. [Default = All].\n");
		fprintf(stderr,"\t-D Only write profile if length is at least <minimum_distance> [0].\n");
		fprintf(stderr,"\t-F Filter the data.  Give full widths of cosine arch filters for xy and z.\n");
		fprintf(stderr,"\t   Defaults are both widths = 0, giving no filtering.\n");
		fprintf(stderr,"\t   Use negative width to highpass.\n");
		fprintf(stderr,"\t-G Do not let profiles have gaps exceeding <gap>. [Default = 10 dist units].\n");
		GMT_explain_option ('H');
		fprintf(stderr,"\t-M Map units TRUE; x,y in degrees, dist units in km.  [Default dist unit = x,y unit].\n");
		fprintf(stderr,"\t-N Write output to separate files named <namestem>.profile#.\n");
		fprintf(stderr,"\t     [Default all to stdout, separated by >].\n");
		fprintf(stderr,"\t-Q Indicate what output you want as one or more of xyzdh in any order;\n");
		fprintf(stderr,"\t  where x,y,z refer to input data locations and optional z-value(s),\n");
		fprintf(stderr,"\t  and d,h are the distance and heading along track.\n");
		fprintf(stderr,"\t  [Default is all fields, i.e., -Qxyzdh (or -Qxydh if -Z is set)].\n");
		fprintf(stderr,"\t-S d,h is supplied.  Input is 5 col x,y,z,d,h with d non-decreasing.\n");
		fprintf(stderr,"\t     [Default input is 3 col x,y,z only and computes d,h from the data].\n");
		fprintf(stderr,"\t-Z No z-values.  Input is 2 col x,y only.\n");
		GMT_explain_option ('V');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf(stderr,"\t     Default input columns is set given -S and -Z options.\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
        }

	if (Ctrl->D.value < 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -D option: Minimum segment distance must be positive\n", GMT_program);
		error++;
	}
	if (Ctrl->C.value <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -C option: Course change tolerance must be positive\n", GMT_program);
		error++;
	}
	if (Ctrl->A.tolerance < 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -A option: Azimuth tolerance must be positive\n", GMT_program);
		error++;
	}
	if (Ctrl->G.value < 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option: Data gap distance must be positive\n", GMT_program);
		error++;
	}
	if (Ctrl->Z.active && Ctrl->S.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Z option: Cannot be used with -S option\n", GMT_program);
		error++;
	}
	if (Ctrl->Z.active && Ctrl->F.z_filter != 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -F option: Cannot specify z-filter while using -Z option\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_OUT] && !Ctrl->N.name) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary output requires a namestem in -N\n", GMT_program);
		error++;
	}
	if (n_outputs > 0 && z_selected && Ctrl->Z.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  -Q:  Cannot request z if -Z have been specified\n", GMT_program);
		error++;
	}
	n_required_fields = (Ctrl->S.active) ? 5 : ((Ctrl->Z.active) ? 2 : 3);
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = n_required_fields;
	if (GMT_io.binary[GMT_IN])
		n_expected_fields = GMT_io.ncol[GMT_IN];
	else
		n_expected_fields = n_required_fields;
	if (GMT_io.binary[GMT_IN] && n_expected_fields < n_required_fields) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data must have at least %ld columns\n", GMT_program, n_required_fields);
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

	for (k = n_outputs = 0; k < SPLITXYZ_N_OUTPUT_CHOICES && Ctrl->Q.col[k]; k++) {
		switch (Ctrl->Q.col[k]) {
			case 'x':
				output_choice[k] = 0;
				break;
			case 'y':
				output_choice[k] = 1;
				break;
			case 'z':
				output_choice[k] = 2;
				z_selected = TRUE;
				break;
			case 'd':
				output_choice[k] = 3;
				break;
			case 'h':
				output_choice[k] = 4;
				break;
			default:
				fprintf (stderr, "%s: GMT SYNTAX ERROR -Q option:  Unrecognized choice %c\n", GMT_program, Ctrl->Q.col[k]);
				error = TRUE;
		}
		n_outputs++;
	}
	if (Ctrl->M.active) {
		GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
		GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
	}
	if (n_outputs == 0) {	/* Generate default -Q setting (all) */
		n_outputs = 5 - Ctrl->Z.active;
		for (i = 0; i < 2; i++) output_choice[i] = (int)i;
		for (i = 2; i < n_outputs; i++) output_choice[i] = (int)i + Ctrl->Z.active;

	}

	Ctrl->A.tolerance *= D2R;
	/* if (Ctrl->A.azimuth > 180.0) Ctrl->A.azimuth -= 180.0; */	/* Put in Easterly strikes  */
	Ctrl->A.azimuth = D2R * (90.0 - Ctrl->A.azimuth);	/* Work in cartesian angle and radians  */
	Ctrl->C.value *= D2R;

	if (fp == NULL) {
		fp = GMT_stdin;
#ifdef SET_IO_MODE
		GMT_setmode (GMT_IN);
#endif
	}
        n_alloc = GMT_CHUNK;
	ndata = n_read = 0;

        data = (struct DATA *) GMT_memory (VNULL, (size_t)n_alloc, sizeof(struct DATA), GMT_program);

	if (!GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) GMT_fgets (buffer, BUFSIZ, fp);

	while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

		n_read++;
		if (GMT_io.status & GMT_IO_MISMATCH) {
			fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields,  n_expected_fields, n_read);
			continue;

		}
                if (ndata == n_alloc) {
                        n_alloc <<= 1;
                        data = (struct DATA *) GMT_memory ((void *)data, (size_t)n_alloc, sizeof(struct DATA), GMT_program);
                }
		if (!Ctrl->S.active) {
			data[ndata].val[0] = in[GMT_X];
			data[ndata].val[1] = in[GMT_Y];
			if (!Ctrl->Z.active) data[ndata].val[2] = in[GMT_Z];
			if (ndata) {
				dy = (data[ndata].val[1] - data[ndata-1].val[1]);
				dx = (data[ndata].val[0] - data[ndata-1].val[0]);
				if (Ctrl->M.active) {
					dy *= project_info.DIST_KM_PR_DEG;
					dx *= (project_info.DIST_KM_PR_DEG * cosd (0.5 * (data[ndata].val[1] + data[ndata-1].val[1])));
				}
				if (dy == 0.0 && dx == 0.0) {
					data[ndata].val[3] = data[ndata-1].val[3];
					data[ndata].val[4] = data[ndata-1].val[4];
				}
				else {
					data[ndata].val[3] = data[ndata-1].val[3] + hypot(dx,dy);
					data[ndata].val[4] = d_atan2(dy,dx);	/* Angles are stored as CCW angles in radians */
				}
			}
			else {
				data[ndata].val[3] = data[ndata].val[4] = 0.0;
			}
			ndata++;
		}
		else {
			data[ndata].val[0] = in[GMT_X];
			data[ndata].val[1] = in[GMT_Y];
			data[ndata].val[2] = in[GMT_Z];
			data[ndata].val[3] = in[3];
			data[ndata].val[4] = D2R * (90.0 - in[4]);	/* Angles are stored as CCW angles in radians */
			ndata++;
		}
	}

	data = (struct DATA *) GMT_memory ((void *)data, (size_t)ndata, sizeof(struct DATA), GMT_program);
	if (!Ctrl->S.active) data[0].val[4] = data[1].val[4];
	if (fp != GMT_stdin) GMT_fclose(fp);

	/* Now we have read the data and can filter z if necessary.  */

	if (Ctrl->F.z_filter < 0.0) {
		hilow = TRUE;
		Ctrl->F.z_filter = fabs(Ctrl->F.z_filter);
		filterz(data, ndata, Ctrl->F.z_filter, fwork, hilow);
	}
	else if (Ctrl->F.z_filter > 0.0) {
		hilow = FALSE;
		filterz(data, ndata, Ctrl->F.z_filter, fwork, hilow);
	}

	if (Ctrl->F.xy_filter < 0.0) {
		hilow = TRUE;
		Ctrl->F.xy_filter = fabs(Ctrl->F.xy_filter);
		filterxy_setup(fwork);
	}
	else if (Ctrl->F.xy_filter > 0.0) {
		hilow = FALSE;
		filterxy_setup(fwork);
	}

	/* Now we are ready to search for segments.  */

	nprofiles = 0;
	begin = 0;
	end = begin;
	while (end < ndata-1) {
		last_d = data[begin].val[3];
		sincos (data[begin].val[4], &last_s, &last_c);
		csum = last_c;
		ssum = last_s;
		ok = TRUE;
		while (ok && end < ndata-1) {
			end++;
			if (data[end].val[3] - last_d > Ctrl->G.value) {
				/* Fails due to too much distance gap  */
				ok = FALSE;
				continue;
			}
			sincos (data[end].val[4], &this_s, &this_c);
			dotprod = this_c * last_c + this_s * last_s;
			if (fabs(dotprod) > 1.0) dotprod = copysign(1.0, dotprod);
			if (d_acos(dotprod) > Ctrl->C.value) {
				/* Fails due to too much change in azimuth  */
				ok = FALSE;
				continue;
			}
			/* Get here when this point belongs with last one:  */
			csum += this_c;
			ssum += this_s;
			last_c = this_c;
			last_s = this_s;
			last_d = data[end].val[3];
		}

		/* Get here when we have found a beginning and end  */

		if (ok) {
			/* Last point in input should be included in this segment  */
			end++;
		}

		if (end - begin - 1) {

			/* There are at least two points in the list.  */

			if ((data[end-1].val[3] - data[begin].val[3]) >= Ctrl->D.value) {

				/* List is long enough.  Check strike. Compute mean_azim in range [-pi/2, pi/2]:  */

				mean_azim = d_atan2 (ssum, csum);
				mean_azim = fabs(mean_azim - Ctrl->A.azimuth);
				if (mean_azim <= Ctrl->A.tolerance) {

					/* List has acceptable strike.  */

					if (Ctrl->F.xy_filter != 0.0) filterxy(data, begin, end, Ctrl->F.xy_filter, fwork, hilow);
					nprofiles++;

					/* If (Ctrl->N.name) we write separate files, else GMT_stdout with > marks:  */

					if (Ctrl->N.name) {
						sprintf(filename, "%s.profile%ld", Ctrl->N.name, nprofiles);
						if ( (fpout = GMT_fopen(filename, GMT_io.w_mode)) == NULL) {
							fprintf(stderr,"%s:  Cannot create %s\n", GMT_program, argv[i]);
							exit (EXIT_FAILURE);
						}
					}
					else {
						sprintf (buffer, "> Start profile # %ld\n", nprofiles);
						GMT_fputs (buffer, GMT_stdout);
					}

					for (i = begin; i < end; i++) {
						for (j = 0; j < n_outputs; j++) {	/* Remember to convert CCW angles back to azimuths */
							out[j] = (output_choice[j] == 4) ? 90.0 - R2D * data[i].val[4] : data[i].val[output_choice[j]];
						}
						GMT_output (fpout, n_outputs, out);
					}

					if (fpout != GMT_stdout) {
						GMT_fclose (fpout);
						if (gmtdefs.verbose) fprintf (stderr,"%s: Wrote %ld points to file %s\n", GMT_program, end-begin, filename);
					}
				}
			}
		}
		begin = end;
	}

	/* Get here when all profiles have been found and written.  */

	if (gmtdefs.verbose) fprintf(stderr,"%s:  Split %ld data into %ld files.\n", GMT_program, ndata, nprofiles);
	GMT_free ((void *)data);

	Free_splitxyz_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}


void filterz (struct DATA *data, GMT_LONG ndata, double z_filter, double *fwork, GMT_LONG hilow)
{

	double	half_width, sum, dt, tmp;
	GMT_LONG	i, j, k, istart, istop;

	half_width = 0.5 * z_filter;
	sum = 0.0;
	dt = SPLITXYZ_F_RES/half_width;
	tmp = M_PI / SPLITXYZ_F_RES;
	for (i = 0; i < SPLITXYZ_F_RES; i++) {
		fwork[i] = 1.0 + cos(i*tmp);
		sum += fwork[i];
	}
	for (i = 1; i < SPLITXYZ_F_RES; i++) {
		fwork[i] /= sum;
	}

	j = 0;
	istart = 0;
	istop = 0;
	while(j < ndata) {
		while(istart < ndata && data[istart].val[3] - data[j].val[3] <= -half_width) istart++;
		while(istop < ndata && data[istop].val[3] - data[j].val[3] < half_width) istop++;
		sum = 0.0;
		data[j].w = 0.0;
		for(i = istart; i < istop; i++) {
			k = (GMT_LONG)floor(dt*fabs(data[i].val[3] - data[j].val[3]));
			sum += fwork[k];
			data[j].w += (data[i].val[2] * fwork[k]);
		}
		data[j].w /= sum;
		j++;
	}
	if (hilow) {
		for (i = 0; i < ndata; i++) data[i].val[2] = data[i].val[2] - data[i].w;
	}
	else {
		for (i = 0; i < ndata; i++) data[i].val[2] = data[i].w;
	}
	return;
}

void filterxy_setup (double *fwork)
{
	GMT_LONG	i;
	double	tmp, sum = 0.0;

	tmp = M_PI / SPLITXYZ_F_RES;
	for (i = 0; i < SPLITXYZ_F_RES; i++) {
		fwork[i] = 1.0 + cos(i*tmp);
		sum += fwork[i];
	}
	for (i = 1; i < SPLITXYZ_F_RES; i++) {
		fwork[i] /= sum;
	}
	return;
}

void filterxy (struct DATA *data, GMT_LONG begin, GMT_LONG end, double xy_filter, double *fwork, GMT_LONG hilow)
{
	GMT_LONG	i, j, k, istart, istop;
	double	half_width, dt, sum;

	half_width = 0.5 * fabs(xy_filter);
	dt = SPLITXYZ_F_RES/half_width;

	j = begin;
	istart = begin;
	istop = begin;
	while(j < end) {
		while(istart < end && data[istart].val[3] - data[j].val[3] <= -half_width) istart++;
		while(istop < end && data[istop].val[3] - data[j].val[3] < half_width) istop++;
		sum = 0.0;
		data[j].w = 0.0;
		for(i = istart; i < istop; i++) {
			k = (GMT_LONG)floor(dt*fabs(data[i].val[3] - data[j].val[3]));
			sum += fwork[k];
			data[j].w += (data[i].val[0] * fwork[k]);
		}
		data[j].w /= sum;
		j++;
	}
	if (hilow) {
		for (i = begin; i < end; i++) data[i].val[0] = data[i].val[0] - data[i].w;
	}
	else {
		for (i = begin; i < end; i++) data[i].val[0] = data[i].w;
	}

	j = begin;
	istart = begin;
	istop = begin;
	while(j < end) {
		while(istart < end && data[istart].val[3] - data[j].val[3] <= -half_width) istart++;
		while(istop < end && data[istop].val[3] - data[j].val[3] < half_width) istop++;
		sum = 0.0;
		data[j].w = 0.0;
		for(i = istart; i < istop; i++) {
			k = (GMT_LONG)floor(dt*fabs(data[i].val[3] - data[j].val[3]));
			sum += fwork[k];
			data[j].w += (data[i].val[1] * fwork[k]);
		}
		data[j].w /= sum;
		j++;
	}
	if (hilow) {
		for (i = begin; i < end; i++) data[i].val[1] = data[i].val[1] - data[i].w;
	}
	else {
		for (i = begin; i < end; i++) data[i].val[1] = data[i].w;
	}

	return;
}

void *New_splitxyz_Ctrl () {	/* Allocate and initialize a new control structure */
	struct SPLITXYZ_CTRL *C;
	
	C = (struct SPLITXYZ_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct SPLITXYZ_CTRL), "New_splitxyz_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
        C->A.azimuth = 90.0;
	C->A.tolerance = 360.0;
	C->G.value = DBL_MAX;	
	return ((void *)C);
}

void Free_splitxyz_Ctrl (struct SPLITXYZ_CTRL *C) {	/* Deallocate control structure */
	if (C->N.name) free ((void *)C->N.name);	
	GMT_free ((void *)C);	
}
