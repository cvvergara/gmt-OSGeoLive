/*--------------------------------------------------------------------
 *	$Id: grdgradient.c 9923 2012-12-18 20:45:53Z pwessel $
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
 *  grdgradient.c
 * read a grid file and compute gradient in azim direction:
 *
 * azim = azimuth clockwise from north in degrees.
 *
 * gradient = -[(dz/dx)sin(azim) + (dz/dy)cos(azim)].
 *
 * the expression in [] is the correct gradient.  We take
 * -[]  in order that data which goes DOWNHILL in the
 * azim direction will give a positive value; this is
 * for image shading purposes.
 *
 *
 * Author:	W.H.F. Smith
 * Date: 	13 Feb 1991
 * Upgraded to v2.0 15-May-1991 Paul Wessel
 *
 * Modified:	1 Mar 94 by WHFS to make -M scale change with j latitude
 *		1 Mar 96 by PW to find gradient direction and magnitude (-S and -D)
 *		13 Mar 96 by WHFS to add exp trans and user-supplied sigma to -N
 *			option, and add optional second azimuth to -A option.
 *		11 Sep 97 by PW now may pass average gradient along with sigma in -N
 *		22 Apr 98 by WHFS to add boundary conditions, switch sense of -S and 
 *			-D, and switch -Da to -Dc, for consistency of args.
 *		6  Sep 05 by J. Luis, added a -E option that allows the Lambertian or
 *			Peucker piecewise linear radiance computations
 * Version:	4
 */
 
#include "gmt.h"

struct GRDGRADIENT_CTRL {
	struct A {	/* -A<azim>[/<azim2>] */
		GMT_LONG active;
		GMT_LONG two;
		double azimuth[2];
	} A;
	struct D {	/* -D[a][o][n] */
		GMT_LONG active;
		GMT_LONG mode;
	} D;
	struct E {	/* -E[s|p]<azim>/<elev[ambient/diffuse/specular/shine]> */
		GMT_LONG active;
		double azimuth, elevation;
		double ambient, diffuse, specular, shine;
		GMT_LONG mode;
	} E;
	struct G {	/* -G<file> */
		GMT_LONG active;
		char *file;
	} G;
	struct L {	/* -L<flag> */
		GMT_LONG active;
		char mode[4];
	} L;
	struct M {	/* -M */
		GMT_LONG active;
	} M;
	struct N {	/* -N[t_or_e][<amp>[/<sigma>[/<offset>]]] */
		GMT_LONG active;
		GMT_LONG mode;
		double norm, sigma, offset;
	} N;
	struct S {	/* -S<slopefile> */
		GMT_LONG active;
		char *file;
	} S;
};

EXTERN_MSC GMT_LONG GMT_parse_f_option (char *arg);	/* Needed for -L setup */

int main (int argc, char **argv)
{
	char	*infile = CNULL, format[BUFSIZ], ptr[BUFSIZ];
	
	GMT_LONG	error = FALSE, sigma_set = FALSE, offset_set = FALSE, bad;
	
	GMT_LONG pos, p[4], entry, mx, my;
	
	GMT_LONG	i, j, ij, k, n, nm, nm2, n_used = 0;
	
	float *data = NULL, *slp = VNULL;

	double	dx_grid, dy_grid, x_factor, y_factor, dzdx, dzdy, ave_gradient;
	double	azim, denom, max_gradient = 0.0, min_gradient = 0.0, rpi, lat;
	double	x_factor2 = 0.0, y_factor2 = 0.0, dzdx2 = 0.0, dzdy2 = 0.0, dzds1, dzds2;
	double	p0 = 0.0, q0 = 0.0, p0q0_cte = 1.0, norm_z, mag, s[3], lim_x, lim_y, lim_z;
	double	k_ads = 0.0, diffuse, spec, r_min = DBL_MAX, r_max = -DBL_MAX, scale;
	
	struct	GRD_HEADER header;
	struct	GMT_EDGEINFO edgeinfo;
	struct GRDGRADIENT_CTRL *Ctrl = NULL;
	
	double specular (double nx, double ny, double nz, double *s);
	void *New_grdgradient_Ctrl (), Free_grdgradient_Ctrl (struct GRDGRADIENT_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDGRADIENT_CTRL *)New_grdgradient_Ctrl ();	/* Allocate and initialize a new control structure */
	
	GMT_boundcond_init (&edgeinfo);
	memset ((void *)s, 0, 3*sizeof(double));

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
					Ctrl->A.active = TRUE;
					j = sscanf(&argv[i][2], "%lf/%lf", &Ctrl->A.azimuth[0], &Ctrl->A.azimuth[1]);
					Ctrl->A.two = (j == 2);
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					j = 2;
					while (argv[i][j]) {
						switch (argv[i][j]) {
							case 'C':
							case 'c':
								Ctrl->D.mode |= 1;
								break;
							case 'O':
							case 'o':
								Ctrl->D.mode |= 2;
								break;
							case 'N':
							case 'n':
								Ctrl->D.mode |= 4;
								break;
							default:
								fprintf (stderr, "%s: GMT SYNTAX ERROR -D option:  Unrecognized modifier\n", GMT_program);
								error++;
								break;
						}
						j++;
					}
					break;
				case 'E':	/* Lambertian family radiance */
					Ctrl->E.active = TRUE;
					switch (argv[i][2]) {
						case 'p':	/* Peucker */
							Ctrl->E.mode = 1;
							break;
						case 's':	/* "simple" Lambertian case */
							Ctrl->E.mode = 2;
							if (sscanf(&argv[i][3], "%lf/%lf", &Ctrl->E.azimuth, &Ctrl->E.elevation) != 2) {
								fprintf(stderr,"%s: GMT SYNTAX ERROR -Es option: Must append azimuth/elevation\n", GMT_program);
								error++;
							}
							break;
						default:
							Ctrl->E.mode = 3;	/* "full" Lambertian case */
							if (sscanf(&argv[i][2], "%lf/%lf", &Ctrl->E.azimuth, &Ctrl->E.elevation) < 2) {
								fprintf(stderr,"%s: GMT SYNTAX ERROR -E option: Must give at least azimuth and elevation\n", GMT_program);
								error++;
							}
							entry = pos = 0;
							while (entry < 6 && (GMT_strtok (&argv[i][2], "/", &pos, ptr))) {
								switch (entry) {
									case 0:
									case 1:
										break;	/* Cases already processed above */
									case 2:
										if (ptr[0] != '=') Ctrl->E.ambient = atof (ptr);
										break;
									case 3:
										if (ptr[0] != '=') Ctrl->E.diffuse = atof (ptr);
										break;
									case 4:
										if (ptr[0] != '=') Ctrl->E.specular = atof (ptr);
										break;
									case 5:
										if (ptr[0] != '=') Ctrl->E.shine = atof (ptr);
										break;
									default:
										break;
								}
								entry++;
							}
							break;
					}
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
				case 'L':
					Ctrl->L.active = TRUE;
					strncpy (Ctrl->L.mode, &argv[i][2], (size_t)4);
					/* We turn on geographic coordinates if -Lg is given by faking -fg */
					if (! strcmp (Ctrl->L.mode, "g")) GMT_parse_f_option ("g");
					break;
				case 'M':
					Ctrl->M.active = TRUE;
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					j = 2;
					if (argv[i][j]) {
						if (argv[i][j] == 't' || argv[i][j] == 'T') {
							Ctrl->N.mode = 1;
							j++;
						}
						else if (argv[i][j] == 'e' || argv[i][j] == 'E') {
							Ctrl->N.mode = 2;
							j++;
						}
						j = sscanf(&argv[i][j], "%lf/%lf/%lf", &Ctrl->N.norm, &Ctrl->N.sigma, &Ctrl->N.offset);
					}
					break;

				case 'S':
					Ctrl->S.active = TRUE;
					Ctrl->S.file = strdup (&argv[i][2]);
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
		fprintf (stderr,"grdgradient %s - Compute directional gradients from grid files\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdgradient <infile> -G<outfile> [-A<azim>[/<azim2>]] [-D[a][o][n]]\n");
		fprintf (stderr, "[-E[s|p]<azim>/<elev[ambient/diffuse/specular/shine]>]\n");
		fprintf (stderr, "[-L<flag>] [-M] [-N[t_or_e][<amp>[/<sigma>[/<offset>]]]] [-S<slopefile>] [-V]\n\n");
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
		fprintf (stderr,"\t<infile> is name of input grid file.\n");
		fprintf (stderr,"\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A sets azimuth (0-360 CW from North (+y)) for directional derivatives.\n");
		fprintf (stderr, "\t  -A<azim>/<azim2> will compute two directions and save the one larger in magnitude.\n");
		fprintf (stderr, "\t-D finds the direction of grad z.\n");
		fprintf (stderr, "\t   Append c to get cartesian angle (0-360 CCW from East (+x)) [Default:  azimuth].\n");
		fprintf (stderr, "\t   Append o to get bidirectional orientations [0-180] rather than directions [0-360].\n");
		fprintf (stderr, "\t   Append n to add 90 degrees to the values from c or o.\n");
		fprintf (stderr, "\t-E Compute Lambertian radiance appropriate to use with grdimage/grdview.\n");
		fprintf (stderr, "\t   -E<azim/elev> sets azimuth and elevation of light vector.\n");
		fprintf (stderr, "\t   -E<azim/elev/ambient/diffuse/specular/shine> sets azim, elev and\n");
		fprintf (stderr, "\t    other parameters that control the reflectance properties of the surface.\n");
		fprintf (stderr, "\t    Default values are: 0.55/0.6/0.4/10.\n");
		fprintf (stderr, "\t    Specify '=' to get the default value (e.g., -E60/30/=/0.5).\n");
		fprintf (stderr, "\t   Append s to use a simpler Lambertian algorithm (note that with this form\n");
		fprintf (stderr, "\t   you only have to provide the azimuth and elevation parameters).\n");
		fprintf (stderr, "\t   Append p to use the Peucker piecewise linear approximation (simpler but faster algorithm).\n");
		fprintf (stderr, "\t   Note that in this case the azimuth and elevation are hardwired to 315 and 45 degrees.\n");
		fprintf (stderr, "\t   This means that even if you provide other values they will be ignored.\n");
		fprintf (stderr, "\t-G output file for results from -A or -D.\n");
		fprintf (stderr, "\t-L sets boundary conditions.  <flag> can be either:\n");
		fprintf (stderr, "\t   g for geographic boundary conditions\n");
		fprintf (stderr, "\t   or one or both of\n");
		fprintf (stderr, "\t   x for periodic boundary conditions on x\n");
		fprintf (stderr, "\t   y for periodic boundary conditions on y\n");
		fprintf (stderr, "\t   [Default:  Natural conditions].\n");
		fprintf (stderr, "\t-M to use map units.  In this case, dx,dy of grid\n");
		fprintf (stderr, "\t   will be converted from degrees lon,lat into meters (Flat-earth approximation).\n");
		fprintf (stderr, "\t   Default computes gradient in units of data/grid_distance.\n");
		fprintf (stderr, "\t-N will normalize gradients so that max |grad| = <amp> [1.0].\n");
		fprintf (stderr, "\t  -Nt will make atan transform, then scale to <amp> [1.0].\n");
		fprintf (stderr, "\t  -Ne will make exp  transform, then scale to <amp> [1.0].\n");
		fprintf (stderr, "\t  -Nt<amp>/<sigma>[/<offset>] or -Ne<amp>/<sigma>[/<offset>] sets sigma\n");
		fprintf (stderr, "\t     (and offset) for transform. [sigma, offset estimated from data].\n");
		fprintf (stderr, "\t-S output file for |grad z|; requires -D.\n");
		GMT_explain_option ('V');
		exit (EXIT_FAILURE);
	}

	if (!(Ctrl->A.active || Ctrl->D.active || Ctrl->E.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -A, -D, or -E\n", GMT_program);
		error++;
	}
	if (Ctrl->S.active && !Ctrl->S.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  Must specify output file\n", GMT_program);
		error++;
	}
	if (!Ctrl->G.file && !Ctrl->S.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option:  Must specify output file\n", GMT_program);
		error++;
	}
	if (!infile) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input file\n", GMT_program);
		error++;
	}
	if (Ctrl->N.norm <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Normalization amplitude must be > 0\n", GMT_program);
		error++;
	}
	if (sigma_set && (Ctrl->N.sigma <= 0.0) ) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Sigma must be > 0\n", GMT_program);
		error++;
	}
	if (Ctrl->E.active && Ctrl->E.mode > 1 && (Ctrl->E.elevation < 0.0 || Ctrl->E.elevation > 90.0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Use 0-90 degree range for elevation\n", GMT_program);
		error++;
	}
	if (Ctrl->E.active && (Ctrl->A.active || Ctrl->D.active || Ctrl->S.active)) {
		fprintf (stderr, "%s: WARNING: -E option overrides -A, -D or -S\n", GMT_program);
		Ctrl->A.active = Ctrl->D.active = Ctrl->S.active = FALSE;
	}
	
	if (Ctrl->L.active && GMT_boundcond_parse (&edgeinfo, Ctrl->L.mode)) error++;

	if (error) exit (EXIT_FAILURE);

	GMT_err_fail (GMT_read_grd_info (infile, &header), infile);

	if (Ctrl->N.active && Ctrl->N.sigma != 0.0) sigma_set = TRUE;
	if (Ctrl->N.active && Ctrl->N.offset != 0.0) offset_set = TRUE;
	if (Ctrl->A.active) {
		while (Ctrl->A.azimuth[0] < 0.0) Ctrl->A.azimuth[0] += 360.0;
		while (Ctrl->A.azimuth[0] > 360.0) Ctrl->A.azimuth[0] -= 360.0;
		if (Ctrl->A.two) {
			while (Ctrl->A.azimuth[1] < 0.0) Ctrl->A.azimuth[1] += 360.0;
			while (Ctrl->A.azimuth[1] > 360.0) Ctrl->A.azimuth[1] -= 360.0;
		}
	}
	if (Ctrl->E.active) {
		while (Ctrl->E.azimuth < 0.0) Ctrl->E.azimuth += 360.0;
		while (Ctrl->E.azimuth > 360.0) Ctrl->E.azimuth -= 360.0;
	}
	if (Ctrl->E.mode == 2) {
		p0 = cosd(90.0 - Ctrl->E.azimuth) * tand(90.0 - Ctrl->E.elevation);
		q0 = sind(90.0 - Ctrl->E.azimuth) * tand(90.0 - Ctrl->E.elevation);
		p0q0_cte = sqrt(1 + p0*p0 + q0*q0);
	}
	if (Ctrl->E.mode == 3) {
		Ctrl->E.elevation = 90 - Ctrl->E.elevation;
		s[0] = sind(Ctrl->E.azimuth) * cosd(Ctrl->E.elevation);
		s[1] = cosd(Ctrl->E.azimuth) * cosd(Ctrl->E.elevation);
		s[2] = sind(Ctrl->E.elevation);
		k_ads = Ctrl->E.ambient + Ctrl->E.diffuse + Ctrl->E.specular;
	}
	GMT_boundcond_param_prep (&header, &edgeinfo);

	GMT_grd_init (&header, argc, argv, TRUE);
	nm = GMT_get_nm (header.nx, header.ny);
	mx = header.nx + 4;
	my = header.ny + 4;
	nm2 = GMT_get_nm (mx, my);
	data = (float *) GMT_memory (VNULL, (size_t)nm2, sizeof (float), GMT_program);
	GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 2;
	if (Ctrl->S.active) slp = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);

	GMT_err_fail (GMT_read_grd (infile, &header, data, header.x_min, header.x_max, header.y_min, header.y_max, GMT_pad, FALSE), infile);

	/* set boundary conditions:  */

	GMT_boundcond_set (&header, &edgeinfo, GMT_pad, data);

	if (Ctrl->M.active) {
		dx_grid = project_info.DIST_M_PR_DEG * header.x_inc * cosd ((header.y_max + header.y_min) / 2.0);
		dy_grid = project_info.DIST_M_PR_DEG * header.y_inc;
	}
	else {
		dx_grid = header.x_inc;
		dy_grid = header.y_inc;
	}
	x_factor = -1.0 / (2.0 * dx_grid);
	y_factor = -1.0 / (2.0 * dy_grid);
	if (Ctrl->A.active) {
		if (Ctrl->A.two) {
			Ctrl->A.azimuth[1] *= (M_PI / 180.0);
			x_factor2 = x_factor * sin (Ctrl->A.azimuth[1]);
			y_factor2 = y_factor * cos( Ctrl->A.azimuth[1]);
		}
		Ctrl->A.azimuth[0] *= (M_PI / 180.0);
		x_factor *= sin (Ctrl->A.azimuth[0]);
		y_factor *= cos (Ctrl->A.azimuth[0]);
	}

	p[0] = 1;	p[1] = -1;	p[2] = mx;	p[3] = -mx;

	min_gradient = DBL_MAX;	max_gradient = -DBL_MAX;
	ave_gradient = 0.0;
	if (Ctrl->E.mode == 3) {
		lim_x = header.x_max - header.x_min;
		lim_y = header.y_max - header.y_min;
		lim_z = header.z_max - header.z_min;
		scale = MAX(lim_z, MAX(lim_x, lim_y));
		lim_x /= scale;	lim_y /= scale;		lim_z /= scale;
		dx_grid /= lim_x;	dy_grid /= lim_y;
		x_factor = -dy_grid / (2 * lim_z);	y_factor = -dx_grid / (2 * lim_z);
	}
	for (j = k = 0; j < header.ny; j++) {
		if (Ctrl->M.active) {
			lat = GMT_j_to_y (j, header.y_min, header.y_max, header.y_inc, 0.5 * header.node_offset, header.ny);
			dx_grid = project_info.DIST_M_PR_DEG * header.x_inc * cosd (lat);
			if (dx_grid > 0.0) x_factor = -1.0 / (2.0 * dx_grid);	/* Use previous value at the poles */
			if (Ctrl->A.active) {
				if (Ctrl->A.two) {
					x_factor2 = x_factor * sin(Ctrl->A.azimuth[1]);
				}
				x_factor *= sin(Ctrl->A.azimuth[0]);
			}
		}
		for (i = 0; i < header.nx; i++, k++) {
			ij = (j + 2) * mx + i + 2;
			for (n = 0, bad = FALSE; !bad && n < 4; n++) if (GMT_is_fnan (data[ij+p[n]])) bad = TRUE;
			if (bad) {	/* One of corners = NaN, skip */
				data[k] = GMT_f_NaN;
				if (Ctrl->S.active) slp[k] = GMT_f_NaN;
				continue;
			}

			dzdx = (data[ij+1] - data[ij-1]) * x_factor;
			dzdy = (data[ij-mx] - data[ij+mx]) * y_factor;
			if (Ctrl->A.two) {
				dzdx2 = (data[ij+1] - data[ij-1]) * x_factor2;
				dzdy2 = (data[ij-mx] - data[ij+mx]) * y_factor2;
			}

			/* Write output to unused NW corner */

			if (Ctrl->A.active) {	/* Directional derivatives */
				if (Ctrl->A.two) {
					dzds1 = dzdx + dzdy;
					dzds2 = dzdx2 + dzdy2;
					data[k] = (float)((fabs(dzds1) > fabs(dzds2)) ? dzds1 : dzds2);
				}
				else {
					data[k] = (float)(dzdx + dzdy);
				}
				ave_gradient += data[k];
				min_gradient = MIN (min_gradient, data[k]);
				max_gradient = MAX (max_gradient, data[k]);
			}
			else if (Ctrl->D.active) {
				azim = (Ctrl->D.mode & 1) ? atan2d (-dzdy, -dzdx) : 90.0 - atan2d (-dzdy, -dzdx);
				if (Ctrl->D.mode & 4) azim += 90.0;
				if (azim < 0.0) azim += 360.0;
				if (azim >= 360.0) azim -= 360.0;
				if (Ctrl->D.mode & 2 && azim >= 180) azim -= 180.0;
				data[k] = (float)azim;
				if (Ctrl->S.active) slp[k] = (float)hypot (dzdx, dzdy);
			}
			else {	/* Ctrl->E.active */
				if (Ctrl->E.mode == 3) {
					norm_z = dx_grid * dy_grid;
					mag = d_sqrt(dzdx*dzdx + dzdy*dzdy + norm_z*norm_z);
					dzdx /= mag;	dzdy /= mag;	norm_z /= mag;
					diffuse = MAX(0,(s[0]*dzdx + s[1]*dzdy + s[2]*norm_z)); 
					spec = specular(dzdx, dzdy, norm_z, s);
					spec = pow(spec, Ctrl->E.shine);
					data[k] = (float)((Ctrl->E.ambient+Ctrl->E.diffuse*diffuse+Ctrl->E.specular*spec) / k_ads);
				}
				else if (Ctrl->E.mode == 2)
					data[k] = (float)( (1 + p0*dzdx + q0*dzdy) / (sqrt(1 + dzdx*dzdx + dzdy*dzdy) * p0q0_cte) );
				else	/* Peucker method */
					data[k] = (float)( -0.4285 * (dzdx - dzdy) - 0.0844 * fabs(dzdx  + dzdy) + 0.6599 );
				r_min = MIN (r_min, (double)data[k]);
				r_max = MAX (r_max, (double)data[k]);
			}
			n_used++;
		}
	}

	if (Ctrl->M.active || GMT_io.in_col_type[GMT_Y] == GMT_IS_LAT) {	/* Data is geographic */
		double sum;
		/* If the N or S poles are included then we only want a single estimate at these repeating points */
		if (header.y_min == -90.0 && header.node_offset == 0) {	/* Average all the multiple N pole estimates */
			for (k = 0, sum = 0.0; k < header.nx; k++) sum += data[k];
			sum /= header.nx;	/* Average NP gradient */
			for (k = 0; k < header.nx; k++) data[k] = (float)sum;
		}
		if (header.y_min == -90.0 && header.node_offset == 0) {	/* Average all the multiple S pole estimates */
			for (i = 0, k = header.nx * (header.ny - 1), sum = 0.0; i < header.nx; i++, k++) sum += data[k];
			sum /= header.nx;	/* Average SP gradient */
			for (i = 0, k = header.nx * (header.ny - 1); i < header.nx; i++, k++) data[k] = (float)sum;
		}
	}
	
	if (Ctrl->E.active) {	/* data must be scaled to the [-1,1] interval, but we'll do it into [-.95, .95] to not get too bright */
		scale = 1.0 / (r_max - r_min);
		for (k = 0; k < nm; k++) {
			if (GMT_is_fnan (data[k])) continue;
			data[k] = (float)((-1.0 + 2.0 * ((data[k] - r_min) * scale)) * 0.95);
		}
	}

	if (offset_set)
		ave_gradient = Ctrl->N.offset;
	else
		ave_gradient /= n_used;

	if (Ctrl->A.active) {	/* Report some statistics */

		if (Ctrl->N.active) {
			if (Ctrl->N.mode == 1) {
				if (sigma_set) {
					denom = 1.0 / Ctrl->N.sigma;
				}
				else {
					denom = 0.0;
					for (k = 0; k < nm; k++) if (!GMT_is_fnan (data[k])) denom += pow(data[k] - ave_gradient, 2.0);
					denom = sqrt( (n_used - 1) / denom);
					Ctrl->N.sigma = 1.0 / denom;
				}
				rpi = 2.0 * Ctrl->N.norm / M_PI;
				for (k = 0; k < nm; k++) if (!GMT_is_fnan (data[k])) data[k] = (float)(rpi * atan((data[k] - ave_gradient)*denom));
				header.z_max = rpi * atan((max_gradient - ave_gradient)*denom);
				header.z_min = rpi * atan((min_gradient - ave_gradient)*denom);
			}
			else if (Ctrl->N.mode == 2) {
				if (!sigma_set) {
					Ctrl->N.sigma = 0.0;
					for (k = 0; k < nm; k++) if (!GMT_is_fnan (data[k])) Ctrl->N.sigma += fabs((double)data[k]);
					Ctrl->N.sigma = M_SQRT2 * Ctrl->N.sigma / n_used;
				}
				denom = M_SQRT2 / Ctrl->N.sigma;
				for (k = 0; k < nm; k++) {
					if (GMT_is_fnan (data[k])) continue;
					if (data[k] < ave_gradient) {
						data[k] = (float)(-Ctrl->N.norm * (1.0 - exp((data[k] - ave_gradient)*denom)));
					}
					else {
						data[k] = (float)(Ctrl->N.norm * (1.0 - exp(-(data[k] - ave_gradient)*denom)));
					}
				}
				header.z_max = Ctrl->N.norm * (1.0 - exp(-(max_gradient - ave_gradient)*denom));
				header.z_min = -Ctrl->N.norm * (1.0 - exp((min_gradient - ave_gradient)*denom));
			}
                	else {
				if ( (max_gradient - ave_gradient) > (ave_gradient - min_gradient) ) {
					denom = Ctrl->N.norm / (max_gradient - ave_gradient);
				}
				else {
					denom = Ctrl->N.norm / (ave_gradient - min_gradient);
				}
				for (k = 0; k < nm; k++) if (!GMT_is_fnan (data[k])) data[k] = (float)((data[k] - ave_gradient) * denom);
				header.z_max = (max_gradient - ave_gradient) * denom;
				header.z_min = (min_gradient - ave_gradient) * denom;
			}
		}
	}

	/* Now we write out: */

	if (Ctrl->A.active) {
		if (Ctrl->N.active) {
			strcpy (header.title, "Normalized directional derivative(s)");
		}
		else {
			strcpy (header.title, "Directional derivative(s)");
		}
		sprintf (format, "\t%s\t%s\t%s\t%s\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		if (gmtdefs.verbose) {
			fprintf (stderr, "%s:  Min Mean Max sigma intensities:", GMT_program);
			fprintf (stderr, format, min_gradient, ave_gradient, max_gradient, Ctrl->N.sigma);
		}
	}
	else {
		if (Ctrl->E.mode > 1)
			strcpy (header.title, "Lambertian radiance");
		else if (Ctrl->E.mode == 1)
			strcpy (header.title, "Peucker piecewise linear radiance");
		else
			strcpy (header.title, "Directions of maximum slopes");
	}

	GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 0;	/* Because of the shift */

	if (Ctrl->G.active)
		GMT_err_fail (GMT_write_grd (Ctrl->G.file, &header, data, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);

	GMT_free ((void *) data);

	if (Ctrl->S.active) {
		strcpy (header.title, "Magnitude of maximum slopes");
		GMT_err_fail (GMT_write_grd (Ctrl->S.file, &header, slp, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->S.file);
		GMT_free ((void *)slp);
	}

	Free_grdgradient_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

double specular (double nx, double ny, double nz, double *s) {
	/* SPECULAR Specular reflectance.
	   R = SPECULAR(Nx,Ny,Nz,S,V) returns the reflectance of a surface with
	   normal vector components [Nx,Ny,Nz].  S and V specify the direction
	   to the light source and to the viewer, respectively. 
	   For the time beeing I'm using V = [azim elev] = [0 90] so the following

	   V[0] =  sind(V[0])*cosd(V[1]);
	   V[1] = -cosd(V[0])*cosd(V[1]);
	   V[2] =  sind(V[1]);

	   Reduces to V[0] = 0;		V[1] = 0;	V[2] = 1 */

	/*r = MAX(0,2*(s[0]*nx+s[1]*ny+s[2]*nz).*(v[0]*nx+v[1]*ny+v[2]*nz) - (v'*s)*ones(m,n)); */

	return (MAX(0, 2 * (s[0]*nx + s[1]*ny + s[2]*nz) * nz - s[2]));
}

void *New_grdgradient_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDGRADIENT_CTRL *C;
	
	C = (struct GRDGRADIENT_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDGRADIENT_CTRL), "New_grdgradient_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->E.ambient = 0.55;
	C->E.diffuse = 0.6;
	C->E.specular = 0.4;
	C->E.shine = 10;
	C->N.norm = 1.0;		
	return ((void *)C);
}

void Free_grdgradient_Ctrl (struct GRDGRADIENT_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	if (C->S.file) free ((void *)C->S.file);	
	GMT_free ((void *)C);	
}
