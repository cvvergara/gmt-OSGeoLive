/*--------------------------------------------------------------------
 *	$Id: grdimage.c 10173 2014-01-01 09:52:34Z pwessel $
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
 * grdimage will read a grid file and image the area using the PostScript
 * image command. If non-linear scaling is chosen, we first resample the data
 * onto the new grid before calling the image command.  The output image
 * will be 1-, 8-, or 24-bit depending on colors used.  As an option, grdimage
 * can read three grid files with red, green, blue components in the 0-255
 * range and use those colors directly (no cpt required).
 *
 * Author:	Paul Wessel
 * Date:	17-SEP-2001
 * Ver:		4
 *
 */

#include "gmt.h"
#include "pslib.h"

struct GRDIMAGE_CTRL {
	struct C {	/* -C<cptfile> */
		GMT_LONG active;
		char *file;
	} C;
#ifdef USE_GDAL
	struct D {	/* -D to read gdal file */
		GMT_LONG active;
		GMT_LONG mode;	/* Use info of -R option to reference image */
	} D;
#endif
	struct E {	/* -Ei|<dpi> */
		GMT_LONG active;
		GMT_LONG device_dpi;
		GMT_LONG dpi;
	} E;
	struct G {	/* -G[f|b]<rgb> */
		GMT_LONG active;
		int f_rgb[3];
		int b_rgb[3];
	} G;
	struct I {	/* -I<intensfile> */
		GMT_LONG active;
		char *file;
	} I;
	struct M {	/* -M */
		GMT_LONG active;
	} M;
	struct N {	/* -N */
		GMT_LONG active;
	} N;
	struct Q {	/* -Q */
		GMT_LONG active;
	} Q;
	struct S {	/* -S[-]b|c|l|n[/<threshold>] */
		GMT_LONG active;
		GMT_LONG antialias;
		GMT_LONG interpolant;
		double threshold;
	} S;
};

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE, do_rgb = FALSE, done, need_to_project, normal_x, normal_y, resampled = FALSE;

	unsigned char *bitimage_8 = NULL, *bitimage_24 = NULL, *rgb_used = NULL;

	char *grdfile[3], *c_method[2] = {"colorimage", "colortiles",};

	GMT_LONG i, j, k, byte, nx, ny, nx_proj = 0, ny_proj = 0, index, grid_type = 0, n_grids = 0;
	GMT_LONG PS_colormask_off = 0, PS_colormask, PS_interpolate, nm, nm2 = 0, node, kk;
	int rgb[3], try;
	
	float *tmp1[3], *tmp2 = NULL, *map[3], *intensity = NULL;

	double  dx, dy, x_side, y_side, x0 = 0.0, y0 = 0.0;
	double west, east, south, north, data_west, data_east, data_south, data_north;

	struct GRD_HEADER g_head[3], r_head[3], i_head, j_head;
	struct GMT_EDGEINFO edgeinfo;
	struct GRDIMAGE_CTRL *Ctrl = NULL;

#ifdef USE_GDAL
	GMT_LONG do_indexed = FALSE;
	unsigned char *r_table = NULL, *g_table = NULL, *b_table = NULL;
	struct GDALREAD_CTRL *to_gdalread = NULL;
	struct GD_CTRL *from_gdalread = NULL;
#endif

	void GMT_set_proj_limits (struct GRD_HEADER *r, struct GRD_HEADER *g, GMT_LONG projected);
	void *New_grdimage_Ctrl (), Free_grdimage_Ctrl (struct GRDIMAGE_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDIMAGE_CTRL *)New_grdimage_Ctrl ();	/* Allocate and initialize a new control structure */

	grdfile[0] = grdfile[1] = grdfile[2] = CNULL;
	west = east = south = north = 0.0;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'B':
				case 'J':
				case 'K':
				case 'O':
				case 'P':
				case 'R':
				case 'U':
				case 'V':
				case 'X':
				case 'x':
				case 'Y':
				case 'y':
				case 'c':
				case 'f':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */

				case 'C':
					Ctrl->C.active = TRUE;
					Ctrl->C.file = strdup (&argv[i][2]);
					break;
#ifdef USE_GDAL
				case 'D':
					Ctrl->D.active = TRUE;
					if (argv[i][2] == 'r') Ctrl->D.mode = TRUE;
					break;
#endif
				case 'E':
					Ctrl->E.active = TRUE;
					if (argv[i][2] == 'i')	/* Interpolate image to device resolution */
						Ctrl->E.device_dpi = TRUE;
					else
						Ctrl->E.dpi = atoi (&argv[i][2]);
					break;
				case 'G':	/* 1-bit fore or background color for transparent masks */
					Ctrl->G.active = TRUE;
					switch (argv[i][2]) {
						case 'F':
						case 'f':
							if (GMT_getrgb (&argv[i][3], Ctrl->G.f_rgb)) {
								GMT_rgb_syntax ('G', " ");
								error++;
							}
							else
								Ctrl->G.b_rgb[0] = -1;
							break;
						case 'B':
						case 'b':
							if (GMT_getrgb (&argv[i][3], Ctrl->G.b_rgb)) {
								GMT_rgb_syntax ('G', " ");
								error++;
							}
							else
								Ctrl->G.f_rgb[0] = -1;
							break;
						default:	/* Same as -Gf */
							if (GMT_getrgb (&argv[i][2], Ctrl->G.f_rgb)) {
								GMT_rgb_syntax ('G', " ");
								error++;
							}
							else
								Ctrl->G.b_rgb[0] = -1;
							break;
					}
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					Ctrl->I.file = strdup (&argv[i][2]);
					break;
				case 'M':
					Ctrl->M.active = TRUE;
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					break;
				case '0':	/* Backwards compatibility */
					gmtdefs.color_image = 0;
					break;
				case '1':	/* Backwards compatibility */
					gmtdefs.color_image = 1;
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
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
				case 'T':
					fprintf (stderr, "%s: Warning: The -T option has become obsolete. Use grdview instead.\n\t-T is similar to -Sn.\n\t-Ts is similar to -Sn -Q.\n\tProcessing continues with these options. Use -E if needed.\n", GMT_program);
					Ctrl->S.active = TRUE;
					Ctrl->S.interpolant = BCR_NEARNEIGHBOR;
					k = 2;
					if (argv[i][2] == 's') Ctrl->Q.active = TRUE, k = 3;
					if (argv[i][k] == 'o') {	/* Want tile outline also */
						fprintf (stderr, "\tThe -To option is no longer supported. Use grdview instead.\n");
						error++;
					}
					break;

				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else {
			if (n_grids == 3) {
				fprintf (stderr, "%s: ERROR, give 1 or 3 grid files\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			grdfile[n_grids++] = argv[i];
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr,"grdimage %s - Plot grid files in 2-D\n\n", GMT_VERSION);
#ifdef USE_GDAL
		fprintf (stderr, "usage: grdimage <grd_z|grd_r grd_g grd_b> %s [%s] [-C<cpt_file>] [-D[r]] [-Ei|<dpi>] [-G[f|b]<rgb>]\n", GMT_J_OPT, GMT_B_OPT);
#else
		fprintf (stderr, "usage: grdimage <grd_z|grd_r grd_g grd_b> %s [%s] [-C<cpt_file>] [-Ei|<dpi>] [-G[f|b]<rgb>]\n", GMT_J_OPT, GMT_B_OPT);
#endif
		fprintf (stderr, "\t[-I<intensity_file>] [-K] [-M] [-N] [-O] [-P] [-Q] [%s] [-S[-]b|c|l|n[/<threshold>]] [-T]\n", GMT_Rgeo_OPT);
		fprintf (stderr, "\t[%s] [-V] [%s] [%s] [%s]\n\n", GMT_U_OPT, GMT_X_OPT, GMT_Y_OPT, GMT_c_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<grd_z> is data set to be plotted.  Its z-values are in user units and will be\n");
		fprintf (stderr, "\t  converted to rgb colors via the cpt file.  Alternatively, give three separate\n");
		fprintf (stderr, "\t  grid files that contain the red, green, and blue components in the 0-255 range.\n");
#ifdef USE_GDAL
		fprintf (stderr, "\t  If -D is used then <grd_z> is instead expected to be an image.\n");
#endif
		GMT_explain_option ('j');
		fprintf (stderr,"\n\tOPTIONS:\n");
		GMT_explain_option ('b');
		fprintf (stderr, "\t-C color palette file to convert z to rgb.\n");
#ifdef USE_GDAL
		fprintf (stderr, "\t-D is used to read an image via GDAL.  Append r to equate image region to -R region.\n");
#endif
		fprintf (stderr, "\t-E sets dpi for the projected grid which must be constructed\n");
		fprintf (stderr, "\t   if -Jx or -Jm is not selected [Default gives same size as input grid].\n");
		fprintf (stderr, "\t   Give i to do the interpolation in PostScript at device resolution.\n");
		GMT_rgb_syntax ('G', "sets transparency color for images that otherwise would result in 1-bit images\n\t  ");
		fprintf (stderr, "\t-I use illumination.  Append name of intensity grid file.\n");
		GMT_explain_option ('K');
		fprintf (stderr, "\t-M force monochrome image.\n");
		fprintf (stderr, "\t-N Do not clip image at the map boundary.\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "\t-Q use PS Level 3 colormasking to make nodes with z = NaN transparent.\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t-S Determines the interpolation mode (b = B-spline, c = bicubic, l = bilinear,\n");
		fprintf (stderr, "\t   n = nearest-neighbor) [Default: bicubic].\n");
		fprintf (stderr, "\t   Optionally, prepend - to switch off antialiasing [Default: on].\n");
		fprintf (stderr, "\t   Append /<threshold> to change the minimum weight in vicinity of NaNs. A threshold of\n");
		fprintf (stderr, "\t   1.0 requires all nodes involved in interpolation to be non-NaN; 0.5 will interpolate\n");
		fprintf (stderr, "\t   about half way from a non-NaN to a NaN node [Default: 0.5].\n");
		fprintf (stderr, "\t-T OBSOLETE: See man pages.\n");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_explain_option ('X');
		GMT_explain_option ('c');
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	do_rgb = (n_grids == 3);

	for (i = 0; i < n_grids; i++) {
		if (!grdfile[i]) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input file\n", GMT_program);
			error++;
		}
	}
#ifdef USE_GDAL
	if (Ctrl->D.active) {} else
#endif
	if (!Ctrl->C.file && !do_rgb) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify color palette table\n", GMT_program);
		error++;
	}
	if (Ctrl->I.active && !Ctrl->I.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I option:  Must specify intensity file\n", GMT_program);
		error++;
	}
	if (Ctrl->E.active && !Ctrl->E.device_dpi && Ctrl->E.dpi <= 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  dpi must be positive\n", GMT_program);
		error++;
	}
	if (Ctrl->M.active && Ctrl->Q.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Q option:  Cannot use -M when doing colormasking\n", GMT_program);
		error++;
	}
	if (Ctrl->G.f_rgb[0] < 0 && Ctrl->G.b_rgb[0] < 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option:  Only one of fore/back-ground can be transparent for 1-bit images\n", GMT_program);
		error++;
	}
	if (Ctrl->S.active && (Ctrl->S.threshold < 0.0 || Ctrl->S.threshold > 1.0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  threshold must be in [0,1] range\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

#ifdef USE_GDAL
	if (Ctrl->D.active) {
		int pad;
		/* One more test though */
		if (Ctrl->I.active) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR cannot use -D and -I options.\n", GMT_program);
			exit (EXIT_FAILURE);
		}
		if (Ctrl->D.mode && !project_info.region_supplied) {
			fprintf (stderr, "%s: GMT WARNING: -Dr without -R makes no sense. Ignoring -Dr.\n", GMT_program);
			Ctrl->D.mode = FALSE;
		}

		/* Allocate new control structures */
		to_gdalread = (struct GDALREAD_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GDALREAD_CTRL), "New_gdalread_Ctrl");
		from_gdalread = (struct GD_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GD_CTRL), "New_gd_Ctrl");
		to_gdalread->F.active = 1;

		if (project_info.region_supplied && !Ctrl->D.mode) {
			char	strR [128]; 
			sprintf(strR, "-R%.10f/%.10f/%.10f/%.10f", west, east, south, north);
			to_gdalread->R.region = strR;
		}

		pad = (project_info.projection != GMT_LINEAR || project_info.xyz_projection[0] == GMT_LOG10 || 
			project_info.xyz_projection[0] == GMT_POW || project_info.xyz_projection[1] == GMT_LOG10 || 
			project_info.xyz_projection[1] == GMT_POW || Ctrl->E.dpi > 0);
		if (pad) {
			to_gdalread->p.active = 1;
			to_gdalread->p.pad = 2;
		}

		if (GMT_gdalread ( grdfile[0], to_gdalread, from_gdalread)) {
			fprintf (stderr, "%s: ERROR reading file with gdalread.\n", GMT_program);
			exit (EXIT_FAILURE);
		}

		if (!from_gdalread->UInt8.active) {
			fprintf(stderr, "%s: Using data type other than byte (unsigned char) is not implemented\n", GMT_program);
			exit (EXIT_FAILURE);
		}

		do_rgb = (from_gdalread->RasterCount == 3);
		if (do_rgb) n_grids = 3;	/* To be compatible with original algo */

		if (gmtdefs.verbose && from_gdalread->ProjectionRefPROJ4 != CNULL)
			fprintf (stderr, "Data projection (Proj4 type)\n\t%s\n", from_gdalread->ProjectionRefPROJ4);

		GMT_grd_init (&g_head[i], argc, argv, FALSE);
		GMT_grd_init (&r_head[i], argc, argv, FALSE);
		if (!Ctrl->D.mode) {
			g_head[0].x_min = g_head[1].x_min = g_head[2].x_min = from_gdalread->hdr[0];
			g_head[0].x_max = g_head[1].x_max = g_head[2].x_max = from_gdalread->hdr[1];
			g_head[0].y_min = g_head[1].y_min = g_head[2].y_min = from_gdalread->hdr[2];
			g_head[0].y_max = g_head[1].y_max = g_head[2].y_max = from_gdalread->hdr[3];
			g_head[0].x_inc = g_head[1].x_inc = g_head[2].x_inc = from_gdalread->hdr[7];
			g_head[0].y_inc = g_head[1].y_inc = g_head[2].y_inc = from_gdalread->hdr[8];
		}
		else {
			g_head[0].x_min = g_head[1].x_min = g_head[2].x_min = west;
			g_head[0].x_max = g_head[1].x_max = g_head[2].x_max = east;
			g_head[0].y_min = g_head[1].y_min = g_head[2].y_min = south;
			g_head[0].y_max = g_head[1].y_max = g_head[2].y_max = north;
			g_head[0].x_inc = g_head[1].x_inc = g_head[2].x_inc = (g_head[0].x_max - g_head[0].x_min) / from_gdalread->RasterXsize;
			g_head[0].y_inc = g_head[1].y_inc = g_head[2].y_inc = (g_head[0].y_max - g_head[0].y_min) / from_gdalread->RasterYsize;
		}
		g_head[0].nx = g_head[1].nx = g_head[2].nx = from_gdalread->RasterXsize;
		g_head[0].ny = g_head[1].ny = g_head[2].ny = from_gdalread->RasterYsize;
		g_head[0].node_offset = g_head[1].node_offset = g_head[2].node_offset = (int)from_gdalread->hdr[6];
	}
#endif

	/* Get/calculate a color palette file */

	if (!do_rgb) {
		if (Ctrl->C.active)
			GMT_read_cpt (Ctrl->C.file);
#ifdef USE_GDAL
		else if (Ctrl->D.active) {
			if (from_gdalread->ColorMap == NULL && !strncmp(from_gdalread->ColorInterp,"Gray",4)) {
    				r_table = (unsigned char *) malloc(256 * sizeof(*r_table));
				for (i = 0; i < 256; i++) r_table[i] = (unsigned char)i;
				GMT_gray = TRUE;
			}
			else if (from_gdalread->ColorMap != NULL) {
    				r_table = (unsigned char *) malloc(256 * sizeof(*r_table));
    				g_table = (unsigned char *) malloc(256 * sizeof(*g_table));
    				b_table = (unsigned char *) malloc(256 * sizeof(*b_table));
				for (i = 0; i < 256; i++) {
					r_table[i] = from_gdalread->ColorMap[i*4];	/* 4 because color table is RGBA */
					g_table[i] = from_gdalread->ColorMap[i*4 + 1];
					b_table[i] = from_gdalread->ColorMap[i*4 + 2];
				}
				do_indexed = TRUE;		/* Now it will be RGB */
				GMT_gray = FALSE;
			}
		}
#endif
	}

	PS_interpolate = (Ctrl->E.device_dpi) ? -1 : +1;
	PS_colormask   = (Ctrl->Q.active) ? -1 : +1;;

	if (gmtdefs.verbose) fprintf (stderr, "%s: Allocates memory and read data file\n", GMT_program);

#ifdef USE_GDAL
	if (!Ctrl->D.active) {
#endif
		for (i = 0; i < n_grids; i++) {
			GMT_grd_init (&g_head[i], argc, argv, FALSE);
			GMT_grd_init (&r_head[i], argc, argv, FALSE);
			GMT_err_fail (GMT_read_grd_info (grdfile[i], &g_head[i]), grdfile[i]);
		}
#ifdef USE_GDAL
	}
#endif
	if (do_rgb) {	/* Must ensure all three grids are coregistered */
		if (!(g_head[0].x_min == g_head[1].x_min && g_head[0].x_min == g_head[2].x_min)) error++;
		if (!(g_head[0].y_min == g_head[1].y_min && g_head[0].y_min == g_head[2].y_min)) error++;
		if (!(g_head[0].x_inc == g_head[1].x_inc && g_head[0].x_inc == g_head[2].x_inc)) error++;
		if (!(g_head[0].nx == g_head[1].nx && g_head[0].nx == g_head[2].nx)) error++;
		if (!(g_head[0].ny == g_head[1].ny && g_head[0].ny == g_head[2].ny)) error++;
		if (!(g_head[0].node_offset == g_head[1].node_offset && g_head[0].node_offset == g_head[2].node_offset)) error++;
		if (error) {
			fprintf (stderr, "%s: The r, g, and b grids are not congruent\n", GMT_program);
			exit (EXIT_FAILURE);
		}
	}

	/* Determine what wesn to pass to map_setup */

	if (!project_info.region_supplied) {
		west = g_head[0].x_min;
		east = g_head[0].x_max;
		south = g_head[0].y_min;
		north = g_head[0].y_max;
	}

	GMT_err_fail (GMT_map_setup (west, east, south, north), "");
	
	/* Determine if grid is to be projected */

	need_to_project = (project_info.projection != GMT_LINEAR || project_info.xyz_projection[0] == GMT_LOG10 || project_info.xyz_projection[0] == GMT_POW || project_info.xyz_projection[1] == GMT_LOG10 || project_info.xyz_projection[1] == GMT_POW || Ctrl->E.dpi > 0);

	/* Determine the wesn to be used to read the grid file; or exit if file is outside -R */

	if (GMT_grd_setregion (&g_head[0], &data_west, &data_east, &data_south, &data_north, need_to_project * Ctrl->S.interpolant)) {
		/* No grid to plot; just do empty map and exit */
		GMT_plotinit (argc, argv);
		GMT_map_basemap ();
		GMT_plotend ();
		Free_grdimage_Ctrl (Ctrl);	/* Deallocate control structure */
		GMT_end (argc, argv);
		exit (EXIT_SUCCESS);
	}

	nx = GMT_get_n (data_west, data_east, g_head[0].x_inc, g_head[0].node_offset);
	ny = GMT_get_n (data_south, data_north, g_head[0].y_inc, g_head[0].node_offset);

#ifdef USE_GDAL
	if (Ctrl->D.active) {	/* Trust more on info from gdal to make it more stable against pixel vs grid registration troubles */
		nx = from_gdalread->RasterXsize;
		ny = from_gdalread->RasterYsize;
	}
#endif

#if 0
	fprintf (stderr, "%g %g %g %g / %g %g\n", g_head[0].x_min, g_head[0].x_max, g_head[0].y_min, g_head[0].y_max, g_head[0].x_inc, g_head[0].y_inc);
	fprintf (stderr, "%g %g %g %g\n", data_west, data_east, data_south, data_north);
	fprintf (stderr, "%li %li\n", nx, ny);
#endif

	GMT_boundcond_init (&edgeinfo);

	if (need_to_project) {
		nm = GMT_get_nm (4 + nx, 4 + ny);
		GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 2;
	}
	else
		nm = GMT_get_nm (nx, ny);

	/* Read data */

#ifdef USE_GDAL
	if (Ctrl->D.active) {
		for (i = 0; i < n_grids; i++) {
			tmp1[i] = (float *) GMT_memory (VNULL, nm, sizeof (float), GMT_program);
			for (j = 0; j < nm; j++) tmp1[i][j] = (float)from_gdalread->UInt8.data[j+i*nm];
		}
	}
	else {
#endif
		for (i = 0; i < n_grids; i++) {
			tmp1[i] = (float *) GMT_memory (VNULL, nm, sizeof (float), GMT_program);
			GMT_err_fail (GMT_read_grd (grdfile[i], &g_head[i], tmp1[i], data_west, data_east, data_south, data_north, GMT_pad, FALSE), grdfile[i]);
		}
#ifdef USE_GDAL
	}
#endif

	/* If given, get intensity file or compute intensities */

	if (Ctrl->I.active) {	/* Illumination wanted */

		if (gmtdefs.verbose) fprintf (stderr, "%s: Allocates memory and read intensity file\n", GMT_program);

		GMT_grd_init (&i_head, argc, argv, FALSE);
		GMT_grd_init (&j_head, argc, argv, FALSE);

		GMT_err_fail (GMT_read_grd_info (Ctrl->I.file, &i_head), Ctrl->I.file);

		tmp2 = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);

		GMT_err_fail (GMT_read_grd (Ctrl->I.file, &i_head, tmp2, data_west, data_east, data_south, data_north, GMT_pad, FALSE), Ctrl->I.file);
		if (i_head.nx != g_head[0].nx || i_head.ny != g_head[0].ny) {
			fprintf (stderr, "%s: Intensity file has improper dimensions!\n", GMT_program);
			exit (EXIT_FAILURE);
		}
	}

	if (need_to_project) {	/* Need to resample the grd file */

		if (gmtdefs.verbose) fprintf (stderr, "%s: project grid files\n", GMT_program);

		if (Ctrl->E.dpi == 0) {	/* Use input # of nodes as # of projected nodes */
			nx_proj = g_head[0].nx;
			ny_proj = g_head[0].ny;
		}
		for (i = 0; i < n_grids; i++) {
			GMT_set_proj_limits (&r_head[i], &g_head[i], need_to_project);
			grid_type = (Ctrl->E.dpi > 0) ? 1 : g_head[i].node_offset;	/* Force pixel if dpi is set */
			GMT_err_fail (GMT_grdproject_init (&r_head[i], 0.0, 0.0, nx_proj, ny_proj, Ctrl->E.dpi, grid_type), grdfile[i]);
			nm2 = ((size_t)r_head[i].nx) * ((size_t)r_head[i].ny);
			map[i] = (float *) GMT_memory (VNULL, nm2, sizeof (float), "grdimage");
			GMT_grd_project (tmp1[i], &g_head[i], map[i], &r_head[i], &edgeinfo, Ctrl->S.antialias, Ctrl->S.interpolant, Ctrl->S.threshold, FALSE);
			GMT_free ((void *)tmp1[i]);
		}
		if (Ctrl->I.active) {
			j_head.x_min = r_head[0].x_min;	j_head.x_max = r_head[0].x_max;
			j_head.y_min = r_head[0].y_min;	j_head.y_max = r_head[0].y_max;

			if (Ctrl->E.dpi == 0) {	/* Use input # of nodes as # of projected nodes */
				nx_proj = i_head.nx;
				ny_proj = i_head.ny;
			}
			GMT_err_fail (GMT_grdproject_init (&j_head, 0.0, 0.0, nx_proj, ny_proj, Ctrl->E.dpi, grid_type), Ctrl->I.file);
			intensity = (float *) GMT_memory (VNULL, nm2, sizeof (float), "grdimage");
			GMT_grd_project (tmp2, &i_head, intensity, &j_head, &edgeinfo, Ctrl->S.antialias, Ctrl->S.interpolant, Ctrl->S.threshold, FALSE);
			GMT_free ((void *)tmp2);
		}
		nm = nm2;
		resampled = TRUE;
	}
	else {	/* Simply copy g_head[0] info to r_head[0] */
		struct GRD_HEADER tmp_header;
		for (i = 0; i < n_grids; i++) {
			map[i] = tmp1[i];
			//r_head[i].nx = g_head[i].nx;		r_head[i].ny = g_head[i].ny;
			//r_head[i].x_inc = g_head[i].x_inc;	r_head[i].y_inc = g_head[i].y_inc;
			memcpy ((void *)&tmp_header, (void *)&g_head[i], sizeof (struct GRD_HEADER));
			r_head[i] = g_head[i];
			GMT_set_proj_limits (&r_head[i], &tmp_header, need_to_project);
		}
		if (Ctrl->I.active) {
			//j_head.nx = i_head.nx;		j_head.ny = i_head.ny;
			//j_head.x_inc = i_head.x_inc;	j_head.y_inc = i_head.y_inc;
			intensity = tmp2;
			memcpy ((void *)&tmp_header, (void *)&i_head, sizeof (struct GRD_HEADER));
			j_head = i_head;
			GMT_set_proj_limits (&j_head, &tmp_header, need_to_project);
		}
		grid_type = g_head[0].node_offset;
	}

	GMT_plotinit (argc, argv);

	if (!Ctrl->N.active) GMT_map_clip_on (GMT_no_rgb, 3);

	if (gmtdefs.verbose) {
		if (GMT_cpt_pattern) fprintf (stderr, "%s: Warning: Patterns in cpt file only apply to -T\n", GMT_program);
		fprintf (stderr, "%s: Evaluate pixel colors\n", GMT_program);
	}

	if (Ctrl->Q.active) {
		if (GMT_gray) {
			if (gmtdefs.verbose) fprintf (stderr, "%s: Your image is grayscale only but -Q requires 24-bit; image will be converted to 24-bit.\n", GMT_program);
			GMT_gray = FALSE;
			GMT_bfn[GMT_NAN].rgb[0] = 255;	GMT_bfn[GMT_NAN].rgb[1] = GMT_bfn[GMT_NAN].rgb[2] = 0;	/* Arbitrarily pick red as the NaN color since image is gray only */
		}
		rgb_used = (unsigned char *) GMT_memory (VNULL, (size_t)(256*256*256), sizeof (unsigned char), GMT_program);
	}
	if (Ctrl->M.active || GMT_gray)
		bitimage_8 = (unsigned char *) GMT_memory (VNULL, (size_t)nm, sizeof (char), GMT_program);
	else {
		if (PS_colormask == -1) PS_colormask_off = 3;
		bitimage_24 = (unsigned char *) GMT_memory (VNULL, (size_t)(3 * nm + PS_colormask_off), sizeof (char), GMT_program);
		if (PS_colormask == -1) {
			bitimage_24[0] = (unsigned char)GMT_bfn[GMT_NAN].rgb[0];
			bitimage_24[1] = (unsigned char)GMT_bfn[GMT_NAN].rgb[1];
			bitimage_24[2] = (unsigned char)GMT_bfn[GMT_NAN].rgb[2];
		}
	}
	normal_x = !(project_info.projection == GMT_LINEAR && !project_info.xyz_pos[0] && !resampled);
	normal_y = !(project_info.projection == GMT_LINEAR && !project_info.xyz_pos[1] && !resampled);
	
	for (try = 0, done = FALSE; !done && try < 2; try++) {	/* Evaluate colors at least once, or twice if -Q and we need to select another NaN color */
		for (j = 0, byte = PS_colormask_off; j < r_head[0].ny; j++) {
			kk = ((size_t)r_head[0].nx) * ((size_t)(normal_y ? j : r_head[0].ny - j - 1));
			for (i = 0; i < r_head[0].nx; i++) {	/* Compute rgb for each pixel */
				node = kk + (normal_x ? i : r_head[0].nx - i - 1);
#ifdef USE_GDAL
				if (Ctrl->D.active && !do_rgb) {
					index = -1;	/* Ensures no illumination done later */
					rgb[0] = r_table[(int)map[0][node]];
					if (do_indexed) {
						rgb[1] = g_table[(int)map[0][node]];
						rgb[2] = b_table[(int)map[0][node]];
					}
				}
				else
#endif
				if (do_rgb) {
					for (k = 0; k < 3; k++) {
						if (GMT_is_fnan (map[k][node])) {	/* If one is NaN they are all assumed to be NaN */
							k = 3;
							memcpy ((void *)rgb, (void *)GMT_bfn[GMT_NAN].rgb, 3 * sizeof (int));
							index = -1;	/* Ensures no illumination done later */
						}
						else {				/* Set color, let index = 0 so illuminate test will work */
							rgb[k] = irint (map[k][node]);	if (rgb[k] < 0) rgb[k] = 0; else if (rgb[k] > 255) rgb[k] = 255;	/* Clip */
							index = 0;
						}
					}
				}
				else
					index = GMT_get_rgb_from_z (map[0][node], rgb);

				if (Ctrl->I.active && index != -1) GMT_illuminate (intensity[node], rgb);

				if (GMT_gray)	/* Color table only has grays, pick r */
					bitimage_8[byte++] = (unsigned char) rgb[0];
				else if (Ctrl->M.active)	/* Convert rgb to gray using the GMT_YIQ transformation */
					bitimage_8[byte++] = (unsigned char) GMT_YIQ (rgb);
				else {
					bitimage_24[byte++] = (unsigned char) rgb[0];
					bitimage_24[byte++] = (unsigned char) rgb[1];
					bitimage_24[byte++] = (unsigned char) rgb[2];
					if (Ctrl->Q.active && index != -1) rgb_used[(rgb[0]*256 + rgb[1])*256+rgb[2]] = TRUE;	/* Keep track of all r/g/b combinations used except for NaN */
				}
			}
		}

		if (Ctrl->Q.active) {	/* Check that colormasking will work OK */
			index = (GMT_bfn[GMT_NAN].rgb[0]*256 + GMT_bfn[GMT_NAN].rgb[1])*256+GMT_bfn[GMT_NAN].rgb[2];
			if (rgb_used[index]) {	/* This r/g/b also appears in the image as a non-NaN color; find a replacement color */
				for (index = 0, k = -1; k == -1 && index < 256*256*256; index++) if (!rgb_used[index]) k = index;
				if (k == -1) {
					fprintf (stderr, "%s: Warning: Colormasking will fail as there is no unused color that can represent transparency\n", GMT_program);
					done = TRUE;
				}
				else {	/* Pick the first unused color and let it play the role of the NaN color for transparency */
					bitimage_24[0] = (unsigned char)(k >> 16);
					bitimage_24[1] = (unsigned char)((k >> 8) & 255);
					bitimage_24[2] = (unsigned char)(k & 255);
					if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: transparency color reset from %d/%d/%d to color %d/%d/%d\n", GMT_program, \
						GMT_bfn[GMT_NAN].rgb[0], GMT_bfn[GMT_NAN].rgb[1], GMT_bfn[GMT_NAN].rgb[2], (int)bitimage_24[0], (int)bitimage_24[1], (int)bitimage_24[2]);
					GMT_bfn[GMT_NAN].rgb[0] = (int)bitimage_24[0];
					GMT_bfn[GMT_NAN].rgb[1] = (int)bitimage_24[1];
					GMT_bfn[GMT_NAN].rgb[2] = (int)bitimage_24[2];
				}	
			}
		}
		else
			done = TRUE;
	}
	if (Ctrl->Q.active) GMT_free ((void *)rgb_used);
	
	for (i = 0; i < n_grids; i++) GMT_free ((void *)map[i]);
	if (Ctrl->I.active) GMT_free ((void *)intensity);

	/* Get actual size of each pixel */

	dx = GMT_get_inc (r_head[0].x_min, r_head[0].x_max, r_head[0].nx, r_head[0].node_offset);
	dy = GMT_get_inc (r_head[0].y_min, r_head[0].y_max, r_head[0].ny, r_head[0].node_offset);

	/* Set lower left position of image on map */

	x0 = r_head[0].x_min;	y0 = r_head[0].y_min;
	if (grid_type == 0) {	/* Grid registration, move 1/2 pixel down/left */
		x0 -= 0.5 * dx;
		y0 -= 0.5 * dy;
	}

	x_side = dx * r_head[0].nx;
	y_side = dy * r_head[0].ny;

	if (gmtdefs.verbose) fprintf (stderr, "%s: Creating PostScript image ", GMT_program);

	if (GMT_gray) for (kk = 0, GMT_b_and_w = TRUE; GMT_b_and_w && kk < nm; kk++) if (!(bitimage_8[kk] == 0 || bitimage_8[kk] == 255)) GMT_b_and_w = FALSE;

	if (GMT_b_and_w) {	/* Can get away with 1 bit image */
		GMT_LONG nx8, shift, byte, b_or_w, nx_pixels, k8;
		unsigned char *bit;

		if (gmtdefs.verbose) fprintf (stderr, "[1-bit B/W image]\n");

		nx8 = (GMT_LONG)ceil (r_head[0].nx / 8.0);
		nx_pixels = nx8 * 8;
		bit = (unsigned char *) GMT_memory (VNULL, (size_t)(nx8 * r_head[0].ny), sizeof (char), GMT_program);

		for (j = k = k8 = 0; j < r_head[0].ny; j++) {
			shift = byte = 0;
			for (i = 0; i < r_head[0].nx; i++, k++) {
				b_or_w = (bitimage_8[k] == 255);
				byte |= b_or_w;
				shift++;
				if (shift == 8) {	/* Time to dump out byte */
					bit[k8++] = (unsigned char) byte;
					byte = shift = 0;
				}
				else
					byte <<= 1;
			}
			if (shift) {
				byte |= 1;
				shift++;
				while (shift < 8) {
					byte <<= 1;
					byte |= 1;
					shift++;
				}
				bit[k8++] = (unsigned char) byte;
			}
		}
		GMT_free ((void *)bitimage_8);

		x_side = nx_pixels * dx;
		ps_bitimage (x0, y0, x_side, y_side, bit, nx_pixels, r_head[0].ny, FALSE, Ctrl->G.f_rgb, Ctrl->G.b_rgb);
		GMT_free ((void *)bit);
	}
	else if (GMT_gray || Ctrl->M.active) {
		if (gmtdefs.verbose) fprintf (stderr, "[8-bit grayshade image]\n");
		ps_colorimage (x0, y0, x_side, y_side, bitimage_8, r_head[0].nx, r_head[0].ny, 8 *PS_interpolate);
		GMT_free ((void *)bitimage_8);
	}
	else {
		if (gmtdefs.verbose) fprintf (stderr, "24-bit [%s]\n", c_method[gmtdefs.color_image]);
		GMT_color_image (x0, y0, x_side, y_side, bitimage_24, r_head[0].nx * PS_colormask, r_head[0].ny, 24 * PS_interpolate);
		GMT_free ((void *)bitimage_24);
	}

	if (!Ctrl->N.active) GMT_map_clip_off();

	GMT_map_basemap ();

	GMT_plotend ();

	Free_grdimage_Ctrl (Ctrl);	/* Deallocate control structure */

#ifdef USE_GDAL
	if (Ctrl->D.active) {
		if (r_table) free((void *) r_table);
		if (g_table) {
			free((void *) g_table);
			free((void *) b_table);
		}
		GMT_free ((void *)to_gdalread);
		GMT_free ((void *)from_gdalread->UInt8.data);
		if (from_gdalread->ColorMap == NULL) GMT_free((void *) from_gdalread->ColorMap);
		GMT_free((void *) from_gdalread);
	}
#endif

	GMT_end (argc, argv);
	exit (EXIT_SUCCESS);

	return(0);
}

void GMT_set_proj_limits (struct GRD_HEADER *r, struct GRD_HEADER *g, GMT_LONG projected)
{
	/* Sets the projected extent of the grid given the map projection
	 * The extreme x/y coordinates are returned in r, and dx/dy, and
	 * nx/ny are set accordingly.  Not that some of these may change
	 * if GMT_grdproject_init is called at a later stage */

	GMT_LONG i, j;
	GMT_LONG all_lats = FALSE, all_lons = FALSE;
	double lon, lat, x, y;

	r->nx = g->nx;
	r->ny = g->ny;
	r->node_offset = g->node_offset;

	if (project_info.projection == GMT_GENPER && project_info.g_width != 0.0) {
		r->x_min = project_info.xmin;	r->x_max = project_info.xmax;
		r->y_min = project_info.ymin;	r->y_max = project_info.ymax;
		return;
	}

	if (GMT_IS_MAPPING) {
		all_lats = GMT_180_RANGE (g->y_max, g->y_min);
		all_lons = GMT_grd_is_global (g);
		if (all_lons && all_lats) {	/* Whole globe, get rectangular box */
			r->x_min = project_info.xmin;	r->x_max = project_info.xmax;
			r->y_min = project_info.ymin;	r->y_max = project_info.ymax;
			return;
		}
	}
	
	/* Must search for extent along perimeter */

	r->x_min = r->y_min = +DBL_MAX;
	r->x_max = r->y_max = -DBL_MAX;

	for (i = j = 0; i < g->nx; i++, j++) {	/* South and north */
		lon = g->x_min + i * g->x_inc;
		GMT_geo_to_xy (lon, g->y_min, &x, &y);
		r->x_min = MIN (r->x_min, x);	r->x_max = MAX (r->x_max, x);
		r->y_min = MIN (r->y_min, y);	r->y_max = MAX (r->y_max, y);
		GMT_geo_to_xy (lon, g->y_max, &x, &y);
		r->x_min = MIN (r->x_min, x);	r->x_max = MAX (r->x_max, x);
		r->y_min = MIN (r->y_min, y);	r->y_max = MAX (r->y_max, y);
	}
	for (i = 0; i < g->ny; j++, i++) {	/* East and west */
		lat = g->y_min + i * g->y_inc;
		GMT_geo_to_xy (g->x_min, lat, &x, &y);
		r->x_min = MIN (r->x_min, x);	r->x_max = MAX (r->x_max, x);
		r->y_min = MIN (r->y_min, y);	r->y_max = MAX (r->y_max, y);
		GMT_geo_to_xy (g->x_max, lat, &x, &y);
		r->x_min = MIN (r->x_min, x);	r->x_max = MAX (r->x_max, x);
		r->y_min = MIN (r->y_min, y);	r->y_max = MAX (r->y_max, y);
	}
	if (projected) {
		if (all_lons) {	/* Full 360, use min/max for x */
			r->x_min = project_info.xmin;	r->x_max = project_info.xmax;
		}
		if (all_lats) {	/* Full -90/+90, use min/max for y */
			r->y_min = project_info.ymin;	r->y_max = project_info.ymax;
		}
	}
}

void *New_grdimage_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDIMAGE_CTRL *C;

	C = (struct GRDIMAGE_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDIMAGE_CTRL), "New_grdimage_Ctrl");

	/* Initialize values whose defaults are not 0/FALSE/NULL */

	C->G.f_rgb[0] = C->G.f_rgb[1] = C->G.f_rgb[2] = 0;
	C->G.b_rgb[0] = C->G.b_rgb[1] = C->G.b_rgb[2] = 255;
	C->S.antialias = TRUE; C->S.interpolant = BCR_BICUBIC; C->S.threshold = 0.5;

	return ((void *)C);
}

void Free_grdimage_Ctrl (struct GRDIMAGE_CTRL *C) {	/* Deallocate control structure */
	if (C->C.file) free ((void *)C->C.file);
	if (C->I.file) free ((void *)C->I.file);
	GMT_free ((void *)C);
}
