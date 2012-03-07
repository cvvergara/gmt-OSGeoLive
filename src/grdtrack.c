/*--------------------------------------------------------------------
 *	$Id: grdtrack.c,v 1.88 2011/07/08 21:27:06 guru Exp $
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
 * grdtrack reads a xyfile, opens the 2d binary gridded grid file, 
 * and samples the dataset at the xy positions with a bilinear or bicubic
 * interpolant.  This new data is added to the input as an extra column
 * and printed to standard output.  In order to evaluate derivatives along
 * the edges of the grid file region, we assume natural bicubic spline
 * boundary conditions (d2z/dn2 = 0, n being the normal to the edge;
 * d2z/dxdy = 0 in the corners).  Rectangles of size x_inc by y_inc are 
 * mapped to [0,1] x [0,1] by affine transformation, and the interpolation
 * done on the normalized rectangle.
 *
 * Author:	Walter H F Smith
 * Date:	23-SEP-1993
 * 
 * Based on the original grdtrack, which had this authorship/date/history:
 *
 * Author:	Paul Wessel
 * Date:	29-JUN-1988
 * Revised:	5-JAN-1990	PW: Updated to v.2.0
 *		4-AUG-1993	PW: Added -Q
 *		14-AUG-1998	PW: GMT 3.1
 *  Modified:	10 Jul 2000 3.3.5  by PW to allow plain -L to indicate geographic coordinates
 *		24 Feb 2006 4.1.1  by PW to allow use of Sandwell/Smith IMG grids directly.
 * Version:	4
 */

#include "gmt.h"

struct GRDTRACK_CTRL {
	struct G {	/* -G<grdfile> */
		GMT_LONG active;
		char *file;
		double scale, lat;
		GMT_LONG mode;
		GMT_LONG type;	/*HIDDEN */
	} G;
	struct L {	/* -L<flag> */
		GMT_LONG active;
		char mode[4];
	} L;
	struct Q {	/* -Q[b|c|l|n][[/]<threshold>] */
		GMT_LONG active;
		GMT_LONG interpolant;
		double threshold;
	} Q;
	struct S {	/* -S */
		GMT_LONG active;
	} S;
	struct Z {	/* -Z */
		GMT_LONG active;
	} Z;
};

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE, pure_ascii = FALSE;

	char line[BUFSIZ], *not_used = NULL;

	GMT_LONG i, j, ix, iy, mx, my, nx, ny, n_fields;
	GMT_LONG n_output = 0, n_expected_fields = 0;
	GMT_LONG nm, n_points = 0, n_read = 0;

	float *f = NULL;

	double value, west, east, south, north, x, y, *in = NULL, *out = NULL;

	FILE *fp = NULL;

	struct GRD_HEADER grd;
	struct GMT_EDGEINFO edgeinfo;
	struct GMT_BCR bcr;
	struct GRDTRACK_CTRL *Ctrl = NULL;

	void *New_grdtrack_Ctrl (), Free_grdtrack_Ctrl (struct GRDTRACK_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDTRACK_CTRL *) New_grdtrack_Ctrl ();	/* Allocate and initialize a new control structure */
	
	west = east = south = north = 0.0;
	memset ((void *)line, 0, (size_t)BUFSIZ);
	out = (double *)NULL;

	GMT_boundcond_init (&edgeinfo);

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
				case 'f':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */

				case 'G':
					Ctrl->G.active = TRUE;
					if (strchr (argv[i], ',') && !strchr (argv[i], '?')) {	/* IMG grid file with required parameters */
						if ((j = sscanf (&argv[i][2], "%[^,],%lf,%" GMT_LL "d,%lf", line, &Ctrl->G.scale, &Ctrl->G.mode, &Ctrl->G.lat)) < 3) {
							fprintf (stderr, "%s: SYNTAX ERROR -G option: Give imgfile, scale, mode [and optionally max_lat]\n", GMT_program);
							error = TRUE;
						}
						else
							Ctrl->G.file = strdup (line);
						Ctrl->G.type = 1;
					}
					else
						Ctrl->G.file = strdup (&argv[i][2]);
					break;
				case 'L':
					if (argv[i][2]) {
						Ctrl->L.active = TRUE;
						strncpy (Ctrl->L.mode, &argv[i][2], (size_t)4);
					}
					else {
						GMT_io.in_col_type[GMT->common->t.active] = GMT_io.out_col_type[GMT->common->t.active] = GMT_IS_LON;
						GMT_io.in_col_type[1-GMT->common->t.active] = GMT_io.out_col_type[1-GMT->common->t.active] = GMT_IS_LAT;
						fprintf (stderr, "%s: Option -L with no argument is obsolete (but is processed correctly).  Please use -f instead\n", GMT_program);
					}
					break;
				case 'N':	/* Backwards compatible */
					Ctrl->Q.interpolant = BCR_NEARNEIGHBOR;
					fprintf (stderr, "%s: Warning: Option -N deprecated. Use -Qn instead.\n", GMT_program);
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
					Ctrl->Q.interpolant = BCR_BILINEAR;
					for (j = 2; j < 5 && argv[i][j]; j++) {
						switch (argv[i][j]) {
							case 'n':
								Ctrl->Q.interpolant = BCR_NEARNEIGHBOR; break;
							case 'l':
								Ctrl->Q.interpolant = BCR_BILINEAR; break;
							case 'b':
								Ctrl->Q.interpolant = BCR_BSPLINE; break;
							case 'c':
								Ctrl->Q.interpolant = BCR_BICUBIC; break;
							case '/':
							default:
								Ctrl->Q.threshold = atof (&argv[i][j]);
								if (j == 2 && Ctrl->Q.threshold < GMT_SMALL) {
									Ctrl->Q.interpolant = BCR_NEARNEIGHBOR;
									fprintf (stderr, "%s: Warning: Option -Q0 deprecated. Use -Qn instead.\n", GMT_program);
								}
								j = 5; break;
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
		else if ((fp = GMT_fopen (argv[i], GMT_io.r_mode)) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[i]);
			exit (EXIT_FAILURE);
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "grdtrack %s - Sampling of a 2-D gridded netCDF grid file along 1-D trackline\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdtrack <xyfile> -G<grdfile> [%s] [-L<flag>] [-Q[b|c|l|n][[/]<threshold>]]\n", GMT_H_OPT); 
		fprintf (stderr, "\t[%s] [-S] [-V] [-Z] [%s] [%s]\n\t[%s] [%s]\n", GMT_Rgeo_OPT, GMT_t_OPT, GMT_b_OPT, GMT_f_OPT, GMT_m_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<xyfile> is an multicolumn ASCII file with (lon,lat) in the first two columns\n");
		fprintf (stderr, "\t-G <grdfile> is the name of the 2-D binary data set to sample\n");
		fprintf (stderr, "\t   If the file is a Sandwell/Smith Mercator grid (IMG format) instead,\n");
		fprintf (stderr, "\t   append comma-separated scale (0.1 or 1), mode, and optionally max latitude [%g].  Modes are\n", GMT_IMG_MAXLAT_80);
		fprintf (stderr, "\t     0 = img file w/ no constraint code, interpolate to get data at track.\n");
                fprintf (stderr, "\t     1 = img file w/ constraints coded, interpolate to get data at track.\n");
                fprintf (stderr, "\t     2 = img file w/ constraints coded, gets data only at constrained points, NaN elsewhere.\n");
                fprintf (stderr, "\t     3 = img file w/ constraints coded, gets 1 at constraints, 0 elsewhere.\n");
                fprintf (stderr, "\t   For mode 2|3 you may want to consider the -Q<threshold> setting.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-L sets boundary conditions.  <flag> can be either\n");
		fprintf (stderr, "\t   g for geographic boundary conditions\n");
		fprintf (stderr, "\t   or one or both of\n");
		fprintf (stderr, "\t   x for periodic boundary conditions on x\n");
		fprintf (stderr, "\t   y for periodic boundary conditions on y\n");
		fprintf (stderr, "\t   [Default is natural conditions for grids and geographic for IMG grids]\n");
		fprintf (stderr, "\t-Q Quick mode, use bilinear rather than bicubic [Default] interpolation.\n");
		fprintf (stderr, "\t   Alternatively, select interpolation mode by adding b = B-spline, c = bicubic,\n");
		fprintf (stderr, "\t   l = bilinear, or n = nearest-neighbor.\n");
		fprintf (stderr, "\t   Optionally, append <threshold> in the range [0,1]. [Default = 1 requires all\n");
		fprintf (stderr, "\t   4 or 16 nodes to be non-NaN.], <threshold> = 0.5 will interpolate about 1/2 way\n");
		fprintf (stderr, "\t   from a non-NaN to a NaN node, while 0.1 will go about 90%% of the way, etc.\n");
		fprintf (stderr, "\t   -Q0 will return the value of the nearest node instead of interpolating (Same as -Qn).\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t-S Suppress output when result equals NaN\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t-Z only output z-values [Default gives all columns]\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2 input columns\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (Ctrl->Q.active && (Ctrl->Q.threshold < 0.0 || Ctrl->Q.threshold > 1.0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Q:  threshold must be in [0,1] range\n", GMT_program);
		error++;
	}
	if (Ctrl->G.mode < 0 || Ctrl->G.mode > 3) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G:  mode must be in 0-3 range\n", GMT_program);
		error++;
	}
	if (Ctrl->G.lat < 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G:  max latitude should be positive\n", GMT_program);
		error++;
	}
	if (!Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G:  Must specify input file\n", GMT_program);
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
	if (Ctrl->L.active && GMT_boundcond_parse (&edgeinfo, Ctrl->L.mode)) error++;

	if (error) exit (EXIT_FAILURE);

	if (edgeinfo.gn) {
		GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
		GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
	}
	
	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

	pure_ascii = !(GMT_io.binary[GMT_IN] || GMT_io.binary[GMT_OUT]);

	if (fp == NULL) {
		fp = GMT_stdin;
		if (gmtdefs.verbose) fprintf (stderr, "%s: Reads from standard input\n", GMT_program);
#ifdef SET_IO_MODE
		GMT_setmode (GMT_IN);
#endif
	}

	GMT_grd_init (&grd, argc, argv, FALSE);
	GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 2;
	if (Ctrl->G.type == 0) {	/* Regular GMT grids */
		GMT_err_fail (GMT_read_grd_info (Ctrl->G.file, &grd), Ctrl->G.file);

		if (west == east) {	/* No subset asked for */
			west = grd.x_min;
			east = grd.x_max;
			south = grd.y_min;
			north = grd.y_max;
		}
		else	/* Make sure the subset fits the grid layout */
			GMT_err_pass (GMT_adjust_loose_wesn (&west, &east, &south, &north, &grd), grd.name);
		
		nx = GMT_get_n (west, east, grd.x_inc, grd.node_offset);
		ny = GMT_get_n (south, north, grd.y_inc, grd.node_offset);
		mx = nx + 4;	my = ny + 4;
		nm = GMT_get_nm (mx, my);

		f = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);

		GMT_err_fail (GMT_read_grd (Ctrl->G.file, &grd, f, west, east, south, north, GMT_pad, FALSE), Ctrl->G.file);

		project_info.w = west;	project_info.e = east;
		project_info.s = south;	project_info.n = north;

	}
	else {	/* Sandwell/Smith Mercator grids */
		GMT_read_img (Ctrl->G.file, &grd, &f, west, east, south, north, Ctrl->G.scale, Ctrl->G.mode, Ctrl->G.lat, TRUE);
		if (GMT_360_RANGE (grd.x_max, grd.x_min)) GMT_boundcond_parse (&edgeinfo, "g");
		GMT_boundcond_parse (&edgeinfo, "g");
		mx = grd.nx + 4;	my = grd.ny + 4;
	}
	
	GMT_boundcond_param_prep (&grd, &edgeinfo);

	/* Initialize bcr structure:  */

	GMT_bcr_init (&grd, GMT_pad, Ctrl->Q.interpolant, Ctrl->Q.threshold, &bcr);

	/* Set boundary conditions  */

	GMT_boundcond_set (&grd, &edgeinfo, GMT_pad, f);
	
	if (GMT_io.io_header[GMT_IN]) {	/* First echo headers, if any */
		for (i = 0; i < GMT_io.n_header_recs - 1; i++) {
			not_used = GMT_fgets (line, BUFSIZ, fp);
			if (!GMT_io.binary[GMT_OUT] && GMT_io.io_header[GMT_OUT]) fprintf (stdout, "%s", line);
		}
		not_used = GMT_fgets (line, BUFSIZ, fp);
		GMT_chop (line);
		if (!GMT_io.binary[GMT_OUT] && GMT_io.io_header[GMT_OUT]) fprintf (stdout, "%s\tsample\n", line);
	}

	ix = (gmtdefs.xy_toggle[GMT_OUT]);	iy = 1 - ix;	/* These are used for output purposes only */
	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : GMT_MAX_COLUMNS;
	memset ((void *)line, 0, (size_t)BUFSIZ);
	if (Ctrl->Z.active) GMT_io.out_col_type[0] = GMT_IS_FLOAT;	/* Since we are outputting z in the longitude column */

	while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {

		n_read++;

		while (GMT_io.status & GMT_IO_SEGMENT_HEADER) {
			GMT_write_segmentheader (GMT_stdout, n_expected_fields);
			n_fields = GMT_input (fp,  &n_expected_fields, &in);
		}
		if ((GMT_io.status & GMT_IO_EOF)) continue;	/* At EOF */

		if ((GMT_io.status & GMT_IO_MISMATCH) && n_fields < 2) {
			fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
			continue;
		}
		if (n_output == 0) n_output = n_expected_fields + 1;

		if (Ctrl->G.type == 1) {	/* Mercator IMG grid - get Mercator coordinates x,y */
			GMT_geo_to_xy (in[GMT_X], in[GMT_Y], &x, &y);
			if (x > grd.x_max) x -= 360.0;
		}
		else {			/* Regular grd, just copy the x,y */
			x = in[GMT_X];
			y = in[GMT_Y];
		}

		/* If point is outside grd area, shift it using periodicity or skip if not periodic. */

		while ( (in[GMT_Y] < grd.y_min) && (edgeinfo.nyp > 0) ) y += (grd.y_inc * edgeinfo.nyp);
		if (y < grd.y_min) continue;

		while ( (y > grd.y_max) && (edgeinfo.nyp > 0) ) y -= (grd.y_inc * edgeinfo.nyp);
		if (y > grd.y_max) continue;

		if (GMT_io.in_col_type[0] == GMT_IS_LON) {
			while (x > grd.x_max) x -= 360.0;
			while (x < grd.x_min) x += 360.0;
		}
		
		while ( (x < grd.x_min) && (edgeinfo.nxp > 0) ) x += (grd.x_inc * edgeinfo.nxp);
		if (x < grd.x_min) continue;

		while ( (x > grd.x_max) && (edgeinfo.nxp > 0) ) x -= (grd.x_inc * edgeinfo.nxp);
		if (x > grd.x_max) continue;

		value = GMT_get_bcr_z (&grd, x, y, f, &edgeinfo, &bcr);

		if (Ctrl->S.active && GMT_is_dnan (value)) continue;

		if (!out) out = (double *) GMT_memory (VNULL, (size_t)n_output, sizeof (double), GMT_program);

		if (Ctrl->Z.active) {	/* Simply print out value */
			GMT_output (GMT_stdout, 1, &value);
		}
		else if (pure_ascii && n_expected_fields > 2) {
			/* Special case: Ascii i/o and at least 3 columns:
			   Columns beyond first two could be text strings */

			/* First get rid of any commas that may cause grief */
			for (i = 0; GMT_io.current_record[i]; i++) if (GMT_io.current_record[i] == ',') GMT_io.current_record[i] = ' ';
			sscanf (GMT_io.current_record, "%*s %*s %[^\n]", line);
			GMT_ascii_output_one (GMT_stdout, in[ix], ix);	GMT_fputs ("\t", GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, in[iy], iy);	GMT_fputs ("\t", GMT_stdout);
			GMT_fputs (line, GMT_stdout);			GMT_fputs ("\t", GMT_stdout); 
			GMT_ascii_output_one (GMT_stdout, value, 2);	GMT_fputs ("\n", GMT_stdout);
		}
		else {	/* Simply copy other columns, append value, and output */
			for (i = 0; i < n_expected_fields; i++) out[i] = in[i];
			out[i] = value;
			GMT_output (GMT_stdout, n_output, out);
		}

		n_points++;
	}
	GMT_fclose (fp);

	if (gmtdefs.verbose) fprintf (stderr, "%s: Sampled %ld points from grid %s (%d x %d)\n", GMT_program,
		n_points, Ctrl->G.file, grd.nx, grd.ny);

	GMT_free ((void *)f);
	GMT_free ((void *)out);

	Free_grdtrack_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_grdtrack_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDTRACK_CTRL *C;
	
	C = (struct GRDTRACK_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDTRACK_CTRL), "New_grdtrack_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->G.scale = 1.0;
	C->Q.interpolant = BCR_BICUBIC; C->Q.threshold = 1.0;
	return ((void *)C);
}

void Free_grdtrack_Ctrl (struct GRDTRACK_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}
