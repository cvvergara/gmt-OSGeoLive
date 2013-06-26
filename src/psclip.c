/*--------------------------------------------------------------------
 *	$Id: psclip.c 9923 2012-12-18 20:45:53Z pwessel $
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
 * psclip reads one or many xy-files and draws polygons just like psxy
 * with the exception that these polygons are set up to act as clip
 * paths for subsequent plotting.  psclip uses the even-odd winding
 * rule to decide what's inside and outside the path.
 *
 * Author:	Paul Wessel
 * Date:	6-FEB-1991-2000
 * Version:	4
 *
 */

#include "gmt.h"
#include "pslib.h"

struct PSCLIP_CTRL {
	struct C {	/* -C */
		GMT_LONG active;
	} C;
	struct E {	/* -Eazim/elev */
		GMT_LONG active;
		double azimuth, elevation;
	} E;
	struct N {	/* -N */
		GMT_LONG active;
	} N;
	struct T {	/* -T */
		GMT_LONG active;
	} T;
	struct Z {	/* -Z<z_level> */
		GMT_LONG active;
		double level;
	} Z;
};

int main (int argc, char **argv)
{
	GMT_LONG i, n, n_alloc;
	GMT_LONG fno, n_args, n_files = 0, n_read, n_fields, n_expected_fields = 0;

	GMT_LONG error = FALSE, nofile, done, first;

	char line[BUFSIZ];

	double west, east, north, south, x0, y0, *in = NULL, *xx = NULL, *yy = NULL;

	FILE *fp = NULL;

	struct PSCLIP_CTRL *Ctrl = NULL;

	void *New_psclip_Ctrl (), Free_psclip_Ctrl (struct PSCLIP_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct PSCLIP_CTRL *)New_psclip_Ctrl ();	/* Allocate and initialize a new control structure */

	west = east = north = south = 0.0;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'B':
				case 'H':
				case 'J':
				case 'K':
				case 'M':
				case 'O':
				case 'P':
				case 'R':
				case 'U':
				case 'V':
				case 'X':
				case 'x':
				case 'Y':
				case 'y':
				case 'b':
				case 'c':
				case 'f':
				case 'm':
				case ':':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */

				case 'S':	/* Undocumented backwards compatible option */
				case 'C':
					Ctrl->C.active = TRUE;
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					sscanf (&argv[i][2], "%lf/%lf", &Ctrl->E.azimuth, &Ctrl->E.elevation);
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					break;
				case 'Z':
					if (argv[i][2]) {
						Ctrl->Z.level = atof (&argv[i][2]);
						Ctrl->Z.active = TRUE;
					}
					else {
						error++;
						fprintf (stderr, "%s: GMT SYNTAX ERROR -Z:  Must append a new z-value [0]\n", GMT_program);
					}
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
		fprintf (stderr,"psclip %s - To set up polygonal clip paths\n\n", GMT_VERSION);
		fprintf (stderr, "usage: psclip -C [-K] [-O]  OR\n");
		fprintf (stderr, "\tpsclip <xy-files> %s %s [%s]\n", GMT_J_OPT, GMT_Rgeo_OPT, GMT_B_OPT);
		fprintf (stderr, "\t[-Eaz/el] [%s] [-K] [-N] [-O] [-P] [-T] [%s] [-V]\n", GMT_H_OPT, GMT_U_OPT);
		fprintf (stderr, "\t[%s] [%s] [-Z<zlevel>] [%s] [%s]\n\t[%s] [%s] [%s]\n\n", GMT_X_OPT, GMT_Y_OPT, GMT_c_OPT, GMT_t_OPT, GMT_bi_OPT, GMT_f_OPT, GMT_mo_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t-C means stop existing clip-path.  -R, -J are not required\n");
		fprintf (stderr, "\t<xy-files> is one or more polygon files.  If none, standard input is read\n");
		GMT_explain_option ('j');
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_explain_option ('b');
		fprintf (stderr, "\t-E set azimuth and elevation of viewpoint for 3-D perspective [180/90]\n");
		GMT_explain_option ('H');
		GMT_explain_option ('K');
		fprintf (stderr, "\t-N uses the outside of the polygons as clip paths\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "\t-T Set of clip path for entire map frame.  No input file is required\n");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_explain_option ('X');
		fprintf (stderr, "\t-Z For 3-D plots: Set the z-level of map [0]\n");
		GMT_explain_option ('c');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2 input columns\n");
		GMT_explain_option ('f');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (!Ctrl->C.active) {
		if (!project_info.region_supplied) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
			error++;
		}
		if (Ctrl->E.elevation <= 0.0 || Ctrl->E.elevation > 90.0) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Elevation must be in 0-90 range\n", GMT_program);
			error++;
		}
	}
	if (Ctrl->T.active) Ctrl->N.active = TRUE;	/* -T implies -N */
	if (Ctrl->T.active && n_files) fprintf (stderr, "%s: Warning:  Option -T ignores all input files\n", GMT_program);

	if (Ctrl->N.active && frame_info.plot) {
		fprintf (stderr, "%s: Warning:  Option -B cannot be used in combination with Options -N or -T. -B is ignored.\n", GMT_program);
		frame_info.plot = FALSE;
	}

	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 2) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Binary input data (-bi) must have at least 2 columns\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

	z_project.view_azimuth = Ctrl->E.azimuth;
	z_project.view_elevation = Ctrl->E.elevation;

	if (!project_info.x_off_supplied && GMT_ps.overlay) GMT_ps.x_origin = 0.0;
	if (!project_info.y_off_supplied && GMT_ps.overlay) GMT_ps.y_origin = 0.0;

	if (Ctrl->C.active)
		GMT_ps.clip = -1;	/* Signal that this program terminates clipping that initiated prior to this process */
	else
		GMT_ps.clip = +1;	/* Signal that this program initiates clipping that wil outlive this process */

	GMT_plotinit (argc, argv);

	if (Ctrl->C.active && !frame_info.plot) {	/* Just undo previous clip-path, no basemap needed */
		ps_clipoff ();
		GMT_plotend ();

		if (gmtdefs.verbose) fprintf (stderr, "%s: Done!\n", GMT_program);

		Free_psclip_Ctrl (Ctrl);	/* Deallocate control structure */
		GMT_end (argc, argv);
		exit (EXIT_SUCCESS);
	}

	GMT_err_fail (GMT_map_setup (west, east, south, north), "");

	if (project_info.three_D) ps_transrotate (-z_project.xmin, -z_project.ymin, 0.0);

	if (Ctrl->Z.active) project_info.z_level = Ctrl->Z.level;

	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 2;

	if (Ctrl->C.active) {	/* Undo previous clip-path and draw basemap */
		ps_clipoff ();
		GMT_map_basemap ();
	}

	else {	/* Start new clip_path */

		GMT_map_basemap ();
		/* Allocate arrays holding the contour xy values */
		xx = (double *) GMT_memory (VNULL, (size_t)GMT_CHUNK, sizeof (double), GMT_program);
		yy = (double *) GMT_memory (VNULL, (size_t)GMT_CHUNK, sizeof (double), GMT_program);
		n_alloc = GMT_CHUNK;

		if (Ctrl->N.active) GMT_map_clip_on (GMT_no_rgb, 1);	/* Must clip map */
		first = !Ctrl->N.active;
		done = Ctrl->T.active;	/* When -T is used, skip all input */
		n_args = (argc > 1) ? argc : 2;
		nofile = (n_files == 0);

		for (fno = 1; !done && fno < n_args; fno++) {	/* Loop over all input files */
			if (!nofile && argv[fno][0] == '-') continue;
			if (nofile) {
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

			if (!nofile && gmtdefs.verbose) {
				fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[fno]);
				sprintf (line, "File: %s", argv[fno]);
				ps_comment (line);
			}
			if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) GMT_fgets (line, BUFSIZ, fp);
			if (GMT_io.multi_segments[GMT_IN]) {
				GMT_fgets (line, BUFSIZ, fp);
				if (gmtdefs.verbose) ps_comment (line);
			}
			n_fields = GMT_input (fp, &n_expected_fields, &in);

			while (! (GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

				while (GMT_io.status & GMT_IO_SEGMENT_HEADER) {
					if (gmtdefs.verbose) ps_comment (GMT_io.segment_header);
					n_fields = GMT_input (fp, &n_expected_fields, &in);
				}
				if ((GMT_io.status & GMT_IO_EOF)) continue;	/* At EOF */

				/* Start a new segment path */

				n_read = n = 0;

				while (! (GMT_io.status & (GMT_IO_SEGMENT_HEADER | GMT_IO_EOF))) {	/* Keep going until FALSE */
					n_read++;
					if (GMT_io.status & GMT_IO_MISMATCH) {
						fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%d) fields near line %ld (skipped)\n", GMT_program, n_fields, 2, n_read);
						continue;
					}

					xx[n] = in[GMT_X];	yy[n] = in[GMT_Y];
					n++;
					if (n == n_alloc) {
						n_alloc <<= 1;
						xx = (double *) GMT_memory ((void *)xx, (size_t)n_alloc, sizeof (double), GMT_program);
						yy = (double *) GMT_memory ((void *)yy, (size_t)n_alloc, sizeof (double), GMT_program);
					}
					n_fields = GMT_input (fp, &n_expected_fields, &in);
				}

				xx = (double *) GMT_memory ((void *)xx, (size_t)n, sizeof (double), GMT_program);
				yy = (double *) GMT_memory ((void *)yy, (size_t)n, sizeof (double), GMT_program);
				n_alloc = n;

				for (i = 0; i < n; i++) {
					GMT_geo_to_xy (xx[i], yy[i], &x0, &y0);
					xx[i] = x0; yy[i] = y0;
				}

				if (project_info.three_D) GMT_2D_to_3D (xx, yy, project_info.z_level, n);

				ps_clipon (xx, yy, n, GMT_no_rgb, first);
				first = FALSE;
			}
			if (fp != GMT_stdin) GMT_fclose (fp);
		}

		/* Finish up the composite clip path */
		ps_clipon (xx, yy, (GMT_LONG)0, GMT_no_rgb, 2 + first);

		GMT_free ((void *)xx);
		GMT_free ((void *)yy);
	}

	if (project_info.three_D) ps_rotatetrans (z_project.xmin, z_project.ymin, 0.0);

	GMT_plotend ();

	if (gmtdefs.verbose) fprintf (stderr, "%s: Done!\n", GMT_program);

	Free_psclip_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_psclip_Ctrl () {	/* Allocate and initialize a new control structure */
	struct PSCLIP_CTRL *C;
	
	C = (struct PSCLIP_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSCLIP_CTRL), "New_psclip_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->E.azimuth = 180.0;
	C->E.elevation = 90.0;
		
	return ((void *)C);
}

void Free_psclip_Ctrl (struct PSCLIP_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);	
}
