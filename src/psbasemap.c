/*--------------------------------------------------------------------
 *	$Id: psbasemap.c 9923 2012-12-18 20:45:53Z pwessel $
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
 * psbasemap plots out a basemap for the given area using the specified map
 * projection.
 *
 * Author:	Paul Wessel
 * Date:	08-JUL-2000
 * Version:	4
 *
 */
 
#include "gmt.h"
#include "pslib.h"

struct PSBASEMAP_CTRL {
	struct E {	/* -Eazim/elev */
		GMT_LONG active;
		double azimuth, elevation;
	} E;
	struct G {	/* -G<fill> */
		GMT_LONG active;
		struct GMT_FILL fill;
	} G;
	struct L {	/* -L */
		GMT_LONG active;
		struct GMT_MAP_SCALE item;
	} L;
	struct T {	/* -L */
		GMT_LONG active;
		struct GMT_MAP_ROSE item;
	} T;
	struct Z {	/* -Z<z_level> */
		GMT_LONG active;
		double level;
	} Z;
};

int main (int argc, char **argv)
{
	GMT_LONG i;

	GMT_LONG error = FALSE;

	double w = 0.0, e = 0.0, s = 0.0, n = 0.0;

	struct PSBASEMAP_CTRL *Ctrl = NULL;

	void *New_psbasemap_Ctrl (), Free_psbasemap_Ctrl (struct PSBASEMAP_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct PSBASEMAP_CTRL *)New_psbasemap_Ctrl ();	/* Allocate and initialize a new control structure */

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
				case '\0':
					error += GMT_parse_common_options (argv[i], &w, &e, &s, &n);
					break;

				/* Supplemental options */

				case 'E':
					Ctrl->E.active = TRUE;
					error += GMT_get_proj3D (&argv[i][2], &Ctrl->E.azimuth, &Ctrl->E.elevation);
					break;

				case 'G':
					Ctrl->G.active = TRUE;
					if (GMT_getfill (&argv[i][2], &Ctrl->G.fill)) {
						GMT_fill_syntax ('G', " ");
						error++;
					}
					break;

				case 'L':
					Ctrl->L.active = TRUE;
					error += GMT_getscale (&argv[i][2], &Ctrl->L.item);
					break;

				case 'T':
					Ctrl->T.active = TRUE;
					error += GMT_getrose (&argv[i][2], &Ctrl->T.item);
					break;

				case 'Z':
					if (argv[i][2]) {
						Ctrl->Z.level = atof (&argv[i][2]);
						Ctrl->Z.active = TRUE;
					}
					else {
						error++;
						fprintf (stderr, "%s: GMT SYNTAX ERROR -Z:  Must append a z-value\n", GMT_program);
					}
					break;

				/* Illegal options */

				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			fprintf (stderr, "%s: Warning: Ignoring filename %s\n", GMT_program, argv[i]);
	}

	if (GMT_give_synopsis_and_exit || argc == 1) {
		fprintf (stderr,"psbasemap %s - To plot PostScript basemaps\n\n", GMT_VERSION);
		fprintf (stderr, "usage: psbasemap %s %s %s [%s] [-G<fill>]\n", GMT_B_OPT, GMT_J_OPT, GMT_Rgeoz_OPT, GMT_E_OPT);
		fprintf (stderr, "\t[-K] [%s] [%s]\n", GMT_Jz_OPT, GMT_SCALE);
		fprintf (stderr, "\t[-O] [-P] [%s] [%s] [-V]\n", GMT_TROSE, GMT_U_OPT);
		fprintf (stderr, "\t[%s] [%s] [-Z<zlevel>] [%s]\n\n", GMT_X_OPT, GMT_Y_OPT, GMT_c_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		GMT_explain_option ('B');
		GMT_explain_option ('J');
		GMT_explain_option ('Z');
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_explain_option ('E');
		GMT_fill_syntax ('G', "Select fill inside of basemap.");
		GMT_explain_option ('K');
		GMT_mapscale_syntax ('L', "Draws a simple map scale centered on <lon0>/<lat0>.");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		GMT_maprose_syntax ('T', "Draws a north-pointing map rose centered on <lon0>/<lat0>.");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_explain_option ('X');
		fprintf (stderr, "\t-Z For 3-D plots: Set the z-level of map [Default is at bottom of z-axis].\n");
		GMT_explain_option ('c');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (!(frame_info.plot || Ctrl->L.active || Ctrl->T.active || Ctrl->G.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify at least one of -B, -G, -L, -T\n", GMT_program);
		error++;
	}
	if (Ctrl->E.elevation <= 0.0 || Ctrl->E.elevation > 90.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Elevation must be in 0-90 range\n", GMT_program);
		error++;
	}
	if (Ctrl->L.active && ! GMT_IS_MAPPING) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  -L applies to geographical data only\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	z_project.view_azimuth = Ctrl->E.azimuth;
	z_project.view_elevation = Ctrl->E.elevation;

	if (gmtdefs.verbose) fprintf (stderr, "psbasemap: Constructing basemap\n");

	GMT_err_fail (GMT_map_setup (w, e, s, n), "");

	GMT_plotinit (argc, argv);

	if (project_info.three_D) ps_transrotate (-z_project.xmin, -z_project.ymin, 0.0);

	if (Ctrl->G.active) {
		double *x, *y;
		GMT_LONG np;
		GMT_LONG donut;
		np = GMT_map_clip_path (&x, &y, &donut);
		GMT_fill (x, y, (1 + donut) * np, &Ctrl->G.fill, FALSE);
		GMT_free ((void *)x);
		GMT_free ((void *)y);
	}

	if (Ctrl->Z.active) project_info.z_level = Ctrl->Z.level;

	GMT_map_basemap ();

	if (Ctrl->L.active) GMT_draw_map_scale (&Ctrl->L.item);

	if (Ctrl->T.active) GMT_draw_map_rose (&Ctrl->T.item);

	if (project_info.three_D) ps_rotatetrans (z_project.xmin, z_project.ymin, 0.0);
	
	GMT_plotend ();

	Free_psbasemap_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_psbasemap_Ctrl () {	/* Allocate and initialize a new control structure */
	struct PSBASEMAP_CTRL *C;
	
	C = (struct PSBASEMAP_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSBASEMAP_CTRL), "New_psbasemap_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->E.azimuth = 180.0;
	C->E.elevation = 90.0;
	GMT_init_fill (&C->G.fill, -1, -1, -1);
	memset ((void *)&C->L.item, 0, sizeof (struct GMT_MAP_SCALE));
	memset ((void *)&C->T.item, 0, sizeof (struct GMT_MAP_ROSE));
		
	return ((void *)C);
}

void Free_psbasemap_Ctrl (struct PSBASEMAP_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);	
}
