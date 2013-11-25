/*--------------------------------------------------------------------
 *	$Id: pscoast.c 9975 2013-01-10 20:16:44Z pwessel $
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
 * pscoast can draw shorelines, paint the landward and/or seaward side of
 * the shoreline, paint lakes, draw political borders, and rivers.  Three
 * data sets are considered: shores, borders, and rivers.  Each of these
 * data sets come in 5 different resolutions.  The lower resolution files
 * are derived from the full resolution data using the Douglas-Peucker (DP)
 * line reduction algorithm.  By giving a tolerance in km, the algorithm
 * will remove points that do not depart more than the tolerance from the
 * straight line that would connect the points if the point in question was
 * removed.  The resolutions are:
 *
 * full		- The complete World Vector Shoreline + CIA data base
 * high		- DP reduced with tolerance = 0.2 km
 * intermediate	- DP reduced with tolerance = 1 km
 * low		- DP reduced with tolerance = 5 km
 * crude	- DP reduced with tolerance = 25 km
 *
 * If selected, pscoast will open the desired binned shoreline file and fill
 * the landmasses and/or oceans/lakes as shaded/colored/textured polygons. The
 * shoreline may also be drawn.  Political boundaries and rivers may be plotted,
 * too.  If only land is filled, then the 'wet' areas remain transparent.
 * If only oceans are filled, then the 'dry' areas remain transparent.  This
 * allows the user to overlay land or ocean without wiping out the plot.
 * For more information about the binned polygon file, see the GMT Technical
 * Reference and Cookbook.
 * Optionally, the user may choose to issue clip paths rather than paint the
 * polygons.  That way one may clip subsequent images to be visible only
 * inside or outside the coastline.
 *
 *
 * Author:	Paul Wessel
 * Date:	24-JUN-2000
 * Version:	4
 *
 */

#include "gmt.h"
#include "pslib.h"

#define LAKE	0
#define RIVER	1

struct PSCOAST_CTRL {
	struct A {	/* -A<min_area>[/<min_level>/<max_level>] */
		GMT_LONG active;
		struct GMT_SHORE_SELECT info;
	} A;
	struct C {	/* -C<fill> */
		GMT_LONG active;
		struct GMT_FILL fill[2];	/* lake and riverlake fill */
	} C;
	struct D {	/* -D<resolution> */
		GMT_LONG active;
		GMT_LONG force;	/* if TRUE, select next highest level if current set is not avaialble */
		char set;	/* One of f, h, i, l, c */
	} D;
	struct E {	/* -Eazim/elev */
		GMT_LONG active;
		double azimuth, elevation;
	} E;
	struct G {	/* -G<fill> */
		GMT_LONG active;
		GMT_LONG clip;
		struct GMT_FILL fill;
	} G;
	struct I {	/* -I<feature>[/<pen>] */
		GMT_LONG active;
		GMT_LONG use[GMT_N_RLEVELS];
		struct GMT_PEN pen[GMT_N_RLEVELS];
	} I;
	struct L {	/* -L */
		GMT_LONG active;
		struct GMT_MAP_SCALE item;
	} L;
	struct N {	/* -N<feature>[/<pen>] */
		GMT_LONG active;
		GMT_LONG use[GMT_N_BLEVELS];
		struct GMT_PEN pen[GMT_N_BLEVELS];
	} N;
	struct Q {	/* -Q */
		GMT_LONG active;
	} Q;
	struct S {	/* -S<fill> */
		GMT_LONG active;
		GMT_LONG clip;
		struct GMT_FILL fill;
	} S;
	struct T {	/* -L */
		GMT_LONG active;
		struct GMT_MAP_ROSE item;
	} T;
	struct W {	/* -W[<feature>/]<pen> */
		GMT_LONG active;
		GMT_LONG use[GMT_MAX_GSHHS_LEVEL];
		struct GMT_PEN pen[GMT_MAX_GSHHS_LEVEL];
	} W;
	struct Z {	/* -Z<z_level> */
		GMT_LONG active;
		double level;
	} Z;
#ifdef DEBUG
	struct DBG {	/* -+<bin> */
		GMT_LONG active;
		int bin;
	} debug;
#endif
};

int main (int argc, char **argv)
{
	int i, np, ind, bin, base, anti_bin = -1, np_new, k, last_k, err;
	int level_to_be_painted = 0, direction, start_direction, stop_direction, last_pen_level;
	int n_blevels = 0, n_rlevels = 0, bin_trouble;

	GMT_LONG n;
	
	GMT_LONG	error = FALSE, shift = FALSE, need_coast_base, recursive;
	GMT_LONG	greenwich = FALSE, possibly_donut_hell = FALSE, fill_in_use = FALSE;
	GMT_LONG clobber_background, paint_polygons = FALSE, donut, dumping = FALSE;
	GMT_LONG donut_hell = FALSE, world_map_save, clipping = FALSE;

	double	west = 0.0, east = 0.0, south = 0.0, north = 0.0, bin_x[5], bin_y[5], out[2];
	double west_border, east_border, anti_lon = 0.0, anti_lat = -90.0, edge = 720.0;
	double *xtmp = NULL, *ytmp = NULL, step, left, right, anti_x, anti_y, x_0, y_0, x_c, y_c, dist;

	char *string = NULL, comment[GMT_LONG_TEXT];
	char *shore_resolution[5] = {"full", "high", "intermediate", "low", "crude"};

	struct GMT_FILL fill[6];	/* Colors for (0) water, (1) land, (2) lakes, (3) islands in lakes, (4) lakes in islands in lakes, (5) riverlakes */
	struct GMT_PEN pen;
	struct GMT_SHORE c;
	struct GMT_BR b, r;
	struct GMT_GSHHS_POL *p = NULL;
	struct PSCOAST_CTRL *Ctrl = NULL;

	void recursive_path (int k0, int np, struct GMT_GSHHS_POL *p, int level, struct GMT_FILL fill[]);
	void *New_pscoast_Ctrl (), Free_pscoast_Ctrl (struct PSCOAST_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct PSCOAST_CTRL *) New_pscoast_Ctrl ();		/* Allocate and initialize defaults in a new control structure */

	/* Check and interpret the command line arguments */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {

				/* Common parameters */

				case 'B':	/* Get tickmark information */
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
				case 'c':
				case 'b':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */

				case 'A':
					Ctrl->A.active = TRUE;
					GMT_set_levels (&argv[i][2], &Ctrl->A.info);
					break;
				case 'C':	/* Lake colors */
					Ctrl->C.active = TRUE;
					if ((argv[i][2] == 'l' || argv[i][2] == 'r') && argv[i][3] == '/') {	/* Specific lake or river-lake fill */
						k = (argv[i][2] == 'l') ? LAKE : RIVER;
						if (argv[i][4] && GMT_getfill (&argv[i][4], &Ctrl->C.fill[k])) {
							GMT_fill_syntax ('C', " ");
							error++;
						}
					}
					else if (argv[i][2]) {
						if (GMT_getfill (&argv[i][2], &Ctrl->C.fill[LAKE])) {
							GMT_fill_syntax ('C', " ");
							error++;
						}
						Ctrl->C.fill[RIVER] = Ctrl->C.fill[LAKE];
					}
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					Ctrl->D.set = argv[i][2];
					if (argv[i][3] == '+') Ctrl->D.force = TRUE;
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					error += GMT_get_proj3D (&argv[i][2], &Ctrl->E.azimuth, &Ctrl->E.elevation);
					break;
				case 'G':		/* Set Gray shade, pattern, or clipping */
					Ctrl->G.active = TRUE;
					if (argv[i][2] == 'c' && argv[i][3] == '\0')
						Ctrl->G.clip = TRUE;
					else if (argv[i][2] && GMT_getfill (&argv[i][2], &Ctrl->G.fill)) {
						GMT_fill_syntax ('G', " ");
						error++;
					}
					break;
				case 'L':
					Ctrl->L.active = TRUE;
					error += GMT_getscale (&argv[i][2], &Ctrl->L.item);
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					if (!argv[i][2]) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR:  -N option takes at least one argument\n", GMT_program);
						error++;
						continue;
					}
					if ((string = strchr (argv[i], '/'))) {		/* Get specified pen */
						if (GMT_getpen (++string, &pen)) {	/* Error decoding pen */
							GMT_pen_syntax ('N', " ");
							error++;
						}
					}
					else	/* No pen specified, use default */
						GMT_init_pen (&pen, GMT_PENWIDTH);

					switch (argv[i][2]) {
						case 'a':
							for (k = 0; k < GMT_N_BLEVELS; k++) Ctrl->N.use[k] = TRUE;
							for (k = 0; k < GMT_N_BLEVELS; k++) Ctrl->N.pen[k] = pen;
							break;
						default:
							k = argv[i][2] - '1';
							if (k < 0 || k >= GMT_N_BLEVELS) {
								fprintf (stderr, "%s: GMT SYNTAX ERROR -N option: Feature not in list!\n", GMT_program);
								error++;
							}
							else {
								Ctrl->N.use[k] = TRUE;
								Ctrl->N.pen[k] = pen;
							}
							break;
					}
					break;

				case 'I':
					Ctrl->I.active = TRUE;
					if (!argv[i][2]) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR:  -I option takes at least one argument\n", GMT_program);
						error++;
						continue;
					}
					if ((string = strchr (argv[i], '/'))) {	/* Get specified pen */
						if (GMT_getpen (++string, &pen)) {	/* Error decoding pen */
							GMT_pen_syntax ('I', " ");
							error++;
						}
					}
					else	/* No pen specified, use default */
						GMT_init_pen (&pen, GMT_PENWIDTH);

					switch (argv[i][2]) {
						case 'a':
							for (k = 0; k < GMT_N_RLEVELS; k++) Ctrl->I.use[k] = TRUE;
							for (k = 0; k < GMT_N_RLEVELS; k++) Ctrl->I.pen[k] = pen;
							break;
						case 'r':
							for (k = 0; k < GMT_RIV_INTERMITTENT; k++) Ctrl->I.use[k] = TRUE;
							for (k = 0; k < GMT_RIV_INTERMITTENT; k++) Ctrl->I.pen[k] = pen;
							break;
						case 'i':
							for (k = GMT_RIV_INTERMITTENT; k < GMT_RIV_CANALS; k++) Ctrl->I.use[k] = TRUE;
							for (k = GMT_RIV_INTERMITTENT; k < GMT_RIV_CANALS; k++) Ctrl->I.pen[k] = pen;
							break;
						case 'c':
							for (k = GMT_RIV_CANALS; k < GMT_N_RLEVELS; k++) Ctrl->I.use[k] = TRUE;
							for (k = GMT_RIV_CANALS; k < GMT_N_RLEVELS; k++) Ctrl->I.pen[k] = pen;
							break;
						default:
							k = atoi (&argv[i][2]);
							if (k < 0 || k >= GMT_N_RLEVELS) {
								fprintf (stderr, "%s: GMT SYNTAX ERROR -I option: Feature not in list!\n", GMT_program);
								error++;
							}
							else {
								Ctrl->I.use[k] = TRUE;
								Ctrl->I.pen[k] = pen;
							}
							break;
					}
					break;

				case 'S':		/* Set ocean color if needed */
					Ctrl->S.active = TRUE;
					if (argv[i][2] == 'c' && argv[i][3] == '\0')
						Ctrl->S.clip = TRUE;
					else if (argv[i][2] && GMT_getfill (&argv[i][2], &Ctrl->S.fill)) {
						GMT_fill_syntax ('S', " ");
						error++;
					}
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					error += GMT_getrose (&argv[i][2], &Ctrl->T.item);
					break;
				case 'W':
					Ctrl->W.active = TRUE;	/* Want to draw shorelines */
					if ((argv[i][2] >= '1' && argv[i][2] <= '4') && argv[i][3] == '/') {	/* Specific pen for this feature */
						k = (GMT_LONG)(argv[i][2] - '1');
						if (argv[i][4] && GMT_getpen (&argv[i][4], &Ctrl->W.pen[k])) {
							GMT_pen_syntax ('W', " ");
							error++;
						}
						Ctrl->W.use[k] = TRUE;
					}
					else if (argv[i][2]) {	/* Same pen for all features */
						Ctrl->W.use[0] = TRUE;
						if (GMT_getpen (&argv[i][2], &Ctrl->W.pen[0])) {
							GMT_pen_syntax ('W', " ");
							error++;
						}
						else {
							for (k = 1; k < GMT_MAX_GSHHS_LEVEL; k++) {
								Ctrl->W.pen[k] = Ctrl->W.pen[0];
								Ctrl->W.use[k] = TRUE;
							}
						}
					}
					else {	/* Accept default pen for all features */
						for (k = 0; k < GMT_MAX_GSHHS_LEVEL; k++) Ctrl->W.use[k] = TRUE;
					}
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
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
#ifdef DEBUG
				case '+':
					Ctrl->debug.active = TRUE;
					Ctrl->debug.bin = atoi (&argv[i][2]);
					break;
#endif
				default:		/* Options not recognized */
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			fprintf (stderr, "%s: Warning: Ignoring filename %s\n", GMT_program, argv[i]);
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {	/* Display usage */
		fprintf (stderr,"pscoast %s - Plot continents, shorelines, rivers, and borders on maps\n\n", GMT_VERSION);
		fprintf (stderr,"usage: pscoast %s %s [%s] [%s]\n", GMT_B_OPT, GMT_J_OPT, GMT_A_OPT, GMT_Rgeoz_OPT);
		fprintf (stderr, "\t[-C[<feature>/]<fill>] [-D<resolution>][+] [%s] [-G[<fill>]] [-I<feature>[/<pen>]]\n", GMT_E_OPT);
		fprintf (stderr, "\t[%s] [-K] [%s]\n", GMT_Jz_OPT, GMT_SCALE);
		fprintf (stderr, "\t[-N<feature>[/<pen>]] [-O] [-P] [-Q] [-S<fill>]\n");
		fprintf (stderr, "\t[%s] [%s]\n", GMT_TROSE, GMT_U_OPT);
		fprintf (stderr, "\t[-V] [-W[<feature>/][<pen>]] [%s] [%s] [-Z<zlevel>]\n", GMT_X_OPT, GMT_Y_OPT);
		fprintf (stderr, "\t[%s] [%s] [%s]\n", GMT_bo_OPT, GMT_c_OPT, GMT_mo_OPT);
#ifdef DEBUG
		fprintf (stderr, "\t[-+<bin>]\n");
#endif

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		GMT_explain_option ('j');
		GMT_explain_option ('Z');
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_GSHHS_syntax ('A', "Place limits on coastline features from the GSHHS data base.");
		GMT_explain_option ('b');
		GMT_fill_syntax ('C', "Sets separate color for lakes and riverlakes [Default is same as ocean].");
		fprintf (stderr, "\t   Alternatively, set custom fills below.  Repeat the -C option as needed.\n");
		fprintf (stderr, "\t      l = Lakes\n");
		fprintf (stderr, "\t      r = River-lakes\n");
		fprintf (stderr, "\t-D Choose one of the following resolutions:\n");
		fprintf (stderr, "\t   f - full resolution (may be very slow for large regions).\n");
		fprintf (stderr, "\t   h - high resolution (may be slow for large regions).\n");
		fprintf (stderr, "\t   i - intermediate resolution.\n");
		fprintf (stderr, "\t   l - low resolution [Default].\n");
		fprintf (stderr, "\t   c - crude resolution, for busy plots that need crude continent outlines only.\n");
		fprintf (stderr, "\t   Append + to use a lower resolution should the chosen one not be available [abort].\n");
		GMT_explain_option ('E');
		GMT_fill_syntax ('G', "Paint or clip \"dry\" areas.");
		fprintf (stderr, "\t   6) c to issue clip paths for land areas.\n");
		GMT_pen_syntax ('I', "draws rivers.  Specify feature and optionally append pen [Default is solid line of unit thickness].");
		fprintf (stderr, "\t   Choose from the features below.  Repeat the -I option as many times as needed.\n");
		fprintf (stderr, "\t      0 = Double-lined rivers (river-lakes).\n");
		fprintf (stderr, "\t      1 = Permanent major rivers.\n");
		fprintf (stderr, "\t      2 = Additional major rivers.\n");
		fprintf (stderr, "\t      3 = Additional rivers.\n");
		fprintf (stderr, "\t      4 = Minor rivers.\n");
		fprintf (stderr, "\t      5 = Intermittent rivers - major.\n");
		fprintf (stderr, "\t      6 = Intermittent rivers - additional.\n");
		fprintf (stderr, "\t      7 = Intermittent rivers - minor.\n");
		fprintf (stderr, "\t      8 = Major canals.\n");
		fprintf (stderr, "\t      9 = Minor canals.\n");
		fprintf (stderr, "\t     10 = Irrigation canals.\n");
		fprintf (stderr, "\t      a = All rivers and canals (0-10).\n");
		fprintf (stderr, "\t      r = All permanent rivers (0-4).\n");
		fprintf (stderr, "\t      i = All intermittent rivers (5-7).\n");
		fprintf (stderr, "\t      c = All canals (8-10).\n");
		GMT_explain_option ('K');
		GMT_mapscale_syntax ('L', "Draws a simple map scale centered on <lon0>/<lat0>.");
		GMT_pen_syntax ('N', "draws boundaries.  Specify feature and optionally append pen [Default is solid line of unit thickness].");
		fprintf (stderr, "\t   Choose from the features below.  Repeat the -N option as many times as needed.\n");
		fprintf (stderr, "\t     1 = National boundaries.\n");
		fprintf (stderr, "\t     2 = State boundaries within the Americas.\n");
		fprintf (stderr, "\t     3 = Marine boundaries.\n");
		fprintf (stderr, "\t     a = All boundaries (1-3).\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr,"\t-Q means terminate previously set clip-paths.\n");
		GMT_fill_syntax ('S', "Paint of clip \"wet\" areas.");
		fprintf (stderr, "\t   6) c to issue clip paths for water areas.\n");
		GMT_maprose_syntax ('T', "Draws a north-pointing map rose centered on <lon0>/<lat0>.");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_pen_syntax ('W', "draw shorelines.  Append pen [Default is solid black line of unit thickness for all levels].");
		fprintf (stderr, "\t   Alternatively, set custom pens below.  Repeat the -W option as many times as needed.\n");
		fprintf (stderr, "\t      1 = Coastline.\n");
		fprintf (stderr, "\t      2 = Lake shores.\n");
		fprintf (stderr, "\t      3 = Island in lakes shores.\n");
		fprintf (stderr, "\t      4 = Lake in island in lake shores.\n");
		fprintf (stderr, "\t   When feature-specific pens are used, those not set are deactivated\n");
		GMT_explain_option ('X');
		fprintf (stderr, "\t-Z For 3-D plots: Set the z-level of map [0].\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('c');
		fprintf (stderr, "\t-m Dump a multisegment ascii (or binary, see -bo) file to standard output.  No plotting occurs\n");
		fprintf (stderr, "\t   Specify any combination of -W, -I, -N.  Optionally, append 1-char\n");
		fprintf (stderr, "\t   segment header <flag> [%c].\n", GMT_io.EOF_flag[GMT_OUT]);
#ifdef DEBUG
		fprintf (stderr, "\t-+ Print only a single bin (debug option).\n");
#endif
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (Ctrl->C.active && !(Ctrl->G.active || Ctrl->S.active || Ctrl->W.active)) {	/* Just lakes, fix -A */
		if (Ctrl->A.info.low < 2) Ctrl->A.info.low = 2;
	}

	/* Check that the options selected are mutually consistent */

	dumping = (GMT_io.multi_segments[GMT_IN] || GMT_io.multi_segments[GMT_OUT]);	/* -M was set */
	clipping = (Ctrl->G.clip || Ctrl->S.clip);
	if (!Ctrl->Q.active) {	/* Need -R -J */
		if (!project_info.region_supplied) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
			error++;
		}
		if (!GMT_IS_MAPPING) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Geographic coastline data require the use of a map projection (use -Jx<scale>d if linear scaling is desired).\n", GMT_program);
			error++;
		}
	}
	for (k = 0; k < GMT_MAX_GSHHS_LEVEL; k++) {
		if (Ctrl->W.pen[k].width < 0.0) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -W option:  Pen thickness for feature %d cannot be negative\n", GMT_program, k);
			error++;
		}
	}
	if (!(Ctrl->G.active || Ctrl->S.active || Ctrl->C.active || Ctrl->W.active || Ctrl->N.active || Ctrl->I.active || Ctrl->Q.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify at least one of -C, -G, -S, -I, -N, and -W\n", GMT_program);
		error++;
	}
	if ((Ctrl->G.active + Ctrl->S.active + Ctrl->C.active) > 1 && clipping > 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Cannot combine -C, -G, -S while clipping\n", GMT_program);
		error++;
	}
	if (Ctrl->G.clip && Ctrl->S.clip) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must choose between clipping land OR water\n", GMT_program);
		error++;
	}
	if (dumping && (Ctrl->G.active || Ctrl->S.active || Ctrl->C.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must choose between dumping and clipping/plotting\n", GMT_program);
		error++;
	}
	if (Ctrl->E.elevation <= 0.0 || Ctrl->E.elevation > 90.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Elevation must be in 0-90 range\n", GMT_program);
		error++;
	}
	if (clipping && project_info.projection == GMT_AZ_EQDIST && fabs (project_info.w - project_info.e) == 360.0 && (project_info.n - project_info.s) == 180.0) {
		fprintf (stderr, "%s: -JE not implemented for global clipping - I quit\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (clipping && (Ctrl->N.active || Ctrl->I.active || Ctrl->W.active)) {
		fprintf (stderr, "%s: Cannot do clipping AND draw coastlines, rivers, or borders\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (dumping && !(Ctrl->N.active || Ctrl->I.active || Ctrl->W.active)) {
		fprintf (stderr, "%s: Must specify at least one of -I, -N, and -W\n", GMT_program);
		exit (EXIT_FAILURE);
	}
		
	if (Ctrl->I.active) {
		for (k = n_rlevels = 0; k < GMT_N_RLEVELS; k++) {
			if (!Ctrl->I.use[k]) continue;
			if (Ctrl->I.pen[k].width < 0.0) {
				fprintf (stderr, "%s: GMT SYNTAX ERROR -I option:  Pen thickness cannot be negative\n", GMT_program);
				error++;
			}
			Ctrl->I.use[n_rlevels] = k;
			Ctrl->I.pen[n_rlevels] = Ctrl->I.pen[k];
			n_rlevels++;
		}
	}

	if (Ctrl->N.active) {
		for (k = n_blevels = 0; k < GMT_N_BLEVELS; k++) {
			if (!Ctrl->N.use[k]) continue;
			if (Ctrl->N.pen[k].width < 0.0) {
				fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Pen thickness cannot be negative\n", GMT_program);
				error++;
			}
			Ctrl->N.use[n_blevels] = k + 1;
			Ctrl->N.pen[n_blevels] = Ctrl->N.pen[k];
			n_blevels++;
		}
	}

	if (error) exit (EXIT_FAILURE);

	if (Ctrl->D.force) Ctrl->D.set = GMT_shore_adjust_res (Ctrl->D.set);
	base = (int)GMT_set_resolution (&Ctrl->D.set, 'D');
	z_project.view_azimuth = Ctrl->E.azimuth;
	z_project.view_elevation = Ctrl->E.elevation;
	fill[0] = Ctrl->S.fill;
	fill[1] = fill[3] = Ctrl->G.fill;
	fill[2] = fill[4] = (Ctrl->C.active) ? Ctrl->C.fill[LAKE] : Ctrl->S.fill;
	fill[5] = (Ctrl->C.active) ? Ctrl->C.fill[RIVER] : fill[2];
	need_coast_base = (Ctrl->G.active || Ctrl->S.active || Ctrl->C.active || Ctrl->W.active);
	if (Ctrl->Q.active) need_coast_base = FALSE;	/* Since we just end clipping */
	clobber_background = (Ctrl->G.active && Ctrl->S.active);
	recursive = (Ctrl->G.active != (Ctrl->S.active || Ctrl->C.active) || clipping);
	paint_polygons = (Ctrl->G.active || Ctrl->S.active || Ctrl->C.active);

	if (east > 360.0) {
		west -= 360.0;
		east -= 360.0;
	}

	if (Ctrl->Q.active && project_info.projection == GMT_NO_PROJ) {	/* Fake area and linear projection */
		GMT_parse_J_option ("x1d");
		east = north = 1.0;
	}
	else if (dumping) {
		GMT_parse_J_option ("x1d");
	}

	GMT_err_fail (GMT_map_setup (west, east, south, north), "");

	world_map_save = GMT_world_map;

	if (need_coast_base && GMT_init_shore (Ctrl->D.set, &c, project_info.w, project_info.e, project_info.s, project_info.n, &Ctrl->A.info))  {
		fprintf (stderr, "%s: %s resolution shoreline data base not installed\n", GMT_program, shore_resolution[base]);
		need_coast_base = FALSE;
	}

	if (Ctrl->N.active && GMT_init_br ('b', Ctrl->D.set, &b, project_info.w, project_info.e, project_info.s, project_info.n)) {
		fprintf (stderr, "%s: %s resolution political boundary data base not installed\n", GMT_program, shore_resolution[base]);
		Ctrl->N.active = FALSE;
	}
	if (need_coast_base && gmtdefs.verbose == 2) fprintf (stderr, "GSHHG version %s\n%s\n%s\n", c.version, c.title, c.source);

	if (Ctrl->I.active && GMT_init_br ('r', Ctrl->D.set, &r, project_info.w, project_info.e, project_info.s, project_info.n)) {
		fprintf (stderr, "%s: %s resolution river data base not installed\n", GMT_program, shore_resolution[base]);
		Ctrl->I.active = FALSE;
	}

	if (! (need_coast_base || Ctrl->N.active || Ctrl->I.active || Ctrl->Q.active)) {
		fprintf (stderr, "%s: No databases available - aborts\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (dumping) {
#ifdef SET_IO_MODE
		GMT_setmode (GMT_OUT);
#endif
		if (!GMT_io.binary[1]) {
			sprintf (comment, "# Data from the %s resolution GMT shoreline, borders, and rivers database\n", shore_resolution[base]);
			GMT_fputs (comment, GMT_stdout);
		}
	}
	else {
		if (Ctrl->Q.active)
			GMT_ps.clip = -1;	/* Signal that this program terminates clipping that initiated prior to this process */
		else if (clipping)
			GMT_ps.clip = +1;	/* Signal that this program initiates clipping that wil outlive this process */

		GMT_plotinit (argc, argv);

		if (Ctrl->Q.active) {  /* Just undo previous clip-path */
	  		ps_clipoff ();

			GMT_map_basemap (); /* Basemap needed */

			GMT_plotend ();

			if (gmtdefs.verbose) fprintf (stderr, "%s: Done!\n", GMT_program);
			
			Free_pscoast_Ctrl (Ctrl);	/* Deallocate control structure */

			GMT_end (argc, argv);
			exit (EXIT_SUCCESS);
		}

		if (project_info.three_D) ps_transrotate (-z_project.xmin, -z_project.ymin, 0.0);
	}

	for (i = 0; i < 5; i++) if (fill[i].use_pattern) fill_in_use = TRUE;

	if (fill_in_use && !clobber_background) {	/* Force clobber when fill is used since our routine cannot deal with clipped fills */
		if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: Pattern fill requires oceans to be painted first\n", GMT_program);
		clobber_background = TRUE;
		recursive = FALSE;
	}

	if (clipping && frame_info.plot) {	/* Want basemap drawn before clipping is activated */
		/* project_info.w = west_border;
		project_info.e = east_border; */
		GMT_map_basemap ();
	}

	if (project_info.projection == GMT_AZ_EQDIST && fabs (project_info.w - project_info.e) == 360.0 && (project_info.n - project_info.s) == 180.0) {
		possibly_donut_hell = TRUE;
		anti_lon = project_info.central_meridian + 180.0;
		if (anti_lon >= 360.0) anti_lon -= 360.0;
		anti_lat = -project_info.pole;
		anti_bin = (int)floor ((90.0 - anti_lat) / c.bsize) * c.bin_nx + (int)floor (anti_lon / c.bsize);
		GMT_geo_to_xy (anti_lon, anti_lat, &anti_x, &anti_y);
		GMT_geo_to_xy (project_info.central_meridian, project_info.pole, &x_0, &y_0);
		if (Ctrl->G.active && gmtdefs.verbose) {
			fprintf (stderr, "%s: Warning: Fill/cliiping continent option (-G) may not work for this projection.\n", GMT_program);
			fprintf (stderr, "If the antipole (%g/%g) is in the ocean then chances are good\n", anti_lon, anti_lat);
			fprintf (stderr, "Else: avoid projection center coordinates that are exact multiples of %g degrees\n", c.bsize);
		}
	}

	if (possibly_donut_hell && paint_polygons && !clobber_background) {	/* Force clobber when donuts may be called for now */
		if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: -JE requires oceans to be painted first\n", GMT_program);
		clobber_background = TRUE;
		recursive = FALSE;
	}

	if (clobber_background) {	/* Paint entire map as ocean first, then lay land on top */
		n = GMT_map_clip_path (&xtmp, &ytmp, &donut);
		if (donut) {
			ps_line (xtmp, ytmp, n, 1, TRUE);
			ps_line (&xtmp[n], &ytmp[n], n, 1, TRUE);
			GMT_fill (xtmp, ytmp, n, &Ctrl->S.fill, -9);
		}
		else
			GMT_fill (xtmp, ytmp, n, &Ctrl->S.fill, FALSE);
		GMT_free ((void *)xtmp);
		GMT_free ((void *)ytmp);
		level_to_be_painted = 1;
	}
	if (recursive) {
		start_direction = stop_direction = (Ctrl->G.active) ? 1 : -1;
		level_to_be_painted = (Ctrl->G.active) ? 1 : 0;
	}
	else {
		start_direction = -1;
		stop_direction = 1;
		level_to_be_painted = (Ctrl->S.active) ? 0 : 1;
	}

	if (west < 0.0 && east > 0.0 && !GMT_IS_LINEAR) greenwich = TRUE;
	if ((360.0 - fabs (project_info.e - project_info.w) ) < c.bsize) {
		GMT_world_map = TRUE;
		if (project_info.projection == GMT_GNOMONIC || project_info.projection == GMT_GENPER) GMT_world_map = FALSE;
		if (GMT_IS_AZIMUTHAL) GMT_world_map = FALSE;
	}
	if (GMT_world_map && greenwich)
		edge = project_info.central_meridian;
	else if (!GMT_world_map && greenwich) {
		shift = TRUE;
		edge = west;
		if (edge < 0.0) edge += 360.0;
		if (edge > 180.0) edge = 180.0;
	}

	if (project_info.w < 0.0 && project_info.e <= 0.0) {	/* Temporarily shift boundaries */
		project_info.w += 360.0;
		project_info.e += 360.0;
		if (project_info.central_meridian < 0.0) project_info.central_meridian += 360.0;
	}
	west_border = floor (project_info.w / c.bsize) * c.bsize;
	east_border = ceil (project_info.e / c.bsize) * c.bsize;

	if (!dumping && Ctrl->W.active) GMT_setpen (&Ctrl->W.pen[0]);

	if (Ctrl->Z.active) project_info.z_level = Ctrl->Z.level;

	/* Maximum step size (in degrees) used for interpolation of line segments along great circles */
	step = gmtdefs.line_step / project_info.x_scale / project_info.DIST_M_PR_DEG;
	last_pen_level = -1;

	if (clipping) ps_clipon ((double *)NULL, (double *)NULL, 0, GMT_no_rgb, 1);	/* Start clip path */

	for (ind = 0; need_coast_base && ind < c.nb; ind++) {	/* Loop over necessary bins only */

		bin = (int)c.bins[ind];
#ifdef DEBUG
		if (Ctrl->debug.active && bin != Ctrl->debug.bin) continue;
#endif
		if ((err = (int)GMT_get_shore_bin (ind, &c))) {
			fprintf (stderr, "%s: %s [%s resolution shoreline]\n", GMT_program, GMT_strerror(err), shore_resolution[base]);
			exit (EXIT_FAILURE);
		}

		if (gmtdefs.verbose) fprintf (stderr, "%s: Working on block # %5d\r", GMT_program, bin);
		if (!dumping) {
			sprintf (comment, "Bin # %d", bin);
			ps_comment (comment);
		}

		if (GMT_world_map && greenwich) {
			left = c.bsize * (bin % c.bin_nx);	right = left + c.bsize;
			shift = ((left - edge) <= 180.0 && (right - edge) > 180.0);
		}

		bin_trouble = (anti_bin == bin) ? anti_bin : -1;

		if (possibly_donut_hell) {
			bin_x[0] = bin_x[3] = bin_x[4] = c.lon_sw;
			bin_x[1] = bin_x[2] = c.lon_sw + c.bsize;
			bin_y[0] = bin_y[1] = bin_y[4] = c.lat_sw;
			bin_y[2] = bin_y[3] = c.lat_sw + c.bsize;
			GMT_geo_to_xy (c.lon_sw + 0.5 * c.bsize, c.lat_sw + 0.5 * c.bsize, &x_c, &y_c);
			dist = hypot (x_c - x_0, y_c - y_0);
			donut_hell = (dist > 0.8 * project_info.r || GMT_non_zero_winding (anti_lon, anti_lat, bin_x, bin_y, 5));
		}

		for (direction = start_direction; paint_polygons && direction <= stop_direction; direction += 2) {

			/* Assemble one or more segments into polygons */

			if ((np = (int)GMT_assemble_shore (&c, direction, TRUE, shift, west_border, east_border, &p)) == 0) continue;

			/* Get clipped polygons in x,y inches that can be plotted */

			np_new = (int)GMT_prep_polygons (&p, np, donut_hell, step, bin_trouble);

			if (clipping) {
				for (k = level_to_be_painted; k < GMT_MAX_GSHHS_LEVEL - 1; k++) recursive_path (-1, np_new, p, k, NULL);

				for (k = 0; k < np_new; k++) {	/* Do any remaining interior polygons */
					if (p[k].n == 0) continue;
					if (p[k].level%2 == level_to_be_painted || p[k].level > 2) {
						ps_line (p[k].lon, p[k].lat, p[k].n, 1, FALSE);
					}
				}
			}
			else if (recursive) {	/* Must avoid pointing anything but the polygons inside */

				for (k = level_to_be_painted; k < GMT_MAX_GSHHS_LEVEL - 1; k++) recursive_path (-1, np_new, p, k, fill);
				for (k = 0; k < np_new; k++) {	/* Do any remaining interior polygons */
					if (p[k].n == 0) continue;
					if (p[k].level%2 == level_to_be_painted || p[k].level > 2)
						GMT_fill (p[k].lon, p[k].lat, p[k].n, &fill[p[k].fid], FALSE);
				}
			}
			else {	/* Simply paints all polygons as is */
				for (k = 0; k < np_new; k++) {
					if (p[k].n == 0 || p[k].level < level_to_be_painted) continue;
					if (donut_hell && GMT_non_zero_winding (anti_x, anti_y, p[k].lon, p[k].lat, p[k].n)) {	/* Antipode inside polygon, must do donut */
						n = GMT_map_clip_path (&xtmp, &ytmp, &donut);
						ps_line (xtmp, ytmp, n, 1, TRUE);
						ps_line (p[k].lon, p[k].lat, p[k].n, 1, TRUE);
						GMT_fill (xtmp, ytmp, n, &fill[p[k].fid], -9);
						GMT_free ((void *)xtmp);
						GMT_free ((void *)ytmp);
					}
					else
						GMT_fill (p[k].lon, p[k].lat, p[k].n, &fill[p[k].fid], FALSE);
				}
			}

			GMT_free_polygons (p, np_new);
			GMT_free ((void *)p);
		}

		if (Ctrl->W.active && c.ns) {	/* Draw or dump shorelines, no need to assemble polygons */
			if ((np = (int)GMT_assemble_shore (&c, direction, FALSE, shift, west_border, east_border, &p)) == 0) continue;

			for (i = 0; i < np; i++) {
				if (dumping) {
					sprintf (GMT_io.segment_header, "%c Shore Bin # %d, Level %ld\n", GMT_io.EOF_flag[GMT_OUT], bin, p[i].level);
					GMT_write_segmentheader (GMT_stdout, 2);
					for (k = 0; k < p[i].n; k++) {
						out[GMT_X] = p[i].lon[k];
						out[GMT_Y] = p[i].lat[k];
						GMT_output (GMT_stdout, 2, out);
					}
				}
				else if (Ctrl->W.use[p[i].level-1]) {
					if (donut_hell) p[i].n = GMT_fix_up_path (&p[i].lon, &p[i].lat, p[i].n, step, 0);
					GMT_n_plot = GMT_geo_to_xy_line (p[i].lon, p[i].lat, p[i].n);
					k = (int)(p[i].level - 1);
					if (k != last_pen_level) {
						GMT_setpen (&Ctrl->W.pen[k]);
						last_pen_level = k;
					}
					if (GMT_n_plot) GMT_plot_line (GMT_x_plot, GMT_y_plot, GMT_pen, GMT_n_plot);
				}
			}

			GMT_free_polygons (p, np);
			GMT_free ((void *)p);
		}

		GMT_free_shore (&c);

	}
	if (need_coast_base) GMT_shore_cleanup (&c);

	if (clipping) ps_clipon ((double *)NULL, (double *)NULL, 0, GMT_no_rgb, 2);	/* Dummy to end path */

	if (gmtdefs.verbose) fprintf (stderr, "\n");

	if (Ctrl->I.active) {	/* Read rivers file and plot as lines */

		if (gmtdefs.verbose) fprintf (stderr, "%s: Adding Rivers...", GMT_program);
		if (!dumping) ps_comment ("Start of River segments");
		last_k = -1;

		for (ind = 0; ind < r.nb; ind++) {	/* Loop over necessary bins only */

			bin = (int)r.bins[ind];
			GMT_get_br_bin (ind, &r, Ctrl->I.use, n_rlevels);

			if (r.ns == 0) continue;

			if (GMT_world_map && greenwich) {
				left = r.bsize * (bin % r.bin_nx);	right = left + r.bsize;
				shift = ((left - edge) <= 180.0 && (right - edge) > 180.0);
			}

			if ((np = (int)GMT_assemble_br (&r, shift, edge, &p)) == 0) {
				GMT_free_br (&r);
				continue;
			}

			for (i = 0; i < np; i++) {
				if (dumping) {
					sprintf (GMT_io.segment_header, "%c River Bin # %d, Level %ld\n", GMT_io.EOF_flag[GMT_OUT], bin, p[i].level);
					GMT_write_segmentheader (GMT_stdout, 2);
					for (k = 0; k < p[i].n; k++) {
						out[GMT_X] = p[i].lon[k];
						out[GMT_Y] = p[i].lat[k];
						GMT_output (GMT_stdout, 2, out);
					}
				}
				else {
					GMT_n_plot = GMT_geo_to_xy_line (p[i].lon, p[i].lat, p[i].n);
					if (!GMT_n_plot) continue;

					k = (int)p[i].level;
					if (k != last_k) {
						GMT_setpen (&Ctrl->I.pen[k]);
						last_k = k;
					}
					GMT_plot_line (GMT_x_plot, GMT_y_plot, GMT_pen, GMT_n_plot);
				}
			}

			/* Free up memory */

			GMT_free_br (&r);
			GMT_free_polygons (p, np);
			GMT_free ((void *)p);
		}
		GMT_br_cleanup (&r);

		if (gmtdefs.verbose) fprintf (stderr, "\n");
	}

	if (Ctrl->N.active) {	/* Read borders file and plot as lines */

		if (gmtdefs.verbose) fprintf (stderr, "%s: Adding Borders...", GMT_program);
		if (!dumping) ps_comment ("Start of Border segments");

		/* Must resample borders because some points may be too far apart and look like 'jumps' */

		step = MAX (fabs(project_info.w - project_info.e), fabs (project_info.n - project_info.s)) / 4.0;

		last_k = -1;

		for (ind = 0; ind < b.nb; ind++) {	/* Loop over necessary bins only */

			bin = (int)b.bins[ind];
			GMT_get_br_bin (ind, &b, Ctrl->N.use, n_blevels);

			if (b.ns == 0) continue;

			if (GMT_world_map && greenwich) {
				left = b.bsize * (bin % b.bin_nx);	right = left + b.bsize;
				shift = ((left - edge) <= 180.0 && (right - edge) > 180.0);
			}

			if ((np = (int)GMT_assemble_br (&b, shift, edge, &p)) == 0) {
				GMT_free_br (&b);
				continue;
			}

			for (i = 0; i < np; i++) {
				if (dumping) {
					sprintf (GMT_io.segment_header, "%c Border Bin # %d, Level %ld\n", GMT_io.EOF_flag[GMT_OUT], bin, p[i].level);
					GMT_write_segmentheader (GMT_stdout, 2);
					for (k = 0; k < p[i].n; k++) {
						out[GMT_X] = p[i].lon[k];
						out[GMT_Y] = p[i].lat[k];
						GMT_output (GMT_stdout, 2, out);
					}
				}
				else {
					p[i].n = GMT_fix_up_path (&p[i].lon, &p[i].lat, p[i].n, step, 0);
					GMT_n_plot = GMT_geo_to_xy_line (p[i].lon, p[i].lat, p[i].n);
					if (!GMT_n_plot) continue;

					k = (int)(p[i].level - 1);
					if (k != last_k) {
						GMT_setpen (&Ctrl->N.pen[k]);
						last_k = k;
					}
					GMT_plot_line (GMT_x_plot, GMT_y_plot, GMT_pen, GMT_n_plot);
				}
			}

			/* Free up memory */

			GMT_free_br (&b);
			GMT_free_polygons (p, np);
			GMT_free ((void *)p);
		}
		GMT_br_cleanup (&b);

		if (gmtdefs.verbose) fprintf (stderr, "\n");
	}

	if (!dumping) {
		if (frame_info.plot) {
			/* project_info.w = west_border;
			project_info.e = east_border; */
			GMT_world_map = world_map_save;
			GMT_map_basemap ();
		}

		if (Ctrl->L.active) GMT_draw_map_scale (&Ctrl->L.item);
		if (Ctrl->T.active) GMT_draw_map_rose (&Ctrl->T.item);

		if (project_info.three_D) ps_rotatetrans (z_project.xmin, z_project.ymin, 0.0);
		GMT_plotend ();
	}

	if (gmtdefs.verbose) fprintf (stderr, "Done\n");
	
	Free_pscoast_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void recursive_path (int k0, int np, struct GMT_GSHHS_POL p[], int level, struct GMT_FILL fill[]) {
	/* To avoid pointing where no paint should go we assemble all the swiss cheese polygons by starting
	 * with the lowest level (0) and looking for polygons inside polygons of the same level.  We do this
	 * using a recursive algorithm */

	int k;
	char text[GMT_TEXT_LEN];
	GMT_LONG add_this_polygon_to_path (int k0, struct GMT_GSHHS_POL *p, int level, int k);
	if (level > GMT_MAX_GSHHS_LEVEL) return;
	for (k = k0 + 1; k < np; k++) {
		if (p[k].n == 0 || p[k].level < level) continue;
		if (add_this_polygon_to_path (k0, p, level, k)) {	/* Add this to the current path */
			sprintf (text, "Polygon %d", k);
			ps_comment (text);
			ps_line (p[k].lon, p[k].lat, p[k].n, 1, FALSE);
#ifdef DEBUGX
			fprintf (stderr, "#\n");
			for (i = 0; i < p[k].n; i++) fprintf (stderr, "%g\t%g\n", p[k].lon[i], p[k].lat[i]);
#endif
			recursive_path (k, np, p, level+1, fill);
			p[k].n = 0;	/* Mark as used */
		}

		if (k0 == -1 && fill) {	/* Done nesting, time to paint the assembled swiss cheese polygon */
			GMT_setfill (&fill[p[k].fid], 0);
			ps_command ("FO");
		}
	}
}

GMT_LONG add_this_polygon_to_path (int k0, struct GMT_GSHHS_POL *p, int level, int k)
{
	/* Determines if we should add the current polygon pol[k] to the growing path we are constructing */

	GMT_LONG j;
	if (p[k].level != level) return (FALSE);	/* Sorry, this polygon does not have the correct level */
	if (k0 == -1) return (TRUE);			/* The first one is always used */
	/* Must make sure the polygon is inside the mother polygon.  Because numerical noise can trick us we try up to two points */
	if (GMT_non_zero_winding (p[k].lon[0], p[k].lat[0], p[k0].lon, p[k0].lat, p[k0].n)) return (TRUE);	/* OK, we're in! */
	/* First point was not inside.  Test another point to see if the first failure was just an accident */
	j = p[k].n / 2;	/* We pick the second point from the other half of the polygon */
	if (GMT_non_zero_winding (p[k].lon[j], p[k].lat[j], p[k0].lon, p[k0].lat, p[k0].n)) return (TRUE);	/* One of those rare occasions when we almost missed a polygon */
	return (FALSE);	/* No, we do not want to use this polygon */
}

void *New_pscoast_Ctrl () {	/* Allocate and initialize a new control structure */
	int k;
	struct PSCOAST_CTRL *C;

	C = (struct PSCOAST_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSCOAST_CTRL), "New_pscoast_Ctrl");

	/* Initialize values whose defaults are not 0/FALSE/NULL */

	C->A.info.high = GMT_MAX_GSHHS_LEVEL;			/* Include all GSHHS levels */
	C->D.set = 'l';						/* Low-resolution coastline data */
	C->E.azimuth = 180.0;
	C->E.elevation = 90.0;
	GMT_init_fill (&C->S.fill, 255, 255, 255);		/* Default Ocean color */
	GMT_init_fill (&C->G.fill, 0, 0, 0);			/* Default Land color */
	GMT_init_fill (&C->C.fill[LAKE], 255, 255, 255);	/* Default Lake color */
	GMT_init_fill (&C->C.fill[RIVER], 255, 255, 255);	/* Default Riverlake color */
	for (k = 0; k < GMT_N_RLEVELS; k++) GMT_init_pen (&C->I.pen[k], GMT_PENWIDTH);		/* Default river pens */
	for (k = 0; k < GMT_N_BLEVELS; k++) GMT_init_pen (&C->N.pen[k], GMT_PENWIDTH);		/* Default border pens */
	for (k = 0; k < GMT_MAX_GSHHS_LEVEL; k++) GMT_init_pen (&C->W.pen[k], GMT_PENWIDTH);	/* Default coastline pens */
	memset ((void *)&C->L.item, 0, sizeof (struct GMT_MAP_SCALE));
	memset ((void *)&C->T.item, 0, sizeof (struct GMT_MAP_ROSE));

	return ((void *)C);
}

void Free_pscoast_Ctrl (struct PSCOAST_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);
}
