/*--------------------------------------------------------------------
 *	$Id: grdinfo.c 17560 2017-02-17 22:05:42Z pwessel $
 *
 *	Copyright (c) 1991-2017 by P. Wessel, W. H. F. Smith, R. Scharroo, J. Luis and F. Wobbe
 *	See LICENSE.TXT file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU Lesser General Public License as published by
 *	the Free Software Foundation; version 3 or any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU Lesser General Public License for more details.
 *
 *	Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/
/*
 * Author:	Paul Wessel
 * Date:	1-JAN-2010
 * Version:	5 API
 *
 * Brief synopsis: grdinfo reads one or more grid file and [optionally] prints
 * out various statistics like mean/standard deviation and median/scale.
 *
 */

#define THIS_MODULE_NAME	"grdinfo"
#define THIS_MODULE_LIB		"core"
#define THIS_MODULE_PURPOSE	"Extract information from grids"
#define THIS_MODULE_KEYS	"<G{+,>T},>DC"

#include "gmt_dev.h"

#define GMT_PROG_OPTIONS "->RVfh"

/* Control structure for grdinfo */

enum Opt_I_modes {
	GRDINFO_GIVE_INCREMENTS = 0,
	GRDINFO_GIVE_REG_ORIG,
	GRDINFO_GIVE_REG_ROUNDED,
	GRDINFO_GIVE_BOUNDBOX};

struct GRDINFO_CTRL {
	struct GRDINFO_C {	/* -C */
		bool active;
	} C;
	struct GRDINFO_F {	/* -F */
		bool active;
	} F;
	struct GRDINFO_I {	/* -Idx[/dy] */
		bool active;
		unsigned int status;
		double inc[2];
	} I;
	struct GRDINFO_M {	/* -M */
		bool active;
	} M;
	struct GRDINFO_L {	/* -L[1|2] */
		bool active;
		unsigned int norm;
	} L;
	struct GRDINFO_T {	/* -T[s]<dz>  -T[<dz>][+s][+a[<alpha>]] */
		bool active;
		unsigned int mode;
		double inc;
		double alpha;
	} T;
	struct GRDINFO_G {	/*  */
		bool active;
		char *opts;
	} G;
};

GMT_LOCAL void *New_Ctrl (struct GMT_CTRL *GMT) {	/* Allocate and initialize a new control structure */
	struct GRDINFO_CTRL *C;

	C = gmt_M_memory (GMT, NULL, 1, struct GRDINFO_CTRL);

	/* Initialize values whose defaults are not 0/false/NULL */
	C->T.alpha = 2.0;	/* 2 % alpha trim is default if selected */
	return (C);
}

GMT_LOCAL void Free_Ctrl (struct GMT_CTRL *GMT, struct GRDINFO_CTRL *C) {	/* Deallocate control structure */
	if (!C) return;
	gmt_M_free (GMT, C);
}

GMT_LOCAL int usage (struct GMTAPI_CTRL *API, int level) {
	gmt_show_name_and_purpose (API, THIS_MODULE_LIB, THIS_MODULE_NAME, THIS_MODULE_PURPOSE);
	if (level == GMT_MODULE_PURPOSE) return (GMT_NOERROR);
	GMT_Message (API, GMT_TIME_NONE, "usage: grdinfo <grid> [-C] [-F] [-I[<dx>[/<dy>]|r|b]] [-L[0|1|2]] [-M]\n");
	GMT_Message (API, GMT_TIME_NONE, "	[%s] [-T[<dz>][+a[<alpha>]][+s]] [%s] [%s]\n\t[%s]\n\n", GMT_Rgeo_OPT, GMT_V_OPT, GMT_f_OPT, GMT_ho_OPT);

	if (level == GMT_SYNOPSIS) return (GMT_MODULE_SYNOPSIS);

	GMT_Message (API, GMT_TIME_NONE, "\t<grid> may be one or more grid files.\n");
	GMT_Message (API, GMT_TIME_NONE, "\n\tOPTIONS:\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-C Format report in fields on a single line using the format\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   file w e s n z0 z1 dx dy n_columns n_rows [x0 y0 x1 y1] [med scale] [mean std rms] [n_nan].\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   (-M gives [x0 y0 x1 y1] and [n_nan]; -L1 gives [med scale]; -L2 gives [mean std rms]).\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-F Report domain in world mapping format [Default is generic].\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-I Return textstring -Rw/e/s/n to nearest multiple of dx/dy.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   If -C is set then rounding off will occur but no -R string is issued.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   If no argument is given then the -I<xinc>/<yinc> string is issued.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   If -Ir is given then the grid's -R string is issued.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   If -Ib is given then the grid's bounding box polygon is issued.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-L Set report mode:\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   -L0 reports range of data by actually reading them (not from header).\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   -L1 reports median and L1-scale of data set.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   -L[2] reports mean, standard deviation, and rms of data set.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-M Search for the global min and max locations (x0,y0) and (x1,y1).\n");
	GMT_Option (API, "R");
	GMT_Message (API, GMT_TIME_NONE, "\t-T Print global -Tzmin/zmax[/dz] (in rounded multiples of dz, if given).\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   Append +a[<alpha>] to trim grid range by excluding the two <alpha>/2 tails [2 %%].\n");
	GMT_Message (API, GMT_TIME_NONE, "\t     Note: +a is limited to a single grid.  Give <alpha> in percent.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   Append +s to force a symmetrical range about zero.\n");
	GMT_Option (API, "V,f,h,.");
	
	return (GMT_MODULE_USAGE);
}

GMT_LOCAL int parse (struct GMT_CTRL *GMT, struct GRDINFO_CTRL *Ctrl, struct GMT_OPTION *options) {

	/* This parses the options provided to grdcut and sets parameters in CTRL.
	 * Any GMT common options will override values set previously by other commands.
	 * It also replaces any file names specified as input or output with the data ID
	 * returned when registering these sources/destinations with the API.
	 */

	unsigned int n_errors = 0, n_files = 0;
char text[GMT_LEN32] = {""};
	struct GMT_OPTION *opt = NULL;

	for (opt = options; opt; opt = opt->next) {	/* Process all the options given */

		switch (opt->option) {
			/* Common parameters */

			case '<':	/* Input files */
				if (gmt_check_filearg (GMT, '<', opt->arg, GMT_IN, GMT_IS_GRID))
					n_files++;
				else
					n_errors++;
				break;

			/* Processes program-specific parameters */

			case 'C':	/* Column format */
				Ctrl->C.active = true;
				break;
			case 'F':	/* World mapping format */
				Ctrl->F.active = true;
				break;
			case 'G':	/* List of GDAL options */
				Ctrl->G.active = true;
				Ctrl->G.opts = strdup(opt->arg);
				break;
			case 'I':	/* Increment rounding */
				Ctrl->I.active = true;
				if (!opt->arg[0])	/* No args given, we want to output the -I string */
					Ctrl->I.status = GRDINFO_GIVE_INCREMENTS;
				else if ((opt->arg[0] == 'r' || opt->arg[0] == '-') && opt->arg[1] == '\0')	/* -Ir: we want to output the actual -R string */
					Ctrl->I.status = GRDINFO_GIVE_REG_ORIG;
				else if (opt->arg[0] == 'b' && opt->arg[1] == '\0')	/* -Ib means return grid perimeter as bounding box */
					Ctrl->I.status = GRDINFO_GIVE_BOUNDBOX;
				else {	/* Report -R to nearest given multiple increment */
					Ctrl->I.status = GRDINFO_GIVE_REG_ROUNDED;
					if (gmt_getinc (GMT, opt->arg, Ctrl->I.inc)) {
						gmt_inc_syntax (GMT, 'I', 1);
						n_errors++;
					}
				}
				break;
			case 'L':	/* Selects norm */
				Ctrl->L.active = true;
				switch (opt->arg[0]) {
					case '\0': case '2':
						Ctrl->L.norm |= 2; break;
					case '1':
						Ctrl->L.norm |= 1; break;
				}
				break;
			case 'M':	/* Global extrema */
				Ctrl->M.active = true;
				break;
			case 'T':	/* CPT range */
				Ctrl->T.active = true;
				if (opt->arg[0] == 's' && gmt_M_compat_check (GMT, 5)) {	/* Old-style format, cast in new syntax */
					GMT_Report (GMT->parent, GMT_MSG_COMPAT, "Warning: -Ts option is deprecated; please use -T[<dz>][+s][+a[<alpha>]] next time.\n");
					sprintf (text, "%s+s", &opt->arg[1]);
				}
				else
					strncpy (text, opt->arg, GMT_LEN32-1);
				if (gmt_validate_modifiers (GMT, text, opt->option, "as")) {
					GMT_Report (GMT->parent, GMT_MSG_COMPAT, "Error -T: Syntax is -T[<dz>][+s][+a[<alpha>]] next time.\n");
					n_errors++;
				}
				else {
					char string[GMT_LEN32] = {""};
					if (text[0] && text[0] != '+')
						Ctrl->T.inc = atof (text);
					if (gmt_get_modifier (text, 's', string))	/* Want symmetrical range about 0, i.e., -3500/3500[/500] */
						Ctrl->T.mode |= 1;
					if (gmt_get_modifier (text, 'a', string)) {	/* Want alpha-trimmed range before determining limits */
						Ctrl->T.mode |= 2;
						if (string[0]) Ctrl->T.alpha = atof (string);
					}
				}
				break;

			default:	/* Report bad options */
				n_errors += gmt_default_error (GMT, opt->option);
				break;
		}
	}

	n_errors += gmt_M_check_condition (GMT, n_files == 0, "Syntax error: Must specify one or more input files\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->T.active && Ctrl->T.inc < 0.0, "Syntax error -T: The optional increment must be positive\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->T.mode & 2 && n_files != 1, "Syntax error -T: The optional alpha-trim value can only work with a single grid file\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->T.active && (Ctrl->T.alpha < 0.0 || Ctrl->T.alpha > 100.0), "Syntax error -T: The optional alpha-trim value must be in the 0 < alpha < 100 %% range\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->I.active && Ctrl->I.status == GRDINFO_GIVE_REG_ROUNDED && (Ctrl->I.inc[GMT_X] <= 0.0 || Ctrl->I.inc[GMT_Y] <= 0.0), "Syntax error -I: Must specify a positive increment(s)\n");
	n_errors += gmt_M_check_condition (GMT, (Ctrl->I.active || Ctrl->T.active) && Ctrl->M.active, "Syntax error -M: Not compatible with -I or -T\n");
	n_errors += gmt_M_check_condition (GMT, (Ctrl->I.active || Ctrl->T.active) && Ctrl->L.active, "Syntax error -L: Not compatible with -I or -T\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->T.active && Ctrl->I.active, "Syntax error: Only one of -I -T can be specified\n");

	return (n_errors ? GMT_PARSE_ERROR : GMT_NOERROR);
}

#define bailout(code) {gmt_M_free_options (mode); return (code);}
#define Return(code) {Free_Ctrl (GMT, Ctrl); gmt_end_module (GMT, GMT_cpy); bailout (code);}

#if defined(HAVE_GDAL) && (GDAL_VERSION_MAJOR >= 2) && (GDAL_VERSION_MINOR >= 1)
#include "gdal_utils.h"
#include "gmt_gdal_librarified.c"
#endif

int GMT_grdinfo (void *V_API, int mode, void *args) {
	int error = 0;
	unsigned int n_grds = 0, o_type = GMT_IS_TEXTSET, n_cols = 0, col;
	bool subset;

	uint64_t ij, n_nan = 0, n = 0;

	double x_min = 0.0, y_min = 0.0, z_min = 0.0, x_max = 0.0, y_max = 0.0, z_max = 0.0, wesn[4];
	double global_xmin, global_xmax, global_ymin, global_ymax, global_zmin, global_zmax;
	double mean = 0.0, median = 0.0, sum2 = 0.0, stdev = 0.0, scale = 0.0, rms = 0.0, x, out[20];

	char format[GMT_BUFSIZ] = {""}, text[GMT_LEN64] = {""}, record[GMT_BUFSIZ] = {""}, grdfile[GMT_LEN256] = {""};
	char *type[2] = { "Gridline", "Pixel"}, *sep = NULL, *projStr = NULL;

	struct GRDINFO_CTRL *Ctrl = NULL;
	struct GMT_GRID *G = NULL;
	struct GMT_OPTION *opt = NULL;
	struct GMT_CTRL *GMT = NULL, *GMT_cpy = NULL;
	struct GMT_OPTION *options = NULL;
	struct GMTAPI_CTRL *API = gmt_get_api_ptr (V_API);	/* Cast from void to GMTAPI_CTRL pointer */

	/*----------------------- Standard module initialization and parsing ----------------------*/

	if (API == NULL) return (GMT_NOT_A_SESSION);
	if (mode == GMT_MODULE_PURPOSE) return (usage (API, GMT_MODULE_PURPOSE));	/* Return the purpose of program */
	options = GMT_Create_Options (API, mode, args);	if (API->error) return (API->error);	/* Set or get option list */

	if (!options || options->option == GMT_OPT_USAGE) bailout (usage (API, GMT_USAGE));	/* Return the usage message */
	if (options->option == GMT_OPT_SYNOPSIS) bailout (usage (API, GMT_SYNOPSIS));	/* Return the synopsis */

	/* Parse the command-line arguments */

	GMT = gmt_begin_module (API, THIS_MODULE_LIB, THIS_MODULE_NAME, &GMT_cpy); /* Save current state */
	if (GMT_Parse_Common (API, GMT_PROG_OPTIONS, options)) Return (API->error);
	Ctrl = New_Ctrl (GMT);	/* Allocate and initialize a new control structure */
	if ((error = parse (GMT, Ctrl, options)) != 0) Return (error);

	/*---------------------------- This is the grdinfo main code ----------------------------*/

	/* OK, done parsing, now process all input grids in a loop */
	
	sep = GMT->current.setting.io_col_separator;
	gmt_M_memcpy (wesn, GMT->common.R.wesn, 4, double);	/* Current -R setting, if any */
	global_xmin = global_ymin = global_zmin = DBL_MAX;
	global_xmax = global_ymax = global_zmax = -DBL_MAX;
	if (Ctrl->C.active) {
		if (API->mode) o_type = GMT_IS_DATASET;	/* With external interface we are returning doubles */
		n_cols = 6;	/* w e s n z0 z1 */
		if (!Ctrl->I.active) {
			n_cols += 4;				/* Add dx dy n_columns n_rows */
			if (Ctrl->M.active) n_cols += 5;	/* Add x0 y0 x1 y1 nnan */
			if (Ctrl->L.norm & 1) n_cols += 2;	/* Add median scale */
			if (Ctrl->L.norm & 2) n_cols += 3;	/* Add mean stdev rms */
		}
	}
	if (GMT_Init_IO (API, o_type, GMT_IS_NONE, GMT_OUT, GMT_ADD_DEFAULT, 0, options) != GMT_NOERROR) {	/* Registers default output destination, unless already set */
		Return (API->error);
	}
	if (GMT_Begin_IO (API, o_type, GMT_OUT, GMT_HEADER_OFF) != GMT_NOERROR) {	/* Enables data output and sets access mode */
		Return (API->error);
	}
	if (GMT_Set_Geometry (API, GMT_OUT, GMT_IS_NONE) != GMT_NOERROR) {	/* Sets output geometry */
		Return (API->error);
	}
	if (n_cols && (error = gmt_set_cols (GMT, GMT_OUT, n_cols)) != 0) Return (error);	/* Set number of output columns */
	
	for (opt = options; opt; opt = opt->next) {	/* Loop over arguments, skip options */ 

		if (opt->option != '<') continue;	/* We are only processing filenames here */

#if defined(HAVE_GDAL) && (GDAL_VERSION_MAJOR >= 2) && (GDAL_VERSION_MINOR >= 1)
	if (Ctrl->G.active)
		grid_gdal_librarified (GMT, opt->arg, Ctrl->G.opts);
#endif

		gmt_set_cartesian (GMT, GMT_IN);	/* Reset since we may get a bunch of files, some geo, some not */

		if ((G = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, opt->arg, NULL)) == NULL) {
			Return (API->error);
		}
		subset = gmt_M_is_subset (GMT, G->header, wesn);	/* Subset requested */
		if (subset) gmt_M_err_fail (GMT, gmt_adjust_loose_wesn (GMT, wesn, G->header), "");	/* Make sure wesn matches header spacing */

		GMT_Report (API, GMT_MSG_VERBOSE, "Processing grid %s\n", G->header->name);

		if (G->header->ProjRefPROJ4 && !Ctrl->C.active && !Ctrl->T.active)
			projStr = strdup(G->header->ProjRefPROJ4);		/* Copy proj string to print at the end */
		else if (G->header->ProjRefWKT && !Ctrl->C.active && !Ctrl->T.active) 
			projStr = strdup(G->header->ProjRefWKT);

		for (n = 0; n < GMT_Z; n++) GMT->current.io.col_type[GMT_OUT][n] = GMT->current.io.col_type[GMT_IN][n];	/* Since grids may differ in types */

		n_grds++;

		if (Ctrl->M.active || Ctrl->L.active || subset || (Ctrl->T.mode & 2)) {	/* Need to read the data (all or subset) */
			if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, wesn, opt->arg, G) == NULL) {
				Return (API->error);
			}
		}

		if (Ctrl->T.mode & 2) strncpy (grdfile, opt->arg, GMT_LEN256-1);
		
		if (Ctrl->M.active || Ctrl->L.active) {	/* Must determine the location of global min and max values */
			uint64_t ij_min, ij_max;
			unsigned int col, row;

			z_min = DBL_MAX;	z_max = -DBL_MAX;
			mean = median = sum2 = 0.0;
			ij_min = ij_max = n = 0;
			gmt_M_grd_loop (GMT, G, row, col, ij) {
				if (gmt_M_is_fnan (G->data[ij])) continue;
				if (G->data[ij] < z_min) {
					z_min = G->data[ij];	ij_min = ij;
				}
				if (G->data[ij] > z_max) {
					z_max = G->data[ij];	ij_max = ij;
				}
				n++;
				if (Ctrl->L.active) {	/* Use Welford (1962) algorithm to compute mean and corrected sum of squares */
					x = G->data[ij] - mean;
					mean += x / n;
					sum2 += x * (G->data[ij] - mean);
				}
			}

			n_nan = G->header->nm - n;
			if (n) {	/* Meaning at least one non-NaN node was found */
				col = (unsigned int)gmt_M_col (G->header, ij_min);
				row = (unsigned int)gmt_M_row (G->header, ij_min);
				x_min = gmt_M_grd_col_to_x (GMT, col, G->header);
				y_min = gmt_M_grd_row_to_y (GMT, row, G->header);
				col = (unsigned int)gmt_M_col (G->header, ij_max);
				row = (unsigned int)gmt_M_row (G->header, ij_max);
				x_max = gmt_M_grd_col_to_x (GMT, col, G->header);
				y_max = gmt_M_grd_row_to_y (GMT, row, G->header);
			}
			else	/* Not a single valid node */
				x_min = x_max = y_min = y_max = GMT->session.d_NaN;
		}

		if (Ctrl->L.norm & 1) {	/* Calculate the median and L1 scale */
			int new_grid;
			struct GMT_GRID *G2 = NULL;

			/* Note that this option rearranges the input grid, so if a memory location is passed then
			 * the grid in the calling program is no longer the original values */
			new_grid = gmt_set_outgrid (GMT, opt->arg, false, G, &G2);	/* true if input is a read-only array */
			gmt_grd_pad_off (GMT, G2);	/* Undo pad if one existed */
			gmt_sort_array (GMT, G2->data, G2->header->nm, GMT_FLOAT);
			median = (n%2) ? G2->data[n/2] : 0.5*(G2->data[n/2-1] + G2->data[n/2]);
			for (ij = 0; ij < n; ij++) G2->data[ij] = (float)fabs (G2->data[ij] - median);
			gmt_sort_array (GMT, G2->data, n, GMT_FLOAT);
			scale = (n%2) ? 1.4826 * G2->data[n/2] : 0.7413 * (G2->data[n/2-1] + G2->data[n/2]);
			if (new_grid) {	/* Free the temporary grid */
				if (GMT_Destroy_Data (API, &G2) != GMT_NOERROR) {
					GMT_Report (API, GMT_MSG_NORMAL, "Failed to free G2\n");
				}
			}
		}
		if (Ctrl->L.norm & 2) {	/* Calculate the mean, standard deviation, and rms */
			x = (double)n;
			stdev = (n > 1) ? sqrt (sum2 / (x-1)) : GMT->session.d_NaN;
			rms = (n > 0) ? sqrt (sum2 / x + mean * mean) : GMT->session.d_NaN;
			mean = (n > 0) ? mean : GMT->session.d_NaN;
		}

		if (gmt_M_is_geographic (GMT, GMT_IN)) {
			if (gmt_M_grd_is_global(GMT, G->header) || (G->header->wesn[XLO] < 0.0 && G->header->wesn[XHI] <= 0.0))
				GMT->current.io.geo.range = GMT_IS_GIVEN_RANGE;
			else if (G->header->wesn[XLO] < 0.0 && G->header->wesn[XHI] >= 0.0)
				GMT->current.io.geo.range = GMT_IS_M180_TO_P180_RANGE;
			else
				GMT->current.io.geo.range = GMT_IS_0_TO_P360_RANGE;
		}

		/* OK, time to report results */

		if (Ctrl->I.active && Ctrl->I.status == GRDINFO_GIVE_REG_ORIG) {
			sprintf (record, "-R");
			gmt_ascii_format_col (GMT, text, G->header->wesn[XLO], GMT_OUT, GMT_X);	strcat (record, text);	strcat (record, "/");
			gmt_ascii_format_col (GMT, text, G->header->wesn[XHI], GMT_OUT, GMT_X);	strcat (record, text);	strcat (record, "/");
			gmt_ascii_format_col (GMT, text, G->header->wesn[YLO], GMT_OUT, GMT_Y);	strcat (record, text);	strcat (record, "/");
			gmt_ascii_format_col (GMT, text, G->header->wesn[YHI], GMT_OUT, GMT_Y);	strcat (record, text);
			GMT_Put_Record (API, GMT_WRITE_TEXT, record);
		} else if (Ctrl->I.active && Ctrl->I.status == GRDINFO_GIVE_INCREMENTS) {
			sprintf (record, "-I");
			gmt_ascii_format_col (GMT, text, G->header->inc[GMT_X], GMT_OUT, GMT_Z);	strcat (record, text);	strcat (record, "/");
			gmt_ascii_format_col (GMT, text, G->header->inc[GMT_Y], GMT_OUT, GMT_Z);	strcat (record, text);
			GMT_Put_Record (API, GMT_WRITE_TEXT, record);
		} else if (Ctrl->I.active && Ctrl->I.status == GRDINFO_GIVE_BOUNDBOX) {
			if (GMT_Set_Geometry (API, GMT_OUT, GMT_IS_POLY) != GMT_NOERROR) {	/* Sets output geometry */
				Return (API->error);
			}
			sprintf (record, "> Bounding box for %s", G->header->name);
			GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			/* LL */
			gmt_ascii_format_col (GMT, record, G->header->wesn[XLO], GMT_OUT, GMT_X);	strcat (record, sep);
			gmt_ascii_format_col (GMT, text,   G->header->wesn[YLO], GMT_OUT, GMT_Y);	strcat (record, text);
			GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			/* LR */
			gmt_ascii_format_col (GMT, record, G->header->wesn[XHI], GMT_OUT, GMT_X);	strcat (record, sep);
			gmt_ascii_format_col (GMT, text,   G->header->wesn[YLO], GMT_OUT, GMT_Y);	strcat (record, text);
			GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			/* UR */
			gmt_ascii_format_col (GMT, record, G->header->wesn[XHI], GMT_OUT, GMT_X);	strcat (record, sep);
			gmt_ascii_format_col (GMT, text,   G->header->wesn[YHI], GMT_OUT, GMT_Y);	strcat (record, text);
			GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			/* UL */
			gmt_ascii_format_col (GMT, record, G->header->wesn[XLO], GMT_OUT, GMT_X);	strcat (record, sep);
			gmt_ascii_format_col (GMT, text,   G->header->wesn[YHI], GMT_OUT, GMT_Y);	strcat (record, text);
			GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			/* LL (repeat to close polygon) */
			gmt_ascii_format_col (GMT, record, G->header->wesn[XLO], GMT_OUT, GMT_X);	strcat (record, sep);
			gmt_ascii_format_col (GMT, text,   G->header->wesn[YLO], GMT_OUT, GMT_Y);	strcat (record, text);
			GMT_Put_Record (API, GMT_WRITE_TEXT, record);
		} else if (Ctrl->C.active && !Ctrl->I.active) {
			if (API->mode) {	/* External interface, return as data with no leading text */
				/* w e s n z0 z1 dx dy n_columns n_rows [x0 y0 x1 y1] [med scale] [mean std rms] [n_nan] */
				gmt_M_memcpy (out, G->header->wesn, 4, double);	/* Place the w/e/s/n limits */
				out[ZLO]   = G->header->z_min;		out[ZHI]   = G->header->z_max;
				out[ZHI+1] = G->header->inc[GMT_X];	out[ZHI+2] = G->header->inc[GMT_Y];
				out[ZHI+3] = G->header->n_columns;		out[ZHI+4] = G->header->n_rows;
				col = ZHI+5;
				if (Ctrl->M.active) {
					out[col++] = x_min;	out[col++] = y_min;
					out[col++] = x_max;	out[col++] = y_max;
				}
				if (Ctrl->L.norm & 1) {
					out[col++] = median;	out[col++] = scale;
				}
				if (Ctrl->L.norm & 2) {
					out[col++] = mean;	out[col++] = stdev;	out[col++] = rms;
				}
				if (Ctrl->M.active) {
					out[col++] = (double)n_nan;
				}
				GMT_Put_Record (API, GMT_WRITE_DATA, out);
			}
			else {	/* Command-line usage */
				sprintf (record, "%s%s", G->header->name, sep);
				gmt_ascii_format_col (GMT, text, G->header->wesn[XLO], GMT_OUT, GMT_X);	strcat (record, text);	strcat (record, sep);
				gmt_ascii_format_col (GMT, text, G->header->wesn[XHI], GMT_OUT, GMT_X);	strcat (record, text);	strcat (record, sep);
				gmt_ascii_format_col (GMT, text, G->header->wesn[YLO], GMT_OUT, GMT_Y);	strcat (record, text);	strcat (record, sep);
				gmt_ascii_format_col (GMT, text, G->header->wesn[YHI], GMT_OUT, GMT_Y);	strcat (record, text);	strcat (record, sep);
				gmt_ascii_format_col (GMT, text, G->header->z_min, GMT_OUT, GMT_Z);	strcat (record, text);	strcat (record, sep);
				gmt_ascii_format_col (GMT, text, G->header->z_max, GMT_OUT, GMT_Z);	strcat (record, text);	strcat (record, sep);
				gmt_ascii_format_col (GMT, text, G->header->inc[GMT_X], GMT_OUT, GMT_X);
				if (isalpha ((int)text[strlen(text)-1])) text[strlen(text)-1] = '\0';	/* Chop of trailing WESN flag here */
				strcat (record, text);	strcat (record, sep);
				gmt_ascii_format_col (GMT, text, G->header->inc[GMT_Y], GMT_OUT, GMT_Y);
				if (isalpha ((int)text[strlen(text)-1])) text[strlen(text)-1] = '\0';	/* Chop of trailing WESN flag here */
				strcat (record, text);	strcat (record, sep);
				gmt_ascii_format_col (GMT, text, (double)G->header->n_columns, GMT_OUT, GMT_Z);	strcat (record, text);	strcat (record, sep);
				gmt_ascii_format_col (GMT, text, (double)G->header->n_rows, GMT_OUT, GMT_Z);	strcat (record, text);

				if (Ctrl->M.active) {
					strcat (record, sep);	gmt_ascii_format_col (GMT, text, x_min, GMT_OUT, GMT_X);	strcat (record, text);
					strcat (record, sep);	gmt_ascii_format_col (GMT, text, y_min, GMT_OUT, GMT_Y);	strcat (record, text);
					strcat (record, sep);	gmt_ascii_format_col (GMT, text, x_max, GMT_OUT, GMT_X);	strcat (record, text);
					strcat (record, sep);	gmt_ascii_format_col (GMT, text, y_max, GMT_OUT, GMT_Y);	strcat (record, text);
				}
				if (Ctrl->L.norm & 1) {
					strcat (record, sep);	gmt_ascii_format_col (GMT, text, median, GMT_OUT, GMT_Z);	strcat (record, text);
					strcat (record, sep);	gmt_ascii_format_col (GMT, text,  scale, GMT_OUT, GMT_Z);	strcat (record, text);
				}
				if (Ctrl->L.norm & 2) {
					strcat (record, sep);	gmt_ascii_format_col (GMT, text,  mean, GMT_OUT, GMT_Z);	strcat (record, text);
					strcat (record, sep);	gmt_ascii_format_col (GMT, text, stdev, GMT_OUT, GMT_Z);	strcat (record, text);
					strcat (record, sep);	gmt_ascii_format_col (GMT, text,   rms, GMT_OUT, GMT_Z);	strcat (record, text);
				}
				if (Ctrl->M.active) { sprintf (text, "%s%" PRIu64, sep, n_nan);	strcat (record, text); }
				GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			}
		}
		else if (!(Ctrl->T.active || (Ctrl->I.active && Ctrl->I.status == GRDINFO_GIVE_REG_ROUNDED))) {
			char *gtype[2] = {"Cartesian grid", "Geographic grid"};
			sprintf (record, "%s: Title: %s", G->header->name, G->header->title);		GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			sprintf (record, "%s: Command: %s", G->header->name, G->header->command);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			sprintf (record, "%s: Remark: %s", G->header->name, G->header->remark);		GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			if (G->header->registration == GMT_GRID_NODE_REG || G->header->registration == GMT_GRID_PIXEL_REG)
				sprintf (record, "%s: %s node registration used [%s]", G->header->name, type[G->header->registration], gtype[gmt_M_is_geographic (GMT, GMT_IN)]);
			else
				sprintf (record, "%s: Unknown registration! Probably not a GMT grid", G->header->name);
			GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			if (G->header->type != k_grd_unknown_fmt)
				sprintf (record, "%s: Grid file format: %s", G->header->name, GMT->session.grdformat[G->header->type]);
			else
				sprintf (record, "%s: Unrecognized grid file format! Probably not a GMT grid", G->header->name);
			GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			if (Ctrl->F.active) {
				if ((fabs (G->header->wesn[XLO]) < 500.0) && (fabs (G->header->wesn[XHI]) < 500.0) && (fabs (G->header->wesn[YLO]) < 500.0) && (fabs (G->header->wesn[YHI]) < 500.0)) {
					sprintf (record, "%s: x_min: %.7f", G->header->name, G->header->wesn[XLO]);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: x_max: %.7f", G->header->name, G->header->wesn[XHI]);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: x_inc: %.7f", G->header->name, G->header->inc[GMT_X]);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: name: %s",    G->header->name, G->header->x_units);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: n_columns: %d",      G->header->name, G->header->n_columns);		GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: y_min: %.7f", G->header->name, G->header->wesn[YLO]);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: y_max: %.7f", G->header->name, G->header->wesn[YHI]);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: y_inc: %.7f", G->header->name, G->header->inc[GMT_Y]);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: name: %s",    G->header->name, G->header->y_units);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: n_rows: %d",      G->header->name, G->header->n_rows);		GMT_Put_Record (API, GMT_WRITE_TEXT, record);
				}
				else {
					sprintf (record, "%s: x_min: %.2f", G->header->name, G->header->wesn[XLO]);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: x_max: %.2f", G->header->name, G->header->wesn[XHI]);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: x_inc: %.2f", G->header->name, G->header->inc[GMT_X]);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: name: %s",    G->header->name, G->header->x_units);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: n_columns: %d",      G->header->name, G->header->n_columns);		GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: y_min: %.2f", G->header->name, G->header->wesn[YLO]);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: y_max: %.2f", G->header->name, G->header->wesn[YHI]);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: y_inc: %.2f", G->header->name, G->header->inc[GMT_Y]);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: name: %s",    G->header->name, G->header->y_units);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
					sprintf (record, "%s: n_rows: %d",      G->header->name, G->header->n_rows);		GMT_Put_Record (API, GMT_WRITE_TEXT, record);
				}
			}
			else {
				sprintf (record, "%s: x_min: ", G->header->name);
				gmt_ascii_format_col (GMT, text, G->header->wesn[XLO], GMT_OUT, GMT_X);	strcat (record, text);
				strcat (record, " x_max: ");
				gmt_ascii_format_col (GMT, text, G->header->wesn[XHI], GMT_OUT, GMT_X);	strcat (record, text);
				gmt_ascii_format_col (GMT, text, G->header->inc[GMT_X], GMT_OUT, GMT_X);
				if (isalpha ((int)text[strlen(text)-1])) text[strlen(text)-1] = '\0';	/* Chop of trailing WESN flag here */
				strcat (record, " x_inc: ");	strcat (record, text);
				strcat (record, " name: ");	strcat (record, G->header->x_units);
				sprintf (text, " n_columns: %d", G->header->n_columns);	strcat (record, text);
				GMT_Put_Record (API, GMT_WRITE_TEXT, record);
				sprintf (record, "%s: y_min: ", G->header->name);
				gmt_ascii_format_col (GMT, text, G->header->wesn[YLO], GMT_OUT, GMT_Y);	strcat (record, text);
				strcat (record, " y_max: ");
				gmt_ascii_format_col (GMT, text, G->header->wesn[YHI], GMT_OUT, GMT_Y);	strcat (record, text);
				gmt_ascii_format_col (GMT, text, G->header->inc[GMT_Y], GMT_OUT, GMT_Y);
				if (isalpha ((int)text[strlen(text)-1])) text[strlen(text)-1] = '\0';	/* Chop of trailing WESN flag here */
				strcat (record, " y_inc: ");	strcat (record, text);
				strcat (record, " name: ");	strcat (record, G->header->y_units);
				sprintf (text, " n_rows: %d", G->header->n_rows);	strcat (record, text);
				GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			}

			if (Ctrl->M.active) {
				if (z_min == DBL_MAX) z_min = GMT->session.d_NaN;
				if (z_max == -DBL_MAX) z_max = GMT->session.d_NaN;
				sprintf (record, "%s: z_min: ", G->header->name);
				gmt_ascii_format_col (GMT, text, z_min, GMT_OUT, GMT_Z);	strcat (record, text);
				strcat (record, " at x = ");
				gmt_ascii_format_col (GMT, text, x_min, GMT_OUT, GMT_X);	strcat (record, text);
				strcat (record, " y = ");
				gmt_ascii_format_col (GMT, text, y_min, GMT_OUT, GMT_Y);	strcat (record, text);
				strcat (record, " z_max: ");
				gmt_ascii_format_col (GMT, text, z_max, GMT_OUT, GMT_Z);	strcat (record, text);
				strcat (record, " at x = ");
				gmt_ascii_format_col (GMT, text, x_max, GMT_OUT, GMT_X);	strcat (record, text);
				strcat (record, " y = ");
				gmt_ascii_format_col (GMT, text, y_max, GMT_OUT, GMT_Y);	strcat (record, text);
				GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			}
			else if (Ctrl->F.active) {
				sprintf (record, "%s: zmin: %g", G->header->name, G->header->z_min);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
				sprintf (record, "%s: zmax: %g", G->header->name, G->header->z_max);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
				sprintf (record, "%s: name: %s", G->header->name, G->header->z_units);	GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			}
			else {
				sprintf (record, "%s: z_min: ", G->header->name);
				gmt_ascii_format_col (GMT, text, G->header->z_min, GMT_OUT, GMT_Z);	strcat (record, text);
				strcat (record, " z_max: ");
				gmt_ascii_format_col (GMT, text, G->header->z_max, GMT_OUT, GMT_Z);	strcat (record, text);
				strcat (record, " name: ");	strcat (record, G->header->z_units);
				GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			}

			/* print scale and offset */
			sprintf (format, "%s: scale_factor: %s add_offset: %s",
			         G->header->name, GMT->current.setting.format_float_out, GMT->current.setting.format_float_out);
			sprintf (record, format, G->header->z_scale_factor, G->header->z_add_offset);
			if (G->header->z_scale_factor != 1.0 || G->header->z_add_offset != 0) {
				/* print packed z-range */
				sprintf (format, "%s packed z-range: [%s,%s]", record,
				         GMT->current.setting.format_float_out, GMT->current.setting.format_float_out);
				sprintf (record, format,
				         (G->header->z_min - G->header->z_add_offset) / G->header->z_scale_factor,
				         (G->header->z_max - G->header->z_add_offset) / G->header->z_scale_factor);
			}
			GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			if (n_nan) {
				double percent = 100.0 * n_nan / G->header->nm;
				sprintf (record, "%s: %" PRIu64 " nodes (%.1f%%) set to NaN", G->header->name, n_nan, percent);
				GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			}
			if (Ctrl->L.norm & 1) {
				sprintf (record, "%s: median: ", G->header->name);
				gmt_ascii_format_col (GMT, text, median, GMT_OUT, GMT_Z);	strcat (record, text);
				strcat (record, " scale: ");
				gmt_ascii_format_col (GMT, text, scale, GMT_OUT, GMT_Z);	strcat (record, text);
				GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			}
			if (Ctrl->L.norm & 2) {
				sprintf (record, "%s: mean: ", G->header->name);
				gmt_ascii_format_col (GMT, text, mean, GMT_OUT, GMT_Z);	strcat (record, text);
				strcat (record, " stdev: ");
				gmt_ascii_format_col (GMT, text, stdev, GMT_OUT, GMT_Z);	strcat (record, text);
				strcat (record, " rms: ");
				gmt_ascii_format_col (GMT, text, rms, GMT_OUT, GMT_Z);	strcat (record, text);
				GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			}
			if (strspn(GMT->session.grdformat[G->header->type], "nc") != 0) {
				/* type is netCDF: report chunk size and deflation level */
				if (G->header->is_netcdf4) {
					sprintf (text, " chunk_size: %" PRIuS ",%" PRIuS " shuffle: %s deflation_level: %u",
					         G->header->z_chunksize[0], G->header->z_chunksize[1],
					         G->header->z_shuffle ? "on" : "off", G->header->z_deflate_level);
				}
				else
					text[0] = '\0';
				sprintf (record, "%s: format: %s%s",
						G->header->name, G->header->is_netcdf4 ? "netCDF-4" : "classic", text);
				GMT_Put_Record (API, GMT_WRITE_TEXT, record);
			}
		} /* !(Ctrl->T.active || (Ctrl->I.active && Ctrl->I.status == GRDINFO_GIVE_REG_ROUNDED))) */
		else {
			if (G->header->z_min < global_zmin) global_zmin = G->header->z_min;
			if (G->header->z_max > global_zmax) global_zmax = G->header->z_max;
			if (G->header->wesn[XLO] < global_xmin) global_xmin = G->header->wesn[XLO];
			if (G->header->wesn[XHI] > global_xmax) global_xmax = G->header->wesn[XHI];
			if (G->header->wesn[YLO] < global_ymin) global_ymin = G->header->wesn[YLO];
			if (G->header->wesn[YHI] > global_ymax) global_ymax = G->header->wesn[YHI];
		}
		if ((Ctrl->T.mode & 2) == 0 && GMT_Destroy_Data (API, &G) != GMT_NOERROR) {
			Return (API->error);
		}
	}

	if (global_zmin == DBL_MAX) global_zmin = GMT->session.d_NaN;	/* Never got set */
	if (global_zmax == -DBL_MAX) global_zmax = GMT->session.d_NaN;

	if (Ctrl->C.active && (Ctrl->I.active && Ctrl->I.status == GRDINFO_GIVE_REG_ROUNDED)) {
		global_xmin = floor (global_xmin / Ctrl->I.inc[GMT_X]) * Ctrl->I.inc[GMT_X];
		global_xmax = ceil  (global_xmax / Ctrl->I.inc[GMT_X]) * Ctrl->I.inc[GMT_X];
		global_ymin = floor (global_ymin / Ctrl->I.inc[GMT_Y]) * Ctrl->I.inc[GMT_Y];
		global_ymax = ceil  (global_ymax / Ctrl->I.inc[GMT_Y]) * Ctrl->I.inc[GMT_Y];
		if (gmt_M_is_geographic (GMT, GMT_IN)) {	/* Must make sure we don't get outside valid bounds */
			if (global_ymin < -90.0) {
				global_ymin = -90.0;
				GMT_Report (API, GMT_MSG_VERBOSE, "Warning: Using -I caused south to become < -90.  Reset to -90.\n");
			}
			if (global_ymax > 90.0) {
				global_ymax = 90.0;
				GMT_Report (API, GMT_MSG_VERBOSE, "Warning: Using -I caused north to become > +90.  Reset to +90.\n");
			}
			if (fabs (global_xmax - global_xmin) > 360.0) {
				GMT_Report (API, GMT_MSG_VERBOSE, "Warning: Using -I caused longitude range to exceed 360.  Reset to a range of 360.\n");
				global_xmin = (global_xmin < 0.0) ? -180.0 : 0.0;
				global_xmax = (global_xmin < 0.0) ? +180.0 : 360.0;
			}
		}
		if (API->mode) {	/* External interface, return as data with no leading text */
			/* w e s n z0 z1 */
			out[XLO] = global_xmin;		out[XHI] = global_xmax;
			out[YLO] = global_ymin;		out[YHI] = global_ymax;
			out[ZLO] = global_zmin;		out[ZHI] = global_zmax;
			GMT_Put_Record (API, GMT_WRITE_DATA, out);
		}
		else {	/* Command-line usage */
			sprintf (record, "%d%s", n_grds, sep);
			gmt_ascii_format_col (GMT, text, global_xmin, GMT_OUT, GMT_X);	strcat (record, text);	strcat (record, sep);
			gmt_ascii_format_col (GMT, text, global_xmax, GMT_OUT, GMT_X);	strcat (record, text);	strcat (record, sep);
			gmt_ascii_format_col (GMT, text, global_ymin, GMT_OUT, GMT_Y);	strcat (record, text);	strcat (record, sep);
			gmt_ascii_format_col (GMT, text, global_ymax, GMT_OUT, GMT_Y);	strcat (record, text);	strcat (record, sep);
			gmt_ascii_format_col (GMT, text, global_zmin, GMT_OUT, GMT_Z);	strcat (record, text);	strcat (record, sep);
			gmt_ascii_format_col (GMT, text, global_zmax, GMT_OUT, GMT_Z);	strcat (record, text);
			GMT_Put_Record (API, GMT_WRITE_TEXT, record);
		}
	}
	else if (Ctrl->T.active) {
		if (Ctrl->T.mode & 2) {	/* Must do alpha trimming first */
			float *tmp_grid = NULL;
			char *file_ptr = grdfile;	/* To avoid a warning */
			if (gmt_M_file_is_memory (file_ptr)) {	/* Must operate on a copy since sorting is required */
				tmp_grid = gmt_M_memory_aligned (GMT, NULL, G->header->size, float);
				gmt_M_memcpy (tmp_grid, G->data, G->header->size, float);
			}
			else
				tmp_grid = G->data;
			gmt_sort_array (GMT, tmp_grid, G->header->size, GMT_FLOAT);	/* Sort so we can find quantiles */
			global_zmin = gmt_quantile_f (GMT, tmp_grid, 0.5 * Ctrl->T.alpha, G->header->size);			/* "Left" quantile */
			global_zmax = gmt_quantile_f (GMT, tmp_grid, 100.0-0.5* Ctrl->T.alpha, G->header->size);	/* "Right" quantile */
			if (GMT_Destroy_Data (API, &G) != GMT_NOERROR) {	/* Delayed destroy due to alpha trimming */
				Return (API->error);
			}
			if (gmt_M_file_is_memory (file_ptr))	/* Now free temp grid */
				gmt_M_free (GMT, tmp_grid);
		}
		if (Ctrl->T.mode & 1) {	/* Get a symmetrical range */
			if (Ctrl->T.inc > 0.0) {	/* Round limits first */
				global_zmin = floor (global_zmin / Ctrl->T.inc) * Ctrl->T.inc;
				global_zmax = ceil  (global_zmax / Ctrl->T.inc) * Ctrl->T.inc;
			}
			global_zmax = MAX (fabs (global_zmin), fabs (global_zmax));
			global_zmin = -global_zmax;
		}
		else {	/* Just use reported min/max values (possibly alpha-trimmed above) */
			if (Ctrl->T.inc > 0.0) {	/* Round limits first */
				global_zmin = floor (global_zmin / Ctrl->T.inc) * Ctrl->T.inc;
				global_zmax = ceil  (global_zmax / Ctrl->T.inc) * Ctrl->T.inc;
			}
		}
		sprintf (record, "-T");
		gmt_ascii_format_col (GMT, text, global_zmin, GMT_OUT, GMT_Z);	strcat (record, text);	strcat (record, "/");
		gmt_ascii_format_col (GMT, text, global_zmax, GMT_OUT, GMT_Z);	strcat (record, text);
		if (Ctrl->T.inc > 0.0) {
			strcat (record, "/");
			gmt_ascii_format_col (GMT, text, Ctrl->T.inc, GMT_OUT, GMT_Z);	strcat (record, text);
		}
		GMT_Put_Record (API, GMT_WRITE_TEXT, record);
	}
	else if ((Ctrl->I.active && Ctrl->I.status == GRDINFO_GIVE_REG_ROUNDED)) {
		global_xmin = floor (global_xmin / Ctrl->I.inc[GMT_X]) * Ctrl->I.inc[GMT_X];
		global_xmax = ceil  (global_xmax / Ctrl->I.inc[GMT_X]) * Ctrl->I.inc[GMT_X];
		global_ymin = floor (global_ymin / Ctrl->I.inc[GMT_Y]) * Ctrl->I.inc[GMT_Y];
		global_ymax = ceil  (global_ymax / Ctrl->I.inc[GMT_Y]) * Ctrl->I.inc[GMT_Y];
		if (gmt_M_is_geographic (GMT, GMT_IN)) {	/* Must make sure we don't get outside valid bounds */
			if (global_ymin < -90.0) {
				global_ymin = -90.0;
				GMT_Report (API, GMT_MSG_VERBOSE, "Warning: Using -I caused south to become < -90.  Reset to -90.\n");
			}
			if (global_ymax > 90.0) {
				global_ymax = 90.0;
				GMT_Report (API, GMT_MSG_VERBOSE, "Warning: Using -I caused north to become > +90.  Reset to +90.\n");
			}
			if (fabs (global_xmax - global_xmin) > 360.0) {
				GMT_Report (API, GMT_MSG_VERBOSE, "Warning: Using -I caused longitude range to exceed 360.  Reset to a range of 360.\n");
				global_xmin = (global_xmin < 0.0) ? -180.0 : 0.0;
				global_xmax = (global_xmin < 0.0) ? +180.0 : 360.0;
			}
		}
		sprintf (record, "-R");
		gmt_ascii_format_col (GMT, text, global_xmin, GMT_OUT, GMT_X);	strcat (record, text);	strcat (record, "/");
		gmt_ascii_format_col (GMT, text, global_xmax, GMT_OUT, GMT_X);	strcat (record, text);	strcat (record, "/");
		gmt_ascii_format_col (GMT, text, global_ymin, GMT_OUT, GMT_Y);	strcat (record, text);	strcat (record, "/");
		gmt_ascii_format_col (GMT, text, global_ymax, GMT_OUT, GMT_Y);	strcat (record, text);
		GMT_Put_Record (API, GMT_WRITE_TEXT, record);
	}

	if (!Ctrl->C.active && !Ctrl->T.active && projStr) {		/* Print the referencing info */
		GMT_Put_Record (API, GMT_WRITE_TEXT, projStr);
		gmt_M_str_free (projStr);
	}

	if (GMT_End_IO (API, GMT_OUT, 0) != GMT_NOERROR) {	/* Disables further data output */
		Return (API->error);
	}

	GMT_Report (API, GMT_MSG_VERBOSE, "Done!\n");
	Return (GMT_NOERROR);
}
