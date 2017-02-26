/*--------------------------------------------------------------------
 *	$Id: grdconvert.c 17449 2017-01-16 21:27:04Z pwessel $
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
 * Brief synopsis: grdconvert.c reads a grid file in one format and outputs it in another
 *
 * Author:	Paul Wessel
 * Date:	1-JAN-2010
 * Version:	5 API
 */

#define THIS_MODULE_NAME	"grdconvert"
#define THIS_MODULE_LIB		"core"
#define THIS_MODULE_PURPOSE	"Convert between different grid formats"
#define THIS_MODULE_KEYS	"<G{,>G}"

#include "gmt_dev.h"

#define GMT_PROG_OPTIONS "-RVf"

struct GRDCONVERT_CTRL {
	struct IO {
		bool active;
		char *file[2];
	} IO;
	struct N {	/* -N */
		bool active;
	} N;
};

GMT_LOCAL void *New_Ctrl (struct GMT_CTRL *GMT) {	/* Allocate and initialize a new control structure */
	struct GRDCONVERT_CTRL *C;
	
	C = gmt_M_memory (GMT, NULL, 1, struct GRDCONVERT_CTRL);
	
	/* Initialize values whose defaults are not 0/false/NULL */
	
	return (C);
}

GMT_LOCAL void Free_Ctrl (struct GMT_CTRL *GMT, struct GRDCONVERT_CTRL *C) {	/* Deallocate control structure */
	if (!C) return;
	gmt_M_str_free (C->IO.file[GMT_IN]);	
	gmt_M_str_free (C->IO.file[GMT_OUT]);	
	gmt_M_free (GMT, C);	
}

GMT_LOCAL int usage (struct GMTAPI_CTRL *API, int level) {
	int i;
	char **grdformats = gmt_grdformats_sorted (API->GMT);

	gmt_show_name_and_purpose (API, THIS_MODULE_LIB, THIS_MODULE_NAME, THIS_MODULE_PURPOSE);
	if (level == GMT_MODULE_PURPOSE) return (GMT_NOERROR);
	GMT_Message (API, GMT_TIME_NONE, "usage: grdconvert <ingrid>[=<id>[/<scale>/<offset>[/<nan>]]]\n\t<outgrid>[=<id>[/<scale>/<offset>[/<nan>]][:<driver>[/<dataType>]]] [-N]\n\t[%s] [%s] [%s]\n\n",
		GMT_Rgeo_OPT, GMT_V_OPT, GMT_f_OPT);

	if (level == GMT_SYNOPSIS) return (GMT_MODULE_SYNOPSIS);

	GMT_Message (API, GMT_TIME_NONE, "\t<ingrid> is the grid file to convert.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t<outgrid> is the new converted grid file.\n");
	GMT_Message (API, GMT_TIME_NONE, "\tscale and offset, if given, will multiply data by scale and add offset.\n");
	GMT_Message (API, GMT_TIME_NONE, "\n\tOPTIONS:\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-N Do NOT write the header (for native grids only - ignored otherwise).\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   Useful when creating files to be used by grdraster.\n");
	GMT_Option (API, "R,V,f,.");

	GMT_Message (API, GMT_TIME_NONE, "\nThe following grid file formats are supported:\n");
	for (i = 1; i < GMT_N_GRD_FORMATS; ++i) {
		if (!strstr (grdformats[i], "not supported"))
			GMT_Message (API, GMT_TIME_NONE, "\t%s\n", grdformats[i]);
	}
#ifdef HAVE_GDAL
	GMT_Message (API, GMT_TIME_NONE, "\n	When <id>=gd on output, the grid will be saved using the GDAL library.\n");
	GMT_Message (API, GMT_TIME_NONE, "	Specify <driver> and optionally <dataType>. Driver names are as in GDAL\n		(e.g., netCDF, GTiFF, etc.)\n");
	GMT_Message (API, GMT_TIME_NONE, "	<dataType> is u8|u16|i16|u32|i32|float32; i|u denote signed|unsigned\n		integer.  Default type is float32.\n");
	GMT_Message (API, GMT_TIME_NONE, "	Both driver names and data types are case insensitive.\n");
#endif
	return (GMT_MODULE_USAGE);
}

GMT_LOCAL int parse (struct GMT_CTRL *GMT, struct GRDCONVERT_CTRL *Ctrl, struct GMT_OPTION *options) {
	/* This parses the options provided to grdconvert and sets parameters in CTRL.
	 * Any GMT common options will override values set previously by other commands.
	 * It also replaces any file names specified as input or output with the data ID
	 * returned when registering these sources/destinations with the API.
	 */

	unsigned int n_errors = 0, n_in = 0;
	struct GMT_OPTION *opt = NULL;
	struct GMTAPI_CTRL *API = GMT->parent;

	for (opt = options; opt; opt = opt->next) {
		switch (opt->option) {

			case '<':	/* Input and Output files */
				if (n_in == 0 && gmt_check_filearg (GMT, '<', opt->arg, GMT_IN, GMT_IS_GRID))
					Ctrl->IO.file[n_in++] = strdup (opt->arg);
				else if (n_in == 1 && gmt_check_filearg (GMT, '>', opt->arg, GMT_OUT, GMT_IS_GRID))
					Ctrl->IO.file[n_in++] = strdup (opt->arg);
				else {
					n_in++;
					GMT_Report (API, GMT_MSG_NORMAL, "Syntax error: Specify only one input and one output file\n");
					n_errors++;
				}
				break;
			case '>':	/* Output file */
				if (gmt_check_filearg (GMT, '>', opt->arg, GMT_OUT, GMT_IS_GRID))
					Ctrl->IO.file[GMT_OUT] = strdup (opt->arg);
				else
					n_errors++;
				n_in++;
				break;

			/* Processes program-specific parameters */

			case 'N':
				Ctrl->N.active = true;
				break;

			default:	/* Report bad options */
				n_errors += gmt_default_error (GMT, opt->option);
				break;
		}
	}

	n_errors += gmt_M_check_condition (GMT, n_in != 2, "Syntax error: Must specify both input and output file names\n");

	return (n_errors ? GMT_PARSE_ERROR : GMT_NOERROR);
}

#define bailout(code) {gmt_M_free_options (mode); return (code);}
#define Return(code) {Free_Ctrl (GMT, Ctrl); gmt_end_module (GMT, GMT_cpy); bailout (code);}

int GMT_grdconvert (void *V_API, int mode, void *args) {
	int error = 0;
	unsigned int hmode, type[2];
	char fname[2][GMT_BUFSIZ];
	char   command[GMT_GRID_COMMAND_LEN320] = {""};
	struct GMT_GRID *Grid = NULL;
	struct GRDCONVERT_CTRL *Ctrl = NULL;
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

	/*---------------------------- This is the grdconvert main code ----------------------------*/

	if ((Grid = gmt_create_grid (API->GMT)) == NULL) Return (API->error);	/* Tmp grid only, no i/o is used */
	gmt_grd_init (GMT, Grid->header, options, false);
	hmode = (Ctrl->N.active) ? GMT_GRID_NO_HEADER : 0;
	gmt_M_err_fail (GMT, gmt_grd_get_format (GMT, Ctrl->IO.file[0], Grid->header, true), Ctrl->IO.file[0]);
	type[0] = Grid->header->type;
	strncpy (fname[0], Grid->header->name, GMT_BUFSIZ);
	gmt_M_err_fail (GMT, gmt_grd_get_format (GMT, Ctrl->IO.file[1], Grid->header, false), Ctrl->IO.file[1]);
	type[1] = Grid->header->type;
	strncpy (fname[1], Grid->header->name, GMT_BUFSIZ);
	gmt_free_grid (GMT, &Grid, true);	/* Free temp grid, Grid is now NULL */

	if (type[1] == GMT_GRID_IS_SD) {
		/* Golden Surfer format 7 is read-only */
		GMT_Report (API, GMT_MSG_NORMAL, "Writing unsupported: %s\n", GMT->session.grdformat[GMT_GRID_IS_SD]);
		Return (GMT_RUNTIME_ERROR);
	}

	if (gmt_M_is_verbose (GMT, GMT_MSG_VERBOSE)) {
		if (Ctrl->IO.file[0][0] == '=') strcpy (fname[0], "<stdin>");
		if (Ctrl->IO.file[1][0] == '=') strcpy (fname[1], "<stdout>");
		GMT_Report (API, GMT_MSG_VERBOSE, "Translating file %s (format %s) to file %s (format %s)\n",
		            fname[0], GMT->session.grdformat[type[0]], fname[1], GMT->session.grdformat[type[1]]);
		if (hmode && GMT->session.grdformat[type[1]][0] != 'c' && GMT->session.grdformat[type[1]][0] != 'n')
			GMT_Report (API, GMT_MSG_NORMAL, "No grd header will be written\n");
	}

	if ((Grid = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, Ctrl->IO.file[0], NULL)) == NULL) {	/* Get header only */
		Return (API->error);
	}

	if (GMT->common.R.active) {	/* Specified a subset */
		bool global = false;
		global = gmt_M_grd_is_global (GMT, Grid->header);
		if (!global && (GMT->common.R.wesn[XLO] < Grid->header->wesn[XLO] || GMT->common.R.wesn[XHI] > Grid->header->wesn[XHI])) error++;
		if (GMT->common.R.wesn[YLO] < Grid->header->wesn[YLO] || GMT->common.R.wesn[YHI] > Grid->header->wesn[YHI]) error++;
		if (error) {
			GMT_Report (API, GMT_MSG_NORMAL, "Subset exceeds data domain!\n");
			Return (GMT_RUNTIME_ERROR);
		}
		if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, GMT->common.R.wesn, Ctrl->IO.file[0], Grid) == NULL) {
			Return (API->error);	/* Get subset */
		}
	}
	else if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, Ctrl->IO.file[0], Grid) == NULL) {
		Return (API->error);	/* Get all */
	}

	Grid->header->type = type[1];

	/* When converting from netcdf to netcdf, we will keep the old command, so we need to make a copy of it now */
	command[0] = '\n';	command[1] = '\t';
	strcat(command, "(old cmd) ");
	strcat(command, Grid->header->command);

	gmt_grd_init (GMT, Grid->header, options, true);

	if (!GMT->common.R.active && ((type[0] >= GMT_GRID_IS_CB && type[0] <= GMT_GRID_IS_CD)  ||	/* That is, from netCDF to netCDF */
	                              (type[0] >= GMT_GRID_IS_NB && type[0] <= GMT_GRID_IS_ND)) &&
	                             ((type[1] >= GMT_GRID_IS_CB && type[1] <= GMT_GRID_IS_CD)  ||
	                              (type[1] >= GMT_GRID_IS_NB && type[1] <= GMT_GRID_IS_ND)) ) {
		/* Do nothing, which means the new grid will keep the command string of the old grid */
		if (GMT_Set_Comment (API, GMT_IS_GRID, GMT_COMMENT_IS_COMMAND, command, Grid))
			Return (API->error);
	}
	else if (GMT_Set_Comment (API, GMT_IS_GRID, GMT_COMMENT_IS_OPTION | GMT_COMMENT_IS_COMMAND, options, Grid))
		Return (API->error);

	if (GMT_Write_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, hmode, NULL, Ctrl->IO.file[1], Grid) != GMT_NOERROR)
		Return (API->error);

	Return (GMT_NOERROR);
}
