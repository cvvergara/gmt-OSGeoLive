/*--------------------------------------------------------------------
 *	$Id: mgd77info.c,v 1.75 2011/07/11 19:22:03 guru Exp $
 *
 *    Copyright (c) 2004-2011 by P. Wessel
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * mgd77info reads one or more MGD77 or MGD77+ files and report on the
 * extent of the file, number of data points etc.  Alternatively, it
 * can echo out the entire MGD77 header section or list columns that
 * are present.
 *
 * Author:	Paul Wessel
 * Date:	26-AUG-2004
 * Version:	1.0 Ideal based on the old gmtinfo.c
 *		2005-SEP-05: Added -P [PW]
 *		2005-OCT-07: Added -C,-I [PW]
 *		2006-MAR-31: Changed -H (header info) to -M (metadata)
 *		2007-JUN-14: Added -Me|h also
 *
 *
 */
 
#include "mgd77.h"
#include "mgd77_codes.h"

int main (int argc, char **argv)
{
	GMT_LONG i, id, rec, argno, col_listing = 3, length, id_col, t_col, x_col, y_col, saved_range, use;
	GMT_LONG header_flag = 0, n_paths, counter[MGD77_MAX_COLS], agency_listing = 1, quad_no, n_quad;
	
	GMT_LONG error = FALSE, greenwich = FALSE, quick_summary = FALSE, dump_raw_header = FALSE;
	GMT_LONG first = TRUE, col_summary = FALSE, read_file, dump_formatted_header = FALSE;
	GMT_LONG quick_agencies = FALSE, dump_e77_header = FALSE, dump_hist_header = FALSE;
	GMT_LONG quad[4] = {FALSE, FALSE, FALSE, FALSE};
	
	double this_dist, this_lon, this_lat, last_lon, last_lat, dx, dy, dlon, ds, lon_w;
	double xmin, xmax, xmin1, xmin2, xmax1, xmax2, ymin, ymax, this_time, tmin, tmax;
	double *dvalue[MGD77_MAX_COLS];
	
	char *tvalue[MGD77_MAX_COLS], buffer[BUFSIZ], **list = NULL;
	
	GMT_cal_rd rata_die;
	
	struct MGD77_CONTROL M, Out;
	struct MGD77_DATASET *D = NULL;
		
	argc = (int)GMT_begin (argc, argv);		/* Initialize GMT Machinery */
	GMT_get_time_system ("unix", &(gmtdefs.time_system));						/* MGD77+ uses GMT's Unix time epoch */
	GMT_init_time_system_structure (&(gmtdefs.time_system));
	
	/* Initialize MGD77 output order and other parameters*/
	
	MGD77_Init (&M);			/* Initialize MGD77 Machinery */
	MGD77_Init (&Out);		/* Initialize MGD77 Machinery */
	Out.fp = GMT_stdout;

	for (i =1; !error && i < argc; i++) {	/* Process input options */
		if (argv[i][0] != '-') continue;

		switch(argv[i][1]) {

			case 'V':
			case '\0':
				error += GMT_parse_common_options (argv[i], NULL, NULL, NULL, NULL);
				break;
						
			case 'C':	/* Get the short list [Default] */
				switch (argv[i][2]) {
					case 'm':
					case 'M':
						col_listing = 1;
						break;
					case 'e':
					case 'E':
						col_listing = 2;
						break;
					default:
						col_listing = 3;
						break;
				}
				col_summary = TRUE;
				break;
	
			case 'M':
				if (argv[i][2] == 'f') {
					dump_formatted_header = TRUE;
					header_flag = MGD77_Select_Header_Item (&M, &argv[i][3]);
					if (header_flag < 0) error++;
				}
				else if (argv[i][2] == 'r') {
					dump_raw_header = TRUE;
				}
				else if (argv[i][2] == 'e') {
					dump_e77_header = TRUE;
				}
				else if (argv[i][2] == 'h') {
					dump_hist_header = TRUE;
				}
				else {
					fprintf (stderr, "%s: Option -M Bad modifier (%c). Use -Mf|r|e|h!\n", GMT_program, argv[i][2]);
					exit (EXIT_FAILURE);
				}
				break;
	
			case 'I':
				MGD77_Process_Ignore (argv[i][1], &argv[i][2]);
				break;

			case 'E':	/* Get the short list [Default] */
				switch (argv[i][2]) {
					case 'm':
					case 'M':
						col_listing = 1;
						break;
					case 'e':
					case 'E':
						col_listing = 2;
						break;
					default:
						col_listing = 3;
						break;
				}
				quick_summary = TRUE;
				break;

			case 'L':	/* Get the list of institutions and vessels  */
				switch (argv[i][2]) {
					case 'a':
						agency_listing = 1;
						break;
					case 'v':
						agency_listing = 2;
						break;
					default:
						col_listing = 1;
						break;
				}
				quick_agencies = TRUE;
				break;
	
			default:		/* Options not recognized */
				error = TRUE;
				break;
		}
	}
	
	if (error) exit (EXIT_FAILURE);

	/* Check that the options selected are mutually consistent */
	
	if (GMT_give_synopsis_and_exit || argc == 1 || header_flag == 1) {	/* Display usage */
		fprintf(stderr,"mgd77info %s - Extract information about MGD77 files\n\n", MGD77_VERSION);
		fprintf(stderr,"usage: mgd77info <cruise(s)> [-C[m|e]] [-E[m|e]] [-I<code>] [-Mf[<item>]|r|e|h] [-L[v]] [-V]\n\n");
         
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
              
		MGD77_Cruise_Explain ();
		fprintf(stderr,"	OPTIONS:\n\n");
		fprintf(stderr,"	-C List abbreviations of all columns present for each cruise\n");
		fprintf(stderr,"	   Append m for listing just the MGD77 columns present\n");
		fprintf(stderr,"	   Append e for listing just any extra columns present\n");
		fprintf(stderr,"	-E Give the information summary of each cruise's geographical/temporal extent\n");
		fprintf(stderr,"	   Append m for counting just the number of non-NaN values for each MGD77 field\n");
		fprintf(stderr,"	   Append e for counting just the of non-NaN values for each extra field\n");
		fprintf(stderr,"	-M Print header items (and MGD77+ history).  Append type of presentation:\n");
		fprintf(stderr,"	     f: Print header items individually, one per line.  Append name of a particular\n");
		fprintf(stderr,"	        item (e.g. Port_of_Departure), all [Default], or - to see a list of items.\n");
		fprintf(stderr,"	        You can also use the number of the item.\n");
		fprintf(stderr,"	     r: Display raw original MGD77 header records.\n");
		fprintf(stderr,"	     e: Display the MGD77+ file's E77 status.\n");
		fprintf(stderr,"	     h: Display the MGD77+ file's history.\n");
		if (header_flag == 1) MGD77_List_Header_Items (&M);
		fprintf(stderr,"	-I Ignore certain data file formats from consideration. Append combination of act to ignore\n");
		fprintf(stderr,"	   (a) MGD77 ASCII, (c) MGD77+ netCDF, or (t) plain table files. [Default ignores none]\n");
		fprintf(stderr,"	-L Just list all the institutions and their 2-character GEODAS codes.  Append v to also\n");
		fprintf(stderr,"	   display the vessels and their 4-character codes for each institution\n");
		fprintf(stderr,"	-V verbose, report progress\n");
		exit (EXIT_FAILURE);
	}

	if (!(dump_raw_header + dump_e77_header + dump_hist_header + quick_summary + col_summary + dump_formatted_header + quick_agencies ) == 1) {
		fprintf(stderr, "%s: ERROR: Specify one of -C, -E, -L, or -M\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (quick_agencies) {	/* Just display the list and exit */
		(agency_listing == 2) ? printf ("CODE\tINSTITUTION/VESSEL\n") : printf ("CODE\tINSTITUTION\n");
		for (id = i = 0; id < MGD77_N_AGENCIES; id++) {
			printf ("%s = %s\n", MGD77_agency[id].code, MGD77_agency[id].name);
			for (; agency_listing == 2 && i < MGD77_N_VESSELS && MGD77_vessel[i].agent == id; i++) {
				printf ("%s\t-> %s\n", MGD77_vessel[i].code, MGD77_vessel[i].name);
			}
		}
		exit (EXIT_SUCCESS);
	}

	n_paths = MGD77_Path_Expand (&M, argv, argc, &list);	/* Get list of requested IDs */
	
	if (n_paths == 0) {
		fprintf(stderr, "%s: ERROR: No cruises given\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	
	read_file = (quick_summary || dump_raw_header);
	
	saved_range = GMT_io.geo.range;	/* We may have to reset thisso keep a copy */
	GMT_io.out_col_type[0] = GMT_IS_LON;	GMT_io.out_col_type[1] = GMT_IS_LAT;
	GMT_io.out_col_type[2] = M.time_format;	
	if (quick_summary) {
		sprintf (buffer, "#Cruise %sID      %sWest    %sEast    %sSouth   %sNorth   %sStartTime%s%sEndTime%s%s%sDist%snRec",
		gmtdefs.field_delimiter, gmtdefs.field_delimiter, gmtdefs.field_delimiter, gmtdefs.field_delimiter, gmtdefs.field_delimiter, gmtdefs.field_delimiter,
		gmtdefs.field_delimiter, gmtdefs.field_delimiter, gmtdefs.field_delimiter, gmtdefs.field_delimiter, gmtdefs.field_delimiter, gmtdefs.field_delimiter);
		GMT_fputs (buffer, GMT_stdout);
	}	

	use = (M.original || M.format != MGD77_FORMAT_CDF) ? MGD77_ORIG : MGD77_REVISED;
	
	for (argno = 0; argno < n_paths; argno++) {		/* Process each ID */
	
		if (MGD77_Open_File (list[argno], &M, MGD77_READ_MODE)) continue;

		if (gmtdefs.verbose) fprintf (stderr, "%s: Now processing cruise %s\n", GMT_program, list[argno]);
		
		D = MGD77_Create_Dataset ();
		
		if (read_file && MGD77_Read_File (list[argno], &M, D)) {
			fprintf (stderr, "%s: Error reading header & data for cruise %s\n", GMT_program, list[argno]);
			exit (EXIT_FAILURE);
		}
		if (!read_file && MGD77_Read_Header_Record (list[argno], &M, &D->H)) {
			fprintf (stderr, "%s: Error reading header sequence for cruise %s\n", GMT_program, list[argno]);
			exit (EXIT_FAILURE);
		}

		if (dump_hist_header) {	/* Dump of MGD77+ history */
			sprintf (buffer, "%s: %s", list[argno], D->H.history);
			GMT_fputs (buffer, GMT_stdout);
			MGD77_Close_File (&M);
			MGD77_Free (D);
			continue;
		}
		if (dump_e77_header) {	/* Dump of e77 header status */
			if (D->H.E77 && strlen(D->H.E77) > 0) {
				sprintf (buffer, "%s: %s\n", list[argno], D->H.E77);
				GMT_fputs (buffer, GMT_stdout);
			}
			else {
				sprintf (buffer, "%s: E77 not applied\n", list[argno]);
				GMT_fputs (buffer, GMT_stdout);
			}
			MGD77_Close_File (&M);
			MGD77_Free (D);
			continue;
		}
		if (dump_formatted_header) {	/* Dump of header items, one per line */
			MGD77_Dump_Header_Params (&M, D->H.mgd77[use]);	
			MGD77_Close_File (&M);
			MGD77_Free (D);
			continue;
		}
		t_col = MGD77_Get_Column ("time", &M);
		x_col = MGD77_Get_Column ("lon", &M);
		y_col = MGD77_Get_Column ("lat", &M);
		id_col = MGD77_Get_Column ("id", &M);
		
		if (first && quick_summary) {	/* Output all column headers */
			for (i = 0; i < M.n_out_columns; i++) {
				if (i == id_col || i == t_col || i == x_col || i == y_col) continue;
				sprintf (buffer,"%s%s", gmtdefs.field_delimiter, D->H.info[M.order[i].set].col[M.order[i].item].abbrev);
				GMT_fputs (buffer, GMT_stdout);
			}
			GMT_fputs ("\n", GMT_stdout);
		}
		
		if (col_summary) {	/* Just list names and info for any extra columns */
			for (i = 0, first = TRUE; i < M.n_out_columns; i++) {
				if (i == id_col || i == t_col || i == x_col || i == y_col) continue;
				if (!first) GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
				if (((col_listing & 1) && M.order[i].set == 0) || ((col_listing & 2) && M.order[i].set == 1)) {
					GMT_fputs (D->H.info[M.order[i].set].col[M.order[i].item].abbrev, GMT_stdout);
					first = FALSE;
				}
			}
			if (first) GMT_fputs ("No columns matching selection found!", GMT_stdout);
			GMT_fputs ("\n", GMT_stdout);
			MGD77_Close_File (&M);
			MGD77_Free (D);
			continue;
		}
			
		if (dump_raw_header) {	/* Write entire MGD77 header */
			GMT_fputs ("-------------------------------", GMT_stdout);
			sprintf (buffer, " Cruise: %8s ", M.NGDC_id);	GMT_fputs (buffer, GMT_stdout);
			GMT_fputs ("-------------------------------\n", GMT_stdout);
			MGD77_Write_Header_Record_m77 ("", &Out, &D->H);
			GMT_fputs ("----------------------------------------", GMT_stdout);
			sprintf (buffer, "----------------------------------------\n");	GMT_fputs (buffer, GMT_stdout);
			if (M.format == MGD77_FORMAT_CDF) {
				sprintf (buffer, "%s\n", D->H.history);
				GMT_fputs (buffer, GMT_stdout);
				for (i = 0; i < M.n_out_columns; i++) {
					if ((M.order[i].set == MGD77_CDF_SET)) {
						sprintf (buffer, "> %s%s%s%s%s%s%s", D->H.info[MGD77_CDF_SET].col[M.order[i].item].abbrev, gmtdefs.field_delimiter,
						D->H.info[MGD77_CDF_SET].col[M.order[i].item].name, gmtdefs.field_delimiter,
						D->H.info[MGD77_CDF_SET].col[M.order[i].item].units, gmtdefs.field_delimiter,
						D->H.info[MGD77_CDF_SET].col[M.order[i].item].comment);
						GMT_fputs ("\n", GMT_stdout);
					}
				}
			}
			GMT_fputs ("\n", GMT_stdout);
		}
		
		GMT_make_dnan (tmin);
		GMT_make_dnan (tmax);
		this_dist = this_lon = this_lat = ds = this_time = 0.0;
		xmin1 = xmin2 = 360.0;
		xmax1 = xmax2 = -360.0;
		ymin = 180.0;
		ymax = -180.0;
		memset ((void *)quad, 0, 4*sizeof(GMT_LONG));	/* Set al to FALSE */
		greenwich = FALSE;
		memset ((void *) counter, 0, (size_t)(MGD77_MAX_COLS * sizeof (GMT_LONG)));
	
		for (i = 0; i < MGD77_MAX_COLS; i++) {
			dvalue[i] = (double *)D->values[i];
			tvalue[i] = (char *)D->values[i];
		}
		
		/* Start processing data */
	
		for (rec = 0; rec < D->H.n_records; rec++) {		/* While able to read a data record */
		
			/* Get min and max time */
			if (t_col >= 0 && !GMT_is_dnan(dvalue[t_col][rec])) {
				if (GMT_is_dnan(tmin) && GMT_is_dnan(tmax)) tmin = tmax = dvalue[t_col][rec];
				this_time = dvalue[t_col][rec];
				tmin = MIN (this_time, tmin);
				tmax = MAX (this_time, tmax);
			}

			/* Compute accumulated distance along track (Flat Earth) */
			last_lon  = this_lon;
			last_lat  = this_lat;
			this_lon  = lon_w = dvalue[x_col][rec];
			this_lat  = dvalue[y_col][rec];
			if (this_lon < 0.0) this_lon += 360.0;	/* Start off with everything in 0-360 range */
			xmin1 = MIN (this_lon, xmin1);
			xmax1 = MAX (this_lon, xmax1);
			quad_no = (GMT_LONG)floor (this_lon/90.0);	/* Yields quadrants 0-3 */
			if (quad_no == 4) quad_no = 0;		/* When this_lon == 360.0 */
			quad[quad_no] = TRUE;
			if (lon_w > 180.0) this_lon -= 360.0;	/* For -180/+180 range */
			xmin2 = MIN (this_lon, xmin2);
			xmax2 = MAX (this_lon, xmax2);
			if (rec > 0) {	/* Need a previous point to calculate distance, speed, and heading */
				dlon = this_lon - last_lon;
				if (fabs (dlon) > 180.0) {
					greenwich = TRUE;
					dlon = copysign ((360.0 - fabs (dlon)), dlon);
				}
				dx = dlon * cosd (0.5 * (this_lat + last_lat));
				dy = this_lat - last_lat;
				ds = project_info.DIST_KM_PR_DEG * hypot (dx, dy);
				this_dist += ds;
			}
			ymin = MIN (this_lat, ymin);
			ymax = MAX (this_lat, ymax);
			
			/* Count the number of non-NaN observations */
			
			for (i = 1; i < M.n_out_columns; i++) {
				if (i == id_col || i == t_col || i == x_col || i == y_col) continue;
				if ((length = D->H.info[M.order[i].set].col[M.order[i].item].text)) {
					if (strncmp (&tvalue[i][rec*length], ALL_NINES, (size_t)length)) counter[i]++;
				}
				else
					if (!GMT_is_dnan (dvalue[i][rec])) counter[i]++;
			}
		}

		GMT_io.geo.range = saved_range;	/* We reset this each time */
		n_quad = quad[0] + quad[1] + quad[2] + quad[3];	/* How many quadrants had data */
		if (quad[0] && quad[3]) {	/* Longitudes on either side of Greenwich only, must use -180/+180 notation */
			xmin = xmin2;
			xmax = xmax2;
			GMT_io.geo.range = 2;	/* Override this setting explicitly */
		}
		else if (quad[1] && quad[2]) {	/* Longitudes on either side of the date line, must user 0/360 notation */
			xmin = xmin1;
			xmax = xmax1;
			GMT_io.geo.range = 0;	/* Override this setting explicitly */
		}
		else if (n_quad == 2 && ((quad[0] && quad[2]) || (quad[1] && quad[3]))) {	/* Funny quadrant gap, pick shortest longitude extent */
			if ((xmax1 - xmin1) < (xmax2 - xmin2)) {	/* 0/360 more compact */
				xmin = xmin1;
				xmax = xmax1;
				GMT_io.geo.range = 0;	/* Override this setting explicitly */
			}
			else {						/* -180/+180 more compact */
				xmin = xmin2;
				xmax = xmax2;
				GMT_io.geo.range = 2;	/* Override this setting explicitly */
			}
		}
		else {						/* Either will do, use default settings */
			xmin = xmin1;
			xmax = xmax1;
		}
		if (xmin > xmax) xmin -= 360.0;
		if (xmin < 0.0 && xmax < 0.0) xmin += 360.0, xmax += 360.0;

		if (GMT_is_dnan(tmin) || GMT_is_dnan(tmax)) {
			int yy[2], mm[2], dd[2];
			if (gmtdefs.verbose) fprintf (stderr, "%s warning: cruise %s no time records.\n",GMT_program, M.NGDC_id);
			yy[0] = (!D->H.mgd77[use]->Survey_Departure_Year[0] || !strncmp (D->H.mgd77[use]->Survey_Departure_Year, ALL_BLANKS, (size_t)4)) ? 0 : atoi (D->H.mgd77[use]->Survey_Departure_Year);
			yy[1] = (!D->H.mgd77[use]->Survey_Arrival_Year[0] || !strncmp (D->H.mgd77[use]->Survey_Arrival_Year, ALL_BLANKS, (size_t)4)) ? 0 : atoi (D->H.mgd77[use]->Survey_Arrival_Year);
			mm[0] = (!D->H.mgd77[use]->Survey_Departure_Month[0] || !strncmp (D->H.mgd77[use]->Survey_Departure_Month, ALL_BLANKS, (size_t)2)) ? 1 : atoi (D->H.mgd77[use]->Survey_Departure_Month);
			mm[1] = (!D->H.mgd77[use]->Survey_Arrival_Month[0] || !strncmp (D->H.mgd77[use]->Survey_Arrival_Month, ALL_BLANKS, (size_t)2)) ? 1 : atoi (D->H.mgd77[use]->Survey_Arrival_Month);
			dd[0] = (!D->H.mgd77[use]->Survey_Departure_Day[0] || !strncmp (D->H.mgd77[use]->Survey_Departure_Day, ALL_BLANKS, (size_t)2)) ? 1 : atoi (D->H.mgd77[use]->Survey_Departure_Day);
			dd[1] = (!D->H.mgd77[use]->Survey_Arrival_Day[0] || !strncmp (D->H.mgd77[use]->Survey_Arrival_Day, ALL_BLANKS, (size_t)2)) ? 1 : atoi (D->H.mgd77[use]->Survey_Arrival_Day);
			if (! (yy[0] == 0 && yy[1] == 0)) {	/* With year we can do something */
				rata_die = GMT_rd_from_gymd (yy[0], mm[0], dd[0]);
				tmin = GMT_rdc2dt (rata_die, 0.0);
				rata_die = GMT_rd_from_gymd (yy[1], mm[1], dd[1]);
				tmax = GMT_rdc2dt (rata_die, 0.0);
			}
		}			
		if (quick_summary) {
			sprintf (buffer,"%8s%s%8s%s", M.NGDC_id, gmtdefs.field_delimiter, D->H.mgd77[use]->Survey_Identifier, gmtdefs.field_delimiter);
			GMT_fputs (buffer, GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, xmin, 0);	GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, xmax, 0);	GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, ymin, 1);	GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, ymax, 1);	GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
			if (!GMT_is_dnan(tmin) && !GMT_is_dnan(tmax)) {
				GMT_ascii_output_one (GMT_stdout, tmin, 2);	GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, tmax, 2);	GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);						
			} else {
				sprintf (buffer, "%4s-%2s-%2s%s%4s-%2s-%2s%s",
				D->H.mgd77[use]->Survey_Departure_Year, D->H.mgd77[use]->Survey_Departure_Month, D->H.mgd77[use]->Survey_Departure_Day, gmtdefs.field_delimiter,
				D->H.mgd77[use]->Survey_Arrival_Year, D->H.mgd77[use]->Survey_Arrival_Month, D->H.mgd77[use]->Survey_Arrival_Day, gmtdefs.field_delimiter);
				GMT_fputs (buffer, GMT_stdout);
			}
			sprintf (buffer, "%ld%s%ld", (GMT_LONG)irint (this_dist), gmtdefs.field_delimiter, D->H.n_records);
			GMT_fputs (buffer, GMT_stdout);
			for (i = 1; i < M.n_out_columns; i++) {
				if (i == id_col || i == t_col || i == x_col || i == y_col) continue;
				if (((col_listing & 1) && M.order[i].set == 0) || ((col_listing & 2) && M.order[i].set == 1)) {
					sprintf (buffer,"%s%ld",	gmtdefs.field_delimiter, counter[i]);
					GMT_fputs (buffer, GMT_stdout);
				}
			}
			GMT_fputs ("\n", GMT_stdout);
		}
		MGD77_Free (D);
	}
		
	MGD77_Path_Free ((int)n_paths, list);
	MGD77_end (&M);

	exit (EXIT_SUCCESS);
}
