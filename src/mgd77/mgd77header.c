/*--------------------------------------------------------------------
 *	$Id: mgd77header.c 10060 2013-06-24 00:04:12Z pwessel $
 *
 *    Copyright (c) 2004-2013 by P. Wessel and Michael Chandler
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * mgd77header.c reads an NGDC A77 file, determines temporal and spatial extents,
 * ten degree boxes, data columns present. Also reads header items from a file (-H<file).
 * Output: Header items determined from data and read from input are output in 
 * H77 or raw header format (-M<r|f>).
 * 
 * Author:	Michael Chandler
 * Date:	23-MAY-2012
 *
 */
 
#include "mgd77.h"
#include "mgd77_codes.h"
#include "math.h"
#include "time.h"

int MGD77_Read_File_nohdr (char *, struct MGD77_CONTROL *, struct MGD77_DATASET *);
int MGD77_Read_File_asc_nohdr (char *, struct MGD77_CONTROL *, struct MGD77_DATASET *);
int MGD77_Read_File_cdf_nohdr (char *, struct MGD77_CONTROL *, struct MGD77_DATASET *);
int MGD77_Read_Header_Record_nohdr (char *, struct MGD77_CONTROL *, struct MGD77_HEADER *);
int MGD77_Read_Header_Record_m77_nohdr (char *, struct MGD77_CONTROL *, struct MGD77_HEADER *);
int MGD77_Read_Header_Record_m77t_nohdr (char *, struct MGD77_CONTROL *, struct MGD77_HEADER *);

EXTERN_MSC int MGD77_Read_Data_asc (char *file, struct MGD77_CONTROL *F, struct MGD77_DATASET *S);
EXTERN_MSC int MGD77_Read_Data_cdf (char *file, struct MGD77_CONTROL *F, struct MGD77_DATASET *S);
EXTERN_MSC int MGD77_Read_Header_Record_cdf (char *file, struct MGD77_CONTROL *F, struct MGD77_HEADER *H);
EXTERN_MSC int MGD77_Decode_Header_m77t (struct MGD77_HEADER_PARAMS *P, char *record);

int main (int argc, char **argv)
{
	GMT_LONG i, id, rec, argno, col_listing = 3, length, id_col, t_col, x_col, y_col, saved_range, use;
	GMT_LONG header_flag = 0, n_paths, counter[MGD77_MAX_COLS], agency_listing = 1, quad_no, n_quad;
	
	GMT_LONG error = FALSE, greenwich = FALSE, quick_summary = FALSE, dump_raw_header = FALSE;
	GMT_LONG first = TRUE, col_summary = FALSE, read_file, dump_formatted_header = FALSE;
	GMT_LONG quick_agencies = FALSE, dump_e77_header = FALSE, dump_hist_header = FALSE;
	GMT_LONG quad[4] = {FALSE, FALSE, FALSE, FALSE}, b_col, m_col, f_col;
	GMT_LONG tendeg[36][18], tenx, teny, nten = 0, tquad;
	
	double this_dist, this_lon, this_lat, last_lon, last_lat, dx, dy, dlon, ds, lon_w;
	double xmin, xmax, xmin1, xmin2, xmax1, xmax2, ymin, ymax, this_time, tmin, tmax;
	double *dvalue[MGD77_MAX_COLS];
    
	FILE *infile = NULL;
    
	time_t tt;
	struct tm *tod = NULL;
	time (&tt);
	tod = localtime(&tt);
	
	char *tvalue[MGD77_MAX_COLS], buffer[BUFSIZ], **list = NULL, name[BUFSIZ], value[BUFSIZ], params[BUFSIZ], hinfile[BUFSIZ], line[BUFSIZ];

	struct MGD77_CONTROL M, Out;
	struct MGD77_DATASET *D = NULL;
		
	argc = (int)GMT_begin (argc, argv);		/* Initialize GMT Machinery */
	
	/* Initialize MGD77 output order and other parameters*/
	
	MGD77_Init (&M);		/* Initialize MGD77 Machinery */
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
	
			case 'H':
				strcpy (hinfile, &argv[i][2]);
				infile = fopen (hinfile, "r");
				if (infile == NULL) {
					fprintf (stderr, "*** Can't open input file (%s) ***\n", hinfile);
					exit(0);
				}
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
		fprintf (stderr, "mgd77header %s - Extract information about MGD77 files\n\n", MGD77_VERSION);
		fprintf (stderr, "usage: mgd77header <cruise(s)> [-C[m|e]] [-E[m|e]] [-H<file> [-I<code>] [-Mf[<item>]|r|e|h] [-L[v]] [-V]\n\n");
         
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
              
		MGD77_Cruise_Explain ();
		fprintf (stderr, "\tOPTIONS:\n\n");
		fprintf (stderr, "\t-C List abbreviations of all columns present for each cruise\n");
		fprintf (stderr, "\t   Append m for listing just the MGD77 columns present\n");
		fprintf (stderr, "\t   Append e for listing just any extra columns present\n");
		fprintf (stderr, "\t-E Give the information summary of each cruise's geographical/temporal extent\n");
		fprintf (stderr, "\t   Append m for counting just the number of non-NaN values for each MGD77 field\n");
		fprintf (stderr, "\t   Append e for counting just the of non-NaN values for each extra field\n");
		fprintf (stderr, "\t-H Read and assign header values from a file. Each input file row gives an exact\n");
		fprintf (stderr, "\t   header_field_name, space or tab, and header value. Values are read according to\n");
		fprintf (stderr, "\t   NGDC's MGD77 header format specification.\n\t\te.g.,\n\t\tSource_Institution Univ. of Hawaii\n");
		fprintf (stderr, "\t\tPort_of_Arrival Honolulu, HAWAII\n\t\t...\n	   See mgd77info -Mf output for recognized header field names.\n");
		fprintf (stderr, "\t-M Print header items (and MGD77+ history).  Append type of presentation:\n");
		fprintf (stderr, "\t     f: Print header items individually, one per line.  Append name of a particular\n");
		fprintf (stderr, "\t        item (e.g., Port_of_Departure), all [Default], or - to see a list of items.\n");
		fprintf (stderr, "\t        You can also use the number of the item.\n");
		fprintf (stderr, "\t     r: Display raw original MGD77 header records.\n");
		fprintf (stderr, "\t     e: Display the MGD77+ file's E77 status.\n");
		fprintf (stderr, "\t     h: Display the MGD77+ file's history.\n");
		if (header_flag == 1) MGD77_List_Header_Items (&M);
		fprintf (stderr, "\t-I Ignore certain data file formats from consideration. Append combination of act to ignore\n");
		fprintf (stderr, "\t   (a) MGD77 ASCII, (c) MGD77+ netCDF, (m) MGD77T ASCII, or (t) plain table files. [Default ignores none]\n");
		fprintf (stderr, "\t-L Just list all the institutions and their 2-character GEODAS codes.  Append v to also\n");
		fprintf (stderr, "\t   display the vessels and their 4-character codes for each institution\n");
		fprintf (stderr, "\t-V verbose, report progress\n");
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
	
	read_file = (quick_summary || dump_raw_header || dump_formatted_header);
	
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
		
		if (read_file && MGD77_Read_File_nohdr (list[argno], &M, D)) {
			fprintf (stderr, "%s: Error reading data for cruise %s\n", GMT_program, list[argno]);
			exit (EXIT_FAILURE);
		}
		if (!read_file && MGD77_Read_Header_Record_nohdr (list[argno], &M, &D->H)) {
			fprintf (stderr, "%s: Error reading header sequence for cruise %s\n", GMT_program, list[argno]);
			exit (EXIT_FAILURE);
		}

		t_col = MGD77_Get_Column ("time", &M);
		x_col = MGD77_Get_Column ("lon", &M);
		y_col = MGD77_Get_Column ("lat", &M);
		b_col = MGD77_Get_Column ("depth", &M);
		m_col = MGD77_Get_Column ("mag", &M);
		f_col = MGD77_Get_Column ("faa", &M);
		id_col = MGD77_Get_Column ("id", &M);
		
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

        	/* Check which ten-degree bins were crossed */
        	for (tenx = 0; tenx < 36; tenx++) {
			for (teny = 0; teny < 18; teny++)
  				tendeg[tenx][teny] = 0;
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
			teny = floor ((this_lat/10)+9);
			tenx = floor (this_lon/10);
			if (!tendeg[tenx][teny]) {
				tendeg[tenx][teny] = 1;
				/* printf ("lon: %.5f lat: %.5f\ttendeg[%d][%d] = %d\n",this_lon,this_lat,tenx,teny,tendeg[tenx][teny]); */
				nten++;
			}
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
			if (gmtdefs.verbose) fprintf (stderr, "%s abort: cruise %s no time records.\n",GMT_program, M.NGDC_id);
			exit (EXIT_FAILURE);
		}

		/* Store data limits in header */
		sprintf (value,"%.0f",floor(ymin));
		sprintf (name,"Bottommost_Latitude");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],value,MGD77_Header_Lookup[id].length);
		sprintf (value,"%.0f",ceil(ymax));
		sprintf (name,"Topmost_Latitude");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],value,MGD77_Header_Lookup[id].length);
		sprintf (value,"%.0f",floor(xmin));
		sprintf (name,"Leftmost_Longitude");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],value,MGD77_Header_Lookup[id].length);
		sprintf (value,"%.0f",ceil(xmax));
		sprintf (name,"Rightmost_Longitude");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],value,MGD77_Header_Lookup[id].length);

		/* Add ten degree identifier string to header */
		sprintf (value,"%02ld",nten);
		sprintf (name,"Number_of_Ten_Degree_Identifiers");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],value,MGD77_Header_Lookup[id].length);
		value[0] = '\0';
		for (tenx = 0; tenx < 36; tenx++) {
			for (teny = 0; teny < 18; teny++) {
				if (tendeg[tenx][teny]) {
					if (tenx >= 18) {
						tquad = 5;
						if (teny >= 9) 
							tquad = 7;
					} else {
						tquad = 3;
						if (teny >= 9)
							tquad = 1;
					}
					sprintf (value,"%s%ld%.0f%02ld, ",value,tquad,fabs(teny-9),tenx);
				}
			}
		}
		sprintf (value, "%s9999", value);
		sprintf (name,"Ten_Degree_Identifier");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],value,MGD77_Header_Lookup[id].length);

		/* Add survey identifier and file creation date */
		sprintf (value,"%.8s",list[argno]);
		sprintf (name,"Survey_Identifier");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],value,MGD77_Header_Lookup[id].length);
		sprintf (value,"MGD77");
		sprintf (name,"Format_Acronym");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],value,MGD77_Header_Lookup[id].length);
		sprintf (value,"11111");
		if (counter[b_col]) value[0] = '5';
		if (counter[m_col]) value[1] = '5';     
		if (counter[f_col]) value[2] = '5';
		strcpy(params,value);
		sprintf (name,"Parameters_Surveyed_Code");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],value,MGD77_Header_Lookup[id].length);
		sprintf (value, "%d",1900+tod->tm_year);
		sprintf (name,"File_Creation_Year");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],value,MGD77_Header_Lookup[id].length);
		sprintf (value, "%02d",1+tod->tm_mon);
		sprintf (name,"File_Creation_Month");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],value,MGD77_Header_Lookup[id].length);
		sprintf (value, "%02d",tod->tm_mday);
		sprintf (name,"File_Creation_Day");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],value,MGD77_Header_Lookup[id].length);
		GMT_ascii_format_one (value, tmin, GMT_IS_ABSTIME);
		sprintf (name,"Survey_Departure_Year");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],value,4);
		sprintf (name,"Survey_Departure_Month");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],&value[5],2);
		sprintf (name,"Survey_Departure_Day");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],&value[8],2);
		GMT_ascii_format_one (value, tmax, GMT_IS_ABSTIME);
		sprintf (name,"Survey_Arrival_Year");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],value,4);
		sprintf (name,"Survey_Arrival_Month");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],&value[5],2);
		sprintf (name,"Survey_Arrival_Day");
		id = MGD77_Get_Header_Item (&M, name);
		strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],&value[8],2);
        
		/* Copy header items from header input file */
		if (infile) {
			while (fgets (line,BUFSIZ,infile)) {
				sscanf (line, "%s %[^\t\n]", name, value);
				if (! strlen(value)) continue;
				if (params[0] == '1' && !strncmp(name,"Bathyme",7)) continue;
				if (params[1] == '1' && !strncmp(name,"Magneti",7)) continue;
				if (params[2] == '1' && !strncmp(name,"Gravity",7)) continue;
				id = MGD77_Get_Header_Item (&M, name);
				strncpy (MGD77_Header_Lookup[id].ptr[MGD77_M77_SET],value,MGD77_Header_Lookup[id].length);
			}
		}
        
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

		if (dump_formatted_header) {	/* Dump of header items, one per line */
			MGD77_Dump_Header_Params (&M, D->H.mgd77[use]);	
			MGD77_Close_File (&M);
			continue;
		}
        
		if (dump_raw_header) {	/* Write entire MGD77 header */
			if (infile == NULL) {
				GMT_fputs ("-------------------------------", GMT_stdout);
				sprintf (buffer, " Cruise: %8s ", M.NGDC_id);	GMT_fputs (buffer, GMT_stdout);
				GMT_fputs ("-------------------------------\n", GMT_stdout);
			}
			MGD77_Write_Header_Record_m77 ("", &Out, &D->H);
			if (infile == NULL) {
				GMT_fputs ("----------------------------------------", GMT_stdout);
				sprintf (buffer, "----------------------------------------\n");	GMT_fputs (buffer, GMT_stdout);
			}
			if (M.format == MGD77_FORMAT_CDF) {
				sprintf (buffer, "%s\n", D->H.history);
				GMT_fputs (buffer, GMT_stdout);
				for (i = 0; i < M.n_out_columns; i++) {
					if (M.order[i].set == MGD77_CDF_SET) {
						sprintf (buffer, "> %s%s%s%s%s%s%s", D->H.info[MGD77_CDF_SET].col[M.order[i].item].abbrev, gmtdefs.field_delimiter,
							D->H.info[MGD77_CDF_SET].col[M.order[i].item].name, gmtdefs.field_delimiter,
							D->H.info[MGD77_CDF_SET].col[M.order[i].item].units, gmtdefs.field_delimiter,
							D->H.info[MGD77_CDF_SET].col[M.order[i].item].comment);
						GMT_fputs ("\n", GMT_stdout);
					}
				}
			}
			if (infile == NULL) GMT_fputs ("\n", GMT_stdout);
		}

		if (quick_summary) {
			sprintf (buffer,"%8s%s%8s%s", M.NGDC_id, gmtdefs.field_delimiter, D->H.mgd77[use]->Survey_Identifier, gmtdefs.field_delimiter);
			GMT_fputs (buffer, GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, xmin, 0);	GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, xmax, 0);	GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, ymin, 1);	GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, ymax, 1);	GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, tmin, 2);	GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, tmax, 2);	GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);						
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
		if (dump_hist_header) {	/* Dump of MGD77+ history */
			sprintf (buffer, "%s: %s", list[argno], D->H.history);
			GMT_fputs (buffer, GMT_stdout);
			MGD77_Close_File (&M);
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
			continue;
		}

		MGD77_Free (D);
	}
		
	MGD77_Path_Free ((int)n_paths, list);
	MGD77_end (&M);

	exit (EXIT_SUCCESS);
}

int MGD77_Read_File_nohdr (char *file, struct MGD77_CONTROL *F, struct MGD77_DATASET *S)
{
	int err = 0;
    
	switch (F->format) {
		case MGD77_FORMAT_M77:	/* Plain MGD77 file */
		case MGD77_FORMAT_M7T:	/* Plain MGD77T file */
		case MGD77_FORMAT_TBL:	/* Plain ascii table */
			err = MGD77_Read_File_asc_nohdr (file, F, S);
			break;
		case MGD77_FORMAT_CDF:	/* netCDF MGD77 file */
			err = MGD77_Read_File_cdf_nohdr (file, F, S);
			break;
		default:
			fprintf (stderr, "%s: Bad format (%d)!\n", GMT_program, F->format);
			err = MGD77_UNKNOWN_FORMAT;
	}
    
	return (err);
}


int MGD77_Read_File_asc_nohdr (char *file, struct MGD77_CONTROL *F, struct MGD77_DATASET *S)	  /* Will read all MGD77 records in current file */
{
	int err;
    
	err = MGD77_Open_File (file, F, MGD77_READ_MODE);
	if (err) return (err);
    
	MGD77_Select_All_Columns (F, &S->H);	/* We know we only deal with items from set 0 here */

	err = MGD77_Read_Header_Record_nohdr (file, F, &S->H);  /* Will read the entire 24-section header structure */
	if (err) return (err);
    
	err = MGD77_Read_Data_asc (file, F, S);	  /* Will read all MGD77 records in current file */
	if (err) return (err);
    
	MGD77_Close_File (F);
    
	return (MGD77_NO_ERROR);
}

int MGD77_Read_File_cdf_nohdr (char *file, struct MGD77_CONTROL *F, struct MGD77_DATASET *S)
{
	int err;
    
	MGD77_Select_All_Columns (F, &S->H);

	err = MGD77_Read_Header_Record_nohdr (file, F, &S->H);  /* Will read the entire 24-section header structure */
	if (err) return (err);
    
	err = MGD77_Read_Data_cdf (file, F, S);
	if (err) return (err);
    
	MGD77_nc_status (nc_close (F->nc_id));
    
	return (MGD77_NO_ERROR);
}

int MGD77_Read_Header_Record_nohdr (char *file, struct MGD77_CONTROL *F, struct MGD77_HEADER *H)
{	/* Reads the header structure form a MGD77[+] file */
	int error;
    
	switch (F->format) {
		case MGD77_FORMAT_M77:	/* Will read MGD77 headers from MGD77 files or ascii tables */
			error = MGD77_Read_Header_Record_m77_nohdr (file, F, H);
			break;
		case MGD77_FORMAT_TBL:	/* Will read MGD77 headers from MGD77 files or ascii tables */
			error = MGD77_Read_Header_Record_m77_nohdr (file, F, H);
			break;
		case MGD77_FORMAT_M7T:
			error = MGD77_Read_Header_Record_m77t_nohdr (file, F, H);
			break;
		case MGD77_FORMAT_CDF:	/* Will read MGD77 headers from a netCDF file */
			error = MGD77_Read_Header_Record_cdf (file, F, H);
			break;
		default:
			error = MGD77_UNKNOWN_FORMAT;
			break;
	}
    
	MGD77_Init_Ptr (MGD77_Header_Lookup, H->mgd77);	/* set pointers */
    
	return (error);
}

int MGD77_Read_Header_Record_m77_nohdr (char *file, struct MGD77_CONTROL *F, struct MGD77_HEADER *H)
{	/* Applies to MGD77 files */
	char *MGD77_header[MGD77_N_HEADER_RECORDS], line[BUFSIZ], *not_used = NULL;
	int i, sequence, err, n_eols, c, n;
	struct GMT_STAT buf;
    
	n_eols = c = n = 0;	/* Also shuts up the boring compiler warnings */
    
	/* argument file is generally ignored since file is already open */
    
	memset ((void *)H, '\0', sizeof (struct MGD77_HEADER));	/* Completely wipe existing header */
	if (F->format == MGD77_FORMAT_M77) {			/* Can compute # records from file size because format is fixed */
		if (GMT_STAT (F->path, &buf)) {	/* Inquiry about file failed somehow */
			fprintf (stderr, "%s: Unable to stat file %s\n", GMT_program, F->path);
			GMT_exit (EXIT_FAILURE);
		}
#ifdef WIN32
		/* Count number of records by counting number of new line characters. The non-Windows solution does not work here
		   because the '\r' characters which are present on Win terminated EOLs are apparently stripped by the stdio and
		   so if we cant find their traces (!!!) */
		while ( (c = fgetc( F->fp )) != EOF ) {
			if (c == '\n') n++;
		}
		H->n_records = n;					/* 0 is the header size */
		rewind (F->fp);					/* Go back to beginning of file */
#else
        
		/* Test if we need to use +2 because of \r\n. We could use the above solution but this one looks more (time) efficient. */
		not_used = fgets (line, BUFSIZ, F->fp);
		rewind (F->fp);					/* Go back to beginning of file */
		n_eols = (line[strlen(line)-1] == '\n' && line[strlen(line)-2] == '\r') ? 2 : 1;
		H->n_records = irint ((double)(buf.st_size) / (double)(MGD77_RECORD_LENGTH + n_eols));
#endif
	}
	else {
		/* Since we do not know the number of records, we must quickly count lines */
		while (fgets (line, BUFSIZ, F->fp)) if (line[0] != '#') H->n_records++;	/* Count every line except comments  */
		rewind (F->fp);					/* Go back to beginning of file */
	}
    
	/* Read Sequences No 01-24: */
    
	for (sequence = 0; sequence < MGD77_N_HEADER_RECORDS; sequence++) {
		MGD77_header[sequence] = (char *)GMT_memory (VNULL, (size_t)(MGD77_HEADER_LENGTH + 2), sizeof (char), GMT_program);
/*		if ((err = MGD77_Read_Header_Sequence (F->fp, MGD77_header[sequence], sequence+1))) return (err);*/
	}
	if (F->format != MGD77_FORMAT_M77) not_used = fgets (line, BUFSIZ, F->fp);	/* Skip the column header for tables */
    
	for (i = 0; i < 2; i++) H->mgd77[i] = (struct MGD77_HEADER_PARAMS *) GMT_memory (VNULL, (size_t)1, sizeof (struct MGD77_HEADER_PARAMS), GMT_program);	/* Allocate parameter header */
    
/*	if ((err = MGD77_Decode_Header_m77 (H->mgd77[MGD77_ORIG], MGD77_header, MGD77_FROM_HEADER))) return (err);	 Decode individual items in the text headers */
	for (sequence = 0; sequence < MGD77_N_HEADER_RECORDS; sequence++) GMT_free ((void *)MGD77_header[sequence]);
    
	/* Fill in info in F */
    
	MGD77_set_plain_mgd77 (H, FALSE);				/* Set the info for the standard 27 data fields in MGD-77 files */
	if ((err = MGD77_Order_Columns (F, H))) return (err);	/* Make sure requested columns are OK; if not given set defaults */
    
	return (MGD77_NO_ERROR);	/* Success, it seems */
}

int MGD77_Read_Header_Record_m77t_nohdr (char *file, struct MGD77_CONTROL *F, struct MGD77_HEADER *H)
{	/* Applies to MGD77T files */
	char *MGD77_header, line[BUFSIZ], *not_used = NULL;
	int i, err, n_eols, c, n;
    
	n_eols = c = n = 0;	/* Also shuts up the boring compiler warnings */
    
	/* argument file is generally ignored since file is already open */
    
	memset ((void *)H, '\0', sizeof (struct MGD77_HEADER));	/* Completely wipe existing header */
	/* Since we do not know the number of records, we must quickly count lines */
	while (fgets (line, BUFSIZ, F->fp)) H->n_records++;	/* Count every line */
	rewind (F->fp);					/* Go back to beginning of file */
    
	not_used = fgets (line, BUFSIZ, F->fp);		/* Skip the column header  */
    
	MGD77_header = (char *)GMT_memory (VNULL, (size_t)MGD77T_HEADER_LENGTH, sizeof (char), GMT_program);
	// not_used = fgets (MGD77_header, BUFSIZ, F->fp);			/* Read the entire header record  */
    
	for (i = 0; i < 2; i++) H->mgd77[i] = (struct MGD77_HEADER_PARAMS *) GMT_memory (VNULL, (size_t)1, sizeof (struct MGD77_HEADER_PARAMS), GMT_program);	/* Allocate parameter header */
    
	if ((err = MGD77_Decode_Header_m77t (H->mgd77[MGD77_ORIG], MGD77_header))) return (err);	/* Decode individual items in the text headers */
	GMT_free ((void *)MGD77_header);
    
	/* Fill in info in F */
    
	MGD77_set_plain_mgd77 (H, TRUE);			/* Set the info for the standard 27 data fields in MGD-77 files */
	if ((err = MGD77_Order_Columns (F, H))) return (err);	/* Make sure requested columns are OK; if not given set defaults */
    
	return (MGD77_NO_ERROR);	/* Success, it seems */
}
