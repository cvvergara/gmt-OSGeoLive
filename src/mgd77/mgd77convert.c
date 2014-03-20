/*--------------------------------------------------------------------
 *	$Id: mgd77convert.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *    Copyright (c) 2005-2014 by P. Wessel
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * mgd77convert allows for format conversions between three file formats:
 * a) The standard MGD-77 ASCII punchcard file format from NGDC
 * c) An enhanced "MGD77+" format based on netCDF that allows extra columns
 * t) A plain ASCII table version of the MGD-77 punch cards
 *
 * Input files are sought from both current directory and the list of data
 * directories given in the mgd77_paths.txt file in $MGD77_HOME.  Output is
 * always written to the current directory.  No file will be overwritten
 * unless this is requested.
 *
 * Author:	Paul Wessel
 * Date:	10-MAR-2006
 * Version:	1.1
 *
 */
 
#include "mgd77.h"

int main (int argc, char **argv)
{
	int code_pos, i, in_out, argno, n_cruises = 0, n_paths, format[2] = {MGD77_NOT_SET, MGD77_NOT_SET};
		
	GMT_LONG error = FALSE, force = FALSE, high_resolution = FALSE, original = FALSE;
	
	char file[BUFSIZ], **list = NULL, *fcode = "actm";
	char *format_name[MGD77_N_FORMATS] = {"MGD77 ASCII", "MGD77+ netCDF", "ASCII table", "MGD77T ASCII"};

	struct MGD77_CONTROL M;
	struct MGD77_DATASET *D = NULL;
	
	argc = (int)GMT_begin (argc, argv);		/* Initialize GMT Machinery */
	
	/* Initialize MGD77 output order and other parameters*/
	
	MGD77_Init (&M);			/* Initialize MGD77 Machinery */
	
	for (argno = 1; !error && argno < argc; argno++) {	/* Process input options */
		if (argv[argno][0] != '-') continue;
		
		in_out   = 0;	/* 0 for input and 1 for output (see -T -F) */
		code_pos = 2;	/* argv[arg_no][code_pos] holds the character code for format */
		switch (argv[argno][1]) {

			case 'V':
			case '\0':
				error += GMT_parse_common_options (argv[argno], NULL, NULL, NULL, NULL);
				break;
			
			case 'L':	/* Determine level of error/warning checking and log destination */
				M.verbose_level = 0;
				for (i = 2; argv[argno][i]; i++) {
					if (argv[argno][i] == 'e') M.verbose_level |= 2;
					if (argv[argno][i] == 'w') M.verbose_level |= 1;
					if (argv[argno][i] == '+') M.verbose_dest = 1;
				}
				break;
			case 'T':
				in_out = 1;	/* The "to" field is 1; the fall-through via a missing "break" is intentional */
				if (argv[argno][code_pos] == '+') force = TRUE, code_pos++;	/* Force overwriting existing files */
			case 'F':
				switch (argv[argno][code_pos]) {									
					case 'a':		/* Standard ascii MGD77 file */
						format[in_out] = MGD77_FORMAT_M77;
						break;
					case 'C':		/* Enhanced MGD77+ netCDF file */
						original = TRUE;	/* Overlook revisions */
					case 'c':
						format[in_out] = MGD77_FORMAT_CDF;
						break;
					case 'm':
						format[in_out] = MGD77_FORMAT_M7T;
						break;
					case 't':		/* Plain ascii dat table */
						format[in_out] = MGD77_FORMAT_TBL;
						break;
					default:
						fprintf (stderr, "%s: Option -%c Bad format (%c)!\n", GMT_program, argv[argno][1], argv[argno][code_pos]);
						exit (EXIT_FAILURE);
						break;
				}
				break;
	
			case '4':	/* Selected high-resolution 4-byte integer MGD77+ format for mag, diur, faa, eot [2-byte integer] */
				high_resolution = TRUE;
				break;
			default:		/* Options not recognized */
				error = TRUE;
				fprintf (stderr, "%s: Option -%c not recognized!\n", GMT_program, argv[argno][1]);
				break;
		}
	}
	
	if (error) exit (EXIT_FAILURE);

	/* Check that the options selected are mutually consistent */
	
	if (GMT_give_synopsis_and_exit || argc == 1) {	/* Display usage */
		fprintf(stderr,"mgd77convert %s - Convert MGD77 data to other file formats\n\n", MGD77_VERSION);
		fprintf (stderr, "usage: mgd77convert <cruise(s)> -Fa|c|m|t -T[+]a|c|m|t [-L[e][w][+]] [-V] [-4]\n\n");
         
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
              
		MGD77_Cruise_Explain ();
		fprintf (stderr, "	[Files are read from data repositories and written to current directory]\n");
		fprintf (stderr, "	-F Convert from a file that is either (a) MGD77 ASCII, (c) MGD77+ netCDF, (m) MGD77T ASCII, or (t) plain table\n");
		fprintf (stderr, "	   Use -FC to recover the original MGD77 setting from the MGD77+ file [Default applies E77 corrections]\n");
		fprintf (stderr, "	-T Convert to a file that is either (a) MGD77 ASCII, (c) MGD77+ netCDF, (m) MGD77T ASCII, or (t) plain table\n");
		fprintf (stderr, "	   By default we will refuse to overwrite existing files.  Prepend + to override this policy.\n");
		fprintf (stderr, "	OPTIONS:\n\n");
		fprintf (stderr, "	-L Log level and destination setting for verification reporting.  Append a combination\n");
		fprintf (stderr, "	   of w for warnings, e for errors, and + to send log to stdout [Default is stderr])\n");
		fprintf (stderr, "	-V verbose, report cruise being processed and error summary\n");
		fprintf (stderr, "	-4 Selects high-resolution, 4-byte storage for mag, diur, faa, eot, and msd with precision\n");
		fprintf (stderr, "	   of 10 fTesla, 1 nGal, 0.01 mm [Default is 2-byte with 0.1 nTesla, 0.1 mGal, m precision]\n");
		exit (EXIT_FAILURE);
	}

	n_paths = MGD77_Path_Expand (&M, argv, argc, &list);	/* Get list of requested IDs */

	if (n_paths == 0) {
		fprintf(stderr, "%s: ERROR: No cruises given\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (format[GMT_IN] == MGD77_NOT_SET) {
		fprintf(stderr, "%s: Must specify format of input files\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	
	if (format[GMT_OUT] == MGD77_NOT_SET) {
		fprintf(stderr, "%s: Must specify format of output files\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	
	if (gmtdefs.verbose && format[GMT_IN] == format[GMT_OUT]) {
		fprintf(stderr, "%s: Warning: The two formats chosen are the same\n", GMT_program);
	}
	
	if (format[GMT_OUT] == MGD77_FORMAT_TBL && !(strcmp (gmtdefs.d_format, "%lg") & strcmp (gmtdefs.d_format, "%g"))) {
		strcpy (gmtdefs.d_format, "%.10g");	/* To avoid loosing precision upon rereading this file */
	}
	
	if (format[GMT_OUT] == MGD77_FORMAT_CDF && high_resolution) MGD77_select_high_resolution ();
	
	for (argno = 0; argno < n_paths; argno++) {		/* Process each ID */
	
		D = MGD77_Create_Dataset ();	/* Get data structure w/header */
		MGD77_Reset (&M);		/* Reset to start fresh for next file */

		M.format  = format[GMT_IN];	/* Set input file's format and read everything into memory */
		M.original = original;
		if (original) M.use_corrections[MGD77_M77_SET] = M.use_corrections[MGD77_CDF_SET] = FALSE;	/* Turn off E77 corrections */
		for (i = 0; i < MGD77_N_FORMATS; i++) MGD77_format_allowed[i] = (M.format == i) ? TRUE : FALSE;	/* Only allow the specified input format */
		if (MGD77_Open_File (list[argno], &M, MGD77_READ_MODE)) continue;
		if (MGD77_Read_Header_Record (list[argno], &M, &D->H)) {
			fprintf (stderr, "%s: Error reading header sequence for cruise %s\n", GMT_program, list[argno]);
			exit (EXIT_FAILURE);
		}
		sprintf (file, "%s.%s", M.NGDC_id, MGD77_suffix[format[GMT_OUT]]);
		if (format[GMT_IN] == format[GMT_OUT] && !(M.path[0] == '/' || M.path[1] == ':')) {
			fprintf(stderr, "%s: Input and Output file have same name! Output file will have extension \".new\" appended\n", GMT_program);
			strcat (file, ".new");	/* To avoid overwriting original file */
		}
		if (!access (file, R_OK)) {	/* File exists */
			if (force) {	/* Must delete the file first */
				if (remove (file)) {	/* Oops, removal failed */
					fprintf(stderr, "%s: Unable to remove existing file %s - skipping the conversion\n", GMT_program, file);
					MGD77_Close_File (&M);
					MGD77_Free (D);	/* Free memory allocated by MGD77_Read_File */
					continue;
				}
			}
			else {	/* Cowardly refuse to do this */
				fprintf(stderr, "\n%s: Output file already exists.  Use -T+%c to force overwriting\n", GMT_program, fcode[format[GMT_OUT]]);
				MGD77_Close_File (&M);
				MGD77_Free (D);	/* Free memory allocated by MGD77_Read_File */
				continue;
			}
		}
		
		/* OK, now we can read the data set */
		
		if (MGD77_Read_Data (list[argno], &M, D)) {
			fprintf (stderr, "%s: Error reading data set for cruise %s\n", GMT_program, list[argno]);
			exit (EXIT_FAILURE);
		}
		MGD77_Close_File (&M);
		
		MGD77_Verify_Prep (&M, D);	/* Get key meta-data derived form data records */

		MGD77_Verify_Header (&M, &(D->H), NULL);	/* Verify the header */
	
		if (format[GMT_IN] == MGD77_FORMAT_CDF && format[GMT_OUT] != MGD77_FORMAT_CDF && (D->H.info[MGD77_CDF_SET].n_col || D->flags[0] || D->flags[1])) {
			fprintf(stderr, "\n%s: Warning: Input file contains enhanced material that the output file format cannot represent\n", GMT_program);
		}

		/* OK, ready to write out converted file */
		
		M.format  = format[GMT_OUT];				/* Change the format to the desired output format and write new file in current directory */
		M.original = TRUE;					/* Always write to original attributes */
		for (i = 0; i < MGD77_N_FORMATS; i++) MGD77_format_allowed[i] = (M.format == i) ? TRUE : FALSE;	/* Only allow the specified output format */
		D->H.author = (char *)GMT_memory (VNULL, strlen (M.user)+1, sizeof (char), GMT_program);	/* Allocate space for author */
		strcpy (D->H.author, M.user);									/* Pass current user login id as author */
		if (D->H.history) GMT_free ((void *)D->H.history);	/* Make sure history is blank so it is reset by MGD77_Write_File */
		if (MGD77_Write_File (file, &M, D)) {
			fprintf (stderr, "%s: Error writing new file for cruise %s\n", GMT_program, list[argno]);
			exit (EXIT_FAILURE);
		}
		if (gmtdefs.verbose) {
			fprintf (stderr, "%s: Converted cruise %s to %s format", GMT_program, list[argno], format_name[format[GMT_OUT]]);
			if (D->H.errors[0]) fprintf (stderr, " [%2.2d header problems (%d warnings + %d errors)]", D->H.errors[0], D->H.errors[1], D->H.errors[2]);
			if (D->errors) fprintf (stderr, " [%d data errors]", D->errors);
			fprintf (stderr, "\n");
		}

		MGD77_Free (D);	/* Free memory allocated by MGD77_Read_File */
		n_cruises++;
	}
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Converted %d MGD77 files\n", GMT_program, n_cruises);
	
	MGD77_Path_Free (n_paths, list);
	MGD77_end (&M);

	exit (EXIT_SUCCESS);
}
