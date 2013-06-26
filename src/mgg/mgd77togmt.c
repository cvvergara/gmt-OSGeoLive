/*--------------------------------------------------------------------
 *	$Id: mgd77togmt.c 9923 2012-12-18 20:45:53Z pwessel $
 *
 *    Copyright (c) 1991-2013 by P. Wessel and W. H. F. Smith
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * mgd77togmt reads a mgd77 tape file from NGDC and creates
 * a gmt file directly.  If no time is provided, a fake time is
 * issued based on a year and a time increment. We set year = 0 for
 * these cruises and estimate the time_inc based on the distance
 * between 2 nav points and an assumed velocity.
 *
 * Author:	Paul Wessel
 * Date:	9-SEP-1989
 *
 * Modified by WHFS 21 July 1993:
 *	Fixed usage message to be consistent with actual run for -H.
 *	Added carter table correction routines for depth.  Now uses
 *	twt if present, and then if not uses corr depth if present.
 * Modified by PW 19 February 1999:
 *	The 3 binary Carter tables have been replaced by a single ASCII
 *	table carter.d which makes it simpler to port this program.
 * Modified by PW 21 June 1999:
 *	The new MGD77 format is now Y2K compliant.  Added support for
 *	the new format, plus Y2K kludge for older files based on the
 *	NGDC-provided fact that the oldest cruise in the database is
 *	from 1939.  Hence, all 2-digit years from 0 to 38 are inter-
 *	preted to mean 20xx.  This fix will fail in 2039 if there will
 *	be fools then that make old-style MGD77 files.
 *	Also, removed the effect of the -H switch since we can check
 *	the first column for what type of record it is.  Left -H in
 * 	for backward compatibility.
 *	If no -A is given, agency is set equal to the legname.
 *	Placed all the Carter table initialization parameters in a
 *	separate include file (carter.h).
 */

#include "gmt.h"
#include "gmt_mgg.h"

#define MPRDEG 111.1949e-3

int main (int argc, char **argv) {
	int n_records, *year = NULL, k;
	int i, rec, n_read, n_files = 0, n_alloc = GMT_CHUNK, leg_year = 0, len;
	int T_INC = 60, fake_running_time = 0;
	int t_flag = FALSE, anom_offset = 0;
	GMT_LONG use_list = FALSE, set_agency = FALSE, mag_rewind = FALSE, greenwich = FALSE, give_synopsis_and_exit = FALSE;
	GMT_LONG error = FALSE;
	double cable_len = 0;
	char file[80], *mfile = NULL, *list = NULL, agency[10], *legfname = NULL, line[BUFSIZ], **mgd77 = NULL, **prefix = NULL;
	struct GMTMGG_TIME *gmt = NULL;
	struct GMTMGG_REC *record = NULL;
	FILE *fp = NULL, *fpo = NULL;

	argc = (int)GMT_begin (argc, argv);

	gmtmggpath_init();

	legfname = list = mfile = CNULL;
	memset ((void *)agency, 0, (size_t)10);
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case '\0':
					error = give_synopsis_and_exit = TRUE;
					break;
				case 'Y':
					leg_year = atoi (&argv[i][2]);
					break;
				case 'I':	/* No time in data file, specify time increment between points */
					T_INC = atoi (&argv[i][2]);
					t_flag = TRUE;	/* fake time starts at 1/1/2000 @ 00:00:00 hr */
					break;
				case 'H':	/* No longer used, kept for nostalgic reasons */
					break;
				case 'F':
					legfname = &argv[i][2];
					break;
				case 'L':
					list = &argv[i][2];
					use_list = TRUE;
					break;
				case 'A':
					strcpy (agency, &argv[i][2]);
					break;
				case 'T':			/* Read Total field instead of anomaly and subtract a cte */
					if (argv[i][2])
						anom_offset = atoi (&argv[i][2]);
					else
						anom_offset = 40000;
					break;
				case 'G':
					greenwich = TRUE;
					break;
				case 'V':
					gmtdefs.verbose = TRUE;
					break;
				case 'W':
					if (argv[i][2])
						cable_len = atof (&argv[i][2]);
					mag_rewind = TRUE;
					break;
				default:
					fprintf (stderr, "SYNTAX ERROR:  Unrecognized option -%c\n", argv[i][1]);
					error = TRUE;
					break;
			}
		}
		else
			mfile = argv[i];
	}
	
	if (!give_synopsis_and_exit) {
		if (use_list && !list) {
			fprintf (stderr, "SYNTAX ERROR -L option:  Must specify list file\n");
			error = TRUE;
		}
		if (use_list && (legfname || leg_year != 0)) {
			fprintf (stderr, "SYNTAX ERROR -L option:  Specify -L or the combination -F, -Y\n");
			error = TRUE;
		}
		if (!use_list && !mfile && !legfname) {
			fprintf (stderr, "SYNTAX ERROR -F option:  When using standard input you must use -F option.\n");
			error = TRUE;
		}
	}
	
	if (argc == 1 || error) {
		fprintf (stderr, "usage: mgd77togmt [mgd77file] [-F<filename>] -Y<leg_year> OR -L<leglist> [-A<10 char agency name>]\n");
		fprintf (stderr, "\t[-G] [NGDC-file -I<time_increment>] [-V] [-T[<offset>]] [-W[<cable_length>]]\n\n");
		
		if (give_synopsis_and_exit) exit (EXIT_FAILURE);
		
		fprintf (stderr, "\t-Y sets start year. If not provided and -L option not used, it tries to get\n");
		fprintf (stderr, "\t   it from header file. The header file must be in the same directory of the\n");
		fprintf (stderr, "\t   main file and must have a name equal to the main but with a .h77 extension.\n");
		fprintf (stderr, "\t-L gives name of a list with several records of <mgd77file> <gmtprefix> <leg_year>\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A sets optional 10 char info for gmt header.\n");
		fprintf (stderr, "\t-F sets gmtfilename prefix (e.g., without path or .gmt). If not given, it\n");
		fprintf (stderr, "\t   will be constructed from the mgd77file name plus the .gmt extension.\n");
		fprintf (stderr, "\t-G force geographical longitudes (-180/+180) [Default is 0-360]\n");
		fprintf (stderr, "\t-I sets fake timeincrement for files without time information\n");
		fprintf (stderr, "\t-T Extracts Total field instead of anomaly. Since F does not hold in a 2 byte int var\n");
		fprintf (stderr, "\t   we subtract a constant [default = 40000] but you can provide another value in <offset>.\n");
		fprintf (stderr, "\t-W Take into account that the magnetometer is not at ship's position.\n");
		fprintf (stderr, "\t   <cable> is magnetometer tow distance [default = 200 meters].\n");
		fprintf (stderr, "\t   If -W only is given (e.g., no <cable_length>) and like with the -Y option, we try\n");
		fprintf (stderr, "\t   to get the tow distance from the header file. Failing, defaults to 200 meters.\n");
		fprintf (stderr, "\t   Note that this option will throw away the first points whose accumulated.\n");
		fprintf (stderr, "\t   distance since the start of magnetic acquisition is less than cable.\n");
		fprintf (stderr, "\t   length, and likewise for the end of the mag profile.\n");
		fprintf (stderr, "\t-V runs in verbose mode.\n");
		
		exit (EXIT_FAILURE);
	}

	if (agency[0]) set_agency = TRUE;

	if (use_list) {
	
		if ((fp = fopen (list, "r")) == NULL) {
			fprintf (stderr, "mgd77togmt: Cannot open file %s!\n", list);
			exit (EXIT_FAILURE);
		}
		mgd77 = (char **) GMT_memory (VNULL, (size_t)n_alloc, sizeof (char *), "mgd77togmt");
		prefix = (char **) GMT_memory (VNULL, (size_t)n_alloc, sizeof (char *), "mgd77togmt");
		year = (int *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (int), "mgd77togmt");
		n_read = 0;
		while (fgets (line, BUFSIZ, fp)) {
			mgd77[n_read]  = (char *) GMT_memory (VNULL, (size_t)80, sizeof (char), "mgd77togmt");
			prefix[n_read] = (char *) GMT_memory (VNULL, (size_t)80, sizeof (char), "mgd77togmt");
			if ((sscanf (line, "%s %s %d", mgd77[n_read], prefix[n_read], &year[n_read]) != (size_t)3)) {
				fprintf (stderr, "mgd77togmt: Trouble reading record # %d in list %s\n", n_read, list);
				exit (EXIT_FAILURE);
			}
			n_read++;
			if (n_read == n_alloc) {
				n_alloc <<= 1;
				mgd77 = (char **) GMT_memory ((void *)mgd77, (size_t)n_alloc, sizeof (char *), "mgd77togmt");
				prefix = (char **) GMT_memory ((void *)prefix, (size_t)n_alloc, sizeof (char *), "mgd77togmt");
				year = (int *) GMT_memory ((void *)year, (size_t)n_alloc, sizeof (int), "mgd77togmt");
			}
		}
		n_files = n_read;
		fclose (fp);
		mgd77  = (char **) GMT_memory ((void *)mgd77, (size_t)n_files, sizeof (char *), "mgd77togmt");
		prefix = (char **) GMT_memory ((void *)prefix, (size_t)n_files, sizeof (char *), "mgd77togmt");
		year = (int *) GMT_memory ((void *)year, (size_t)n_files, sizeof (int), "mgd77togmt");
	}
	else {
		int	free_legfname = 0;
		char	*ext;
		if (!legfname) {
			legfname = strdup(mfile);
			free_legfname = 1;
			GMT_chop_ext(legfname);
		}
		else {		/* Make sure that the legname does not have the .gmt extension */
			ext = GMT_chop_ext(legfname);
			if (ext && strcmp(ext,".gmt"))		/* Hoops, removed a wrong extension */
				strcat(legfname,ext);
		}
		mgd77  = (char **) GMT_memory (VNULL, (size_t)1, sizeof (char *), "mgd77togmt");
		prefix = (char **) GMT_memory (VNULL, (size_t)1, sizeof (char *), "mgd77togmt");
		year = (int *) GMT_memory (VNULL, (size_t)1, sizeof (int), "mgd77togmt");
		mgd77[0]  = (char *) GMT_memory (VNULL, (size_t)80, sizeof (char), "mgd77togmt");
		prefix[0] = (char *) GMT_memory (VNULL, (size_t)80, sizeof (char), "mgd77togmt");
		strcpy (mgd77[0], mfile);
		strcpy (prefix[0], legfname);
		year[0] = leg_year;
		n_files = 1;
		if (free_legfname) free ((void *)legfname);
	}

	n_alloc = GMT_CHUNK;

	for (i = 0; i < n_files; i++) {
	
		if (gmtdefs.verbose) fprintf (stderr, "mgd77togmt: Processing file %s\n", mgd77[i]);
		
		if (mgd77[i][0] == 0) {
			fp = stdin;
		}
		else if ((fp = fopen (mgd77[i], "r")) == NULL) {
			fprintf (stderr, "mgd77togmt: Unable to open file %s - skipping\n", mgd77[i]);
			GMT_free ((void *)mgd77[i]);
			GMT_free ((void *)prefix[i]);
			continue;
		}

		if (!set_agency) {	/* Use legname as agency string */
			strncpy ((void *)agency, prefix[i], (size_t)10);
			if (gmtdefs.verbose) fprintf (stderr, "mgd77togmt: Agency string set to %s\n", prefix[i]);
		}

		/* See if we have the header file arround and if yes try to get year and mag tow distance (if needed) */
		if (!leg_year || (mag_rewind && !cable_len)) {
			char	*hdr, s_yr[5];
			FILE	*fph;
			hdr = strdup(&mgd77[i][0]);
			GMT_chop_ext(hdr);
			strcat(hdr,".h77");
			k = 0;
			if ((fph = fopen (hdr, "r")) != NULL) {
				while (fgets (line, BUFSIZ, fph) && k < 13) {
					k++;
					if (k == 4 && !leg_year) {
						strncpy (s_yr, &line[0], (size_t) 4);	s_yr[4] = 0;
						year[i] = atoi (s_yr);
						if (gmtdefs.verbose) fprintf (stderr, "\tGot year %d from header file %s\n", year[i], hdr);
						if (!mag_rewind || cable_len > 0) k = 14;	/* End reading header */
					}
					else if (k == 13) {	/* Seek for the magnetometer tow distance */
						strncpy (s_yr, &line[5], (size_t) 4);	s_yr[4] = 0;
						if (s_yr[0] == ' ' && s_yr[1] == ' ' && s_yr[2] == ' ' && s_yr[3] == ' ')
							cable_len = 200;	/* The default tow distance */
						else
							cable_len = atof (s_yr);
						if (cable_len > 500) cable_len = 200;	/* This accounts also for 999 non-values in header */
						if (gmtdefs.verbose) fprintf (stderr, "\tUsing magnetometer tow distance = %.0f\n", cable_len);
					}
				}
				fclose(fph);
			}
			else {
				if (!leg_year) {
					fprintf (stderr, "WARNING: Year not provided and companion header file not found. Jumping this file\n");
					continue;
				}
				if (mag_rewind) cable_len = 200;
			}

			free ((void *) hdr);
		}

		gmt = gmtmgg_init (year[i]);

		record = (struct GMTMGG_REC *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (struct GMTMGG_REC), "mgd77togmt");

		rec = n_read = 0;
		while (fgets (line, BUFSIZ, fp)) {
			n_read++;
			if (!(line[0] == '3' || line[0] == '5')) continue;	/* Only process data records */
			if ((len = (int)strlen(line)) != 121) {
				fprintf (stderr, "mgd77togmt: Record # %d has incorrect length (%d), skipped\n", n_read, len);
				continue;
			}
			if (!gmtmgg_decode_MGD77 (line, t_flag, &record[rec], &gmt, anom_offset)) {
				if (t_flag) record[rec].time = (fake_running_time += T_INC);
				if (greenwich && record[rec].lon > 180000000) record[rec].lon -= 360000000;
				rec++;
			}
			else
				fprintf (stderr, "mgd77togmt: Trouble decoding record # %d (skipped)\n", n_read);


			if (rec == n_alloc) {
				n_alloc <<= 1;
				record = (struct GMTMGG_REC *) GMT_memory ((void *)record, (size_t)n_alloc, sizeof (struct GMTMGG_REC), "mgd77togmt");
			}
		}
		if (fp != stdin) fclose (fp);
		GMT_free ((void *)gmt);
		n_records = rec;
	
		sprintf (file, "%s.gmt", prefix[i]);
		if ((fpo = fopen (file, "wb")) == NULL) {
			fprintf (stderr, "mgd77togmt: Could not create file %s\n", file);
			exit (EXIT_FAILURE);
		}
	
		if (fwrite ((void *)(&year[i]), (size_t)4, (size_t)1, fpo) != (size_t)1) {
			fprintf (stderr,"mgd77togmt: Error while writing first year\n");
			exit (EXIT_FAILURE);
		}
		if (fwrite ((void *)(&n_records), (size_t)4, (size_t)1, fpo) != (size_t)1) {
			fprintf (stderr,"mgd77togmt: Error while writing no of records\n");
			exit (EXIT_FAILURE);
		}
		if (fwrite ((void *)agency, (size_t)10, (size_t)1, fpo) != (size_t)1) {

			fprintf (stderr,"mgd77togmt: Error while writing info-header\n");
			exit (EXIT_FAILURE);
		}

		if (!mag_rewind) {
			for (rec = 0; rec < n_records; rec++) {
				if (fwrite ((void *)(&record[rec]), (size_t)18, (size_t)1, fpo) != (size_t)1) {
					fprintf (stderr,"mgd77togmt: Error writing data record no %d\n",rec);
					exit (EXIT_FAILURE);
				}
			}
		}
		else {			/* At ship's position we'll get the mag reading down the track that is closest to cable length */
			int	n, dlon, last_lon = 0, last_lat = 0, itmp;
			double	*ds, dds, dx, dy;
			ds = (double *) GMT_memory (VNULL, (size_t)n_records, sizeof (double), "mgd77togmt");
			for (rec = 0; rec < n_records; rec++) {
				if (rec == 0) {
					last_lon = record[0].lon;
					last_lat = record[0].lat;
					ds[0] = 0.0;
				}
				else {
					dlon = record[rec].lon - last_lon;
					dx = (double) dlon * cosd (0.5e-06*(double)(record[rec].lat+last_lat));
					dy = (double) (record[rec].lat - last_lat);
					ds[rec] = ds[rec-1] + MPRDEG * hypot (dx, dy);
					last_lon = record[rec].lon;
					last_lat = record[rec].lat;
				}
			}

			for (rec = 0; rec < n_records; rec++) {
				dds = ds[rec] - cable_len;
				n = rec;
				if (dds < 0) {			/* First points (of distance < cable_len) are lost */
					record[rec].gmt[1] = GMTMGG_NODATA;
				}
				else {
					while ((ds[n] - dds) > 0) n--;
				}
				itmp = record[rec].gmt[1];
				record[rec].gmt[1] = record[n].gmt[1];
				if (fwrite ((void *)(&record[rec]), (size_t)18, (size_t)1, fpo) != (size_t)1) {
					fprintf (stderr,"mgd77togmt: Error writing data record no %d\n",rec);
					exit (EXIT_FAILURE);
				}
				record[rec].gmt[1] = itmp;	/* Reset to original to be used when its turn arrives */
			}
		}
		fclose (fpo);
	
		GMT_free ((void *)record);
		
		GMT_free ((void *)mgd77[i]);
		GMT_free ((void *)prefix[i]);
	}
	
	GMT_free ((void *)mgd77);
	GMT_free ((void *)prefix);
	GMT_free ((void *)year);

	gmtmgg_end ();
	GMT_end (argc, argv);
	
	exit (EXIT_SUCCESS);
}
