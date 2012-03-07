/*-----------------------------------------------------------------
 *	$Id: x2sys_init.c,v 1.46 2011/07/11 19:22:07 guru Exp $
 *
 *      Copyright (c) 1999-2011 by P. Wessel
 *      See LICENSE.TXT file for copying and redistribution conditions.
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; version 2 or any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/
/* x2sys_init will accept command line options to specify track data
 * file formats, region of interest, bin spacing, etc, and a unique
 * system identifier.  It will then initialize the databases associated
 * with this system where information about the tracks and their
 * intersections will be stored.
 *
 * Author:	Paul Wessel
 * Date:	14-JUN-2004
 * Version:	1.1, based on the spirit of the old mgg/s_system code
 *
 */

#include "x2sys.h"
extern void x2sys_set_home (void);

int main (int argc, char **argv)
{
	char *sfile = CNULL, *TAG = CNULL, *suffix = NULL;
	time_t right_now;
#ifndef WIN32
	struct passwd *pw = NULL;
#endif
	char tag_file[BUFSIZ], track_file[BUFSIZ], bin_file[BUFSIZ], def_file[BUFSIZ];
	char path_file[BUFSIZ], path[BUFSIZ], line[BUFSIZ];

	GMT_LONG error = FALSE, geographic = FALSE, force = FALSE;

	FILE *fp = NULL, *fp_def = NULL;

	int n_tags = 0, i, n_found = 0, geodetic = 0, r_entry, i_entry, d_entry, d_start;
	int e_entry,  g_entry, m_entry, wd_entry, wt_entry, n_entry[2], c_entry, n_n = 0;

	double x_min = 0.0, x_max = 0.0, y_min = 0.0, y_max = 0.0, bin_x = 0.0, bin_y = 0.0;

	c_entry = n_entry[0] = n_entry[1] = r_entry = i_entry = d_entry = e_entry = g_entry = m_entry = wd_entry = wt_entry = 0;

	argc = (int)GMT_begin (argc, argv);
	
	for (i = (int)strlen(argv[0]); i >= 0 && argv[0][i] != '/'; i--);
	X2SYS_program = &argv[0][i+1];	/* Name without full path */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'R':
					r_entry = i;
				case 'H':
				case 'V':
				case '\0':
					error += GMT_parse_common_options (argv[i], &x_min, &x_max, &y_min, &y_max);
					break;

				/* Supplemental parameters */

				case 'C':	/* Distance calculation flag */
					if (!strchr ("cefg", (int)argv[i][2])) {
						fprintf(stderr, "%s: ERROR -C: Flag must be c, f, g, or e\n", GMT_program);
						error++;
					}
					c_entry = i;
					break;
				case 'D':
					sfile = &argv[i][2];
					d_entry = i;
					break;
				case 'E':
					suffix = &argv[i][2];
					e_entry = i;
					break;
				case 'G':	/* Geographical coordinates, set discontinuity */
					geographic = TRUE;
					geodetic = 0;
					if (argv[i][2] == 'd') geodetic = 2;
					g_entry = i;
					break;
				case 'F':
					force = TRUE;
					break;
				case 'I':
					if (argv[i][2]) GMT_getinc (&argv[i][2], &bin_x, &bin_y);
					i_entry = i;
					break;
				case 'M':
				case 'm':
					m_entry = i;
					break;
				case 'N':	/* Distance and speed unit selection */
					switch (argv[i][2]) {
						case 'd':	/* Distance unit selection */
						case 's':	/* Speed unit selection */
							if (!strchr ("cekmn", (int)argv[i][3])) {
								fprintf(stderr, "%s: ERROR -N%c: Unit must be c, e, k, m, or n\n", GMT_program, argv[i][2]);
								error++;
							}
							break;
						default:
							fprintf(stderr, "%s: ERROR -N: Choose from -Nd and -Ns\n", GMT_program);
							error++;
							break;
					}
					if (n_n >= 2) {
						fprintf(stderr, "%s: ERROR -N: Can only be given once each for -Nd -Ns\n", GMT_program);
						error++;
					}
					else
						n_entry[n_n++] = i;
					break;
				case 'W':
					switch (argv[i][2]) {
						case 't':	/* Get new timegap */
							wt_entry = i;
							break;
						case 'd':	/* Get new distgap */
							wd_entry = i;
							break;
						default:
							fprintf (stderr, "%s: Syntax Error: -Wt|d<width>\n", GMT_program);
							error++;
							break;
					}
					break;
				default:
					error = TRUE;
					break;
			}
		}
		else {
			TAG = argv[i];
			n_tags++;
		}
	}

	if (argc == 1 || error || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "x2sys_init %s - Initialize a new track system database\n\n", X2SYS_VERSION);
		fprintf (stderr, "usage: x2sys_init <TAG> [-Cc|f|g|e] [-D<deffile>] [-E<suffix>] [-F] [-G[d/g]] [-I[<binsize>]]\n");
		fprintf (stderr, "\t[-N[d|s][c|e|k|m|n]]] [%s] [-V] [-Wt|d|n<gap>] [%s]\n\n", GMT_Rgeo_OPT, GMT_m_OPT);
		fprintf (stderr, "\t<TAG> is the unique system identifier.  Files created will be placed in\n");
		fprintf (stderr, "\t   the directory X2SYS_HOME/<TAG>.\n");

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-C Select procedure for along-track distance and azimuth calculations:\n");
		fprintf (stderr, "\t   c Plain Cartesian\n");
		fprintf (stderr, "\t   f Flat Earth\n");
		fprintf (stderr, "\t   g Great circle [Default]\n");
		fprintf (stderr, "\t   e Ellipsoidal (geodesic) using current ellipsoid\n");
		fprintf (stderr, "\t-D definition file for the track data set [<TAG>.def]\n");
		fprintf (stderr, "\t-E Extension (suffix) for these data files\n");
		fprintf (stderr, "\t   [Default equals the prefix for the definition file]\n");
		fprintf (stderr, "\t-F Force creating new files if old ones are present [Default will abort if old files are found]\n");
		fprintf (stderr, "\t-G for geographical coordinates.  Append g for discontinuity at Greenwich (output 0/360 [Default])\n");
		fprintf (stderr, "\t   and append d for discontinuity at Dateline (output -180/+180)\n");
		fprintf (stderr, "\t-I sets bin size for track bin index output [1/1]\n");
		fprintf (stderr, "\t-N Append (d)istances or (s)peed, and your choice for unit. Choose among:\n");
		fprintf (stderr, "\t   c Cartesian distance (user-dist-units, user user-dist-units/user-time-units)\n");
		fprintf (stderr, "\t   e Metric units I (meters, m/s)\n");
		fprintf (stderr, "\t   k Metric units II (km, km/hr)\n");
		fprintf (stderr, "\t   m British/US units (miles, miles/hr)\n");
		fprintf (stderr, "\t   n Nautical units (nautical miles, knots)\n");
		fprintf (stderr, "\t   [Default is -Ndk -Nse]\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t   [Default region is 0/360/-90/90]\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t-W Set maximum gaps allowed at crossover.  Option may be repeated.\n");
		fprintf (stderr, "\t   -Wt sets maximum time gap (in user units) [Default is infinite]\n");
		fprintf (stderr, "\t   -Wd sets maximum distance gap (in user units) [Default is infinite]\n");
		GMT_explain_option ('m');
		exit (EXIT_FAILURE);
	}

	if (n_tags == 0) {
		fprintf (stderr, "%s: No system tag given!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (n_tags > 1) {
		fprintf (stderr, "%s: Only give one system tag!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (r_entry && (x_min >= x_max || y_min >= y_max)) {
		fprintf (stderr, "%s: -R given inconsistent values!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (i_entry && (bin_x <= 0.0 || bin_y <= 0.0)) {
		fprintf (stderr, "%s: -Idx/dy must be positive!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (!sfile) sfile = TAG;	/* Default def file */
	sprintf (def_file, "%s.def", sfile);
	if (access (def_file, R_OK)) {	/* No such local *.def file */
		fprintf (stderr,"%s: Unable to find local definition file : %s\n", GMT_program, def_file);
		exit (EXIT_FAILURE);
	}
	else if ((fp_def = fopen (def_file, "r")) == NULL) {
		fprintf (stderr,"%s: Unable to open local definition file : %s\n", GMT_program, def_file);
		exit (EXIT_FAILURE);
	}
	for (d_start = (int)strlen (sfile)-1; d_start >= 0 && sfile[d_start] != '/'; d_start--);	/* Find pos of last slash */
	d_start++;		/* Find start of file name */
	
	/* Determine the TAG directory */
	
	x2sys_set_home ();
	x2sys_path (TAG, path);
	if (x2sys_access (TAG, R_OK)) {	/* No such dir */
		if (mkdir (path, (mode_t)0777)) {
			fprintf (stderr,"%s: Unable to create TAG directory : %s\n", GMT_program, path);
			exit (EXIT_FAILURE);
		}
	}
	else if (!force) {	/* Directory exists but -F not on */
		fprintf (stderr,"%s: TAG directory already exists: %s\n", GMT_program, path);
		exit (EXIT_FAILURE);
	}
	
	/* Initialize the system TAG files in X2SYS_HOME/TAG */

	sprintf (tag_file, "%s/%s.tag", TAG, TAG);
	sprintf (def_file, "%s/%s.def", TAG, &sfile[d_start]);
	sprintf (path_file, "%s/%s_paths.txt", TAG, TAG);
	sprintf (track_file, "%s/%s_tracks.d", TAG, TAG);
	sprintf (bin_file, "%s/%s_index.b", TAG, TAG);

	if (!x2sys_access (tag_file, R_OK)) {
		fprintf (stderr,"%s: File exists: %s\n", GMT_program, tag_file);
		x2sys_path (tag_file, path);
		if (force) {
			if (remove (path))
				fprintf (stderr,"%s: Unable to remove %s\n", GMT_program, path);
			else
				fprintf (stderr,"%s: Removed file %s\n", GMT_program, path);
		}
		else 
			n_found++;
	}
	if (!x2sys_access (def_file, R_OK)) {
		fprintf (stderr,"%s: File exists: %s\n", GMT_program, def_file);
		x2sys_path (def_file, path);
		if (force) {
			if (remove (path))
				fprintf (stderr,"%s: Unable to remove %s\n", GMT_program, path);
			else
				fprintf (stderr,"%s: Removed file %s\n", GMT_program, path);
		}
		else 
			n_found++;
	}
	if (!x2sys_access (track_file, R_OK)) {
		fprintf (stderr,"%s: File exists: %s\n", GMT_program, track_file);
		x2sys_path (track_file, path);
		if (force) {
			if (remove (path))
				fprintf (stderr,"%s: Unable to remove %s\n", GMT_program, path);
			else
				fprintf (stderr,"%s: Removed file %s\n", GMT_program, path);
		}
		else 
			n_found++;
	}
	if (!x2sys_access (path_file, R_OK)) {
		fprintf (stderr,"%s: File exists: %s\n", GMT_program, path_file);
		x2sys_path (path_file, path);
		if (force) {
			if (remove (path))
				fprintf (stderr,"%s: Unable to remove %s\n", GMT_program, path);
			else
				fprintf (stderr,"%s: Removed file %s\n", GMT_program, path);
		}
		else 
			n_found++;
	}
	if (!x2sys_access (bin_file, R_OK)) {
		fprintf (stderr,"%s: File exists: %s\n", GMT_program, bin_file);
		x2sys_path (bin_file, path);
		if (force) {
			if (remove (path))
				fprintf (stderr,"%s: Unable to remove %s\n", GMT_program, path);
			else
				fprintf (stderr,"%s: Removed file %s\n", GMT_program, path);
		}
		else 
			n_found++;
	}
	if (n_found) {
		fprintf (stderr,"%s: Remove/rename old files or use -F to overwrite\n", GMT_program);
		exit (EXIT_FAILURE);
	}


	if (gmtdefs.verbose) fprintf (stderr, "%s: Initialize %s\n", GMT_program, tag_file);
	if ((fp = x2sys_fopen (tag_file, "w")) == NULL) {
		fprintf (stderr,"%s: Could not create file %s\n", GMT_program, tag_file);
		exit (EXIT_FAILURE);
	}

        right_now = time ((time_t *)0);
	fprintf (fp, "# TAG file for system: %s\n", TAG);
	fprintf (fp, "#\n# Initialized on: %s", ctime (&right_now));
#ifndef WIN32
	if ((pw = getpwuid (getuid ())) != NULL)
		fprintf (fp, "# Initialized by: %s\n#\n", pw->pw_name);
	else
#endif
		fprintf (fp, "# Initialized by: unknown\n#\n");
	fprintf (fp, "-D%s", &sfile[d_start]);	/* Now a local *.def file in the TAG directory */
	if (c_entry) fprintf (fp, " %s", argv[c_entry]);
	if (e_entry) fprintf (fp, " %s", argv[e_entry]);
	if (g_entry) fprintf (fp, " %s", argv[g_entry]);
	if (m_entry) fprintf (fp, " -m%s", &argv[m_entry][2]);
	if (n_entry[0]) fprintf (fp, " %s", argv[n_entry[0]]);
	if (n_entry[1]) fprintf (fp, " %s", argv[n_entry[1]]);
	if (wt_entry) fprintf (fp, " %s", argv[wt_entry]);
	if (wd_entry) fprintf (fp, " %s", argv[wd_entry]);
	(i_entry) ? fprintf (fp, " %s", argv[i_entry]) : fprintf (fp, " -I1/1");
	(r_entry) ? fprintf (fp, " %s", argv[r_entry]) : fprintf (fp, " -R0/360/-90/90");
	fprintf (fp, "\n");
	x2sys_err_fail (x2sys_fclose (tag_file, fp), tag_file);

	/* Initialize the system's definition file  */

	if (gmtdefs.verbose) fprintf (stderr, "%s: Initialize %s\n", GMT_program, def_file);
	if ((fp = x2sys_fopen (def_file, "w")) == NULL) {
		fprintf (stderr, "%s: Could not create %s\n", GMT_program, def_file);
		exit (EXIT_FAILURE);
	}
	while (fgets (line, BUFSIZ, fp_def)) fprintf (fp, "%s", line);
	x2sys_err_fail (x2sys_fclose (def_file, fp), def_file);
	fclose (fp_def);	/* Close local def file */

	/* Initialize the system's tracks data base  */

	if (gmtdefs.verbose) fprintf (stderr, "%s: Initialize %s\n", GMT_program, track_file);
	if ((fp = x2sys_fopen (track_file, "w")) == NULL) {
		fprintf (stderr, "%s: Could not create %s\n", GMT_program, track_file);
		exit (EXIT_FAILURE);
	}
	x2sys_err_fail (x2sys_fclose (track_file, fp), track_file);

	/* Initialize the system's index data base  */

	if (gmtdefs.verbose) fprintf (stderr, "%s: Initialize %s\n", GMT_program, bin_file);
	if ((fp = x2sys_fopen (bin_file, "wb")) == NULL) {
		fprintf (stderr,"%s: Could not create %s\n", GMT_program, bin_file);
		exit (EXIT_FAILURE);
	}
	x2sys_err_fail (x2sys_fclose (bin_file, fp), bin_file);

	/* Initialize the system's track path file  */

	if (gmtdefs.verbose) fprintf (stderr, "%s: Initialize %s\n", GMT_program, path_file);
	if ((fp = x2sys_fopen (path_file, "wb")) == NULL) {
		fprintf (stderr,"%s: Could not create %s\n", GMT_program, path_file);
		exit (EXIT_FAILURE);
	}
	fprintf (fp, "# Directories with data files for TAG %s\n", TAG);
	fprintf (fp, "# The current directory is always searched first.\n");
	fprintf (fp, "# Add full paths to search additional directories\n");
	x2sys_err_fail (x2sys_fclose (path_file, fp), path_file);
	
	if (gmtdefs.verbose) fprintf (stderr, "%s completed successfully\n", GMT_program);

	GMT_free ((void *)X2SYS_HOME);
	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}
