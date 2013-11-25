/*--------------------------------------------------------------------
 *	$Id: mgd77list.c 10110 2013-10-31 22:02:27Z pwessel $
 *
 *    Copyright (c) 2004-2013 by P. Wessel
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * mgd77list produces ASCII listings of <ngdc-id>.mgd77 files. The *.mgd77
 * files distributed from NGDC contain along-track geophysical observations
 * such as bathymetry, gravity, and magnetics, and the user may extract
 * any combination of these parameters as well as four generated quantities
 * such as distance (in km), heading, velocity (m/s), and weight by using
 * the -F option.  The order of the choices given to the -F option is used
 * to determine the sequence in which the parameters will be printed out.
 * If -F is not specified, the default will output all the data.  E.g. to
 * create an input file for surface, use -Flon,lat,depth (for bathymetry).
 *
 * To select a sub-section of the track, specify the start/endpoints by:
 *	1) Start-time (yyyy-mm-ddT[hh:mm:ss]) OR start-distance (km)
 *	2) Stop-time  (yyyy-mm-ddT[hh:mm:ss]) OR stop-distance (km)
 * To select data inside an area, use the -R option.
 * To start output with a header string, use -H.
 * To separate each data set with a multisegment header string, use -M.
 *
 * Author:	Paul Wessel
 * Date:	19-JUN-2004
 * Version:	1.0 Based somewhat on the old gmtlist.c
 *		31-MAR-2006: Changed -X to -L to avoid GMT collision
 *		23-MAY-2006: Added -Q for limits on speed/azimuths
 *		21-FEB-2008: Added -Ga|b<rec> for limits on rec range
 *		08-MAR-2012: Deal with MGD77T format files and add more aux
 *				columns to recreate orig format (date, tz, hhmm, dmin)
 *
 */
 
#include "mgd77.h"

#define MGD77_FMT  "drt,id,tz,year,month,day,hour,dmin,lat,lon,ptc,twt,depth,bcc,btc,mtf1,mtf2,mag,msens,diur,msd,gobs,eot,faa,nqc,sln,sspn"
#define MGD77_ALL  "drt,id,time,lat,lon,ptc,twt,depth,bcc,btc,mtf1,mtf2,mag,msens,diur,msd,gobs,eot,faa,nqc,sln,sspn"
#define MGD77T_FMT "id,tz,date,hhmm,lat,lon,ptc,nqc,twt,depth,bcc,btc,bqc,mtf1,mtf2,mag,msens,diur,msd,mqc,gobs,eot,faa,gqc,sln,sspn"
#define MGD77T_ALL "id,time,lat,lon,ptc,nqc,twt,depth,bcc,btc,bqc,mtf1,mtf2,mag,msens,diur,msd,mqc,gobs,eot,faa,gqc,sln,sspn"
#define MGD77_GEO  "time,lat,lon,twt,depth,mtf1,mtf2,mag,gobs,faa"
#define MGD77_AUX  "dist,azim,vel,weight"

#define ADJ_CT	0
#define ADJ_DP	1
#define ADJ_GR	2
#define ADJ_MG	3

int main (int argc, char **argv)
{
	int i, c, id, k, kx, pos, argno, code, n_cruises = 0, n_paths, use, n_items = 0, t_pos = MGD77_NOT_SET;
	int t_col, x_col, y_col, z_col, e_col = 0, m_col = 0, f_col = 0, g_col = 0, m1_col = 0, m2_col = 0, ms_col = 0, twt_col = 0;
	int n_sub, n_out_columns, n_cols_to_process, n_aux, select_option, adj_code[4];
	int time_column, lon_column, lat_column, dist_flag = 2, GF_version = MGD77_NOT_SET;
	
	GMT_LONG rec, prevrec, n_out = 0, start_rec = 0, stop_rec;
	
	GMT_LONG error = FALSE, exact = FALSE, string_output = FALSE, limit_institution = FALSE, need_depth = FALSE, fake_times = FALSE;
	GMT_LONG negative_depth = FALSE,  negative_msd = FALSE, need_distances, need_time, limit_on_time = FALSE, limit_on_recs = FALSE;
	GMT_LONG limit_on_dist = FALSE, limit_on_wesn = FALSE, need_lonlat = FALSE, first_cruise = TRUE, need_twt = FALSE, this_limit_on_time;
	GMT_LONG need_date, need_sound = FALSE, apply_corrections = FALSE, adj_force = FALSE, PDR_wrap, skip_if_time_is_NaN = FALSE;
	GMT_LONG limit_on_azimuth = FALSE, limit_on_velocity = FALSE, lonlat_not_NaN, first_warning = TRUE, has_prev_twt = FALSE;
	
	char *start_date = NULL, *stop_date = NULL, f_setting[BUFSIZ], fx_setting[BUFSIZ], **list = NULL, **item_names = NULL;
	char *tvalue[MGD77_MAX_COLS], institution[3], s_unit[2], d_unit[2], *aux_tvalue[N_MGD77_AUX], *correction_table = NULL;
	char buffer[BUFSIZ];
	
	double west, east, south, north, start_time, stop_time, start_dist, IGRF[7], correction, min_vel = 0.0, max_vel = DBL_MAX;
	double stop_dist, dist_scale, vel_scale, ds, dt, cumulative_dist, aux_dvalue[N_MGD77_AUX], *out = NULL;
	double sound_speed = 0.0, i_sound_speed = 0.0, date = 0.0, g, m, z, v, twt, *dvalue[MGD77_MAX_COLS];
	double min_az = 0.0, max_az = 360.0, prev_twt = 0, d_twt, twt_pdrwrap_corr;
	
	struct MGD77_CONTROL M;
	struct MGD77_DATASET *D = NULL;
	struct MGD77_AUX_INFO aux[N_MGD77_AUX];
	struct GMT_gcal cal;
	struct MGD77_CARTER Carter;
	struct MGD77_CORRTABLE **CORR = NULL;
	struct MGD77_AUXLIST auxlist[N_MGD77_AUX] = {
		{ "dist",    MGD77_AUX_DS, 0, 0, "d(km)"},
		{ "azim",    MGD77_AUX_AZ, 0, 0, "azimuth"},
		{ "vel",     MGD77_AUX_SP, 0, 0, "v(m/s)"},
		{ "year",    MGD77_AUX_YR, 0, 0, "year"},
		{ "month",   MGD77_AUX_MO, 0, 0, "month"},
		{ "day",     MGD77_AUX_DY, 0, 0, "day"},
		{ "hour",    MGD77_AUX_HR, 0, 0, "hour"},
		{ "min",     MGD77_AUX_MI, 0, 0, "minute"},
		{ "dmin",    MGD77_AUX_DM, 0, 0, "dec-minute"},
		{ "sec",     MGD77_AUX_SC, 0, 0, "second"},
		{ "date",    MGD77_AUX_DA, 1, 0, "date"},
		{ "hhmm",    MGD77_AUX_HM, 0, 0, "hourmin"},
		{ "weight",  MGD77_AUX_WT, 0, 0, "weight"},
		{ "drt",     MGD77_AUX_RT, 0, 0, "rectype"},
		{ "igrf",    MGD77_AUX_MG, 0, 0, "IGRF"},
		{ "carter",  MGD77_AUX_CT, 0, 0, "Carter"},
		{ "ngrav",   MGD77_AUX_GR, 0, 0, "IGF"},
		{ "ngdcid",  MGD77_AUX_ID, 1, 0, "NGDC-ID"}
#ifdef USE_CM4
	, { "cm4",  MGD77_AUX_CM, 0, 0, "CM4"}
#endif
	};
#ifdef USE_CM4
	struct MGD77_CM4 CM4;
#endif
	PFD GMT_azimuth_func;
	
	GMT_LONG separate_aux_columns (struct MGD77_CONTROL *F, char *fx_setting, struct MGD77_AUX_INFO *aux, struct MGD77_AUXLIST *A);
	GMT_LONG augment_aux_columns (int n_items, char **item_name, struct MGD77_AUX_INFO *aux, struct MGD77_AUXLIST *auxlist, int n_aux);

	argc = (int)GMT_begin (argc, argv);		/* Initialize GMT Machinery */
	
	/* Initialize MGD77 output order and other parameters*/
	
	MGD77_Init (&M);			/* Initialize MGD77 Machinery */
	west = 0.0;	east = 360.0;	south = -90.0;	north = 90.0;
	start_time = start_dist = 0.0;	/* Default covers all times and distances */
	stop_time = stop_dist = DBL_MAX;
	stop_rec = INT_MAX;
	aux_dvalue[MGD77_AUX_WT] = 1.0;	/* Default weight */
	aux_dvalue[MGD77_AUX_RT] = 5.0;	/* Default record type */
	memset ((void *)adj_code, 0, 4*sizeof (int));
	memset ((void *)institution, 0, (size_t)3);
	d_unit[0] = 'k';	d_unit[1] = '\0';	/* Default is -Ndk */
	s_unit[0] = 'e';	s_unit[1] = '\0';	/* Default is -Nse */
	strcpy (f_setting, MGD77_ALL);	/* Default is a full MGD77 record */

	for (i = 1; !error && i < argc; i++) {	/* Process input options */
		if (argv[i][0] != '-') continue;
		
		switch(argv[i][1]) {

			case 'H':
			case 'M':
			case 'R':
			case 'V':
			case 'b':
			case 'm':
			case '\0':
				error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
				break;
						
			case 'A':	/* Adjustment flags */
				k = 2;
				if (argv[i][2] == '+') {	/* Recalculate anomalies even if original anomaly == NaN [Default leaves NaNs unchanged] */
					adj_force = TRUE;
					k = 3;
				}
				switch (argv[i][k]) {
					case 'c':	/* Carter correction adjustment */
						code = argv[i][k+1] - '0';
						if (code < 1 || code > 11) {
							fprintf (stderr, "%s: ERROR -Ac<code>.  <code> must be 1,2,4,8 or binary combination.\n", GMT_program);
							error++;
						}
						if (argv[i][k+2] == ',') {
							sound_speed = atof (&argv[i][k+3]);
							if (sound_speed < 1400.0 || sound_speed > 1600.0) {
								fprintf (stderr, "%s: ERROR -Ac<code>,<speed>.  <speed> in m/s in the 1400-1600 range.\n", GMT_program);
							error++;
							}
						}
						adj_code[ADJ_CT] |= code;
						break;
					case 'd':	/* depth adjustment */
						code = argv[i][k+1] - '0';
						if (code < 1 || code > 7) {
							fprintf (stderr, "%s: ERROR -Ad<code>.  <code> must be 1,2,4 or binary combination.\n", GMT_program);
							error++;
						}
						if (argv[i][k+2] == ',') {
							sound_speed = atof (&argv[i][k+3]);
							if (sound_speed < 1400.0 || sound_speed > 1600.0) {
								fprintf (stderr, "%s: ERROR -Ad<code>,<speed>.  <speed> in m/s in the 1400-1600 range.\n", GMT_program);
							error++;
							}
						}
						adj_code[ADJ_DP] |= code;
						break;
					case 'f':	/* faa adjustment */
						code = argv[i][k+1] - '0';
						if (code < 1 || code > 3) {
							fprintf (stderr, "%s: ERROR -Af<code>.  <code> must be 1-3.\n", GMT_program);
							error++;
						}
						if (argv[i][k+2] == ',') {
							GF_version = atoi (&argv[i][k+3]);
							if (GF_version < MGD77_IGF_HEISKANEN || GF_version > MGD77_IGF_1980) {
								fprintf (stderr, "%s: ERROR -Af<code>,<field>.  Select <field> is 1-4 range.\n", GMT_program);
							error++;
							}
						}
						adj_code[ADJ_GR] |= code;
						break;
					case 'm':	/* mag adjustment */
						code = atoi (&argv[i][k+1]);
						if (code < 1 || code > 31) {
							fprintf (stderr, "%s: ERROR -Am<code>.  <code> must be 1,2,4 or binary combination.\n", GMT_program);
							error++;
						}
						adj_code[ADJ_MG] |= code;
						break;
					case 't':	/* fake time requires */
						fake_times = TRUE;
						break;
					default:
						fprintf (stderr, "%s: ERROR -A<flag>.  <flag> must be c, d, g, m, or t.\n", GMT_program);
						error++;
						break;
				}
				break;
			
			case 'C':	/* Distance calculation flag */
				if (argv[i][2] == 'f') dist_flag = 1;
				if (argv[i][2] == 'g') dist_flag = 2;
				if (argv[i][2] == 'e') dist_flag = 3;
				if (dist_flag < 1 || dist_flag > 3) {
					fprintf(stderr, "%s: ERROR -C: Flag must be f, g, or e\n", GMT_program);
					error++;
				}
				break;

			case 'D':		/* Assign start/stop times for sub-section */
				limit_on_time = TRUE;
				switch (argv[i][2]) {
				 	case 'A':		/* Start date, skip records with time = NaN */
						skip_if_time_is_NaN = TRUE;
				 	case 'a':		/* Start date */
						start_date = &argv[i][3];
						break;
					case 'B':		/* Stop date, skip records with time = NaN */
						skip_if_time_is_NaN = TRUE;
					case 'b':		/* Stop date */
						stop_date = &argv[i][3];
						break;
					default:
						error = TRUE;
						break;
				}
				break;

			case 'E':	/* Exact parameter match */
				exact = TRUE;
				break;
	
			case 'F':	/* Selected output fields */
				strcpy (f_setting, &argv[i][2]);
				if (!strcmp (f_setting, "mgd77")) strcpy (f_setting, MGD77_FMT);
				if (!strcmp (f_setting, "mgd77+")) {
					strcpy (f_setting, MGD77_FMT);
					strcat (f_setting, ",");
					strcat (f_setting, MGD77_AUX);
				}
				if (!strcmp (f_setting, "mgd77t")) strcpy (f_setting, MGD77T_FMT);
				if (!strcmp (f_setting, "mgd77t+")) {
					strcpy (f_setting, MGD77T_FMT);
					strcat (f_setting, ",");
					strcat (f_setting, MGD77_AUX);
				}
				if (!strcmp (f_setting, "all")) strcpy (f_setting, MGD77_ALL);
				if (!strcmp (f_setting, "all+")) {
					strcpy (f_setting, MGD77_ALL);
					strcat (f_setting, ",");
					strcat (f_setting, MGD77_AUX);
				}
				if (!strcmp (f_setting, "allt")) strcpy (f_setting, MGD77T_ALL);
				if (!strcmp (f_setting, "allt+")) {
					strcpy (f_setting, MGD77T_ALL);
					strcat (f_setting, ",");
					strcat (f_setting, MGD77_AUX);
				}
				if (!strcmp (f_setting, "geo")) strcpy (f_setting, MGD77_GEO);
				if (!strcmp (f_setting, "geo+")) {
					strcpy (f_setting, MGD77_GEO);
					strcat (f_setting, ",");
					strcat (f_setting, MGD77_AUX);
				}
				break;
	
			case 'G':		/* Assign start/stop records for sub-section */
				limit_on_recs = TRUE;
				if (argv[i][2] == 'a') {		/* Start rec */
					start_rec = atoi (&argv[i][3]);
				}
				else if (argv[i][2] == 'b')	 {	/* Stop rec */
					stop_rec = atoi (&argv[i][3]);
				}
				else
					error = TRUE;
				break;
	
			case 'I':
				MGD77_Process_Ignore (argv[i][1], &argv[i][2]);
				break;

			case 'L':	/* Crossover correction table */
				correction_table = &argv[i][2];
				apply_corrections = TRUE;
				break;

			case 'N':	/* Nautical units (knots, nautical miles) */
				switch (argv[i][2]) {
					case 'd':	/* Distance unit selection */
					case 'D':
						d_unit[0] = argv[i][3];
						if (!strchr ("ekmn", (int)d_unit[0])) {
							fprintf(stderr, "%s: ERROR -Nd: Unit must be e, k, m, or n\n", GMT_program);
							error++;
						}
						break;
					case 's':	/* Speed unit selection */
					case 'S':
						s_unit[0] = argv[i][3];
						if (!strchr ("ekmn", (int)s_unit[0])) {
							fprintf(stderr, "%s: ERROR -Nd: Unit must be e, k, m, or n\n", GMT_program);
							error++;
						}
						break;
				}
				break;

			case 'P':	/* Restrict output to this institution */
				/* 01/31/2006: OBSOLETE - Left for backwards compatibility */
				limit_institution = TRUE;
				institution[0] = argv[i][2];
				institution[1] = argv[i][3];
				break;
			case 'Q':		/* Assign min/max values for speeds or azimuth */
				switch (argv[i][2]) {
					case 'a':	/* Azimuth min/max */
						if (sscanf (&argv[i][3], "%lf/%lf", &min_az, &max_az) != 2) {
							fprintf(stderr, "%s: ERROR -Qa: append min/max azimuth limits [0/360]\n", GMT_program);
							error++;
						}
						limit_on_azimuth = TRUE;
						break;
					case 'v':	/* Velocity min/max */
						code = sscanf (&argv[i][3], "%lf/%lf", &min_vel, &max_vel);
						if (code == 1)
							max_vel = DBL_MAX;
						else if (code <= 0) {
							fprintf(stderr, "%s: ERROR -Qv: append min[/max] velocity limits [0]\n", GMT_program);
							error++;
						}
						limit_on_velocity = TRUE;
						break;
					default:
						fprintf(stderr, "%s: ERROR -Q: Syntax is -Qa|v<min>/<max>\n", GMT_program);
						error++;
						break;
				}
				break;
				
			case 'S':		/* Assign start/stop position for sub-section (converted to meters) */
				limit_on_dist = TRUE;
				if (argv[i][2] == 'a') {		/* Start position */
					MGD77_Set_Unit (&argv[i][3], &dist_scale, 1);
					start_dist = atof (&argv[i][3]) * dist_scale;
				}
				else if (argv[i][2] == 'b') {	/* Stop position */
					MGD77_Set_Unit (&argv[i][3], &dist_scale, 1);
					stop_dist = atof (&argv[i][3]) * dist_scale;
				}
				else
					error = TRUE;
				break;
	
			case 'T':	/* Disable automatic corrections */
				switch (argv[i][2]) {
					case '\0':	/* Both sets */
						M.use_corrections[MGD77_M77_SET] = FALSE;
						M.use_corrections[MGD77_CDF_SET] = FALSE;
						break;
					case 'm':	/* MGD77 set */
						M.use_corrections[MGD77_M77_SET] = FALSE;
						break;
					case 'e':	/* extra CDF set */
						M.use_corrections[MGD77_CDF_SET] = FALSE;
						break;
					default:
						fprintf(stderr, "%s: ERROR -T: append m, e, or neither\n", GMT_program);
						error++;
						break;
				}
				break;
			case 'W':		/* Assign a weight to these data */
				aux_dvalue[MGD77_AUX_WT] = (!strcmp (&argv[i][2], "NaN")) ? GMT_d_NaN : atof (&argv[i][2]);
				break;
	
			case 'Z':		/* -Z- is negative down for depths */
				negative_msd = negative_depth = (argv[i][2] == '-') ? TRUE : FALSE;
				break;
	
			default:		/* Options not recognized */
				fprintf(stderr, "%s: ERROR -%s: option not recognized\n", GMT_program, &argv[i][1]);
				error = TRUE;
				break;
		}
	}
	
	/* Check that the options selected are mutually consistent */
	
	if (start_date && start_dist > 0.0) {
		fprintf(stderr, "%s: ERROR: Cannot specify both start time AND start distance\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (stop_date && stop_dist < DBL_MAX) {
		fprintf(stderr, "%s: ERROR: Cannot specify both stop time AND stop distance\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (east < west || south > north) {
		fprintf(stderr, "%s: ERROR: Region set incorrectly\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (aux_dvalue[MGD77_AUX_WT] <= 0.0) {
		fprintf(stderr, "%s: ERROR: -W weight must be nonzero\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (start_date && GMT_verify_expectations (GMT_IS_ABSTIME, GMT_scanf (start_date, GMT_IS_ABSTIME, &start_time), start_date)) {
		fprintf (stderr, "%s: ERROR -Da: Start time (%s) in wrong format\n", GMT_program, start_date);
		exit (EXIT_FAILURE);
	}
	if (stop_date && GMT_verify_expectations (GMT_IS_ABSTIME, GMT_scanf (stop_date, GMT_IS_ABSTIME, &stop_time), stop_date)) {
		fprintf (stderr, "%s: ERROR -Db : Stop time (%s) in wrong format\n", GMT_program, stop_date);
		exit (EXIT_FAILURE);
	}
	if (start_time > stop_time) {
		fprintf(stderr, "%s: ERROR -D: Start time exceeds stop time!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (start_dist > stop_dist) {
		fprintf (stderr, "%s: ERROR -S: Start distance exceeds stop distance!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (limit_on_azimuth && min_az >= max_az) {
		fprintf(stderr, "%s: ERROR -Qa: Minimum azimuth equals or exceeds maximum azimuth\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (limit_on_velocity && (min_vel >= max_vel || min_vel < 0.0)) {
		fprintf(stderr, "%s: ERROR -Qv: Minimum velocity equals or exceeds maximum velocity or is negative\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (GMT_give_synopsis_and_exit || argc == 1) {	/* Display usage */
		fprintf(stderr,"mgd77list %s - Extract data from MGD77 files\n\n", MGD77_VERSION);
		fprintf (stderr,"usage: mgd77list <cruise(s)> -F<dataflags>[,<tests>] [-A[+]c|d|f|m|t[code]] [-Cf|g|e] [-Da<startdate>] [-Db<stopdate>] [-E]\n");
		fprintf (stderr, "\t[-Ga<startrec>] [-Gb<stoprec>] [-H] [-I<code>] [-L[<corrtable.txt>]] [-N[s|p][e|k|n|m]]] [-Qa|v<min>/<max>] [%s]\n", GMT_Rgeo_OPT);
		fprintf (stderr, "\t[-Sa<startdist>[unit]] [-Sb<stopdist>[unit]] [-T[m|e]] [-V] [-W<Weight>] [-Z[+|-] [%s] [%s]\n\n", GMT_bo_OPT, GMT_mo_OPT);
         
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
              
		MGD77_Cruise_Explain ();
		fprintf (stderr, "\t-F <dataflags> is a comma-separated string made up of one or more of these abbreviations\n");
		fprintf (stderr, "\t   (for standard MGD77 files - use mgd77info to probe for other columns in MGD77+ files):\n");
		fprintf (stderr, "\t   >Track information:\n");
		fprintf (stderr, "\t     time:    Choose between Absolute time [default], Relative time, or fractional year:\n");
		fprintf (stderr, "\t       atime: Absolute time (formatted according to OUTPUT_DATE_FORMAT, OUTPUT_CLOCK_FORMAT)\n");
		fprintf (stderr, "\t       rtime: Relative time (formatted according to D_FORMAT and TIME_SYSTEM (or TIME_EPOCH, TIME_UNIT))\n");
		fprintf (stderr, "\t       ytime: Absolute time as decimal year (formatted according to D_FORMAT)\n");
		fprintf (stderr, "\t       year:  Record year\n");
		fprintf (stderr, "\t       month: Record month (1-12)\n");
		fprintf (stderr, "\t       day :  Record day of month (1-31)\n");
		fprintf (stderr, "\t       hour:  Record hour(0-23)\n");
		fprintf (stderr, "\t       min:   Record minute (0-59)\n");
		fprintf (stderr, "\t       sec:   Record second (0-60)\n");
		fprintf (stderr, "\t       dmin:  Decimal minute (0-59.xxxx)\n");
		fprintf (stderr, "\t       hhmm:  Clock hhmm.xxxx (0-2359.xxxx)\n");
		fprintf (stderr, "\t       date:  yyyymmdd string\n");
		fprintf (stderr, "\t       tz :   Time zone adjustment in hours (-13 to +12)\n");
		fprintf (stderr, "\t     lon:     Longitude (formatted according to OUTPUT_DEGREE_FORMAT)\n");
		fprintf (stderr, "\t     lat:     Latitude (formatted according to OUTPUT_DEGREE_FORMAT)\n");
		fprintf (stderr, "\t     id:      Survey leg ID [TEXTSTRING]\n");
		fprintf (stderr, "\t     ngdcid:  NGDC ID [TEXTSTRING]\n");
		fprintf (stderr, "\t     dist:    Along-track distances (see -C for method and -N for units)\n");
		fprintf (stderr, "\t     azim:    Track azimuth (Degrees east from north)\n");
		fprintf (stderr, "\t     vel:     Ship velocity (m/s)\n");
		fprintf (stderr, "\t   >Geophysical Observations:\n");
		fprintf (stderr, "\t     twt:     Two-way traveltime (s)\n");
		fprintf (stderr, "\t     depth:   Corrected bathymetry (m) [Also see -Z]\n");
		fprintf (stderr, "\t     mtf1:    Magnetic Total Field Sensor 1 (gamma, nTesla)\n");
		fprintf (stderr, "\t     mtf2:    Magnetic Total Field Sensor 2 (gamma, nTesla)\n");
		fprintf (stderr, "\t     mag:     Magnetic residual anomaly (gamma, nTesla)\n");
		fprintf (stderr, "\t     gobs:    Observed gravity (mGal)\n");
		fprintf (stderr, "\t     faa:     Free-air gravity anomaly (mGal)\n");
		fprintf (stderr, "\t   >Codes, Corrections, and Information:\n");
		fprintf (stderr, "\t     drt:     Data record type [5]\n");
		fprintf (stderr, "\t     ptc:     Position type code\n");
		fprintf (stderr, "\t     bcc:     Bathymetric correction code\n");
		fprintf (stderr, "\t     btc:     Bathymetric type code\n");
		fprintf (stderr, "\t     carter:  Carter correction from twt\n");
		fprintf (stderr, "\t     msens:   Magnetic sensor for residual field\n");
		fprintf (stderr, "\t     msd:     Magnetic sensor depth/altitude (m)\n");
		fprintf (stderr, "\t     diur:    Magnetic diurnal correction (gamma, nTesla)\n");
		fprintf (stderr, "\t     igrf:    International Geomagnetic Reference Field (gamma, nTesla)\n");
#ifdef USE_CM4
		fprintf (stderr, "\t     cm4:     Comprehensive Model CM4 Geomagnetic Reference Field (gamma, nTesla)\n");
#endif
		fprintf (stderr, "\t     eot:     Eotvos correction (mGal)\n");
		fprintf (stderr, "\t     ngrav:   IGF, or Theoretical (Normal) Gravity Field (mGal)\n");
		fprintf (stderr, "\t     sln:     Seismic line number string [TEXTSTRING]\n");
		fprintf (stderr, "\t     sspn:    Seismic shot point number string [TEXTSTRING]\n");
		fprintf (stderr, "\t     weight:  Give weight specified in -W\n");
		fprintf (stderr, "\t     nqc:     Navigation quality code\n");
		fprintf (stderr, "\t     bqc:     Bathymetric quality code, if available\n");
		fprintf (stderr, "\t     mqc:     Magnetics quality code, if available\n");
		fprintf (stderr, "\t     gqc:     Gravity quality code, if available\n");
		fprintf (stderr, "\t  The data are written in the order specified in <dataflags>\n");
		fprintf (stderr, "\t  Shortcut flags are:\n");
		fprintf (stderr, "\t     mgd77:   The full set of all 27 fields in the MGD77 specification\n");
		fprintf (stderr, "\t     mgd77t:  The full set of all 26 columns in the MGD77T specification\n");
		fprintf (stderr, "\t     geo:     time,lon,lat + the 7 geophysical observations\n");
		fprintf (stderr, "\t     all:     As mgd77 but with time items written as a date-time string\n");
		fprintf (stderr, "\t     allt:    As mgd77t but with time items written as a date-time string\n");
		fprintf (stderr, "\t    Append + to include the 4 derived quantities dist, azim, vel, and weight [see -W]\n");
		fprintf (stderr, "\t    [Default is all]\n");
		fprintf (stderr, "\t  Abbreviations in UPPER CASE will suppress records where any such column is NaN.\n");
		fprintf (stderr, "\t  (Note that -E is a shorthand to set all abbreviations to upper case).\n");
		fprintf (stderr, "\t  Optionally, append comma-separated logical tests that columns must pass to be output.\n");
		fprintf (stderr, "\t  Format is <flag><OP><value>, where flag is any of the dataflags above, and <OP> is\n");
		fprintf (stderr, "\t  one of the operators <, <=, =, >=, >, |, and !=.  <value> is the limit you are testing,\n");
		fprintf (stderr, "\t  including NaN (with = and != only).  If <flag> is UPPERCASE the test MUST be passed;\n");
		fprintf (stderr, "\t  else at least ONE of the tests must pass for output to take place.  When using operators\n");
		fprintf (stderr, "\t  involving characters <, >, and |, put entire argument to -F in single quotes.\n");
		fprintf (stderr, "\t  Finally, for MGD77+ files you may optionally append : followed by one or more comma-\n");
		fprintf (stderr, "\t  separated -+|-<col> terms.  This compares specific bitflags for each listed column\n");
		fprintf (stderr, "\t  + means bit must be 1, - means it must be 0.  All bit tests given must be passed.\n");
		fprintf (stderr, "\t  By default, MGD77+ files with error bit flags will use the flags to suppress bad data.\n");
		fprintf (stderr, "\t  Turn this behavior off by append : with no arguments.\n");
		fprintf (stderr, "\tOPTIONS:\n\n");
		fprintf (stderr, "\t-A Adjust some data values before output. Append c|d|f|m|t to select field:\n");
		fprintf (stderr, "\t   c<code>[,<v>] Adjust field carter. <v>, the sound velocity in water, is taken from\n");
		fprintf (stderr, "\t     the MGD77 header (or 1500 if invalid); optionally append your <v> (in m/s)\n");
		fprintf (stderr, "\t     Here, C(twt) is Carter correction, U(twt,v) is uncorrected depth (given <v>).\n");
		fprintf (stderr, "\t     TC(z) is twt from inverse Carter correction, TU(z,v) is twt from uncorrected depth.\n");
		fprintf (stderr, "\t       c1 return difference between U(twt,v) and depth [Default].\n");
		fprintf (stderr, "\t       c2 return difference between U(twt,v) and Carter(twt).\n");
		fprintf (stderr, "\t       c4 return difference between (uncorrected) depth and Carter (TU(depth,v)).\n");
		fprintf (stderr, "\t       c8 return difference between U(TC(depth),v) and depth.\n");
		fprintf (stderr, "\t   d<code>[,<v>] Adjust field depth. <v> is optional sound speed in water (m/s)\n");
		fprintf (stderr, "\t       d1 return depth as stored in file [Default].\n");
		fprintf (stderr, "\t       d2 return calculated uncorrected depth U(twt,v).\n");
		fprintf (stderr, "\t       d4 return calculated corrected depth Carter (twt,v).\n");
		fprintf (stderr, "\t   f<code>[,<field>] Adjust field faa. <field>, the IGF reference field, is taken\n");
		fprintf (stderr, "\t     from the MGD77 header (or 4 if invalid); optionally append your <field> from\n");
		fprintf (stderr, "\t     1 = Heiskanen 1924 formula:\n\t       ");
		MGD77_IGF_text (stderr, 1);
		fprintf (stderr, "\t     2 = International 1930 formula:\n\t       ");
		MGD77_IGF_text (stderr, 2);
		fprintf (stderr, "\t     3 = International 1967 formula:\n\t       ");
		MGD77_IGF_text (stderr, 3);
		fprintf (stderr, "\t     4 = International 1980 formula:\n\t       ");
		MGD77_IGF_text (stderr, 4);
		fprintf (stderr, "\t       f1 return faa as stored in file [Default].\n");
		fprintf (stderr, "\t       f2 return difference gobs - ngrav.\n");
		fprintf (stderr, "\t       f4 return difference gobs + eot - ngrav.\n");
		fprintf (stderr, "\t   m<code> Adjust field mag.\n");
		fprintf (stderr, "\t       m1 return mag as stored in file [Default].\n");
		fprintf (stderr, "\t       m2 return difference mtfx - igrf, where x = msens (or 1 if undefined).\n");
		fprintf (stderr, "\t       m4 return difference mtfx - igrf, where x != msens (or 2 if undefined).\n");
#ifdef USE_CM4
		fprintf (stderr, "\t       m8 return difference mtfx - cm4, where x = msens (or 1 if undefined).\n");
		fprintf (stderr, "\t       m16 return difference mtfx - cm4, where x != msens (or 2 if undefined).\n");
#endif
		fprintf (stderr, "\t   t will compute fake times for cruises with known duration but lacking record times\n");
		fprintf (stderr, "\t   The optional -A+ means selected anomalies will be recalculated even when the original\n");
		fprintf (stderr, "\t   anomaly is NaN [Default honors NaNs in existing anomalies]\n");
		fprintf (stderr, "\t-C Select procedure for along-track distance and azimuth calculations:\n");
		fprintf (stderr, "\t   f Flat Earth\n");
		fprintf (stderr, "\t   g Great circle [Default]\n");
		fprintf (stderr, "\t   e Ellipsoidal (geodesic) using current ellipsoid\n");
		fprintf (stderr, "\t-Da<date> lists from date (given as yyyy-mm-ddT[hh:mm:ss]) [Start of cruise]\n");
		fprintf (stderr, "\t  b<date> lists up to date (given as yyyy-mm-ddT[hh:mm:ss]) [End of cruise]\n");
		fprintf (stderr, "\t  If A|B is used instead or a|b then records with no time are excluded from output\n");
		fprintf (stderr, "\t-E Only records that exactly matches the requested geophysical information in -F will be used.\n");
		fprintf (stderr, "\t   [Default will output all record that matches at least one column]\n");
		fprintf (stderr, "\t-Ga<rec> lists from given record [Start of cruise]\n");
		fprintf (stderr, "\t  b<rec> lists up to given record [End of cruise]\n");
		fprintf (stderr, "\t-H Write one header record with column names\n");
		fprintf (stderr, "\t-I Ignore certain data file formats from consideration. Append combination of act to ignore\n");
		fprintf (stderr, "\t   (a) MGD77 ASCII, (c) MGD77+ netCDF, (m) MGD77T ASCII, or (t) plain table files. [Default ignores none]\n");
		fprintf (stderr, "\t-L Subtract systematic corrections from the data. If no correction file is given,\n");
		fprintf (stderr, "\t   the default file mgd77_corrections.txt in $MGD77_HOME is assumed.\n");
		fprintf (stderr, "\t-N Append (d)istances or (s)peed, and your choice for unit. Choose among:\n");
		fprintf (stderr, "\t   e Metric units I (meters, m/s)\n");
		fprintf (stderr, "\t   k Metric units II (km, km/hr)\n");
		fprintf (stderr, "\t   m British/US units (miles, miles/hr)\n");
		fprintf (stderr, "\t   n Nautical units (nautical miles, knots)\n");
		fprintf (stderr, "\t   [Default is -Ndk -Nse]\n");
		fprintf (stderr, "\t-Q Only return data whose azimuth (-Qa) or velocity (-Qv) fall inside specified range:\n");
		fprintf (stderr, "\t   -Qa<min_az>/<max_az>, where <min_az> < <max_az> [all azimuths, i.e., 0/360]\n");
		fprintf (stderr, "\t   -Qv<min_vel>[/<max_vel>], where <max_vel> is optional [all velocities, i.e., 0/infinity]\n");
		fprintf (stderr, "\t      Velocities are given in m/s unless changed by -Ns\n");
		fprintf (stderr, "\t-R Only return data inside the specified region [0/360/-90/90]\n");
		fprintf (stderr, "\t-Sa<dist> lists from dist  (in m; append k, m, or n) [Start of the cruise]\n");
		fprintf (stderr, "\t-Sb<dist> lists up to dist (in m; append k, m, or n) [End of the cruise]\n");
		fprintf (stderr, "\t-T turns OFF the otherwise automatic adjustment of values based on correction terms\n");
		fprintf (stderr, "\t   stored in the mgd77+ file (option has no effect on plain MGD77 ASCII files).\n");
		fprintf (stderr, "\t   Append m or e to indicate the MGD77 data set or the extended columns set. [Default is both]\n");
		fprintf (stderr, "\t-V Verbose, report progress\n");
		fprintf (stderr, "\t-W Sets weight for these data [1]\n");
		fprintf (stderr, "\t-Z Append - to report bathymetry & msd as negative depths [Default is positive -Z+]\n");
		GMT_explain_option ('o');
		fprintf (stderr, "\t-m Write multisegment header records for each cruise\n");
		exit (EXIT_FAILURE);
	}

	n_paths = MGD77_Path_Expand (&M, argv, argc, &list);	/* Get list of requested IDs */

	if (n_paths == 0) {
		fprintf(stderr, "%s: ERROR: No cruises given\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (M.adjust_time) start_time = MGD77_time2utime (&M, start_time);	/* Convert to Unix time if need be */
	if (M.adjust_time) stop_time = MGD77_time2utime (&M, stop_time);
	if (apply_corrections) {	/* Scan the ephemeral correction table for needed auxilliary columns */
		char path[BUFSIZ];
		if (!correction_table) {	/* Try default correction table */
			sprintf (path, "%s/mgd77_corrections.txt", M.MGD77_HOME);
			if (access (path, R_OK)) {
				fprintf (stderr, "%s: No default MGD77 Correction table (%s) found!\n", GMT_program, path);
				exit (EXIT_FAILURE);
			}
			correction_table = path;
		}
		n_items = MGD77_Scan_Corrtable (correction_table, list, n_paths, M.n_out_columns, (char **)M.desired_column, &item_names, 2);
	}
	
	limit_on_wesn = (!(west == 0.0 && east == 360.0 && south == -90.0 && north == 90.0));	/* User specified a sub region */
	select_option = MGD77_RESET_CONSTRAINT | MGD77_RESET_EXACT;	/* Make sure these start at zero */
	if (exact) select_option |= MGD77_SET_ALLEXACT;			/* Sets all columns listed as "must be present" */
	MGD77_Select_Columns (f_setting, &M, select_option);		/* This is the list of columns the user ultimately wants output */
	n_out_columns = M.n_out_columns;				/* This is the total number of columns in the final output */
	if (MGD77_Get_Column ("depth", &M) == MGD77_NOT_SET) negative_depth = FALSE;	/* Just so we don't accidently access dvalue[z_col] further down in the loop */
	if (MGD77_Get_Column ("msd", &M) == MGD77_NOT_SET) negative_msd = FALSE;	/* Just so we don't accidently access dvalue[m_col] further down in the loop */
	n_aux = (int)separate_aux_columns (&M, fx_setting, aux, auxlist);				/* Determine which auxillary columns are requested (if any) */
	if (apply_corrections) {
		n_aux = (int)augment_aux_columns (n_items, item_names, aux, auxlist, n_aux);		/* Determine which auxillary columns are needed by -L */
		for (i = 0; i < n_items; i++) GMT_free ((void *)item_names[i]);
		if (n_items) GMT_free ((void *)item_names);
	}
	aux_tvalue[MGD77_AUX_ID] = (char *)GMT_memory (VNULL, (size_t)GMT_TEXT_LEN, sizeof (char), GMT_program);	/* Just in case */
	aux_tvalue[MGD77_AUX_DA] = (char *)GMT_memory (VNULL, (size_t)GMT_TEXT_LEN, sizeof (char), GMT_program);	/* Just in case */
	use = (M.original) ? MGD77_ORIG : MGD77_REVISED;
	
	/* Most auxillary columns depend on values in the data columns.  If the user did not specify the required data columns
	 * then we must append them to make sure we have access to the values we need to calculate the auxillary values.
	 * Also, so limit tests on data records (e.g., distances, region, or time) also implies the need for certain data
	 * columns such as time, lon, and lat.
	 */
	 
	need_distances = (limit_on_dist || auxlist[MGD77_AUX_SP].requested || auxlist[MGD77_AUX_DS].requested || auxlist[MGD77_AUX_AZ].requested);	/* Distance is requested */
	need_lonlat = (auxlist[MGD77_AUX_MG].requested || auxlist[MGD77_AUX_GR].requested || auxlist[MGD77_AUX_CT].requested || adj_code[ADJ_MG] > 1 || adj_code[ADJ_DP] & 4 || adj_code[ADJ_CT] >= 2 || adj_code[ADJ_GR] > 1 || fake_times);	/* Need lon, lat to calculate reference fields or Carter correction */
	need_time = (auxlist[MGD77_AUX_YR].requested || auxlist[MGD77_AUX_MO].requested || auxlist[MGD77_AUX_DY].requested || auxlist[MGD77_AUX_HR].requested || auxlist[MGD77_AUX_MI].requested || auxlist[MGD77_AUX_SC].requested || auxlist[MGD77_AUX_DM].requested || auxlist[MGD77_AUX_HM].requested || auxlist[MGD77_AUX_DA].requested || auxlist[MGD77_AUX_MG].requested || (adj_code[ADJ_MG] > 1));
#ifdef USE_CM4
	if (auxlist[MGD77_AUX_CM].requested) need_lonlat = need_time = TRUE;
#endif

	n_sub = 0;	/* This value will hold the number of columns that we will NOT printout (they are only needed to calculate auxillary values) */
	if (need_distances || need_lonlat) {	/* Must make sure we get lon,lat if they are not already requested */
		 if (MGD77_Get_Column ("lat", &M) == MGD77_NOT_SET) strcat (fx_setting, ",lat"), n_sub++;	/* Append lat to requested list */
		 if (MGD77_Get_Column ("lon", &M) == MGD77_NOT_SET) strcat (fx_setting, ",lon"), n_sub++;	/* Append lon to requested list */
	}
	if ((limit_on_time || need_time || auxlist[MGD77_AUX_SP].requested) && MGD77_Get_Column ("time", &M) == MGD77_NOT_SET) strcat (fx_setting, ",time"), n_sub++;	/* Append time to requested list */
	need_twt = (auxlist[MGD77_AUX_CT].requested || (adj_code[ADJ_CT] > 0 && adj_code[ADJ_CT] < 3) || (adj_code[ADJ_DP] > 1));
	if (need_twt) {	/* Want to estimate Carter corrections */
		 if (MGD77_Get_Column ("twt", &M) == MGD77_NOT_SET) strcat (fx_setting, ",twt"), n_sub++;	/* Must append twt to requested list */
		MGD77_carter_init (&Carter);	/* Initialize Carter machinery */
	}
	need_depth = ((adj_code[ADJ_CT] & (1 | 3 | 8)) || (adj_code[ADJ_DP] & 1));
	if (need_depth) {	/* Need depth*/
		 if (MGD77_Get_Column ("depth", &M) == MGD77_NOT_SET) strcat (fx_setting, ",depth"), n_sub++;	/* Must append depth to requested list */
	}
	if (adj_code[ADJ_GR] > 1) {	/* Need gobs */
		 if (MGD77_Get_Column ("gobs", &M) == MGD77_NOT_SET) strcat (fx_setting, ",gobs"), n_sub++;	/* Must append gobs to requested list */
	}
	if (adj_code[ADJ_GR] == 3) {	/* Need eot */
		 if (MGD77_Get_Column ("eot", &M) == MGD77_NOT_SET) strcat (fx_setting, ",eot"), n_sub++;	/* Must append eot to requested list */
	}
	if (adj_code[ADJ_MG] > 1) {	/* Need mtf1,2, and msens */
		 if (MGD77_Get_Column ("mtf1", &M) == MGD77_NOT_SET) strcat (fx_setting, ",mtf1"), n_sub++;	/* Must append mtf1 to requested list */
		 if (MGD77_Get_Column ("mtf2", &M) == MGD77_NOT_SET) strcat (fx_setting, ",mtf2"), n_sub++;	/* Must append mtf2 to requested list */
		 if (MGD77_Get_Column ("msens", &M) == MGD77_NOT_SET) strcat (fx_setting, ",msens"), n_sub++;	/* Must append msens to requested list */
	}
	/* If logical tests are specified we must make sure the required columns are included as auxillary */
	for (i = 0; i < M.n_constraints; i++) {
		if (MGD77_Get_Column (M.Constraint[i].name, &M) != MGD77_NOT_SET) continue;	/* OK, already included */
		strcat (fx_setting, ",");
		strcat (fx_setting, M.Constraint[i].name);	/* Must add to our list */
		n_sub++;
	}
	need_sound = (((adj_code[ADJ_CT] & (1 | 2 | 8)) || adj_code[ADJ_DP] & 2) && sound_speed == 0.0);
	sound_speed *= 0.5;	/* Takes care of the 2 in 2-way travel time */
	MGD77_Select_Columns (fx_setting, &M, 0);	/* Only deal with col names - leave constraints/exacts unchanged from last call */
	n_cols_to_process = M.n_out_columns - n_sub;
	
	MGD77_Set_Unit (d_unit, &dist_scale, -1);	/* Gets scale which multiplies meters to chosen distance unit */
	MGD77_Set_Unit (s_unit, &vel_scale,  -1);	/* Sets output scale for distances using in velocities */
	switch (s_unit[0]) {
		case 'e':
			vel_scale /= dist_scale;			/* Must counteract any distance scaling to get meters. dt is in sec so we get  m/s */
			strcpy (auxlist[MGD77_AUX_SP].header, "v(m/s)");
			break;
		case 'k':
			vel_scale *= (3600.0 / dist_scale);		/* Must counteract any distance scaling to get km. dt is in sec so 3600 gives  km/hr */
			strcpy (auxlist[MGD77_AUX_SP].header, "v(km/hr)");
			break;
		case 'm':
			vel_scale *= (3600.0 / dist_scale);		/* Must counteract any distance scaling to get miles. dt is in sec so 3600 gives  miles/hr */
			strcpy (auxlist[MGD77_AUX_SP].header, "v(mi/hr)");
			break;
		case 'n':
			vel_scale *= (3600.0 / dist_scale);		/* Must counteract any distance scaling to get miles. dt is in sec so 3600 gives  miles/hr */
			strcpy (auxlist[MGD77_AUX_SP].header, "v(kts)");
			break;
	}
	switch (d_unit[0]) {
		case 'e':
			strcpy (auxlist[MGD77_AUX_SP].header, "d(m)");
			break;
		case 'k':
			strcpy (auxlist[MGD77_AUX_SP].header, "d(km)");
			break;
		case 'm':
			strcpy (auxlist[MGD77_AUX_SP].header, "d(miles)");
			break;
		case 'n':
			strcpy (auxlist[MGD77_AUX_SP].header, "d(nm)");
			break;
	}

	switch (dist_flag) {
		case 1:	/* Flat earth */
			GMT_distance_func = GMT_flatearth_dist_meter;
			GMT_azimuth_func  = GMT_az_backaz_flatearth;
			break;
		case 2:	/* Great circle */
			GMT_distance_func = GMT_great_circle_dist_meter;
			GMT_azimuth_func  = GMT_az_backaz_sphere;
			break;
		default:	/* Geodesic */
			GMT_distance_func = GMT_geodesic_dist_meter;
			GMT_azimuth_func  = GMT_az_backaz_geodesic;
			break;
	}

	start_dist *= dist_scale;	stop_dist *= dist_scale;	/* Convert the meters to the same units used for cumulative distances */

#ifdef USE_CM4
	if (auxlist[MGD77_AUX_CM].requested) MGD77_CM4_init (&M, &CM4);	/* Initialize CM4 structure */
#endif

	if (apply_corrections) {	/* Load an ephemeral correction table */
		char path[BUFSIZ];
		if (!correction_table) {	/* Try default correction table */
			sprintf (path, "%s/mgd77_corrections.txt", M.MGD77_HOME);
			if (access (path, R_OK)) {
				fprintf (stderr, "%s: No default MGD77 Correction table (%s) found!\n", GMT_program, path);
				exit (EXIT_FAILURE);
			}
			correction_table = path;
		}
		MGD77_Parse_Corrtable (correction_table, list, n_paths, M.n_out_columns, (char **)M.desired_column, 2, &CORR);
	}

	for (argno = 0; argno < n_paths; argno++) {		/* Process each ID */
	
		if (MGD77_Open_File (list[argno], &M, MGD77_READ_MODE)) continue;

		if (gmtdefs.verbose) fprintf (stderr, "%s: Now processing cruise %s\n", GMT_program, list[argno]);
		
		if (limit_institution && strncmp (M.NGDC_id, institution, (size_t)2)) {	/* Institution does not match specification */
			MGD77_Close_File (&M);
			continue;
		}
			
		D = MGD77_Create_Dataset ();

		error = MGD77_Read_Header_Record (list[argno], &M, &D->H);
		if (error) {
			if (error == MGD77_ERROR_NOSUCHCOLUMN)
				fprintf (stderr, "%s: One or more requested columns not present in cruise %s - skipping\n", GMT_program, list[argno]);
			else
				fprintf (stderr, "%s: Error reading header sequence for cruise %s - skipping\n", GMT_program, list[argno]);
			MGD77_Free (D);
			continue;
		}

		/* Having the header we can process -F and assign indices that refers to this particular data set */
		
		
		if (first_cruise) {
			for (i = 0, string_output = FALSE ; i < n_cols_to_process; i++) {	/* Prepare GMT output formatting machinery */
				if (D->H.info[M.order[i].set].col[M.order[i].item].text) string_output = TRUE;
			}
			if (auxlist[MGD77_AUX_ID].requested || auxlist[MGD77_AUX_DA].requested) string_output = TRUE;
			if (string_output && GMT_io.binary[1]) {
				fprintf(stderr, "%s: ERROR: Cannot specify binary output with text fields\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			first_cruise = FALSE;
			if (!string_output) out = (double *) GMT_memory (VNULL, (size_t)n_out_columns, sizeof(double), GMT_program);

		}
		
		if (MGD77_Read_Data (list[argno], &M, D)) {
			fprintf (stderr, "%s: Error reading data set for cruise %s\n", GMT_program, list[argno]);
			exit (EXIT_FAILURE);
		}
		MGD77_Close_File (&M);
		
		/* The 1*, 2*, 3* below is just there to ensure we dont end up with multiple cases all == MGD77_NOT_SET */
		time_column = ((i = MGD77_Get_Column ("time", &M)) != MGD77_NOT_SET && M.order[i].set == MGD77_M77_SET) ? M.order[i].item : 1 * MGD77_NOT_SET;
		lon_column  = ((i = MGD77_Get_Column ("lon",  &M)) != MGD77_NOT_SET && M.order[i].set == MGD77_M77_SET) ? M.order[i].item : 2 * MGD77_NOT_SET;
		lat_column  = ((i = MGD77_Get_Column ("lat",  &M)) != MGD77_NOT_SET && M.order[i].set == MGD77_M77_SET) ? M.order[i].item : 3 * MGD77_NOT_SET;
		
		if (time_column != MGD77_NOT_SET && GMT_io.binary[GMT_OUT] && gmtdefs.verbose && first_warning) {	/* Warn that binary time output is in Unix secs */
			if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: For binary output, time is stored as seconds since 1970 (Use TIME_SYSTEM=Unix to decode)\n", GMT_program);
			first_warning = FALSE;
		}
		for (i = kx = pos = 0; pos < n_out_columns; i++, pos++) {	/* Prepare GMT output formatting machinery */
			while (kx < n_aux && aux[kx].pos == i) {	/* Insert formatting for auxillary column (none are special) */
				GMT_io.out_col_type[pos] = GMT_IS_FLOAT;
				pos++, kx++;
			}
			if (i >= n_cols_to_process) continue;	/* Dont worry about helper columns that wont be printed */
			c  = M.order[i].set;
			id = M.order[i].item;
			if (c == 0 && id == time_column)	{	/* Special time formatting; time can only be in original set 0 */
				GMT_io.out_col_type[pos] = M.time_format;
				t_pos = pos;	/* Output order of time */
			}
			else if (id == lon_column)	/* Special lon formatting */
				GMT_io.out_col_type[pos] = (c == 0) ? GMT_IS_LON : GMT_IS_FLOAT;
			else if (id == lat_column)	/* Special lat formatting */
				GMT_io.out_col_type[pos] = (c == 0) ? GMT_IS_LAT : GMT_IS_FLOAT;
			else 		/* Everything else is float (not true for the 3 strings though) */
				GMT_io.out_col_type[pos] = GMT_IS_FLOAT;
		}
		
		if (!GMT_io.binary[GMT_OUT] && GMT_io.io_header[GMT_OUT]) {	/* Write out header record */
			GMT_fputs ("# ", GMT_stdout);
			for (i = kx = pos = 0; pos < n_out_columns; i++, pos++) {
				while (kx < n_aux && aux[kx].pos == i) {	/* Insert auxillary column */
					GMT_fputs (auxlist[aux[kx].type].header, GMT_stdout);
					if ((pos+1) < n_out_columns) GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
					pos++, kx++;
				}
				if (i >= n_cols_to_process) continue;
				c  = M.order[i].set;
				id = M.order[i].item;
				sprintf (buffer, "%7s", D->H.info[c].col[id].abbrev);
				GMT_fputs (buffer, GMT_stdout);
				if ((pos+1) < n_out_columns) GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
			}
			GMT_fputs ("\n", GMT_stdout);
		}

		if (GMT_io.multi_segments[GMT_OUT]) {	/* Write multisegment header between each cruise */
			sprintf (GMT_io.segment_header, "%c %s\n", GMT_io.EOF_flag[GMT_OUT], list[argno]);
			GMT_write_segmentheader (GMT_stdout, n_out_columns);
		}
		aux_dvalue[MGD77_AUX_DS] = cumulative_dist = ds = 0.0;
		if (auxlist[MGD77_AUX_ID].requested) strcpy (aux_tvalue[MGD77_AUX_ID], M.NGDC_id);
	
		t_col = MGD77_Get_Column ("time",   &M);
		x_col = MGD77_Get_Column ("lon",    &M);
		y_col = MGD77_Get_Column ("lat",    &M);
		z_col = MGD77_Get_Column ("depth",  &M);
		if (need_twt) twt_col = MGD77_Get_Column ("twt",  &M);
		if (adj_code[ADJ_GR]) f_col = MGD77_Get_Column ("faa",  &M);
		if (adj_code[ADJ_GR] > 1) g_col = MGD77_Get_Column ("gobs",  &M);
		if (adj_code[ADJ_GR] == 3) e_col = MGD77_Get_Column ("eot",  &M);
		if (adj_code[ADJ_MG]) m_col = MGD77_Get_Column ("mag",  &M);
		if (adj_code[ADJ_MG] > 1) {	/* Need more magnetics items */
			m1_col = MGD77_Get_Column ("mtf1",  &M);
			m2_col = MGD77_Get_Column ("mtf2",  &M);
			ms_col = MGD77_Get_Column ("msens",  &M);
		}
		if ((auxlist[MGD77_AUX_GR].requested || (adj_code[ADJ_GR] > 1 )) && GF_version == MGD77_NOT_SET) {
			GF_version = D->H.mgd77[use]->Gravity_Theoretical_Formula_Code - '0';
			if (GF_version < MGD77_IGF_HEISKANEN || GF_version > MGD77_IGF_1980) {
				fprintf (stderr, "%s: Invalid Gravity Theoretical Formula Code (%c) - default to %d\n", GMT_program, D->H.mgd77[use]->Gravity_Theoretical_Formula_Code, MGD77_IGF_1980);
				GF_version = MGD77_IGF_1980;
			}
		}
		for (i = 0; i < M.n_out_columns; i++) {
			dvalue[i] = (double *)D->values[i];
			tvalue[i] = (char *)D->values[i];
		}

		this_limit_on_time = limit_on_time;	/* Since we might change it below */
		if (time_column != MGD77_NOT_SET && D->H.no_time) {	/* Cannot know if ASCII MGD77 dont have time until after reading */
			GMT_LONG faked = FALSE;
			if (fake_times) {	/* Try to make fake times based on duration and distances */
				faked = MGD77_fake_times (&M, &(D->H), dvalue[x_col], dvalue[y_col], dvalue[t_col], D->H.n_records);
				if (faked && gmtdefs.verbose) fprintf (stderr, "%s: Warning: Time column for cruise %s created from distances and duration\n", GMT_program, list[argno]);
			}
			if (!faked && gmtdefs.verbose) {
				fprintf (stderr, "%s: Warning: Time column not present in cruise %s - set to NaN\n", GMT_program, list[argno]);
				if (this_limit_on_time) fprintf (stderr, "%s: Warning: -D limits cannot be used for cruise %s\n", GMT_program, list[argno]);
			}
			if (!faked && !skip_if_time_is_NaN) this_limit_on_time = FALSE;	/* To avoid pointless tests against NaN in loop */
		}
		if (need_sound) {	/* We opted to go with the value in the header [or 1500] */
			v = atof (D->H.mgd77[use]->Bathymetry_Assumed_Sound_Velocity) * 0.1;
			sound_speed = 0.5 * ((v < 1400.0 || v > 1600.0) ? 1500.0 : v);
		}
		
		if (sound_speed > 0.0) i_sound_speed = 1.0 / sound_speed;
		
		if (apply_corrections) MGD77_Init_Correction (CORR[argno], dvalue);	/* Initialize origins if needed */
		
		has_prev_twt = PDR_wrap = FALSE;
		twt_pdrwrap_corr = 0.0;
		
		/* Start processing records  */
		
		for (rec = 0, prevrec = -1; rec < D->H.n_records; rec++) {
		
			/* Compute accumulated distance along track (Great circles or Flat Earth) */
		
			if (need_distances) {
				lonlat_not_NaN = !( GMT_is_dnan (dvalue[x_col][rec]) || GMT_is_dnan (dvalue[y_col][rec]));
				if (rec == 0) {	/* Azimuth at 1st point set to azimuth of 2nd point since there is no previous point */
					if (auxlist[MGD77_AUX_AZ].requested) aux_dvalue[MGD77_AUX_AZ] = GMT_azimuth_func (dvalue[x_col][1], dvalue[y_col][1], dvalue[x_col][0], dvalue[y_col][0], TRUE);
				}
				else {		/* Need a previous point to calculate distance and heading */
					if (lonlat_not_NaN && prevrec >= 0) {	/* We have to records with OK lon,lat and can compute a distance from the previous OK point */
						ds = dist_scale * GMT_distance_func (dvalue[x_col][rec], dvalue[y_col][rec], dvalue[x_col][prevrec], dvalue[y_col][prevrec]);
						if (auxlist[MGD77_AUX_AZ].requested) aux_dvalue[MGD77_AUX_AZ] = GMT_azimuth_func (dvalue[x_col][rec], dvalue[y_col][rec], dvalue[x_col][prevrec], dvalue[y_col][prevrec], TRUE);
						cumulative_dist += ds;
						aux_dvalue[MGD77_AUX_DS] = cumulative_dist;
					}
					else {
						aux_dvalue[MGD77_AUX_DS] = GMT_d_NaN;
						if (auxlist[MGD77_AUX_AZ].requested) aux_dvalue[MGD77_AUX_AZ] = GMT_d_NaN;
					}
				}
				if (auxlist[MGD77_AUX_SP].requested) {
					if (rec == 0 || prevrec < 0) {	/* Initialize various counters */
						dt = dvalue[t_col][1] - dvalue[t_col][0];
						if (auxlist[MGD77_AUX_SP].requested) aux_dvalue[MGD77_AUX_SP] = (GMT_is_dnan (dt) || dt == 0.0) ? GMT_d_NaN : vel_scale * ds / dt;
					}
					else {		/* Need a previous point to calculate speed */
						dt = dvalue[t_col][rec] - dvalue[t_col][prevrec];
						if (auxlist[MGD77_AUX_SP].requested) aux_dvalue[MGD77_AUX_SP] = (GMT_is_dnan (dt) || dt == 0.0) ? GMT_d_NaN : vel_scale * ds / dt;
					}
				}
				if (lonlat_not_NaN) prevrec = rec;	/* This was a record with OK lon,lat; make it the previous point for distance calculations */
			}
			
			/* Check if rec no, time or distance falls outside specified ranges */
		
			if (limit_on_recs && (rec < start_rec || rec > stop_rec)) continue;
			if (limit_on_dist && (cumulative_dist < start_dist || cumulative_dist >= stop_dist)) continue;
			if (skip_if_time_is_NaN && GMT_is_dnan (dvalue[t_col][rec])) continue;
			if (this_limit_on_time && (dvalue[t_col][rec] < start_time || dvalue[t_col][rec] >= stop_time)) continue;
			if (limit_on_wesn) {	/* Check is lat/lon is outside specified area */
				if (dvalue[y_col][rec] < south || dvalue[y_col][rec] > north) continue;
				while (dvalue[x_col][rec] > east) dvalue[x_col][rec] -= 360.0;
				while (dvalue[x_col][rec] < west) dvalue[x_col][rec] += 360.0;
				if (dvalue[x_col][rec] > east) continue;
			}
			
			if (limit_on_velocity) {	/* Check if we are outside velocity range */
				if (aux_dvalue[MGD77_AUX_SP] < min_vel || aux_dvalue[MGD77_AUX_SP] > max_vel) continue;
			}
			
			if (limit_on_azimuth) {	/* Check if we are outside azimuth range */
				while (aux_dvalue[MGD77_AUX_AZ] > min_az) aux_dvalue[MGD77_AUX_AZ] -= 360.0;	/* Wind down to be sure az < min azimuth */
				while (aux_dvalue[MGD77_AUX_AZ] < min_az) aux_dvalue[MGD77_AUX_AZ] += 360.0;	/* Now add 360 until we pass min azimuth */	
				if (aux_dvalue[MGD77_AUX_AZ] > max_az) continue;				/* Outside azimuth range */
			}
			/* Check if it passes any given column data constraints */
			
			if (!MGD77_Pass_Record (&M, D, rec)) continue;	/* Failed the test */

			/* This record will now be printed out */
		
			if (need_time) {	/* Need auxillary time columns such as year, days etc, hence we get the calendar first, then use MGD77_cal_to_fyear */
				MGD77_gcal_from_dt (&M, dvalue[t_col][rec], &cal);	/* No adjust for TZ; this is GMT UTC time */
				aux_dvalue[MGD77_AUX_YR] = (double)cal.year;
				aux_dvalue[MGD77_AUX_MO] = (double)cal.month;
				aux_dvalue[MGD77_AUX_DY] = (double)cal.day_m;
				aux_dvalue[MGD77_AUX_HR] = (double)cal.hour;
				aux_dvalue[MGD77_AUX_MI] = (double)cal.min;
				aux_dvalue[MGD77_AUX_SC] = cal.sec;
				aux_dvalue[MGD77_AUX_DM] = cal.min + cal.sec / 60.0;
				aux_dvalue[MGD77_AUX_HM] = 100.0 * cal.hour + aux_dvalue[MGD77_AUX_DM];
				date = MGD77_cal_to_fyear (&cal);	/* Get date as decimal year */
				if (auxlist[MGD77_AUX_DA].requested) sprintf (aux_tvalue[MGD77_AUX_DA], "%4.4ld%2.2ld%2.2ld", cal.year, cal.month, cal.day_m);
				need_date = FALSE;
			}
			else
				need_date = TRUE;
			
			if (auxlist[MGD77_AUX_MG].requested) {	/* Evaluate IGRF */
				double date = 0.0;
				date = MGD77_cal_to_fyear (&cal);	/* Get date as decimal year */
				aux_dvalue[MGD77_AUX_MG] = (MGD77_igrf10syn (0, date, 1, 0.0, dvalue[x_col][rec], dvalue[y_col][rec], IGRF)) ? GMT_d_NaN : IGRF[MGD77_IGRF_F];
			}
#ifdef USE_CM4
			if (auxlist[MGD77_AUX_CM].requested) {	/* Evaluate CM4 */
				double date;
				date = MGD77_cal_to_fyear (&cal);	/* Get date as decimal year */
/* Change this --> */		aux_dvalue[MGD77_AUX_MG] = (MGD77_igrf10syn (0, date, 1, 0.0, dvalue[x_col][rec], dvalue[y_col][rec], IGRF)) ? GMT_d_NaN : IGRF[MGD77_IGRF_F];
			}
#endif

			if (auxlist[MGD77_AUX_GR].requested)	/* Evaluate Theoretical Gravity Model */
				aux_dvalue[MGD77_AUX_GR] = MGD77_Theoretical_Gravity (dvalue[x_col][rec], dvalue[y_col][rec], GF_version);

			if (auxlist[MGD77_AUX_CT].requested) {	/* Carter is one of the output columns */
				if (adj_code[ADJ_CT]) {	/* We have requested some adjustment to the carter value */
					aux_dvalue[MGD77_AUX_CT] = GMT_d_NaN;
					if (adj_code[ADJ_CT] & 1)	/* Try uncorr. depth - obs. depth */
						aux_dvalue[MGD77_AUX_CT] = dvalue[twt_col][rec] * sound_speed - dvalue[z_col][rec];	/* Factor of 2 dealt with earlier */
					if (adj_code[ADJ_CT] & 2 && GMT_is_dnan (aux_dvalue[MGD77_AUX_CT])) {	/* Try uncorr. depth - Carter depth */
						MGD77_carter_depth_from_xytwt (dvalue[x_col][rec], dvalue[y_col][rec], 1000.0 * dvalue[twt_col][rec], &Carter, &z);
						aux_dvalue[MGD77_AUX_CT] = dvalue[twt_col][rec] * i_sound_speed - z;
					}
					if (adj_code[ADJ_CT] & 4 && GMT_is_dnan (aux_dvalue[MGD77_AUX_CT])) {	/* Try uncorr. depth - inferred Carter depth */
						twt = dvalue[z_col][rec] * i_sound_speed;	/* Factor of 2 dealt with earlier */
						MGD77_carter_depth_from_xytwt (dvalue[x_col][rec], dvalue[y_col][rec], twt, &Carter, &z);
						aux_dvalue[MGD77_AUX_CT] = dvalue[z_col][rec] - z;
					}
					if (adj_code[ADJ_CT] & 8 && GMT_is_dnan (aux_dvalue[MGD77_AUX_CT])) {	/* Try inferred uncorr. depth - obs. depth */
						MGD77_carter_twt_from_xydepth (dvalue[x_col][rec], dvalue[y_col][rec], dvalue[z_col][rec], &Carter, &twt);
						z = twt * sound_speed;
						aux_dvalue[MGD77_AUX_CT] = z - dvalue[z_col][rec];
					}
				}
				else {
					twt = 1000.0 * dvalue[twt_col][rec];
					aux_dvalue[MGD77_AUX_CT] = MGD77_carter_correction (dvalue[x_col][rec], dvalue[y_col][rec], twt, &Carter);
				}
				if (negative_depth) aux_dvalue[MGD77_AUX_CT] = -aux_dvalue[MGD77_AUX_CT];
			}

			if (z_col != MGD77_NOT_SET && adj_code[ADJ_DP]) {	/* We have requested some adjustment to the depth value */
				z = GMT_d_NaN;
				if (adj_code[ADJ_DP] & 1)	/* Try obs. depth */
					z = dvalue[z_col][rec];
				if (adj_code[ADJ_DP] & 2 && GMT_is_dnan (z))	/* Try uncorr. depth */
					z = dvalue[twt_col][rec] * i_sound_speed;
				if (adj_code[ADJ_DP] & 4 && GMT_is_dnan (z)) {	/* Try Carter depth */
					twt = dvalue[twt_col][rec];
					if (!GMT_is_dnan (twt)) {	/* OK, valid twt */
						if (has_prev_twt) {	/* OK, may look at change in twt */
							d_twt = twt - prev_twt;
							if (fabs (d_twt) > TWT_PDR_WRAP_TRIGGER) {
								twt_pdrwrap_corr += copysign (TWT_PDR_WRAP, -d_twt);
								if (!PDR_wrap) fprintf (stderr, "%s: PDR travel time wrap detected for cruise %s\n", GMT_program, list[argno]);
								PDR_wrap = TRUE;
							}
						}
						has_prev_twt = TRUE;
						prev_twt = twt;
					}
					twt += twt_pdrwrap_corr;
					MGD77_carter_depth_from_xytwt (dvalue[x_col][rec], dvalue[y_col][rec], 1000.0 * twt, &Carter, &z);
				}
				if (adj_force || !GMT_is_dnan(dvalue[z_col][rec])) dvalue[z_col][rec] = z;
			}
			
			if (f_col != MGD77_NOT_SET && adj_code[ADJ_GR]) {	/* We have requested some adjustment to the faa value */
				g = GMT_d_NaN;
				if (adj_code[ADJ_GR] == 1)	/* Try faa */
					g = dvalue[f_col][rec];
				if (adj_code[ADJ_GR] == 2 && GMT_is_dnan (g))	/* Try gobs - ngrav */
					g = dvalue[g_col][rec] - MGD77_Theoretical_Gravity (dvalue[x_col][rec], dvalue[y_col][rec], GF_version);
				if (adj_code[ADJ_GR] == 3 && GMT_is_dnan (g))	/* Try gobs + eot - ngrav */
					g = dvalue[g_col][rec] + dvalue[e_col][rec] - MGD77_Theoretical_Gravity (dvalue[x_col][rec], dvalue[y_col][rec], GF_version);
				if (adj_force || !GMT_is_dnan(dvalue[f_col][rec])) dvalue[f_col][rec] = g;
			}
			
			if (m_col != MGD77_NOT_SET && adj_code[ADJ_MG]) {	/* We have requested some adjustment to the mag value */
				m = GMT_d_NaN;
				if (adj_code[ADJ_MG] & 1)	/* Try mag */
					m = dvalue[m_col][rec];
				if (adj_code[ADJ_MG] & 2 && GMT_is_dnan (m)) {	/* Try mtf 1st - igrf */
					if (need_date) {	/* Did not get computed already */
						date = MGD77_time_to_fyear (&M, dvalue[t_col][rec]);
						need_date = FALSE;
					}
					i = irint (dvalue[ms_col][rec]);
					k = (i == 2) ? m2_col : m1_col;
					m = MGD77_Recalc_Mag_Anomaly_IGRF (&M, date, dvalue[x_col][rec], dvalue[y_col][rec], dvalue[k][rec], FALSE);
				}
				if (adj_code[ADJ_MG] & 4 && GMT_is_dnan (m)) {	/* Try mtf 2nd - igrf */
					if (need_date) {	/* Did not get computed already */
						date = MGD77_time_to_fyear (&M, dvalue[t_col][rec]);
						need_date = FALSE;
					}
					i = irint (dvalue[ms_col][rec]);
					k = (i == 2) ? m1_col : m2_col;
					m = MGD77_Recalc_Mag_Anomaly_IGRF (&M, date, dvalue[x_col][rec], dvalue[y_col][rec], dvalue[k][rec], FALSE);
				}
#ifdef USE_CM4
				if (adj_code[ADJ_MG] & 8 && GMT_is_dnan (m)) {	/* Try mtf 1st - cm4 */
					if (need_date) {	/* Did not get computed already */
						date = MGD77_time_to_fyear (&M, dvalue[t_col][rec]);
						need_date = FALSE;
					}
					i = irint (dvalue[ms_col][rec]);
					k = (i == 2) ? m2_col : m1_col;
					m = MGD77_Recalc_Mag_Anomaly_CM4 (&M, date, dvalue[x_col][rec], dvalue[y_col][rec], dvalue[k][rec], FALSE, &CM4);
				}
				if (adj_code[ADJ_MG] & 16 && GMT_is_dnan (m)) {	/* Try mtf 2nd - cm4 */
					if (need_date) {	/* Did not get computed already */
						date = MGD77_time_to_fyear (&M, dvalue[t_col][rec]);
						need_date = FALSE;
					}
					i = irint (dvalue[ms_col][rec]);
					k = (i == 2) ? m1_col : m2_col;
					m = MGD77_Recalc_Mag_Anomaly_CM4 (&M, date, dvalue[x_col][rec], dvalue[y_col][rec], dvalue[k][rec], FALSE, &CM4);
				}
#endif
				if (adj_force || !GMT_is_dnan(dvalue[m_col][rec])) dvalue[m_col][rec] = m;
			}
			
			if (negative_depth) dvalue[z_col][rec] = -dvalue[z_col][rec];
			if (negative_msd) dvalue[m_col][rec] = -dvalue[m_col][rec];
			
			if (string_output) {	/* Must do it col by col and deal with the requested string(s) */
				for (i = kx = pos = 0; pos < n_out_columns; i++, pos++) {
					while (kx < n_aux && aux[kx].pos == i) {	/* Insert auxillary column */
						if (aux[kx].text)
							GMT_fputs (aux_tvalue[aux[kx].type], GMT_stdout);
						else
							GMT_ascii_output_one (GMT_stdout, aux_dvalue[aux[kx].type], pos);
						if ((pos+1) < n_out_columns) GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
						kx++, pos++;
					}
					if (i >= n_cols_to_process) continue;
					c  = M.order[i].set;
					id = M.order[i].item;
					if (D->H.info[c].col[id].text) {
						if (tvalue[i][rec*D->H.info[c].col[id].text] == 0)
							GMT_fputs ("NaN", GMT_stdout);
						else {
							for (k = 0; k < D->H.info[c].col[id].text && tvalue[i][rec*D->H.info[c].col[id].text+k]; k++) GMT_fputc ((int)tvalue[i][rec*D->H.info[c].col[id].text+k], GMT_stdout);
						}
					}
					else if (id == time_column) {	/* Time */
						if (GMT_io.out_col_type[pos] == GMT_IS_FLOAT) {	/* fractional year */
							if (need_date) {	/* Did not get computed already */
								date = MGD77_time_to_fyear (&M, dvalue[t_col][rec]);
								need_date = FALSE;
							}
						}
						else if (M.adjust_time)
							date = MGD77_utime2time (&M, dvalue[t_col][rec]);
						else
							date = dvalue[t_col][rec];
						GMT_ascii_output_one (GMT_stdout, date, pos);
					}
					else {
						correction = (apply_corrections) ? MGD77_Correction (CORR[argno][i].term, dvalue, aux_dvalue, rec) : 0.0;
						GMT_ascii_output_one (GMT_stdout, dvalue[i][rec] - correction, pos);
					}
					if ((pos+1) < n_out_columns) GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
				}
				GMT_fputs ("\n", GMT_stdout);
			}
			else {	/* Use GMT output machinery which can handle binary output, if requested */
				for (i = kx = pos = 0; pos < n_out_columns; i++, pos++) {
					while (kx < n_aux && aux[kx].pos == i) {	/* Insert auxillary column */
						out[pos] = aux_dvalue[aux[kx].type];
						pos++, kx++;
					}
					if (i >= n_cols_to_process) continue;
					if (pos == t_pos) {	/* This is the time column */
						if (GMT_io.out_col_type[pos] == GMT_IS_FLOAT) {	/* fractional year */
							if (need_date) {	/* Did not get computed already */
								date = MGD77_time_to_fyear (&M, dvalue[t_col][rec]);
								need_date = FALSE;
							}
							out[pos] = date;
						}
						else if (M.adjust_time)
							out[pos] = MGD77_utime2time (&M, dvalue[t_col][rec]);
						else
							out[pos] = dvalue[t_col][rec];
					}
					else {
						correction = (apply_corrections) ? MGD77_Correction (CORR[argno][i].term, dvalue, aux_dvalue, rec) : 0.0;
						out[pos] = dvalue[i][rec] - correction;
					}
				}
				GMT_output (GMT_stdout, n_out_columns, out);
			}
			n_out++;
		}
		MGD77_Free (D);
		n_cruises++;
	}
	
	if (!string_output) GMT_free ((void *)out);
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Returned %ld output records from %d cruises\n", GMT_program, n_out, n_cruises);
	
	MGD77_Path_Free (n_paths, list);
	if (apply_corrections) MGD77_Free_Correction (CORR, n_paths);
#ifdef USE_CM4
	if (auxlist[MGD77_AUX_CM].requested) MGD77_CM4_end (&CM4);	/* Free up CM4 structure */
#endif
	MGD77_end (&M);

	exit (EXIT_SUCCESS);
}

GMT_LONG separate_aux_columns (struct MGD77_CONTROL *F, char *fx_setting, struct MGD77_AUX_INFO *aux, struct MGD77_AUXLIST *auxlist)
{
	GMT_LONG i, j, k, this_aux, n_aux;
	
	fx_setting[0] = '\0';
	for (i = k = n_aux = 0; i < F->n_out_columns; i++) {
		for (j = 0, this_aux = MGD77_NOT_SET; j < N_MGD77_AUX && this_aux == MGD77_NOT_SET; j++) if (!strcmp (auxlist[j].name, F->desired_column[i])) this_aux = j;
		if (this_aux == MGD77_NOT_SET) {	/* Just pass other columns through */
			if (k) strcat (fx_setting, ",");
			strcat (fx_setting, F->desired_column[i]);
			k++;
		}
		else
		{	/* Found a request for an auxillary column  */
			aux[n_aux].type = auxlist[this_aux].type;
			aux[n_aux].text = auxlist[this_aux].text;
			aux[n_aux].pos = k;
			auxlist[this_aux].requested = TRUE;
			n_aux++;
		}
	}
	return (n_aux);
}

GMT_LONG augment_aux_columns (int n_items, char **item_name, struct MGD77_AUX_INFO *aux, struct MGD77_AUXLIST *auxlist, int n_aux)
{
	/* This adds additional aux colums that are required by the correction table and not already requested by other means (e.g. -F) */
	GMT_LONG i, j, k, this_aux, n;
	
	for (i = k = 0, n = n_aux; i < n_items; i++) {
		for (j = 0, this_aux = MGD77_NOT_SET; j < N_MGD77_AUX && this_aux == MGD77_NOT_SET; j++) if (!strcmp (auxlist[j].name, item_name[i])) this_aux = j;
		if (this_aux != MGD77_NOT_SET && !auxlist[this_aux].requested) {	/* Found a request for an auxillary column not yet requested  */
			aux[n].type = auxlist[this_aux].type;
			aux[n].text = auxlist[this_aux].text;
			aux[n].pos = k;
			auxlist[this_aux].requested = TRUE;
			n++;
		}
	}
	return (n);
}
