/*--------------------------------------------------------------------
 *	$Id: mgd77manage.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *    Copyright (c) 2005-2014 by P. Wessel
 * mgd77manage is used to (1) remove data columns from mgd77+ files
 * or (2) add a new data column to mgd77+ files.  Data can be added
 * from data tables, created from reference field formulas, or by
 * sampling grids along track.  Data from tables may be numbers or
 * text strings.
 *
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 *
 * Author:	Paul Wessel
 * Date:	18-OCT-2005
 * Version:	1.0
 *
 */
 
#include "mgd77.h"
#include "mgd77_e77.h"	/* E77 Header Errata Codes */

#define N_PAR		7
#define COL_SCALE	0
#define COL_OFFSET	1
#define IMG_SPACING	2
#define IMG_SCALE	3
#define IMG_MODE	4
#define IMG_LAT		5
#define COL_TYPE	6

#define ADD_IGRF	2
#define ADD_CARTER	3
#define ADD_GRAV	4
#define ADD_RMAG	5
#ifdef USE_CM4
#define ADD_CM4		6
#define ADD_RMAG4	7
#endif

#define N_E77_MODES	5
#define E77_HEADER_MODE	0
#define E77_TREND_MODE	1
#define E77_NAV_MODE	2
#define E77_VALUES_MODE	3
#define E77_SLOPES_MODE	4

#define set_bit(k) (1 << (k))

int main (int argc, char **argv)
{
	int cdf_var_id, n_dims = 0, dims[2];		/* netCDF variables should be declared as int */
	size_t start[2] = {0, 0}, count[2] = {0, 0};	/* NetCDF offset variables are size_t */
	
	GMT_LONG i, j, k = 0, ii, jj, argno, n_cruises = 0, n_paths = 0, mx = 0, my, column, result;
	GMT_LONG width, n_delete = 0, GF_version = MGD77_NOT_SET, n_expected_fields, n_fields = 0;
	GMT_LONG kind = GMT_IS_FLOAT, MTF_col = 1, set, dist_flag = 2, n_bad, n_alloc = GMT_CHUNK;
	GMT_LONG interpolant = BCR_BICUBIC, check, n_sampled = 0, n_changed = 0, n = 0, pos;
	
	time_t now;
	
	nc_type c_nc_type;
	
	GMT_LONG got_c = FALSE, got_a = FALSE, got_d = FALSE, got_g = FALSE, got_i = FALSE, got_n = FALSE, got_t = FALSE;
	GMT_LONG error = FALSE, replace = FALSE, delete = FALSE, interpolate = FALSE, transform, verified, force = FALSE;
	GMT_LONG strings = FALSE, got_grid, got_table, two_cols = FALSE, constant, ok_to_read = TRUE, got_e = FALSE;
	GMT_LONG ignore_verify = FALSE, e77_skip_mode[N_E77_MODES];
	
	char line[BUFSIZ], file[BUFSIZ], *d_list = NULL, p[BUFSIZ], history[BUFSIZ], c_size = 0, **list = NULL;
	char c_abbrev[GMT_TEXT_LEN], c_units[GMT_TEXT_LEN], c_name[MGD77_COL_NAME_LEN], c_comment[MGD77_COL_COMMENT_LEN];
	char no_char = '9', not_given[GMT_TEXT_LEN], word[BUFSIZ], **tmp_string = NULL, *text = NULL, d_unit[2];
	signed char LEN = 0, OLDLEN = 0;
	
	float *f = NULL;

	double threshold = 1.0, x, y, match_value, single_val, dist_scale = 1.0;
	double parameters[N_PAR], *xtmp, *coldnt = NULL, *colvalue = VNULL, *in = NULL, limits[2];

	struct MGD77_CONTROL In;
	struct MGD77_DATASET *D = NULL;
	struct GRD_HEADER grd;
	struct GMT_EDGEINFO edgeinfo;
	struct GMT_BCR bcr;
	struct MGD77_CARTER Carter;
	
	FILE *fp = NULL, *fp_err = NULL;

	GMT_LONG decode_A_options (GMT_LONG mode, char *line, char *file, double parameters[]);
	GMT_LONG decode_I_options (char *line, char *abbrev, char *name, char *units, char *size, char *comment, double parameters[]);
	GMT_LONG skip_if_missing (char *name, char *file, struct MGD77_CONTROL *F, struct MGD77_DATASET *D);
	GMT_LONG got_default_answer (char *line, char *answer);

	argc = (int)GMT_begin (argc, argv);		/* Initialize GMT Machinery */
	
	GMT_boundcond_init (&edgeinfo);
	parameters[COL_SCALE]   = 1.0;	/* Output column scaling */
	parameters[COL_OFFSET]  = 0.0;	/* Output column offset */
	parameters[IMG_SPACING] = 0.0;	/* IMG data minute spacing */
	parameters[IMG_SCALE]   = 1.0;	/* IMG data scaling */
	parameters[IMG_MODE]    = 0.0;	/* IMG data mode */
	parameters[IMG_LAT ]    = 0.0;	/* IMG lat not set */
	parameters[COL_TYPE]    = 0;	/* netCDF type */
	GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 2;
	d_unit[0] = d_unit[1] = '\0';
	memset ((void *)e77_skip_mode, 0, (size_t)(N_E77_MODES * sizeof (GMT_LONG)));
	/* Default e77_skip_mode will apply header and fix corrections if prefix is Y and set all data bits */
	
	MGD77_Init (&In);			/* Initialize MGD77 Machinery */

	for (i = 1; !error && i < argc; i++) {	/* Process input options */
		k = 0;
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {
			
				case 'V':
				case '\0':
					error += GMT_parse_common_options (argv[i], NULL, NULL, NULL, NULL);
					break;
					
				case 'A':	/* Adding a new column */
					k = 2;
					if (argv[i][k] == '+') {
						replace = TRUE;
						k++;
					}
					switch (argv[i][k]) {
						case 'a':	/* Plain column data file of exact same # of records */
							got_a = TRUE;
							error += decode_A_options (0, &argv[i][k+1], file, parameters);
							break;
						case 'c':	/* Add reference field or correction term */
							got_c = TRUE;
							error += decode_A_options (0, &argv[i][k+1], file, parameters);
							break;
						case 'D':	/* dist,val data file - interpolate to get values at all records */
							interpolate = TRUE;
						case 'd':	/* dist,val data file - only update records with matching distances */
							got_d = TRUE;
							error += decode_A_options (0, &argv[i][k+1], file, parameters);
							break;
						case 'E':	/* Plain E77 error flag file from mgd77sniffer */
							ignore_verify = TRUE;	/* Process raw e77 files that have not been verified */
						case 'e':	/* Plain E77 error flag file from mgd77sniffer */
							got_e = TRUE;
							while (argv[i][++k]) {
								switch (argv[i][k]) {
									case 'h':	/* Ignore all header recommendations regardless of Y/N prefix */
										e77_skip_mode[E77_HEADER_MODE] = TRUE;
										break;
									case 'f':	/* Ignore all systematic trend recommendations regardless of Y/N prefix */
										e77_skip_mode[E77_TREND_MODE] = TRUE;
										break;
									case 'n':	/* Ignore all NAV flags */
										e77_skip_mode[E77_NAV_MODE] = TRUE;
										break;
									case 'v':	/* Ignore all VALUE flags */
										e77_skip_mode[E77_VALUES_MODE] = TRUE;
										break;
									case 's':	/* Ignore all SLOPE flags */
										e77_skip_mode[E77_SLOPES_MODE] = TRUE;
										break;
									default:
										fprintf (stderr, "%s: ERROR: -Ae modifiers must be combination of hfnvs\n", GMT_program);
										break;
								}
							}
							break;
						case 'g':	/* Sample along track from this GMT grid file */
							got_g = TRUE;
							error += decode_A_options (0, &argv[i][k+1], file, parameters);
							break;
						case 'i':	/* Sample along track from this *.img grid file */
							got_i = TRUE;
							error += decode_A_options (1, &argv[i][k+1], file, parameters);
							break;
						case 'n':	/* recno,val data file - only update records with matching rec number */
							got_n = TRUE;
							error += decode_A_options (0, &argv[i][k+1], file, parameters);
							break;
						case 'T':	/* abstime,val data file - interpolate to get values at all records */
							interpolate = TRUE;
						case 't':	/* abstime,val data file - only update records with matching times */
							got_t = TRUE;
							kind = GMT_IS_ABSTIME;
							error += decode_A_options (0, &argv[i][k+1], file, parameters);
							break;
						default:
							fprintf (stderr, "%s: ERROR: -A modifier must be a|c|d|D|e|g|i|n|t|T\n", GMT_program);
							error++;
							break;
					}
					break;

				case 'C':	/* Distance calculation method */
					if (argv[i][2] == 'f') dist_flag = 1;
					if (argv[i][2] == 'g') dist_flag = 2;
					if (argv[i][2] == 'e') dist_flag = 3;
					if (dist_flag < 1 || dist_flag > 3) {
						fprintf(stderr, "%s: ERROR -C: Flag must be f, g, or e\n", GMT_program);
						error++;
					}
					break;
				case 'D':	/* Columns to delete */
					d_list = &argv[i][2];
					delete = TRUE;
					
				case 'E':	/* character to generate no-string value */
					no_char = argv[i][2];
					break;

				case 'F':	/* Force mode */
					force = TRUE;
					break;

				case 'I':	/* Column attribute information */
					error += decode_I_options (&argv[i][2], c_abbrev, c_name, c_units, &c_size, c_comment, parameters);
					break;
					
				case 'N':	/* Set distance units */
					d_unit[0] = argv[i][2];	d_unit[1] = 0;
					if (!strchr ("ekmn", (int)d_unit[0])) {
						fprintf(stderr, "%s: ERROR -N: Unit must be e, k, m, or n\n", GMT_program);
						error++;
					}
					break;
					
				case 'Q':	/* Interpolation parameters */
					interpolant = BCR_BILINEAR;
					for (j = 2; j < 5 && argv[i][j]; j++) {
						switch (argv[i][j]) {
							case 'n':
								interpolant = BCR_NEARNEIGHBOR; break;
							case 'l':
								interpolant = BCR_BILINEAR; break;
							case 'b':
								interpolant = BCR_BSPLINE; break;
							case 'c':
								interpolant = BCR_BICUBIC; break;
							case '/':
							default:
								threshold = atof (&argv[i][j]);
								if (j == 2 && threshold < GMT_SMALL) interpolant = BCR_NEARNEIGHBOR;
								j = 5; break;
						}
					}
					break;

				default:		/* Options not recognized */
					error = TRUE;
					break;
			}
		}
		else
			n_cruises++;
	}
	
	/* Check that the options selected are mutually consistent */
	
	if (GMT_give_synopsis_and_exit || argc == 1) {	/* Display usage */
		fprintf(stderr,"mgd77manage %s - Manage the content of MGD77+ files\n\n", MGD77_VERSION);
		fprintf(stderr,"usage: mgd77manage <cruise(s)> [-A[+]a|c|d|D|e|E|g|i|n|t|T<info>] [-Cf|g|e] [-Dname1,name2,...]\n");
		fprintf(stderr,"\t[-E<no_char>] [-F] [-I<abbrev>/<name>/<units>/<size>/<scale>/<offset>/\"comment\"]\n");
		fprintf(stderr,"\t[-Ne|k|m|n[+|-]] [-Q[b|c|l|n][[/]<threshold>]] [-V] [%s]\n\n", GMT_bi_OPT);
         
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
              
		fprintf (stderr, "\t<cruises> is one or more MGD77+ legnames, e.g., 01010083.\n");
		fprintf (stderr, "\n\tOPTIONS:\n\n");
		fprintf (stderr, "\t-A Append a new data column to the given files.  The code letters are:\n");
		fprintf (stderr, "\t   +: Optional.  Will overwrite an existing column with same name with new data.\n");
		fprintf (stderr, "\t   [Default will refuse if an existing column has the same abbreviation as the new data].\n");
		fprintf (stderr, "\t   a: Give filename with a new column to add.  We expect a single-column file\n");
		fprintf (stderr, "\t      with the same number of records as the MGD77 file.  Only one cruise can be set.\n");
		fprintf (stderr, "\t      If filename is - we read from stdin.\n");
		fprintf (stderr, "\t   c: Create a new column to be calculated from existing columns.  Add code:\n");
#ifdef USE_CM4
		fprintf (stderr, "\t        4 = CM4 field, m = IGRF total field, c = Carter correction, g = IGF (\"normal gravity\")\n");
		fprintf (stderr, "\t        R = recomputed magnetic anomaly rmag = mtfx - CM4 total field\n");
		fprintf (stderr, "\t        r = recomputed magnetic anomaly rmag = mtfx - IGRF total field\n");
#else
		fprintf (stderr, "\t        m = IGRF total field, c = Carter correction, g = IGF (\"normal gravity\")\n");
		fprintf (stderr, "\t        r = recomputed magnetic anomaly rmag = mtfx - IGRF total field\n");
#endif
		fprintf (stderr, "\t        Append x for which mtfx field to use (1 or 2) [1]\n");
		fprintf (stderr, "\t        For g, optionally append 1-4 to select the gravity formula to use:\n");
		fprintf (stderr, "\t        1 = Heiskanen 1924, 2 = International 1930, 3 = IGF1967, 4 = IGF1980.\n");
		fprintf (stderr, "\t        [Default uses formula specified in the MGD77 header, or 4 if not valid]\n");
		fprintf (stderr, "\t   d: Give filename with (dist [see -N], data) for a new column.  We expect a two-column file\n");
		fprintf (stderr, "\t      with distances (in km) in first column and data values in 2nd.  Only one cruise can be set.\n");
		fprintf (stderr, "\t      If filename is - we read from stdin.  Only records with mathcing distance will have data assigned.\n");
		fprintf (stderr, "\t   D: Same as d but we interpolate between the dist,data pairs to fill in all data records.\n");
		fprintf (stderr, "\t   e: Ingest MGD77 error/correction information (e77) produced by mgd77sniffer.  We will look\n");
		fprintf (stderr, "\t      for the <cruise>.e77 file in the current directory or in MGD77_HOME/E77 [%s/E77]\n", In.MGD77_HOME);
		fprintf (stderr, "\t      By default we will apply recommended header (h) and systematic fixes (f) and set all data bit flags.\n");
		fprintf (stderr, "\t      Append a combination of these flags to change the default accordingly:\n");
		fprintf (stderr, "\t        h = Ignore all header recommendations\n");
		fprintf (stderr, "\t        f = Ignore all systematic fixes recommendations\n");
		fprintf (stderr, "\t        n = Ignore data record bitflags pertaining to navigation (time, lon, lat).\n");
		fprintf (stderr, "\t        v = Ignore data record bitflags pertaining to data values.\n");
		fprintf (stderr, "\t        s = Ignore data record bitflags pertaining to data slopes (gradients).\n");
		fprintf (stderr, "\t      Use -DE to ignore the verification status of the e77 file [Default requires verification to be Y]\n");
		fprintf (stderr, "\t      NOTE:  Previous E77 information will be removed prior to processing this E77 information.\n");
		fprintf (stderr, "\t   g: Sample a GMT grid along track. (also see -Q).\n");
		fprintf (stderr, "\t      Append filename of the GMT grid.\n");
		fprintf (stderr, "\t   i: Sample a Sandwell/Smith *.img Mercator grid along track (also see -Q).\n");
		fprintf (stderr, "\t      Give filename and append comma-separated scale, mode, and optionally max latitude [%g].\n", GMT_IMG_MAXLAT_80);
		fprintf (stderr, "\t      The scale (0.1 or 1) is used to multiply after read; give mode as follows:\n");
		fprintf (stderr, "\t        0 = img file w/ no constraint code, interpolate to get data at track.\n");
                fprintf (stderr, "\t        1 = img file w/ constraints coded, interpolate to get data at track.\n");
                fprintf (stderr, "\t        2 = img file w/ constraints coded, gets data only at constrained points, NaN elsewhere.\n");
                fprintf (stderr, "\t        3 = img file w/ constraints coded, gets 1 at constraints, 0 elsewhere.\n");
                fprintf (stderr, "\t        For mode 2|3 you may want to consider the -Q<value> setting.\n");
		fprintf (stderr, "\t   n: Give filename with (rec_no, data) for a new column.  We expect a two-column file\n");
		fprintf (stderr, "\t      with record numbers (0 means 1st row) in first column and data values in 2nd.  Only one cruise can be set.\n");
		fprintf (stderr, "\t      If filename is - we read from stdin.  Only records with matching record numbers will have data assigned.\n");
		fprintf (stderr, "\t   t: Give filename with (abstime, data) for a new column.  We expect a two-column file\n");
		fprintf (stderr, "\t      with dateTclock strings in first column and data values in 2nd.  Only one cruise can be set.\n");
		fprintf (stderr, "\t      If filename is - we read from stdin.  Only records with matching times will have data assigned.\n");
		fprintf (stderr, "\t   T: Same as t but we interpolate between the time, data pairs to fill in all data records.\n");
		fprintf (stderr, "\t-C Append code for distance calculation procedure (when -Ad|D is set):\n");
		fprintf (stderr, "\t     f Flat Earth\n");
		fprintf (stderr, "\t     g Great circle [Default]\n");
		fprintf (stderr, "\t     e Ellipsoidal (geodesic) using current GMT ellipsoid\n");
		fprintf (stderr, "\t-D Delete the columns listed from all the cruise data files.\n");
		fprintf (stderr, "\t   The columns are removed before any data are added.  It is not a substitute for -A+.\n");
		fprintf (stderr, "\t   However, sometimes the shape of new data demands the old to be deleted first (you will be told)\n");
		fprintf (stderr, "\t-E Give character used to fill empty/missing string columns [9]\n");
		fprintf (stderr, "\t-F Force mode.  This allows you to even replace the standard MGD77 columns [only extended columns can be changed]\n");
		fprintf (stderr, "\t-I In addition to the file information above, you must also specify column information:\n");
		fprintf (stderr, "\t      abbrev:  Short, abbreviated word (lower case only), like satfaa (%d char max)\n", MGD77_COL_ABBREV_LEN);
		fprintf (stderr, "\t      name:    Descriptive name, like \"Geosat/ERS-1 Free-air gravity\" (%d char max)\n", MGD77_COL_NAME_LEN);
		fprintf (stderr, "\t      units:   Units for the column (e.g., mGal, gamma, km) (%d char max)\n", MGD77_COL_NAME_LEN);
		fprintf (stderr, "\t      size:    Either t(ext), b(yte), s(hort), f(loat), i(nt), or d(ouble)\n");
		fprintf (stderr, "\t      scale:   Multiply data by this scale before writing to mgd77+ file\n");
		fprintf (stderr, "\t      offset:  Add after scaling before writing to mgd77+ file\n");
		fprintf (stderr, "\t      comment: Any text (in double quotes) for information about column (%d char max)\n", MGD77_COL_COMMENT_LEN);
		fprintf (stderr, "\t      -I is ignored by -Ae\n");
		fprintf (stderr, "\t   Note for text: Interpolation is not allowed, and \"not-a-string\" is created from -E.\n");
		fprintf (stderr, "\t-N Append your choice for distance unit (if -Ad|D are set). Choose among:\n");
		fprintf (stderr, "\t   (e) meter, (k) km, (m) miles, or (n) nautical miles [Default is -Nk]\n");
		fprintf (stderr, "\t    See -C for selecting distance calculation procedure.\n");
		fprintf (stderr, "\t-Q Quick mode, use bilinear rather than bicubic [Default] interpolation.\n");
		fprintf (stderr, "\t   Alternatively, select interpolation mode by adding b = B-spline, c = bicubic,\n");
		fprintf (stderr, "\t   l = bilinear, or n = nearest-neighbor.\n");
		fprintf (stderr, "\t   Optionally, append <threshold> in the range [0,1]. [Default = 1 requires all\n");
		fprintf (stderr, "\t   4 or 16 nodes to be non-NaN.], <threshold> = 0.5 will interpolate about 1/2 way\n");
		fprintf (stderr, "\t   from a non-NaN to a NaN node, while 0.1 will go about 90%% of the way, etc.\n");
		fprintf (stderr, "\t   -Q0 will return the value of the nearest node instead of interpolating (Same as -Qn).\n");
		GMT_explain_option ('V');
		GMT_explain_option ('i');
		exit (EXIT_FAILURE);
	}

	got_table = (got_a || got_d || got_n || got_t);	/* Got a table to read */
	got_grid = (got_g || got_i);					/* Got a grid to read */
	c_nc_type = (nc_type) irint (parameters[COL_TYPE]);		/* NC data type */
	strings = (c_nc_type == NC_CHAR);				/* TRUE if our new column contains strings */
	
	n_paths = MGD77_Path_Expand (&In, argv, argc, &list);	/* Get list of requested IDs */

	if (n_paths == 0) {
		fprintf(stderr, "%s: ERROR: No cruises given\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if ((got_table + got_grid) > 1) {
		fprintf (stderr, "%s: ERROR: You must select one, and only one, of the -A options.\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if ((interpolate + strings) > 1) {
		fprintf (stderr, "%s: ERROR: Cannot interpolate column if data are strings\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (got_table && got_c) {
		fprintf (stderr, "%s: ERROR: Only one -A option can be specified\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (got_table && n_paths != 1) {
		fprintf (stderr, "%s: ERROR: With -Aa|d|D|n|t|T you can only select one cruise at the time.\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (!got_grid && interpolant != BCR_BICUBIC) {
		fprintf (stderr, "%s: ERROR -Q:  Requires -Ag|i.\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (threshold < 0.0 || threshold > 1.0) {
		fprintf (stderr, "%s: ERROR -Q:  threshold must be in [0,1] range\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (!(delete || got_e)) {
		if (strlen (c_abbrev) > MGD77_COL_ABBREV_LEN) {
			fprintf(stderr, "%s: ERROR: Column abbreviation too long - %d characters is maximum!\n", GMT_program, MGD77_COL_ABBREV_LEN);
			exit (EXIT_FAILURE);
		}
		if (strlen (c_name) > MGD77_COL_NAME_LEN) {
			fprintf(stderr, "%s: ERROR: Column name too long - %d characters is maximum!\n", GMT_program, MGD77_COL_NAME_LEN);
			exit (EXIT_FAILURE);
		}
		if (strlen (c_comment) > MGD77_COL_COMMENT_LEN) {
			fprintf(stderr, "%s: ERROR: Column comment too long - %d characters is maximum!\n", GMT_program, MGD77_COL_COMMENT_LEN);
			exit (EXIT_FAILURE);
		}
	}
	MGD77_Set_Unit (d_unit, &dist_scale, -1);	/* Gets scale which multiplies meters to chosen distance unit */

	memset ((void *)not_given, (int)no_char, (size_t)GMT_TEXT_LEN);	/* Text representing "no text value" */
	not_given[GMT_TEXT_LEN-1] = '\0';
	fp_err = (In.verbose_dest == 1) ? GMT_stdout : stderr;
	
	if (got_c) {	/* Calculate values to be stored */
		GMT_LONG version = MGD77_NOT_SET, mfield = 1;
		/* "file" is here either m, c, or g[1-4] */
		if (file[0] == 'm' && file[1] == '\0') {
			got_c = ADD_IGRF;
		}
#ifdef USE_CM4
		else if (file[0] == '4' && file[1] == '\0') {
			got_c = ADD_CM4;
		}
#endif
		else if (file[0] == 'c' && file[1] == '\0') {
			got_c = ADD_CARTER;
			MGD77_carter_init (&Carter);	/* Initialize Carter machinery */
		}
		else if (file[0] == 'g' && (file[1] == '\0' || ((version = (file[1] - '0')) >= MGD77_IGF_HEISKANEN && version <= MGD77_IGF_1980)) ) {
			got_c = ADD_GRAV;
			GF_version = version;
		}
		else if (file[0] == 'r' && (file[1] == '\0' || ((mfield = (file[1] - '0')) >= 1 && mfield <= 2)) ) {
			got_c = ADD_RMAG;
			MTF_col = mfield;
		}
#ifdef USE_CM4
		else if (file[0] == 'R' && (file[1] == '\0' || ((mfield = (file[1] - '0')) >= 1 && mfield <= 2)) ) {
			got_c = ADD_RMAG4;
			MTF_col = mfield;
		}
#endif
		else {
#ifdef USE_CM4
			fprintf(stderr, "%s: ERROR: -Ac expects 4, m, c, or g[1-4]\n", GMT_program);
#else
			fprintf(stderr, "%s: ERROR: -Ac expects m, c, or g[1-4]\n", GMT_program);
#endif
			exit (EXIT_FAILURE);
		}
	}
	else if (got_e) {	/* Do E77 work by ignoring previous E77 settings */
		In.use_flags[MGD77_M77_SET] = In.use_flags[MGD77_CDF_SET] = FALSE;		/* Turn use of flag bits OFF */
		In.use_corrections[MGD77_M77_SET] = In.use_corrections[MGD77_CDF_SET] = FALSE;	/* Turn use of systematic corrections OFF */
	}
	else if (got_g) {	/* Read regular GMT grid */

		GMT_err_fail (GMT_read_grd_info (file, &grd), file);
		if (GMT_360_RANGE (grd.x_max, grd.x_min)) GMT_boundcond_parse (&edgeinfo, "g");
	
		GMT_boundcond_param_prep (&grd, &edgeinfo);
	
		/* Initialize bcr structure with 2 rows/cols boundaries:  */

		GMT_bcr_init (&grd, GMT_pad, interpolant, threshold, &bcr);
		
		mx = grd.nx + GMT_pad[0] + GMT_pad[2];	my = grd.ny + GMT_pad[1] + GMT_pad[3];

		f = (float *) GMT_memory (VNULL, (size_t)(mx * my), sizeof (float), GMT_program);

		GMT_err_fail (GMT_read_grd (file, &grd, f, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), file);
		interpolate = (threshold > 0.0);
	}
	else if (got_i) {	/* Read Sandwell/Smith IMG file */
		
		GMT_read_img (file, &grd, &f, 0.0, 0.0, 0.0, 0.0, parameters[IMG_SCALE], (GMT_LONG)irint(parameters[IMG_MODE]), parameters[IMG_LAT], TRUE);
		if (GMT_360_RANGE (grd.x_max, grd.x_min)) GMT_boundcond_parse (&edgeinfo, "g");
		GMT_boundcond_param_prep (&grd, &edgeinfo);
		mx = grd.nx + GMT_pad[0] + GMT_pad[2];	my = grd.ny + GMT_pad[1] + GMT_pad[3];
		interpolate = (threshold > 0.0);
	}
	else if (got_table) {	/* Got a one- or two-column table to read */
		GMT_LONG n_ave = 0;
		double last_dnt = -DBL_MAX, sum_z = 0.0;
		char *not_used = NULL;
		
		if (file[0] == '-') {   /* Just read from standard input */
			fp = GMT_stdin;
#ifdef SET_IO_MODE
			GMT_setmode (GMT_IN);
#endif
		}
		else {
			if ((fp = GMT_fopen (file, GMT_io.r_mode)) == NULL) {
				fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, file);
				exit (EXIT_FAILURE);
			}
		}

		/* Skip any header records */
		if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) not_used = GMT_fgets (line, BUFSIZ, fp);

		two_cols = (got_d || got_n || got_t);
		n = (two_cols) ? -1 : 0;
		n_alloc = GMT_CHUNK;
		n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : GMT_MAX_COLUMNS;
		colvalue = (double *) GMT_memory (VNULL, n_alloc, sizeof (double), GMT_program);
		if (two_cols) {	/* Got an abscissae column as well (dnt: d = dist, n = rec number, t = time) */
			coldnt = (double *) GMT_memory (VNULL, n_alloc, sizeof (double), GMT_program);
			GMT_io.in_col_type[0] = kind;
		}
		if (strings && !two_cols) {	/* Must read strings directly from file since GMT_input would barf */
			ok_to_read = FALSE;
			tmp_string = (char **) GMT_memory (VNULL, n_alloc, sizeof (char *), GMT_program);
			while (GMT_fgets (word, BUFSIZ, fp)) {
				if (word[0] == '#') continue;
				width = strlen (word);
				tmp_string[n] = (char *) GMT_memory (VNULL, (size_t)(width+1), sizeof (char), GMT_program);
				strcpy (tmp_string[n], word);
				if (width > LEN) LEN = (signed char)width;
				n++;
			}
		}
		else if (strings) {		/* Pretend to read one column and get the text string form the text record */
			tmp_string = (char **) GMT_memory (VNULL, n_alloc, sizeof (char *), GMT_program);
			n_expected_fields = 1;
		}
		
		if (ok_to_read) n_fields = GMT_input (fp, &n_expected_fields, &in);

		while (ok_to_read && ! (GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

			while ((GMT_io.status & GMT_IO_SEGMENT_HEADER) && !(GMT_io.status & GMT_IO_EOF)) {
				n_fields = GMT_input (fp, &n_expected_fields, &in);
			}
			if ((GMT_io.status & GMT_IO_EOF)) continue;	/* At EOF */

			if (GMT_io.status & GMT_IO_MISMATCH) {
				fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld\n", GMT_program, n_fields, n_expected_fields, n);
				exit (EXIT_FAILURE);
			}

			if (strings) {	/* number in col1, string in col2 */
				coldnt[n]   = in[0];
				sscanf (GMT_io.current_record, "%*s %s", word);
				tmp_string[n] = (char *) GMT_memory (VNULL, (size_t)(strlen(word)+1), sizeof (char), GMT_program);
				strcpy (tmp_string[n], word);
			}
			else if (two_cols) {
				if (in[0] > last_dnt) {	/* Start of new averaging scheme (if needed) */
					if (n_ave) {
						colvalue[n] = sum_z / n_ave;
						coldnt[n]   = last_dnt;
					}
					n_ave = 0;
					sum_z = 0.0;
					n++;
					last_dnt = in[0];
				}
				sum_z += in[1];	/* Add up possibly multiple values for same dnt */
				n_ave++;
			}
			else
				colvalue[n] = in[0];
			if (!two_cols) n++;
			if (n == n_alloc) {
				n_alloc <<= 1;
				if (strings)
					tmp_string = (char **) GMT_memory ((void *)tmp_string, n_alloc, sizeof (char *), GMT_program);
				else
					colvalue = (double *) GMT_memory ((void *)colvalue, n_alloc, sizeof (double), GMT_program);
				if (two_cols) coldnt = (double *) GMT_memory ((void *)coldnt, n_alloc, sizeof (double), GMT_program);
			}

			n_fields = GMT_input (fp, &n_expected_fields, &in);
		}
		if (fp != GMT_stdin) GMT_fclose (fp);
		if (two_cols && n_ave) { colvalue[n] = sum_z / n_ave; coldnt[n++] = last_dnt;}
		if (!strings) colvalue = (double *) GMT_memory ((void *)colvalue, (size_t)n, sizeof (double), GMT_program);
		if (two_cols) coldnt = (double *) GMT_memory ((void *)coldnt, (size_t)n, sizeof (double), GMT_program);
	}
	
	MGD77_Ignore_Format (MGD77_FORMAT_ANY);	/* Reset to all formats OK, then ... */
	MGD77_Ignore_Format (MGD77_FORMAT_M77);	/* disallow ASCII MGD77 files */
	MGD77_Ignore_Format (MGD77_FORMAT_TBL);	/* and ASCII tables */
	
	In.format = MGD77_FORMAT_CDF;	/* Only file type allowed as input */
	
	for (argno = 0; argno < n_paths; argno++) {		/* Process each ID */
	
		if (MGD77_Open_File (list[argno], &In, MGD77_READ_MODE)) continue;
				
		if (gmtdefs.verbose) fprintf (stderr, "%s: Now processing cruise %s\n", GMT_program, list[argno]);
		
		D = MGD77_Create_Dataset ();
		In.n_out_columns = 0;

		if (MGD77_Read_File (list[argno], &In, D)) {
			fprintf (stderr, "%s: Error reading data set for cruise %s\n", GMT_program, list[argno]);
			exit (EXIT_FAILURE);
		}

		/* Start reading data from file */
	
		column = MGD77_Get_Column (c_abbrev, &In);
		set    = MGD77_Get_Set (c_abbrev);
		
		if (!got_e && column != MGD77_NOT_SET) {	/* A column with same abbreviation is already present in the file */
			if (set == MGD77_M77_SET && !force) {
				fprintf (stderr, "%s: column %s is part of the standard MGD77 set and cannot be removed unless you use -F!\n", GMT_program, c_abbrev);
				exit (EXIT_FAILURE);
			}
			if (!replace) {
				fprintf (stderr, "%s: A columned named %s is already present in %s.  use -A+ to overwrite [default is to skip]\n", GMT_program, c_abbrev, list[argno]);
				MGD77_Free (D);	/* Free memory already allocated by MGD77_Read_File for this aborted effort */
				continue;
			}
			n_dims = (D->H.info[In.order[column].set].col[In.order[column].item].constant) ? 0 : 1;
			if (D->H.info[In.order[column].set].col[In.order[column].item].text) n_dims++;
		}

		if (delete) {	/* Must create a new file with everything except the fields to be deleted */
			GMT_LONG id, c, reset_column = FALSE;
			char oldfile[BUFSIZ];
			
			if (column != MGD77_NOT_SET) {	/* Get info about this existing column to see if it is compatible with new data */
				n_dims = (D->H.info[In.order[column].set].col[In.order[column].item].constant) ? 0 : 1;
				if (D->H.info[In.order[column].set].col[In.order[column].item].text) n_dims++;
			}
			
			pos = n_delete = 0;
			(void) time (&now);
			sprintf (history, "%s [%s] removed columns", ctime(&now), In.user);
			for (i = 0; history[i]; i++) if (history[i] == '\n') history[i] = ' ';	/* Remove the \n returned by ctime() */
			while ((GMT_strtok (d_list, ",", &pos, p))) {	/* For each named column */
				k = MGD77_Get_Column (p, &In);
				if (k == MGD77_NOT_SET) {
					fprintf (stderr, "%s: No column named %s in %s - cannot delete it. \n", GMT_program, p, list[argno]);
					continue;
				}
				c = In.order[k].set;
				id = In.order[k].item;
				D->H.info[c].col[id].abbrev[0] = D->H.info[c].col[id].name[0] = D->H.info[c].col[id].units[0] = D->H.info[c].col[id].comment[0] = '\0';
				D->H.info[c].col[id].pos = D->H.info[c].col[id].var_id = MGD77_NOT_SET;
				D->H.info[c].bit_pattern = 0;
				D->H.info[c].col[id].present = FALSE;
				D->H.info[c].n_col--;
				D->H.n_fields--;
				D->n_fields--;
				In.n_out_columns--;
				for (i = k; i < In.n_out_columns; i++) {	/* Move remaining columns over */
					D->values[i] = D->values[i+1];
					strcpy (In.desired_column[i], In.desired_column[i+1]);
					In.order[i].set = In.order[i+1].set;
					In.order[i].item = In.order[i+1].item;
				}
				strcat (history, " ");
				strcat (history, p);
				n_delete++;
				if (k == column && c == set) {	/* Just removed the old column by the same name, must unset column */
					reset_column = TRUE;
				}
			}
			
			/* Rename the old file for now */
			
			sprintf (oldfile, "%s.old", In.path);
			if (rename (In.path, oldfile)) {
				fprintf (stderr, "%s: Unable to rename %s to %s\n", GMT_program, In.path, oldfile);
				exit (EXIT_FAILURE);
			}
			
			/* Update header history */

			k = strlen (history);
			for (i = 0; i < k; i++) if (history[i] == '\n') history[i] = ' ';	/* Remove the \n returned by ctime() */
			history[k++] = '\n';	history[k] = '\0';				/* Add LF at end of line */
			k += strlen (D->H.history);
			D->H.history = (char *)GMT_memory ((void *)D->H.history, (size_t)k, sizeof (char), GMT_program);
			strcat (D->H.history, history);		/* MGD77_Write_FILE_cdf will use this to create the history attribute, thus preserving earlier history */

			if (MGD77_Write_File (In.path, &In, D)) {	/* Create the new, slimmer file */
				fprintf (stderr, "%s: Error writing slimmer version of %s\n", GMT_program, list[argno]);
				exit (EXIT_FAILURE);
			}

			/* Now we can safely remove the old file */
			
			if (remove (oldfile)) {
				fprintf (stderr, "%s: Error removing the old version of %s\n", GMT_program, list[argno]);
				exit (EXIT_FAILURE);
			}
			
			MGD77_Free (D);
			if (column == MGD77_NOT_SET) continue;	/* Nothing more to do for this file */
			
			/* Now reread header etc since things have changed in the file */
			
			In.n_out_columns = 0;
			D = MGD77_Create_Dataset ();
			if (MGD77_Read_File (list[argno], &In, D)) {
				fprintf (stderr, "%s: Error reading data set for cruise %s\n", GMT_program, list[argno]);
				exit (EXIT_FAILURE);
			}
			if (reset_column)
				column = MGD77_NOT_SET;
			else
				n_changed++;
		}

		if (got_c == ADD_IGRF) {	/* Append IGRF column */
			GMT_LONG ix, iy, it;
			double date, *xvar, *yvar, *tvar, IGRF[7];
			
			if ((ix = skip_if_missing ("lon",  list[argno], &In, D)) == MGD77_NOT_SET) continue;
			if ((iy = skip_if_missing ("lat",  list[argno], &In, D)) == MGD77_NOT_SET) continue;
			if ((it = skip_if_missing ("time", list[argno], &In, D)) == MGD77_NOT_SET) continue;

			xvar = (double *)D->values[ix];
			yvar = (double *)D->values[iy];
			tvar = (double *)D->values[it];
			colvalue = (double *)GMT_memory (VNULL, (size_t)D->H.n_records, sizeof (double), GMT_program);
			
			for (i = n_sampled = 0; i < D->H.n_records; i++) {
				date = MGD77_time_to_fyear (&In, tvar[i]);	/* Get date as decimal year */
				colvalue[i] = (MGD77_igrf10syn (0, date, 1, 0.0, xvar[i], yvar[i], IGRF)) ? GMT_d_NaN : IGRF[MGD77_IGRF_F];
				n_sampled++;
			}
			if (gmtdefs.verbose) fprintf (stderr, "%s: Estimated IGRF at %ld locations out of %ld for cruise %s\n", GMT_program, n_sampled, D->H.n_records, list[argno]);
		}
#ifdef USE_CM4	
		else if (got_c == ADD_CM4) {	/* Append CM4 column */
			GMT_LONG ix, iy, it;
			double date, *xvar, *yvar, *tvar;
			struct MGD77_CM4 CM4;
			
			if ((ix = skip_if_missing ("lon",  list[argno], &In, D)) == MGD77_NOT_SET) continue;
			if ((iy = skip_if_missing ("lat",  list[argno], &In, D)) == MGD77_NOT_SET) continue;
			if ((it = skip_if_missing ("time", list[argno], &In, D)) == MGD77_NOT_SET) continue;

			xvar = (double *)D->values[ix];
			yvar = (double *)D->values[iy];
			tvar = (double *)D->values[it];
			colvalue = (double *)GMT_memory (VNULL, (size_t)D->H.n_records, sizeof (double), GMT_program);
			MGD77_CM4_init (&In, &CM4);
			
			for (i = n_sampled = 0; i < D->H.n_records; i++) {
				date = MGD77_time_to_fyear (&In, tvar[i]);	/* Get date as decimal year */
				colvalue[i] = MGD77_Calc_CM4 (date, xvar[i], yvar[i], FALSE, &CM4);
				n_sampled++;
			}
			MGD77_CM4_end (&CM4);
			if (gmtdefs.verbose) fprintf (stderr, "%s: Estimated CM4 at %ld locations out of %ld for cruise %s\n", GMT_program, n_sampled, D->H.n_records, list[argno]);
		}
		else if (got_c == ADD_RMAG4) {	/* Append recomputed residual mag column */
			GMT_LONG ix, iy, it, im;
			double date, *xvar, *yvar, *tvar, *mvar, IGRF[7];
			char field[5];
			
			if ((ix = skip_if_missing ("lon",  list[argno], &In, D)) == MGD77_NOT_SET) continue;
			if ((iy = skip_if_missing ("lat",  list[argno], &In, D)) == MGD77_NOT_SET) continue;
			if ((it = skip_if_missing ("time", list[argno], &In, D)) == MGD77_NOT_SET) continue;
			sprintf (field, "mtf%d", MTF_col);
			if ((im = skip_if_missing (field, list[argno], &In, D)) == MGD77_NOT_SET) continue;

			xvar = (double *)D->values[ix];
			yvar = (double *)D->values[iy];
			tvar = (double *)D->values[it];
			mvar = (double *)D->values[im];
			colvalue = (double *)GMT_memory (VNULL, (size_t)D->H.n_records, sizeof (double), GMT_program);
			
			for (i = n_sampled = 0; i < D->H.n_records; i++) {
				date = MGD77_time_to_fyear (&In, tvar[i]);	/* Get date as decimal year */
/* Change this--> */		check = MGD77_igrf10syn (0, date, 1, 0.0, xvar[i], yvar[i], IGRF);
				colvalue[i] = (check) ? GMT_d_NaN : mvar[i] - IGRF[MGD77_IGRF_F];
				n_sampled++;
			}
			if (gmtdefs.verbose) fprintf (stderr, "%s: Estimated recomputed magnetic anomaly at %ld locations out of %ld for cruise %s\n", GMT_program, n_sampled, D->H.n_records, list[argno]);
		}
#endif
		else if (got_c == ADD_GRAV) {	/* Append IGF column */
			GMT_LONG ix, iy, use;
			double *xvar, *yvar;
			
			if ((ix = skip_if_missing ("lon", list[argno], &In, D)) == MGD77_NOT_SET) continue;
			if ((iy = skip_if_missing ("lat", list[argno], &In, D)) == MGD77_NOT_SET) continue;

			if (GF_version == MGD77_NOT_SET) {
				use = (In.original) ? MGD77_ORIG : MGD77_REVISED;
				GF_version = D->H.mgd77[use]->Gravity_Theoretical_Formula_Code - '0';
				if (GF_version < MGD77_IGF_HEISKANEN || GF_version > MGD77_IGF_1980) {
					fprintf (stderr, "%s: Invalid Gravity Theoretical Formula Code (%c) - default to %d\n", GMT_program, D->H.mgd77[use]->Gravity_Theoretical_Formula_Code, MGD77_IGF_1980);
					GF_version = MGD77_IGF_1980;
				}
			}
			xvar = (double *)D->values[ix];
			yvar = (double *)D->values[iy];
			colvalue = (double *)GMT_memory (VNULL, (size_t)D->H.n_records, sizeof (double), GMT_program);
			
			for (i = 0; i < D->H.n_records; i++) colvalue[i] = MGD77_Theoretical_Gravity (xvar[i], yvar[i], (int)GF_version);
			if (gmtdefs.verbose) fprintf (stderr, "%s: Estimated IGRF at %ld locations out of %ld for cruise %s\n", GMT_program, D->H.n_records, D->H.n_records, list[argno]);
		}
		else if (got_c == ADD_CARTER) {	/* Append Carter correction column */
			GMT_LONG ix, iy, it;
			double *xvar, *yvar, *tvar;
			
			if ((ix = skip_if_missing ("lon", list[argno], &In, D)) == MGD77_NOT_SET) continue;
			if ((iy = skip_if_missing ("lat", list[argno], &In, D)) == MGD77_NOT_SET) continue;
			if ((it = skip_if_missing ("twt", list[argno], &In, D)) == MGD77_NOT_SET) continue;

			xvar = (double *)D->values[ix];
			yvar = (double *)D->values[iy];
			tvar = (double *)D->values[it];
			colvalue = (double *)GMT_memory (VNULL, (size_t)D->H.n_records, sizeof (double), GMT_program);
			
			for (i = 0; i < D->H.n_records; i++) colvalue[i] = MGD77_carter_correction (xvar[i], yvar[i], 1000.0 * tvar[i], &Carter);
			if (gmtdefs.verbose) fprintf (stderr, "%s: Estimated IGRF at %ld locations out of %ld for cruise %s\n", GMT_program, D->H.n_records, D->H.n_records, list[argno]);
		}
		else if (got_c == ADD_RMAG) {	/* Append recomputed residual mag column */
			GMT_LONG ix, iy, it, im;
			double date, *xvar, *yvar, *tvar, *mvar, IGRF[7];
			char field[5];
			
			if ((ix = skip_if_missing ("lon",  list[argno], &In, D)) == MGD77_NOT_SET) continue;
			if ((iy = skip_if_missing ("lat",  list[argno], &In, D)) == MGD77_NOT_SET) continue;
			if ((it = skip_if_missing ("time", list[argno], &In, D)) == MGD77_NOT_SET) continue;
			sprintf (field, "mtf%ld", MTF_col);
			if ((im = skip_if_missing (field, list[argno], &In, D)) == MGD77_NOT_SET) continue;

			xvar = (double *)D->values[ix];
			yvar = (double *)D->values[iy];
			tvar = (double *)D->values[it];
			mvar = (double *)D->values[im];
			colvalue = (double *)GMT_memory (VNULL, (size_t)D->H.n_records, sizeof (double), GMT_program);
			
			for (i = n_sampled = 0; i < D->H.n_records; i++) {
				date = MGD77_time_to_fyear (&In, tvar[i]);	/* Get date as decimal year */
				check = MGD77_igrf10syn (0, date, 1, 0.0, xvar[i], yvar[i], IGRF);
				colvalue[i] = (check) ? GMT_d_NaN : mvar[i] - IGRF[MGD77_IGRF_F];
				n_sampled++;
			}
			if (gmtdefs.verbose) fprintf (stderr, "%s: Estimated recomputed magnetic anomaly at %ld locations out of %ld for cruise %s\n", GMT_program, n_sampled, D->H.n_records, list[argno]);
		}
		else if (got_grid) {	/* Sample grid along track (or Mercator-projected) track */
			GMT_LONG ix, iy;
			double *xvar, *yvar;
			
			if ((ix = skip_if_missing ("lon", list[argno], &In, D)) == MGD77_NOT_SET) continue;
			if ((iy = skip_if_missing ("lat", list[argno], &In, D)) == MGD77_NOT_SET) continue;

			xvar = (double *)D->values[ix];
			yvar = (double *)D->values[iy];
			colvalue = (double *)GMT_memory (VNULL, (size_t)D->H.n_records, sizeof (double), GMT_program);
			
			for (i = n_sampled = 0; i < D->H.n_records; i++) {
				colvalue[i] = GMT_d_NaN;	/* In case we are outside grid */
	
				/* If point is outside grd area, shift it using periodicity or skip if not periodic. */

				if (got_i)	/* Mercator IMG grid */
					GMT_geo_to_xy (xvar[i], yvar[i], &x, &y);
				else {		/* Regular geographic grd */
					x = xvar[i];
					y = yvar[i];
				}
				if (y < grd.y_min || y > grd.y_max) continue;

				while ( (x < grd.x_min) && (edgeinfo.nxp > 0) ) x += (grd.x_inc * edgeinfo.nxp);
				if (x < grd.x_min) continue;

				while ( (x > grd.x_max) && (edgeinfo.nxp > 0) ) x -= (grd.x_inc * edgeinfo.nxp);
				if (x > grd.x_max) continue;

				if (interpolate) {	/* IMG has been corrected, and GRD is good to go */
					colvalue[i] = GMT_get_bcr_z (&grd, x, y, f, &edgeinfo, &bcr);
				}
				else {	/* Take IMG nearest node and do special stuff (values already set during read) */
					ii = GMT_x_to_i (x, grd.x_min, grd.x_inc, grd.xy_off, grd.nx);
					jj = GMT_y_to_j (y, grd.y_min, grd.y_inc, grd.xy_off, grd.ny);
					colvalue[i] = f[(jj+GMT_pad[3])*mx+ii+GMT_pad[0]];
				}
				n_sampled++;
			}
			if (gmtdefs.verbose) fprintf (stderr, "%s: Sampled grid at %ld locations out of %ld for cruise %s\n", GMT_program, n_sampled, D->H.n_records, list[argno]);
		}
		else if (got_a) {	/* Just got a single column to paste in, assuming the row numbers match */
			if (n != D->H.n_records) {
				fprintf (stderr, "%s: Extra column data records (%ld) do not match # of cruise records (%ld) for %s\n", GMT_program, n, D->H.n_records, list[argno]);
				exit (EXIT_FAILURE);
			}
			if (gmtdefs.verbose) fprintf (stderr, "%s: Appended column data for all %ld records for cruise %s\n", GMT_program, D->H.n_records, list[argno]);
		}
		else if (got_d || got_n || got_t) {	/* Got either (time,data) or (dist,data) */
			GMT_LONG ix, iy, it;
			double *x = VNULL, *y = VNULL, *d = VNULL;
			colvalue = (double *) GMT_memory ((void *)colvalue, (size_t)D->H.n_records, sizeof (double), GMT_program);
			if (got_d) {	/* Must create distances in user's units */
				if ((ix = skip_if_missing ("lon", list[argno], &In, D)) == MGD77_NOT_SET) continue;
				if ((iy = skip_if_missing ("lat", list[argno], &In, D)) == MGD77_NOT_SET) continue;
				x = (double *)D->values[ix];
				y = (double *)D->values[iy];
				GMT_err_fail (GMT_distances (x, y, D->H.n_records, dist_scale, dist_flag, &d), "");
				x = d;
			}
			else if (got_t) {	/* Time */
				if ((it = skip_if_missing ("time", list[argno], &In, D)) == MGD77_NOT_SET) continue;
				x = (double *)D->values[it];
			}
			if (interpolate) {	/* Using given table to interpolate the values at all mgd77 records */
				y = (double *) GMT_memory (VNULL, (size_t)D->H.n_records, sizeof (double), GMT_program);
				result = GMT_intpol (coldnt, colvalue, n, D->H.n_records, x, y, gmtdefs.interpolant);
				if (result != 0) {
					fprintf (stderr, "%s: Error from GMT_intpol near row %ld!\n", GMT_program, result+1);
					exit (EXIT_FAILURE);
				}
				memcpy ((void *)colvalue, (void *)y, (size_t)(D->H.n_records * sizeof (double)));
				GMT_free ((void *)y);
			}
			else if (strings && n < D->H.n_records) {	/* Only update the exact matching records */
				text = (char *) GMT_memory (VNULL, (size_t)(D->H.n_records * LEN), sizeof (char), GMT_program);
				for (i = j = n_sampled = 0; i < D->H.n_records && j < n; i++) {
					match_value = (got_n) ? i+1 : x[i];
					strncpy (&text[i*LEN], not_given, (size_t)LEN);	/* In case we have no data at this time */
					while (coldnt[j] < match_value && j < n) j++;
					if (coldnt[j] == match_value) {
						strncpy (&text[i*LEN], tmp_string[j], (size_t)LEN);
						n_sampled++;
					}
				}
				GMT_free ((void *)tmp_string);
				if (gmtdefs.verbose) fprintf (stderr, "%s: Appended column data for %ld locations out of %ld for cruise %s\n", GMT_program, n_sampled, D->H.n_records, list[argno]);
			}
			else if (strings) {	/* One to one match */
				text = (char *) GMT_memory (VNULL, (size_t)(D->H.n_records * LEN), sizeof (char), GMT_program);
				for (i = 0; i < n; i++) strncpy (&text[i*LEN], tmp_string[i], (size_t)LEN);
				GMT_free ((void *)tmp_string);
				if (gmtdefs.verbose) fprintf (stderr, "%s: Appended column data for %ld locations out of %ld for cruise %s\n", GMT_program, n_sampled, D->H.n_records, list[argno]);
			}
			else {	/* Only update the exact matching records */
				y = (double *) GMT_memory (VNULL, (size_t)D->H.n_records, sizeof (double), GMT_program);
				for (i = j = n_sampled = 0; i < D->H.n_records && j < n; i++) {
					match_value = (got_n) ? i+1 : x[i];
					y[i] = GMT_d_NaN;	/* In case we have no data at this time */
					while (coldnt[j] < match_value && j < n) j++;
					if (coldnt[j] == match_value) {	/* Found our guy */
						y[i] = colvalue[j];
						n_sampled++;
					}
				}
				memcpy ((void *)colvalue, (void *)y, (size_t)(D->H.n_records * sizeof (double)));
				if (gmtdefs.verbose) fprintf (stderr, "%s: Appended column data for %ld locations out of %ld for cruise %s\n", GMT_program, n_sampled, D->H.n_records, list[argno]);
				GMT_free ((void *)y);
			}
			if (got_d) GMT_free ((void *)d);
		}
		else if (got_e)
		{
			/* Read any header scale/offset factors to be set.
			 * Decode error flags to give bit flag, store in colvalue
			 */
			FILE *fp_e;
			int cdf_var_id, cdf_adjust;
			char ID[16], date[16], field[GMT_TEXT_LEN], efile[BUFSIZ], E77[256], timestamp[GMT_TEXT_LEN], answer[BUFSIZ], code[BUFSIZ], YorN, kind;
			GMT_LONG n_recs, rec, number, tz;
			GMT_LONG type, it, id, key, n_E77_flags, from, to, day, month, year, item;
			GMT_LONG n_E77_headers, n_E77_scales, n_E77_offsets, n_E77_recalcs, n_tz_corr, n_unprocessed, e_error = 0;
			unsigned int *flags, pattern;
			size_t length;
			GMT_LONG has_time, old_flags, tz_errors = FALSE;
			short *tz_corr = NULL;
			struct MGD77_HEADER_PARAMS *P;
			double rec_time, del_t, value, *tvar = VNULL;
			
			if (D->H.E77 && strlen(D->H.E77) > 0 && !replace) {
				fprintf (stderr, "%s: E77 corrections are already present in %s.  use -A+e to overwrite with new corrections\n", GMT_program, list[argno]);
				MGD77_Free (D);	/* Free memory allocated by MGD77_Read_File for this aborted effort */
				continue;
			}
			
			sprintf (efile, "%s.e77", list[argno]);
			if ((fp_e = GMT_fopen (efile, "r")) == NULL) {	/* Not in current directory, try MGD77_HOME/E77 */
				sprintf (efile, "%s/E77/%s.e77", In.MGD77_HOME, list[argno]);
				if ((fp_e = GMT_fopen (efile, "r")) == NULL) {	/* Not here either */
					fprintf(stderr, "%s: ERROR: The file %s.e77 could not be found in current directory or in MGD77_HOME/E77 - skipped\n", GMT_program, list[argno]);
					MGD77_Free (D);	/* Free memory allocated by MGD77_Read_File for this aborted effort */
					continue;
				}
			
			}
			/* We will first do many checks to make sure this E77 file goes with the specified cruise and that
			 * all the verification steps has been taken to make this E77 a correction file
			 */
			 
			P = D->H.mgd77[MGD77_ORIG];	/* Because E77 is absolute and not incremental we start from original settings */
			if (!GMT_fgets (line, BUFSIZ, fp_e)) {
				fprintf(stderr, "%s: ERROR: Could not read record #1 from %s.e77 - aborting\n", GMT_program, list[argno]);
				e_error++;
			}
			sscanf (&line[1], "%*s %s %*s %*s %*s %*s %*s %s %*s %" GMT_LL "d", ID, date, &n_recs);
			if (strcmp (In.NGDC_id, ID)) {
				fprintf(stderr, "%s: ERROR: E77 Conflict %s : ID = %s versus %s - aborting\n", GMT_program, efile, ID, In.NGDC_id);
				e_error++;
			}
			/* Make sure the File creation dates from the data file and the E77 match */
			day = atoi (&date[6]);
			date[6] = 0;
			month = atoi (&date[4]);
			date[4] = 0;
			year = atoi (date);
			
			if (!(year == atoi (P->File_Creation_Year) && month == atoi (P->File_Creation_Month) && day == atoi (P->File_Creation_Day))) {
				fprintf(stderr, "%s: ERROR: E77 Conflict %s: File Creation Date: %s versus %s%s%s - aborting\n", GMT_program, efile, date,
					P->File_Creation_Year, P->File_Creation_Month, P->File_Creation_Day);
				e_error++;
			}
			if (n_recs != D->H.n_records) {
				fprintf(stderr, "%s: ERROR: E77 Conflict %s: n_recs = %ld versus %ld = aborting\n", GMT_program, efile, n_recs, D->H.n_records);
				e_error++;
			}
			verified = FALSE;
			while (GMT_fgets (line, BUFSIZ, fp_e) && strncmp (line, "# Errata: Header", (size_t)14)) {	/* Read until we get to Header record section */
				if (line[0] == '#') continue;	/* Skip comments */
				GMT_chop (line);		/* Rid the world of CR/LF */
				if (!strncmp (line, "Y Errata table verification status", (size_t)34)) verified = TRUE;
			}
			if (!verified && !ignore_verify) {
				fprintf(stderr, "%s: ERROR: E77 file %s not yet verified.  E77 not applied\n", GMT_program, efile);
				e_error++;
			}
			
			if (e_error) {
				fprintf(stderr, "%s: ERROR: The file %s has too many errors.  E77 not applied\n", GMT_program, efile);
				MGD77_Free (D);	/* Free memory allocated by MGD77_Read_File for this aborted effort */
				continue;
			}
			
			/* OK, we got this far so the meta data for this E77 file seems OK */
			
			/* Quickly scan through file to make sure there are no unprocessed recommendations or bad records before making changes */
			
			e_error = n_unprocessed = 0;
			while (GMT_fgets (line, BUFSIZ, fp_e)) {
				GMT_chop (line);		/* Rid the world of CR/LF */
				if (line[0] == '#' || line[0] == '\0') continue;	/* Comments or blank lines are OK */
				if (line[1] == '-') {		/* Header record */
					if (!(line[0] == 'Y' || line[0] == 'N') && !ignore_verify) {		/* Unprocessed recommendation? */
						fprintf (stderr, "%s: UNDECIDED: %s\n", list[argno], line);
						if (line[0] == '?') n_unprocessed++;
						e_error++;
					}
					sscanf (line, "%c-%c-%[^-]-%[^-]-%" GMT_LL "d", &YorN, &kind, ID, field, &item);
				}
				else				/* Data record */
					sscanf (line, "%c %s %s %" GMT_LL "d %s", &YorN, ID, timestamp, &rec, code);
				if (strcmp (In.NGDC_id, ID)) {
					fprintf(stderr, "%s: ERROR: E77 Conflict %s : ID = %s versus %s in header records!\n", GMT_program, efile, ID, In.NGDC_id);
					e_error++;
				}
			}
			
			if (e_error) {
				fprintf(stderr, "%s: ERROR: The file %s has too many errors.  E77 not applied\n", GMT_program, efile);
				GMT_fclose (fp_e);
				MGD77_Free (D);	/* Free memory allocated by MGD77_Read_File for this aborted effort */
				continue;
			}
			if (n_unprocessed) {
				fprintf(stderr, "%s: ERROR: The file %s has unprocessed E77 recommendations.  E77 not applied\n", GMT_program, efile);
				GMT_fclose (fp_e);
				MGD77_Free (D);	/* Free memory allocated by MGD77_Read_File for this aborted effort */
				continue;
			}
			
			/* OK, here we believe the E77 file contains the correct information for this cruise. Rewind and start from top */
			
			GMT_rewind (fp_e);
			while (GMT_fgets (line, BUFSIZ, fp_e) && strncmp (line, "# Errata: Header", (size_t)14));	/* Read until we get to Header record section */
			
			flags = (unsigned int *)GMT_memory (VNULL, (size_t)D->H.n_records, sizeof (unsigned int), GMT_program);
			n_E77_flags = n_E77_headers = n_E77_scales = n_E77_offsets = n_E77_recalcs = n_tz_corr = 0;

			MGD77_nc_status (nc_open (In.path, NC_WRITE, &In.nc_id));	/* Open the file */
			MGD77_nc_status (nc_redef (In.nc_id));				/* Enter define mode */
			old_flags = MGD77_Remove_E77 (&In);				/* Remove any previously revised header parameters */
			while (GMT_fgets (line, BUFSIZ, fp_e) && strncmp (line, "# Errata: Data", (size_t)14)) {	/* Read until we get to data record section */
				GMT_chop (line);					/* Rid the world of CR/LF */
				if (line[0] == '#' || line[0] == '\0') continue;	/* Skip comments */
				/* Example of expected line:  
				   Y-E-06050010-H15-01: Invalid Gravity Departure Base Station Value: (0000000) [1000009]
				*/
				sscanf (line, "%c-%c-%[^-]-%[^-]-%" GMT_LL "d", &YorN, &kind, ID, field, &item);
				if (strcmp (In.NGDC_id, ID)) {
					fprintf(stderr, "%s: ERROR: E77 Conflict %s : ID = %s versus %s in header records - skipped\n", GMT_program, efile, ID, In.NGDC_id);
					e_error++;
					continue;
				}
				if (field[0] == 'H') {
					type = E77_HEADER_MODE;
					number = (GMT_LONG)atoi (&field[1]);
				}
				else {
					type = 1;
					number = (GMT_LONG)item;
				}
				if (e77_skip_mode[type]) continue;
				if (!e77_skip_mode[type] && YorN == 'N') continue;
				if (kind == 'W') {	/* Output the warning (if Y) and goto next line*/
					if (gmtdefs.verbose == 2 && (YorN == 'Y' || (ignore_verify && YorN == '?'))) fprintf (stderr, "%s: WARNING: %s\n", list[argno], line);
					continue;
				}
				if (!got_default_answer (line, answer)) continue;
				
				/* Here we must do something */
				
				if (type == E77_HEADER_MODE) {	/* Header meta data fixes */
				
					key = MGD77_Param_Key ((int)number, (int)item);	/* Returns -ve if sequence not found or item not found, >=0 otherwise */
					
					switch (key) {
						case MGD77_BAD_HEADER_RECNO:
							fprintf (stderr, "%s: Warning: Sequence number %ld is outside range - skipped\n", GMT_program, number);
							break;
						case MGD77_BAD_HEADER_ITEM:
							fprintf (stderr, "%s: Warning: Sequence number %ld, Item %ld is not supported - skipped\n", GMT_program, number, item);
							break;
						default:	/* Update internal structure as well as netCDF file attributes */
							length = (MGD77_Header_Lookup[key].length == 1) ? 1 : strlen (answer);
							strncpy (MGD77_Header_Lookup[key].ptr[MGD77_REVISED], answer, (size_t)length);
							MGD77_Put_Param (&In, MGD77_Header_Lookup[key].name, length, MGD77_Header_Lookup[key].ptr[MGD77_ORIG], length, MGD77_Header_Lookup[key].ptr[MGD77_REVISED], 2);
							n_E77_headers++;
							break;
					}
				}
				else {			/* Systematic fixes */
					if ((id = MGD77_Get_Column (field, &In)) == MGD77_NOT_SET) {
						fprintf (stderr, "%s: Warning: Correction found for %s which is not in this cruise?\n", GMT_program, field);
					}
					else {
						k = MGD77_Info_from_Abbrev (field, &(D->H), &set, &item);
						value = atof (answer);
						switch (number) {
							case E77_HDR_PDR:	/* Must deal with undetected Precision Depth Recorder wrap-arounds - this also force recalc of depth when data is read*/
								MGD77_nc_status (nc_put_att_double (In.nc_id, NC_GLOBAL, "PDR_wrap", NC_DOUBLE, (size_t)1, &value));
								cdf_adjust = MGD77_COL_ADJ_TWT;
								MGD77_nc_status (nc_put_att_int (In.nc_id, D->H.info[set].col[id].var_id, "adjust", NC_INT, (size_t)1, &cdf_adjust));
								n_E77_recalcs++;
								if ((id = MGD77_Get_Column ("depth", &In)) == MGD77_NOT_SET) {
									fprintf (stderr, "%s: Warning: Correction implied for %s which is not in this cruise?\n", GMT_program, field);
									break;
								}
								/* no break - we want to fall through and also set depth adjustment */
							case E77_HDR_CARTER:	/* Recalculate Carter depth from twt */
								cdf_adjust = MGD77_COL_ADJ_DEPTH;
								MGD77_nc_status (nc_put_att_int (In.nc_id, D->H.info[set].col[id].var_id, "adjust", NC_INT, (size_t)1, &cdf_adjust));
								n_E77_recalcs++;
								break;
							case E77_HDR_ANOM_MAG:	/* Recalculate anomaly mag as mtf1 - igrf */
								cdf_adjust = MGD77_COL_ADJ_MAG;
								MGD77_nc_status (nc_put_att_int (In.nc_id, D->H.info[set].col[id].var_id, "adjust", NC_INT, (size_t)1, &cdf_adjust));
								n_E77_recalcs++;
								break;
							case E77_HDR_ANOM_FAA:	/* Recalculate anomaly faa as gobs - igf */
								cdf_adjust = MGD77_COL_ADJ_FAA;
								MGD77_nc_status (nc_put_att_int (In.nc_id, D->H.info[set].col[id].var_id, "adjust", NC_INT, (size_t)1, &cdf_adjust));
								n_E77_recalcs++;
								break;
							case E77_HDR_ANOM_FAA_EOT:	/* Recalculate anomaly faa as gobs - igf + eot */
								cdf_adjust = MGD77_COL_ADJ_FAA_EOT;
								MGD77_nc_status (nc_put_att_int (In.nc_id, D->H.info[set].col[id].var_id, "adjust", NC_INT, (size_t)1, &cdf_adjust));
								n_E77_recalcs++;
								break;
							case E77_HDR_SCALE:	/* Correction scale factor */
								if (D->H.info[set].col[id].corr_factor == 1.0) {	/* Must add a new attribute to the file */
									D->H.info[set].col[id].corr_factor = value;
									MGD77_nc_status (nc_put_att_double (In.nc_id, D->H.info[set].col[id].var_id, "corr_factor", NC_DOUBLE, (size_t)1, &D->H.info[set].col[id].corr_factor));
								}
								n_E77_scales++;
								break;
							case E77_HDR_DCSHIFT:	/* Correction offset */
								if (D->H.info[set].col[id].corr_offset == 0.0) {	/* Must add a new attribute to the file */
									D->H.info[set].col[id].corr_offset = value;
									MGD77_nc_status (nc_put_att_double (In.nc_id, D->H.info[set].col[id].var_id, "corr_offset", NC_DOUBLE, (size_t)1, &D->H.info[set].col[id].corr_offset));
								}
								n_E77_offsets++;
								break;
							case E77_HDR_GRID_OFFSET:	/* Range of bad values - set flags to BAD  */
							case E77_HDR_FLAGRANGE:		/* Range of bad values - set flags to BAD  */
								sscanf (answer, "%" GMT_LL "d-%" GMT_LL "d", &from, &to);
								if (from < 1 || from > D->H.n_records || to < 1 || to > D->H.n_records || to < from) {
									fprintf (stderr, "%s: Error: Record range %s is invalid.  Correction skipped\n", GMT_program, answer);
									break;
								}
								pattern = set_bit (id);
 								for (rec = from-1; rec < to; rec++, n_E77_flags++) flags[rec] |= pattern;	/* -1 to get C indices */
								break;
							default:
								break;
						}
					}
				}
			}
			/* Now start on data record section */
			has_time = TRUE;
			if ((it = skip_if_missing ("time", list[argno], &In, D)) == MGD77_NOT_SET)
				has_time = FALSE;
			else {	/* See if we really have time or if they are all NaN */
				tvar = (double *)D->values[it];
				for (rec = 0, has_time = FALSE; !has_time && rec < D->H.n_records; rec++) if (!GMT_is_dnan (tvar[rec])) has_time = TRUE;
			}
			while (GMT_fgets (line, BUFSIZ, fp_e)) {	/* Read until EOF */
				sscanf (line, "%c %s %s %" GMT_LL "d %s", &YorN, ID, timestamp, &rec, code);
				if (strcmp (In.NGDC_id, ID)) {
					fprintf(stderr, "%s: ERROR: E77 Conflict %s : ID = %s versus %s in data records - skipped\n", GMT_program, efile, ID, In.NGDC_id);
					e_error++;
					continue;
				}
				if (YorN == 'N') continue;			/* Already decided NOT to use this correction */
				if (YorN == '?' && !ignore_verify) {		/* Undecided: Output the warning and goto next line unless we ignore verification */
					fprintf (stderr, "%s: UNDECIDED: %s\n", list[argno], line);
					continue;
				}
				/* Here, YorN is 'Y' (or '?' if ignore_verify is TRUE) */
				rec--;	/* E77 starts with rec = 1 for first data record */
				if (has_time) {
					if (!strcmp(timestamp,"NaN")) {
						if (gmtdefs.verbose == 2) fprintf (stderr, "%s: WARNING: %s: E77 time stamp %s, using recno\n", GMT_program, ID, timestamp);
					}
					else {	/* Must try to interpret the timestamp */
						if (GMT_verify_expectations (GMT_IS_ABSTIME, GMT_scanf (timestamp, GMT_IS_ABSTIME, &rec_time), timestamp)) {
							fprintf (stderr, "%s: ERROR: %s: E77 time stamp (%s) in wrong format? - skipped\n", GMT_program, ID, timestamp);
							continue;
						}
						del_t = fabs (tvar[rec] - rec_time);
						if (del_t > (0.06 + GMT_CONV_LIMIT)) {	/* 0.06 is finest time step in MGD77 file so we allow that much slop */
							fprintf (stderr, "%s: ERROR: %s: E77 time stamp and record number do not match record time (del_t = %g s) - skipped\n", GMT_program, ID, del_t);
							continue;
						}
					}
				}
				pos = 0;
				item = -1;	/* So we increment to 0 inside the loop */
				while (GMT_strtok (code, "-", &pos, p)) {
					item++;
					if (e77_skip_mode[item+2]) continue;	/* Ignore this sort of code */
					if (p[0] == '0') continue;
					for (k = 0; k < (int)strlen(p); k++) {	/* Loop over one or more codes */
					
						if (item == 0) {	/* NAV */
							switch (p[k]) {
								case 'A':	/* Time out of range */
									flags[rec] |= set_bit(NCPOS_TIME);
									n_E77_flags++;
									break;
								case 'B':
									if (gmtdefs.verbose == 2) fprintf (stderr, "%s: Decreasing time %s - Source Institution need to sort records\n", list[argno], timestamp);
									break;
								case 'C':	/* Excessive speed - flag time, lon, lat */
									flags[rec] |= set_bit(NCPOS_TIME);
									flags[rec] |= set_bit(NCPOS_LON);
									flags[rec] |= set_bit(NCPOS_LAT);
									n_E77_flags++;
									break;
								case 'D':	/* On land */ 
									flags[rec] |= set_bit(NCPOS_LON);
									flags[rec] |= set_bit(NCPOS_LAT);
									n_E77_flags++;
									break;
								case 'E':	/* Undefined nav - flag time, lon, lat */
									flags[rec] |= set_bit(NCPOS_TIME);
									flags[rec] |= set_bit(NCPOS_LON);
									flags[rec] |= set_bit(NCPOS_LAT);
									n_E77_flags++;
									break;
								case 'F':	/* Time-zone error */
									if (!tz_errors) {	/* First time we must allocate tz_corr array */
										tz_corr = (short *)GMT_memory (VNULL, (size_t)D->H.n_records, sizeof (short), GMT_program);
										tz_errors = TRUE;
									}
									sscanf (line, "%*c %*s %*s %*s %*s %*s %*s %*s %" GMT_LL "d", &tz);
									tz_corr[rec] = (short)tz;
									n_tz_corr++;
									break;
								default:
									fprintf (stderr, "%s: Unrecognized NAV code %c - skipped\n", list[argno], p[k]);
									break;
							}
						}
						else if (p[k] < 'A' || p[k] > 'X') {
							fprintf (stderr, "%s: Unrecognized error field %c - skipped\n", list[argno], p[k]);
						}
						else {			/* EO, RANGE, or SLOPE */
							if (p[k] >= 'A' && p[k] <= 'X')	{ /* Valid codes */
								if (p[k] > 'G')
 									key = p[k] - 'A' - 4;	/* H (lat) = 3, J (ptc) = 5, etc (only expect J-X though) */
								else if (p[k] < 'C')
									key = p[k] - 'A' + 1;	/* A (rectype) = 1, B (TZ) = 2 */
								else
  									key = 0;	/* C-G (yyyy,mm,dd,hh,mi) all map to time 0 */
								flags[rec] |= set_bit(key);
								n_E77_flags++;
							}
							else {
								fprintf (stderr, "%s: Unrecognized error field %c - skipped\n", list[argno], p[k]);
							}
						}
					}		
				}
			}
			GMT_fclose (fp_e);

			/* Update E77 history */

			(void) time (&now);
			sprintf (E77, "%s [%s] E77 corrections applied to header: %ld scale: %ld offset: %ld recalc: %ld flags: %ld tz_corr: %ld", ctime(&now), In.user, n_E77_headers, n_E77_scales, n_E77_offsets, n_E77_recalcs, n_E77_flags, n_tz_corr);
			for (i = 0; E77[i]; i++) if (E77[i] == '\n') E77[i] = ' ';	/* Remove the \n returned by ctime() */
			k = strlen (E77);
			D->H.E77 = (char *)GMT_memory ((void *)D->H.E77, (size_t)(k+1), sizeof (char), GMT_program);
			strcpy (D->H.E77, E77);
			MGD77_nc_status (nc_put_att_text (In.nc_id, NC_GLOBAL, "E77", (size_t)k, D->H.E77));
		
			if (tz_errors) {	/* Use time var_id */
				MGD77_nc_status (nc_put_att_short (In.nc_id, D->H.info[MGD77_M77_SET].col[0].var_id, "tz_corr", NC_SHORT, (size_t)D->H.n_records, tz_corr));
				memset ((void *)answer, 0, (size_t)BUFSIZ);	/* No default answer */
				strcpy (answer, "tz_corr in hours to subtract from UTC");
				MGD77_nc_status (nc_put_att_text (In.nc_id, D->H.info[MGD77_M77_SET].col[0].var_id, "comment", strlen (answer), answer));
				GMT_free ((void *)tz_corr);
			}
			old_flags =  (nc_inq_varid (In.nc_id, "MGD77_flags", &cdf_var_id) == NC_NOERR);	/* TRUE if flag variable exists already */
			
			if (n_E77_flags) {	/* Add flags to netCDF file */
				if (old_flags) {	/* Flag variable exists already - simply replace existing flags with the new ones */
					if (D->flags[0])	/* Was allocated and read */
						memcpy ((void *)D->flags[0], (void *)flags, (size_t)(D->H.n_records * sizeof (int)));
					else	/* Was not allcoated */
						D->flags[0] = flags;
				}
				else {	/* We need to define the flags for the first time */
					dims[0] = In.nc_recid;
					MGD77_nc_status (nc_def_var (In.nc_id, "MGD77_flags", NC_INT, 1, dims, &cdf_var_id));	/* Define an array variable */
					memset ((void *)answer, 0, (size_t)BUFSIZ);	/* No default answer */
					strcpy (answer, "MGD77 flags (ON = Bad, OFF = Good) derived from E77 errata");
					MGD77_nc_status (nc_put_att_text (In.nc_id, cdf_var_id, "comment", strlen (answer), answer));
					D->flags[0] = flags;
				}
				MGD77_nc_status (nc_enddef (In.nc_id));	/* End define mode. */
				start[0] = 0;
				count[0] = D->H.n_records;
				MGD77_nc_status (nc_put_vara_int (In.nc_id, cdf_var_id, start, count, (int *)D->flags[0]));
			}
			else if (old_flags) {	/* Had flags from before which we cannot delete */
				MGD77_nc_status (nc_enddef (In.nc_id));	/* End define mode. */
				fprintf (stderr, "%s: File %s contains flags from an earlier E77 but this E77 do not contain any flags.\n", GMT_program, list[argno]);
				fprintf (stderr, "%s: The flags in the file %s will all be set to zero but cannot be removed.\n", GMT_program, list[argno]);
				fprintf (stderr, "%s: If possible, recreate the MGD77+ file %s from the MGD77 original, then reapply E77.\n", GMT_program, list[argno]);
				start[0] = 0;
				count[0] = D->H.n_records;
				memset ((void *)D->flags[0], 0, (size_t)(D->H.n_records * sizeof (int)));	/* Reset all flags to 0 (GOOD) */
				MGD77_nc_status (nc_put_vara_int (In.nc_id, cdf_var_id, start, count, (int *)D->flags[0]));
			}
			
			MGD77_Free (D);
			MGD77_Close_File (&In);
			n_changed++;
			continue;	/* Nothing more to do for this file */
		}
		
		/* Specify the information for the extra column. */
		
		constant = (LEN == 0) ? MGD77_dbl_are_constant (colvalue, D->H.n_records, limits) : MGD77_txt_are_constant (text, D->H.n_records, (int)LEN);	/* Do we need to store 1 or n values? */

		if (column != MGD77_NOT_SET) {	/* Is it possible just to replace the existing column? */
			error = 0;
			if (LEN) {
				if (OLDLEN != LEN) {
					fprintf (stderr, "%s: Revised text column %s differs in width (%d) from the old values (%d).\n", GMT_program, c_abbrev, (int)LEN, (int)OLDLEN);
					error = TRUE;
				}
				if (constant && n_dims == 2) {
					fprintf (stderr, "%s: Revised text column %s is constant whereas old values were in an array\n", GMT_program, c_abbrev);
					error = TRUE;
				}
				if (!constant && n_dims == 1) {
					fprintf (stderr, "%s: Revised text column %s is an array whereas old values is a constant\n", GMT_program, c_abbrev);
					error = TRUE;
				}
			}
			else {
				if (constant && n_dims == 1) {
					fprintf (stderr, "%s: Revised data column %s is constant whereas old values were in an array\n", GMT_program, c_abbrev);
					error = TRUE;
				}
				if (!constant && n_dims == 0) {
					fprintf (stderr, "%s: Revised data column %s is an array whereas old values is a constant\n", GMT_program, c_abbrev);
					error = TRUE;
				}
			}
			if (error) {
				fprintf (stderr, "%s: You must use -D to delete the old information before adding the new information\n", GMT_program);
				continue;
			}
		}
		
		/* OK, here we may either replace an exiting column or add a new one */
		
		if (MGD77_Open_File (list[argno], &In, MGD77_WRITE_MODE)) return (-1);	/* Only creates the full path to the new file */
	
		MGD77_nc_status (nc_open (In.path, NC_WRITE, &In.nc_id));	/* Open the file */
		MGD77_nc_status (nc_redef (In.nc_id));				/* Enter define mode */
		
		dims[0] = In.nc_recid;	dims[1] = LEN;
		start[0] = start[1] = 0;
		count[0] = D->H.n_records;	count[1] = LEN;
		
		if (column == MGD77_NOT_SET) {	/*Adding a new column */
			if (constant) {	/* Simply store one value */
				if (LEN)	/* Text variable */
					MGD77_nc_status (nc_def_var (In.nc_id, c_abbrev, c_nc_type, 1, &dims[1], &cdf_var_id));	/* Define a single text variable */
				else
					MGD77_nc_status (nc_def_var (In.nc_id, c_abbrev, c_nc_type, 0, NULL, &cdf_var_id));		/* Define a single variable */
			}
			else {	/* Must store array */
				if (LEN)	/* Text variable */
					MGD77_nc_status (nc_def_var (In.nc_id, c_abbrev, c_nc_type, 2, dims, &cdf_var_id));		/* Define a 2-D text variable */
				else
					MGD77_nc_status (nc_def_var (In.nc_id, c_abbrev, c_nc_type, 1, dims, &cdf_var_id));		/* Define a number array variable */
			}
		}
		else	/* Reuse, get id */
			MGD77_nc_status (nc_inq_varid (In.nc_id, c_abbrev, &cdf_var_id));
		
		if (c_name[0]) MGD77_nc_status (nc_put_att_text   (In.nc_id, cdf_var_id, "long_name", strlen (c_name), c_name));
		if (c_units[0]) MGD77_nc_status (nc_put_att_text   (In.nc_id, cdf_var_id, "units", strlen (c_units), c_units));
		MGD77_nc_status (nc_put_att_double   (In.nc_id, cdf_var_id, "actual_range", NC_DOUBLE, (size_t)2, limits));
		if (c_comment[0]) MGD77_nc_status (nc_put_att_text   (In.nc_id, cdf_var_id, "comment", strlen (c_comment), c_comment));
		MGD77_nc_status (nc_put_att_double (In.nc_id, cdf_var_id, "_FillValue", c_nc_type, (size_t)1, &MGD77_NaN_val[c_nc_type]));
		MGD77_nc_status (nc_put_att_double (In.nc_id, cdf_var_id, "missing_value", c_nc_type, (size_t)1, &MGD77_NaN_val[c_nc_type]));
		if (parameters[COL_SCALE]  != 1.0) MGD77_nc_status (nc_put_att_double (In.nc_id, cdf_var_id, "scale_factor", NC_DOUBLE, (size_t)1, &parameters[COL_SCALE]));
		if (parameters[COL_OFFSET] != 0.0) MGD77_nc_status (nc_put_att_double (In.nc_id, cdf_var_id, "add_offset",   NC_DOUBLE, (size_t)1, &parameters[COL_OFFSET]));
					
		/* Update history */

		(void) time (&now);
		sprintf (history, "%s [%s] Column %s added", ctime(&now), In.user, c_abbrev);
		k = strlen (history);
		for (i = 0; i < k; i++) if (history[i] == '\n') history[i] = ' ';	/* Remove the \n returned by ctime() */
		history[k++] = '\n';	history[k] = '\0';				/* Add LF at end of line */
		k += strlen (D->H.history);
		D->H.history = (char *)GMT_memory ((void *)D->H.history, (size_t)k, sizeof (char), GMT_program);
		strcat (D->H.history, history);
		MGD77_nc_status (nc_put_att_text (In.nc_id, NC_GLOBAL, "history", strlen (D->H.history), D->H.history));
		
		MGD77_nc_status (nc_enddef (In.nc_id));	/* End define mode.  Now we can write/update data */

		transform = (! (parameters[COL_SCALE] == 1.0 && parameters[COL_OFFSET] == 0.0));	/* TRUE if we must transform before writing */
		n_bad = 0;
		if (constant) {	/* Simply store one value */
			if (LEN)
				MGD77_nc_status (nc_put_vara_schar (In.nc_id, cdf_var_id, start, &count[1], (signed char *)text));	/* Just write one text string */
			else {
				n_bad = MGD77_do_scale_offset_before_write (&single_val, colvalue, (size_t)1, parameters[COL_SCALE], parameters[COL_OFFSET], c_nc_type);
				MGD77_nc_status (nc_put_var1_double (In.nc_id, cdf_var_id, start, &single_val));
			}
		}
		else {	/* Must store array */
			if (LEN)
				MGD77_nc_status (nc_put_vara_schar (In.nc_id, cdf_var_id, start, count, (signed char *)text));
			else if (transform) {
				xtmp = (double *) GMT_memory (VNULL, count[0], sizeof (double), GMT_program);
				n_bad = MGD77_do_scale_offset_before_write (xtmp, colvalue, D->H.n_records, parameters[COL_SCALE], parameters[COL_OFFSET], c_nc_type);
				MGD77_nc_status (nc_put_vara_double (In.nc_id, cdf_var_id, start, count, xtmp));
				GMT_free ((void *)xtmp);
			}
			else 
				MGD77_nc_status (nc_put_vara_double (In.nc_id, cdf_var_id, start, count, colvalue));
		}
		if (n_bad) {	/* Report what we found */
			if (In.verbose_level | 1) fprintf (fp_err, "%s: %s [%s] had %ld values outside valid range <%g,%g> for the chosen type (set to NaN = %g)\n",
				GMT_program, In.NGDC_id, c_abbrev, n_bad, MGD77_Low_val[c_nc_type], MGD77_High_val[c_nc_type], MGD77_NaN_val[c_nc_type]);
		}
		
		MGD77_Close_File (&In);
		MGD77_Free (D);
		n_changed++;
	}
	
	if (got_g || got_i) GMT_free ((void *)f);
	if (got_table) GMT_free ((void *)colvalue);
	if (two_cols) GMT_free ((void *)coldnt);
	
	if (gmtdefs.verbose) {
		if (delete)
			fprintf (stderr, "%s: Removed %ld data columns from %ld MGD77 files\n", GMT_program, n_delete, n_changed);
		else if (got_e)
			fprintf (stderr, "%s: E77 corrections applied to %ld MGD77 files\n", GMT_program, n_changed);
		else
			fprintf (stderr, "%s: Sampled data for %ld MGD77 files\n", GMT_program, n_changed);
	}
	
	MGD77_Path_Free ((int)n_paths, list);
	MGD77_end (&In);

	exit (EXIT_SUCCESS);
}

/* Help functions to decode the -A and -I options */

GMT_LONG decode_A_options (GMT_LONG mode, char *line, char *file, double parameters[])
{
	GMT_LONG error = 0, n;
	
	if (mode == 1) {	/* For *.img files since we need to know data scale and grid mode */
		/* -A[+]i<filename>,<scale>/<mode>[/<lat>] */
		n = sscanf (line, "%[^,],%lf,%lf,%lf", file, &parameters[IMG_SCALE], &parameters[IMG_MODE], &parameters[IMG_LAT]);
		if (n < 3) error = 1;
	}
	else {	/* GMT grid or table: No data scale and mode to worry about */
		/* -A[+]a|c|d|D|e|g|n|t|T<filename> */
		strcpy (file, line);
	}
	
	return (error);
}

GMT_LONG decode_I_options (char *line, char *abbrev, char *name, char *units, char *size, char *comment, double parameters[])
{	/* -I<abbrev>/<name>/<units>/<size>/<scale>/<offset>/\"comment\" */
	GMT_LONG i = 0, k, error;
	GMT_LONG pos = 0;
	char p[BUFSIZ];
	
	while (i < 7 && GMT_strtok (line, "/", &pos, p)) {	/* Process the 7 items */
		switch (i) {
			case 0:
				strcpy (abbrev, p);
				/* Check abbrev for COARDS compliance as well as being lower case */
				for (k = error = 0; abbrev[k]; k++) {
					if (isupper ((int)abbrev[k])) error++;
					if (isalpha ((int)abbrev[k])) continue;
					if (isdigit ((int)abbrev[k]) && k > 0) continue;
					if (abbrev[k] == '_' && k > 0) continue;
					error++;
				}
				if (error) {
					fprintf (stderr, "%s: Abbreviation name should only contain lower case letters, digits, and underscores\n", GMT_program);
					return (TRUE);
				}
				break;
			case 1:
				strcpy (name, p);
				break;
			case 2:
				strcpy (units, p);
				break;
			case 3:
				*size = p[0];
				break;
			case 4:
				parameters[COL_SCALE]  = atof (p);
				break;
			case 5:
				parameters[COL_OFFSET]  = atof (p);
				break;
			case 6:
				strcpy (comment, p);
				break;
		}
		i++;
	}
	
	switch (*size) {	/* Given size, set the NC type */
		case 'b':
			parameters[COL_TYPE] = NC_BYTE;
			break;
		case 'd':
			parameters[COL_TYPE] = NC_DOUBLE;
			break;
		case 'f':
			parameters[COL_TYPE] = NC_FLOAT;
			break;
		case 'i':
			parameters[COL_TYPE] = NC_INT;
			break;
		case 's':
			parameters[COL_TYPE] = NC_SHORT;
			break;
		case 't':
			parameters[COL_TYPE] = NC_CHAR;
			break;
		default:
			fprintf (stderr, "%s: Unknown data type flag %c\n", GMT_program, *size);
			parameters[COL_TYPE] = MGD77_NOT_SET;
			break;
	}
	return ((GMT_LONG)(irint (parameters[COL_TYPE]) == MGD77_NOT_SET) || (i != 7));
}

GMT_LONG skip_if_missing (char *name, char *file, struct MGD77_CONTROL *F, struct MGD77_DATASET *D)
{	/* Used when a needed column is not present and we must free memory and goto next file */
	GMT_LONG id;

	if ((id = MGD77_Get_Column (name, F)) == MGD77_NOT_SET) {
		fprintf (stderr, "%s: Cruise %s is missing column %s which is required for selected operation - skipping\n", GMT_program, file, name);
		MGD77_Free (D);	/* Free memory already allocated by MGD77_Read_File for this aborted effort */
	}
	return (id);
}

GMT_LONG got_default_answer (char *line, char *answer)
{
	GMT_LONG i, k, len;
	
	len = strlen(line) - 1;
	memset ((void *)answer, 0, (size_t)BUFSIZ);	/* No default answer */
	if (line[len] == ']') {	/* Got a default answer for this item */
		for (k = i = len; i && line[i] != '['; i--);
		strncpy (answer, &line[i+1], (size_t)(k - i - 1));
	}
	return (answer[0] != '\0');
}
