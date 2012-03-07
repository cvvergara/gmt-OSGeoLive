/*-----------------------------------------------------------------
 *	$Id: x2sys_datalist.c,v 1.65 2011/07/11 19:22:07 guru Exp $
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
 *      Contact info: www.soest.hawaii.edu/pwessel
 *--------------------------------------------------------------------*/
/* x2sys_datalist will read one or several data files and dump their
 * contents to stdout in ascii or binary (double precision) mode.
 * Input data file formats are determined by the definition file
 * given by the -D option.
 *
 * Author:	Paul Wessel
 * Date:	15-JUN-2004
 * Version:	1.1, based on the spirit of the old xsystem code
 *
 */

#include "x2sys.h"

struct X2SYS_ADJUST {
	int n;
	double *d, *c;
};

GMT_LONG x2sys_load_adjustments (char *DIR, char *TAG, char *track, char *column, struct X2SYS_ADJUST **A);

int main (int argc, char **argv)
{
	char *sfile, *def = "x2sys", *fflags = CNULL, *TAG = CNULL, *correction_table = CNULL;
	char **trk_name = NULL, buffer[BUFSIZ];

	GMT_LONG i, j, k, bad, trk_no, n_data_col_out = 0;
	int n_tracks;

	GMT_LONG error = FALSE, suppress = FALSE, apply_corrections = FALSE, t_given, set_r = FALSE;
	GMT_LONG adjust = FALSE, cmdline_files, special_formatting = FALSE, xy_check, *adj_col = NULL;

	double **data = NULL, west, east, south, north, *out = NULL, correction = 0.0, aux_dvalue[N_GENERIC_AUX];
	double ds = 0.0, cumulative_dist, dist_scale = 1.0, dt, vel_scale = 1.0, adj_amount;

	struct X2SYS_INFO *s = NULL;
	struct X2SYS_FILE_INFO p;		/* File information */
	struct X2SYS_BIX B;
	struct MGD77_CORRTABLE **CORR = NULL;
	struct MGD77_AUX_INFO aux[N_MGD77_AUX];
	struct MGD77_AUXLIST auxlist[N_GENERIC_AUX] = {
		{ "dist",    MGD77_AUX_DS, 0, 0, "d(km)"},
		{ "azim",    MGD77_AUX_AZ, 0, 0, "azimuth"},
		{ "vel",     MGD77_AUX_SP, 0, 0, "v(m/s)"}
	};
	struct X2SYS_ADJUST **A = NULL;
	PFD GMT_azimuth_func;

	argc = (int)GMT_begin (argc, argv);
	for (i = strlen(argv[0]); i >= 0 && argv[0][i] != '/'; i--);
	X2SYS_program = &argv[0][i+1];	/* Name without full path */
	sfile = def;
	west = east = south = north = 0.0;
	t_given = FALSE;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'H':
				case 'M':
				case 'R':
				case 'V':
				case 'b':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */

				case 'A':
					adjust = TRUE;
					break;
				case 'F':
					fflags = &argv[i][2];
					break;
				case 'L':	/* Crossover correction table */
					correction_table = &argv[i][2];
					apply_corrections = TRUE;
					break;
				case 'S':
					suppress = TRUE;
					break;
				case 'T':
					TAG = &argv[i][2];
					t_given = TRUE;
					break;
				default:
					error = TRUE;
					break;
			}
		}
	}

	if (argc == 1 || error) {
		fprintf (stderr, "x2sys_datalist %s - listing of datafiles to stdout\n\n", X2SYS_VERSION);
		fprintf (stderr, "usage: x2sys_datalist <files> -T<TAG> [-A] [-F<fields>] [-L[<corrtable.txt>]]\n");
		fprintf (stderr, "\t[%s] [-S] [-V] [%s] [%s]\n\n", GMT_Rgeo_OPT, GMT_bo_OPT, GMT_mo_OPT);
		fprintf (stderr, "\t<files> is one or more datafiles, or give =<files.lis> for a file with a list of datafiles\n");
		fprintf (stderr, "\t-T <TAG> is the system tag for the data set.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A Use any adjustment splines per track to redistribute COEs between tracks\n");
		fprintf (stderr, "\t   according to their relative weight [no adjustments].\n");
		fprintf (stderr, "\t-F is comma-separated list of column names to output [Default are all fields]\n");
		fprintf (stderr, "\t-L Subtract systematic corrections from the data. If no correction file is given,\n");
		fprintf (stderr, "\t   the default file <TAG>_corrections.txt in $X2SYS_HOME/<TAG> is assumed.\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t-S Suppress output records where all data columns are NaN [Output all records]\n");
		GMT_explain_option ('V');
		GMT_explain_option ('o');
		fprintf (stderr, "\t-m will write a multi-segment header between the output from each file\n");
		exit (EXIT_FAILURE);
	}

	if (!t_given) {
		fprintf (stderr, "%s: ERROR: Must specify -T\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if ((n_tracks = x2sys_get_tracknames (argc, argv, &trk_name, &cmdline_files)) == 0) {
		fprintf (stderr, "%s: No datafiles given!\n", GMT_program);
		exit (EXIT_FAILURE);		
	}

	if (!project_info.region_supplied) set_r = TRUE;

	x2sys_err_fail (x2sys_set_system (TAG, &s, &B, &GMT_io), TAG);

	if (fflags) x2sys_err_fail (x2sys_pick_fields (fflags, s), "-F");

	s->ascii_out = !GMT_io.binary[1];
	xy_check = (s->use_column[s->x_col[GMT_IN]] && s->use_column[s->y_col[GMT_IN]]);	/* TRUE means we requested lon,lat so we can check these for inside/outside */

	if (set_r) {
		west = B.x_min;		east = B.x_max;
		south = B.y_min;	north = B.y_max;
		if (!xy_check) {
			fprintf (stderr, "%s: Region specified but lon/lat (or x/y) columns not selected; -R ignored\n", GMT_program);
		}
	}
	if (project_info.region_supplied) {
		/* Supply dummy linear proj */
		project_info.projection = project_info.xyz_projection[0] = project_info.xyz_projection[1] = GMT_LINEAR;
		project_info.pars[0] = project_info.pars[1] = 1.0;
		if (s->geographic) project_info.degree[0] = project_info.degree[1] = TRUE;
		if (west < 0.0 && east < 0.0) {
			west += 360.0;
			east += 360.0;
		}
		GMT_err_fail (GMT_map_setup (west, east, south, north), "");
	}

	out = (double *) GMT_memory (VNULL, (size_t)s->n_fields, sizeof (double), "x2sys_datalist");

	for (i = 0; i < s->n_out_columns; i++) {	/* Set output formats */
		if (i == s->t_col[GMT_OUT])
			GMT_io.out_col_type[i] = GMT_IS_ABSTIME;
		else if (i == s->x_col[GMT_OUT])
			GMT_io.out_col_type[i] = (!strcmp (s->info[s->out_order[i]].name, "lon")) ? GMT_IS_LON : GMT_IS_FLOAT;
		else if (i == s->y_col[GMT_OUT])
			GMT_io.out_col_type[i] = (!strcmp (s->info[s->out_order[i]].name, "lat")) ? GMT_IS_LAT : GMT_IS_FLOAT;
		else
			GMT_io.out_col_type[i] = GMT_IS_FLOAT;

		if (s->info[s->out_order[i]].format[0] != '-') special_formatting = TRUE;
	}
	if (GMT_io.binary[GMT_OUT]) special_formatting = FALSE;

	if (suppress) {	/* Must count output data columns (except t, x, y) */
		for (i = n_data_col_out = 0; i < s->n_out_columns; i++) {
			if (i == s->t_col[GMT_OUT]) continue;
			if (i == s->x_col[GMT_OUT]) continue;
			if (i == s->y_col[GMT_OUT]) continue;
			n_data_col_out++;
		}
	}

	MGD77_Set_Unit (s->unit[X2SYS_DIST_SELECTION], &dist_scale, -1);	/* Gets scale which multiplies meters to chosen distance unit */
	MGD77_Set_Unit (s->unit[X2SYS_SPEED_SELECTION], &vel_scale,  -1);	/* Sets output scale for distances using in velocities */
	switch (s->unit[X2SYS_SPEED_SELECTION][0]) {
		case 'c':
			vel_scale = 1.0;
			break;
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
	switch (s->unit[X2SYS_DIST_SELECTION][0]) {
		case 'c':
			strcpy (auxlist[MGD77_AUX_SP].header, "d(user)");
			break;
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

	switch (s->dist_flag) {
		case 0:	/* Cartesian */
			GMT_distance_func = GMT_cartesian_dist;
			GMT_azimuth_func  = GMT_az_backaz_cartesian;
			break;
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
	
	if (apply_corrections) {	/* Load an ephemeral correction table */
		x2sys_get_corrtable (s, correction_table, n_tracks, trk_name, NULL, aux, auxlist, &CORR);
		if (auxlist[MGD77_AUX_SP].requested && s->t_col[GMT_IN] == -1) {
			fprintf (stderr, "%s: Selected correction table requires velocity which implies time (not selected)\n", GMT_program);
			MGD77_Free_Correction (CORR, n_tracks);
			x2sys_free_list (trk_name, n_tracks);
			exit (EXIT_FAILURE);
		}
	}

	if (adjust) {
		A = (struct X2SYS_ADJUST **) GMT_memory (VNULL, s->n_out_columns, sizeof (struct X2SYS_ADJUST *), GMT_program);
		adj_col = (GMT_LONG *) GMT_memory (VNULL, s->n_out_columns, sizeof (GMT_LONG), GMT_program);
	}
	
	for (trk_no = 0; trk_no < n_tracks; trk_no++) {

		if (gmtdefs.verbose) fprintf (stderr, "x2sys_datalist: Reading track %s\n", trk_name[trk_no]);

		x2sys_err_fail ((s->read_file) (trk_name[trk_no], &data, s, &p, &GMT_io, &k), trk_name[trk_no]);

		if (apply_corrections && s->t_col[GMT_IN] >= 0) MGD77_Init_Correction (CORR[trk_no], data);	/* Initialize origins if needed */

		if (adjust) {
			for (k = 0; k < s->n_out_columns; k++) adj_col[k] = x2sys_load_adjustments (X2SYS_HOME, TAG, trk_name[trk_no], s->info[s->out_order[k]].name, &A[k]);
		}

		if (GMT_io.multi_segments[GMT_OUT]) GMT_write_segmentheader (GMT_stdout, s->n_fields);

		cumulative_dist = 0.0;
		for (j = 0; j < p.n_rows; j++) {
			if (project_info.region_supplied && xy_check && GMT_map_outside (data[s->x_col[GMT_OUT]][j], data[s->y_col[GMT_OUT]][j])) continue;
			if (suppress) {
				for (k = bad = 0; k < s->n_out_columns; k++) {
					if (k == s->t_col[GMT_OUT]) continue;
					if (k == s->x_col[GMT_OUT]) continue;
					if (k == s->y_col[GMT_OUT]) continue;
					if (GMT_is_dnan (data[k][j])) bad++;
				}
				if (bad == n_data_col_out) continue;
			}
			if (auxlist[MGD77_AUX_AZ].requested) {
				if (j == 0)	/* Look forward at first point to get an azimuth */
					aux_dvalue[MGD77_AUX_AZ] = GMT_azimuth_func (data[s->x_col[GMT_OUT]][1], data[s->y_col[GMT_OUT]][1], data[s->x_col[GMT_OUT]][0], data[s->y_col[GMT_OUT]][0], FALSE);
				else		/* else go from previous to this point */
					aux_dvalue[MGD77_AUX_AZ] = GMT_azimuth_func (data[s->x_col[GMT_OUT]][j], data[s->y_col[GMT_OUT]][j], data[s->x_col[GMT_OUT]][j-1], data[s->y_col[GMT_OUT]][j-1], FALSE);
			}
			if (auxlist[MGD77_AUX_DS].requested) {
				ds = (j == 0) ? 0.0 : dist_scale * GMT_distance_func (data[s->x_col[GMT_OUT]][j], data[s->y_col[GMT_OUT]][j], data[s->x_col[GMT_OUT]][j-1], data[s->y_col[GMT_OUT]][j-1]);
				cumulative_dist += ds;
				aux_dvalue[MGD77_AUX_DS] = cumulative_dist;
			}
			if (auxlist[MGD77_AUX_SP].requested) {
				dt =  (j == 0) ? data[s->t_col[GMT_OUT]][1] - data[s->t_col[GMT_OUT]][0] : data[s->t_col[GMT_OUT]][j] - data[s->t_col[GMT_OUT]][j-1];
				aux_dvalue[MGD77_AUX_SP] = (GMT_is_dnan (dt) || dt == 0.0) ? GMT_d_NaN : vel_scale * ds / dt;
			}
			for (k = 0; k < s->n_out_columns; k++) {	/* Load output record */
				correction = (apply_corrections) ? MGD77_Correction (CORR[trk_no][k].term, data, aux_dvalue, j) : 0.0;
				if (adjust && adj_col[k]) {
					if (GMT_intpol (A[k]->d, A[k]->c, A[k]->n, 1, &aux_dvalue[MGD77_AUX_DS], &adj_amount, gmtdefs.interpolant)) {
						fprintf (stderr, "%s: Error interpolating adjustment for %s near row %ld - no adjustment made!\n", GMT_program, s->info[s->out_order[k]].name, j);
						adj_amount = 0.0;
					}
					correction -= adj_amount;
				}
				out[k] = data[k][j] - correction;	/* This loads out in the correct output order */
			}
			if (special_formatting)  {	/* use the specified formats */
				for (k = 0; k < s->n_out_columns; k++) {
					if (s->info[s->out_order[k]].format[0] == '-')
						GMT_ascii_output_one (GMT_stdout, out[k], k);
					else {
						if (!GMT_is_dnan (out[k])) {
							sprintf (buffer, s->info[s->out_order[k]].format, out[k]);
							GMT_fputs (buffer, GMT_stdout);
						}
						else
							GMT_fputs ("NaN", GMT_stdout);
					}
					(k == (s->n_out_columns - 1)) ? GMT_fputs ("\n", GMT_stdout) : GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
				}
			}
			else {
				GMT_output (GMT_stdout, s->n_out_columns, out);
			}
		}

		x2sys_free_data (data, s->n_out_columns, &p);
		for (k = 0; k < s->n_fields; k++) if (adjust && adj_col[k]) {
			GMT_free ((void *)A[k]->d);
			GMT_free ((void *)A[k]->c);
			GMT_free ((void *)A[k]);
		}
	}

	if (apply_corrections) MGD77_Free_Correction (CORR, n_tracks);

	x2sys_end (s);
	GMT_free ((void *)out);
	if (adjust) {
		GMT_free ((void *)A);
		GMT_free ((void *)adj_col);
	}
	x2sys_free_list (trk_name, n_tracks);
	
	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

GMT_LONG x2sys_load_adjustments (char *DIR, char *TAG, char *track, char *column, struct X2SYS_ADJUST **A)
{
	GMT_LONG n_fields, n_expected_fields = 2, n = 0, k, n_alloc = GMT_CHUNK;
	GMT_LONG type[2] = {GMT_IS_FLOAT, GMT_IS_FLOAT};
	double *in;
	char file[BUFSIZ];
	FILE *fp;
	struct X2SYS_ADJUST *adj;
	
	sprintf (file, "%s/%s/%s.%s.adj", DIR, TAG, track, column);
	if ((fp = GMT_fopen (file, "r")) == NULL) return FALSE;	/* Nuthin' to read */
	
	adj = (struct X2SYS_ADJUST *)GMT_memory (VNULL, 1, sizeof (struct X2SYS_ADJUST), GMT_program);
	adj->d = (double *)GMT_memory (VNULL, n_alloc, sizeof (double), GMT_program);
	adj->c = (double *)GMT_memory (VNULL, n_alloc, sizeof (double), GMT_program);
	for (k = 0; k < 2; k++) l_swap (type[k], GMT_io.in_col_type[k]);	/* Save original input type setting */
	while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */
		adj->d[n] = in[0];
		adj->c[n] = in[1];
		n++;
		if (n == n_alloc) {
			n_alloc <<= 1;
			adj->d = (double *)GMT_memory ((void *)adj->d, n_alloc, sizeof (double), GMT_program);
			adj->c = (double *)GMT_memory ((void *)adj->c, n_alloc, sizeof (double), GMT_program);
		}
	}
	GMT_fclose (fp);
	adj->d = (double *)GMT_memory ((void *)adj->d, n, sizeof (double), GMT_program);
	adj->c = (double *)GMT_memory ((void *)adj->c, n, sizeof (double), GMT_program);
	adj->n = (int)n;
	*A = adj;
	for (k = 0; k < 2; k++) l_swap (GMT_io.in_col_type[k], type[k]);	/* Restore original input type setting */
	return (TRUE);
}
