/*
 *	$Id: gmtstitch.c 9923 2012-12-18 20:45:53Z pwessel $
 *	Copyright (c) 1991-2013 by P. Wessel and W. H. F. Smith
 */
/* gmtstitch will combine pieces of coastlines or similar segments into a
 * continuous line, polygon, or group of lines/polygons so that the jump
 * between segment endpoints exceeds a specified threshold.
 *
 * Paul Wessel, March, 2006, derived from our earlier GSHHG processing tools July, 1994.
 */

#include "gmt.h"

#define SEG_I	0
#define SEG_J	1
#define END_A	0
#define END_B	1

#define CLOSED	0
#define OPEN	1

struct BUDDY {
	GMT_LONG id;
	GMT_LONG orig_id;
	GMT_LONG end_order;
	double dist, next_dist;
};

struct LINK {
	GMT_LONG id;
	GMT_LONG orig_id;
	GMT_LONG group;
	GMT_LONG pos;
	GMT_LONG n;
	GMT_LONG used;
	double x_end[2];
	double y_end[2];
	struct BUDDY buddy[2];
};


int main (int argc, char **argv)
{
	struct LINK *seg = NULL;
	GMT_LONG nearest_end[2][2], ii, end, n_open;
	GMT_LONG i, j, k, np, distance_flag = 0, n_files = 0, ns, n_args, id, pos, start_id, done, end_order;
	GMT_LONG n_new, n, chain = 0, n_islands = 0, n_trouble = 0, n_closed = 0, id2, L, G, error = 0, mode = 0;
	GMT_LONG n_alloc = GMT_SMALL_CHUNK, n_id_alloc = GMT_CHUNK, out_seg, match = 0, n_steps, n_points = 0;	
	GMT_LONG individual_file = FALSE, save_type = FALSE, write_links = FALSE, do_lists = FALSE, first, *skip = NULL;
	GMT_LONG nn_check = FALSE, separate = FALSE;
	struct GMT_DATASET *D = NULL;
	double dd[2][2], p_dummy_x, p_dummy_y, p_last_x, p_last_y, p_first_x, p_first_y, distance;
	double cutoff = 0.0, nn_dist = 0.0, scl = 1.0, closed_dist = 0.0;
	char format[BUFSIZ], filename[BUFSIZ], linkfile[BUFSIZ], lists[BUFSIZ], A[GMT_TEXT_LEN], B[GMT_TEXT_LEN];
	char *BE = "BE", closedfile[BUFSIZ];
	FILE *fp = NULL, *fp3 = NULL, *fpl[2] = {NULL, NULL}, *fp_closed = NULL;

	void Write_This_Segment (FILE *fp, struct GMT_LINE_SEGMENT *line, GMT_LONG start, GMT_LONG end);
	GMT_LONG connect (struct LINK *S, int id, int order, double cutoff, GMT_LONG nn_check, double nn_dist);
	
	argc = (int)GMT_begin (argc, argv);
	memset ((void *)format, 0, BUFSIZ);
	memset ((void *)linkfile, 0, BUFSIZ);
	memset ((void *)lists, 0, BUFSIZ);
	memset ((void *)closedfile, 0, BUFSIZ);
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'H':
				case 'M':
				case 'V':
				case ':':
				case 'b':
				case 'f':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], NULL, NULL, NULL, NULL);
					break;

				/* Supplemental parameters */

				case 'C':               /* Separate closed from open segments  */
					separate = TRUE;
					if (argv[i][2]) strcpy (closedfile, &argv[i][2]);
					break;
				case 'D':               /* Write each segment to a separate output file */
					if (argv[i][2]) strcpy (format, &argv[i][2]);
					individual_file = TRUE;
					break;
				case 'L':               /* Write link information to file */
					if (argv[i][2]) strcpy (linkfile, &argv[i][2]);
					write_links = TRUE;
					break;
				case 'Q':               /* Write names of individual files to list(s) */
					if (argv[i][2]) strcpy (lists, &argv[i][2]);
					do_lists = TRUE;
					break;
				case 'T':
					n = sscanf (&argv[i][2], "%[^/]/%s", A, B);
					cutoff = GMT_getradius (A);
					if (A[strlen(A)-1] == 'k') distance_flag = 1;
					if (A[strlen(A)-1] == 'e') {distance_flag = 1, scl = 0.001;}	/* Convert to km */
					if (A[strlen(A)-1] == 'K') distance_flag = 2;
					if (A[strlen(A)-1] == 'E') {distance_flag = 2, scl = 0.001;}	/* Convert to km */
					cutoff *= scl;
					if (n == 2) {
						nn_dist = scl * GMT_getradius (B);
						nn_check = TRUE;
					}
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			n_files++;
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "gmtstitch %s - Join individual lines whose end points match within tolerance\n\n", GMT_VERSION);
		fprintf (stderr, "usage: gmtstitch [<infiles>] [%s] [-C<closedfile>] [-D[<template>]] [-L[<linkfile>]] [-Q<list>]\n", GMT_H_OPT);
		fprintf (stderr, "\t-T<cutoff>[m|c|e|E|k|K][/<nn_dist>] [-V[l]] [%s] [%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_b_OPT, GMT_f_OPT, GMT_m_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\tinfiles (in ASCII or binary) have 2 or more columns with (x,y) or (y,x) in first columns.\n");
		fprintf (stderr, "\t  If no file(s) is given, standard input is read.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-C will write the already-closed polygons to <closedfile> [gmtstitch_closed.d]\n");
		fprintf (stderr, "\t   and write all other segments unchanged to stdout in multisegment format.\n");
		fprintf (stderr, "\t-D writes individual segments to separate files [Default writes one multisegment file to stdout].\n");
		fprintf (stderr, "\t   Append file name template which MUST contain a C-format specified for an integer (e.g., %%d).\n");
		fprintf (stderr, "\t   If the format also includes a %%c string BEFORE the %%d part we replace it with C(losed) or O(pen).\n");
		fprintf (stderr, "\t   [Default uses gmtstitch_segment_%%d.d]\n");
		fprintf (stderr, "\t-L write link information (seg id, begin/end nearest seg id, end, and distance) to file [link.d].\n");
		fprintf (stderr, "\t   Link output excludes duplicates and segments already forming a closed polygon.\n");
		GMT_explain_option ('H');
		GMT_explain_option ('V');
		fprintf (stderr, "\t-T sets cutoff distance in data units; append m or c for minutes or seconds [0].\n");
		fprintf (stderr, "\t   Append e for meters or k for km (implies -fg), use flat Earth approximation.\n");
		fprintf (stderr, "\t   Append E for meters or K for km (implies -fg), use exact geodesic distances.\n");
		fprintf (stderr, "\t   If the current ELLIPSOID is Sphere then spherical great circle distances are used.\n");
		fprintf (stderr, "\t   If two lines has endpoints that are closer than this cutoff they will be joined.\n");
		fprintf (stderr, "\t   Optionally, append <nn_dist> which adds the requirement that the second closest\n");
		fprintf (stderr, "\t   match must exceed <nn_dist> (assumed to be in the same units as <cutoff>)\n");
		fprintf (stderr, "\t-Q Used with -D to write names of files to a list.  Optionally give list name [gmtstitch_list.d].\n");
		fprintf (stderr, "\t   Embed %%c in the list name to write two separate lists: one for C(losed) and one for O(pen).\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2 input columns\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (n_files == 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  No input files specified\n", GMT_program);
		error++;
	}
	if (cutoff < 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  -T cutoff must be >= 0!\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 2) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least %d columns\n", GMT_program, 2);
		error++;
	}
	if (separate && individual_file) {
               fprintf (stderr, "%s: GMT SYNTAX ERROR.  Option -C cannot be used with -D!\n", GMT_program);
		error++;
	}
	if (do_lists && !individual_file) {
               fprintf (stderr, "%s: GMT SYNTAX ERROR.  Option -Q requires -D!\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	/* Now we are ready to take on some input values */

	if (individual_file) {
		if (!format[0]) strcpy (format, "gmtstitch_segment_%ld.d");
		if (!lists[0]) strcpy (lists, "gmtstitch_list.d");
		if (strstr (format, "%c")) save_type = TRUE;
		if (do_lists) {
			if (strstr (lists, "%c")) {
				sprintf (filename, lists, 'C');
				if ((fpl[CLOSED] = fopen (filename, GMT_io.w_mode)) == NULL ) {
					fprintf (stderr, "%s: Error creating file %s\n", GMT_program, filename);
					exit (EXIT_FAILURE);
				}
				sprintf (filename, lists, 'O');
				if ((fpl[OPEN] = fopen (filename, GMT_io.w_mode)) == NULL ) {
					fprintf (stderr, "%s: Error creating file %s\n", GMT_program, filename);
					exit (EXIT_FAILURE);
				}
			}
			else {
				if ((fpl[CLOSED] = fopen (lists, GMT_io.w_mode)) == NULL ) {
					fprintf (stderr, "%s: Error creating file %s\n", GMT_program, lists);
					exit (EXIT_FAILURE);
				}
				fpl[OPEN] = fpl[CLOSED];	/* Same file */
			}
			
		}
	}
	
	n_args = (argc > 1) ? argc : 2;

	switch (distance_flag) {	/* Take different action depending on how we want distances calculated */
		case 0:		/* Cartesian distance (or spherical degrees if -fg was used) */
			GMT_distance_func = (GMT_io.in_col_type[0] == GMT_IS_LON) ? GMT_great_circle_dist_km : GMT_cartesian_dist;
			break;
		case 1:		/* Flat Earth Approximation */
			GMT_distance_func = GMT_flatearth_dist_km;
			break;
		case 2:		/* Full spherical calculation */
			GMT_distance_func = GMT_great_circle_dist_km;
			break;
		case 3:		/* Full Ellipsoidal calculation */
			GMT_distance_func = GMT_geodesic_dist_km;
		break;
	}

	/* Allocate memory and read in all the files; each file can have many lines (-M) */
	
	D = (struct GMT_DATASET *) GMT_memory (VNULL, (size_t)1, sizeof (struct GMT_DATASET), GMT_program);
	if (n_files == 0) {	/* Read stdin */
		D->table = (struct GMT_TABLE **) GMT_memory (VNULL, (size_t)1, sizeof (struct GMT_TABLE *), GMT_program);
		D->n_tables = 1;
		if (gmtdefs.verbose) fprintf (stderr, "%s: Reading standard input\n", GMT_program);
		GMT_import_table ((void *)GMT_stdin, GMT_IS_STREAM, &D->table[0], 0.0, FALSE, FALSE, TRUE);
		n_points = D->table[0]->n_records;
	}
	else {
		D->table = (struct GMT_TABLE **) GMT_memory (VNULL, (size_t)n_alloc, sizeof (struct GMT_TABLE *), GMT_program);
		for (k = 1, i = 0; k < argc; k++) {
			if (argv[k][0] == '-') continue;

			if (gmtdefs.verbose) fprintf (stderr, "%s: Reading file %s\n", GMT_program, argv[k]);

			if (GMT_import_table ((void *)argv[k], GMT_IS_FILE, &D->table[i], 0.0, FALSE, FALSE, TRUE) == GMT_IO_EOF) continue;
			n_points += D->table[i]->n_records;
			i++;
			if (i == n_alloc) {
				n_alloc <<= 1;
				D->table = (struct GMT_TABLE **) GMT_memory ((void *)D->table, (size_t)n_alloc, sizeof(struct GMT_TABLE *), GMT_program);
			}
		}
		if (i < n_alloc) D->table = (struct GMT_TABLE **) GMT_memory ((void *)D->table, (size_t)i, sizeof(struct GMT_TABLE *), GMT_program);
		D->n_tables = i;
	}
	if (n_points == 0) {	/* Empty files, nothing to do */
		GMT_free_dataset (D);
		GMT_exit (EXIT_SUCCESS);
	}
	
	seg = (struct LINK *) GMT_memory (VNULL, (size_t)n_id_alloc, sizeof(struct LINK), GMT_program);
	id = pos = ns = out_seg = 0;
	if (gmtdefs.verbose) fprintf (stderr, "%s: Check for closed polygons\n", GMT_program);
	
	/* Closed polygons are already finished - just identify, write out, and move on */
	
	if (separate) {
		if (!closedfile[0]) strcpy (closedfile, "gmtstitch_closed.d");
		if ((fp_closed = GMT_fopen (closedfile, GMT_io.w_mode)) == NULL ) {
			fprintf (stderr, "%s: Error creating file %s\n", GMT_program, closedfile);
			exit (EXIT_FAILURE);
		}
	}
	
	ns = id = -1;
	n_open = n_closed = 0;
	closed_dist = (separate) ? cutoff : 0.0;
	for (k = 0; k < D->n_tables; k++) {
		for (j = 0; j < D->table[k]->n_segments; j++) {
			np = D->table[k]->segment[j]->n_rows;
			ns++;
			distance = (GMT_distance_func) (D->table[k]->segment[j]->coord[GMT_X][0], D->table[k]->segment[j]->coord[GMT_Y][0], D->table[k]->segment[j]->coord[GMT_X][np-1], D->table[k]->segment[j]->coord[GMT_Y][np-1]);
			if (np > 2 && distance <= closed_dist) {	/* Already closed, just write out and forget in the rest of the program */
				if (individual_file) {
					(save_type) ? sprintf (filename, format, 'C', out_seg) : sprintf (filename, format, out_seg);
					if ((fp = GMT_fopen (filename, GMT_io.w_mode)) == NULL ) {
						fprintf (stderr, "%s: Error creating file %s\n", GMT_program, filename);
						exit (EXIT_FAILURE);
					}
					if (do_lists) fprintf (fpl[CLOSED], "%s\n", filename);
				}
				else if (separate)
					fp = fp_closed;
				else
					fp = GMT_stdout;
				if (D->table[k]->segment[j]->header) strcpy (GMT_io.segment_header, D->table[k]->segment[j]->header);
				GMT_write_segmentheader (fp, D->table[k]->segment[j]->n_columns);
				Write_This_Segment (fp, D->table[k]->segment[j], (GMT_LONG)0, np-1);
				if (separate && distance > 0.0) Write_This_Segment (fp, D->table[k]->segment[j], (GMT_LONG)0, 0);	/* Close polygon */
				if (individual_file) GMT_fclose (fp);
				n_islands++;
				out_seg++;
				n_closed++;
			}
			else if (separate) {	/* Write open segment to stdout */
				if (D->table[k]->segment[j]->header) strcpy (GMT_io.segment_header, D->table[k]->segment[j]->header);
				GMT_write_segmentheader (GMT_stdout, D->table[k]->segment[j]->n_columns);
				Write_This_Segment (GMT_stdout, D->table[k]->segment[j], (GMT_LONG)0, np-1);
				n_open++;
			}
			else { /* Here we have a segment that is not closed.  Store refs to D->table and copy end points */
				if (id == -1) id = 0;
				if (np == 1 && gmtdefs.verbose) fprintf (stderr, "%s: Segment %ld only has a single point.  Stitching may require additional stitching.\n", GMT_program, id);
				seg[id].id = id;
				seg[id].orig_id = ns;
				seg[id].group = k;
				seg[id].pos = j;
				seg[id].n = np;
				seg[id].x_end[0] = D->table[k]->segment[j]->coord[GMT_X][0];
				seg[id].y_end[0] = D->table[k]->segment[j]->coord[GMT_Y][0];
				seg[id].x_end[1] = D->table[k]->segment[j]->coord[GMT_X][np-1];
				seg[id].y_end[1] = D->table[k]->segment[j]->coord[GMT_Y][np-1];
				seg[id].buddy[0].dist = seg[id].buddy[1].dist = seg[id].buddy[0].next_dist = seg[id].buddy[1].next_dist = DBL_MAX;
				id++;
				if (id == n_id_alloc) {
					n_id_alloc <<= 1;
					seg = (struct LINK *) GMT_memory ((void *)seg, (size_t)n_id_alloc, sizeof(struct LINK), GMT_program);
				}
			}
		}
	}
	if (separate) {	/* Report and exit */
		if (gmtdefs.verbose) fprintf (stderr, "%s: Separated %ld closed and %ld open segments\n", GMT_program, n_closed, n_open);
		GMT_free_dataset (D);
		exit (EXIT_SUCCESS);
	}
	if (id == -1) {	/* All are closed or we wrote them out as separate */
		if (gmtdefs.verbose) fprintf (stderr, "%s: All %ld segments already form closed polygons\n", GMT_program, n_closed);
		GMT_free_dataset (D);
		exit (EXIT_SUCCESS);
	}
	ns = id;
	
	if (ns < n_id_alloc) seg = (struct LINK *) GMT_memory ((void *)seg, (size_t)ns, sizeof(struct LINK), GMT_program);
	skip = (GMT_LONG *) GMT_memory (VNULL, (size_t)ns, sizeof(GMT_LONG), GMT_program);
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Found %ld closed polygons\n", GMT_program, n_islands);
	
	/* The algorithm will be confused if there are identical duplicates of segments - thus we check */
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Check for duplicate lines\n", GMT_program);
	for (i = 0; i < ns; i++) {
		if (skip[i]) continue;	/* Skip segment that has been determined to be a duplicate segment */
		for (j = i + 1; j < ns; j++) {
			if (skip[j]) continue;	/* Skip segment that has been determined to be a duplicate segment */
			if ((seg[i].x_end[0] == seg[j].x_end[0] && seg[i].y_end[0] == seg[j].y_end[0]) ||
			    (seg[i].x_end[0] == seg[j].x_end[1] && seg[i].y_end[0] == seg[j].y_end[1]) ||
			    (seg[i].x_end[1] == seg[j].x_end[0] && seg[i].y_end[1] == seg[j].y_end[0]) ||
			    (seg[i].x_end[1] == seg[j].x_end[1] && seg[i].y_end[1] == seg[j].y_end[1])) {
			    	if (seg[i].n == seg[j].n) {
					for (k = match = 0; k < seg[i].n && k == match; k++) {
						match += (D->table[seg[i].group]->segment[seg[i].pos]->coord[GMT_X][k] == D->table[seg[j].group]->segment[seg[j].pos]->coord[GMT_X][k] && 
						          D->table[seg[i].group]->segment[seg[i].pos]->coord[GMT_Y][k] == D->table[seg[j].group]->segment[seg[j].pos]->coord[GMT_Y][k]);
					}
					match = (match == seg[i].n) ? 1 : 0;
					if (match) {
						if (gmtdefs.verbose) fprintf (stderr, "%s: Segments %ld and %ld are duplicates - Segment %ld will be ignored\n", GMT_program, i, j, j);
						skip[j] = TRUE;
					}
				}
			}
		}
	}
	
	/* Eliminate the duplicate segments from consideration */
	
	for (i = j = 0; i < ns; i++) {
		if (skip[i]) continue;
		if (i > j) seg[j] = seg[i];
		j++;
	}
	if (j < ns && gmtdefs.verbose) fprintf (stderr, "%s: %ld duplicate segment removed\n", GMT_program, ns - j);
	ns = j;
	GMT_free ((void *)skip);
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Calculate and rank end point separations [cutoff = %g nn_dist = %g]\n", GMT_program, cutoff, nn_dist);
	
	/* We determine the distance from each segments two endpoints to the two endpoints on every other
	 * segment; this is four distances per segment.  We then assign the nearest endpoint to each end
	 * of a segment to the buddy structure which keeps the id of the nearest segment so far.
	 */
	 
	for (i = 0; i < ns; i++) {
		for (j = i; j < ns; j++) {
			/* nearest_end indicates which end is closest to this end */
			if (i == j) {	/* Store offset between the endpoints of a single segment (should be 0 if closed) */
				dd[SEG_I][END_A] = dd[SEG_J][END_B] = DBL_MAX;
				dd[SEG_I][END_B] = dd[SEG_J][END_A] = (GMT_distance_func) (seg[i].x_end[END_A], seg[i].y_end[END_A], seg[i].x_end[END_B], seg[i].y_end[END_B]);
    				nearest_end[SEG_I][END_A] = nearest_end[SEG_J][END_A] = END_B;
    				nearest_end[SEG_J][END_B] = nearest_end[SEG_I][END_B] = END_A;
			}
			else {	/* Store the distances between the 4 possible end-to-end configurations */
				dd[SEG_I][END_A] = (GMT_distance_func) (seg[i].x_end[END_A], seg[i].y_end[END_A], seg[j].x_end[END_A], seg[j].y_end[END_A]);
				dd[SEG_I][END_B] = (GMT_distance_func) (seg[i].x_end[END_A], seg[i].y_end[END_A], seg[j].x_end[END_B], seg[j].y_end[END_B]);
				dd[SEG_J][END_A] = (GMT_distance_func) (seg[i].x_end[END_B], seg[i].y_end[END_B], seg[j].x_end[END_A], seg[j].y_end[END_A]);
				dd[SEG_J][END_B] = (GMT_distance_func) (seg[i].x_end[END_B], seg[i].y_end[END_B], seg[j].x_end[END_B], seg[j].y_end[END_B]);
    				for (end = 0; end < 2; end++) nearest_end[SEG_I][end] = (dd[end][END_A] < dd[end][END_B]) ? END_A : END_B;
    				for (end = 0; end < 2; end++) nearest_end[SEG_J][end] = (dd[END_A][end] < dd[END_B][end]) ? END_A : END_B;
    			}
    			/* Update list of closest matches for both ends */
    			for (ii = 0; ii < 2; ii++) {	/* For each end of the segment */
    				end = nearest_end[SEG_I][ii];	/* The end of segment j that was closest to segment i's end ii */
    				if (dd[ii][end] < seg[i].buddy[ii].dist) {	/* This distance is shorter than the previous shortest distance */
					seg[i].buddy[ii].next_dist = seg[i].buddy[ii].dist;	/* Previous closest distance */
					seg[i].buddy[ii].orig_id = seg[j].orig_id;
					seg[i].buddy[ii].id = j;
					seg[i].buddy[ii].dist = dd[ii][end];
					seg[i].buddy[ii].end_order = end;
    				}
    				end = nearest_end[SEG_J][ii];	/* The end of segment i that was closest to segment j's end ii */
    				if (dd[end][ii] < seg[j].buddy[ii].dist) {	/* This distance is shorter than the previous shortest distance */
 					seg[j].buddy[ii].next_dist = seg[j].buddy[ii].dist;	/* Previous closest distance */
					seg[j].buddy[ii].orig_id = seg[i].orig_id;
 					seg[j].buddy[ii].id = i;
					seg[j].buddy[ii].dist = dd[end][ii];
					seg[j].buddy[ii].end_order = end;
    				}
    			}
		}
	}
	if (write_links) {
		char name[BUFSIZ], name0[BUFSIZ], name1[BUFSIZ], *pp = NULL;
		if (!linkfile[0]) strcpy (linkfile, "links.d");
		if ((fp3 = fopen (linkfile, "w")) == NULL) {
			fprintf (stderr, "%s: Error creating file %s\n", GMT_program, linkfile);
			exit (EXIT_FAILURE);
		}
		fprintf (fp3, "# segid\tbegin_id\tb_pt\tb_dist\tb_nndist\tend_id\te_pt\te_dist\te_nndist\n");
		for (i = 0; i < ns; i++) {
			G = seg[i].group;	L = seg[i].pos;
			if (D->table[G]->segment[L]->header && (pp = strstr (D->table[G]->segment[L]->header, "-L"))) {
				strcpy (name, &pp[2]);
				for (j = 0; name[j]; j++) if (name[j] == ' ') name[j] = '\0';		/* Just truncate after 1st word */
			} else sprintf (name, "%ld", seg[i].orig_id);
			G = seg[seg[i].buddy[0].id].group;	L = seg[seg[i].buddy[0].id].pos;
			if (D->table[G]->segment[L]->header && (pp = strstr (D->table[G]->segment[L]->header, "-L"))) {
				strcpy (name0, &pp[2]);
				for (j = 0; name0[j]; j++) if (name0[j] == ' ') name0[j] = '\0';	/* Just truncate after 1st word */
			} else sprintf (name0, "%ld", seg[i].buddy[0].orig_id);
			G = seg[seg[i].buddy[1].id].group;	L = seg[seg[i].buddy[1].id].pos;
			if (D->table[G]->segment[L]->header && (pp = strstr (D->table[G]->segment[L]->header, "-L"))) {
				strcpy (name1, &pp[2]);
				for (j = 0; name1[j]; j++) if (name1[j] == ' ') name1[j] = '\0';	/* Just truncate after 1st word */
			} else sprintf (name1, "%ld", seg[i].buddy[1].orig_id);
			fprintf (fp3, "%s\t%s\t%c\t%g\t%g\t%s\t%c\t%g\t%g\n", name, name0, BE[seg[i].buddy[0].end_order], seg[i].buddy[0].dist, seg[i].buddy[0].next_dist, name1, \
				BE[seg[i].buddy[1].end_order], seg[i].buddy[1].dist, seg[i].buddy[1].next_dist);
		}
		fclose (fp3);
	}
	start_id = done = n_closed = 0;
	p_dummy_x = p_dummy_y = DBL_MAX;
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Assemble new segments\n", GMT_program);
	while (!done) {
	
		/* Find the 'beginning' of the chain that this segment belongs to by tracing the connections
		 * until we either reappear at the starting point (a closed loop) or we reach an end (i.e.,
		 * the nearest next endpoint is beyond the separation threshold. */
		
		done = FALSE;
		id = start_id;
		end_order = n_steps = 0;
#if DEBUG2
		if (gmtdefs.verbose) fprintf (stderr, "%ld\n", seg[id].orig_id);
#endif
		while (!done && connect (seg, (int)id, (int)end_order, cutoff, nn_check, nn_dist)) {
			id2 = seg[id].buddy[end_order].id;
#if DEBUG2
			if (gmtdefs.verbose) fprintf (stderr, "%ld\n", seg[id2].orig_id);
#endif
			if (id2 == start_id)	/* Closed polygon, start here */
				done = TRUE;
			if (id2 == id || n_steps > ns) {	/* Not good... */
				done = TRUE;
				n_trouble++;
			}
			else {	/* Trace the connection to the next segment */
				end_order = !seg[id].buddy[end_order].end_order;
				id = id2;
			}
			n_steps++;
		}
				
		/* This id should be the beginning of a segment.  Now trace forward and dump out the chain */
		
		/* First dump start segment */
		
		start_id = id;
		
		memset (GMT_io.segment_header, 0, (size_t)BUFSIZ);
		if (individual_file) {
			mode = OPEN;
			(save_type) ? sprintf (filename, format, 'O', out_seg) : sprintf (filename, format, out_seg);
			if ((fp = GMT_fopen (filename, GMT_io.w_mode)) == NULL ) {
				fprintf (stderr, "%s: Error creating file %s\n", GMT_program, filename);
				exit (EXIT_FAILURE);
			}
		}
		else
			fp = GMT_stdout;
		sprintf (GMT_io.segment_header, "%c Possibly a composite segment; see comments for individual segment headers\n", GMT_io.EOF_flag[GMT_OUT]);
		GMT_write_segmentheader (fp, D->table[seg[id].group]->segment[seg[id].pos]->n_columns);
		
		p_first_x = p_last_x = p_dummy_x;
		p_first_y = p_last_y = p_dummy_y;
		n_new = 0;
		k = 0;
		done = FALSE;
		first = TRUE;
		do {
			G = seg[id].group;
			L = seg[id].pos;
			np = seg[id].n;
			if (end_order == 0) {	/* Already in the right order */
				if (D->table[G]->segment[L]->coord[GMT_X][0] == p_last_x && D->table[G]->segment[L]->coord[GMT_Y][0] == p_last_y) {	/* Skip duplicate anchor point */
					j = 1;
					n = np - 1;
				}
				else {	/* We need all the points */
					j = 0;
					n = np;
				}
				if (D->table[G]->segment[L]->header)
					sprintf (GMT_io.segment_header, "# %s", D->table[G]->segment[L]->header);
				else
					sprintf (GMT_io.segment_header, "#\n");
				GMT_write_segmentheader (fp, D->table[G]->segment[L]->n_columns);
				Write_This_Segment (fp, D->table[G]->segment[L], j, np-1);
				p_last_x = D->table[G]->segment[L]->coord[GMT_X][np-1];
				p_last_y = D->table[G]->segment[L]->coord[GMT_Y][np-1];
				if (first) p_first_x = D->table[G]->segment[L]->coord[GMT_X][0], p_first_y = D->table[G]->segment[L]->coord[GMT_Y][0];
			}
			else {	/* Must reverse the segment's order of points */
				if (D->table[G]->segment[L]->coord[GMT_X][np-1] == p_last_x && D->table[G]->segment[L]->coord[GMT_Y][np-1] == p_last_y) {	/* Skip duplicate anchor point */
					j = 1;
					n = np - 1;
				}
				else {	/* We need all the points */
					j = 0;
					n = np;
				}
				if (D->table[G]->segment[L]->header)
					sprintf (GMT_io.segment_header, "# Reversed %s", D->table[G]->segment[L]->header);
				else
					sprintf (GMT_io.segment_header, "#\n");
				GMT_write_segmentheader (fp, D->table[G]->segment[L]->n_columns);
				Write_This_Segment (fp, D->table[G]->segment[L], np-1-j, (GMT_LONG)0);
				p_last_x = D->table[G]->segment[L]->coord[GMT_X][0];
				p_last_y = D->table[G]->segment[L]->coord[GMT_Y][0];
				if (first) p_first_x = D->table[G]->segment[L]->coord[GMT_X][np-1], p_first_y = D->table[G]->segment[L]->coord[GMT_Y][np-1];
			}
			first = FALSE;
			n_new += n;
			end_order = !end_order;
			seg[id].used = TRUE;
			if (seg[id].buddy[end_order].dist <= cutoff && !seg[seg[id].buddy[end_order].id].used) {
				/* Not done, trace into the next connecting segment */
				id2 = seg[id].buddy[end_order].id;
				end_order = seg[id].buddy[end_order].end_order;
				done = (id2 == start_id || id2 == id);
				id = id2;
			}
			else	/* End of the D->table for this segment */
				done = TRUE;
			k++;
		} while (!done);
		if (individual_file) GMT_fclose (fp);
		if (gmtdefs.verbose) fprintf (stderr, "%s: Segment %ld made from %ld pieces\n", GMT_program, out_seg, k);
		
		if (p_first_x == p_last_x && p_first_y == p_last_y) {
			if (individual_file && save_type) {	/* Ended up closed, rename with the C type */
				char filename2[BUFSIZ];
				sprintf (filename2, format, 'C', out_seg);
				rename (filename, filename2);
				strcpy (filename, filename2);
				mode = CLOSED;
			}
			n_closed++;
		}
		if (do_lists) fprintf (fpl[mode], "%s\n", filename);
		
		chain++;
		out_seg++;
		
		/* Wind to the next unused segments to start the connection search again */
		start_id = 0;
		while (start_id < ns && seg[start_id].used) start_id++;
		done = (start_id == ns);	/* No more unused segments */
	}
	if (do_lists) {
		fclose (fpl[CLOSED]);
		if (fpl[CLOSED] != fpl[OPEN]) fclose (fpl[OPEN]);
	}

	if (gmtdefs.verbose) {
		fprintf (stderr, "%s: Segments in: %ld Segments out: %ld\n", GMT_program, ns + n_islands, chain + n_islands);
		if (n_trouble) fprintf (stderr, "%s: %ld trouble spots\n", GMT_program, n_trouble);
		if (n_closed) fprintf (stderr, "%s: %ld new closed segments\n", GMT_program, n_closed);
		if (n_islands) fprintf (stderr, "%s: %ld were already closed\n", GMT_program, n_islands);
	}
	
	GMT_free ((void *)seg);
	GMT_free_dataset (D);

	exit (EXIT_SUCCESS);
}

GMT_LONG connect (struct LINK *S, int id, int order, double cutoff, GMT_LONG nn_check, double nn_dist)
{	/* Checks if OK to connect this segment to its nearest neighbor and returns TRUE if OK */

	if (S[S[id].buddy[order].id].used) return (FALSE);		/* Segment has been used already */
	if (S[id].buddy[order].dist > cutoff) return (FALSE);		/* Exceeds minimum gap */
	if (!nn_check) return (TRUE);					/* Passed all requirements */
	if (S[id].buddy[order].next_dist > nn_dist) return (TRUE);	/* Next neighboor is far enough away */
	return (FALSE);							/* Failed all tests */
}

void Write_This_Segment (FILE *fp, struct GMT_LINE_SEGMENT *line, GMT_LONG start, GMT_LONG end)
{
	GMT_LONG i, j, inc, done = FALSE;
	static double out[GMT_MAX_COLUMNS];
	
	inc = (start < end) ? +1 : -1;
	for (i = start; !done; i += inc) {
		for (j = 0; j < line->n_columns; j++) out[j] = line->coord[j][i];
		GMT_output (fp, line->n_columns, out);
		done = (i == end);
	}
}
