/*--------------------------------------------------------------------
 *	$Id: grdspotter.c,v 1.34 2011/07/11 19:22:06 guru Exp $
 *
 *   Copyright (c) 1999-2011 by P. Wessel
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; version 2 or any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   Contact info: www.soest.hawaii.edu/pwessel
 *--------------------------------------------------------------------*/
/*
 * GRDSPOTTER will (1) read a grid with topo/grav/whatever of seamounts,
 * (2) read an ascii file with stage or finite rotations, and (3)
 * convolve the flowline of each node with its prism volume, and
 * (4) build a cumulative volcano amplitude (CVA) grid.  Optionally, it
 * can also calculate the Data Importance (DI) associated with each node's
 * flowline and the predicted age (PA) for each node.  Finally, errors in
 * the location of the CVA maximum can be assessed by running a bootstrap
 * estimation.  The grids are all written out in GMT format and can be
 * processed and plotted with GMT.
 * GRDSPOTTER is part of the SPOTTER supplemental GMT package and should
 * be installed accordingly, see README.spotter.
 *
 * Author:	Paul Wessel, SOEST, Univ. of Hawaii, Honolulu, HI, USA
 * Date:	20-OCT-2007
 * Version:	2.0
 *
 *-------------------------------------------------------------------------
 * The Euler file must have following format:
 *
 * 1. Any number of comment lines starting with # in first column
 * 2. Any number of blank lines (just carriage return, no spaces)
 * 2. Any number of stage pole records which each have the format:
 *    lon(deg)  lat(deg)  tstart(Ma)  tstop(Ma)  ccw-angle(deg)
 * 3. stage records must go from oldest to youngest rotation
 * 4. Note tstart is larger (older) that tstop for each record
 * 5. No gaps allowed: tstart must equal the previous records tstop
 *
 * Example: Duncan & Clague [1985] Pacific-Hotspot rotations:
 *
 * # Time in Ma, angles in degrees
 * # lon  lat	tstart	tend	ccw-angle
 * 165     85	150	100	24.0
 * 284     36	100	74	15.0
 * 265     22	74	65	7.5
 * 253     17	65	42	14.0
 * 285     68	42	0	34.0
 *
 * If finite reconstruction poles are preferred, use the -C option and
 * give a file like DC85 poles in finite rotation format:
 *
 * #longitude	latitude	time(My)	angle(deg)
 * 285.00000	 68.00000	 42.0000	 34.0000
 * 275.66205	 53.05082	 65.0000	 43.5361
 * 276.02501	 48.34232	 74.0000	 50.0405
 * 279.86436	 46.30610	100.0000	 64.7066
 * 265.37800	 55.69932	150.0000	 82.9957
 *
 * Note that finite rotation files have one less column since each rotation
 * implicitly goes to time 0.
 *
 * Seamount data grid is a standard GMT grid file, presumably
 * representing topography or gravity of residual seamounts
 * (background residual should be removed first).
 *
 *------------------------------------------------------------------------------
 * REFERENCES:
 *
 * -> The hotspotting technique:
 *
 * Aslanian, D., L. Geli, and J.-L. Olivet, 1998, Hotspotting called into
 *    question, Nature, 396, 127.
 * Wessel, P., and L. W. Kroenke, 1997, A geometric technique for
 *    relocating hotspots and refining absolute plate motions, Nature,
 *    387, 365-369.
 * Wessel, P., and L. W. Kroenke, 1998a, Factors influencing the locations
 *    of hot spots determined by the hot-spotting technique, Geophys. Res.
 *    Lett., 25, 55-58.
 * Wessel, P., and L. W. Kroenke, 1998b, The geometric relationship between
 *    hot spots and seamounts: implications for Pacific hot spots, Earth
 *    Planet. Sci. Lett., 158, 1-18.
 * Wessel, P., and L. W. Kroenke, 1998c, Hotspotting called into question
 *    - Reply, Nature, 396, 127-128.
 *
 * -> Suitable Seamount data set:
 *
 * Sandwell, D, and W. H. F. Smith, 1995, JGR, ...
 * Smith, W. H. F., and D. Sandwell, 1997, Science, ...
 *
 * -> Plate motion models (stage poles):
 *
 * Duncan, R.A., and D. Clague, 1985, Pacific plate motion recorded by linear
 *    volcanic chains, in: A.E.M. Nairn, F. G. Stehli, S. Uyeda (eds.), The
 *    Ocean Basins and Margins, Vol. 7A, Plenum, New York, pp. 89-121.
 * Wessel and Kroenke, 1997, (see above)
 *
 */
 
#include "spotter.h"

#define B_TO_MB	(1.0 / 1048576.0)
#define OUTSIDE -1

GMT_LONG get_flowline (double xx, double yy, double tt, struct EULER *p, GMT_LONG n_stages, double d_km, GMT_LONG step, int forth_flag, double wesn[], double **flow);
GMT_LONG set_age (double *t_smt, float *age, GMT_LONG node, double upper_age, GMT_LONG truncate, GMT_LONG no_ages);

int main (int argc, char **argv)
{

	GMT_LONG n_nodes;			/* Number of nodes processed */
	GMT_LONG n_stages;			/* Number of stage rotations (poles) */
	GMT_LONG node;			/* The current node index */
	GMT_LONG one_or_zero;		/* Will be 1 or 0 depending on grid format (-F) */
	GMT_LONG n_try = 0;			/* Number of bootstrap estimates required */
	GMT_LONG try;			/* Number of current bootstrap estimate */
	GMT_LONG nm;			/* Number of nodes in the output CVA grid */
	GMT_LONG node_offset = 0;
	size_t n_alloc = 0, inc_alloc = BIG_CHUNK;
	
	unsigned short pa = 0;		/* Placeholder for PA along track */

	/* Misc. counters */

	GMT_LONG i, j, k, ij, m, row, col, k_step, np, max_ij = 0, n_flow, mem = 0, n_unique_nodes = 0, Qid;

	GMT_LONG error = FALSE;		/* TRUE when arguments are wrong */
	GMT_LONG no_ages = TRUE;		/* TRUE when no age grid is given */
	GMT_LONG truncate_ages = FALSE;	/* Truncate all nodes ages > upper_age */
	GMT_LONG finite = FALSE;		/* TRUE if stage pole file contains finite rotation poles instead */
	GMT_LONG get_DI = FALSE;		/* TRUE if we want to obtain Data Importance using max CVA along flowlines */
	GMT_LONG get_PA = FALSE;	/* TRUE if we want to use CVA flowlines to predict node age */
	GMT_LONG blabber = FALSE;	/* TRUE if we want excessive verbosity */
	GMT_LONG bootstrap = FALSE;	/* TRUE if we want bootstrap estimates */
	GMT_LONG keep_flowlines = FALSE;	/* TRUE if get_DI, get_PA, or bootstrap is TRUE */
	GMT_LONG multi_slices = FALSE;	/* TRUE if we want CVAs for many different z-sizes */
	GMT_LONG check_IDs = FALSE;	/* TRUE if we want to compute CVS for certain ID points only */
	GMT_LONG forth_flag;		/* Holds the do_time + 10 flag passed to forthtrack */
	GMT_LONG normalize = FALSE;	/* TRUE if we want CVA in percents */
	GMT_LONG *ID = NULL;			/* Optional array with IDs for each node */
	GMT_LONG slow_and_painful = FALSE;	/* When we run out of flowline memory */
	GMT_LONG set_fixed = FALSE;	/* Optional TRUE if we want to set all z[ij]'s to a fixed value */
	GMT_LONG *processed_node;		/* Pointer to array with TRUE/FALSE values for each grid node */
	
	float *CVA = NULL;			/* The CVA surface we seek to calculate */
	float *z = NULL;			/* The topo or grav input surface */
	float *age = NULL;			/* The age input surface */
	float *p_ages = NULL;			/* The age output surface [optional] */
	float *DI = NULL;			/* The Data Importance output surface [optional] */
	float zz, *tmp = NULL;

	double sampling_int_in_km;	/* Sampling interval along flowline (in km) */
	double upper_age = 180.0;	/* Upper age assigned to nodes on undated seafloor */
	double *x_smt = NULL;		/* node longitude (input degrees, stored as radians) */
	double *y_smt = NULL;		/* node latitude (input degrees, stored as radians) */
	double t_smt;			/* node upper age (up to age of seafloor) */
	double *c = NULL;		/* Array with one flowline */
	double CVA_west_rad, CVA_east_rad, CVA_south_rad, CVA_north_rad;	/* Region in radians */
	double CVA_xinc_rad, CVA_yinc_rad;		/* Grid spacing in radians */
	double offset, half;		/* 0.0 or 0.5, depending on -F */
	double z_min = 0.0, z_max = DBL_MAX;		/* Only deal with z-values inside this range [0/Inf] */
	double CVA_max, wesn[4], cva_contribution;
	double out[3], scale, z_inc, sdist = 0.0, area, yg;
	double *lat_area = NULL;	/* Area of each dx by dy note in km as function of latitude */
	double CVA_scale;		/* Used to normalize CVAs to percent */
	double this_wesn[4];
	double this_pa, pa_val = 0.0, n_more_than_once = 0.0;
	double *x_cva = NULL, *y_cva = NULL;		/* Coordinates on the CVA grid */
	double fixed_value = 0.0;

	char *CVA_file = CNULL;		/* Output name of CVA grid */
	char *Z_file = CNULL;		/* Input name of z grid */
	char *T_file = CNULL;		/* Input name of age grid */
	char *ID_file = CNULL;		/* Input name of ID grid */
	char *PA_file = CNULL;		/* Output name of predicted age grid */
	char *DI_file = CNULL;		/* Output name of optional improved z grid */
	char *Q_file = CNULL;		/* Input name of optional ID info file */
	char *euler_file = CNULL;	/* Input name of file with APM rotations */
	
	struct GRD_HEADER CVA_map;	/* Header structure for output CVA grid */
	struct GRD_HEADER zgrid;	/* Header structure for input topo/grav grid */
	struct GRD_HEADER tgrid;	/* Header structure for input age grid */
	struct GRD_HEADER idgrid;	/* Header structure for input ID grid */
	
	struct EULER *p = NULL;		/* Array of structures with Euler stage rotations */

	struct FLOWLINE *flowline = NULL;	/* Array with flowline structures */
	
	struct ID {			/* Information regarding one chain ID */
		double wesn[4];		/* Do not calculate flowlines outside this box */
		GMT_LONG ok;		/* TRUE if we want to calculate this CVA */
		GMT_LONG check_region;	/* TRUE if w, e, s, n is more restrictive than command line -R */
	} *ID_info = NULL;

	/* ------------------- END OF DECLARATIONS ------------------------------------------------------------*/

	argc = (int)GMT_begin (argc, argv);			/* Initialize GMT environment */
		
	GMT_grd_init (&CVA_map, argc, argv, FALSE);	/* Initialize grid structure */
	
	/* Check command line arguments */
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			
				/* Common parameters */
			
				case 'V':
				case 'R':
				case ':':
				case '\0':
					error += GMT_get_common_args (argv[i], &CVA_map.x_min, &CVA_map.x_max, &CVA_map.y_min, &CVA_map.y_max);
					break;
				
				/* Supplemental parameters */
				
				case 'A':
					T_file = &argv[i][2];
					break;
				case 'B':
					n_try = atoi (&argv[i][2]);
					bootstrap = TRUE;
					break;
				case 'C':	/* Use finite rotation poles */
					finite = TRUE;
					break;
				case 'E':
					euler_file = &argv[i][2];
					break;
				case 'F':
					CVA_map.node_offset = 1;
					break;
				case 'G':
					CVA_file = &argv[i][2];
					break;
				case 'D':
					DI_file = &argv[i][2];
					get_DI = TRUE;
					break;
				case 'P':
					PA_file = &argv[i][2];
					get_PA = get_DI = TRUE;
					break;
				case 'I':
					GMT_getinc (&argv[i][2], &CVA_map.x_inc, &CVA_map.y_inc);
					break;
				case 'L':
					ID_file = &argv[i][2];
					break;
				case 's':
					sdist = atof(&argv[i][2]);
					break;
				case 'M':
					slow_and_painful = TRUE;
					break;
				case 'N':
					upper_age = atof (&argv[i][2]);
					break;
				case 'Z':
					m = sscanf (&argv[i][2], "%lf/%lf/%lf", &z_min, &z_max, &z_inc);
					if (m == 1) z_max = 1.0e300;	/* Max not specified */
					if (m == 3) multi_slices = TRUE;	/* Want several slices */
					break;
				case 'Q':
					if (!access (&argv[i][2], R_OK)) {	/* The file exists */
						check_IDs = 2;
						Q_file = &argv[i][2];
					}
					else if ((Qid = atoi (&argv[i][2])) > 0) {	/* Got OK id value */
						check_IDs = 1;
					}
					else {
						fprintf (stderr, "%s: Error -Q: Must give valid file or ID value\n", GMT_program);
						error++;
					}
					break;
				case 'S':
					normalize = TRUE;
					break;
				case 'T':
					truncate_ages = TRUE;
					break;
				case 'U':
					fixed_value = atof (&argv[i][2]);;
					set_fixed = TRUE;
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			Z_file = argv[i];
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "%s %s - Create CVA image from surface grids\n\n", GMT_program, SPOTTER_VERSION);
		fprintf (stderr, "usage: %s zgrdfile -E<stage_file> -G<CVAgrid> -I<dx[m|c]>[/<dy>[m|c]]\n", GMT_program);
		fprintf (stderr, "\t-R<west/east/south/north> [-A<agegrid>] [-B<n_try] [-C] [-D<newzgrid] [-F] [-L<IDgrid>]\n");
		fprintf (stderr, "\t[-M] [-N<upper_age>] [-P<predict_agegrid>] [-Q<IDinfo>] [-S] [-T] [-U<val>] [-V] [-Z<z_min>[/<z_max>[/<z_inc>]]]\n\n");
		
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
		
		fprintf (stderr, "\t<zgrdfile> is the grid with topo or gravity\n");
		fprintf (stderr, "\t-E specifies the rotations to be used (see man page for format)\n\n");
		fprintf (stderr, "\t-G Specify file name for output CVA convolution grid.\n");
		fprintf (stderr, "\t-I specifies grid interval(s); Append m [or c] to <dx> and/or <dy> for minutes [or seconds].\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A co-registered grid with upper ages to use [Default is flowlines for all ages]\n");
		fprintf (stderr, "\t-B Get <n_try> bootstrap estimates of maximum CVA location [Default is no bootstrapping].\n");
		fprintf (stderr, "\t-C The file given with -E contains finite rotation poles [Default is stage poles]\n");
		fprintf (stderr, "\t-D Use flowlines to estimate data importance DI\n");
		fprintf (stderr, "\t-F Force pixel-registration [Default is gridline registration]\n");
		fprintf (stderr, "\t-L co-registered grid with chain ID for each node [Default ignores IDs]\n");
		fprintf (stderr, "\t-M Do flowline calculations as needed rather than storing in memory.\n");
		fprintf (stderr, "\t   You may have to use this option if -R is too large. Cannot be used with -B or -Z-slicing\n");
		fprintf (stderr, "\t-N sets upper age in m.y. for nodes whose plate age is NaN [180]\n");
		fprintf (stderr, "\t-P Use flowlines to estimate predicted ages at node locations\n");
		fprintf (stderr, "\t-Q Either single ID to use or file with list of IDs [Default uses all IDs]\n");
		fprintf (stderr, "\t   Each line would be TAG ID [w e s n] with optional zoom box\n");
		fprintf (stderr, "\t-S Normalize CVA grid to percentages of the CVA maximum\n");
		fprintf (stderr, "\t-T truncate all ages to max age in stage pole model [Default extrapolates]\n");
		fprintf (stderr, "\t-U After a node passes the -Z test, use this fixed value instead in CVA calculations\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t-Z ignore nodes with z-value lower than z_min [0] and optionally larger than z_max [Inf]\n");
		fprintf (stderr, "\t   Give z_min/z_max/z_inc to make CVA grids for each z-slice {Default makes 1 CVA grid]\n");
		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&CVA_map.x_inc, &CVA_map.y_inc, &node_offset, NULL);
	CVA_map.node_offset = (int)node_offset;

	if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (CVA_map.x_inc <= 0.0 || CVA_map.y_inc <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
		error = TRUE;
	}
	if (!CVA_file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G:  Must specify output file\n", GMT_program);
		error = TRUE;
	}
	if (!Z_file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Must give name of topo gridfile\n", GMT_program);
		error++;
	}
	if (ID_file && !check_IDs) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Must specify both -L and -Q if one is present\n", GMT_program);
		error++;
	}
	if (slow_and_painful && (bootstrap || multi_slices)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Cannot use -M with -B or -Z (slicing)\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	/* ------------------- END OF PROCESSING COMMAND LINE ARGUMENTS  --------------------------------------*/

	GMT_lat_swap_init ();	/* Initialize auxiliary latitude machinery */

	/* Load in the Euler stage poles */

	n_stages = spotter_init (euler_file, &p, TRUE, finite, FALSE, &upper_age, gmtdefs.verbose);

	/* Initialize the CVA grid and structure */

	CVA_xinc_rad = CVA_map.x_inc * D2R;
	CVA_yinc_rad = CVA_map.y_inc * D2R;

	if (CVA_map.node_offset) {	/* Pixel grid */
		one_or_zero = 0;
		offset = 0.0;
	}
	else {
		one_or_zero = 1;
		offset = 0.5;
	}

	CVA_map.nx = (int)(irint ((CVA_map.x_max - CVA_map.x_min) / CVA_map.x_inc) + one_or_zero);
	CVA_map.ny = (int)(irint ((CVA_map.y_max - CVA_map.y_min) / CVA_map.y_inc) + one_or_zero);
	nm = GMT_get_nm (CVA_map.nx, CVA_map.ny);
	CVA = (float *) GMT_memory (VNULL, (size_t)nm, sizeof(float), GMT_program);

	/* Assign grid-region variables in radians to avoid conversions inside convolution loop */

	CVA_west_rad   = CVA_map.x_min * D2R;
	CVA_east_rad   = CVA_map.x_max * D2R;
	CVA_south_rad  = CVA_map.y_min * D2R;
	CVA_north_rad  = CVA_map.y_max * D2R;
	wesn[0] = CVA_west_rad;		wesn[1] = CVA_east_rad;
	wesn[2] = CVA_south_rad;	wesn[3] = CVA_north_rad;

	/* Allocate T/F array */

	processed_node = (GMT_LONG *) GMT_memory (VNULL, nm, sizeof (GMT_LONG), GMT_program);

	/* Set flowline sampling interval to 1/2 of the shortest distance between x-nodes */

	/* sampling_int_in_km = 0.5 * CVA_xinc_rad * EQ_RAD * ((fabs (CVA_north_rad) > fabs (CVA_south_rad)) ? cos (CVA_north_rad) : cos (CVA_south_rad)); */
	sampling_int_in_km = CVA_xinc_rad * EQ_RAD * ((fabs (CVA_north_rad) > fabs (CVA_south_rad)) ? cos (CVA_north_rad) : cos (CVA_south_rad));
	if (sdist != 0.0) sampling_int_in_km = sdist;
	if (gmtdefs.verbose) fprintf (stderr, "%s: Flowline sampling interval = %.3f km\n", GMT_program, sampling_int_in_km);

	if (truncate_ages && gmtdefs.verbose) fprintf (stderr, "%s: Ages truncated to %g\n", GMT_program, upper_age);

	/* Start to read input data */
	
	if (GMT_read_grd_info (Z_file, &zgrid)) {
		fprintf (stderr, "%s: Error opening file %s\n", GMT_program, Z_file);
		exit (EXIT_FAILURE);
	}
	half = (zgrid.node_offset) ? 0.5 : 0.0;
	z = (float *) GMT_memory (VNULL, (size_t)(zgrid.nx * zgrid.ny), sizeof(float), GMT_program);

	if (GMT_read_grd (Z_file, &zgrid, z, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE)) {
		fprintf (stderr, "%s: Error reading file %s\n", GMT_program, Z_file);
		exit (EXIT_FAILURE);
	}
	area = 111.195 * zgrid.y_inc * 111.195 * zgrid.x_inc;	/* In km^2 at Equator */
	x_smt = (double *) GMT_memory (VNULL, (size_t)zgrid.nx, sizeof(double), GMT_program);
	for (i = 0; i < zgrid.nx; i++) x_smt[i] = D2R * GMT_i_to_x (i, zgrid.x_min, zgrid.x_max, zgrid.x_inc, half, zgrid.nx);
	y_smt = (double *) GMT_memory (VNULL, (size_t)zgrid.ny, sizeof(double), GMT_program);
	for (j = 0; j < zgrid.ny; j++) y_smt[j] = D2R * GMT_lat_swap (GMT_j_to_y (j, zgrid.y_min, zgrid.y_max, zgrid.y_inc, half, zgrid.ny), GMT_LATSWAP_G2O);	/* Convert to geocentric */
	lat_area = (double *) GMT_memory (VNULL, (size_t)zgrid.ny, sizeof(double), GMT_program);

	for (j = 0; j < zgrid.ny; j++) lat_area[j] = area * cos (y_smt[j]);
	
	x_cva = (double *) GMT_memory (VNULL, (size_t)CVA_map.nx, sizeof(double), GMT_program);
	for (i = 0; i < CVA_map.nx; i++) x_cva[i] = GMT_i_to_x (i, CVA_map.x_min, CVA_map.x_max, CVA_map.x_inc, half, CVA_map.nx);
	y_cva = (double *) GMT_memory (VNULL, (size_t)CVA_map.ny, sizeof(double), GMT_program);
	for (j = 0; j < CVA_map.ny; j++) y_cva[j] = GMT_j_to_y (j, CVA_map.y_min, CVA_map.y_max, CVA_map.y_inc, half, CVA_map.ny);
	if (T_file) {
		if (GMT_read_grd_info (T_file, &tgrid)) {
			fprintf (stderr, "%s: Error opening file %s\n", GMT_program, T_file);
			exit (EXIT_FAILURE);
		}
		if (!(tgrid.nx == zgrid.nx && tgrid.ny == zgrid.ny && tgrid.x_min == zgrid.x_min && tgrid.y_min == zgrid.y_min)) {
			fprintf (stderr, "%s: GMT ERROR:  topo grid and age grid must coregister\n", GMT_program);
			exit (EXIT_FAILURE);
		}
		
		age = (float *) GMT_memory (VNULL, (size_t)(tgrid.nx * tgrid.ny), sizeof(float), GMT_program);

		if (GMT_read_grd (T_file, &tgrid, age, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE)) {
			fprintf (stderr, "%s: Error reading file %s\n", GMT_program, T_file);
			exit (EXIT_FAILURE);
		}
		no_ages = FALSE;
	}
	if (ID_file) {
		if (GMT_read_grd_info (ID_file, &idgrid)) {
			fprintf (stderr, "%s: Error opening file %s\n", GMT_program, ID_file);
			exit (EXIT_FAILURE);
		}
		if (!(idgrid.nx == zgrid.nx && idgrid.ny == zgrid.ny && idgrid.x_min == zgrid.x_min && idgrid.y_min == zgrid.y_min)) {
			fprintf (stderr, "%s: GMT ERROR:  topo grid and ID grid must coregister\n", GMT_program);
			exit (EXIT_FAILURE);
		}
		
		ID = (GMT_LONG *) GMT_memory (VNULL, (size_t)(idgrid.nx * idgrid.ny), sizeof (GMT_LONG), GMT_program);
		tmp = (float *) GMT_memory (VNULL, (size_t)(idgrid.nx * idgrid.ny), sizeof(float), GMT_program);

		if (GMT_read_grd (ID_file, &idgrid, tmp, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE)) {
			fprintf (stderr, "%s: Error reading file %s\n", GMT_program, ID_file);
			exit (EXIT_FAILURE);
		}
		for (i = 0; i < idgrid.nx * idgrid.ny; i++) ID[i] = (GMT_LONG) rint ((double)tmp[i]);
		GMT_free ((void *)tmp);
		ID_info = (struct ID *) GMT_memory (VNULL, (size_t)(rint (idgrid.z_max) + 1), sizeof (struct ID), GMT_program);
		if (check_IDs == 1) {	/* Only doing one CVA with no extra restrictions */
			ID_info[Qid].ok = TRUE;	/* Every other info in struct array is NULL or 0 */
		}
		else {	/* Must read from file */
			FILE *fp;
			double wq, eq, sq, nq;
			char line[BUFSIZ];
			
			if ((fp = fopen (Q_file, "r")) == NULL) {	/* Oh, oh... */
				fprintf (stderr, "%s: GMT ERROR:  -Q info file unreadable/nonexistent\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			while (fgets (line, BUFSIZ, fp)) {
				GMT_chop (line);					/* Rid the world of CR/LF */
				if (line[0] == '#' || line[0] == '\0') continue;
				k = sscanf (line, "%*s %" GMT_LL "d %lf %lf %lf %lf", &Qid, &wq, &eq, &sq, &nq);
				ID_info[Qid].ok = TRUE;
				if (k == 5) {	/* Got restricted wesn also */
					ID_info[Qid].check_region = TRUE;
					ID_info[Qid].wesn[0] = wq * D2R;
					ID_info[Qid].wesn[1] = eq * D2R;
					ID_info[Qid].wesn[2] = sq * D2R;
					ID_info[Qid].wesn[3] = nq * D2R;
				}
			}
			fclose (fp);
		}
		check_IDs = TRUE;
	}

	if (slow_and_painful) {
		keep_flowlines = FALSE;	/* Do it the hard way to save memory */
		k_step = 2;		/* Reset back to 3 if needed later */
		forth_flag = 10;	/* The 10 is used to limit the flowline calculation to only resample track within the rectangular box of interest */
	}
	else {
		keep_flowlines = (get_PA || get_DI || bootstrap);
		forth_flag = get_PA + 10;	/* The 10 is used to limit the flowline calculation to only resample track within the rectangular box of interest */
		k_step = (get_PA) ? 3 : 2;
	}
	if (keep_flowlines) {
		n_alloc = inc_alloc;
		flowline = (struct FLOWLINE *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (struct FLOWLINE), GMT_program);
		if (gmtdefs.verbose) {
			fprintf (stderr, "%s: Will attempt to keep all flowlines in memory.  However, should this not be possible\n", GMT_program);
			fprintf (stderr, "the program might crash.  If so consider using the -M option\n");
		}
	}
	
	for (row = ij = n_flow = n_nodes = 0; row < zgrid.ny; row++) {	/* Loop over all input rows */
		for (col = 0; col < zgrid.nx; col++, ij++) {	/* Loop over all input columns */
			/* if (! (col == 1267 && row == 840)) continue; */
			/* STEP 1: Determine if z exceeds threshold and if so assign age */
			if (GMT_is_fnan (z[ij]) || z[ij] < z_min || z[ij] > z_max) continue;	/* Skip node since it is NaN or outside the z_min < z < z_max range */
			if (check_IDs && !ID_info[(GMT_LONG)ID[ij]].ok) continue;			/* Skip because of wrong ID */
			if (!set_age (&t_smt, age, ij, upper_age, truncate_ages, no_ages)) continue;
			/* STEP 2: Calculate this node's flowline */
			if (check_IDs && ID_info[(GMT_LONG)ID[ij]].check_region) /* Set up a box-limited flowline sampling */
				memcpy ((void *)this_wesn, (void *)ID_info[(GMT_LONG)ID[ij]].wesn, 4*sizeof(double));
			else
				memcpy ((void *)this_wesn, (void *)wesn, 4*sizeof(double));
			np = get_flowline (x_smt[col], y_smt[row], t_smt, p, n_stages, sampling_int_in_km, k_step, (int)forth_flag, this_wesn, &c);
			if (np == 0) continue;	/* No flowline inside this wesn */

			/* STEP 3: Convolve this flowline with node shape and add to CVA grid */

			/* Our convolution is approximate:  We sample the flowline frequently and use
			 * one of the points on the flowline that are closest to the node.  Ideally we
			 * want the nearest distance from each node to the flowline. Later versions may
			 * improve on this situation */

			memset ((void *)processed_node, FALSE, nm*sizeof(GMT_LONG));	/* Fresh start for this flowline convolution */
			
			if (keep_flowlines) {
				flowline[n_nodes].node = (GMT_LONG *) GMT_memory (VNULL, (size_t)np, sizeof (GMT_LONG), GMT_program);
				if (get_PA) flowline[n_nodes].PA = (unsigned short *) GMT_memory (VNULL, (size_t)np, sizeof (unsigned short), GMT_program);
			}
			
			/* Keep in mind that although first and last are entry/exit into the CVA grid, we must
			 * expect the flowline between those points to also exit/reentry the CVA grid; hence
			 * we must still check if points are in/out of the box.  Here, we do not skip points but
			 * set the node index to OUTSIDE */
			 
			cva_contribution = lat_area[row] * (set_fixed ? fixed_value : z[ij]);	/* This node's contribution to the convolution */

#ifdef DEBUG2
			printf ("> %ld %ld %ld %ld %g\n", n_nodes, np, row, col, cva_contribution);
#endif
			for (m = 0, k = 1; m < np; m++) {	/* Store nearest node indices only */
				i = GMT_x_to_i (c[k], CVA_west_rad,  CVA_xinc_rad, offset, CVA_map.nx);	k++;
				yg = GMT_lat_swap (R2D*c[k], GMT_LATSWAP_O2G);	/* Convert back to geodetic */
				j = GMT_y_to_j (yg, CVA_map.y_min, CVA_map.y_inc, offset, CVA_map.ny);	k++;
				if (i < 0 || i >= CVA_map.nx || j < 0 || j >= CVA_map.ny)	/* Outside the CVA box, flag as outside */
					node = OUTSIDE;
				else								/* Inside the CVA box, assign node ij */
					node = j * CVA_map.nx + i;
				if (keep_flowlines) {
					flowline[n_nodes].node[m] = node;
					if (get_PA) flowline[n_nodes].PA[m] = (unsigned short) irint (c[k++] * T_2_PA);
				}
				/* If we did not keep flowlines then there is no PA to skip over (hence no k++) */
				if (node != OUTSIDE) {
					if (!processed_node[node]) {	/* Have not added to the CVA at this node yet */
						CVA[node] += (float)cva_contribution;
						processed_node[node] = TRUE;		/* Now we have visited this node */
						n_unique_nodes++;
#ifdef DEBUG2
						printf ("%g\t%g\n", x_cva[i],y_cva[j]);
#endif
					}
					n_more_than_once += 1.0;
				}
			}
			if (keep_flowlines) {
				flowline[n_nodes].n = np;	/* Number of points in flowline */
				flowline[n_nodes].ij = ij;	/* Originating node in topo grid */
				mem += sizeof (struct FLOWLINE) + np * sizeof (GMT_LONG) + ((get_PA) ? np * sizeof (unsigned short) : 0);
			}

			GMT_free ((void *)c);	/* Free the flowline vector */

			n_nodes++;	/* Go to next node */

			if (keep_flowlines && n_nodes == (GMT_LONG)n_alloc) {
				inc_alloc *= 2;
				n_alloc += inc_alloc;
				flowline = (struct FLOWLINE *) GMT_memory ((void *)flowline, (size_t)n_alloc, sizeof (struct FLOWLINE), GMT_program);
			}
			
			if (gmtdefs.verbose && !(n_nodes%100)) fprintf (stderr, "%s: Row %5ld Processed %5ld nodes [%5ld/%.1f]\r", GMT_program, row, n_nodes, n_flow, mem * B_TO_MB);
		}
		if (gmtdefs.verbose) fprintf (stderr, "%s: Row %5ld Processed %5ld nodes [%5ld/%.1f]\r", GMT_program, row, n_nodes, n_flow, mem * B_TO_MB);
	}
	if (gmtdefs.verbose) fprintf (stderr, "%s: Row %5ld Processed %5ld nodes [%5ld/%.1f]\n", GMT_program, row, n_nodes, n_flow, mem * B_TO_MB);
	if (gmtdefs.verbose) fprintf (stderr, "%s: On average, each node was visited %g times\n", GMT_program, n_more_than_once / n_unique_nodes);

	if (keep_flowlines && n_nodes != (GMT_LONG)n_alloc) {
		flowline = (struct FLOWLINE *) GMT_memory ((void *)flowline, (size_t)n_nodes, sizeof (struct FLOWLINE), GMT_program);
	}
	
	/* OK, Done processing, time to write out */

	if (normalize) {	/* Convert CVA values to percent of CVA maximum */		
		if (gmtdefs.verbose) fprintf (stderr, "%s: Normalize CVS grid to percentages of max CVA\n", GMT_program);
		CVA_map.z_min = +DBL_MAX;
		CVA_map.z_max = -DBL_MAX;
		for (node = 0; node < (GMT_LONG)nm; node++) {
			if (GMT_is_fnan (CVA[node])) continue;
			if (CVA[node] < CVA_map.z_min) CVA_map.z_min = CVA[node];
			if (CVA[node] > CVA_map.z_max) CVA_map.z_max = CVA[node];
		}
		if (gmtdefs.verbose) fprintf (stderr, "%s: CVA min/max: %g %g -> ", GMT_program, CVA_map.z_min, CVA_map.z_max);
		CVA_scale = 100.0 / CVA_map.z_max;
		for (node = 0; node < (GMT_LONG)nm; node++) CVA[node] *= (float)CVA_scale;
		CVA_map.z_min *= CVA_scale;
		CVA_map.z_max *= CVA_scale;
		if (gmtdefs.verbose) fprintf (stderr, "%g %g\n", CVA_map.z_min, CVA_map.z_max);
	}
		
	if (gmtdefs.verbose) fprintf (stderr, "%s: Write CVA grid %s\n", GMT_program, CVA_file);

	if (GMT_write_grd (CVA_file, &CVA_map, CVA, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE)) {
		fprintf (stderr, "%s: Error writing file %s\n", GMT_program, CVA_file);
		exit (EXIT_FAILURE);
	}

	if (multi_slices) {	/* Do CVA calculations for each z-slice using stored flowlines */
		GMT_LONG nz;
		char file[256], format[256];
		double z0, z1;
		float *CVA_inc;
		
		if (gmtdefs.verbose) fprintf (stderr, "%s: Start z-slice CVA calculations\n", GMT_program);
		for (i = strlen (CVA_file); i >= 0 && CVA_file[i] != '.'; i--);
		if (CVA_file[i] == '.') {	/* Make a filename template from the CVA filename using the period as delimeter */
			strncpy (format, CVA_file, (size_t)i);		/* Should keep the prefix from a file called prefix.ext */
			strcat (format, "_%%ld");		/* Make filenames like prefix_#.ext */
			strcat (format, &CVA_file[i]);		/* Should add the extension from said file */
		}
		CVA_inc = (float *) GMT_memory (VNULL, nm, sizeof(float), GMT_program);
		nz = irint ((z_max - z_min) / z_inc);
		for (i = 0; i < nz; i++) {
			z0 = z_min + i * z_inc;
			z1 = z0 + z_inc;
			if (gmtdefs.verbose) fprintf (stderr, "%s: Start z-slice %g - %g\n", GMT_program, z0, z1);
			memset ((void *)CVA_inc, 0, nm * sizeof (float));	/* Fresh start for this z-slice */
			for (m = 0; m < n_nodes; m++) {				/* Loop over all active flowlines */
				ij = flowline[m].ij;
				if (z[ij] <= z0 || z[ij] > z1) continue;	/* z outside current slice */
				cva_contribution = lat_area[ij/zgrid.nx] * (set_fixed ? fixed_value : z[ij]);	/* This node's contribution to the convolution */
				memset ((void *)processed_node, FALSE, nm);			/* Fresh start for this flowline convolution */
				for (k = 0; k < flowline[m].n; k++) {			/* For each point along this flowline */
					node = flowline[m].node[k];
					if (node != OUTSIDE && !processed_node[node]) {	/* Have not added to the CVA at this node yet */
						CVA_inc[node] += (float)cva_contribution;
						processed_node[node] = TRUE;		/* Now we have visited this node */
					}
				}
				if (blabber && !(m%10000)) fprintf (stderr, "%s: Processed %5ld flowlines\r", GMT_program, m);
			}
			if (blabber) fprintf (stderr, "%s: Processed %5ld flowlines\n", GMT_program, n_nodes);
			
			/* Time to write out this z-slice grid */
			if (normalize) {	/* Convert CVA values to percent of CVA maximum */		
				if (gmtdefs.verbose) fprintf (stderr, "%s: Normalize CVS grid to percentages of max CVA\n", GMT_program);
				CVA_map.z_min = +DBL_MAX;
				CVA_map.z_max = -DBL_MAX;
				for (node = 0; node < (GMT_LONG)nm; node++) {
					if (GMT_is_fnan (CVA_inc[node])) continue;
					if (CVA_inc[node] < CVA_map.z_min) CVA_map.z_min = CVA_inc[node];
					if (CVA_inc[node] > CVA_map.z_max) CVA_map.z_max = CVA_inc[node];
				}
				if (gmtdefs.verbose) fprintf (stderr, "%s: CVA min/max: %g %g -> ", GMT_program, CVA_map.z_min, CVA_map.z_max);
				CVA_scale = 100.0 / CVA_map.z_max;
				for (node = 0; node < (GMT_LONG)nm; node++) CVA_inc[node] *= (float)CVA_scale;
				CVA_map.z_min *= CVA_scale;
				CVA_map.z_max *= CVA_scale;
			}
			sprintf (CVA_map.remark, "CVA for z-range %g - %g only", z0, z1);
			sprintf (file, format, i);
			if (gmtdefs.verbose) fprintf (stderr, "%s: Save z-slice CVA to file %s\n", GMT_program, file);
			if (GMT_write_grd (file, &CVA_map, CVA_inc, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE)) {
				fprintf (stderr, "%s: Error writing file %s\n", GMT_program, file);
				exit (EXIT_FAILURE);
			}
		}
		GMT_free ((void *)CVA_inc);
	}
			
	if (get_DI || get_PA) {	/* Must determine max CVA along each flowline */
		if (get_DI) DI = (float *) GMT_memory (VNULL, (size_t)(zgrid.nx * zgrid.ny), sizeof(float), GMT_program);
		if (get_PA) p_ages = (float *) GMT_memory (VNULL, (size_t)(zgrid.nx * zgrid.ny), sizeof(float), GMT_program);
		if (gmtdefs.verbose) fprintf (stderr, "%s: Compute DI and/or PA grids\n", GMT_program);

		if (keep_flowlines) {
			for (m = 0; m < n_nodes; m++) {	/* Loop over all active flowlines */
				CVA_max = 0.0;						/* Fresh start for this flowline convolution */
				for (k = 0; k < flowline[m].n; k++) {			/* For each point along this flowline */
					node = flowline[m].node[k];
					if (node != OUTSIDE && CVA[node] > CVA_max) {	/* Found a higher CVA value */
						CVA_max = CVA[node];
						if (get_PA) pa = flowline[m].PA[k];	/* Keep the estimate of age at highest CVA */
					}
				}
				if (get_DI) DI[flowline[m].ij] = (float)CVA_max;	/* Store the maximum CVA associated with this node's flowline */
				if (get_PA) p_ages[flowline[m].ij] = (float) (pa * PA_2_T);
				if (blabber && !(m%10000)) fprintf (stderr, "%s: Processed %5ld flowlines\r", GMT_program, m);
			}
			if (blabber && gmtdefs.verbose) fprintf (stderr, "%s: Processed %5ld flowlines\n", GMT_program, n_nodes);
		}
		else {	/* Must recreate flowlines */
			k_step = 3;	/* FLowlines have (x,y,t) here */
			forth_flag = get_PA + 10;	/* The 10 is used to limit the flowline calculation to only resample track within the rectangular box of interest */
			for (row = ij = n_flow = n_nodes = 0; row < zgrid.ny; row++) {	/* Loop over all input rows */
				for (col = 0; col < zgrid.nx; col++, ij++) {	/* Loop over all input columns */
					/* STEP 1: Determine if z exceeds threshold and if so assign age */
					if (GMT_is_fnan (z[ij]) || z[ij] < z_min || z[ij] > z_max) continue;	/* Skip node since it is NaN or outside the z_min < z < z_max range */
					if (check_IDs && !ID_info[(GMT_LONG)ID[ij]].ok) continue;			/* Skip because of wrong ID */
					if (!set_age (&t_smt, age, ij, upper_age, truncate_ages, no_ages)) continue;
					/* STEP 2: Calculate this node's flowline */
					if (check_IDs && ID_info[(GMT_LONG)ID[ij]].check_region) /* Set up a box-limited flowline sampling */
						memcpy ((void *)this_wesn, (void *)ID_info[(GMT_LONG)ID[ij]].wesn, 4*sizeof(double));
					else
						memcpy ((void *)this_wesn, (void *)wesn, 4*sizeof(double));
					np = get_flowline (x_smt[col], y_smt[row], t_smt, p, n_stages, sampling_int_in_km, k_step, (int)forth_flag, this_wesn, &c);
					if (np == 0) continue;	/* No flowline inside this wesn */
			 		n_nodes++;
					/* Fresh start for this flowline convolution */
					CVA_max = 0.0;
					this_pa = GMT_d_NaN;
					for (m = 0, k = 1; m < np; m++) {	/* Store nearest node indices only */
						i = GMT_x_to_i (c[k], CVA_west_rad,  CVA_xinc_rad, offset, CVA_map.nx);	k++;
						yg = GMT_lat_swap (R2D*c[k], GMT_LATSWAP_O2G);	/* Convert back to geodetic */
						j = GMT_y_to_j (yg, CVA_map.y_min, CVA_map.y_inc, offset, CVA_map.ny);	k++;
						if (get_PA) pa_val = c[k++];
						if (i < 0 || i >= CVA_map.nx || j < 0 || j >= CVA_map.ny) continue;	/* Outside the CVA box, flag as outside */
						node = j * CVA_map.nx + i;
						if (CVA[node] <= CVA_max) continue;	/* Already seen higher CVA values */
						CVA_max = CVA[node];
						if (get_PA) this_pa = pa_val;
					}
					if (get_DI) DI[ij] = (float)CVA_max;	/* Store the maximum CVA associated with this node's flowline */
					if (get_PA) p_ages[ij] = (float) this_pa;
					GMT_free ((void *)c);
				}
				if (blabber) fprintf (stderr, "%s: Row %5ld: Processed %5ld flowlines\r", GMT_program, row, n_nodes);
			}
			if (blabber) fprintf (stderr, "%s: Row %5ld: Processed %5ld flowlines\n", GMT_program, row, n_nodes);
		}
		

		if (get_DI) {
			if (gmtdefs.verbose) fprintf (stderr, "%s: Write DI grid %s\n", GMT_program, DI_file);
			GMT_grd_init (&zgrid, argc, argv, TRUE);
			sprintf (zgrid.remark, "CVA maxima along flowlines from each node");
			if (GMT_write_grd (DI_file, &zgrid, DI, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE)) {
				fprintf (stderr, "%s: Error writing file %s\n", GMT_program, DI_file);
				exit (EXIT_FAILURE);
			}
			GMT_free ((void *)DI);
		}
		if (get_PA) {
			if (gmtdefs.verbose) fprintf (stderr, "%s: Write PA grid %s\n", GMT_program, PA_file);
			GMT_grd_init (&zgrid, argc, argv, TRUE);
			sprintf (zgrid.remark, "Predicted age for each node");
			if (GMT_write_grd (PA_file, &zgrid, p_ages, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE)) {
				fprintf (stderr, "%s: Error writing file %s\n", GMT_program, PA_file);
				exit (EXIT_FAILURE);
			}
			GMT_free ((void *)p_ages);
		}
	}
	
	if (bootstrap) {	/* Use bootstrapping to estimate confidence region for CVA maxima */

		if (gmtdefs.verbose) {
			fprintf (stderr, "%s: Preprocessed %5ld flowlines\n", GMT_program, n_nodes);
			fprintf (stderr, "%s: %ld of %ld total flowlines entered CVA region\n", GMT_program, n_nodes, n_flow);
			fprintf (stderr, "%s: Flowlines consumed %ld Mb of memory\n", GMT_program, (GMT_LONG)irint (mem * B_TO_MB));
			fprintf (stderr, "%s: Estimate %ld CVA max locations using bootstrapping\n", GMT_program, n_try);
		}

		/* Now do bootstrap sampling of flowlines */
	
		try = 0;
		srand ((unsigned int)time(NULL));	/* Initialize random number generator */
		scale = (double)n_nodes / (double)RAND_MAX;
		for (try = 1; try <= n_try; try++) {
			if (gmtdefs.verbose) fprintf (stderr, "%s: Bootstrap try %ld\r", GMT_program, try);
		
			memset ((void *)CVA, 0, nm * sizeof (float));	/* Start with fresh grid */
			for (m = 0; m < n_nodes; m++) {	/* Loop over all indices */
				ij = (GMT_LONG)floor (rand() * scale);		/* Get a random integer in 0 to n_nodes-1 range */
				memset ((void *)processed_node, FALSE, nm);		/* Fresh start for this flowline convolution */
				zz = z[flowline[ij].ij];
				for (k = 0; k < flowline[ij].n; k++) {		/* For each point along this flowline */
					node = flowline[ij].node[k];
					if (node != OUTSIDE && !processed_node[node]) {	/* Have not added to the CVA at this node yet */
						CVA[node] += zz;
						processed_node[node] = TRUE;		/* Now we have visited this node; flag it */
					}
				}
				if (blabber && !(m%10000)) fprintf (stderr, "%s: Processed %5ld flowlines\r", GMT_program, m);
			}
			if (blabber) fprintf (stderr, "%s: Processed %5ld flowlines\n", GMT_program, n_nodes);
		
			/* Find max CVA location */
		
			CVA_max = 0.0;
			for (ij = 0; ij < (GMT_LONG)nm; ij++) {		/* Loop over all CVA nodes */
				if (CVA[ij] > CVA_max) {	/* Update new max location */
					CVA_max = CVA[ij];
					max_ij = ij;
				}
			}
			i = max_ij % CVA_map.nx;
			j = max_ij / CVA_map.nx;
			out[0] = GMT_i_to_x (i, CVA_map.x_min, CVA_map.x_max, CVA_map.x_inc, half, CVA_map.nx);
			out[1] = GMT_j_to_y (j, CVA_map.y_min, CVA_map.y_max, CVA_map.y_inc, half, CVA_map.ny);
			out[2] = CVA_max;
			
			GMT_output (GMT_stdout, 3, out);
		}
		if (gmtdefs.verbose) fprintf (stderr, "%s: Bootstrap try %ld\n", GMT_program, n_try);
	}
	
	/* Clean up memory */

	GMT_free ((void *)CVA);
	GMT_free ((void *)processed_node);
	GMT_free ((void *)z);
	for (i = 0; keep_flowlines && i < n_nodes; i++) {
		GMT_free ((void *)flowline[i].node);
		if (get_PA) GMT_free ((void *)flowline[i].PA);
	}
	GMT_free ((void *)p);
	GMT_free ((void *)x_smt);
	GMT_free ((void *)y_smt);
	if (!no_ages) GMT_free ((void *)age);
	if (keep_flowlines) GMT_free ((void *)flowline);
	if (check_IDs) GMT_free ((void *)ID_info);

	if (gmtdefs.verbose) fprintf (stderr, "%s: Done\n", GMT_program);

	GMT_end (argc, argv);
	
	exit (EXIT_SUCCESS);
}

GMT_LONG get_flowline (double xx, double yy, double tt, struct EULER *p, GMT_LONG n_stages, double d_km, GMT_LONG step, int flag, double wesn[], double **flow)
{
	GMT_LONG n_chunk, n_track, m, kx, ky, first, last, np;
	double *c, *f, y_min, y_max;
	
	yy = D2R*GMT_lat_swap (R2D*yy, GMT_LATSWAP_G2O);		/* Convert to geocentric */
	y_min = D2R*GMT_lat_swap (R2D*wesn[2], GMT_LATSWAP_G2O);	/* Convert to geocentric */
	y_max = D2R*GMT_lat_swap (R2D*wesn[3], GMT_LATSWAP_G2O);	/* Convert to geocentric */
	/* Get the flowline from this point back to time tt, restricted to the given wesn box */
	n_chunk = spotter_forthtrack (&xx, &yy, &tt, (GMT_LONG)1, p, n_stages, d_km, 0.0, flag, wesn, &c);

	n_track = irint (c[0]);				/* Number of point pairs making up this flowline */

	/* Find the first point on the flowline inside the desired CVA region */
			
	for (m = 0, ky = 2, first = -1; m < n_track && first == -1; m++, ky += step) {	/* For each point along flowline */
		if (c[ky] < y_min || c[ky] > y_max) continue;		/* Latitude outside region */
		kx = ky - 1;						/* Index for the x-coordinate */
		while (c[kx] > wesn[1]) c[kx] -= TWO_PI;		/* Elaborate W/E test because of 360 periodicity */
		while (c[kx] < wesn[0]) c[kx] += TWO_PI;
		if (c[kx] > wesn[1]) continue;				/* Longitude outside region */
		first = kx;						/* We are inside, this terminates the for loop */
	}
			
	if (first == -1) { 	/* Was never inside the grid, skip the entire flowline and move on */
		GMT_free ((void *)c);	/* Free the flowline vector */
		return 0;
	}
			
	/* Here we know searching from the end will land inside the grid eventually so last can never exit as -1 */
			
	for (m = n_track - 1, ky = step * m + 2, last = -1; m >= 0 && last == -1; m--, ky -= step) {	/* For each point along flowline */
		if (c[ky] < y_min || c[ky] > y_max) continue;		/* Latitude outside region */
		kx = ky - 1;						/* Index for the x-coordinate */
		while (c[kx] > wesn[1]) c[kx] -= TWO_PI;		/* Elaborate W/E test because of 360 periodicity */
		while (c[kx] < wesn[0]) c[kx] += TWO_PI;
		if (c[kx] > wesn[1]) continue;				/* Longitude outside region */
		last = kx;						/* We are inside, this terminates the for loop */
	}
			
	np = (last - first) / step + 1;			/* Number of (x,y[,t]) points on this flowline inside the region */
	if (np < n_track) {	/* Just copy out the subset of points we want */
		GMT_LONG n_alloc;
		n_alloc = np * step;	/* Number of (x,y[,t]) to copy */
		f = (double *) GMT_memory (VNULL, (size_t)(n_alloc+1), sizeof (double), GMT_program);
		f[0] = (double)np;	/* Number of points found */
		memcpy ((void *)&f[1], (void *)&c[first], n_alloc * sizeof (double));
		GMT_free ((void *)c);	/* Free the old flowline vector */
		*flow = f;		/* Return pointer to trimmed flowline */
	}
	else
		*flow = c;		/* Return the entire flowline as is */
	return (np);
}

GMT_LONG set_age (double *t_smt, float *age, GMT_LONG node, double upper_age, GMT_LONG truncate, GMT_LONG no_ages)
{
	/* Returns the age of this node based on either a given seafloor age grid
	 * or the upper age, truncated if necessary */
	 
	if (no_ages || GMT_is_fnan (age[node]))		/* Age is NaN, assign upper value */
		*t_smt = upper_age;
	else {	/* Assign given value */
		*t_smt = age[node];
		if (*t_smt > upper_age) {	/* Exceeding our limit */
			if (truncate)		/* Allowed to truncate to max age */
				*t_smt = upper_age;
			else {			/* Consider this an error or just skip */
				fprintf (stderr, "%s: Node %ld has age (%g) > oldest stage (%g) (skipped)\n", GMT_program, node, *t_smt, upper_age);
				return (FALSE);
			}
		}
	}
	return (TRUE);	/* We are returning a useful age */
}
