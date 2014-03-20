/*--------------------------------------------------------------------
 *	$Id: hotspotter.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *   Copyright (c) 1999-2014 by P. Wessel
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
 * HOTSPOTTER will (1) read ascii file(s) with records for each seamount
 * (2) read an ascii file with stage (Euler) rotations, and (3)
 * convolve the flowline of each seamount with its gravimetric shape, and
 * (4) build a cumulative volcano amplitude (CVA) grid.  The grid is
 * written out in GMT format and can be processed and plotted with GMT.
 * HOTSPOTTER is part of the SPOTTER supplemental GMT package and should
 * be installed accordingly, see README.spotter.
 *
 * Author:	Paul Wessel, SOEST, Univ. of Hawaii, Honolulu, HI, USA
 * Date:	22-JUN-1999
 * Version:	1.0
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
 *
 * Seamount data file(s) must have the following format:
 * 1. Any number of comment lines starting with # in first column
 * 2. Any number of blank lines (just carriage return, no spaces)
 * 3. ANy number of seamount records which each have the format:
 *    lon(deg)  lat(deg)  amplitude radius(km) age(Ma)
 * 4. The amplitude is in user units (m, mGal, km^3, whatever)
 * 5. Age represents the upper possible age for seamount, which
 *    is usually the age of the seafloor beneath it.  If the
 *    crustal age is not known, set it to NaN in the file.
 *    NaN values will be replaced with -N (see usage)
 *
 * Example: Wessel & Lyons [1997] Pacific seamounts (just a few records):
 * # Pacific seamounts > 100 Eotvos amplitude in Vertical Gravity Gradient
 * # From Wessel & Lyons [1997]
 * #LON		  FAA	 	 VGG	RADIUS	CRUST_AGE
 * 134.38333	0.9436415	120.5	 22.97	37.606796
 * 136.05	7.6325042	102.7	 18.67	NaN
 * 131.28333	1.1423035	129.0	 17.16	NaN
 * .....
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
 * -> Seamount data set:
 *
 * Wessel, P., and S. Lyons, 1997, Distribution of large Pacific seamounts
 *    from Geosat/ERS-1: Implications for the history of intraplate volcanism,
 *    J. Geophys. Res., 102,22,459-22,475.
 * Wessel, P., 1997, Sizes and ages of seamounts using remote sensing:
 *    Implications for intraplate volcanism, Science, 277, 802-805.
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

int main (int argc, char **argv)
{

	GMT_LONG n_smts;		/* Number of seamounts read */
	GMT_LONG n_stages;		/* Number of stage rotations (poles) */
	GMT_LONG n_chunk;		/* Number of path values returned by libspotter functions */
	GMT_LONG n_track;		/* Number of points along a single flowline */
	GMT_LONG node_x_width;		/* Number of x-nodes covered by the seamount in question (y-dependent) */
	GMT_LONG node_y_width;		/* Number of y-nodes covered by the seamount */
	GMT_LONG n_args;		/* Command line argument counter */
	GMT_LONG n_files = 0;		/* Number of input files on the command line */
	GMT_LONG node;			/* The current node index */
	GMT_LONG one_or_zero;		/* Will be 1 or 0 depending on grid format (-F) */
	GMT_LONG n_fields, n_expected_fields;
	GMT_LONG n_read = 0;		/* Number of records read */
	GMT_LONG n_nodes;		/* Number of nodes in the output CVA grid */
	GMT_LONG node_offset = 0;

	/* Misc. counters */

	GMT_LONG i, j, kx, ky, m, ii, jj, i0, j0, k0, dummy[4], fno;

	GMT_LONG	error = FALSE;		/* TRUE when arguments are wrong */
	GMT_LONG nofile = TRUE;		/* TRUE when data are to be read from stdin */
	GMT_LONG truncate_ages = FALSE;	/* Truncate all seamounts ages > upper_age */
	GMT_LONG finite = FALSE;		/* TRUE if stage pole file contains finite rotation poles instead */
	GMT_LONG normalize = FALSE;	/* TRUE if we want CVA in percents */
	GMT_LONG done;

	float *CVA = NULL;		/* The CVA surface we seek to calculate */

	double sampling_int_in_km;	/* Sampling interval along flowline (in km) */
	double sampling_factor = 0.5;	/* Sets how many points along flowline per node */
	double upper_age = 180.0;	/* Upper age assigned to seamounts on undated seafloor */
	double x_smt;			/* Seamount longitude (input degrees, stored as radians) */
	double y_smt;			/* Seamount latitude (input degrees, stored as radians) */
	double z_smt;			/* Seamount amplitude (in user units) */
	double r_smt;			/* Seamount radius (input km, stored as radians) */
	double t_smt;			/* Seamount upper age (up to age of seafloor) */
	double norm;			/* Normalization factor based on r_smt */
	double *xpos = NULL, *ypos = NULL;	/* Coordinates of the output grid (in radians) */
	double *c = NULL;		/* Array with one flowline */
	double dx, dy;			/* x,y distance from projected seamount center to nearest node */
	double *latfactor = NULL, *ilatfactor = NULL;	/* Arrays of latitudinal-dependent x-scales (cos(lat)-stuff) */
	double x_part, y_part;		/* Components of radius from projected seamount center to current node */
	double y_part2;			/* y_part squared */
	double r2;			/* Radius squared from projected seamount center to current node */
	double west, east, south, north;	/* Region in radians */
	double xinc_r, yinc_r;		/* Grid spacing in radians */
	double i_xinc_r, i_yinc_r;	/* Inverse grid spacing in 1/radians */
	double offset;			/* 0.0 or 0.5, depending onf -F */
	double half_dx, half_dy;	/* Half grid sizes or 0.0, depending on -F */
	double *in = NULL;			/* GMT read array */
	double south_c, north_c, yg;	/* Same but geocentric */

	char *CVA_file = CNULL;		/* Output name of CVA grid */
	char *euler_file = CNULL;	/* Input name of file with stage (Euler) rotations */
	char line[BUFSIZ];		/* Read buffer */
	char *processed_node = NULL;	/* Pointer to array with TRUE/FALSE values for each grid node */
	char *not_used = NULL;

	FILE *fp = NULL;		/* Input data file pointer */

	struct GRD_HEADER CVA_map;	/* Header structure for output CVA grid */

	struct EULER *p = NULL;		/* Array of structures with Euler stage rotations */

	/* ------------------- END OF DECLARATIONS ------------------------------------------------------------*/

	argc = (int)GMT_begin (argc, argv);			/* Initialize GMT environment */

	GMT_grd_init (&CVA_map, argc, argv, FALSE);	/* Initialize grid structure */

	GMT_io.in_col_type[0] = GMT_IS_LON;
	GMT_io.in_col_type[1] = GMT_IS_LAT;

	/* Check command line arguments */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'V':
				case 'H':
				case 'R':
				case ':':
				case '\0':
					error += GMT_parse_common_options (argv[i], &CVA_map.x_min, &CVA_map.x_max, &CVA_map.y_min, &CVA_map.y_max);
					break;

				/* Supplemental parameters */

				case 'b':
					error += GMT_parse_b_option (&argv[i][2]);
					break;

				case 'C':	/* Use finite rotation poles */
					finite = TRUE;
					break;

				case 'D':
					sampling_factor = atof (&argv[i][2]);
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
				case 'I':
					GMT_getinc (&argv[i][2], &CVA_map.x_inc, &CVA_map.y_inc);
					break;
				case 'N':
					upper_age = atof (&argv[i][2]);
					break;
				case 'S':
					normalize = TRUE;
					break;
				case 'T':
					truncate_ages = TRUE;
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
		fprintf (stderr, "%s %s - Create CVA image from seamount locations\n\n", GMT_program, SPOTTER_VERSION);
		fprintf (stderr, "usage: %s [<files>] -E<stage_file> -G<CVAgrid> %s\n", GMT_program, GMT_Id_OPT);
		fprintf (stderr, "\t%s [-C] [-D<factor>] [-F] [%s] [-N<upper_age>]\n", GMT_Rgeo_OPT, GMT_H_OPT);
		fprintf (stderr, "\t[-S] [-T] [-V] [%s] [%s]\n\n", GMT_t_OPT, GMT_bi_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<files> is one or more seamount (x,y,z,r,t) files\n");
		fprintf (stderr, "\t-E specifies the rotations to be used (see man page for format)\n\n");
		fprintf (stderr, "\t-G Specify file name for output CVA grid.\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I specifies grid interval(s); Append m [or c] to <dx> and/or <dy> for minutes [or seconds].\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-C The file given with -E contains finite rotation poles [Default is stage poles]\n");
		fprintf (stderr, "\t-D Scale affecting distance between points along flowline [0.5]\n");
		fprintf (stderr, "\t-F Force pixel-registration [Default is gridline registration]\n");
		fprintf (stderr, "\t-N sets upper age in m.y. for seamounts whose plate age is NaN [180]\n");
		fprintf (stderr, "\t-S Normalize CVA grid to percentages of the CVA maximum\n");
		fprintf (stderr, "\t-T truncate all ages to max age in stage pole model [Default extrapolates]\n");
		GMT_explain_option ('V');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "	   Default is 5 input columns\n");
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
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 5;
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 5) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 5 columns\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	/* ------------------- END OF PROCESSING COMMAND LINE ARGUMENTS  --------------------------------------*/

	GMT_lat_swap_init ();	/* Initialize auxiliary latitude machinery */

	/* Load in the Euler stage poles */

	n_stages = spotter_init (euler_file, &p, TRUE, finite, FALSE, &upper_age, gmtdefs.verbose);

	/* Initialize the CVA grid and structure */

	xinc_r = CVA_map.x_inc * D2R;
	yinc_r = CVA_map.y_inc * D2R;

	if (CVA_map.node_offset) {	/* Pixel grid */
		one_or_zero = 0;
		offset = 0.0;
		half_dx = 0.5 * xinc_r;
		half_dy = 0.5 * yinc_r;
	}
	else {
		one_or_zero = 1;
		offset = 0.5;
		half_dx = half_dy = 0.0;
	}

	CVA_map.nx = (int)(irint ((CVA_map.x_max - CVA_map.x_min) / CVA_map.x_inc) + one_or_zero);
	CVA_map.ny = (int)(irint ((CVA_map.y_max - CVA_map.y_min) / CVA_map.y_inc) + one_or_zero);
	n_nodes = GMT_get_nm (CVA_map.nx, CVA_map.ny);
	CVA = (float *) GMT_memory (VNULL, (size_t)n_nodes, sizeof(float), GMT_program);

	/* Assign grid-region variables in radians to avoid conversions inside convolution loop */

	west   = CVA_map.x_min * D2R;
	east   = CVA_map.x_max * D2R;
	south  = CVA_map.y_min * D2R;
	north  = CVA_map.y_max * D2R;
	south_c  = GMT_lat_swap (CVA_map.y_min, GMT_LATSWAP_G2O) * D2R;	/* Convert to geocentric */ 
	north_c  = GMT_lat_swap (CVA_map.y_max, GMT_LATSWAP_G2O) * D2R;	/* Convert to geocentric */ 
	i_xinc_r = 1.0 / xinc_r;
	i_yinc_r = 1.0 / yinc_r;

	/* Precalculate coordinates xpos[], ypos[] and scale factors(lat) on the grid */

	xpos = (double *) GMT_memory (VNULL, (size_t)CVA_map.nx, sizeof(double), GMT_program);
	ypos = (double *) GMT_memory (VNULL, (size_t)CVA_map.ny, sizeof(double), GMT_program);
	latfactor  = (double *) GMT_memory (VNULL, (size_t)CVA_map.ny, sizeof(double), GMT_program);
	ilatfactor = (double *) GMT_memory (VNULL, (size_t)CVA_map.ny, sizeof(double), GMT_program);

	for (i = 0; i < CVA_map.nx; i++) xpos[i] = west + i * xinc_r + half_dx;
	for (j = 0; j < CVA_map.ny; j++) {
		ypos[j] = north - j * yinc_r - half_dy;
		latfactor[j] = xinc_r * cos (ypos[j]);
		ilatfactor[j] = 1.0 / latfactor[j];
	}

	/* Allocate T/F array */

	processed_node = (char *) GMT_memory (VNULL, n_nodes, sizeof (char), GMT_program);

	/* Set flowline sampling interval to 1/2 of the shortest distance between x-nodes */

	sampling_int_in_km = sampling_factor * xinc_r * EQ_RAD * ((fabs (north) > fabs (south)) ? cos (north) : cos (south));
	if (gmtdefs.verbose) fprintf (stderr, "%s: Flowline sampling interval = %.3f km\n", GMT_program, sampling_int_in_km);

	if (truncate_ages && gmtdefs.verbose) fprintf (stderr, "%s: Seamount ages truncated to %g\n", GMT_program, upper_age);

	/* Start to read input data */

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;
	n_args = (argc > 1) ? argc : 2;

	done = FALSE;

	n_smts = 0;

	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 5;

	for (fno = 1; !done && fno < n_args; fno++) {	/* Loop over all input files */

		if (!nofile && argv[fno][0] == '-') continue;	/* Skip argument */
		if (nofile) {
			fp = stdin;
			done = TRUE;
		}
		else if ((fp = GMT_fopen (argv[fno], "r")) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[fno]);
			continue;
		}

		if (gmtdefs.verbose) (nofile) ? fprintf (stderr, "%s: Reading stdin\n", GMT_program) : fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[fno]);

		if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) not_used = GMT_fgets (line, BUFSIZ, fp);

		while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */
			n_read++;
			while ((GMT_io.status & GMT_IO_SEGMENT_HEADER) && !(GMT_io.status & GMT_IO_EOF)) {
				GMT_write_segmentheader (GMT_stdout, n_expected_fields);
				n_fields = GMT_input (fp, &n_expected_fields, &in);
				n_read++;
			}
			if (GMT_io.status & GMT_IO_EOF) continue;

			if (GMT_io.status & GMT_IO_MISMATCH) {
				fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
				continue;
			}

			in[GMT_Y] = GMT_lat_swap (in[GMT_Y], GMT_LATSWAP_G2O);	/* Convert to geocentric */

			/* STEP 1: Read information about a single seamount from input record */

			x_smt = in[0] * D2R;	/* Seamount positions in RADIANS */
			y_smt = in[1] * D2R;
			z_smt = in[2];
			r_smt = in[3];
			if (GMT_is_dnan (in[4])) {	/* Age is NaN, assign value */
				t_smt = upper_age;
			}
			else {			/* Assign given value, truncate if necessary */
				t_smt = in[4];
				if (t_smt > upper_age) {
					if (truncate_ages) {
						t_smt = upper_age;
					}
					else {
						fprintf (stderr, "%s: Seamounts near line %ld has age (%g) > oldest stage (%g) (skipped)\n", GMT_program, n_read, t_smt, upper_age);
						continue;
					}
				}
			}

			/* Do some normalizations here to save processing inside convolution later */

			r_smt /= EQ_RAD;				/* Converts radius in km to radians */
			norm = -4.5 / (r_smt * r_smt);			/* Gaussian normalization */
			node_y_width = (GMT_LONG)ceil (i_yinc_r * r_smt);	/* y-node coverage */

			/* STEP 2: Calculate this seamount's flowline */

			n_chunk = spotter_forthtrack (&x_smt, &y_smt, &t_smt, (GMT_LONG)1, p, n_stages, sampling_int_in_km, 0.0, FALSE, NULL, &c);

			/* STEP 3: Convolve this flowline with seamount shape and add to CVA grid */

			/* Our convolution is approximate:  We sample the flowline frequently and use
			 * one of the points on the flowline that are closest to the node.  Ideally we
			 * want the nearest distance from each node to the flowline. Later versions may
			 * improve on this situation */

			n_track = irint (c[0]);				/* Number of point pairs making up this flowline */

			memset ((void *)processed_node, 0, n_nodes);	/* Fresh start for this flowline convolution */

			for (m = 0, kx = 1; m < n_track; m++, kx += 2) {		/* For each point along flowline */

				ky = kx + 1;	/* Index for the y-coordinate */

				/* First throw out points outside specified grid region */

				if (c[ky] < south_c || c[ky] > north_c) continue;	/* Latitude outside region */

				if (c[kx] > east) c[kx] -= TWO_PI;
				while (c[kx] < west) c[kx] += TWO_PI;
				if (c[kx] > east) continue;			/* Longitude outside region */

				/* OK, this point is within our region, get node index */

				yg = GMT_lat_swap (c[ky], GMT_LATSWAP_O2G);	/* Convert back to geodetic */
				i = (GMT_LONG) floor (((c[kx] - west)  * i_xinc_r) + offset);
				j = (GMT_LONG) floor (((north - yg) * i_yinc_r) + offset);
				node = j * CVA_map.nx + i;

				if (!processed_node[node]) {	/* Have not added to the CVA at this node yet */

					/* Shape is z_smt * exp (r^2 * norm) */

					node_x_width = (GMT_LONG) ceil (r_smt * ilatfactor[j]);
					dx = c[kx] - xpos[i];
					dy = c[ky] - ypos[j];

					/* Loop over a square that circumscribes this seamounts basal outline */

					for (jj = -node_y_width, j0 = j - node_y_width; jj <= node_y_width; jj++, j0++) {

						if (j0 < 0 || j0 >= CVA_map.ny) continue;	/* Outside grid */

						y_part = jj * yinc_r - dy;
						y_part2 = y_part * y_part;
						k0 = j0 * CVA_map.nx;

						for (ii = -node_x_width, i0 = i - node_x_width; ii <= node_x_width; ii++, i0++) {

							if (i0 < 0 || i0 >= CVA_map.nx) continue;	/* Outside grid */

							x_part = ii * latfactor[j] - dx;
							r2 = (x_part * x_part + y_part2) * norm;
							CVA[k0+i0] += (float)(z_smt * exp (r2));
						}
					}
					processed_node[node] = 1;	/* Now we have visited this node */
				}
			}

			GMT_free ((void *)c);	/* Free the flowline vector */

			n_smts++;	/* Go to next seamount */

			if (gmtdefs.verbose && !(n_smts%100)) fprintf (stderr, "%s: Processed %5ld seamounts\r", GMT_program, n_smts);
		}

		if (fp != stdin) GMT_fclose (fp);
	}

	/* OK, Done processing, time to write out */

	if (gmtdefs.verbose) fprintf (stderr, "%s: Processed %5ld seamounts\n", GMT_program, n_smts);

	if (normalize) {	/* Convert CVA values to percent of CVA maximum */
		size_t k;
		double scale;
		
		if (gmtdefs.verbose) fprintf (stderr, "%s: Normalize CVS grid to percentages of max CVA\n", GMT_program);
		CVA_map.z_min = +DBL_MAX;
		CVA_map.z_max = -DBL_MAX;
		for (k = 0; k < (size_t)n_nodes; k++) {
			if (GMT_is_fnan (CVA[k])) continue;
			if (CVA[k] < CVA_map.z_min) CVA_map.z_min = CVA[k];
			if (CVA[k] > CVA_map.z_max) CVA_map.z_max = CVA[k];
		}
		if (gmtdefs.verbose) fprintf (stderr, "%s: CVA min/max: %g %g -> ", GMT_program, CVA_map.z_min, CVA_map.z_max);
		scale = 100.0 / CVA_map.z_max;
		for (k = 0; k < (size_t)n_nodes; k++) CVA[k] *= (float)scale;
		CVA_map.z_min *= scale;
		CVA_map.z_max *= scale;
		if (gmtdefs.verbose) fprintf (stderr, "%g %g\n", CVA_map.z_min, CVA_map.z_max);
	}
	if (gmtdefs.verbose) fprintf (stderr, "%s: Write CVA grid %s\n", GMT_program, CVA_file);

	dummy[0] = dummy[1] = dummy[2] = dummy[3] = 0;	/* No grid boundary padding today */

	GMT_err_fail (GMT_write_grd (CVA_file, &CVA_map, CVA, 0.0, 0.0, 0.0, 0.0, dummy, FALSE), CVA_file);

	/* Clean up memory */

	GMT_free ((void *)CVA);
	GMT_free ((void *)processed_node);
	GMT_free ((void *)latfactor);
	GMT_free ((void *)ilatfactor);
	GMT_free ((void *)xpos);
	GMT_free ((void *)ypos);
	GMT_free ((void *)p);

	if (gmtdefs.verbose) fprintf (stderr, "%s: Done\n", GMT_program);

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}
