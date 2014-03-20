/*--------------------------------------------------------------------
 *	$Id: mgd77magref.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *    Copyright (c) 2009-2014 by J. Luis and P. Wessel
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * mgd77magref produces output derived from input locations and time and
 * the CM4 or IGRF magnetic field models.
 *
 * Author:	Joaquim Luis and Paul Wessel
 * Date:	1-MAY-2009
 * Version:	1.0
 *
 */
 
#include "mgd77.h"

int main (int argc, char **argv) {

	GMT_LONG i, j, s, pos_slash = 0, nval = 0, nfval = 0, fno, n_files = 0, n_args, t_col = 3, n_out = 0, n_in;
	GMT_LONG lval = 0, lfval = 0, n_field_components, this_i, pos, error = 0;
	GMT_LONG nofile = TRUE, copy_input = FALSE, time_in_years = FALSE, fixed_alt = FALSE, fixed_time = FALSE;
	GMT_LONG do_CM4 = TRUE, joint_IGRF_CM4 = FALSE, do_IGRF = FALSE, do_CM4core = FALSE, cm4_igrf_T = FALSE, done;
	size_t n_alloc = 0, need = 0;
	double	the_altitude, the_time, *time_array = NULL, *alt_array = NULL, *time_years = NULL, IGRF[7], out[GMT_MAX_COLUMNS];
	double	*igrf_xyz = NULL;	/* Temporary storage for the joint_IGRF_CM4 case */
	char p[GMT_LONG_TEXT], txt[GMT_LONG_TEXT], tfixed[GMT_LONG_TEXT];
	FILE *fp = NULL;
	struct MGD77_CONTROL M;
	struct	MGD77_CM4 *Ctrl = NULL;
	struct GMT_TABLE *T = NULL;
	void *New_CM4_Ctrl (), Free_CM4_Ctrl (struct MGD77_CM4 *C);

	Ctrl = (struct MGD77_CM4 *)New_CM4_Ctrl ();	/* Allocate and initialize a new control structure */

	argc = (int) GMT_begin (argc, argv);		/* Initialize GMT Machinery */
	
	MGD77_Init (&M);			/* Initialize MGD77 Machinery */
	MGD77_CM4_init (&M, Ctrl);		/* Presets path using strdup */

	Ctrl->D.dst = (double *) calloc((size_t)(1), sizeof(double));	/* We need at least a size of one in case a value is given in input */
	GMT_io.in_col_type[t_col] = GMT_io.out_col_type[t_col] = GMT_IS_ABSTIME;	/* By default, time is in 4th input column */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'H':
				case 'V':
				case 'b':
				case 'f':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'A':
					strcpy (txt, &argv[i][2]);
					pos = 0;
					while ((GMT_strtok (txt, "+", &pos, p))) {
						switch (p[0]) {
							case 'a':
								fixed_alt = TRUE;
								the_altitude = atof (&p[1]);
								t_col = 2;	/* Since we are missing the altitude column */
								break;
							case 't':
								fixed_time = TRUE;
								the_time = atof (&p[1]);
								strcpy (tfixed, &p[1]);
								GMT_io.out_col_type[3] = GMT_IS_FLOAT;
								break;
							case 'y':
								time_in_years = TRUE;
								GMT_io.in_col_type[2] = GMT_io.out_col_type[2] = GMT_io.in_col_type[3] = GMT_io.out_col_type[3] = GMT_IS_FLOAT;
								break;
							default:
								break;
						}
					}
					if (fixed_time) {
						if (time_in_years)
							the_time = atof (tfixed);
						else
							GMT_scanf_arg (tfixed, GMT_IS_ABSTIME, &the_time);
					}
					break;
				case 'C':	/* Alternate CM4 coefficient file */
					free ((void *)Ctrl->M.path);
					Ctrl->M.path = strdup (&argv[i][2]);
					break;
				case 'D':
					this_i = 2;
					if (argv[i][2] == '-') this_i++;
					if ((argv[i][this_i] > 47) && (argv[i][this_i] < 58)) {	/* arg is numeric -> Dst Index */
						Ctrl->D.dst[0] = atof(&argv[i][2]);
						Ctrl->D.index = FALSE; 
					}
					else {
						free ((void *)Ctrl->D.path);
						Ctrl->D.path = strdup (&argv[i][2]);
						Ctrl->D.load = TRUE;
					}
					break;
				case 'E':
					if ((argv[i][2] > 47) && (argv[i][2] < 58)) {	/* arg is numeric -> Dst Index */
						Ctrl->I.F107 = atof(&argv[i][2]);
						Ctrl->I.index = FALSE;
					}
					else {
						free ((void *)Ctrl->I.path);
						Ctrl->I.path = strdup (&argv[i][2]);
						Ctrl->I.load = TRUE;
					}
					break;
				case 'F':
					Ctrl->F.active = TRUE;

					pos_slash = 0;
					for (j = 2; argv[i][j]; j++) {
						if (argv[i][j] == '/') {
							pos_slash = j;
							break;
						}
						switch (argv[i][j]) {
							case 'r':		/* Echo input record */
								copy_input = TRUE; 
								break;
							case 't':		/* Total field is requested */
								Ctrl->F.field_components[nval++] = 0; 
								break;
							case 'h':		/* Horizontal field is requested */
								Ctrl->F.field_components[nval++] = 1; 
								break;
							case 'x':		/* X component is requested */
								Ctrl->F.field_components[nval++] = 2; 
								break;
							case 'y':		/* Y component is requested */
								Ctrl->F.field_components[nval++] = 3; 
								break;
							case 'z':		/* Z component is requested */
								Ctrl->F.field_components[nval++] = 4; 
								break;
							case 'd':		/* Declination is requested */
								Ctrl->F.field_components[nval++] = 5; 
								break;
							case 'i':		/* Inclination is requested */
								Ctrl->F.field_components[nval++] = 6; 
								break;
						}
					}
					Ctrl->F.n_field_components = (int)nval; 

					if (pos_slash) {
						for (j = pos_slash; argv[i][j]; j++) {
							switch (argv[i][j]) {
								case '0':		/* IGRF */
									do_IGRF = TRUE;
									do_CM4  = FALSE;
									joint_IGRF_CM4 = FALSE;
									break;
								case '1':		/* Main field 1 */
									Ctrl->F.field_sources[nfval++] = 0;
									do_CM4core = TRUE; 
									break;
								case '2':		/* Main field 2 */
									Ctrl->F.field_sources[nfval++] = 1; 
									break;
								case '3':		/* Primary Magnetospheric field */
									Ctrl->F.field_sources[nfval++] = 2; 
									break;
								case '4':		/* Induced Magnetospheric field */
									Ctrl->F.field_sources[nfval++] = 3; 
									break;
								case '5':		/* Primary ionospheric field */
									Ctrl->F.field_sources[nfval++] = 4; 
									break;
								case '6':		/* Induced ionospheric field */
									Ctrl->F.field_sources[nfval++] = 5; 
									break;
								case '7':		/* Toroidal field */
									Ctrl->F.field_sources[nfval++] = 6; 
									break;
								case '9':		/* Main field is computed with the IGRF */
									do_IGRF = FALSE;/* No contradiction, the way will be through joint_IGRF_CM4 */
									do_CM4  = TRUE;	/* Well maybe, if some other source is selected also */
									joint_IGRF_CM4 = TRUE;
									break;
							}
						}
						Ctrl->F.n_field_sources = (int)nfval;
					}
					break;
				case 'G':
					Ctrl->G.geodetic = FALSE;
					break;
				case 'L':
					Ctrl->L.curr = TRUE;

					pos_slash = 0;
					for (j = 2; argv[i][j]; j++) {
						if (argv[i][j] == '/') {
							pos_slash = j;
							break;
						}
						switch (argv[i][j]) {
							case 'r':		/* Echo input record */
								copy_input = TRUE; 
								break;
							case 't':		/* Total field is requested */
								Ctrl->L.curr_components[lval++] = 0; 
								break;
							case 'x':		/* X component is requested */
								Ctrl->L.curr_components[lval++] = 1; 
								break;
							case 'y':		/* Y component is requested */
								Ctrl->L.curr_components[lval++] = 2; 
								break;
							case 'z':		/* Z component is requested */
								Ctrl->L.curr_components[lval++] = 3; 
								break;
						}
					}
					Ctrl->L.n_curr_components = (int)lval; 

					if (pos_slash) {
						for (j = pos_slash; argv[i][j]; j++) {
							switch (argv[i][j]) {
								case '1':		/* Induced Magnetospheric field */
									Ctrl->L.curr_sources[lfval++] = 0; 
									break;
								case '2':		/* Primary ionospheric field */
									Ctrl->L.curr_sources[lfval++] = 1; 
									break;
								case '3':		/* Induced ionospheric field */
									Ctrl->L.curr_sources[lfval++] = 2; 
									break;
								case '4':		/* Poloidal field */
									Ctrl->L.curr_sources[lfval++] = 3; 
									break;
							}
						}
						Ctrl->L.n_curr_sources = (int)lfval; 
					}
					break;
				case 'S':
					if (argv[i][2] == 'c') {
						j = sscanf (&argv[i][3], "%d/%d", &Ctrl->S.nlmf[0], &Ctrl->S.nhmf[0]);
						if (j != 2) {
							fprintf (stderr, "%s: ERROR: -Sc option usage is -Sc<low/high>\n", GMT_program);
							error++;
						}
					}
					if (argv[i][2] == 'l') {
						j = sscanf (&argv[i][3], "%d/%d", &Ctrl->S.nlmf[1], &Ctrl->S.nhmf[1]);
						if (j != 2) {
							fprintf (stderr, "%s: ERROR: -Sl option usage is -Sl<low/high>\n", GMT_program);
							error++;
						}
					}
					break;
				default:
					fprintf (stderr, "%s: ERROR: Unrecognized option %s\n", GMT_program, argv[i]);
					error++;
					break;
			}
		}
		else
			n_files++;
	}

	if (error || GMT_give_synopsis_and_exit) {	/* Display usage */
		fprintf (stderr, "mgd77magref %s - Evaluating the IGRF or CM4 magnetic field models\n\n", GMT_VERSION);
		fprintf (stderr, "usage: mgd77magref [<inputfile>] [-A+y+a<alt>+t<date>] [-C<cm4file>] [-D<dstfile>] [-E<f107file>]\n");
		fprintf (stderr, "\t[-F<rthxyzdi[/[0|9]1234567]>] [-G] [%s] [-L<rtxyz[/1234]>] [-Sc|l<low/high>]\n", GMT_H_OPT);
		fprintf (stderr, "\t[-V] [%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_b_OPT, GMT_m_OPT);
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t<inputfile> contains records that must contain lon, lat, alt, time[, other cols]\n");
		fprintf (stderr, "\t   longitude and latitude is the geocentric position on the ellipsoid [but see -G].\n");
		fprintf (stderr, "\t   alt is the altitude in km positive above the ellipsoid.\n");
		fprintf (stderr, "\t   time is the time of data aquisition, in <date>T<clock> format (but see -A+y).\n");
		fprintf (stderr, "\t   We read <stdin> if no input file is given.\n");
		fprintf (stderr, "\t-A Adjusts how the input records are interpreted. Append\n");
		fprintf (stderr, "\t   +a<alt> to indicate a constant altitude [Default is 3rd column].\n");
		fprintf (stderr, "\t   +t<time> to indicate a constant time [Default is 4th column].\n");
		fprintf (stderr, "\t   +y to indicate times are given in decimal years [Default is ISO <date>T<clock> format].\n");
		fprintf (stderr, "\t-C Selects an alternate file with coefficients for the CM4 model [%s/umdl.CM4].\n", GMT_SHAREDIR);
		fprintf (stderr, "\t-D Selects an alternate file with hourly means of the Dst index for CM4 [%s/Dst_all.wdc],\n", GMT_SHAREDIR);
		fprintf (stderr, "\t   OR a single Dst index to apply for all records.\n");
		fprintf (stderr, "\t-E Selects an alternate file with monthly means of absolute F10.7 solar radio flux for CM4 [%s/F107_mon.plt],\n", GMT_SHAREDIR);
		fprintf (stderr, "\t   OR a single solar radio flux to apply for all records.\n");
		fprintf (stderr, "\t-F Dataflags is a string made up of 1 or more of these characters:\n");
		fprintf (stderr, "\t	 r means output all input columns before adding the items below (all in nTesla).\n");
		fprintf (stderr, "\t	 t means list total field.\n");
		fprintf (stderr, "\t	 h means list horizontal field.\n");
		fprintf (stderr, "\t	 x means list X component.\n");
		fprintf (stderr, "\t	 y means list Y component.\n");
		fprintf (stderr, "\t	 z means list Z component.\n");
		fprintf (stderr, "\t	 d means list declination.\n");
		fprintf (stderr, "\t	 i means list inclination.\n");
		fprintf (stderr, "\t   Append a number to indicate the requested field contribution(s)\n");
		fprintf (stderr, "\t	 0 means Core field from IGRF only (no CM4 evalution).\n");
		fprintf (stderr, "\t	 1 means Core field.\n");
		fprintf (stderr, "\t	 2 means Lithospheric field.\n");
		fprintf (stderr, "\t	 3 Primary Magnetospheric field.\n");
		fprintf (stderr, "\t	 4 Induced Magnetospheric field.\n");
		fprintf (stderr, "\t	 5 Primary ionospheric field.\n");
		fprintf (stderr, "\t	 6 Induced ionospheric field.\n");
		fprintf (stderr, "\t	 7 Toroidal field.\n");
		fprintf (stderr, "\t	 9 means Core field from IGRF and other contributions from CM4. DO NOT USE BOTH 1 AND 9.\n");
		fprintf (stderr, "\t   Append several numbers to add up the different contributions. For example,\n");
		fprintf (stderr, "\t     -Ft/12 computes the total field due to CM4 Core and Lithospheric sources.\n");
		fprintf (stderr, "\t     Two special cases are allowed which mix which Core field from IGRF and other sources from CM4.\n");
		fprintf (stderr, "\t     -Ft/934 computes Core field due to IGRF plus terms 3 and 4 from CM4.\n");
		fprintf (stderr, "\t     -Fxyz/934 the same as above but output the field components.\n");
		fprintf (stderr, "\t	 The data is written out in the order specified in <dataflags>\n");
		fprintf (stderr, "\t	 [Default is -Frthxyzdi/1]\n");
		fprintf (stderr, "\t-G Specifies that coordinates are geocentric [geodetic].\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-L Computes J field vectors from certain external sources.\n");
		fprintf (stderr, "\t   Dataflags is a string made up of 1 or more of these characters:\n");
		fprintf (stderr, "\t	 r means output all input columns before adding the items below (all in Ampers/m).\n");
		fprintf (stderr, "\t	 t means list magnitude field.\n");
		fprintf (stderr, "\t	 x means list X component.\n");
		fprintf (stderr, "\t	 y means list Y component.\n");
		fprintf (stderr, "\t	 z means list Z or current function Psi.\n");
		fprintf (stderr, "\t   Append a number to indicate the requested J contribution(s)\n");
		fprintf (stderr, "\t	 1 means Induced Magnetospheric field.\n");
		fprintf (stderr, "\t	 2 means Primary ionospheric field.\n");
		fprintf (stderr, "\t	 3 means Induced ionospheric field.\n");
		fprintf (stderr, "\t	 4 means Poloidal field.\n");
		fprintf (stderr, "\t-S limits the CM4 contributions from core and lithosphere to certain harmonic degree bands.\n");
		fprintf (stderr, "\t   Append c(ore) or l(ithosphere) and the low and high degrees to use [-Sc1/13 -Sl14/65].\n");
		GMT_explain_option ('V');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 4 input columns (unless -A is used).  Note for binary input, absolute time must\n");
		fprintf (stderr, "\t   be in the unix time-system (unless -A+y is used).\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	n_out = 4 - (int)(fixed_alt + fixed_time);	/* Minimum input columns (could be more) */
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = n_out;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < n_out) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least %ld columns\n", GMT_program, n_out);
		error++;
	}
	n_in = n_out;

	if (Ctrl->F.active && Ctrl->L.curr) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR. You cannot select both -F and -L options.\n", GMT_program);
		error++;
	}
						
	if ((do_CM4core && do_IGRF) || (do_CM4core && joint_IGRF_CM4)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR: You cannot select both CM4 core (1) and IGRF as they are both core fields.\n", GMT_program);
			error++;
	}

	if ( nfval && do_IGRF )
		fprintf (stderr, "%s: WARNING. Source fields other than IGRF will be ignored. It's in the manual\n", GMT_program);

	/* ------------- Test, and take measures, if mix mode IGRF/CM4 is to be used -------------------------- */
	if (joint_IGRF_CM4) {
		if (nfval == 0) {	/* It means we had a -F.../9 which is exactly iqual to -F.../0 */
			do_IGRF = TRUE;
			do_CM4  = FALSE;
			joint_IGRF_CM4 = FALSE;
		}
		else {
			cm4_igrf_T = FALSE;
			if ( (nval == 1) && (Ctrl->F.field_components[0] == 0) ) {
				for (i = 0; i < 3; i++) Ctrl->F.field_components[i] = (int)(i+2);	/* Force the x,y,z request */
				cm4_igrf_T = TRUE;
			}
			else if ( !((nval == 3) && 
				(Ctrl->F.field_components[0] == 2) && (Ctrl->F.field_components[1] == 3) && (Ctrl->F.field_components[2] == 4)) ) {
				fprintf (stderr, "%s: GMT ERROR. In mix CM4/IGRF mode -F option can oly be -Ft[r]/... or -Fxyz[r]/...\n", GMT_program);
				error++;
			}

			nval = 3;
			Ctrl->F.n_field_components = (int)nval; 
		}
	}
	/* ----------------------------------------------------------------------------------------------------- */

	if (error) exit (EXIT_FAILURE);

	if (!Ctrl->F.active && !Ctrl->L.curr) Ctrl->F.active = TRUE;

	/* Sort the order in which the parameters appear */
	if (Ctrl->F.active) {
		if (nval == 0) {		/* No components selected, default used */
			copy_input = TRUE;
			Ctrl->F.n_field_components = 7;
			for (i = 0; i < 7; i++) Ctrl->F.field_components[i] = (int)i;
		}
		if (nfval == 0 && !joint_IGRF_CM4) {		/* No sources selected, retain only the main field */
			Ctrl->F.field_sources[0] = 0; 
			Ctrl->F.n_field_sources = 1;
		}
		n_field_components = Ctrl->F.n_field_components;
	}
	else {
		if (lval == 0) {		/* Nothing selected, default used */
			copy_input = TRUE;
			Ctrl->L.n_curr_components = 1;
			Ctrl->L.curr_components[0] = 0;
		}
		if (lfval == 0) {		/* Nothing selected, retain only the induced magnetospheric field */
			Ctrl->L.curr_sources[0] = 0; 
			Ctrl->L.n_curr_sources = 1;
		}
		n_field_components = Ctrl->L.n_curr_components;
	}

	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;
	n_args = (argc > 1) ? argc : 2;

	if (fixed_alt) {	/* A single altitude should apply to all records; set array to point to this single altitude */
		alt_array = &the_altitude;
		Ctrl->DATA.n_altitudes = 1;
	}
	if (fixed_time) {	/* A single time should apply to all records; set array to point to this single time */
		if (!time_in_years) {
			if (M.adjust_time) the_time = MGD77_time2utime (&M, the_time);	/* Convert to Unix time if need be */
			the_time = MGD77_time_to_fyear (&M, the_time);			/* Get decimal year */
		}
		time_array = &the_time;
		Ctrl->DATA.n_times = 1;
	}
	else	/* Make sure input time columns are encoded/decoded properly since here we know t_col is set. */
		GMT_io.in_col_type[t_col] = GMT_io.out_col_type[t_col] = (time_in_years) ? GMT_IS_FLOAT : GMT_IS_ABSTIME;

	GMT_io.in_col_type[t_col+1] = GMT_io.out_col_type[t_col+1] = GMT_IS_FLOAT;		/* Override any previous t_col = 3 settings */
	if (!copy_input) GMT_io.out_col_type[2] = GMT_io.out_col_type[3] = GMT_IS_FLOAT;	/* No time on output */
	
	/* Process all input files */
	
	for (fno = 1, done = FALSE; !done && fno < n_args; fno++) {	/* Loop over all input files */
		if (!nofile && argv[fno][0] == '-') continue;
		if (nofile) {
			if (gmtdefs.verbose) fprintf (stderr, "%s: Reading from standard input\n", GMT_program);
			fp = GMT_stdin;
			done = TRUE;
#ifdef SET_IO_MODE
			GMT_setmode (GMT_IN);
#endif
		}
		else if ((fp = GMT_fopen (argv[fno], GMT_io.r_mode)) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[fno]);
			continue;
		}
		if (!nofile && gmtdefs.verbose) {
			fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[fno]);
		}
		
		/* Read in an entire file */
		GMT_import_table ((void *)fp, GMT_IS_STREAM, &T, 0.0, FALSE, FALSE, TRUE);
		
		n_out = n_field_components + ((copy_input) ? T->segment[0]->n_columns : 0);
		if (GMT_io.binary[GMT_OUT] && GMT_io.ncol[GMT_OUT] > 0 && n_out > GMT_io.ncol[GMT_OUT]) {
			fprintf (stderr, "%s: Binary output must have at least %ld columns (your -bo option only set %ld)\n", 
				GMT_program, n_out, GMT_io.ncol[GMT_OUT]);
			exit (EXIT_FAILURE);
		}

		if (T->n_columns < n_in) {
			fprintf (stderr, "%s: File %ld has %ld columns, but from the used options we expect %ld\n",
				GMT_program, fno, T->n_columns, n_in);
			continue;
		}

		for (s = 0; s < T->n_segments; s++) {	/* Process each file segment separately */
			if (GMT_io.multi_segments[GMT_OUT]) {
				if (T->segment[s]->header)
					sprintf (GMT_io.segment_header, "%s", T->segment[s]->header);
				else
					GMT_io.segment_header[0] = '\0';
				GMT_write_segmentheader (GMT_stdout, n_out);
			}
			need = T->segment[s]->n_rows;	/* Size of output array needed in MGD77_cm4field */
			if (need > n_alloc) {	/* Need to reallocate */
				n_alloc = need;
				Ctrl->DATA.out_field = (double *) GMT_memory ((void *)Ctrl->DATA.out_field, n_alloc * n_field_components, sizeof(double), GMT_program);
				if (!(time_in_years || fixed_time)) 
					time_years = (double *) GMT_memory ((void *)time_years, n_alloc, sizeof(double), GMT_program);

				if (joint_IGRF_CM4)
					igrf_xyz = (double *) GMT_memory ((void *)igrf_xyz, n_alloc * 3, sizeof(double), GMT_program);
			}

			if (!fixed_alt) {	/* Assign the alt_array to the provided altitude array */
				alt_array = T->segment[s]->coord[GMT_Z];
				Ctrl->DATA.n_altitudes = (int)T->segment[s]->n_rows;
			}

			if (!fixed_time) {	/* Assign the time_array to the provided time array */
				Ctrl->DATA.n_times = (int)T->segment[s]->n_rows;
				if (time_in_years)
					time_array = T->segment[s]->coord[t_col];
				else {	/* Must convert internal GMT time to decimal years first */
					for (i = 0; i < T->segment[s]->n_rows; i++) 
						time_years[i] = MGD77_time_to_fyear (&M, T->segment[s]->coord[t_col][i]);
					time_array = time_years;
				}
			}

			if (!(do_IGRF || joint_IGRF_CM4 ) && !s && time_array[0] > 2002.7) {	/* Only atmospheric terms may be reliable */
				fprintf (stderr, "%s: WARNING. Time is outside the CM4 strict validity domain [1960.0-2002.7].\n", GMT_program);
				fprintf (stderr, "\tThe secular variation estimation will be unreliable. In this case\n"
							"\tyou really should use the IGRF to estimate the core contribution.\n");
			}

			Ctrl->DATA.n_pts = (int)T->segment[s]->n_rows;
			if (do_IGRF || joint_IGRF_CM4) {
				int type;
				type = (Ctrl->G.geodetic) ? 1 : 2;
				for (i = 0; i < T->segment[s]->n_rows; i++) {
					the_altitude = (fixed_alt) ? alt_array[0] : alt_array[i];
					the_time = (fixed_time) ? time_array[0] : time_array[i];
					if (type == 2) the_altitude += 6371.2;
					MGD77_igrf10syn (0, the_time, type, the_altitude, T->segment[s]->coord[GMT_X][i], 
							T->segment[s]->coord[GMT_Y][i], IGRF);
					if (!joint_IGRF_CM4) {		/* IGRF only */
						for (j = 0; j < Ctrl->F.n_field_components; j++) 
							Ctrl->DATA.out_field[i*n_field_components+j] = IGRF[Ctrl->F.field_components[j]];
					}
					else {				/* Store the IGRF x,y,z components for later use */
						for (j = 0; j < 3; j++) 
							igrf_xyz[i*3+j] = IGRF[Ctrl->F.field_components[j]];
					}
				}
			}

			if (do_CM4) {				/* DO CM4 only. Eval CM4 at all points */
				if ((error = MGD77_cm4field (Ctrl, T->segment[s]->coord[GMT_X], T->segment[s]->coord[GMT_Y],
					alt_array, time_array))) {
					fprintf(stderr, "Error: this segment has a record generating error.\n"
						"Unfortunately, this means all other eventually good records\n"
						"are also ignored. Fix the bad record and rerun the command.\n");
					continue;
				}
			}

			if ((do_CM4 || do_IGRF) && !joint_IGRF_CM4) {	/* DID CM4 or (exclusive) IGRF only. */
				for (i = 0; i < T->segment[s]->n_rows; i++) {	/* Output the requested columns */
					n_out = 0;
					if (copy_input) for (j = 0; j < T->segment[s]->n_columns; j++) out[n_out++] = T->segment[s]->coord[j][i];
					for (j = 0; j < n_field_components; j++)
						out[n_out++] = Ctrl->DATA.out_field[i*n_field_components+j];

					GMT_output (GMT_stdout, n_out, out);
				}
			}
			else {					/* DID CM4 and IGRF */
				double x, y, z;
				for (i = 0; i < T->segment[s]->n_rows; i++) {	/* Output the requested columns */
					n_out = 0;
					if (copy_input) for (j = 0; j < T->segment[s]->n_columns; j++) out[n_out++] = T->segment[s]->coord[j][i];
					if (cm4_igrf_T) {
						x = Ctrl->DATA.out_field[i*3  ] + igrf_xyz[i*3  ];
						y = Ctrl->DATA.out_field[i*3+1] + igrf_xyz[i*3+1];
						z = Ctrl->DATA.out_field[i*3+2] + igrf_xyz[i*3+2];
						out[n_out++] = sqrt(x*x + y*y + z*z);
					}
					else {
						for (j = 0; j < 3; j++)
							out[n_out++] = Ctrl->DATA.out_field[i*3+j] + igrf_xyz[i*3+j];
					}

					GMT_output (GMT_stdout, n_out, out);
				}
			}

		}
		GMT_free_table (T);
		if (fp != GMT_stdin) GMT_fclose (fp);
	}
	
	GMT_free ((void *)Ctrl->DATA.out_field);
	if (!(time_in_years || fixed_time)) GMT_free ((void *)time_years);
	if (joint_IGRF_CM4) GMT_free ((void *)igrf_xyz);
	
	Free_CM4_Ctrl (Ctrl);	/* Deallocate control structure */
	
	MGD77_end (&M);
	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_CM4_Ctrl () {	/* Allocate and initialize a new control structure */
	struct MGD77_CM4 *C;
	
	C = (struct MGD77_CM4 *) calloc ((size_t)1, sizeof (struct MGD77_CM4));
	return ((void *)C);
}

void Free_CM4_Ctrl (struct MGD77_CM4 *C) {	/* Deallocate control structure */
	free ((void *)C);	
}
