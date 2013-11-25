/*--------------------------------------------------------------------
 *	$Id: xyz2grd.c 9987 2013-01-30 21:35:50Z pwessel $
 *
 *	Copyright (c) 1991-2013 by P. Wessel and W. H. F. Smith
 *	See LICENSE.TXT file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation; version 2 or any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/
/*
 * xyz2grd reads a xyz file from standard input and creates the
 * corresponding grd-file. The input file does not have to be
 * sorted. xyz2grd will report if some nodes are missing.
 *
 * Author:	Paul Wessel
 * Date:	3-JAN-1991-2000
 * Version:	2.0	Based on old v1.x
 * Version:	3.1	Based on old 3.0
 * Version:	3.2	Added -S option for swapping only
 *		3.3.1	-S option need to set fopen modes to rb an wb
 *		3.3.4	-L option finally fixed
 *		3.3.5	-L option finally fixed, part II.  Also added -A[z|n]
 * Version:	3.3.6
 *		3.4	Now can handle multiple input files (except with -S)
 *		4
 */
 
#include "gmt.h"

#include "gmt_parse_z_io.h"	/* To define the Z structure used for parsing */

struct XYZ2GRD_CTRL {
	struct A {	/* -A[n|z|u|l] */
		GMT_LONG active;
		char mode;
	} A;
	struct D {	/* -D<xname>/<yname>/<zname>/<scale>/<offset>/<title>/<remark> */
		GMT_LONG active;
		char *information;
	} D;
	struct E {	/* -E[<nodata>] */
		GMT_LONG active;
		GMT_LONG set;
		double nodata;
	} E;
	struct F {	/* -F */
		GMT_LONG active;
	} F;
	struct G {	/* -G<output_grdfile> */
		GMT_LONG active;
		char *file;
	} G;
	struct I {	/* -Idx[/dy] */
		GMT_LONG active;
		double xinc, yinc;
	} I;
	struct N {	/* -N<nodata> */
		GMT_LONG active;
		double value;
	} N;
	struct S {	/* -S */
		GMT_LONG active;
		char *file;
	} S;
	struct Z Z;
};

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE, do_grid = FALSE, count = FALSE, average = TRUE;
	GMT_LONG b_only = FALSE, done = FALSE, nofile = TRUE;
	
	GMT_LONG z = 2, fno, n_args, *flag = VNULL;
	GMT_LONG n_expected_fields, n_fields, n_files = 0, high_low = 0;
	GMT_LONG i, ii, jj , gmt_ij, nm = 0, n_read = 0, n_filled = 0, n_used = 0;
	GMT_LONG ij = -1;	/* Will be incremented to 0 or set first time around */
	GMT_LONG n_empty = 0, n_stuffed = 0, n_bad = 0, n_confused = 0;

	double *in, d_value, wesn[4];

	float no_data_f, *a = VNULL;

	char line[BUFSIZ];

	FILE *fp = NULL;

	struct GRD_HEADER grd;
	struct GMT_Z_IO io;
	struct XYZ2GRD_CTRL *Ctrl = NULL;
	
	void *New_xyz2grd_Ctrl (), Free_xyz2grd_Ctrl (struct XYZ2GRD_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct XYZ2GRD_CTRL *)New_xyz2grd_Ctrl ();	/* Allocate and initialize a new control structure */

	GMT_grd_init (&grd, argc, argv, FALSE);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'b':
					b_only = TRUE;
				case 'H':
				case 'R':
				case 'V':
				case ':':
				case 'f':
				case '\0':
					error += (GMT_LONG)GMT_parse_common_options (argv[i], &grd.x_min, &grd.x_max, &grd.y_min, &grd.y_max);
					break;

				case 'A':
					if (argv[i][2] && !strchr ("nluz", argv[i][2])) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -A option:  Select -An, -Al, -Au, or -A[z]\n", GMT_program);
						error++;
					}
					else {
						Ctrl->A.active = TRUE;
						Ctrl->A.mode = argv[i][2];
					}
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					Ctrl->D.information = strdup (&argv[i][2]);
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					if (argv[i][2]) {
						Ctrl->E.nodata = atof (&argv[i][2]);
						Ctrl->E.set = TRUE;
					}
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					if (GMT_getinc (&argv[i][2], &Ctrl->I.xinc, &Ctrl->I.yinc)) {
						GMT_inc_syntax ('I', 1);
						error = TRUE;
					}
					break;
				case 'L':	/* Obsolete, but backward compatibility prevails [use -f instead] */
					GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
					GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
					fprintf (stderr, "%s: Option -L is obsolete (but is processed correctly).  Please use -f instead\n", GMT_program);
					break;
				case 'N':
					if (!argv[i][2]) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Must specify value or NaN\n", GMT_program);
						error++;
					}
					else {
						Ctrl->N.value = (argv[i][2] == 'N' || argv[i][2] == 'n') ? GMT_d_NaN : atof (&argv[i][2]);
						Ctrl->N.active = TRUE;
					}
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					if (argv[i][2]) Ctrl->S.file = strdup (&argv[i][2]);
					break;
				case 'Z':
					Ctrl->Z.active = TRUE;
					error += (GMT_LONG)GMT_parse_z_io (&argv[i][2], &Ctrl->Z);
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
		fprintf (stderr, "xyz2grd %s - Converting [xy]z data to a GMT grid file\n\n", GMT_VERSION);
		fprintf (stderr, "usage: xyz2grd [<[xy]zfile(s)>] -G<grdfile> %s %s\n", GMT_I_OPT, GMT_Rgeo_OPT);
		fprintf (stderr, "\t[-A[n|z|u|l]] [%s] [E[<nodata>]] [-F] [%s]\n", GMT_GRDEDIT, GMT_Ho_OPT);
		fprintf (stderr, "\t[-N<nodata>] [-S[<zfile]] [-V] [-Z[<flags>]] [%s] [%s] [%s]\n", GMT_t_OPT, GMT_bi_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\txyzfile(s) is a 1- or 3-column file [or standard input] with [x,y,]z values.\n");
		fprintf (stderr, "\t-G to name the output grid file.\n");
		GMT_inc_syntax ('I', 0);
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A (or -Az): Sum up multiple entries at the same node.\n");
		fprintf (stderr, "\t   Append n (-An): Count number of multiple entries per node instead.\n");
		fprintf (stderr, "\t   Append u (-Au): Keep maximum value if multiple entries per node.\n");
		fprintf (stderr, "\t   Append l (-Al): Keep minimum value if multiple entries per node.\n");
		fprintf (stderr, "\t   [Default (no -A option) will compute mean values].\n");
		fprintf (stderr, "\t-D to enter header information.  Specify '=' to get default value.\n");
		fprintf (stderr, "\t-E input is an ESRI ArcInfo interchange ASCII grid.  Optionally,\n");
		fprintf (stderr, "\t   append the nodata value that represents nodes that should be NaN\n");
		fprintf (stderr, "\t   [Normally, a nodata_value declaration may be found in the ESRI grid].\n");
		fprintf (stderr, "\t   Information for -R , -I, and -F are obtained from the ESRI file.\n");
		fprintf (stderr, "\t-F will force pixel registration [Default is grid registration].\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-N set value for nodes without input xyz triplet [Default is NaN].\n");
		fprintf (stderr, "\t   Z-table entries that equal <nodata> are replaced by NaN.\n");
		fprintf (stderr, "\t-S Requires -Z and will swap the byte-order of the input data and write\n");
		fprintf (stderr, "\t   the result to <zfile> (or stdout if no file given).  No grid file created!\n");
		fprintf (stderr, "\t   For this option, only one input file (or stdin) is allowed.\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t-Z sets exact specification of incoming 1-column z-table.\n");
		fprintf (stderr, "\t   If data is in row format, state if first row is at T(op) or B(ottom).\n");
		fprintf (stderr, "\t     Then, append L or R to indicate starting point in row.\n");
		fprintf (stderr, "\t   If data is in column format, state if first columns is L(left) or R(ight).\n");
		fprintf (stderr, "\t     Then, append T or B to indicate starting point in column.\n");
		fprintf (stderr, "\t   To skip a header of size <n> bytes, append s<n> [<n> = 0].\n");
		fprintf (stderr, "\t   To swap byte-order in 2-byte words, append w.\n");
		fprintf (stderr, "\t   Append x if gridline-registered, periodic data in x without repeating column at xmax.\n");
		fprintf (stderr, "\t   Append y if gridline-registered, periodic data in y without repeating row at ymax.\n");
		fprintf (stderr, "\t   Specify one of the following data types (all binary except a):\n");
		fprintf (stderr, "\t     A  Ascii (multiple floating point values per record).\n");
		fprintf (stderr, "\t     a  Ascii (one value per record).\n");
		fprintf (stderr, "\t     c  signed 1-byte character.\n");
		fprintf (stderr, "\t     u  unsigned 1-byte character.\n");
		fprintf (stderr, "\t     h  signed short 2-byte integer.\n");
		fprintf (stderr, "\t     H  unsigned short 2-byte integer.\n");
		fprintf (stderr, "\t     i  signed 4-byte integer.\n");
		fprintf (stderr, "\t     I  unsigned 4-byte integer.\n");
		fprintf (stderr, "\t     l  signed long (4- or 8-byte) integer.\n");
		fprintf (stderr, "\t     f  4-byte floating point single precision.\n");
		fprintf (stderr, "\t     d  8-byte floating point double precision.\n");
		fprintf (stderr, "\t   [Default format is scanline orientation in ASCII representation: -ZTLa].\n");
		fprintf (stderr, "\t   This option assumes all nodes have data values.\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 3 input columns.\n");
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&Ctrl->I.xinc, &Ctrl->I.yinc, &Ctrl->F.active, &Ctrl->I.active);

	if (Ctrl->S.active && !Ctrl->Z.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR option -S:  Must also specify -Z\n", GMT_program);
		error++;
	}
	if (Ctrl->S.active) {	/* Reading and writing binary file */
		if (n_files > 1) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  -S can only handle one input file\n", GMT_program);
			error++;
		}
		else {
			strcpy (GMT_io.r_mode, "rb");
			strcpy (GMT_io.w_mode, "wb");
		}
	}

	GMT_init_z_io (Ctrl->Z.format, Ctrl->Z.repeat, Ctrl->Z.swab, Ctrl->Z.skip, Ctrl->Z.type, &io);
	if (b_only && Ctrl->Z.active) GMT_io.binary[GMT_IN] = FALSE;

	do_grid = !(Ctrl->S.active || Ctrl->E.active);
	if (do_grid) {
		if (!project_info.region_supplied) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
			error++;
		}
		if (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
			error++;
		}
	}
	if (!Ctrl->S.active && !(Ctrl->G.active || Ctrl->G.file)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR option -G:  Must specify output file\n", GMT_program);
		error++;
	}

	if ((GMT_io.binary[GMT_IN] || io.binary) && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 3;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 3) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 3 columns\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	no_data_f = (float)Ctrl->N.value;
	
	if (do_grid) {
		grd.x_inc = Ctrl->I.xinc;
		grd.y_inc = Ctrl->I.yinc;
		grd.node_offset = (int)Ctrl->F.active;
		GMT_RI_prepare (&grd);	/* Ensure -R -I consistency and set nx, ny */
	}
	
	if (Ctrl->A.active) {
		switch (Ctrl->A.mode) {
			case 'n':	/* Count the number of values for each node */
				count = TRUE;
				average = FALSE;
				break;
			case 'l':	/* Return the lowest (minimum) value at each node */
				high_low = -1;
				break;
			case 'u':	/* Return the upper (maximum)  value at each node */
				high_low = +1;
				break;
			case 'z':	/* Return the sum of multiple values at each node */
			case '\0':
				average = count = FALSE;
				break;
			default:
				fprintf (stderr, "%s: GMT SYNTAX ERROR -A option:  Select -An, -Al, -Au, or -A[z]\n", GMT_program);
				error++;
				break;
		}
	}
	if (b_only && Ctrl->Z.active) fprintf (stderr, "%s: GMT Warning.  -Z overrides -bi\n", GMT_program);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

	if (GMT_io.binary[GMT_IN] && GMT_io.in_col_type[2] & GMT_IS_RATIME && GMT_io.single_precision[GMT_IN]) {
		fprintf (stderr, "%s: Warning: Your single precision binary input data are unlikely to hold absolute time coordinates without serious truncation.\n", GMT_program);
		fprintf (stderr, "%s: Warning: You must use double precision when storing absolute time coordinates in binary data tables.\n", GMT_program);
	}

	if (Ctrl->D.active) GMT_decode_grd_h_info (Ctrl->D.information, &grd);

	if (Ctrl->S.active)
		io.swab = TRUE;
	else if (!Ctrl->E.active) {
		grd.xy_off = 0.5 * grd.node_offset;

		GMT_err_fail (GMT_grd_RI_verify (&grd, 1), Ctrl->G.file);

		if (gmtdefs.verbose) fprintf (stderr, "%s: nx = %d  ny = %d\n", GMT_program, grd.nx, grd.ny);

		nm = GMT_get_nm (grd.nx, grd.ny);

		a = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
		flag = (GMT_LONG *) GMT_memory (VNULL, (size_t)nm, sizeof (GMT_LONG), GMT_program);

		n_expected_fields = (Ctrl->Z.active) ? 1 : Ctrl->A.mode == 'n' ? 2 : 3;
		z = (Ctrl->Z.active) ? 0 : 2;
	}

	GMT_err_fail (GMT_set_z_io (&io, &grd), Ctrl->G.file);

	GMT_set_xy_domain (wesn, &grd);	/* May include some padding if gridline-registered */

	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;
	n_args = (argc > 1) ? argc : 2;

	for (fno = 1; !done && fno < n_args; fno++) {	/* Loop over input files, if any */
		if (!nofile && argv[fno][0] == '-') continue;

		if (nofile) {	/* Just read standard input */
			fp = GMT_stdin;
			done = TRUE;
			if (gmtdefs.verbose) fprintf (stderr, "%s: Reading from standard input\n", GMT_program);
#ifdef SET_IO_MODE
			GMT_setmode (GMT_IN);
#endif
		}
#ifdef WIN32
		else if (Ctrl->E.active && (fp = fopen (argv[fno], GMT_io.r_mode)) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[fno]);
			continue;
		}
		else if ((fp = GMT_fopen (argv[fno], GMT_io.r_mode)) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[fno]);
			continue;
		}
#else
		else if ((fp = GMT_fopen (argv[fno], GMT_io.r_mode)) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[fno]);
			continue;
		}
#endif

		if (!nofile && gmtdefs.verbose) fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[fno]);

		if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) GMT_fgets (line, BUFSIZ, fp);

		if (Ctrl->S.active) {	/* Just swap data and exit */
			FILE *fps = NULL;

			if (gmtdefs.verbose) fprintf (stderr, "%s: Swapping data bytes only\n", GMT_program);

			if (!Ctrl->S.file) {
				fps = GMT_stdout;
#ifdef SET_IO_MODE
				GMT_setmode (GMT_OUT);
#endif
			}
			else if ((fps = GMT_fopen (Ctrl->S.file, GMT_io.w_mode)) == NULL) {
				fprintf (stderr, "%s: Cannot create file %s\n", GMT_program, Ctrl->S.file);
				exit (EXIT_FAILURE);
			}

			if (io.skip) GMT_fseek (fp, (long)io.skip, SEEK_CUR);
			while ((io.read_item) (fp, &d_value)) (io.write_item) (fps, d_value);

			if (fp != GMT_stdin) GMT_fclose (fp);
			if (fps != GMT_stdout) GMT_fclose (fps);

			exit (EXIT_SUCCESS);	/* We are done here */
		}

		if (Ctrl->Z.active) {	/* Read separately because of all the possible formats */
			if (Ctrl->N.active && GMT_is_dnan (Ctrl->N.value)) Ctrl->N.active = FALSE;	/* No point testing */
			if (io.skip) GMT_fseek (fp, (long)io.skip, SEEK_CUR);
			while ((io.read_item) (fp, &d_value)) {
				ij++;
				if (ij == io.n_expected) {
					fprintf (stderr, "%s: More than %ld records, only %ld was expected (aborting)!\n", GMT_program, ij, io.n_expected);
					exit (EXIT_FAILURE);
				}
				(io.get_gmt_ij) (&io, ij, &gmt_ij);
				if (Ctrl->N.active && d_value == Ctrl->N.value)
					a[gmt_ij] = GMT_f_NaN;
				else
					a[gmt_ij] = (float)d_value;
			}
			ij++;
			if (ij != io.n_expected) {	/* Input amount doesn't match expectations */
				fprintf (stderr, "%s: Found %ld records, but %ld was expected (aborting)!\n", GMT_program, ij, io.n_expected);
				exit (EXIT_FAILURE);
			}

			GMT_check_z_io (&io, a);	/* This fills in missing periodic row or column */
		}
		else if (Ctrl->E.active) {	/* Read an ESRI Arc Interchange grid format in ASCII */
			GMT_LONG n_left, j;
			double value;
			
			grd.node_offset = 0;
			(void)fgets (line, BUFSIZ, fp);
			GMT_chop (line);					/* Rid the world of CR/LF */
			if (sscanf (line, "%*s %d", &grd.nx) != 1) {
				fprintf (stderr, "%s: Error decoding ncols record\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			(void)fgets (line, BUFSIZ, fp);
			GMT_chop (line);					/* Rid the world of CR/LF */
			if (sscanf (line, "%*s %d", &grd.ny) != 1) {
				fprintf (stderr, "%s: Error decoding ncols record\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			(void)fgets (line, BUFSIZ, fp);	GMT_str_tolower (line);
			GMT_chop (line);					/* Rid the world of CR/LF */
			if (sscanf (line, "%*s %lf", &grd.x_min) != 1) {
				fprintf (stderr, "%s: Error decoding xll record\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			if (!strncmp (line, "xllcorner", (size_t)9)) grd.node_offset = 1;	/* Pixel grid */
			(void)fgets (line, BUFSIZ, fp);	GMT_str_tolower (line);
			GMT_chop (line);					/* Rid the world of CR/LF */
			if (sscanf (line, "%*s %lf", &grd.y_min) != 1) {
				fprintf (stderr, "%s: Error decoding yll record\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			if (!strncmp (line, "yllcorner", (size_t)9)) grd.node_offset = 1;	/* Pixel grid */
			(void)fgets (line, BUFSIZ, fp);
			GMT_chop (line);					/* Rid the world of CR/LF */
			if (sscanf (line, "%*s %lf", &grd.x_inc) != 1) {
				fprintf (stderr, "%s: Error decoding cellsize record\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			grd.y_inc = grd.x_inc;
			grd.xy_off = 0.5 * grd.node_offset;
			grd.x_max = grd.x_min + (grd.nx - 1 + grd.node_offset) * grd.x_inc;
			grd.y_max = grd.y_min + (grd.ny - 1 + grd.node_offset) * grd.y_inc;
			
			GMT_err_fail (GMT_grd_RI_verify (&grd, 1), Ctrl->G.file);

			if (gmtdefs.verbose) fprintf (stderr, "%s: nx = %d  ny = %d\n", GMT_program, grd.nx, grd.ny);
			nm = n_left = GMT_get_nm (grd.nx, grd.ny);

			a = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
			if (!Ctrl->E.set) {	/* Only read if not specified in -E */
				(void)fgets (line, BUFSIZ, fp);
				GMT_chop (line);					/* Rid the world of CR/LF */
				if (sscanf (line, "%*s %lf", &Ctrl->E.nodata) != 1) {
					fprintf (stderr, "%s: Error decoding nodata_value record\n", GMT_program);
					exit (EXIT_FAILURE);
				}
			}
			i = j = 0;
			fscanf (fp, "%s", line);	GMT_str_tolower (line);
			if (!strcmp (line, "nodata_value")) {	/* Found the optional nodata word */
				fscanf (fp, "%lf", &value);
				if (Ctrl->E.set && !GMT_IS_ZERO (value - Ctrl->E.nodata)) {
					fprintf (stderr, "%s: Your -E%g overrides the nodata_value of %g found in the ESRI file\n", GMT_program, Ctrl->E.nodata, value);
				}
				else
					Ctrl->E.nodata = value;
			}
			else {	/* Instead got the very first data value */
				value = atof (line);
				a[0] = (value == Ctrl->E.nodata) ? GMT_f_NaN : (float)value;
				if (++i == grd.nx) i = 0, j++;
				n_left--;
			}
			/* ESRI grids are scanline oriented (top to bottom) as are GMT grids */
			while (fscanf (fp, "%lf", &value) == 1 && n_left) {
				ij = j * grd.nx + i;
				a[ij] = (value == Ctrl->E.nodata) ? GMT_f_NaN : (float) value;
				if (++i == grd.nx) i = 0, j++;
				n_left--;
			}
			if (n_left) {
				fprintf (stderr, "%s: Expected %ld points, found only %ld\n", GMT_program, nm, nm - n_left);
				exit (EXIT_FAILURE);
			}
		}
		else {	/* Get x, y, z */
			gmtdefs.nan_is_gap = FALSE;	/* Here, all x,y matters */
			while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

				n_read++;
				if (GMT_io.status & GMT_IO_MISMATCH) {
					fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
					continue;
				}

				if (GMT_y_is_outside (in[GMT_Y],  wesn[2], wesn[3])) continue;	/* Outside y-range */
				if (GMT_x_is_outside (&in[GMT_X], wesn[0], wesn[1])) continue;	/* Outside x-range */

				/* Ok, we are inside the region - process data */

				ii = GMT_x_to_i (in[GMT_X], grd.x_min, grd.x_inc, grd.xy_off, grd.nx);
				if (ii == -1) ii++, n_confused++;
				if (ii == grd.nx) ii--, n_confused++;
				jj = GMT_y_to_j (in[GMT_Y], grd.y_min, grd.y_inc, grd.xy_off, grd.ny);
				if (jj == -1) jj++, n_confused++;
				if (jj == grd.ny) jj--, n_confused++;
				ij = GMT_IJ (jj, ii, grd.nx);
				if (high_low) {	/* Always come here if looking for extreme values */
					if (flag[ij]) {	/* Already assigned the first value */
						if ((high_low == -1 && (in[z] < (double)a[ij])) || (high_low == +1 && (in[z] > (double)a[ij]))) a[ij] = (float)in[z];
					}
					else {	/* First time, just assign the current value */
						a[ij] = (float)in[z];
						flag[ij] = 1;
					}
				}
				else { 	/* Add up incase we must average */
					a[ij] += (float)in[z];
					flag[ij]++;
				}
				n_used++;
			}
		}
#ifdef WIN32
		if (fp != GMT_stdin) {
			(Ctrl->E.active) ? fclose (fp) : GMT_fclose (fp);
		}
#else
		if (fp != GMT_stdin) GMT_fclose (fp);
#endif
	}

	if (!(Ctrl->Z.active || Ctrl->E.active)) {	/* xyz data could have resulted in duplicates */
		if (GMT_io.in_col_type[0] == GMT_IS_LON && GMT_360_RANGE (grd.x_max, grd.x_min) && !grd.node_offset) {	/* Make sure longitudes got replicated */
			GMT_LONG j, ij_west, ij_east;

			for (j = 0; j < grd.ny; j++) {	/* For each row, look at west and east bin */
				ij_west = j * grd.nx;
				ij_east = ij_west + grd.nx - 1;
				
				if (flag[ij_west] && !flag[ij_east]) {		/* Nothing in east bin, just copy from west */
					a[ij_east] = a[ij_west];
					flag[ij_east] = flag[ij_west];
				}
				else if (flag[ij_east] && !flag[ij_west]) {	/* Nothing in west bin, just copy from east */
					a[ij_west] = a[ij_east];
					flag[ij_west] = flag[ij_east];
				}
				else {	/* Both have some stuff, consolidate combined value into the west bin, then replicate to the east */
					if (high_low) {	/* Always come here if looking for extreme values */
						if ((high_low == -1 && (a[ij_east] < a[ij_west])) || (high_low == +1 && (a[ij_east] > a[ij_west]))) a[ij_west] = a[ij_east];
					}
					else { 	/* Add up incase we must average */
						a[ij_west] += a[ij_east];
						flag[ij_west] += flag[ij_east];
					}
					/* Replicate: */
					a[ij_east] = a[ij_west];
					flag[ij_east] = flag[ij_west];
				}

			}
		}

		for (ij = 0; ij < nm; ij++) {	/* Check if all nodes got one value only */
			if (flag[ij] == 1) {	/* This catches nodes with one value or the -Al|u single values */
				if (count) a[ij] = 1.0;
				n_filled++;
			}
			else if (flag[ij] == 0) {
				n_empty++;
				a[ij] = no_data_f;
			}
			else {	/* More than 1 value went to this node */
				if (count)
					a[ij] = (float)flag[ij];
				else if (average)
					a[ij] /= (float)flag[ij];
				/* implicit else means return the sum of the values */
				n_filled++;
				n_stuffed++;
			}
		}
		
		if (gmtdefs.verbose) {
			sprintf (line, "%s\n", gmtdefs.d_format);
			fprintf (stderr, "%s:  n_read: %ld  n_used: %ld  n_filled: %ld  n_empty: %ld set to ", GMT_program,
				n_read, n_used, n_filled, n_empty);
			(GMT_is_dnan (Ctrl->N.value)) ? fprintf (stderr, "NaN\n") : fprintf (stderr, line, Ctrl->N.value);
			if (n_bad) fprintf (stderr, "%s: %ld records unreadable\n", GMT_program, n_bad);
			if (n_stuffed && !count) fprintf (stderr, "%s: Warning - %ld nodes had multiple entries that were processed\n", GMT_program, n_stuffed);
			if (n_confused) fprintf (stderr, "%s: Warning - %ld values gave bad indices: Pixel vs gridline confusion?\n", GMT_program, n_confused);
		}
	}

	GMT_err_fail (GMT_write_grd (Ctrl->G.file, &grd, a, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);

	GMT_free ((void *)a);
	if (!Ctrl->E.active) GMT_free ((void *)flag);

	Free_xyz2grd_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_xyz2grd_Ctrl () {	/* Allocate and initialize a new control structure */
	struct XYZ2GRD_CTRL *C;
	
	C = (struct XYZ2GRD_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct XYZ2GRD_CTRL), "New_xyz2grd_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->N.value = GMT_d_NaN;
	C->Z.type = 'a';
	C->Z.format[0] = 'T';	C->Z.format[1] = 'L';
	return ((void *)C);
}

void Free_xyz2grd_Ctrl (struct XYZ2GRD_CTRL *C) {	/* Deallocate control structure */
	if (C->D.information) free ((void *)C->D.information);	
	if (C->G.file) free ((void *)C->G.file);	
	if (C->S.file) free ((void *)C->S.file);	
	GMT_free ((void *)C);	
}
