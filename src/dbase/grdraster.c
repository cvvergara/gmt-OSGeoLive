/*--------------------------------------------------------------------
 *	$Id: grdraster.c 9923 2012-12-18 20:45:53Z pwessel $
 *
 *    Copyright (c) 1991-2013 by P. Wessel and W. H. F. Smith
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/* grdraster.c -- read a rasterfile and extract a region as a grid file.
 *
 * This is a complete rewrite for GMT version 3.0.  This is based on
 * the earlier versions written by me for installations at Scripps and
 * NOAA, and does not resemble the version grdraster.c_supplied which
 * was used at Hawaii.
 *
 * Author:	Walter H. F. Smith
 * Date:	20 January, 1995
 * Update:	10 March 1999 PW: Now use new default dir $GMT_SHAREDIR/dbase
 *		05 April 1999 PW: Now deals with DOS/UNIX directory signs
 *		07 April 1999 PW: ...and DOS drive letters.
 *		27 April 1999 PW: Added option to swap via #define GMTSWAP
 *				  which can manually be set in the makefile.
 *		15 June 1999  PW: Added features to add a last column in
 *				  grdraster.info which may hold the character
 *				  L (for Littleendian) or B (for Bigendian).
 *				  If present, and different from the endianness
 *				  of the current machine, byte-swapping will occur.
 *		18 Oct 1999   PW: Replaced GMTSWAP with WORDS_BIGENDIAN as provided
 *				  by the new configure script.
 *		05 May 2000   PW: Added option -Z [-bo[s]] to write xyz to stdout instead
 *		12 Jan 2001   PW: Dynamically allocate info structures; no longer hardwired.
 *		27 Jun 2005   PW: Allow modifiers G(ographic) and C(artesian) in registration flag.
 *		19 Sep 2005   PW: Allow optional H<bytes> at end of line for skipping headers.
 *		07 Sep 2006   PW: Only do 360 wrap if geographic grid.
 *		27 Apr 2009   PW: Added sanity check on file size (actual vs computed from grdraster.info)
 */

#include "gmt.h"

struct GRDRASTER_CTRL {
	struct G {	/* -G<output_grdfile> */
		GMT_LONG active;
		char *file;
	} G;
	struct I {	/* -Idx[/dy] */
		GMT_LONG active;
		double xinc, yinc;
	} I;
};

struct GRDRASTER_INFO {
	struct GRD_HEADER h;
	GMT_LONG	id;		/* File number  */
	GMT_LONG	nglobal;	/* If not 0, ras is global and i%nglobal makes it periodic  */
	GMT_LONG	nanflag;
	GMT_LONG	nanset;		/* True if raster uses nanflag to signal NaN  */
	GMT_LONG	skip;		/* Skip this number of header bytes when opening file  */
	char	type;
	GMT_LONG swap_me;	/* TRUE if data set need to be swapped */
	GMT_LONG geo;		/* TRUE if we believe x/y is lon/lat, FALSE otherwise */
};

GMT_LONG get_byte_size (char type);

int main (int argc, char **argv)
{
	GMT_LONG xyz_out = FALSE;

#ifndef WORDS_BIGENDIAN
	char my_endian = 'L';	/* This machine is Little endian */
#else
	char my_endian = 'B';	/* This machine is Big endian */
#endif
	char	*buffer = NULL, *r_opt = NULL, *tselect = CNULL, match[GRD_REMARK_LEN];
	unsigned char *ubuffer = NULL;
	static unsigned char maskset[8] = {128, 64, 32, 16, 8, 4, 2, 1};

	GMT_LONG	i, j, k, ksize = 0, iselect, nselected = 0, imult, jmult, nrasters, ij_offset;
	GMT_LONG	irasstart, jrasstart, n_nan, iras, jras, ij, ijras, jseek, nmask = 0;
	GMT_LONG	error = FALSE, firstread, nm;

	float	*grd = NULL, *floatrasrow = NULL;

	double	tol;
	double	grdlatorigin, grdlonorigin, raslatorigin, raslonorigin;
	double *x = NULL, y, out[3];

	struct GRDRASTER_INFO myras;
	struct GRD_HEADER h;
	struct GRDRASTER_INFO *rasinfo = NULL;

	FILE *fp = NULL;

	struct GRDRASTER_CTRL *Ctrl = NULL;

	void *New_grdraster_Ctrl (), Free_grdraster_Ctrl (struct GRDRASTER_CTRL *C);
	void convert_u_row (struct GRDRASTER_INFO ras, float *row, unsigned char *buffer);
	void convert_c_row (struct GRDRASTER_INFO ras, float *row, char *buffer);
	void convert_d_row (struct GRDRASTER_INFO ras, float *row, short unsigned int *buffer);
	void convert_i_row (struct GRDRASTER_INFO ras, float *row, short int *buffer);
	void convert_l_row (struct GRDRASTER_INFO ras, float *row, int *buffer);
	GMT_LONG load_rasinfo (struct GRDRASTER_INFO **rasinfo, char endian);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDRASTER_CTRL *)New_grdraster_Ctrl ();	/* Allocate and initialize a new control structure */

	GMT_grd_init (&h, argc, argv, FALSE);

	if (!(nrasters = load_rasinfo(&rasinfo, my_endian))) {
		fprintf(stderr, "%s:  ERROR reading grdraster.info file.\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	project_info.region_supplied = FALSE;	/* Because () sets it to TRUE by using GMT_parse_common_options to process file regions */

	for (i = 1; !error && i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				case 'R':
					r_opt = argv[i];
				case 'J':
				case 'V':
				case '\0':
					error += GMT_parse_common_options (argv[i], &h.x_min, &h.x_max, &h.y_min, &h.y_max);
					break;

				case 'b':
					error += GMT_parse_b_option (&argv[i][2]);
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
				default:
					error = TRUE;
					break;
			}
		}
		else {
			nselected++;
			tselect = argv[i];
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "grdraster %s - Extract a region from a raster and save in a grid file\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdraster <file number>|<text> %s [-G<grdfilename>] [%s]\n", GMT_Rgeo_OPT, GMT_Id_OPT);
		fprintf (stderr, " [%s]\n", GMT_bo_OPT);

		fprintf (stderr, "\t<file number> (#) or <text> corresponds to one of the datasets listed.\n");
		fprintf (stderr, "\t[<text> can be a unique substring of the description].\n\n");
		fprintf (stderr, "#	Data Description	Unit	Coverage		Spacing	Registration\n");
		fprintf (stderr, "------------------------------------------------------------------------------------\n");
		for (i = 0; i < nrasters; i++) fprintf (stderr, "%s\n", rasinfo[i].h.command);
		fprintf (stderr, "------------------------------------------------------------------------------------\n\n");
#ifndef WORDS_BIGENDIAN
		fprintf (stderr, "grdraster default binary byte order is Little-endian\n");
#else
		fprintf (stderr, "grdraster default binary byte order is Big-endian\n");
#endif
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t-R specifies the west, east, south, and north edges of the area.\n");
		fprintf (stderr, "\t   Use dd:mm format for regions given in degrees and minutes.\n");
		fprintf (stderr, "\t   Append r if -R specifies the longitudes/latitudes of the lower left\n");
		fprintf (stderr, "\t   and upper right corners of a rectangular area.  If r is appended\n");
		fprintf (stderr, "\t   you must also specify a projection with -J (set scale = 1).\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-G sets the filename for output grid.  If no file is given then grdraster\n");
		fprintf (stderr, "\t   will instead write ASCII (or binary, see -b) xyz triplets to stdout.\n");
		fprintf (stderr, "\t-I specifies the sampling interval of the grid [Default is raster spacing].\n");
		fprintf (stderr, "\t   Give -Idx or -Idx/dy if dy not equal dx.  Append m for minutes.\n");
		fprintf (stderr, "\t   (-I does not do any filtering; it just sub-samples the raster.)\n");
		GMT_explain_option('j');
		GMT_explain_option('V');
		GMT_explain_option ('o');
		fprintf (stderr, "\t   This option only applies if no gridfile is given (see -G)\n");
		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&Ctrl->I.xinc, &Ctrl->I.yinc, NULL, &Ctrl->I.active);

	/* Check that arguments were valid:  */
	if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (Ctrl->I.active && (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0) ) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
                error++;
        }
        if (!Ctrl->G.file) {
                if (gmtdefs.verbose) fprintf (stderr, "%s: No grid file given - will write xyz to stdout\n", GMT_program);
                xyz_out = TRUE;
        }
	if (nselected != 1) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  You must specify one and only one raster file number.\n", GMT_program);
		error++;
	}
	else {
		/* Check if given argument is an integer ID.  If so, assign iselect, else set it to -1 */
		GMT_str_toupper (tselect);	/* Make it upper case - which wont affect integers */
		for (j = i = 0; tselect[j] && i == 0; j++) if (!isdigit ((int)tselect[j])) i = 1;
		iselect = (i == 0) ? atoi (tselect) : -1;
		j = -1;
		for (i = 0; !error && i < nrasters; i++) {
			if (iselect != -1) {	/* We gave an integer ID */
				if (rasinfo[i].id == iselect) {
					if (j == -1)
						j = i;
					else {
						fprintf (stderr, "%s ERROR:  At least two rasters have the same file number in grdraster.info\n", GMT_program);
						error++;
					}
				}
			}
			else {	/* We gave a text snippet to match in command */
				strcpy (match, rasinfo[i].h.command);
				GMT_str_toupper (match);	/* Make it upper case  */
				if (strstr (match, tselect)) {	/* Found a matching text string */
					if (j == -1)
						j = i;
					else {
						fprintf (stderr, "%s ERROR:  At least two rasters have the same text [%s] in grdraster.info\n", GMT_program, tselect);
						error++;
					}
				}
			}
		}
		if (j == -1) {
			if (iselect != -1)
				fprintf (stderr, "%s ERROR:  No raster with file number %ld in grdraster.info\n", GMT_program, iselect);
			else
				fprintf (stderr, "%s ERROR:  No raster with text %s in grdraster.info\n", GMT_program, tselect);
			error++;
		}
		else {
			myras = rasinfo[j];
		}
	}
	if (error) exit (EXIT_FAILURE);

	GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = (myras.geo) ? GMT_IS_LON : GMT_IS_FLOAT;
	GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = (myras.geo) ? GMT_IS_LAT : GMT_IS_FLOAT;

	/* Everything looks OK so far.  If (Ctrl->I.active) verify that it will work, else set it.  */
	if (Ctrl->I.active) {
		h.x_inc = Ctrl->I.xinc;
		h.y_inc = Ctrl->I.yinc;
		tol = 0.01 * myras.h.x_inc;
		imult = irint(h.x_inc / myras.h.x_inc);
		if (imult < 1 || fabs(h.x_inc - imult * myras.h.x_inc) > tol) error++;
		tol = 0.01 * myras.h.y_inc;
		jmult = irint(h.y_inc / myras.h.y_inc);
		if (jmult < 1 || fabs(h.y_inc - jmult * myras.h.y_inc) > tol) error++;
		if (error) {
			fprintf(stderr, "%s ERROR:  Your -I option does not create a grid which fits the selected raster.\n", GMT_program);
			fprintf(stderr, "\t%s\n", myras.h.command);
			exit (EXIT_FAILURE);
		}
	}
	else {
		h.x_inc = myras.h.x_inc;
		h.y_inc = myras.h.y_inc;
		imult = jmult = 1;
	}

	if (!project_info.region && project_info.projection != GMT_NO_PROJ) {

		GMT_err_fail (GMT_map_setup (h.x_min, h.x_max, h.y_min, h.y_max), "");

		h.x_min = floor (project_info.w / h.x_inc) * h.x_inc;
		h.x_max = ceil (project_info.e / h.x_inc) * h.x_inc;
		h.y_min = floor (project_info.s / h.y_inc) * h.y_inc;
		h.y_max = ceil (project_info.n / h.y_inc) * h.y_inc;

		if (gmtdefs.verbose && rint (h.x_inc * 60.0) == (h.x_inc * 60.0)) {	/* Spacing in even minutes */
			GMT_LONG w, e, s, n, wm, em, sm, nm;

			w = (GMT_LONG) floor (h.x_min);	wm = (GMT_LONG) irint ((h.x_min - w) * 60.0);
			e = (GMT_LONG) floor (h.x_max);	em = (GMT_LONG) irint ((h.x_max - e) * 60.0);
			s = (GMT_LONG) floor (h.y_min);	sm = (GMT_LONG) irint ((h.y_min - s) * 60.0);
			n = (GMT_LONG) floor (h.y_max);	nm = (GMT_LONG) irint ((h.y_max - n) * 60.0);
			fprintf (stderr, "%s: %s -> -R%ld:%2.2ld/%ld:%2.2ld/%ld:%2.2ld/%ld:%2.2ld\n", GMT_program, r_opt, w, wm, e, em, s, sm, n, nm);
		}
		else if (gmtdefs.verbose)
			fprintf (stderr, "%s: %s -> -R%g/%g/%g/%g\n", GMT_program, r_opt, h.x_min, h.x_max, h.y_min, h.y_max);
	}

	/* Now Enforce that wesn will fit x_inc, y_inc.  Set nx, ny but reset later based on G or P  */
	tol = 0.01 * h.x_inc;
	h.nx = irint((h.x_max - h.x_min)/h.x_inc);
	if (fabs( (h.x_max - h.x_min) - h.x_inc * h.nx) > tol) error++;
	tol = 0.01 * h.y_inc;
	h.ny = irint((h.y_max - h.y_min)/h.y_inc);
	if (fabs( (h.y_max - h.y_min) - h.y_inc * h.ny) > tol) error++;
	if (error) {	/* Must cleanup and give warning */
		h.x_min = floor (h.x_min / h.x_inc) * h.x_inc;
		h.x_max =  ceil (h.x_max / h.x_inc) * h.x_inc;
		h.y_min = floor (h.y_min / h.y_inc) * h.y_inc;
		h.y_max =  ceil (h.y_max / h.y_inc) * h.y_inc;
		h.nx = irint ((h.x_max - h.x_min) / h.x_inc);
		h.ny = irint ((h.y_max - h.y_min) / h.y_inc);
		fprintf(stderr, "%s WARNING:  Your -R option does not create a region divisible by x_inc, y_inc.\n", GMT_program);
		if (GMT_IS_ZERO (rint (h.x_inc * 60.0) - h.x_inc * 60.0)) {	/* Spacing in even minutes */
			GMT_LONG w, e, s, n, wm, em, sm, nm;
			w = (GMT_LONG) floor (h.x_min);	wm = (GMT_LONG) irint ((h.x_min - w) * 60.0);
			e = (GMT_LONG) floor (h.x_max);	em = (GMT_LONG) irint ((h.x_max - e) * 60.0);
			s = (GMT_LONG) floor (h.y_min);	sm = (GMT_LONG) irint ((h.y_min - s) * 60.0);
			n = (GMT_LONG) floor (h.y_max);	nm = (GMT_LONG) irint ((h.y_max - n) * 60.0);
			if (project_info.region)
				fprintf(stderr, "%s WARNING:  Region reset to -R%ld:%2.2ld/%ld:%2.2ld/%ld:%2.2ld/%ld:%2.2ld.\n", GMT_program, w, wm, e, em, s, sm, n, nm);
			else
				fprintf(stderr, "%s WARNING:  Region reset to -R%ld:%2.2ld/%ld:%2.2ld/%ld:%2.2ld/%ld:%2.2ldr\n", GMT_program, w, wm, s, sm, e, em, n, nm);
		}
		else {
			if (project_info.region)
				fprintf(stderr, "%s WARNING:  Region reset to -R%g/%g/%g/%g.\n", GMT_program, h.x_min, h.x_max, h.y_min, h.y_max);
			else
				fprintf(stderr, "%s WARNING:  Region reset to -R%g/%g/%g/%gr.\n", GMT_program, h.x_min, h.y_min, h.x_max, h.y_max);
		}
		error = 0;
	}

	/* Now we are ready to go:  */
	if (!myras.h.node_offset) {
		h.nx++;
		h.ny++;
	}
	strcpy(h.title, myras.h.title);
	strcpy(h.z_units, myras.h.z_units);
	strcpy(h.remark, myras.h.remark);
	if (myras.geo) {
		strcpy(h.x_units, "Longitude [degrees_east]");
		strcpy(h.y_units, "Latitude [degrees_north]");
	}
	else {
		strcpy(h.x_units, "x");
		strcpy(h.y_units, "y");
	}
	h.node_offset = myras.h.node_offset;
	h.z_min = DBL_MAX;
	h.z_max = -DBL_MAX;
	h.xy_off = 0.5 * h.node_offset;
	myras.h.xy_off = 0.5 * myras.h.node_offset;

	grdlatorigin = GMT_j_to_y (0, h.y_min, h.y_max, h.y_inc, h.xy_off, h.ny);
	grdlonorigin = GMT_i_to_x (0, h.x_min, h.x_max, h.x_inc, h.xy_off, h.nx);
	raslatorigin = GMT_j_to_y (0, myras.h.y_min, myras.h.y_max, myras.h.y_inc, myras.h.xy_off, myras.h.ny);
	raslonorigin = GMT_i_to_x (0, myras.h.x_min, myras.h.x_max, myras.h.x_inc, myras.h.xy_off, myras.h.nx);
	irasstart = irint( (grdlonorigin - raslonorigin) / myras.h.x_inc);
	jrasstart = irint( (raslatorigin - grdlatorigin) / myras.h.y_inc);
	if (myras.nglobal) while (irasstart < 0) irasstart += myras.nglobal;
	n_nan = 0;

	/* Get space:  */
	if (xyz_out) {	/* Need just space for one row */
		grd = (float *)GMT_memory(VNULL, (size_t)(h.nx), sizeof(float), GMT_program);
		x = (double *)GMT_memory(VNULL, (size_t)(h.nx), sizeof(double), GMT_program);
		for (i = 0; i < h.nx; i++) x[i] = GMT_i_to_x (i, h.x_min, h.x_max, h.x_inc, h.xy_off, h.nx);
		ij_offset = 0;
	} else {	/* Need entire grid */
		nm = GMT_get_nm (h.nx, h.ny);
		grd = (float *)GMT_memory(VNULL, (size_t)nm, sizeof(float), GMT_program);
		ij_offset = h.nx;
	}

	ksize = get_byte_size (myras.type);
	if (ksize == 0) {	/* Bits; Need to read the whole thing:  */
		nmask = (GMT_LONG)ceil (myras.h.nx * myras.h.ny * 0.125);
		ubuffer = (unsigned char *)GMT_memory(VNULL, (size_t)nmask, (size_t)1, GMT_program);
	}
	else {	/* Need to read by rows, and convert each row to float:  */
		buffer = GMT_memory (VNULL, (size_t)myras.h.nx, (size_t)ksize, GMT_program);
		floatrasrow = (float *)GMT_memory(VNULL, (size_t)myras.h.nx, sizeof(float), GMT_program);
	}

	/* Now open file and do it:  */

	if ( (fp = fopen(myras.h.remark, "rb") ) == NULL) {
		fprintf(stderr,"%s ERROR opening %s for read.\n", GMT_program, myras.h.remark);
		exit (EXIT_FAILURE);
	}
	if (myras.skip && fseek (fp, (long) (myras.skip), SEEK_CUR) ) {
		fprintf(stderr,"%s ERROR skipping %ld bytes in %s.\n", GMT_program, myras.skip, myras.h.remark);
		exit (EXIT_FAILURE);
	}
	if (gmtdefs.verbose) fprintf(stderr, "%s:  Reading from raster %s\n", GMT_program, myras.h.remark);
	if (gmtdefs.verbose && myras.swap_me) fprintf (stderr, "%s:  Data from %s will be byte-swapped\n", GMT_program, myras.h.remark);

	if (myras.type == 'b') {
		if ( (fread((void *)ubuffer, sizeof (unsigned char), (size_t)nmask, fp)) != (size_t)nmask) {
			fprintf(stderr,"%s ERROR:  Failure to read a bitmap raster from %s.\n", GMT_program, myras.h.remark);
			GMT_free ((void *)ubuffer);
			GMT_free ((void *)grd);
			fclose(fp);
			exit (EXIT_FAILURE);
		}
		for (j = 0, jras = jrasstart; j < h.ny; j++, jras += jmult) {
			y = GMT_j_to_y (j, h.y_min, h.y_max, h.y_inc, h.xy_off, h.ny);
			if (jras < 0 || jras > myras.h.ny) {
				/* This entire row is outside the raster:  */
				for (i = 0, ij = j * ij_offset; i < h.nx; i++, ij++) grd[ij] = GMT_f_NaN;
				n_nan += h.nx;
			}
			else {
				iras = irasstart;
				ijras = jras * myras.h.nx;
				for (i = 0, ij = j * ij_offset; i < h.nx; i++, ij++) {
					if (myras.nglobal && iras >= myras.nglobal) iras = iras%myras.nglobal;
					if (iras < 0 || iras >= myras.h.nx) {
						grd[ij] = GMT_f_NaN;
						n_nan++;
					}
					else {
						k = ijras + iras;
						grd[ij] = (float)((ubuffer[k/8] & maskset[k%8]) ? 1.0 : 0.0);
						if (grd[ij] > h.z_max) h.z_max = grd[ij];
						if (grd[ij] < h.z_min) h.z_min = grd[ij];
					}
					iras += imult;
				}
			}
			if (xyz_out) {
				out[1] = y;
				for (i = 0; i < h.nx; i++) {
					out[0] = x[i];
					out[2] = grd[i];
					GMT_output (GMT_stdout, 3, out);
				}
			}
		}
		GMT_free ((void *)ubuffer);
	}
	else {
		firstread = TRUE;
		for (j = 0, jras = jrasstart; j < h.ny; j++, jras += jmult) {
			y = GMT_j_to_y (j, h.y_min, h.y_max, h.y_inc, h.xy_off, h.ny);
			if (jras < 0 || jras > myras.h.ny) {
				/* This entire row is outside the raster:  */
				for (i = 0, ij = j * ij_offset; i < h.nx; i++, ij++) grd[ij] = GMT_f_NaN;
				n_nan += h.nx;
			}
			else {
				if (firstread) {
					jseek = (jras != 0) ? jras : 0;
					firstread = FALSE;
				}
				else if (jmult > 1)
					jseek = jmult - 1;
				else
					jseek = 0;
				/* This will be slow on SGI because seek is broken there */
				if (jseek && fseek (fp, (long) (jseek * ksize * myras.h.nx), SEEK_CUR) ) {
					fprintf(stderr,"%s: ERROR seeking in %s\n", GMT_program, myras.h.remark);
					fclose(fp);
					GMT_free((void *)buffer);
					GMT_free ((void *)grd);
					exit (EXIT_FAILURE);
				}
				if ( (fread((void *)buffer, (size_t)ksize, (size_t)myras.h.nx, fp)) != (size_t)myras.h.nx) {
					fprintf(stderr,"%s: ERROR reading in %s\n", GMT_program, myras.h.remark);
					fclose(fp);
					GMT_free((void *)buffer);
					GMT_free ((void *)grd);
					exit (EXIT_FAILURE);
				}
#ifdef DEBUG
				fprintf (stderr, "%s: Doing line %6.6ld\r", GMT_program, j);
#endif
				switch (myras.type) {
					case 'u':
						convert_u_row(myras, floatrasrow, (unsigned char *)buffer);
						break;
					case 'c':
						convert_c_row(myras, floatrasrow, buffer);
						break;
					case 'd':
						convert_d_row(myras, floatrasrow, (unsigned short int *)buffer);
						break;
					case 'i':
						convert_i_row(myras, floatrasrow, (short int *)buffer);
						break;
					case 'l':
						convert_l_row(myras, floatrasrow, (int *)buffer);
						break;
				}
				iras = irasstart;
				for (i = 0, ij = j * ij_offset; i < h.nx; i++, ij++) {
					if (myras.nglobal && iras >= myras.nglobal) iras = iras%myras.nglobal;
					if (iras < 0 || iras >= myras.h.nx) {
						grd[ij] = GMT_f_NaN;
						n_nan++;
					}
					else {
						grd[ij] = floatrasrow[iras];
						if (grd[ij] > h.z_max) h.z_max = grd[ij];
						if (grd[ij] < h.z_min) h.z_min = grd[ij];
					}
					iras += imult;
				}
			}
			if (xyz_out) {
				out[1] = y;
				for (i = 0; i < h.nx; i++) {
					out[0] = x[i];
					out[2] = grd[i];
					GMT_output (GMT_stdout, 3, out);
				}
			}
		}
		GMT_free ((void *)buffer);
		GMT_free ((void *)floatrasrow);
	}
	fclose(fp);

	if (gmtdefs.verbose) fprintf (stderr, "%s:  Finished reading from %s\n", GMT_program, myras.h.remark);
	if (gmtdefs.verbose) fprintf (stderr, "%s:  min max and # NaN found: %g %g %ld\n", GMT_program, h.z_min, h.z_max, n_nan);

	if (n_nan == h.nx * h.ny) fprintf(stderr,"%s: WARNING - Your grid file is entirely full of NaNs.\n", GMT_program);

	if (xyz_out) {
		GMT_free ((void *)x);
	}
	else
		GMT_err_fail (GMT_write_grd (Ctrl->G.file, &h, grd, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);

	GMT_free ((void *)grd);
	GMT_free ((void *)rasinfo);

	Free_grdraster_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void convert_u_row(struct GRDRASTER_INFO ras, float *row, unsigned char *buffer)
{
	GMT_LONG	i, tempval;
	for (i = 0; i < ras.h.nx; i++) {
		tempval = (GMT_LONG)buffer[i];
		if (ras.nanset && tempval == ras.nanflag) {
			row[i] = GMT_f_NaN;
		}
		else {
			row[i] = (float)tempval;
			if (ras.h.z_scale_factor != 1.0) row[i] *= (float)ras.h.z_scale_factor;
			if (ras.h.z_add_offset != 0.0) row[i] += (float)ras.h.z_add_offset;
		}
	}
	return;
}

void convert_c_row(struct GRDRASTER_INFO ras, float *row, char *buffer)
{
	GMT_LONG	i, tempval;
	for (i = 0; i < ras.h.nx; i++) {
		tempval = (GMT_LONG)buffer[i];
		if (ras.nanset && tempval == ras.nanflag) {
			row[i] = GMT_f_NaN;
		}
		else {
			row[i] = (float)tempval;
			if (ras.h.z_scale_factor != 1.0) row[i] *= (float)ras.h.z_scale_factor;
			if (ras.h.z_add_offset != 0.0) row[i] += (float)ras.h.z_add_offset;
		}
	}
	return;
}

void convert_d_row(struct GRDRASTER_INFO ras, float *row, short unsigned int *buffer)
{
	GMT_LONG	i, tempval;
	for (i = 0; i < ras.h.nx; i++) {
		if (ras.swap_me) buffer[i] = GMT_swab2 (buffer[i]);

		tempval = buffer[i];
		if (ras.nanset && tempval == ras.nanflag) {
			row[i] = GMT_f_NaN;
		}
		else {
			row[i] = (float)tempval;
			if (ras.h.z_scale_factor != 1.0) row[i] *= (float)ras.h.z_scale_factor;
			if (ras.h.z_add_offset != 0.0) row[i] += (float)ras.h.z_add_offset;
		}
	}
	return;
}

void convert_i_row(struct GRDRASTER_INFO ras, float *row, short int *buffer)
{
	GMT_LONG	i, tempval;
	for (i = 0; i < ras.h.nx; i++) {
		if (ras.swap_me) buffer[i] = GMT_swab2 (buffer[i]);

		tempval = buffer[i];
		if (ras.nanset && tempval == ras.nanflag) {
			row[i] = GMT_f_NaN;
		}
		else {
			row[i] = (float)tempval;
			if (ras.h.z_scale_factor != 1.0) row[i] *= (float)ras.h.z_scale_factor;
			if (ras.h.z_add_offset != 0.0) row[i] += (float)ras.h.z_add_offset;
		}
	}
	return;
}

void convert_l_row(struct GRDRASTER_INFO ras, float *row, int *buffer)
{
	GMT_LONG	i, tempval;
	for (i = 0; i < ras.h.nx; i++) {
		if (ras.swap_me) buffer[i] = GMT_swab4 (buffer[i]);

		tempval = buffer[i];
		if (ras.nanset && tempval == ras.nanflag) {
			row[i] = GMT_f_NaN;
		}
		else {
			row[i] = (float)tempval;
			if (ras.h.z_scale_factor != 1.0) row[i] *= (float)ras.h.z_scale_factor;
			if (ras.h.z_add_offset != 0.0) row[i] += (float)ras.h.z_add_offset;
		}
	}
	return;
}

GMT_LONG	load_rasinfo(struct GRDRASTER_INFO **ras, char endian)
{
	/* Read the file grdraster.info
		Store the i'th row of the file in rasinfo[i].h.command.
		Store the filename in rasinfo[i].h.remark.
		Store the description in rasinfo[i].h.title.
		Store the units in rasinfo[i].h.z_units.
		After all has parsed correctly, truncate rasinfo[i].h.command
			so it can be printed out as an abbreviated description
			for the user.
		Figure out if file is global, and set nglobal.
		Set nx and ny also.

	Return 0 if cannot read files correctly, or nrasters if successful.  */

	GMT_LONG	i, j, length, stop_point, nfound = 0, ksize = 0;
	size_t	n_alloc;
	off_t expected_size;
	double	global_lon, lon_tol;
	char	buf[GRD_REMARK_LEN], dir[GRD_REMARK_LEN], *l = NULL;
	FILE	*fp = NULL;
	struct GRDRASTER_INFO *rasinfo = NULL;
	struct GMT_STAT F;

	/* Find and open the file grdraster.info */

	if (!(GMT_getdatapath("grdraster.info", dir) || GMT_getsharepath("dbase", "grdraster", ".info", dir))) {
		fprintf(stderr, "%s: ERROR cannot find file grdraster.info\n", GMT_program);
		return (0);
	}

	if ( (fp = fopen(dir, "r")) == NULL) {
		fprintf(stderr, "%s: ERROR cannot open file %s\n", GMT_program, dir);
		return(0);
	}

	/* Truncate the pathname of grdraster.info to just the directory name */

	if ((l = strstr (dir, "grdraster.info"))) *l = '\0';

	n_alloc = GMT_SMALL_CHUNK;
	rasinfo = (struct GRDRASTER_INFO *) GMT_memory (VNULL, n_alloc, sizeof (struct GRDRASTER_INFO), GMT_program);

	while (fgets(rasinfo[nfound].h.command, GRD_COMMAND_LEN, fp) ) {
		if (rasinfo[nfound].h.command[0] == '#') continue;
		GMT_chop (rasinfo[nfound].h.command);
		if (rasinfo[nfound].h.command[0] == '\0') continue;	 /* Blank line */

		length = strlen(rasinfo[nfound].h.command);

		/* Find the integer file name first:  */
		i = 0;
		while (i < length && (rasinfo[nfound].h.command[i] == ' ' ||  rasinfo[nfound].h.command[i] == '\t') ) i++;
		if (i == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tFile number conversion error.\n");
			return(0);
		}
		j = i+1;
		while (j < length && !(rasinfo[nfound].h.command[j] == ' ' ||  rasinfo[nfound].h.command[j] == '\t') ) j++;
		if (j == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tFile number conversion error.\n");
			return(0);
		}

		strncpy(buf, &rasinfo[nfound].h.command[i], (size_t)j-i);
		buf[j-i] = '\0';
		if ( (sscanf(buf, "%" GMT_LL "d", &rasinfo[nfound].id) ) != 1) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tFile number conversion error.\n");
			return(0);
		}

		/* Now find the title string:  */
		i = j+1;
		while (i < length && (rasinfo[nfound].h.command[i] != '"') ) i++;
		if (i == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tTitle string conversion error.\n");
			return(0);
		}
		j = i+1;
		while (j < length && (rasinfo[nfound].h.command[j] != '"') ) j++;
		if (j == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tTitle string conversion error.\n");
			return(0);
		}
		i++;
		if (i == j) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tTitle string conversion error.\n");
			return(0);
		}
		strncpy(rasinfo[nfound].h.title, &rasinfo[nfound].h.command[i], (size_t)j-i);
		rasinfo[nfound].h.title[j-i] = '\0';

		/* Now find the z_unit string:  */
		i = j+1;
		while (i < length && (rasinfo[nfound].h.command[i] != '"') ) i++;
		if (i == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tUnits string conversion error.\n");
			return(0);
		}
		j = i+1;
		while (j < length && (rasinfo[nfound].h.command[j] != '"') ) j++;
		if (j == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tUnits string conversion error.\n");
			return(0);
		}
		i++;
		if (i == j) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tUnits string conversion error.\n");
			return(0);
		}
		strncpy(rasinfo[nfound].h.z_units, &rasinfo[nfound].h.command[i], (size_t)j-i);
		rasinfo[nfound].h.z_units[j-i] = '\0';

		/* Now find the -R string:  */
		i = j+1;
		while (i < length && (rasinfo[nfound].h.command[i] != '-') ) i++;
		if (i == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\t-R string conversion error.\n");
			return(0);
		}
		j = i+1;
		while (j < length && !(rasinfo[nfound].h.command[j] == ' ' ||  rasinfo[nfound].h.command[j] == '\t') ) j++;
		if (j == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\t-R string conversion error.\n");
			return(0);
		}
		strncpy(buf, &rasinfo[nfound].h.command[i], (size_t)j-i);
		buf[j-i]='\0';
		if (strchr (buf, ':') || strchr (buf, 'W') || strchr (buf, 'E') || strchr (buf, 'S') || strchr (buf, 'N')) {
			GMT_io.in_col_type[0] = GMT_IS_LON;
			GMT_io.in_col_type[1] = GMT_IS_LAT;
		}
		else {
			GMT_io.in_col_type[0] = GMT_IS_FLOAT;
			GMT_io.in_col_type[1] = GMT_IS_FLOAT;
		}

		if (GMT_parse_common_options(buf, &rasinfo[nfound].h.x_min, &rasinfo[nfound].h.x_max, &rasinfo[nfound].h.y_min, &rasinfo[nfound].h.y_max) ) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\t-R string conversion error.\n");
			return(0);
		}
		rasinfo[nfound].geo = (fabs (rasinfo[nfound].h.x_min) > 360.0 || fabs (rasinfo[nfound].h.x_max) > 360.0 || fabs (rasinfo[nfound].h.y_min) > 90.0 || fabs (rasinfo[nfound].h.y_max) > 90.0) ? FALSE : TRUE;
		/* Now find the -I string:  */
		i = j+1;
		while (i < length && (rasinfo[nfound].h.command[i] != '-') ) i++;
		if (i == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\t-I string conversion error.\n");
			return(0);
		}
		j = i+1;
		while (j < length && !(rasinfo[nfound].h.command[j] == ' ' ||  rasinfo[nfound].h.command[j] == '\t') ) j++;
		if (j == length || i+2 >= j) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\t-I string conversion error.\n");
			return(0);
		}
		i += 2;
		strncpy(buf, &rasinfo[nfound].h.command[i], (size_t)j-i);
		buf[j-i]='\0';
		if (GMT_getinc(buf, &rasinfo[nfound].h.x_inc, &rasinfo[nfound].h.y_inc) ) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\t-I string conversion error.\n");
			return(0);
		}

		/* Get P or G:  */
		i = j+1;
		while(i < length && !(rasinfo[nfound].h.command[i] == 'P' || rasinfo[nfound].h.command[i] == 'G') ) i++;
		if (i == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tP or G not found.\n");
			return(0);
		}
		rasinfo[nfound].h.node_offset = (rasinfo[nfound].h.command[i] == 'P') ? 1 : 0;

		/* Check if we have optional G (geographic) or C (Cartesian) that should override the auto test above */

		if (rasinfo[nfound].h.command[i+1] == 'G') {	/* Explicit geographic grid */
			rasinfo[nfound].geo = TRUE;
			i++;
		}
		else if (rasinfo[nfound].h.command[i+1] == 'C') {	/* Explicit Cartesian grid */
			rasinfo[nfound].geo = FALSE;
			i++;
		}

		stop_point = i + 1;

		/* Get type  */
		j = i + 1;
		while (j < length && (rasinfo[nfound].h.command[j] == ' ' || rasinfo[nfound].h.command[j] == '\t') ) j++;
		if (j == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tType conversion error.\n");
			return(0);
		}
		switch (rasinfo[nfound].h.command[j]) {
			case 'b':
			case 'c':
			case 'd':
			case 'i':
			case 'l':
			case 'u':
				rasinfo[nfound].type = rasinfo[nfound].h.command[j];
				break;
			default:
				fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
				fprintf(stderr,"\tInvalid type\n");
				return(0);
		}

		/* Get scale factor  */
		i = j + 1;
		while (i < length && (rasinfo[nfound].h.command[i] == ' ' || rasinfo[nfound].h.command[i] == '\t') ) i++;
		if (i == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tScale factor conversion error.\n");
			return(0);
		}
		j = i + 1;
		while (j < length && !(rasinfo[nfound].h.command[j] == ' ' || rasinfo[nfound].h.command[j] == '\t') ) j++;
		if (j == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tScale factor conversion error.\n");
			return(0);
		}
		strncpy(buf, &rasinfo[nfound].h.command[i], (size_t)j-i);
		buf[j-i] = '\0';
		if ( (sscanf(buf, "%lf", &rasinfo[nfound].h.z_scale_factor) ) != 1) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tScale factor conversion error.\n");
			return(0);
		}

		/* Get offset  */
		i = j+1;
		while (i < length && (rasinfo[nfound].h.command[i] == ' ' || rasinfo[nfound].h.command[i] == '\t') ) i++;
		if (i == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tOffset conversion error.\n");
			return(0);
		}
		j = i + 1;
		while (j < length && !(rasinfo[nfound].h.command[j] == ' ' || rasinfo[nfound].h.command[j] == '\t') ) j++;
		if (j == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tOffset conversion error.\n");
			return(0);
		}
		strncpy(buf, &rasinfo[nfound].h.command[i], (size_t)j-i);
		buf[j-i] = '\0';
		if ( (sscanf(buf, "%lf", &rasinfo[nfound].h.z_add_offset) ) != 1) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tOffset conversion error.\n");
			return(0);
		}

		/* Get NaNflag  */
		i = j+1;
		while (i < length && (rasinfo[nfound].h.command[i] == ' ' || rasinfo[nfound].h.command[i] == '\t') ) i++;
		if (i == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tNaN flag conversion error.\n");
			return(0);
		}
		j = i + 1;
		while (j < length && !(rasinfo[nfound].h.command[j] == ' ' || rasinfo[nfound].h.command[j] == '\t') ) j++;
		if (j == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tNaN flag conversion error.\n");
			return(0);
		}
		strncpy(buf, &rasinfo[nfound].h.command[i], (size_t)j-i);
		buf[j-i] = '\0';
		if (buf[0] == 'n' || buf[0] == 'N') {
			rasinfo[nfound].nanset = 0;
		}
		else {
			rasinfo[nfound].nanset = 1;
			if ( (sscanf(buf, "%" GMT_LL "d", &rasinfo[nfound].nanflag) ) != 1) {
				fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
				fprintf(stderr,"\tNaN flag conversion error.\n");
				return(0);
			}
		}


		/* Get filename:  */
		i = j+1;
		while (i < length && (rasinfo[nfound].h.command[i] == ' ' || rasinfo[nfound].h.command[i] == '\t') ) i++;
		if (i == length) {
			fprintf(stderr,"%s: error reading grdraster.info\n", GMT_program);
			fprintf(stderr,"\tFile name conversion error.\n");
			return(0);
		}
		j = i + 1;
		while (j < length && !(rasinfo[nfound].h.command[j] == ' ' || rasinfo[nfound].h.command[j] == '\t') ) j++;
		strncpy(buf, &rasinfo[nfound].h.command[i], (size_t)j-i);
		buf[j-i] = '\0';

#if WIN32
		if (buf[0] == '/' || buf[1] == ':') {
#else
		if (buf[0] == '/') {
#endif
			strcpy(rasinfo[nfound].h.remark, buf);
		}
		else {
			sprintf(rasinfo[nfound].h.remark, "%s%s", dir, buf);
		}

		/* Decode SWAP flag or SKIP command, if present  */

		i = j + 1;
		while (i < length && (rasinfo[nfound].h.command[i] == ' ' || rasinfo[nfound].h.command[i] == '\t') ) i++;
		if (i < length) {	/* Swap or skip flag set*/
			switch (rasinfo[nfound].h.command[i]) {
				case 'L':
				case 'l':	/* Little endian byte order */
					rasinfo[nfound].swap_me = (endian != 'L');	/* Must swap */
					break;
				case 'B':
				case 'b':	/* Big endian byte order */
					rasinfo[nfound].swap_me = (endian != 'B');	/* Must swap */
					break;
				case 'H':
				case 'h':	/* Give header size for skipping */
					rasinfo[nfound].skip = atoi (&rasinfo[nfound].h.command[i+1]);	/* Must skip header */
					break;
				default:
					fprintf (stderr,"%s: error reading grdraster.info\n", GMT_program);
					fprintf(stderr,"\tByte order or skip conversion error.\n");
					return (0);
			}
		}

		j = i + 1;
		while (j < length && (rasinfo[nfound].h.command[j] == ' ' || rasinfo[nfound].h.command[j] == '\t') ) j++;

		if (j < length) {	/* Skip of swap flag following the first skip or Swap flag */
			switch (rasinfo[nfound].h.command[j]) {
				case 'L':
				case 'l':	/* Little endian byte order */
					rasinfo[nfound].swap_me = (endian != 'L');	/* Must swap */
					break;
				case 'B':
				case 'b':	/* Big endian byte order */
					rasinfo[nfound].swap_me = (endian != 'B');	/* Must swap */
					break;
				case 'H':
				case 'h':	/* Give header size for skipping */
					rasinfo[nfound].skip = atoi (&rasinfo[nfound].h.command[j+1]);	/* Must skip header */
					break;
				default:
					fprintf (stderr,"%s: error reading grdraster.info\n", GMT_program);
					fprintf(stderr,"\tByte order or skip conversion error.\n");
					return (0);
			}
		}

		/* Get here when all is OK for this line:  */
		global_lon = 360.0 - (1 - rasinfo[nfound].h.node_offset)*rasinfo[nfound].h.x_inc;
		lon_tol = 0.01 * rasinfo[nfound].h.x_inc;
		global_lon -= lon_tol;	/* make sure we don't fail to find a truly global file  */
		if (rasinfo[nfound].geo && rasinfo[nfound].h.x_max - rasinfo[nfound].h.x_min >= global_lon) {
			rasinfo[nfound].nglobal = irint(360.0/rasinfo[nfound].h.x_inc);
		}
		else {
			rasinfo[nfound].nglobal = 0;
		}

		rasinfo[nfound].h.command[stop_point] = '\0';

		i = irint( (rasinfo[nfound].h.x_max - rasinfo[nfound].h.x_min)/rasinfo[nfound].h.x_inc);
		rasinfo[nfound].h.nx = (int)((rasinfo[nfound].h.node_offset) ? i : i + 1);
		j = irint( (rasinfo[nfound].h.y_max - rasinfo[nfound].h.y_min)/rasinfo[nfound].h.y_inc);
		rasinfo[nfound].h.ny = (int)((rasinfo[nfound].h.node_offset) ? j : j + 1);

		if ((ksize = get_byte_size (rasinfo[nfound].type)) == 0)
			expected_size = (off_t)(ceil (GMT_get_nm (rasinfo[nfound].h.nx, rasinfo[nfound].h.ny) * 0.125) + rasinfo[nfound].skip);
		else
			expected_size = (off_t)(GMT_get_nm (rasinfo[nfound].h.nx, rasinfo[nfound].h.ny) * ksize + rasinfo[nfound].skip);
		if (GMT_STAT (rasinfo[nfound].h.remark, &F)) {	/* Inquiry about file failed somehow */
			fprintf (stderr, "%s: Warning: Unable to stat file %s\n", GMT_program, rasinfo[nfound].h.remark);
		}
		else if (F.st_size != expected_size) {
			fprintf (stderr, "%s: Actual size of file %s [%ld] differs from expected [%ld]. Verify file and its grdraster.info details.\n", GMT_program, rasinfo[nfound].h.remark, (GMT_LONG)F.st_size, (GMT_LONG)expected_size);
			exit (EXIT_FAILURE);
		}
		nfound++;

		if ((size_t)nfound == n_alloc) {
			n_alloc <<= 1;
			rasinfo = (struct GRDRASTER_INFO *) GMT_memory ((void *)rasinfo, (size_t)n_alloc, sizeof (struct GRDRASTER_INFO), GMT_program);
		}
	}
	fclose (fp);
	if (nfound > 0) rasinfo = (struct GRDRASTER_INFO *) GMT_memory ((void *)rasinfo, (size_t)nfound, sizeof (struct GRDRASTER_INFO), GMT_program);

	*ras = rasinfo;
	return (nfound);
}

GMT_LONG get_byte_size (char type) {
	/* Return byte size of each item, or 0 if bits */
	int ksize;
	switch (type) {
		case 'b':
			ksize = 0;
			break;
		case 'c':
		case 'u':
			ksize = 1;
			break;
		case 'd':
		case 'h':
		case 'i':
			ksize = 2;
			break;
		case 'l':
			ksize = 4;
			break;
		default:
			fprintf (stderr, "%s: ERROR: Invalid data type [%c]\n", GMT_program, type);
			exit (EXIT_FAILURE);
			break;
	}
	return (ksize);
}

void *New_grdraster_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDRASTER_CTRL *C;

	C = (struct GRDRASTER_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDRASTER_CTRL), "New_grdraster_Ctrl");

	return ((void *)C);
}

void Free_grdraster_Ctrl (struct GRDRASTER_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);
	GMT_free ((void *)C);
}
