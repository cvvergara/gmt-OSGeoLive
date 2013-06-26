/* $Id: img2mercgrd.c 9923 2012-12-18 20:45:53Z pwessel $
 *
 * Copyright (c) 1991-2013 by P. Wessel and W. H. F. Smith
 * See LICENSE.TXT file for copying and redistribution conditions.
 *
 * img2mercgrd.c
 *
 * img2mercgrd reads an "img" format file (used by Sandwell and Smith
 * in their estimates of gravity and depth from satellite altimetry),
 * extracts a Region [with optional N by N averaging], and writes out
 * the data [or Track control] as a pixel-registered GMT "grd" file, 
 * preserving the Mercator projection (gmtdefaults ELLIPSOID = Sphere) 
 * inherent in the img file.
 *
 * img file is read from dir GMT_IMGDIR if this is defined as an
 * environment variable; else img file name is opened directly.
 *
 * The coordinates in the img file are determined by the latitude and
 * longitude span of the img file and the width of an img pixel in
 * minutes of longitude.  These may be changed from their default values
 * using the lower-case [-W<maxlon>] [-D<minlat>/<maxlat>] [-m<minutes>]
 * but must be set to match the coverage of the input img file.
 *
 * The user-supplied -R<w>/<e>/<s>/<n> will be adjusted by rounding
 * outward so that the actual range spanned by the file falls exactly on
 * the edges of the input pixels.  If the user uses the option to 
 * average the input pixels into N by N square block averages, then the
 * -R will be adjusted so that the range spans the block averages.
 *
 * The coordinates in the output file are in spherical mercator projected
 * units ranging from 0,0 in the lower left corner of the output to
 * xmax, ymax in the upper right corner, where xmax,ymax are the width
 * and height of the (spherical) Mercator map with -Rww/ee/ss/nn and -Jm1.
 * Here ww/ee/ss/nn are the outward-rounded range values.
 *
 * Further details on the coordinate system can be obtained from the
 * comments in gmt_imgsubs.c and gmt_imgsubs.h
 *
 * This is a complete rebuild (for v3.1) of the old program by this name.
 * New functionality added here is the averaging option [-N] and the
 * options to define -m, -W, -y, which permit handling very early versions
 * of these files and anticipate higher-resolutions to come down the road.
 * Also added is the option to look for the img file in GMT_IMGDIR
 *
 * Author:  Walter H F Smith
 * Date:    8 October, 1998
 * 
 **WHFS 27 Nov 1998 fixed bug so that -T0 gives correct results
 **  PW 18 Oct 1999 Use WORDS_BIGENDIAN macro set by configure.
 *   PW 12 Apr 2006 Replaced -x -y with -W -D
 *   PW 28 Nov 2006 Added -C for setting origin to main img origin (0,0)
 *
 */


#include "gmt_imgsubs.h"

struct IMG2MERCGRD_CTRL {
	char *file;	/* Input file name */
	struct C {	/* -C */
		GMT_LONG active;
	} C;
	struct D {	/* -D[<minlat>/<maxlat>] */
		GMT_LONG active;
		double min, max;
	} D;
	struct G {	/* -G<output_grdfile> */
		GMT_LONG active;
		char *file;
	} G;
	struct N {	/* -N<ave> */
		GMT_LONG active;
		int value;
	} N;
	struct S {	/* -S<scale> 	(use 0.1 for grav,curv,north,east files.) */
		GMT_LONG active;
		double value;
	} S;
	struct T {	/* -T<type> */
		GMT_LONG active;
		int value;
	} T;
	struct W {	/* -W<maxlon> */
		GMT_LONG active;
		double value;
	} W;
	struct m {	/* -m<minutes> */
		GMT_LONG active;
		double value;
	} m;
};

int main (int argc, char **argv)
{
	GMT_LONG	i, navgsq, error = FALSE, irarg = -1, tempint, n_files = 0;
	GMT_LONG	navg;	/* navg by navg pixels are averaged if navg > 1; else if navg == 1 do nothing */
	GMT_LONG	iout, jout, iinstart, iinstop, jinstart, jinstop, k, kk, ion, jin, jj, iin, ii, kstart;
	GMT_LONG nm;
	GMT_LONG ij;
	GMT_LONG	*ix = NULL;
	double	west, east, south, north, toplat, botlat, dx, rnavgsq, csum, dsum, left, bottom;
	struct GRD_HEADER h;
	struct GMT_IMG_COORD imgcoord;
	float	*a = NULL, empty_val;
	short int *row = NULL;
	char	infile[BUFSIZ];
	FILE	*fp = NULL;

	struct GMT_IMG_RANGE imgrange = { GMT_IMG_MAXLON, GMT_IMG_MINLAT, GMT_IMG_MAXLAT, GMT_IMG_MPIXEL };

	struct IMG2MERCGRD_CTRL *Ctrl;
	
	void *New_img2mercgrd_Ctrl (), Free_img2mercgrd_Ctrl (struct IMG2MERCGRD_CTRL *C);

	argc = (int)GMT_begin (argc, argv);
	
	Ctrl = (struct IMG2MERCGRD_CTRL *)New_img2mercgrd_Ctrl ();	/* Allocate and initialize a new control structure */

	GMT_grd_init (&h, argc, argv, FALSE);
	
	empty_val = GMT_f_NaN;

	west = east = south = north = 0.0;

	GMT_io.in_col_type[0] = GMT_IS_LON;	GMT_io.in_col_type[1] = GMT_IS_LAT;
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters:  */

				case 'R':
				case 'V':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					if (argv[i][1] == 'R') irarg = i;
					break;

				/* Supplemental parameters:  */

				case 'C':
					Ctrl->C.active = TRUE;
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					if (argv[i][2] && (sscanf(&argv[i][2], "%lf/%lf", &Ctrl->D.min, &Ctrl->D.max)) != 2) {
						error++;
						fprintf(stderr,"%s:  GMT SYNTAX ERROR -D: Failed to decode <minlat>/<maxlat>.\n", GMT_program);
					}
					else {
						Ctrl->D.min = GMT_IMG_MINLAT_80;
						Ctrl->D.max = GMT_IMG_MAXLAT_80;
					}
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					if ((sscanf(&argv[i][2], "%d", &Ctrl->N.value)) != 1 || Ctrl->N.value < 1) {
						error++;
						fprintf(stderr,"%s:  GMT SYNTAX ERROR -N requires an integer > 1.\n", GMT_program);
					}
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					if ((sscanf(&argv[i][2], "%lf", &Ctrl->S.value)) != 1) {
						error++;
						fprintf(stderr,"%s:  GMT SYNTAX ERROR -S requires a scale factor.\n", GMT_program);
					}
					break;				
				case 'T':
					Ctrl->T.active = TRUE;
					if ((sscanf(&argv[i][2], "%d", &Ctrl->T.value)) != 1) {
						error++;
						fprintf(stderr,"%s:  GMT SYNTAX ERROR -T requires an output type 0-3.\n", GMT_program);
					}
					break;
				case 'W':
					Ctrl->W.active = TRUE;
					if ((sscanf(&argv[i][2], "%lf", &Ctrl->W.value)) != 1) {
						error++;
						fprintf(stderr,"%s:  GMT SYNTAX ERROR -W requires a longitude >= 360.0.\n", GMT_program);
					}
					break;
				case 'm':
					Ctrl->m.active = TRUE;
					if ((sscanf(&argv[i][2], "%lf", &Ctrl->m.value)) != 1) {
						error++;
						fprintf(stderr,"%s:  GMT SYNTAX ERROR -m requires a positive value.\n", GMT_program);
					}
					break;
				default:
					error++;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else if (n_files == 0) {
			Ctrl->file = strdup (argv[i]);
			n_files++;
		}
		else {
			fprintf(stderr,"%s:  GMT SYNTAX ERROR: More than one world image file name given.\n", GMT_program);
			error++;
		}
	}

	if (GMT_give_synopsis_and_exit || argc == 1) {
		fprintf(stderr,"img2mercgrd %s - Extract a grid file from an img file, preserving Mercator projection.\n\n", GMT_VERSION);
		fprintf(stderr,"usage:  img2mercgrd <world_image_filename> %s -G<grdfile> -T<type>\n", GMT_Rgeo_OPT);
		fprintf(stderr,"\t\t[-C] [-D[<minlat>/<maxlat>]] [-N<navg>] [-S<scale>] [-V] [-W<maxlon>] [-m<minutes>]\n\n");

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
		
		fprintf (stderr, "\tREQUIRED ARGUMENTS:\n");
		fprintf(stderr,"\t<world_image_filename> gives location of img file.\n");
		fprintf(stderr,"\t-G sets filename for the output .grd format file.\n");
		fprintf(stderr,"\t-R specifies the region in degrees or degrees:minutes.\n");
		fprintf(stderr,"\t-T selects output type:\n");
		fprintf(stderr,"\t\t-T0 for obsolete img files w/ no constraint code, gets data.\n");
		fprintf(stderr,"\t\t-T1 for new img file w/ constraints coded, gets data at all points.\n");
		fprintf(stderr,"\t\t-T2 for new img file w/ constraints coded, gets data only at constrained points, NaN elsewhere.\n");
		fprintf(stderr,"\t\t-T3 for new img file w/ constraints coded, gets 1 at constraints, 0 elsewhere.\n\n");
		
		fprintf (stderr, "\tOPTIONAL ARGUMENTS:\n");
		fprintf(stderr,"\t-C Refer Mercator coordinates to img source origin [Default sets lower left to 0,0].\n");
		fprintf(stderr,"\t-N<navg> will ouput averages of input in navg by navg squares.  [no averaging]\n");
		fprintf(stderr,"\t-S<scale> will multiply img integer values by scale for output [1]\n");
		GMT_explain_option ('V');
		
		fprintf (stderr, "\n\tOPTIONAL ARGUMENTS WHICH DEPEND ON THE img FILE VERSION:\n");
		fprintf (stderr, "\t-m<minutes> input img pixels are <minutes> minutes of longitude wide. [2.0]\n");
		fprintf (stderr, "\t-D[<minlat>/<maxlat>] input img file bottom and top latitudes. [%.3f/%.3f]\n", GMT_IMG_MINLAT, GMT_IMG_MAXLAT);
		fprintf (stderr, "\t   If no latitudes are given it is taken to mean %.3f/%.3f\n", GMT_IMG_MINLAT_80, GMT_IMG_MAXLAT_80);
		fprintf (stderr, "\t-W<maxlon> input img file runs from 0 to <maxlon> longitude. [360.0]\n");
		exit(EXIT_FAILURE);
	}

	if (irarg < 0 || west >= east || south >= north) {
		fprintf(stderr, "%s: GMT SYNTAX ERROR:  Must specify -R with west < east and south < north.\n", GMT_program);
		error++;
	}
	if (!Ctrl->G.active || Ctrl->G.file == NULL) {
		fprintf(stderr, "%s: GMT SYNTAX ERROR:  Must specify output grid file name with -G.\n", GMT_program);
		error++;
	}
	if (Ctrl->file == NULL) {
		fprintf(stderr, "%s: GMT SYNTAX ERROR:  Must specify input imgfile name.\n", GMT_program);
		error++;
	}
	if (Ctrl->D.active && (Ctrl->D.min <= -90 || Ctrl->D.max >= 90.0 || Ctrl->D.max <= Ctrl->D.min)) {
		fprintf(stderr, "%s: GMT SYNTAX ERROR:  Min/max latitudes are invalid\n", GMT_program);
		error++;
	}
	if (Ctrl->T.value < 0 || Ctrl->T.value > 3) {
		fprintf(stderr, "%s: GMT SYNTAX ERROR:  Must specify output type in the range 0-3 with -T.\n", GMT_program);
		error++;
	}
	if (Ctrl->W.active && Ctrl->W.value < 360.0) {
		fprintf(stderr, "%s: GMT SYNTAX ERROR:  Requires a maximum longitude >= 360.0 with -W.\n", GMT_program);
		error++;
	}
	if (Ctrl->m.active && Ctrl->m.value <= 0.0) {
		fprintf(stderr, "%s: GMT SYNTAX ERROR:  Requires a positive value with -m.\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	if (Ctrl->W.active) imgrange.maxlon = Ctrl->W.value;
	if (Ctrl->m.active) imgrange.mpixel = Ctrl->m.value;
	if (Ctrl->D.active) {
		imgrange.minlat = Ctrl->D.min;
		imgrange.maxlat = Ctrl->D.max;
	}
	
	if (GMT_img_setup_coord ( &imgrange, &imgcoord) ) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Error in img coordinate specification [-m -W or -D]\n", GMT_program);
		error++;
	}
	else if (Ctrl->N.value && (imgcoord.nx360%Ctrl->N.value != 0 || imgcoord.nyrow%Ctrl->N.value != 0) ) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Bad choice of navg in -N.  Must divide %ld and %ld\n", GMT_program,
			imgcoord.nx360, imgcoord.nyrow);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	GMT_getdatapath (Ctrl->file, infile);
	if ((fp = fopen(infile, "rb")) == NULL) {
		fprintf(stderr, "%s: GMT SYNTAX ERROR:  Cannot open %s for binary read.\n", GMT_program, infile);
		error++;
	}
	if (error) exit (EXIT_FAILURE);
	
	GMT_io.out_col_type[0] = GMT_io.out_col_type[1] = GMT_IS_FLOAT;		/* Since output is no longer lon/lat */
	navg = Ctrl->N.value;
	
	/* Expected edges of input image based on coordinate initialization (might not exactly match user spec):  */
	toplat = GMT_img_ypix_to_lat (0.0, &imgcoord);
	botlat = GMT_img_ypix_to_lat ( (double)imgcoord.nyrow, &imgcoord);
	dx = 1.0 / ( (double)imgcoord.nx360 / 360.0);
	if (gmtdefs.verbose)
		fprintf (stderr, "%s expects %s to be %ld by %ld pixels spanning 0/%5.1f/%.8g/%.8g\n",
			GMT_program, infile, imgcoord.nxcol, imgcoord.nyrow, dx*imgcoord.nxcol, botlat, toplat);

	if (toplat < north) {
		fprintf(stderr, "%s:  WARNING:  Your top latitude (%.12g) lies outside top latitude of input (%.12g) - now truncated.\n",
			GMT_program, north, toplat);
		north = toplat - GMT_CONV_LIMIT;	/* To ensure proper round-off in calculating ny */
	}
	if (botlat > south) {
		fprintf(stderr, "%s:  WARNING:  Your bottom latitude (%.12g) lies outside bottom latitude of input (%.12g) - now truncated.\n",
			GMT_program, south, botlat);
		south = botlat + GMT_CONV_LIMIT;	/* To ensure proper round-off in calculating ny */
	}
	
	/* Re-adjust user's -R so that it falls on pixel coordinate boundaries:  */
	
	jinstart = navg * (int)floor (GMT_img_lat_to_ypix (north, &imgcoord) / navg);
	jinstop  = navg * (int)ceil  (GMT_img_lat_to_ypix (south, &imgcoord) / navg);
	/* jinstart <= jinputrow < jinstop  */
	h.ny = (int)((jinstop - jinstart) / navg);
	north = GMT_img_ypix_to_lat ((double)jinstart, &imgcoord);
	south = GMT_img_ypix_to_lat ((double)jinstop,  &imgcoord);

	iinstart = navg * (GMT_LONG)floor (west/(dx*navg));
	iinstop  = navg * (GMT_LONG)ceil  (east/(dx*navg));
	/* iinstart <= ipixelcol < iinstop, but modulo all with imgcoord.nx360  */
	/* Reset left and right edges of user area:  */
	west = iinstart * dx;
	east = iinstop  * dx;
	h.nx = (int)((iinstop - iinstart) / navg);

	if (gmtdefs.verbose) {
		fprintf(stderr, "%s:  To fit [averaged] input, your %s is adjusted to -R%.12g/%.12g/%.12g/%.12g\n",
			GMT_program, argv[irarg], west, east, south, north);
		fprintf (stderr, "%s:  The output will be %d by %d pixels.\n", GMT_program, h.nx, h.ny);
	}

	/* Set iinstart so that it is non-negative, for use to index pixels.  */
	while (iinstart < 0) iinstart += imgcoord.nx360;
	
	/* Set navgsq, rnavgsq, for the averaging:  */
	navgsq = navg * navg;
	rnavgsq = 1.0 / navgsq;

	/* Set up header with Mercatorized dimensions assuming -Jm1  */
	if (Ctrl->C.active) {
		GMT_LONG equator;
		equator = irint (GMT_img_lat_to_ypix (0.0, &imgcoord));
		h.x_min = iinstart * navg * dx;
		h.x_max = h.x_min + h.nx * navg * dx;
		h.y_max = (imgcoord.nyrow - jinstart - equator) * navg * dx;
		h.y_min = h.y_max - h.ny * navg * dx;
		left = bottom = 0.0;
		if (h.x_max > 360.0) {
			h.x_max -= 360.0;
			h.x_min -= 360.0;
		}
	}
	else {
		h.x_min = 0.0;
		h.x_max = h.nx * navg * dx;
		h.y_min = 0.0;
		h.y_max = h.ny * navg * dx;
		left = west;
		bottom = south;
	}
	sprintf (h.x_units, "Spherical Mercator projected Longitude, -Jm1, length from %.12g", left);
	sprintf (h.y_units, "Spherical Mercator projected Latitude, -Jm1, length from %.12g", bottom);
	h.x_inc = navg * dx;
	h.y_inc = navg * dx;
	h.node_offset = 1;
	h.z_scale_factor = 1.0;
	h.z_add_offset = 0.0;
	if (Ctrl->T.value < 3)
		strcpy (h.z_units, "meters, mGal, Eotvos, or micro-radians, depending on img file and -S.");
	else
		strcpy (h.z_units, "T/F, one or more constraints fell in this pixel.");
	strcpy (h.title, "Data from Altimetry");
	sprintf (h.remark, "Spherical Mercator Projected with -Jm1 -R%.12g/%.12g/%.12g/%.12g", west, east, south, north);
	h.z_min = DBL_MAX;
	h.z_max = -DBL_MAX;

	GMT_err_fail (GMT_grd_RI_verify (&h, 1), Ctrl->G.file);

	/* Now malloc some space for float grd array, integer pixel index, and short integer data buffer.  */

	row = (short int *) GMT_memory (VNULL, (size_t)(navg * imgcoord.nxcol), sizeof (short int), GMT_program);
	ix = (GMT_LONG *) GMT_memory (VNULL, (size_t)(navgsq * h.nx), sizeof (GMT_LONG), GMT_program);
	nm = GMT_get_nm (h.nx, h.ny);
	a = (float *) GMT_memory (VNULL, nm, sizeof (float), GMT_program);

	/* Load ix with the index to the correct column, for each output desired.  This helps for Greenwich, 
		also faster averaging of the file, etc.  Note for averaging each n by n block is looped in turn. */
	
	if (navg > 1) {
		k = 0;
		for (iout = 0; iout < h.nx; iout++) {
			ion = iout * navg;
			for (jin = 0; jin < navg; jin++) {
				jj = jin * imgcoord.nxcol;
				for (iin = 0; iin < navg; iin++) {
					ii = (iin + iinstart + ion) % imgcoord.nx360;
					ix[k] = ii + jj;
					k++;
				}
			}
		}
	}
	else {
		for (iout = 0; iout < h.nx; iout++) {
			ix[iout] = (iout + iinstart) % imgcoord.nx360;
		}
	}


	/* Now before beginning data loop, fseek if needed.  */
	if (jinstart > 0 && jinstart < imgcoord.nyrow) {
		fseek (fp, (long)(2 * imgcoord.nxcol * jinstart), SEEK_SET);
	}
	
	/* Now loop over output points, reading and handling data as needed */

	for (ij = 0, jout = 0; jout < h.ny; jout++) {
		jin = jinstart + navg * jout;
		if (jin < 0 || jin >= imgcoord.nyrow) {
			for (iout = 0; iout < h.nx; iout++, ij++) {
				a[ij] = empty_val;
			}
		}
		else {
			if ( (fread((void *)row, sizeof (short int), (size_t)(navg * imgcoord.nxcol), fp) ) != (size_t)(navg * imgcoord.nxcol) ) {
				fprintf(stderr,"%s:  ERROR:  Read failure at jin = %ld.\n", GMT_program, jin);
				exit (EXIT_FAILURE);
			}

#if defined(_WIN32) || !defined(WORDS_BIGENDIAN)
			for (iout = 0; iout < navg * imgcoord.nxcol; iout++) row[iout] = GMT_swab2 (row[iout]);
#endif

			for (iout = 0, kstart = 0; iout < h.nx; iout++, ij++, kstart += navgsq) {
				if (navg) {
					csum = 0.0;
					dsum = 0.0;
					for (k = 0, kk = kstart; k < navgsq; k++, kk++) {
						tempint = (int)row[ix[kk]];
						if (Ctrl->T.value) {
							if ( ( (GMT_abs(tempint))%2) != 0) {
								csum += 1.0;
								tempint--;
							}
						}
						dsum += (double) tempint;
					}
					csum *= rnavgsq;
					dsum *= rnavgsq;
				}
				else {
					tempint = (int)row[ix[iout]];
					if (Ctrl->T.value && GMT_abs(tempint)%2 != 0) {
						csum = 1.0;
						tempint--;
					}
					else
						csum = 0.0;
					dsum = (double) tempint;
				}
				
				if (Ctrl->S.active) dsum *= Ctrl->S.value;
				
				switch (Ctrl->T.value) {
					case 0:
					case 1:
						a[ij] = (float) dsum;
						break;
					case 2:
						a[ij] = (float)((csum >= 0.5) ? dsum : GMT_f_NaN);
						break;
					case 3:
						a[ij] = (float)csum;
						break;
				}
				

				if (Ctrl->T.value != 2 || csum >= 0.5) {
					if (h.z_min > a[ij]) h.z_min = a[ij];
					if (h.z_max < a[ij]) h.z_max = a[ij];
				}
			}
		}
	}
	
	fclose (fp);
	if (gmtdefs.verbose)
		fprintf(stderr,"Created %d by %d Mercatorized grid file.  Min, Max values are %.8g  %.8g\n", h.nx, h.ny, h.z_min, h.z_max);

	GMT_err_fail (GMT_write_grd (Ctrl->G.file, &h, a, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);

	GMT_free((void *)a);
	GMT_free((void *)ix);
	GMT_free((void *)row);
	
	Free_img2mercgrd_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_img2mercgrd_Ctrl () {	/* Allocate and initialize a new control structure */
	struct IMG2MERCGRD_CTRL *C;
	
	C = (struct IMG2MERCGRD_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct IMG2MERCGRD_CTRL), "New_img2mercgrd_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->D.min = GMT_IMG_MINLAT;
	C->D.max = GMT_IMG_MAXLAT;
	C->N.value = 1;		/* No averaging */
	C->T.value = -1;	/* Force user to change this */
	C->m.value = GMT_IMG_MPIXEL;
	
	return ((void *)C);
}

void Free_img2mercgrd_Ctrl (struct IMG2MERCGRD_CTRL *C) {	/* Deallocate control structure */
	if (C->file) free ((void *)C->file);	
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}
