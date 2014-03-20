/*--------------------------------------------------------------------
 *    $Id: psvelo.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *    Copyright (c) 1996-2014 by G. Patau
 *    Distributed under the GNU Public Licence
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*

psvelo will read <x,y> pairs (or <lon,lat>) from inputfile and
plot symbols on a map. Velocity ellipses, strain
crosses, or strain wedges, may be specified, some of which require
additional columns of data.  Only one symbol may be plotted at a time.
PostScript code is written to stdout.


 Author:     Kurt Feigl
 Date:	      7 July 1998
 Version:     4
 Roots:       based on psxy.c
 Adapted to version 3.3 by Genevieve Patau (25 June 1999)
 Last modified : 18 February 2000

 */

#include "gmt.h"	/* to have gmt environment */
#include "pslib.h"	/* to have pslib environment */
#include "utilvelo.h"
#include "utilstrain.h"

#define EPSIL 0.0001

#define CINE 1
#define ANISO 2
#define WEDGE 3
#define CROSS 4


#define veclen(x, y) sqrt((x) * (x) + (y) * (y))

     /* parameters for writing text */
#define ANGLE 0.0
#define FORM 0
#define TIRET_WIDTH 3
#define POINTSIZE 0.005

int main (int argc, char **argv)
{
	GMT_LONG 	i, symbol = 0, n, ix = 0, iy = 1, n_files = 0, fno;
	GMT_LONG	n_args;
	GMT_LONG     n_rec = 0;

	GMT_LONG	outline = FALSE;
	GMT_LONG	error = FALSE, nofile = TRUE, polygon = FALSE;
	GMT_LONG shade_uncert = FALSE;
	GMT_LONG rescale_sigma = FALSE;
	GMT_LONG done, no_size_needed, greenwich;
	GMT_LONG read_ellipse = FALSE, read_rotated_ellipse = FALSE;
	GMT_LONG des_ellipse = TRUE, des_arrow = TRUE;
	GMT_LONG read_anisotropy = FALSE;
	GMT_LONG read_wedge = FALSE;
	GMT_LONG read_cross = FALSE;
	GMT_LONG old_GMT_world_map, skip_if_outside = TRUE;

	double xy[2], west = 0.0, east = 0.0, south = 0.0, north = 0.0;
	double plot_x, plot_y, scale = 0.0;
	double vxy[2], plot_vx, plot_vy;
	double eps1 = 0.0, eps2 = 0.0, spin = 0.0, spinsig = 0.0, theta = 0.0;
	double direction, small_axis, great_axis; 
	double sigma_x, sigma_y, corr_xy; 
	double confidence = 0., sigma_scale = 1.0, conrad = 1.0;
	double v_width = 0.01, h_length = 0.12, h_width = 0.03, vector_shape = 0.4;
	double wedge_amp = 1.e7;
	double t11 = 1.0, t12 = 0.0, t21 = 0.0, t22 = 1.0;
	double hl,hw,vw, fontsize = 0.0;

	char *station_name, *not_used = NULL;
	char txt_a[GMT_TEXT_LEN], txt_b[GMT_TEXT_LEN], txt_c[GMT_TEXT_LEN];
	char line[BUFSIZ], symbol_type, col[12][GMT_TEXT_LEN];

	FILE *fp = NULL;


	struct GMT_PEN pen;
	struct GMT_FILL fill, efill; /* efill is for uncertainty wedge */

	GMT_LONG justify;

	argc = (int)GMT_begin (argc, argv);

	GMT_init_pen (&pen, GMT_PENWIDTH);
	GMT_init_fill (&fill, 0, 0, 0);
	GMT_init_fill (&efill, 255, 255, 255);



	/* Check and interpret the command line arguments */

	for (i = 1; !error && i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {

				/* Common parameters */

				case 'B':
				case 'H':
				case 'J':
				case 'K':
				case 'O':
				case 'P':
				case 'R':
				case 'U':
				case 'V':
				case 'X':
				case 'x':
				case 'Y':
				case 'y':
				case 'c':
				case ':':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */

                                case 'A':       /* Change size of arrow head */
                                        strcpy(txt_a, &argv[i][2]);
                                        strcpy(txt_b, strchr(txt_a+1, '/')+1);
                                        strcpy(txt_c, strchr(txt_b+1, '/')+1);
                                        n=0; while (txt_a[n] && txt_a[n] != '/') n++; txt_a[n]=0;
                                        n=0; while (txt_b[n] && txt_b[n] != '/') n++; txt_b[n]=0;
                                        v_width = GMT_convert_units (txt_a, GMT_INCH);
                                        h_length = GMT_convert_units (txt_b, GMT_INCH);
                                        h_width = GMT_convert_units (txt_c, GMT_INCH);
                                        break; 
                                case 'D':       /* Rescale Sigmas */
                                        rescale_sigma = TRUE;
					sscanf (&argv[i][2], "%lf",&sigma_scale);
                                        break;
				case 'E':		/* Set color for error ellipse  */
					GMT_getfill (&argv[i][2], &efill);
					shade_uncert = TRUE;
					break;
				case 'G':		/* Set Gray shade for polygon */
					GMT_getfill (&argv[i][2], &fill);
					polygon = TRUE;
					break;
				case 'L':		/* Draw the outline */
					outline = TRUE;
					break;
				case 'N':		/* Do not skip points outside border */
					skip_if_outside = FALSE;
					break;
				case 'S':		/* Get symbol [and size] */
                                        symbol_type = argv[i][2];
                                        if(symbol_type == 'e' || symbol_type == 'r') {
                                            strcpy(txt_a, &argv[i][3]);
                                            n=0; while (txt_a[n] && txt_a[n] != '/') n++; txt_a[n]=0;
                                            scale = GMT_convert_units (txt_a, GMT_INCH);
                                            sscanf(strchr(&argv[i][3],'/')+1, "%lf/%lf", &confidence, &fontsize);
					    /* confidence scaling */
					    conrad = sqrt( -2.0 * log(1.0 - confidence)); 
					  }
                                        if(symbol_type == 'n' || symbol_type == 'x' )
                                            scale = GMT_convert_units (&argv[i][3], GMT_INCH);
                                        if(symbol_type == 'w' && strlen(argv[i]) > 3) {
                                            strcpy(txt_a, &argv[i][3]);
                                            n=0; while (txt_a[n] && txt_a[n] != '/') n++; txt_a[n]=0;
                                            scale = GMT_convert_units (txt_a, GMT_INCH);
                                            sscanf(strchr(&argv[i][3],'/')+1, "%lf", &wedge_amp);
                                        }
					switch (symbol_type) {
                                                case 'e':
                                                        symbol = CINE;
                                                        read_ellipse = TRUE;
                                                        break;
                                                case 'r':
                                                        symbol = CINE;
                                                        read_rotated_ellipse = TRUE;
                                                        break;
                                                case 'n':
                                                        symbol = ANISO;
                                                        read_anisotropy = TRUE;
                                                        break;
                                                case 'w':
                                                        symbol = WEDGE;
                                                        read_wedge = TRUE;
                                                        break;
                                                case 'x':
                                                        symbol = CROSS;
                                                        read_cross = TRUE;
                                                        break;
						default:
							error = TRUE;
							break;
					}
					break;
				case 'W':		/* Set line attributes */
					GMT_getpen (&argv[i][2], &pen);
					break;

				/* Illegal options */

				default:		/* Options not recognized */
					error = TRUE;
					break;
			}
		}
		else
			n_files++;
	}

	/* Check that the options selected are mutually consistent */

	no_size_needed = (read_ellipse || read_rotated_ellipse || read_anisotropy || read_cross || read_wedge );
	error += GMT_check_rgb (pen.rgb)
	      + GMT_check_rgb (fill.rgb)
	      + GMT_check_rgb (efill.rgb)
	      + GMT_check_rgb (gmtdefs.basemap_frame_rgb);
        /* Only one allowed */
	if ((read_ellipse + read_rotated_ellipse + read_anisotropy + read_cross + read_wedge ) > 1) error = TRUE;
	if (!no_size_needed && (symbol > 1 && scale <= 0.0)) error = TRUE;
        if (rescale_sigma && ! (read_ellipse || read_wedge)) error = TRUE;

	if (argc == 1 || GMT_give_synopsis_and_exit || error) {	/* Display usage */
		fprintf (stderr,"psvelo %s - Plot symbols on maps\n\n", GMT_VERSION);
		fprintf (stderr,"usage: psvelo <infiles> %s %s [%s] [-A<awidth>/<alength>/<hwidth>]\n", GMT_J_OPT, GMT_Rgeo_OPT, GMT_B_OPT);
		fprintf (stderr, "	[-G<fill>] [%s] [-K] [-L] [-N] [-O]\n", GMT_H_OPT);
		fprintf (stderr, "	[-P] [-S<symbol><scale><fontsize>] [%s] [-V] [-W<pen>] [%s]\n", GMT_U_OPT, GMT_c_OPT);
		fprintf (stderr, "	[%s] [%s] [%s]\n\n", GMT_X_OPT, GMT_Y_OPT, GMT_t_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "	<infiles> is one or more files.  If no, read standard input\n");
		GMT_explain_option ('j');
		GMT_explain_option ('R');
		GMT_explain_option ('b');
                fprintf (stderr, "      -A Change the size of arrow head; specify arrow_width, head_length, head_width;");
                fprintf (stderr, "         default is 0.01/0.12/0.03;");
                fprintf (stderr, "      -D Multiply uncertainties by sigscale. (Se and Sw only) \n");
		fprintf (stderr, "      -E Set color used for uncertainty wedges in -Sw option.\n");
		fprintf (stderr, "	-G Specify color (for symbols/polygons) or pattern (for polygons). fill can be either\n");
		fprintf (stderr, "	   1) <r/g/b> (each 0-255) for color or <gray> (0-255) for gray-shade [0].\n");
		fprintf (stderr, "	   2) p[or P]<iconsize>/<pattern> for predefined patterns (0-31).\n");
		GMT_explain_option ('H');
		GMT_explain_option ('K');
		fprintf (stderr, "	-L draw line or symbol outline using the current pen (see -W).\n");
		fprintf (stderr, "	-N Do Not skip/clip symbols that fall outside map border [Default will ignore those outside]\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "	-S to select symbol type and scale. Choose between\n");
                fprintf (stderr, "         (e) Velocity ellipses: in X,Y,Vx,Vy,SigX,SigY,CorXY,name format \n");
                fprintf (stderr, "         (r) Velocity ellipses: in X,Y,Vx,Vy,a,b,theta,name format \n");
                fprintf (stderr, "         (n) Anisotropy : in X,Y,Vx,Vy \n");
                fprintf (stderr, "         (w) Rotational wedges: in X,Y,Spin,Spinsig \n");
                fprintf (stderr, "         (x) Strain crosses : in X,Y,Eps1,Eps2,Theta \n");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		fprintf (stderr, "	-W sets pen attributes [width = %gp, color = (%d/%d/%d), texture = solid line].\n", 
			pen.width, pen.rgb[0], pen.rgb[1], pen.rgb[2]);
		GMT_explain_option ('X');
		GMT_explain_option ('c');
		GMT_explain_option (':');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;
	n_args = (argc > 1) ? argc : 2;

	greenwich = (project_info.w < 0.0 || project_info.e <= 0.0);

	GMT_err_fail (GMT_map_setup (west, east, south, north), "");

	GMT_plotinit (argc, argv);

	GMT_setpen (&pen);
	ps_setfont (gmtdefs.annot_font[0]);
	if (shade_uncert) outline = TRUE;

	if (skip_if_outside) GMT_map_clip_on (GMT_no_rgb, 3);

	old_GMT_world_map = GMT_world_map;

	ix = (gmtdefs.xy_toggle[0]);	iy = 1 - ix;

	done = FALSE;

	station_name = (char *) GMT_memory (VNULL, (size_t)64, sizeof (char), GMT_program);

	for (fno = 1; !done && fno < n_args; fno++) {	/* Loop over all input files */
		if (!nofile && argv[fno][0] == '-') continue;
		if (nofile) {
			fp = GMT_stdin;
			done = TRUE;
		}
		else if ((fp = GMT_fopen (argv[fno], "r")) == NULL) {
			fprintf (stderr, "psvelo: Cannot open file %s\n", argv[fno]);
			continue;
		}

		if (!nofile && gmtdefs.verbose) {
			fprintf (stderr, "psvelo: Working on file %s\n", argv[fno]);
			sprintf (line, "File: %s", argv[fno]);
			ps_comment (line);
		}
		if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) not_used = GMT_fgets (line, 512, fp);

		if (gmtdefs.verbose && (read_ellipse || read_rotated_ellipse))	 
		  fprintf (stderr, "psvelo: 2-D confidence interval and scaling factor %f %f\n",confidence, conrad);

		while (GMT_fgets (line, 512, fp)) {
                        n_rec++;
			memset ((void *)col, 0, 12 * GMT_TEXT_LEN * sizeof (char));
			memset ((void *)station_name, 0, 64 * sizeof (char));
                       if (read_ellipse || read_rotated_ellipse) {
			  sscanf (line, "%s %s %s %s %s %s %s %[^\n]\n", 
				  col[0], col[1], col[2], col[3], col[4], col[5], col[6], station_name);
			  if (strlen(station_name) <= 0) sprintf(station_name,"\n");
			}
                        else {
                          sscanf (line, "%s %s %s %s %s %s %s %s %s %[^\n]\n",
				  col[0], col[1], col[2], col[3], col[4], col[5], col[6], col[7], col[8], station_name);
			  if (strlen(station_name) <= 0) sprintf(station_name,"\n");
			}

			xy[ix] = atof (col[0]);
			xy[iy] = atof (col[1]);

                        if(read_ellipse) {
                                    vxy[ix] = atof(col[2]);
                                    vxy[iy] = atof(col[3]);
                                    sigma_x = atof(col[4]);
                                    sigma_y = atof(col[5]);
                                    corr_xy = atof(col[6]);
                                    /* rescale uncertainties if necessary */
                                    if (rescale_sigma) {
				      sigma_x = sigma_scale * sigma_x;
				      sigma_y = sigma_scale * sigma_y;
				         } 
                                    if(fabs(sigma_x) < EPSIL && fabs(sigma_y) < EPSIL)
                                        des_ellipse = FALSE;
                                    else { 
                                        des_ellipse = TRUE;
                                        ellipse_convert(sigma_x, sigma_y, corr_xy, conrad, &small_axis, &great_axis, &direction);

                                    /* convert to degrees */
                                       direction = direction * R2D;
                                    }
                        }
                        else if(read_rotated_ellipse) {
                                    vxy[ix] = atof(col[2]);
                                    vxy[iy] = atof(col[3]);
                                    great_axis = conrad*atof(col[4]);
                                    small_axis = conrad*atof(col[5]);
                                    direction = atof(col[6]);
                                    if(fabs(great_axis) < EPSIL && fabs(small_axis) < EPSIL)
                                        des_ellipse = FALSE;
                                    else
                                        des_ellipse = TRUE;
                        }
                        else if(read_anisotropy) {
                                    vxy[ix] = atof(col[2]);
                                    vxy[iy] = atof(col[3]);
                        }
                        else if(read_cross) {
                                    eps1 = atof(col[2]);
                                    eps2 = atof(col[3]);
                                    theta = atof(col[4]);
			}
                        else if(read_wedge) {
                                    spin = atof(col[2]);
                                    spinsig = atof(col[3]);
                                    if (rescale_sigma) spinsig = spinsig * sigma_scale;
                        }

			if (skip_if_outside) {
				GMT_map_outside (xy[0], xy[1]);
				if (GMT_abs (GMT_x_status_new) > 1 || GMT_abs (GMT_y_status_new) > 1) continue;
			}

			GMT_geo_to_xy (xy[0], xy[1], &plot_x, &plot_y);


			switch (symbol) {
                                case CINE:
                                        des_arrow = veclen((vxy[0]), (vxy[1])) < 1.e-8 ? FALSE : TRUE; 
                                        trace_arrow(xy[0], xy[1], vxy[0], vxy[1], scale, &plot_x, &plot_y, &plot_vx, &plot_vy);
                                        get_trans (xy[0], xy[1], &t11, &t12, &t21, &t22);
                                        if(des_ellipse) {
                                            if (shade_uncert) 
                                               paint_ellipse (plot_vx, plot_vy, direction, great_axis, small_axis, scale,
                                               t11,t12,t21,t22, shade_uncert, efill.rgb, outline);
                                            else
                                               paint_ellipse (plot_vx, plot_vy, direction, great_axis, small_axis, scale,
                                               t11,t12,t21,t22, shade_uncert, fill.rgb, outline);
					}
                                        if(des_arrow) {
					  /* verify that arrow is not ridiculously small */
					  if (veclen (plot_x-plot_vx,plot_y-plot_vy) <= 1.5 * h_length) {
					    hl = veclen (plot_x-plot_vx,plot_y-plot_vy) * 0.6;
					    hw = hl * h_width/h_length;
					    vw = hl * v_width/h_length; 
					    if (vw < 2./(double)gmtdefs.dpi) vw = 2./(double)gmtdefs.dpi;
					  }
					  else {
					    hw = h_width;
					    hl = h_length;
					    vw = v_width;
					  }
					  ps_vector(plot_x, plot_y, plot_vx, plot_vy, vw, hl, hw, vector_shape, 
						    fill.rgb, outline);
					  justify = plot_vx - plot_x > 0. ? 7 : 5;
                                          if(fontsize > 0.0 && strlen(station_name) > 0) {
						ps_text(plot_x + (6 - justify) / 25.4 , plot_y, fontsize,
					     	station_name, ANGLE, justify,FORM); /* 1 inch = 2.54 cm */
                                          }
                                        }
                                        else {
                                            ps_circle(plot_x, plot_y, POINTSIZE, fill.rgb, 1);
                                            justify = 10;
                                            if(fontsize > 0.0 && strlen(station_name) > 0) {
                                                 ps_text(plot_x, plot_y - 1. / 25.4, fontsize, station_name, ANGLE, justify, FORM);
                                            }
                                            /*  1 inch = 2.54 cm */
                                        }
                                        i = 0;
                                        while(col[7][i] != '\0') {
                                            col[7][i] = ' ';
                                            i++;
                                        }
                                        i = 0;
                                        while(station_name[i] != '\0') {
                                            station_name[i] = ' ';
                                            i++;
                                        }
                                        break;
                                case ANISO:
                                        trace_arrow(xy[0], xy[1], vxy[0], vxy[1], scale, &plot_x, &plot_y, &plot_vx, &plot_vy);
                                        ps_plot(plot_x, plot_y, PSL_PEN_MOVE);
                                        ps_plot(plot_vx, plot_vy, PSL_PEN_DRAW_AND_STROKE);
                                        break;
                                case CROSS:
                                        vector_shape = 0.1; /* triangular arrowheads */
                                        trace_cross(xy[0],xy[1],eps1,eps2,theta,scale,v_width,h_length,h_width,vector_shape,outline,pen);
                                        break;
                                case WEDGE:
                                        sprintf (line, "begin wedge number %li",n_rec);
					ps_comment (line);
                                        GMT_geo_to_xy(xy[0], xy[1], &plot_x, &plot_y);
                                        get_trans (xy[0], xy[1], &t11, &t12, &t21, &t22);
                                        paint_wedge (plot_x, plot_y, spin, spinsig, scale, wedge_amp, t11,t12,t21,t22,
						     polygon, fill.rgb, 
						     shade_uncert, efill.rgb,outline);
                                        break;
			}
		}
	if (fp != stdin) GMT_fclose (fp);
	}

	GMT_free((void *) station_name);

	if (gmtdefs.verbose) 
			fprintf (stderr, "psvelo: Number of records read: %li\n", n_rec);

	if (gmtdefs.verbose && rescale_sigma) 
			fprintf (stderr, "psvelo: Rescaling uncertainties by a factor of %f\n", sigma_scale);


	if (skip_if_outside) GMT_map_clip_off ();

	if (pen.texture) ps_setdash (CNULL, 0);

	GMT_map_basemap ();

	GMT_plotend ();

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}
