/*--------------------------------------------------------------------
 *    $Id: pspolar.c 9923 2012-12-18 20:45:53Z pwessel $ 
 *
 *    Copyright (c) 1996-2013 by G. Patau
 *    Distributed under the GNU Public Licence
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * pspolar will read azimuth, take-off angle of seismic ray and polarities
 * in stations and plots polarities on the inferior focal half-sphere.
 * A variety of symbols may be specified.
 * Only one event may be plotted at a time.
 * PostScript code is written to stdout.
 *
 * Author:         Genevieve Patau
 * Date:           19-OCT-1995 (psprojstations)
 * Version:        4
 * Roots:          heavily based on psxy.c
 * Last modified : 02-APR-2001
 *
 */

#include "gmt.h"	/* to have gmt environment */
#include "pslib.h"	/* to have pslib environment */
#include <string.h>

#define CROSS       2
#define POINT       3
#define CIRCLE      4
#define SQUARE      5
#define TRIANGLE    6
#define DIAMOND     7
#define STAR        8
#define HEXAGON     9
#define ITRIANGLE  10

#define POINTSIZE 0.005

int main (int argc, char **argv)
{
    GMT_LONG     i, symbol = 0, n, n_files = 0, fno;
    GMT_LONG     n_args;
    GMT_LONG     form_s = 0, justify_s = 5;
    
    GMT_LONG error = FALSE, nofile = TRUE;
    GMT_LONG done, greenwich, label_s = FALSE;
    GMT_LONG change_position = FALSE;
    GMT_LONG get_position = FALSE, get_size = FALSE, get_symbol = FALSE;
    GMT_LONG outline_E = FALSE, outline_G = FALSE, outline_F = FALSE;
    GMT_LONG old_GMT_world_map, skip_if_outside = TRUE;
    GMT_LONG plot_polS = FALSE, vecS = FALSE, scolor = FALSE, outline_s = FALSE;
    GMT_LONG def_cpen = FALSE, def_fpen = FALSE, def_tpen = FALSE;
    GMT_LONG def_gpen = FALSE, def_epen = FALSE, def_spen = FALSE;
    GMT_LONG hypo = FALSE;
    
    double west = 0.0, east = 0.0, south = 0.0, north = 0.0;
    double plot_x, plot_y, symbol_size = 0.0, symbol_size2;
    double lon, lat, plot_x0, plot_y0;
    double new_lon, new_lat, new_plot_x0, new_plot_y0;
    double radius, ech = 0., azimut, ih, plongement;
    double c_pointsize = 0.015, fontsize_s = 12.0;
    double angle_s = 0.0;
    double azS, sizeS = GMT_d_NaN;
    double v_width = GMT_d_NaN, h_length = GMT_d_NaN, h_width = GMT_d_NaN, shape;
    double si, co;

    char line[BUFSIZ], symbol_type, col[4][GMT_TEXT_LEN];
    char pol, *not_used = NULL;
    char txt_a[GMT_TEXT_LEN],txt_b[GMT_TEXT_LEN],txt_c[GMT_TEXT_LEN], txt_d[GMT_TEXT_LEN];
    char stacode[GMT_TEXT_LEN];
    
    FILE *fp = NULL;
    
    struct GMT_PEN pen, cpen, fpen, gpen, epen, spen, tpen;
    struct GMT_FILL fill, ffill, gfill, efill, sfill;
    struct GMT_FILL black, nofill;
    
    argc = (int)GMT_begin (argc, argv);
    
    GMT_init_pen (&pen, GMT_PENWIDTH);
    GMT_init_fill (&nofill, -1, -1, -1); 
    GMT_init_fill (&black, 0, 0, 0);     
    GMT_init_fill (&efill, 250, 250, 250);
    
    fill = nofill;
    ffill = nofill;
    sfill = nofill;
    gfill = black;

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
                case '\0':
                    error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
                    break;
                
                /* Supplemental parameters */
            
                case 'C':       /* New coordinates */
                         change_position = TRUE;
                         sscanf(&argv[i][2], "%lf/%lf", &new_lon, &new_lat);
                         if(strchr(argv[i], 'W')) {
                             GMT_getpen (strchr(argv[i]+1, 'W')+1, &cpen);
                             def_cpen = TRUE;
                         }
                         if(strchr(argv[i], 'P')) {
                             sscanf(strchr(argv[i]+1, 'P')+1, "%lf", &c_pointsize);
                         }
                         break;
                case 'D':       /* Coordinates */
                         get_position = TRUE;
                         sscanf(&argv[i][2], "%lf/%lf", &lon, &lat);
                         break;
                case 'E':        /* Set color for station in extensive part */
                         GMT_getfill (&argv[i][2], &efill);
                         break;
                case 'e':        /* Outline station symbol in extensive part */
                         outline_E = TRUE;
                         if(strlen(argv[i]) > 2) {
                             GMT_getpen (&argv[i][2], &epen);
                             def_epen = TRUE;
                         }
                         break;
                case 'F':        /* Set background color of beach ball */
                         GMT_getfill (&argv[i][2], &ffill);
                         break;
                case 'f':        /* Outline beach ball */
                         outline_F = TRUE;
                         if(strlen(argv[i]) > 2) {
                             GMT_getpen (&argv[i][2], &fpen);
                             def_fpen = TRUE;
                         }
                         break;
                case 'G':        /* Set color for station in compressive part */
                         GMT_getfill (&argv[i][2], &gfill);
                         break;
                case 'g':        /* Outline station symbol in compressive part */
                         outline_G = TRUE;
                         if(strlen(argv[i]) > 2) {
                             GMT_getpen (&argv[i][2], &gpen);
                             def_gpen = TRUE;
                         }
                         break;
                case 'h':    /* Use HYPO71 format */
                         hypo = TRUE;
                         break;
                case 'M':    /* Focal sphere size */
                         get_size = TRUE;
                         sscanf(&argv[i][2], "%s", txt_a);
                         ech = GMT_convert_units (txt_a, GMT_INCH);
                         break;
                case 'N':        /* Do not skip points outside border */
                         skip_if_outside = FALSE;
                         break;
                case 'S':        /* Get symbol [and size] */
                         symbol_type = argv[i][2];
                         symbol_size = GMT_convert_units (&argv[i][3], GMT_INCH);
                         get_symbol = TRUE;
                         switch (symbol_type) {
                             case 'a':
                                 symbol = STAR;
                                 break;
                             case 'c':
                                 symbol = CIRCLE;
                                 break;
                             case 'd':
                                 symbol = DIAMOND;
                                 break;
                             case 'h':
                                 symbol = HEXAGON;
                                 break;
                             case 'i':
                                 symbol = ITRIANGLE;
                                 break;
                             case 'p':
                                 symbol = POINT;
                                 break;
                             case 's':
                                 symbol = SQUARE;
                                 break;
                             case 't':
                                 symbol = TRIANGLE;
                                 break;
                             case 'x':
                                 symbol = CROSS;
                                 break;
                             default:
                                 error = TRUE;
                                 fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  Unrecognized symbol type %c\n", argv[0], symbol_type);
                                 break;
                         }
                         break;
                case 's':        /* Get S polarity */
                         plot_polS = TRUE;
                         strcpy(txt_a, &argv[i][3]);
                         n=0; while (txt_a[n] && txt_a[n] != '/' && txt_a[n] != 'V' && txt_a[n] != 'G' && txt_a[n] != 'L') n++; txt_a[n]=0;
                         sizeS = GMT_convert_units (txt_a, GMT_INCH);

                         if(strchr(argv[i], 'V')) {
                             vecS = TRUE;
                             strcpy(txt_a,strchr(argv[i], 'V'));
                             if(strncmp(txt_a,"VG",(size_t)2) == 0 || strncmp(txt_a,"VL",(size_t)2) == 0 || strlen(txt_a) == 1) {
                                 v_width = 0.03; h_length = 0.12; h_width = 0.1; shape = gmtdefs.vector_shape;
                                 if (!gmtdefs.measure_unit) {
                                     v_width = 0.075; h_length = 0.3; h_width = 0.25; shape = gmtdefs.vector_shape;
                                 }
                             }
                             else {
                                strcpy(txt_a, strchr(argv[i], 'V')+1);
                                strcpy(txt_b, strchr(txt_a+1, '/')+1);
                                strcpy(txt_c, strchr(txt_b+1, '/')+1);
                                strcpy(txt_d, strchr(txt_c+1, '/')+1);
                                n=0; while (txt_a[n] && txt_a[n] != '/') n++; txt_a[n]=0;
                                n=0; while (txt_b[n] && txt_b[n] != '/') n++; txt_b[n]=0;
                                n=0; while (txt_c[n] && txt_c[n] != '/') n++; txt_c[n]=0;
                                n=0; while (txt_d[n] && txt_d[n] != '/' && txt_d[n] != 'L' && txt_d[n] != 'G') n++; txt_d[n]=0;
                                v_width = GMT_convert_units (txt_a, GMT_INCH);
                                h_length = GMT_convert_units (txt_b, GMT_INCH);
                                h_width = GMT_convert_units (txt_c, GMT_INCH);
                                shape = atof(txt_d);
                             }
                         }
                         if(strchr(argv[i], 'G')) {
                             sscanf (strchr(argv[i]+1,'G')+1,"%d/%d/%d",&sfill.rgb[0],&sfill.rgb[1],&sfill.rgb[2]);
                             sprintf(txt_a, "%d/%d/%d",sfill.rgb[0],sfill.rgb[1],sfill.rgb[2]);
                             GMT_getfill (txt_a, &sfill);
                             scolor = TRUE;
                         }
                         if(strchr(argv[i], 'L')) outline_s = TRUE;
                         break;
                case 'T':       /* Information about label printing */
                        label_s = TRUE;
                        if (strlen(argv[i]) > 2) {
                             sscanf (&argv[i][2], "%lf/%" GMT_LL "d/%" GMT_LL "d/%lf/", &angle_s,
                                     &form_s, &justify_s, &fontsize_s);
                        }
                        break;
                case 't':       /* Set color for station label */
                         GMT_getpen (&argv[i][2], &tpen);
                         def_tpen = TRUE;
                         break;

                case 'W':        /* Set line attributes */
                         GMT_getpen (&argv[i][2], &pen);
                         break;
                    
                /* Illegal options */
            
                default:        /* Options not recognized */
                         error = TRUE;
                         GMT_default_error (argv[i][1]);
                         break;
            }
        }
        else
            n_files++;
    }

    /* Check that the options selected are consistent */
    
    if (!project_info.region_supplied) {
        error++;
    }
    if(ech <= 0.) {
        error++;
    }
    if(get_position + get_size + get_symbol < 3) {
        error++;
    }
    
    if (argc == 1 || GMT_give_synopsis_and_exit || error) {    /* Display usage */
        fprintf (stderr,"%s %s - Plot polarities on the inferior focal half-sphere on maps\n\n",argv[0], GMT_VERSION);
        fprintf (stderr,"usage: argv[0] <infiles> %s %s\n", GMT_J_OPT, GMT_Rgeo_OPT);
        fprintf (stderr, " -Dlongitude/latitude -Msize[i/c] -S<symbol><size>[i/c]\n");
        fprintf (stderr, " [-A] [%s] [-Clongitude/latitude[W<pen>][Ppointsize]] [-E<fill>]\n", GMT_B_OPT);
        fprintf (stderr, " [-e[<pen>]] [-F<fill>] [-f[<pen>]] [-G<fill>] [-g[<pen>]] [%s] [-K] [-N] [-O] [-P]\n", GMT_Ho_OPT);
        fprintf (stderr, " [-s<half-size>/[V[<v_width/h_length/h_width/shape]][G<r/g/b>][L]\n");
        fprintf (stderr, " [-T[<labelinfo>]] [-t<pen>] [%s] [-V] [-W<pen>]\n", GMT_U_OPT);
        fprintf (stderr, " [%s] [%s] [%s]\n", GMT_X_OPT, GMT_Y_OPT, GMT_c_OPT);
        
        if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
        
        fprintf (stderr, "    <infiles> is one or more files.  If no, read standard input\n");
        GMT_explain_option ('j');
        GMT_explain_option ('R');
        fprintf (stderr, "        -D Set longitude/latitude\n");
        fprintf (stderr, "        -M Set size of beach ball in %s\n", GMT_unit_names[gmtdefs.measure_unit]);
        fprintf (stderr, "        -S to select symbol type and symbol size (in %s).  Choose between\n", GMT_unit_names[gmtdefs.measure_unit]);
        fprintf (stderr, "           st(a)r, (c)ircle, (d)iamond, (h)exagon, (i)nvtriangle\n");
        fprintf (stderr, "           (p)oint, (s)quare, (t)riangle, (x)cross\n");
        fprintf (stderr, "\n\tOPTIONS:\n");
        GMT_explain_option ('B');
        GMT_explain_option ('b');
        fprintf (stderr, "        -C Set new_longitude/new_latitude[W<pen>][Ppointsize]\n");
        fprintf (stderr, "           A line will be plotted between both positions\n");
        fprintf (stderr, "           Default is width = 3, color = current pen and pointsize = 0.015\n");
        fprintf (stderr, "        -E Specify color symbol for station in extensive part.\n");
        fprintf (stderr, "           Fill can be either <r/g/b> (each 0-255) for color \n");
        fprintf (stderr, "           or <gray> (0-255) for gray-shade [0].\n");
        fprintf (stderr, "           Default is light gray.\n");
        fprintf (stderr, "        -e Outline of station symbol in extensive part.\n");
        fprintf (stderr, "           Default is current pen.\n");
        fprintf (stderr, "        -F Specify background color of beach ball. It can be\n");
        fprintf (stderr, "           <r/g/b> (each 0-255) for color or <gray> (0-255) for gray-shade [0].\n");
        fprintf (stderr, "           Default is no fill\n");
        fprintf (stderr, "        -f Outline beach ball\n");
        fprintf (stderr, "           Add <pen attributes> if not current pen.\n");
        fprintf (stderr, "        -G Specify color symbol for station in compressive part. Fill can be either\n");
        fprintf (stderr, "           Fill can be either <r/g/b> (each 0-255) for color\n");
        fprintf (stderr, "           or <gray> (0-255) for gray-shade [0].\n");
        fprintf (stderr, "           Add L[<pen>] to outline\n");
        fprintf (stderr, "           Default is black.\n");
        fprintf (stderr, "        -g Outline of station symbol in compressive part.\n");
        fprintf (stderr, "           Add <pen attributes> if not current pen.\n");
        fprintf (stderr, "        -h Use special format derived from HYPO71 output.\n");
        GMT_explain_option ('H');
        GMT_explain_option ('K');
        fprintf (stderr, "        -N Do Not skip/clip symbols that fall outside map border\n");
        fprintf (stderr, "           [Default will ignore those outside]\n");
        GMT_explain_option ('O');
        GMT_explain_option ('P');
        fprintf (stderr, "        -s to plot S polarity azimuth.\n");
        fprintf (stderr, "           Azimuth of S polarity is in last column.\n");
        fprintf (stderr, "           It may be a vector (V option) or a segment. Give half-size in cm.\n");
        fprintf (stderr, "           L option is for outline\n");
        fprintf (stderr, "           -s<half-size>/[V[<v_width/h_length/h_width/shape>]][G<r/g/b>][L]\n");
        fprintf (stderr, "           Default definition of v is 0.075/0.3/0.25/1\n");
        fprintf (stderr, "           Outline is current pen\n");
        fprintf (stderr, "        -T[<info about labal printing>] to write station code.\n");
        fprintf (stderr, "           <angle/form/justify/fontsize in points>\n");
        fprintf (stderr, "           Default is 0.0/0/5/12\n");
        fprintf (stderr, "        -t sets pen attributes to write station codes [default is current pen]\n");
        GMT_explain_option ('U');
        GMT_explain_option ('V');
        fprintf (stderr, "        -W sets current pen attributes [width = %gp, color = (%d/%d/%d), texture = solid line].\n", pen.width, pen.rgb[0], pen.rgb[1], pen.rgb[2]);
        GMT_explain_option ('X');
        GMT_explain_option ('c');
        GMT_explain_option ('.');
        exit (EXIT_FAILURE);
    }
    
    if(!def_cpen) {
        cpen = pen; cpen.width = 3;
    }                                         /* pen for change position */
    if(!def_fpen) fpen = pen;                 /* outline beach ball */
    if(!def_gpen) gpen = pen;                 /* outline compressive stations */
    if(!def_epen) epen = pen;                 /* outline extensive stations */
    if(!def_spen) spen = pen;                 /* outline S_pol segment */
    if(!def_tpen) tpen = pen;                 /* pen to print station name */

    if (n_files > 0)
        nofile = FALSE;
    else
        n_files = 1;
    n_args = (argc > 1) ? argc : 2;
    
    greenwich = (west < 0.0 || east <= 0.0);
    
    GMT_err_fail (GMT_map_setup (west, east, south, north), "");

    GMT_plotinit (argc, argv);
    
    if (label_s) ps_setfont (gmtdefs.annot_font[0]);

    if (symbol > 0 && skip_if_outside) GMT_map_clip_on (GMT_no_rgb, 3);
    
    old_GMT_world_map = GMT_world_map;
    
    done = FALSE;


    for (fno = 1; !done && fno < n_args; fno++) {    /* Loop over all input files */
        if (!nofile && argv[fno][0] == '-') continue;
        if (nofile) {
            fp = GMT_stdin;
            done = TRUE;
        }
        else if ((fp = GMT_fopen (argv[fno], "r")) == NULL) {
            fprintf (stderr, "%s: Cannot open file %s\n", argv[0], argv[fno]);
            continue;
        }

        if (!nofile && gmtdefs.verbose) {
            fprintf (stderr, "%s: Working on file %s\n", argv[0], argv[fno]);
            sprintf (line, "File: %s", argv[fno]);
            ps_comment (line);
        }
        if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) not_used = GMT_fgets (line, 512, fp);
        
        GMT_world_map = TRUE;
        
        GMT_geo_to_xy(lon, lat, &plot_x0, &plot_y0);
        if(change_position) {
            GMT_setpen (&cpen);
            GMT_geo_to_xy(new_lon, new_lat, &new_plot_x0, &new_plot_y0);
            ps_circle(plot_x0, plot_y0, c_pointsize, cpen.rgb, 1);
            ps_plot(plot_x0, plot_y0, PSL_PEN_MOVE);
            ps_plot(new_plot_x0, new_plot_y0, PSL_PEN_DRAW_AND_STROKE);
            plot_x0 = new_plot_x0;
            plot_y0 = new_plot_y0;
        }
        if (skip_if_outside) {
            GMT_map_outside (lon, lat);
            if (GMT_abs (GMT_x_status_new) > 1 || GMT_abs (GMT_y_status_new) > 1) continue;
        }

        GMT_setpen (&fpen);
        ps_circle (plot_x0, plot_y0, ech, ffill.rgb, outline_F);

        while (GMT_fgets (line, BUFSIZ, fp)) {
            switch (hypo) {
                case 0 :
                    if(!plot_polS) {
                        sscanf (line, "%s %lf %lf %c", stacode, &azimut, &ih, &pol);
                    }
                    else {
                        n = sscanf (line, "%s %lf %lf %c %lf", stacode, &azimut, &ih, &pol, &azS);
                        if(n == 4) azS = -1.;
                    }
                    break;
                case 1 :
 	    	    memset ((void *)col, 0, 4 * GMT_TEXT_LEN * sizeof (char));
                   if(!plot_polS) {
                        sscanf (line, "%s %s %s %s %lf %lf %c", col[0], col[1], col[2], stacode, &azimut, &ih, col[3]);
                        pol = col[3][2];
                    }
                    else {
                        n = sscanf (line, "%s %s %s %s %lf %lf %c %lf", col[0], col[1], col[2], stacode, &azimut, &ih, col[3], &azS);
                        pol = col[3][2];
                        if(n == 7) azS = -1.;
                    }
                    break;
            }
        
            if(strcmp(col[0],"000000")!=0) {
                plongement = (ih - 90.) * M_PI / 180.;
                if(plongement  < 0.) {
                    plongement = -plongement;
                    azimut += 180 ;
                    symbol_size2 = symbol_size * 0.8;
                }
                else symbol_size2 = symbol_size;
                radius = sqrt(1. - sin(plongement));
                if(radius >= 0.97) radius = 0.97;
                azimut += 180;
                azimut *= M_PI / 180.;
                sincos (azimut, &si, &co);
                plot_x = radius * si * ech / 2. + plot_x0;
                plot_y = radius * co * ech / 2. + plot_y0;
                if (symbol == CROSS || symbol == POINT) ps_setpaint (fill.rgb);
            
                if(label_s) {
                    GMT_setpen (&tpen);
                    switch (justify_s) {
                        case 11 :
                            ps_text(plot_x-symbol_size2-0.1, plot_y-symbol_size2-0.1, fontsize_s, stacode, angle_s, 11, form_s);
                            break;
                        case 10 :
                            ps_text(plot_x, plot_y-symbol_size2-0.1, fontsize_s, stacode, angle_s, 10, form_s);
                            break;
                        case 9 :
                            ps_text(plot_x+symbol_size2+0.1, plot_y-symbol_size2-0.1, fontsize_s, stacode, angle_s, 9, form_s);
                            break;
                        case 7:
                            ps_text(plot_x-symbol_size2-0.1, plot_y, fontsize_s, stacode, angle_s, 7, form_s);
                            break;
                        case 6:
                            ps_text(plot_x, plot_y, fontsize_s, stacode, angle_s, 6, form_s);
                            break;
                        case 5:
                            ps_text(plot_x+symbol_size2+0.1, plot_y, fontsize_s, stacode, angle_s, 5, form_s);
                            break;
                        case 3:
                            ps_text(plot_x-symbol_size2-0.1, plot_y+symbol_size2+0.1, fontsize_s, stacode, angle_s, 3, form_s);
                            break;
                        case 2:
                            ps_text(plot_x, plot_y+symbol_size2+0.1, fontsize_s, col[3], angle_s, 2, form_s);
                            break;
                        case 1:
                            ps_text(plot_x+symbol_size2+0.1, plot_y+symbol_size2+0.1, fontsize_s, stacode, angle_s, 1, form_s);
                            break;
                    }
                }
    
                switch (symbol) {
                    case STAR:
                              if(pol == 'u' || pol == 'U' || pol == 'c' || pol == 'C' || pol == '+') {
                                  GMT_setpen (&gpen);
                                  ps_star (plot_x, plot_y, symbol_size2, gfill.rgb, outline_G);
                              }
                              else if(pol == 'r' || pol == 'R' || pol == 'd' || pol == 'D' || pol == '-') {
                                  GMT_setpen (&epen);
                                  ps_star (plot_x, plot_y, symbol_size2, efill.rgb, outline_E);
                              }
                              else {
                                  GMT_setpen (&pen);
                                  ps_cross (plot_x, plot_y, symbol_size2);
                              }
                              break;
                    case CROSS:
                              GMT_setpen (&pen);
                              ps_cross (plot_x, plot_y, symbol_size2);
                              break;
                    case POINT:
                              GMT_setpen (&pen);
                              ps_cross (plot_x, plot_y, POINTSIZE);
                              break;
                    case CIRCLE:
                              if(pol == 'u' || pol == 'U' || pol == 'c' || pol == 'C' || pol == '+') {
                                  GMT_setpen (&gpen);
                                  ps_circle (plot_x, plot_y, symbol_size2, gfill.rgb, outline_G);
                              }
                              else if(pol == 'r' || pol == 'R' || pol == 'd' || pol == 'D' || pol == '-') {
                                  GMT_setpen (&epen);
                                  ps_circle (plot_x, plot_y, symbol_size2, efill.rgb, outline_E);
                              }
                              else {
                                  GMT_setpen (&pen);
                                  ps_cross (plot_x, plot_y, symbol_size2);
                              }
                              break;
                    case SQUARE:
                              if(pol == 'u' || pol == 'U' || pol == 'c' || pol == 'C' || pol == '+') {
                                  GMT_setpen (&gpen);
                                  ps_square (plot_x, plot_y, symbol_size2, gfill.rgb, outline_G);
                              }
                              else if(pol == 'r' || pol == 'R' || pol == 'd' || pol == 'D' || pol == '-') {
                                  GMT_setpen (&epen);
                                  ps_square (plot_x, plot_y, symbol_size2, efill.rgb, outline_E);
                              }
                              else {
                                  GMT_setpen (&pen);
                                  ps_cross (plot_x, plot_y, symbol_size2);
                              }
                              break;
                    case HEXAGON:
                              if(pol == 'u' || pol == 'U' || pol == 'c' || pol == 'C' || pol == '+') {
                                  GMT_setpen (&gpen);
                                  ps_hexagon (plot_x, plot_y, symbol_size2, gfill.rgb, outline_G);
                              }
                              else if(pol == 'r' || pol == 'R' || pol == 'd' || pol == 'D' || pol == '-') {
                                  GMT_setpen (&epen);
                                  ps_hexagon (plot_x, plot_y, symbol_size2, efill.rgb, outline_E);
                              }
                              else {
                                  GMT_setpen (&pen);
                                  ps_cross (plot_x, plot_y, symbol_size2);
                              }
                              break;
                    case TRIANGLE:
                              if(pol == 'u' || pol == 'U' || pol == 'c' || pol == 'C' || pol == '+') {
                                  GMT_setpen (&gpen);
                                  ps_triangle (plot_x, plot_y, symbol_size2, gfill.rgb, outline_G);
                              }
                              else if(pol == 'r' || pol == 'R' || pol == 'd' || pol == 'D' || pol == '-') {
                                  GMT_setpen (&epen);
                                  ps_triangle (plot_x, plot_y, symbol_size2, efill.rgb, outline_E);
                              }
                              else {
                                  GMT_setpen (&pen);
                                  ps_cross (plot_x, plot_y, symbol_size2);
                              }
                              break;
                    case ITRIANGLE:
                              if(pol == 'u' || pol == 'U' || pol == 'c' || pol == 'C' || pol == '+') {
                                  GMT_setpen (&gpen);
                                  ps_itriangle (plot_x, plot_y, symbol_size2, gfill.rgb, outline_G);
                              }
                              else if(pol == 'r' || pol == 'R' || pol == 'd' || pol == 'D' || pol == '-') {
                                  GMT_setpen (&epen);
                                  ps_itriangle (plot_x, plot_y, symbol_size2, efill.rgb, outline_E);
                              }
                              else {
                                  GMT_setpen (&pen);
                                  ps_cross (plot_x, plot_y, symbol_size2);
                              }
                              break;
                    case DIAMOND:
                              if(pol == 'u' || pol == 'U' || pol == 'c' || pol == 'C' || pol == '+') {
                                  GMT_setpen (&gpen);
                                  ps_diamond (plot_x, plot_y, symbol_size2, gfill.rgb, outline_G);
                              }
                              else if(pol == 'r' || pol == 'R' || pol == 'd' || pol == 'D' || pol == '-') {
                                  GMT_setpen (&epen);
                                  ps_diamond (plot_x, plot_y, symbol_size2, efill.rgb, outline_E);
                              }
                              else {
                                  GMT_setpen (&pen);
                                  ps_cross (plot_x, plot_y, symbol_size2);
                              }
                              break;
                }
                if(plot_polS && azS >= 0.) {
                    GMT_setpen (&spen);
                    sincos (azS*M_PI/180., &si, &co);
                    if(vecS) {
                        ps_vector(plot_x - sizeS*si, plot_y - sizeS*co, 
                            plot_x + sizeS*si, plot_y + sizeS*co, v_width, 
                            h_length, h_width, gmtdefs.vector_shape, sfill.rgb, outline_s);
                    }
                    else { 
                        if(scolor) ps_setpaint (sfill.rgb);
                        else ps_setpaint (pen.rgb);
                        ps_plot(plot_x - sizeS*si, plot_y - sizeS*co, PSL_PEN_MOVE);
                        ps_plot(plot_x + sizeS*si, plot_y + sizeS*co, PSL_PEN_DRAW_AND_STROKE); 
                    }
                }
            }
        }
        if (fp != stdin) GMT_fclose (fp);
    }
    
    if (skip_if_outside) GMT_map_clip_off ();

    GMT_world_map = old_GMT_world_map;
    
    if (pen.texture) ps_setdash (CNULL, 0);

    GMT_map_basemap ();

    GMT_plotend ();
    
    GMT_end (argc, argv);

    exit (EXIT_SUCCESS);
}
