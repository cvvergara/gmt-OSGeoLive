/*--------------------------------------------------------------------
 *    $Id: pscoupe.c 10201 2014-02-07 19:07:39Z pwessel $
 *
 *    Copyright (c) 1996-2014 by G. Patau
 *    Distributed under the GNU Public Licence
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*

pscoupe will read <x,y> pairs (or <lon,lat>) from inputfile and
plot symbols on a cross-section. Focal mechanisms may be specified 
(double couple or moment tensor) and require additional columns of data.
PostScript code is written to stdout.


 Author:       Genevieve Patau
 Date:         9 September 1992
 Last change : 02 April 2001
 Version:      4
 Roots:        based on psxy.c version 3.0

 */

#include "gmt.h"	/* to have gmt environment */
#include "pslib.h"	/* to have pslib environment */

#include "meca.h"
#include "utilmeca.h"
#include "submeca.h"

     /* symbols for plotting seismicity and/or axis */
#define CROSS       1
#define CIRCLE      2
#define SQUARE      3
#define TRIANGLE    4
#define DIAMOND     5
#define STAR        6
#define HEXAGON     7
#define ITRIANGLE   8

int main (int argc, char **argv)
{

    GMT_LONG     i, j, symbol = 0, n, ix = 0, iy = 1, n_files = 0, fno;
    GMT_LONG     n_args;
    GMT_LONG     n_rec = 0;

    GMT_LONG    outline = FALSE, p_outline = FALSE, t_outline = FALSE;
    GMT_LONG    error = FALSE, nofile = TRUE;
    GMT_LONG    done, get_rgb = FALSE, no_size_needed;
    GMT_LONG    read_size = FALSE, read_proj = FALSE;
    GMT_LONG    read_cmt = FALSE, read_aki = FALSE, read_planes = FALSE;
    GMT_LONG    read_tensor = FALSE, read_axis = FALSE, read_point = FALSE;
    GMT_LONG    plot_dc = FALSE, plot_tensor = FALSE, plot_axis = FALSE;
    GMT_LONG    plot_zerotrace = FALSE, tr_zerotrace = FALSE;
    GMT_LONG    transparence = FALSE, one_size = FALSE, no_label = FALSE;
    GMT_LONG    frame_coupe = FALSE;
    GMT_LONG    skip_if_outside = TRUE;
    GMT_LONG    transparence_old = FALSE, not_defined = FALSE;
    
    double xy[2], west = 0.0, east = 0.0, south = 0.0, north = 0.0;
    double plot_x, plot_y, scale = 0.0;

    char event_title[BUFSIZ], txt_a[GMT_TEXT_LEN];

    char line[BUFSIZ], symbol_type, col[15][GMT_TEXT_LEN], *cpt = CNULL, *not_used = NULL;
    char newfile[GMT_LONG_TEXT], extracted_file[GMT_LONG_TEXT];
    
    FILE *fp = NULL;
    FILE *pnew = NULL, *pextract = NULL; 

    struct GMT_PEN pen, lpen, tr_pen, tz_pen;
    struct GMT_PEN ppen, tpen;
    struct GMT_FILL fill, efill, nofill;
    struct GMT_FILL pfill, tfill;

    struct nodal_plane NP1, NP2, PREF;
    st_me meca, mecar;
    struct MOMENT moment;
    struct M_TENSOR mt, mtr;
    struct AXIS T, N, P, Tr, Nr, Pr;

    GMT_LONG n_plane = 0, n_plane_old = 0;
    GMT_LONG default_justify = 2;
    GMT_LONG justify = default_justify, form = 0;

    double default_fontsize = 9.0;
    double fontsize = default_fontsize;
    double default_offset = (gmtdefs.measure_unit == GMT_CM ? 0.1 : 0.04);
    double offset = default_offset, angle = 0.0;
    double size = 0.0;
    double tmp1, tmp2, tmp3, tmp4, tmp5;
    double n_dep, distance;
    double xlonref, ylatref;
    double fault, depth, dmin, dmax;
    double p_length, p_width, lon1, lat1, lon2, lat2;
    char project_type;
    GMT_LONG syscoord, fuseau = 0;

    double P_x, P_y, T_x, T_y, a_size = GMT_d_NaN;
    char P_sym_type = 0, T_sym_type = 0;
    GMT_LONG P_sym = 0, T_sym = 0;
	void distaz (double lat1,double lon1,double lat2,double lon2,double *distrad,double *distdeg,double *distkm,double *az12rad,double *az12deg,double *az21rad,double *az21deg,GMT_LONG syscoord);

    argc = (int)GMT_begin (argc, argv);
    
    GMT_init_pen (&pen, GMT_PENWIDTH);
    GMT_init_fill (&fill, 0, 0, 0);
    GMT_init_fill (&efill, 255, 255, 255);
    GMT_init_fill (&nofill, -1, -1, -1);
    
    lpen = pen;
    ppen = pen;
    tpen = pen;
    tr_pen = pen;
    tz_pen = pen;
    pfill = fill;
    tfill = efill;
    memset ((void *)&meca, 0, sizeof (meca));
    memset ((void *)&moment, 0, sizeof (moment));
    
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
            
                case 'A':   /* Cross-section definition */
                    read_proj = TRUE;
                    project_type = argv[i][2];
                    j = strlen(argv[i]) - 1;
                    if(argv[i][j] == 'f') frame_coupe = TRUE;

                    if(project_type == 'a' || project_type == 'c') {
                        sscanf (&argv[i][3], "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
                            &lon1, &lat1, &lon2, &lat2, &PREF.dip, &p_width, &dmin, &dmax);
                        syscoord = project_type == 'a' ? 0 : 2;
                        distaz(lat1, lon1, lat2, lon2, &tmp1, &tmp2, &p_length, &tmp3, &PREF.str,
                            &tmp4, &tmp5, syscoord);
                        sprintf(newfile, "A%c%.1f_%.1f_%.1f_%.1f_%.0f_%.0f_%.0f_%.0f",
                            project_type, lon1, lat1, lon2, lat2, PREF.dip, p_width, dmin, dmax);
                        sprintf(extracted_file,"A%c%.1f_%.1f_%.1f_%.1f_%.0f_%.0f_%.0f_%.0f_map",
                            project_type, lon1, lat1, lon2, lat2, PREF.dip, p_width, dmin, dmax);
                    } else if(project_type == 'b' || project_type == 'd') {
                        sscanf (&argv[i][3], "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
                            &lon1, &lat1, &PREF.str, &p_length, &PREF.dip, &p_width, &dmin, &dmax);
                        sprintf(newfile, "A%c%.1f_%.1f_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f",
                            project_type, lon1, lat1, PREF.str, p_length, PREF.dip, p_width, dmin, dmax);
                        sprintf(extracted_file,"A%c%.1f_%.1f_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_map",
                            project_type, lon1, lat1, PREF.str, p_length, PREF.dip, p_width, dmin, dmax);
                    }
                    PREF.rake = 0.;
                    pnew = fopen(newfile, "w");
                    pextract = fopen(extracted_file, "w");
                    if(project_type == 'a' || project_type == 'b')
                        fuseau = gutm(lon1, lat1, &xlonref, &ylatref, fuseau);
                    else {
                        fuseau = - 1;
                        xlonref = lon1;
                        ylatref = lat1;
                    }
                    break;
                case 'E':    /* Set color for extensive parts  */
                    if (!argv[i][2] || (argv[i][2] && GMT_getfill (&argv[i][2], &efill))) {
                        GMT_fill_syntax ('G', " ");
                        error++;
                    }
                    break;
                case 'G':    /* Set color for compressive parts */
                    if (!argv[i][2] || (argv[i][2] && GMT_getfill (&argv[i][2], &fill))) {
                        GMT_fill_syntax ('G', " ");
                        error++;
                    }
                    break;
                case 'L':    /* Draw outline [set outline attributes] */
                    outline = TRUE;
                    if (argv[i][2] && GMT_getpen (&argv[i][2], &lpen)) {
                        GMT_pen_syntax ('L', " ");
                        error++;
                    }
                    break;
                case 'M':    /* Same size for any magnitude */
                    one_size = TRUE;
                    break;
                case 'N':    /* Do not skip points outside border */
                    skip_if_outside = FALSE;
                    break;
                case 'S':    /* Mechanisms : get format [and size] */
                    symbol_type = argv[i][2];
                    strcpy(txt_a, &argv[i][3]);
                    n=0; while (txt_a[n] && txt_a[n] != '/') n++; txt_a[n]=0;
                    scale = GMT_convert_units (txt_a, GMT_INCH);

                    if(strstr(argv[i], "/\0") != NULL) {
                        sscanf(strstr(argv[i], "/\0"), "/%lf/%lf",
                           &fontsize, &offset);
                        if (GMT_IS_ZERO (fontsize)) fontsize = default_fontsize;
                        if (fontsize < 0.0) no_label = TRUE;
                        if (GMT_IS_ZERO (offset)) offset = default_offset;
                        if (argv[i][strlen(argv[i]) - 1] == 'u') justify = 10;
                    }

                    switch (symbol_type) {
                        case 'c':
                            read_cmt = TRUE;
                            plot_dc = TRUE;
                            break;
                        case 'a':
                            read_aki = TRUE;
                            plot_dc = TRUE;
                            break;
                        case 'p':
                            read_planes = TRUE;
                            plot_dc = TRUE;
                            break;
                        case 'm':
                            read_tensor = TRUE;
                            plot_tensor = TRUE;
                            break;
                        case 'd':
                            read_tensor = TRUE;
                            plot_dc = TRUE;
                            break;
                        case 'z':
                            read_tensor = TRUE;
                            plot_zerotrace = TRUE;
                            break;
                        case 'x':
                            read_axis = TRUE;
                            plot_tensor = TRUE;
                            break;
                        case 'y':
                            read_axis = TRUE;
                            plot_dc = TRUE;
                            break;
                        case 't':
                            read_axis = TRUE;
                            plot_zerotrace = TRUE;
                            break;
                        default:
                            error = TRUE;
                            break;
                    }
                    break;
                case 's':    /* Only points : get symbol [and size] */
                    read_point = TRUE;
                    symbol_type = argv[i][2];
                    if(strlen(argv[i]) >3) {
                        strcpy(txt_a, &argv[i][3]);
                        n=0; while (txt_a[n] && txt_a[n] != '/') n++; txt_a[n]=0;
                        scale = GMT_convert_units (txt_a, GMT_INCH);
                    }

                    if(strstr(argv[i], "/\0") != NULL) {
                        sscanf(strstr(argv[i], "/\0"), "/%lf/%lf",
                           &fontsize, &offset);
                        if (GMT_IS_ZERO (fontsize)) fontsize = default_fontsize;
                        if (fontsize < 0.0) no_label = TRUE;
                        if (GMT_IS_ZERO (offset)) offset = default_offset;
                        if (argv[i][strlen(argv[i]) - 1] == 'u') justify = 10;
                    }

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
                            fprintf (stderr, "%s: GMT SYNTAX ERROR -s option:  Unrecognized symbol type %c\n", 
                                argv[0], symbol_type);
                        break;
                    }
                    if (GMT_IS_ZERO (scale)) {
                        read_size = TRUE;
                    }
                    size = scale;

                    break;
                case 'T':
                    transparence = TRUE;
                    sscanf (&argv[i][2], "%" GMT_LL "d",&n_plane);
                    if(strlen(argv[i]) > 3) { /* Set transparent attributes */
                        if (argv[i][2] && GMT_getpen (&argv[i][2], &tr_pen)) {
                            GMT_pen_syntax ('T', " ");
                            error++;
                        }
                    }
                    break;
                case 'W':    /* Set line attributes */
                    if (argv[i][2] && GMT_getpen (&argv[i][2], &pen)) {
                        GMT_pen_syntax ('W', " ");
                        error++;
                    }
                    break;
                case 'Z':    /* Vary symbol color with z */
                    cpt = &argv[i][2];
                    get_rgb = TRUE;
                    break;
                case 'a':    /* plot axis */
                    plot_axis = TRUE;
                    if(strlen(argv[i]) == 2) {
                        strcpy(txt_a,"0.08i");
                        a_size = GMT_convert_units (txt_a, GMT_INCH);
                        P_sym_type = 'c';
                        T_sym_type = 'c';
                    }
                    else {
                        a_size = GMT_convert_units (txt_a, GMT_INCH);
                        strcpy(txt_a, &argv[i][3]);
                        n=0; while (txt_a[n] && txt_a[n] != '/') n++; txt_a[n]=0;
                        a_size = GMT_convert_units (txt_a, GMT_INCH);

                        if(strstr(argv[i], "/\0") != NULL) {
                            strcpy(txt_a,strstr(argv[i], "/\0"));
                            switch (strlen(txt_a)) {
                                case 1:
                                    P_sym_type = 'c';
                                    T_sym_type = 'c';
                                    break;
                                case 2:
                                    P_sym_type = txt_a[1];
                                    T_sym_type = txt_a[1];
                                    break;
                                case 3:
                                    P_sym_type = txt_a[1];
                                    T_sym_type = txt_a[2];
                                    break;
                            }
                        }
                    }

                    switch (P_sym_type) {
                        case 'a':
                            P_sym = STAR;
                            break;
                        case 'c':
                            P_sym = CIRCLE;
                            break;
                        case 'd':
                            P_sym = DIAMOND;
                            break;
                        case 'h':
                            P_sym = HEXAGON;
                            break;
                        case 'i':
                            P_sym = ITRIANGLE;
                            break;
                        case 's':
                            P_sym = SQUARE;
                            break;
                        case 't':
                            P_sym = TRIANGLE;
                            break;
                        case 'x':
                            P_sym = CROSS;
                            break;
                    }
                    switch (T_sym_type) {
                        case 'a':
                            T_sym = STAR;
                            break;
                        case 'c':
                            T_sym = CIRCLE;
                            break;
                        case 'd':
                            T_sym = DIAMOND;
                            break;
                        case 'h':
                            T_sym = HEXAGON;
                            break;
                        case 'i':
                            T_sym = ITRIANGLE;
                            break;
                        case 's':
                            T_sym = SQUARE;
                            break;
                        case 't':
                            T_sym = TRIANGLE;
                            break;
                        case 'x':
                            T_sym = CROSS;
                            break;
                    }
                    break;
                case 'e':    /* Set color for T axis symbol */
                    if (!argv[i][2] || (argv[i][2] && GMT_getfill (&argv[i][2], &tfill))) {
                        GMT_fill_syntax ('e', " ");
                        error++;
                    }
                    break;
                case 'g':    /* Set color for P axis symbol */
                    if (!argv[i][2] || (argv[i][2] && GMT_getfill (&argv[i][2], &pfill))) {
                        GMT_fill_syntax ('g', " ");
                        error++;
                    }
                    GMT_getfill (&argv[i][2], &pfill);
                    break;
                case 'p':    /* Draw outline of P axis symbol [set outline attributes] */
                    p_outline = TRUE;
                    if(strlen(argv[i]) > 2) {
                        if (argv[i][2] && GMT_getpen (&argv[i][2], &ppen)) {
                            GMT_pen_syntax ('p', " ");
                            error++;
                        }
                    }
                    break;
                case 't':    /* Draw outline of T axis symbol [set outline attributes] */
                    t_outline = TRUE;
                    if(strlen(argv[i]) > 2){
                        if (argv[i][2] && GMT_getpen (&argv[i][2], &tpen)) {
                            GMT_pen_syntax ('t', " ");
                            error++;
                        }
                    }
                    break;

                /* Illegal options */
            
                default:    /* Options not recognized */
                    error = TRUE;
                    fprintf(stderr, "%s, -%c option not recognized.\n", argv[0], argv[i][1]);
                    break;
            }
        }
        else
            n_files++;
    }
    
    /* Check that the options selected are mutually consistent */
    
    no_size_needed = (read_cmt || read_planes || read_aki || read_axis || read_tensor);
    error += GMT_check_rgb (pen.rgb)
          + GMT_check_rgb (lpen.rgb)
          + GMT_check_rgb (ppen.rgb)
          + GMT_check_rgb (tpen.rgb)
          + GMT_check_rgb (fill.rgb)
          + GMT_check_rgb (efill.rgb)
          + GMT_check_rgb (pfill.rgb)
          + GMT_check_rgb (tfill.rgb)
          + GMT_check_rgb (gmtdefs.basemap_frame_rgb);

    /* Only one allowed */
    if ((read_cmt +  read_aki + read_planes + read_tensor + read_axis + read_point) > 1) 
        error++;    
    
    if(!read_proj) error++;
    if(no_size_needed && symbol > 0) error++;

    if (!no_size_needed && (symbol > 0 && scale < 0.0)) error++;

    if (argc == 1 || GMT_give_synopsis_and_exit || error) {    /* Display usage */
        fprintf (stderr,"pscoupe %s - Plot seismological symbols on cross-sections\n\n", GMT_VERSION);
        fprintf (stderr,"usage: pscoupe <infiles> -A<params> %s %s\n", GMT_J_OPT, GMT_Rgeo_OPT);
        fprintf (stderr, "   [%s] [-E<fill>] [-G<fill>]\n", GMT_B_OPT);
        fprintf (stderr, "   [%s] [-K] [-L<pen>] [-M] [-N] [-O] [-P]\n", GMT_Ho_OPT);
        fprintf (stderr, "   [-S<format><scale>[/fontsize[/justify/offset/angle/form]]]\n");
        fprintf (stderr, "   [-s<symbol><scale>[/fontsize[/justify/offset/angle/form]]]\n");
        fprintf (stderr, "   [-Tnplane[/<pen>]] [%s]\n", GMT_U_OPT);
        fprintf (stderr, "   [-V] [-W<pen>] [%s] [%s]\n", GMT_X_OPT, GMT_Y_OPT);
        fprintf (stderr, "   [-Z<cpt>] [-a[size[Psymbol[Tsymbol]]]\n\n");
        fprintf (stderr, "   [-p<pen>] [-t<pen>] [-e<fill>] -g<fill>]\n\n");
        fprintf (stderr, "   [%s] [%s]\n", GMT_c_OPT, GMT_t_OPT);
        
        if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
        
        fprintf (stderr, "    <infiles> is one or more files. If no files are given, read standard input\n");
        fprintf (stderr, "        -A Specify cross-section parameters. Choose between\n");
        fprintf (stderr, "           -Aa<lon1/lat1/lon2/lat2/dip/p_width/dmin/dmax>[f]\n");
        fprintf (stderr, "           -Ab<lon1/lat1/strike/p_length/dip/p_width/dmin/dmax>[f]\n");
        fprintf (stderr, "           -Ac<x1/y1/x2/y2/dip/p_width/dmin/dmax>[f]\n");
        fprintf (stderr, "           -Ad<x1/y1/strike/p_length/dip/p_width/dmin/max>[f]\n");
        fprintf (stderr, "           Add f to get the frame from the cross-section parameters.\n");
        GMT_explain_option ('j');
        GMT_explain_option ('R');
        GMT_explain_option ('b');
        fprintf (stderr, "        -E Set color used for extensive parts. [default is white]\n");
        fprintf (stderr, "        -G Set color used for compressive parts. [default is black]\n");
        fprintf (stderr, "           <r/g/b> (each 0-255) for color or <gray> (0-255) for gray-shade [0].\n");
        GMT_explain_option ('H');
        GMT_explain_option ('K');
        fprintf (stderr, "        -L draw line or symbol outline using the current pen (see -W) or sets pen attribute for outline.\n");
        fprintf (stderr, "        -M Same size for any magnitude. Size is given with -S.\n");
        fprintf (stderr, "        -N Do Not skip/clip symbols that fall outside map border [Default will ignore those outside]\n");
        GMT_explain_option ('O');
        GMT_explain_option ('P');
        fprintf (stderr, "        -S Select format type and symbol size (in measure_unit).\n");
        fprintf (stderr, "           Choose format between\n");
        fprintf (stderr, "         (c) Focal mechanisms in Harvard CMT convention\n");
        fprintf (stderr, "             X, Y, depth, strike1, dip1, rake1, strike2, dip2, rake2, moment, event_title\n");
        fprintf (stderr, "             with moment in 2 columns : mantiss and exponent corresponding to seismic moment in dynes-cm\n");
        fprintf (stderr, "         (a) Focal mechanism in Aki & Richard's convention:\n");
        fprintf (stderr, "             X, Y, depth, strike, dip, rake, mag, event_title\n");
        fprintf (stderr, "         (p) Focal mechanism defined with\n");
        fprintf (stderr, "             X, Y, depth, strike1, dip1, strike2, fault, mag, event_title\n");
        fprintf (stderr, "             fault = -1/+1 for a normal/inverse fault\n");
        fprintf (stderr, "         (m) Seismic moment tensor (Harvard CMT, with zero trace)\n");
        fprintf (stderr, "             X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp, event_title\n");
        fprintf (stderr, "         (z) Anisotropic part of seismic moment tensor (Harvard CMT, with zero trace)\n");
        fprintf (stderr, "             X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp, event_title\n");
        fprintf (stderr, "         (d) Best double couple defined from seismic moment tensor (Harvard CMT, with zero trace)\n");
        fprintf (stderr, "             X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp, event_title\n");
        fprintf (stderr, "         (x) Principal axis\n");
        fprintf (stderr, "             X,Y,depth,T_value,T_azimuth,T_plunge,N_value,N_azimuth,N_plunge,\n");
        fprintf (stderr, "             P_value,P_azimuth,P_plunge,exp,event_title\n");
        fprintf (stderr, "         (t) Zero trace moment tensor defined from principal axis\n");
        fprintf (stderr, "             X, Y, depth, T_value, T_azim, T_plunge, N_value, N_azim, N_plunge\n");
        fprintf (stderr, "             P_value, P_azim, P_plunge, exp, newX, newY, event_title\n");
        fprintf (stderr, "         (y) Best double couple defined from principal axis\n");
        fprintf (stderr, "             X,Y,depth,T_value,T_azimuth,T_plunge,N_value,N_azimuth,N_plunge,\n");
        fprintf (stderr, "             P_value,P_azimuth,P_plunge,exp,event_title\n");
        fprintf (stderr, "         Optionally add /fontsize[/offset][u]\n");
        fprintf (stderr, "      Default values are /%g/%f\n", default_fontsize, default_offset);
        fprintf (stderr, "      fontsize < 0 : no label written;\n");
        fprintf (stderr, "      offset is from the limit of the beach ball.\n");
        fprintf (stderr, "      By default label is above the beach ball. Add u to plot it under.\n");
        fprintf (stderr, "        -s to select symbol type and symbol size (in user_unit)\n");
        fprintf (stderr, "          Choose between\n");
        fprintf (stderr, "      st(a)r, (c)ircle, (d)iamond, (h)exagon, (i)nvtriangle, (s)quare, (t)riangle.\n");
        fprintf (stderr, "        -Tn[/<pen>] draw nodal planes and circumference only to provide a transparent beach ball using the current pen (see -W) or sets pen attribute. \n");
        fprintf (stderr, "         n = 1 the only first nodal plane is plotted\n");
        fprintf (stderr, "         n = 2 the only second nodal plane is plotted\n");
        fprintf (stderr, "         n = 0 both nodal planes are plotted.\n");
        fprintf (stderr, "         If moment tensor is required, nodal planes overlay moment tensor.\n");
        GMT_explain_option ('U');
        GMT_explain_option ('V');
        fprintf (stderr, "        -W sets pen attributes [width = %gp, color = (%d/%d/%d), texture = solid line].\n", 
            pen.width, pen.rgb[0], pen.rgb[1], pen.rgb[2]);
        fprintf (stderr, "        -Z Use cpt-file to assign colors based on depth-value in 3rd column\n");
        fprintf (stderr, "        -a plots axis. Default symbols are circles.\n");
        fprintf (stderr, "        -p draws P_symbol outline using the current pen (see -W) or sets pen attribute for outline.\n");
        fprintf (stderr, "        -t draws T_symbol outline using the current pen (see -W) or sets pen attribute for outline.\n");
        fprintf (stderr, "        -g Sets color used for P_symbol. [default is compressive parts color]\n");
        fprintf (stderr, "        -e Sets color used for T_symbol. [default is extensive parts color]\n");
        GMT_explain_option ('X');
        GMT_explain_option ('c');
        GMT_explain_option (':');
        GMT_explain_option ('.');
        exit (EXIT_FAILURE);
    }

    if (get_rgb) GMT_read_cpt (cpt);

    if (n_files > 0)
        nofile = FALSE;
    else
        n_files = 1;
    n_args = (argc > 1) ? argc : 2;
    
    if(frame_coupe) {
        west = 0.;
        east = p_length;
        south = dmin;
        north = dmax;
        if (GMT_IS_ZERO (PREF.dip)) PREF.dip = 1.0;
    }

    GMT_err_fail (GMT_map_setup (west, east, south, north), "");

    GMT_plotinit (argc, argv);
    
    GMT_setpen (&pen);
    ps_setfont (gmtdefs.annot_font[0]);
    
    if (skip_if_outside) GMT_map_clip_on (GMT_no_rgb, 3);
        
    ix = (gmtdefs.xy_toggle[0]);    iy = 1 - ix;

    done = FALSE;


    for (fno = 1; !done && fno < n_args; fno++) {    /* Loop over all input files */
        if (!nofile && argv[fno][0] == '-') continue;
        if (nofile) {
            fp = GMT_stdin;
            done = TRUE;
        }
        else if ((fp = GMT_fopen (argv[fno], "r")) == NULL) {
            fprintf (stderr, "pscoupe: Cannot open file %s\n", argv[fno]);
            continue;
        }

        if (!nofile && gmtdefs.verbose) {
            fprintf (stderr, "pscoupe: Working on file %s\n", argv[fno]);
            sprintf (line, "File: %s", argv[fno]);
            ps_comment (line);
        }
        if (GMT_io.io_header[GMT_IN]) 
            for (i = 0; i < GMT_io.n_header_recs; i++) not_used = GMT_fgets (line, BUFSIZ, fp);
        
        while (GMT_fgets (line, BUFSIZ, fp)) {
            n_rec++;
 	    memset ((void *)col, 0, 15 * GMT_TEXT_LEN * sizeof (char));
	    memset ((void *)event_title, 0, BUFSIZ * sizeof (char));
           if (read_cmt) {
                sscanf (line, "%s %s %s %s %s %s %s %s %s %s %s %s %s %[^\n]\n",
                    col[0], col[1], col[2], col[3], col[4], col[5], col[6], 
                    col[7], col[8], col[9], col[10], col[11], col[12], event_title);
                if (strlen(event_title) <= 0) sprintf(event_title,"\n");
            }
            else if (read_aki) {
                sscanf (line, "%s %s %s %s %s %s %s %s %s %[^\n]\n",
                    col[0], col[1], col[2], col[3], col[4], col[5], col[6], 
                    col[7], col[8], event_title);
                if (strlen(event_title) <= 0) sprintf(event_title,"\n");
            }
            else if (read_planes) {
                sscanf (line, "%s %s %s %s %s %s %s %s %s %s %[^\n]\n",
                    col[0], col[1], col[2], col[3], col[4], col[5], col[6], 
                    col[7], col[8], col[9], event_title);
                if (strlen(event_title) <= 0) sprintf(event_title,"\n");
            }
            else if (read_axis) {
                sscanf (line, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %[^\n]\n", 
                    col[0], col[1], col[2], col[3], col[4], col[5], col[6], col[7], 
                    col[8], col[9], col[10], col[11], col[12], col[13], col[14], event_title);
                if (strlen(event_title) <= 0) sprintf(event_title,"\n");
            }   
            else if (read_tensor) {
                sscanf (line, "%s %s %s %s %s %s %s %s %s %s %s %s %[^\n]\n", 
                    col[0], col[1], col[2], col[3], col[4], col[5], col[6], col[7], 
                    col[8], col[9], col[10], col[11], event_title);
                if (strlen(event_title) <= 0) sprintf(event_title,"\n");
            }
            else if(symbol > 0 && read_size) {
                sscanf (line, "%s %s %s %s %[^\n]\n", col[0], col[1], col[2], col[3], event_title);
            }
            else {
                sscanf (line, "%s %s %s %[^\n]\n", col[0], col[1], col[2], event_title);
            }

            xy[ix] = atof (col[0]);
            xy[iy] = atof (col[1]);
            depth = atof (col[2]);
            
            if (dans_coupe(xy[0], xy[1], depth, xlonref, ylatref, fuseau, PREF.str, 
                    PREF.dip, p_length, p_width, &distance, &n_dep)) {
                xy[0] = distance;
                xy[1] = n_dep;
            }
            else {
                xy[0] = -1;
                xy[1] = -1;
            }

            if (skip_if_outside) {
                GMT_map_outside (xy[0], xy[1]);
                if (GMT_abs (GMT_x_status_new) > 1 || GMT_abs (GMT_y_status_new) > 1) continue;
            }

            if (get_rgb) {
                GMT_get_rgb_from_z (depth, fill.rgb);
            }

            GMT_geo_to_xy (xy[0], xy[1], &plot_x, &plot_y);

            if(symbol > 0) {
                if (read_size) {
                    size = GMT_convert_units (col[3], GMT_INCH);
                }
                fprintf(pnew, "%f %f %f %s\n", distance, n_dep, depth, col[3]);
                fprintf(pextract, "%s", line);

                switch (symbol) {
                    case STAR:
                        ps_star (plot_x, plot_y, size, fill.rgb, outline);
                        break;
                    case CROSS:
                        ps_cross (plot_x, plot_y, size);
                        break;
                    case CIRCLE:
                        ps_circle (plot_x, plot_y, size, fill.rgb, outline);
                        break;
                    case SQUARE:
                        ps_square (plot_x, plot_y, size, fill.rgb, outline);
                        break;
                    case HEXAGON:
                        ps_hexagon (plot_x, plot_y, size, fill.rgb, outline);
                        break;
                    case TRIANGLE:
                        ps_triangle (plot_x, plot_y, size, fill.rgb, outline);
                        break;
                    case DIAMOND:
                        ps_diamond (plot_x, plot_y, size, fill.rgb, outline);
                        break;
                    case ITRIANGLE:
                        ps_itriangle (plot_x, plot_y, size, fill.rgb, outline);
                        break;
                }

            } else if (read_cmt) {
                meca.NP1.str = atof(col[3]);
                meca.NP1.dip = atof(col[4]);
                meca.NP1.rake = atof(col[5]);
                meca.NP2.str = atof(col[6]);
                meca.NP2.dip = atof(col[7]);
                meca.NP2.rake = atof(col[8]);
                moment.mant = atof(col[9]);
                moment.exponent = atoi(col[10]);
                if(moment.exponent == 0)
                meca.magms = atof(col[9]);

                rot_meca(meca, PREF, &mecar);

            } else if(read_aki) {
                NP1.str = atof(col[3]);
                NP1.dip = atof(col[4]);
                NP1.rake = atof(col[5]);
                meca.magms = atof(col[6]);    

                moment.exponent = 0;
                moment.mant = meca.magms;
                define_second_plane(NP1, &NP2);

                meca.NP1.str = NP1.str;
                meca.NP1.dip = NP1.dip;
                meca.NP1.rake = NP1.rake;
                meca.NP2.str = NP2.str;
                meca.NP2.dip = NP2.dip;
                meca.NP2.rake = NP2.rake;

                rot_meca(meca, PREF, &mecar);

            } else if(read_planes) {
                meca.NP1.str = atof(col[3]);
                meca.NP1.dip = atof(col[4]);
                meca.NP2.str = atof(col[5]);
                fault = atof(col[6]);
                meca.magms = atof(col[7]);

                moment.exponent = 0;
                moment.mant = meca.magms;
                meca.NP2.dip = computed_dip2(meca.NP1.str, meca.NP1.dip, meca.NP2.str);
                if(meca.NP2.dip == 1000.) {
                    not_defined = TRUE;
                    transparence_old = transparence;
                    n_plane_old = n_plane;
                    transparence = TRUE;
                    n_plane = 1;
                    meca.NP1.rake = 1000.;
                    if(gmtdefs.verbose) {
                        fprintf(stderr, "WARNING : second plane is not defined for event %s only first plane is plotted.\n", line);
                    }
                }
                else {
                    meca.NP1.rake = computed_rake2(meca.NP2.str, meca.NP2.dip, meca.NP1.str, meca.NP1.dip, fault);
                }
                meca.NP2.rake = computed_rake2(meca.NP1.str, meca.NP1.dip, meca.
NP2.str, meca.NP2.dip, fault);

                rot_meca(meca, PREF, &mecar);

            } else if(read_axis) {
                T.val = atof(col[3]);
                T.str = atof(col[4]);
                T.dip = atof(col[5]);
                T.e = atoi(col[12]);

                N.val = atof(col[6]);
                N.str = atof(col[7]);
                N.dip = atof(col[8]);
                N.e = atoi(col[12]);

                P.val = atof(col[9]);
                P.str = atof(col[10]);
                P.dip = atof(col[11]);
                P.e = atoi(col[12]);

/*
F. A. Dahlen and Jeroen Tromp, Theoretical Seismology, Princeton, 1998, p.167.
Definition of scalar moment.
*/
                meca.moment.exponent = T.e;
                meca.moment.mant = sqrt(squared(T.val) + squared(N.val) + squared(P.val)) / M_SQRT2;
                meca.magms = 0.;

/* normalization by M0 */
                T.val /= meca.moment.mant;
                N.val /= meca.moment.mant;
                P.val /= meca.moment.mant;

                rot_axis(T, PREF, &Tr);
                rot_axis(N, PREF, &Nr);
                rot_axis(P, PREF, &Pr);
                Tr.val = T.val;
                Nr.val = N.val;
                Pr.val = P.val;
                Tr.e = T.e;
                Nr.e = N.e;
                Pr.e = P.e;

                if(plot_dc || transparence) {
                    axe2dc(Tr, Pr, &NP1, &NP2);
                    meca.NP1.str = NP1.str;
                    meca.NP1.dip = NP1.dip;
                    meca.NP1.rake = NP1.rake;
                    meca.NP2.str = NP2.str;
                    meca.NP2.dip = NP2.dip;
                    meca.NP2.rake = NP2.rake;
                }

            } else if(read_tensor) {
                for(i=3;i<9;i++) mt.f[i-3] = atof(col[i]);
                mt.expo = atoi(col[i]);

                moment.exponent = mt.expo;
/*
F. A. Dahlen and Jeroen Tromp, Theoretical Seismology, Princeton, 1998, p.167.
Definition of scalar moment.
*/
                moment.mant = sqrt(squared(mt.f[0]) + squared(mt.f[1]) + squared(mt.f[2]) + 2.*(squared(mt.f[3]) + squared(mt.f[4]) + squared(mt.f[5]))) / M_SQRT2;
                meca.magms = 0.;

/* normalization by M0 */
                for(i=0;i<=5;i++) mt.f[i] /= moment.mant;

                rot_tensor(mt, PREF, &mtr);
                GMT_momten2axe(mtr, &T, &N, &P);

                if(plot_dc || transparence) {
                    axe2dc(T, P, &NP1, &NP2);
                    meca.NP1.str = NP1.str;
                    meca.NP1.dip = NP1.dip;
                    meca.NP1.rake = NP1.rake;
                    meca.NP2.str = NP2.str;
                    meca.NP2.dip = NP2.dip;
                    meca.NP2.rake = NP2.rake;
                }
            }
                
            if(no_size_needed) {
                if(one_size) {
                    moment.mant = 4.;
                    moment.exponent = 23;
                }

                size = (computed_mw(moment, meca.magms) / 5.) * scale;

                fprintf(pextract, "%s", line);
                if(read_axis) {
                    fprintf(pnew, "%f %f %f %f %f %f %f %f %f %f %f %f %ld 0 0 %s\n", 
                        xy[0], xy[1], depth, Tr.val, Tr.str, Tr.dip, Nr.val, Nr.str, Nr.dip, 
                        Pr.val, Pr.str, Pr.dip, moment.exponent, event_title);
                    T = Tr;
                    N = Nr;
                    P = Pr;
                }
                else if(read_tensor) {
                    fprintf(pnew, "%f %f %f %f %f %f %f %f %f %ld 0 0 %s\n",
                        xy[0], xy[1], depth, mtr.f[0], mtr.f[1], mtr.f[2], mtr.f[3], mtr.f[4], mtr.f[5],
                        moment.exponent, event_title);
                    mt = mtr;
                }
                else {
                    fprintf(pnew, "%f %f %f %f %f %f %f %f %f %f %ld 0 0 %s\n",
                        xy[0], xy[1], depth, mecar.NP1.str, mecar.NP1.dip, mecar.NP1.rake, 
                        mecar.NP2.str, mecar.NP2.dip, mecar.NP2.rake,
                        moment.mant, moment.exponent, event_title);
                    meca = mecar;
                }

                if(plot_tensor) {
                    GMT_setpen(&lpen);
                    ps_tensor(plot_x,plot_y,size,T,N,P,fill.rgb,efill.rgb,outline,plot_zerotrace);
                }

                if(tr_zerotrace) {
                    GMT_setpen(&tz_pen);
                    ps_tensor(plot_x,plot_y,size,T,N,P,nofill.rgb,nofill.rgb,tr_zerotrace,tr_zerotrace);
                }

                if(transparence) {
                    GMT_setpen(&tr_pen);
                    if(n_plane == 0) ps_meca(plot_x, plot_y, meca, size);
                    else ps_plan(plot_x, plot_y, meca, size, n_plane);
                    if(not_defined) {
                        not_defined = FALSE;
                        transparence = transparence_old;
                        n_plane = n_plane_old;
                    }
                } else if(plot_dc) {
                    GMT_setpen(&lpen);
                    ps_mechanism(plot_x,plot_y,meca,size,fill.rgb,efill.rgb,outline);
                }                                        

                if(plot_axis) {
                    if(read_tensor || read_axis)
                        axis2xy(plot_x, plot_y, size, P.str, P.dip, T.str, T.dip, &P_x, &P_y, &T_x, &T_y);
                    else {
                        dc_to_axe(meca, &T, &N, &P);
                        axis2xy(plot_x, plot_y, size, P.str, P.dip, T.str, T.dip, &P_x, &P_y, &T_x, &T_y);
                    }
                    GMT_setpen(&ppen);
                    switch (P_sym) {
                        case STAR:
                            ps_star (P_x, P_y, a_size, pfill.rgb, p_outline);
                            break;
                        case CROSS:
                            ps_cross (P_x, P_y, a_size);
                            break;
                        case CIRCLE:
                            ps_circle (P_x, P_y, a_size, pfill.rgb, p_outline);
                            break;
                        case SQUARE:
                            ps_square (P_x, P_y, a_size, pfill.rgb, p_outline);
                            break;
                        case HEXAGON:
                            ps_hexagon (P_x, P_y, a_size, pfill.rgb, p_outline);
                            break;
                        case TRIANGLE:
                            ps_triangle (P_x, P_y, a_size, pfill.rgb, p_outline);
                            break;
                        case ITRIANGLE:
                            ps_itriangle (P_x, P_y, a_size, pfill.rgb, p_outline);
                            break;
                        case DIAMOND:
                            ps_diamond (P_x, P_y, a_size, pfill.rgb, p_outline);
                            break;
                    }
                    GMT_setpen(&tpen);
                    switch (T_sym) {
                        case STAR:
                            ps_star (T_x, T_y, a_size, tfill.rgb, t_outline);
                            break;
                        case CROSS:
                            ps_cross (T_x, T_y, a_size);
                            break;
                        case CIRCLE:
                            ps_circle (T_x, T_y, a_size, tfill.rgb, t_outline);
                            break;
                        case SQUARE:
                            ps_square (T_x, T_y, a_size, tfill.rgb, t_outline);
                            break;
                        case HEXAGON:
                            ps_hexagon (T_x, T_y, a_size, tfill.rgb, t_outline);
                            break;
                        case TRIANGLE:
                            ps_triangle (T_x, T_y, a_size, tfill.rgb, t_outline);
                            break;
                        case ITRIANGLE:
                            ps_itriangle (T_x, T_y, a_size, tfill.rgb, t_outline);
                            break;
                        case DIAMOND:
                            ps_diamond (T_x, T_y, a_size, tfill.rgb, t_outline);
                            break;
                    }
                }
            }

            if(! no_label) {
                GMT_setpen(&pen);
                switch (justify) {
                    case 2 :
                        ps_text(plot_x, plot_y + size * 0.5 + offset, fontsize, event_title, angle, justify, form);
                        break;
                    case 10 :
                        ps_text(plot_x, plot_y - size * 0.5 - offset, fontsize, event_title, angle, justify, form);
                        break;
                }
            }

        }
        if (fp != stdin) GMT_fclose (fp);
    }
    
    if (gmtdefs.verbose) 
        fprintf (stderr, "pscoupe: Number of records read: %li\n", n_rec);


    if (skip_if_outside) GMT_map_clip_off ();
    
    ps_setpaint (gmtdefs.basemap_frame_rgb);

    if (pen.texture) ps_setdash (CNULL, 0);

    GMT_map_basemap ();

    GMT_plotend ();
    
    GMT_end (argc, argv);

    exit (EXIT_SUCCESS);
}
