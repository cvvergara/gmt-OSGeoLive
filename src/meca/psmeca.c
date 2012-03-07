/*--------------------------------------------------------------------
 *    $Id: psmeca.c,v 1.61 2011/04/28 16:23:53 remko Exp $
 *
 *    Copyright (c) 1996-2011 by G. Patau
 *    Distributed under the GNU Public Licence
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*

psmeca will read <x,y> pairs (or <lon,lat>) from inputfile and
plot symbols on a map. Focal mechanisms are specified (double couple or
moment tensor).
PostScript code is written to stdout.


 Author:     Genevieve Patau
 Date:         7 July 1998
 Version:    4
 Roots:      based on psxy.c    

 Last change  18 February 2000
 */

#include "gmt.h"        /* to have gmt environment */
#include "pslib.h"      /* to have pslib environment */

#include "meca.h"
#include "utilmeca.h"

     /* symbols for plotting axis */
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

    GMT_LONG     i, symbol = 0, n, ix = 0, iy = 1, n_files = 0, fno, last = 0;
    GMT_LONG     n_args;
    GMT_LONG     n_rec = 0;

    GMT_LONG    outline = FALSE, p_outline = FALSE, t_outline = FALSE;
    GMT_LONG    error = FALSE, nofile = TRUE;
    GMT_LONG    done, get_rgb = FALSE, no_size_needed;
    GMT_LONG    read_cmt = FALSE, read_aki = FALSE, read_planes = FALSE;
    GMT_LONG    read_tensor = FALSE, plot_dc = FALSE, plot_zerotrace = FALSE;
    GMT_LONG    read_axis = FALSE;
    GMT_LONG    read_change_position = FALSE;
    GMT_LONG    transparence = FALSE, one_size = FALSE, no_label = FALSE;
    GMT_LONG    plot_axis = FALSE, tr_zerotrace = FALSE;
    GMT_LONG    skip_if_outside = TRUE;
    GMT_LONG    def_lpen = FALSE, def_cpen = FALSE, def_tr_pen = FALSE;
    GMT_LONG    def_ppen = FALSE, def_tpen = FALSE, def_tz_pen = FALSE;
    GMT_LONG    def_pfill = FALSE, def_tfill = FALSE;
    GMT_LONG    draw_box = FALSE;
    GMT_LONG    transparence_old = FALSE, not_defined = FALSE;
    
    double xy[2], xynew[2], west = 0.0, east = 0.0, south = 0.0, north = 0.0;
    double plot_x, plot_y, scale = 0.0;
    double plot_xnew, plot_ynew;
    double t11 = 1.0, t12 = 0.0, t21 = 0.0, t22 = 1.0;
    double delaz;

    char string[BUFSIZ], event_title[BUFSIZ], txt_a[GMT_TEXT_LEN];

    char line[BUFSIZ], symbol_type, col[15][GMT_TEXT_LEN], *cpt = CNULL, *not_used = NULL;
    
    FILE *fp = NULL;
    

    struct GMT_PEN pen, lpen, cpen, tr_pen, tz_pen;
    struct GMT_PEN ppen, tpen;
    struct GMT_FILL fill, efill, nofill;
    struct GMT_FILL pfill, tfill, bfill;

    st_me meca;
    struct MOMENT moment;
    struct M_TENSOR mt;
    struct AXIS T, N, P;

    GMT_LONG n_plane = 0, n_plane_old = 0;
    GMT_LONG default_justify = 2;
    GMT_LONG justify = default_justify, form = 0;
    GMT_LONG new = 1;

    double default_fontsize = 9.0;
    double fontsize = default_fontsize;
    double default_offset = (gmtdefs.measure_unit == GMT_CM ? 0.1 : 0.04);
    double default_pointsize = 0.005;
    double pointsize = default_pointsize, offset = default_offset, angle = 0.0;
    double fault, depth, depmin = 0.0, depmax = 900.0;
    double size, a_size = GMT_d_NaN;
    double P_x, P_y, T_x, T_y;
    char P_sym_type = 0, T_sym_type = 0;
    GMT_LONG P_sym = 0, T_sym = 0;

    argc = (int)GMT_begin (argc, argv);
    
    GMT_init_pen (&pen, GMT_PENWIDTH);
    GMT_init_pen (&lpen, GMT_PENWIDTH);
    GMT_init_pen (&cpen, GMT_PENWIDTH);
    GMT_init_pen (&tr_pen, GMT_PENWIDTH);
    GMT_init_pen (&tz_pen, GMT_PENWIDTH);
    GMT_init_pen (&ppen, GMT_PENWIDTH);
    GMT_init_pen (&tpen, GMT_PENWIDTH);

    GMT_init_fill (&fill, 0, 0, 0);
    GMT_init_fill (&efill, 255, 255, 255);
    GMT_init_fill (&nofill, -1, -1, -1);
    GMT_init_fill (&bfill, 255, 255, 255);
    GMT_init_fill (&pfill, 255, 255, 255);
    GMT_init_fill (&tfill, 255, 255, 255);
    event_title[0] = 0;
    memset ((void *)&meca, 0, sizeof (meca));
    
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
            
                case 'C':   /* Change position [set line attributes] */
                    read_change_position = TRUE;
                    if (strlen(argv[i]) > 2) {
                        if (strchr(argv[i], 'P')) {
                            strcpy(txt_a, strchr(argv[i]+1, 'P')+1);
                            pointsize = GMT_convert_units (txt_a, GMT_INCH);
                        }
                        if (argv[i][2] != 'P') {
                            def_cpen = TRUE;
                            strcpy(txt_a, &argv[i][2]);
                            n=0; while (txt_a[n] != 'P') n++; txt_a[n]=0;
                            if (GMT_getpen (txt_a, &cpen)) {
                                GMT_pen_syntax ('C', " ");
                                error++;
                            }
                        }
                    }
                    break;
                case 'D':    /* Plot events between depmin and depmax deep */
                    sscanf(&argv[i][2], "%lf/%lf", &depmin, &depmax);
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
                    if (argv[i][2]) {
                        def_lpen = TRUE;
                        if (GMT_getpen (&argv[i][2], &lpen)) {
                            GMT_pen_syntax ('L', " ");
                            error++;
                        }
                    }
                    break;
                case 'M':    /* Same size for any magnitude */
                    one_size = TRUE;
                    break;
                case 'N':    /* Do not skip points outside border */
                    skip_if_outside = FALSE;
                    break;
                case 'r':    /* draw box around text */
                    draw_box = TRUE;
                    if (strlen(argv[i]) > 2 ) GMT_getfill (&argv[i][2], &bfill);
                    break;
                case 'S':    /* Get symbol [and size] */
                    symbol_type = argv[i][2];
                    strcpy(txt_a, &argv[i][3]);
                    n=0;
		    while (txt_a[n] && txt_a[n] != '/') n++;
		    txt_a[n]=0;
                    scale = GMT_convert_units (txt_a, GMT_INCH);

                    if (strstr(argv[i], "/") != NULL) {
                        sscanf(strstr(argv[i], "/"), "/%lf/%lf", &fontsize, &offset);
                        if (GMT_IS_ZERO (fontsize)) fontsize = default_fontsize;
                        if (fontsize < 0.0) no_label = TRUE;
                        if (GMT_IS_ZERO (offset)) offset = default_offset;
                        if (argv[i][strlen(argv[i]) - 1] == 'u') justify = 10;
                    }

                    switch (symbol_type) {
                        case 'c':
                            read_cmt = TRUE;
                            break;
                        case 'a':
                            read_aki = TRUE;
                            break;
                        case 'p':
                            read_planes = TRUE;
                            break;
                        case 'x':
                            read_axis = TRUE;
                            break;
                        case 'y':
                            read_axis = plot_dc = TRUE;
                            break;
                        case 't':
                            read_axis = plot_zerotrace = TRUE;
                            break;
                        case 'm':
                            read_tensor = TRUE;
                            break;
                        case 'd':
                            read_tensor = plot_dc = TRUE;
                            break;
                        case 'z':
                            read_tensor = plot_zerotrace = TRUE;
                            break;
                        default:
                            error = TRUE;
                            break;
                    }
                    break;
                case 'T':
                    transparence = TRUE;
                    sscanf (&argv[i][2], "%" GMT_LL "d",&n_plane);
                    if (strlen(argv[i]) > 4) { /* Set transparent attributes */
                        def_tr_pen = TRUE;
                        if (GMT_getpen (&argv[i][4], &tr_pen)) {
                            GMT_pen_syntax ('T', " ");
                            error++;
                        }
                    }
                    break;
                case 'z': /* overlay zerotrace moment tensor */
                    tr_zerotrace = TRUE;
                    if (strlen(argv[i]) > 2) { /* Set transparent attributes */
                        def_tz_pen = TRUE;
                        if (GMT_getpen (&argv[i][2], &tz_pen)) {
                            GMT_pen_syntax ('z', " ");
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
                    if (strlen(argv[i]) == 2) {
                        strcpy(txt_a,"0.08i");
                        P_sym_type = 'c';
                        T_sym_type = 'c';
                    }
                    else {
                        strcpy(txt_a, &argv[i][3]);
                        n=0;
			while (txt_a[n] && txt_a[n] != '/') n++;
			txt_a[n]=0;
                        a_size = GMT_convert_units (txt_a, GMT_INCH);

                        if (strstr(argv[i], "/") != NULL) {
                            strcpy(txt_a,strstr(argv[i], "/"));
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
                    if (strlen(argv[i]) > 2) {
                        def_tfill = TRUE;
                        if (GMT_getfill (&argv[i][2], &tfill)) {
                            GMT_fill_syntax ('e', " ");
                            error++;
                        }
                    }
                    break;
                case 'g':    /* Set color for P axis symbol */
                    if (strlen(argv[i]) > 2) {
                        def_pfill = TRUE;
                        if (GMT_getfill (&argv[i][2], &pfill)) {
                            GMT_fill_syntax ('g', " ");
                            error++;
                        }
                    }
                    break;
                case 'p':    /* Draw outline of P axis symbol [set outline attributes] */
                    p_outline = TRUE;
                    if (strlen(argv[i]) > 2) {
                        def_ppen = TRUE;
                        if (GMT_getpen (&argv[i][2], &ppen)) {
                            GMT_pen_syntax ('p', " ");
                            error++;
                        }
                    }
                    break;
                case 't':    /* Draw outline of T axis symbol [set outline attributes] */
                    t_outline = TRUE;
                    if (strlen(argv[i]) > 2){
                        def_tpen = TRUE;
                        if (GMT_getpen (&argv[i][2], &tpen)) {
                            GMT_pen_syntax ('t', " ");
                            error++;
                        }
                    } 
                    break;
                case 'o':   /* use psvelomeca format (without depth in 3rd column) */
                    new = 0;
                    break;

                /* Illegal options */
            
                default:    /* Options not recognized */
                    error = TRUE;
                    break;
            }
        }
        else
            n_files++;
    }
    
    /* Check that the options selected are mutually consistent */
    
    no_size_needed = (read_cmt || read_planes || read_aki || read_tensor || read_axis);
    error += GMT_check_rgb (pen.rgb)
          + GMT_check_rgb (cpen.rgb)
          + GMT_check_rgb (lpen.rgb)
          + GMT_check_rgb (ppen.rgb)
          + GMT_check_rgb (tpen.rgb)
          + GMT_check_rgb (fill.rgb)
          + GMT_check_rgb (efill.rgb)
          + GMT_check_rgb (pfill.rgb)
          + GMT_check_rgb (tfill.rgb)
          + GMT_check_rgb (bfill.rgb)
          + GMT_check_rgb (gmtdefs.basemap_frame_rgb);

    /* Only one allowed */
    if (no_size_needed > 1) error++;    
    if (!no_size_needed && (symbol > 1 && scale <= 0.0)) error++;

    if (get_rgb && !new) error++;

    if (argc == 1 || GMT_give_synopsis_and_exit || error) {   /* Display usage */
        fprintf (stderr,"psmeca %s - Plot seismological symbols on maps\n\n", GMT_VERSION);
        fprintf (stderr,"usage: psmeca <infiles> %s %s\n", GMT_J_OPT, GMT_Rgeo_OPT);
        fprintf (stderr, "   -S<format><scale>[/fontsize[/justify/offset/angle/form]]\n");
        fprintf (stderr, "   [%s] [-C[<pen>][P<pointsize>]]\n", GMT_B_OPT);
        fprintf (stderr, "   [-Ddepmin/depmax] [-E<fill>] [-G<fill>]\n");
        fprintf (stderr, "   [%s] [-K] [-L<pen>] [-M] [-N] [-O] [-P] [-r]\n", GMT_Ho_OPT);
        fprintf (stderr, "   [-Tnplane[/<pen>]] [%s]\n", GMT_U_OPT);
        fprintf (stderr, "   [-V] [-W<pen>] [%s] [%s]\n", GMT_X_OPT, GMT_Y_OPT);
        fprintf (stderr, "   [-Z<cpt>] [-z] [-a[size[/Psymbol[Tsymbol]]]\n\n");
        fprintf (stderr, "   [-p<pen>] [-t<pen>] [-e<fill>] -g<fill>]\n\n");
        fprintf (stderr, "   [-o] [%s] [%s]\n", GMT_c_OPT, GMT_t_OPT);
        
        if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
        
        fprintf (stderr, "    <infiles> is one or more files. If no files are given, read standard input\n");
        GMT_explain_option ('j');
        GMT_explain_option ('R');
        GMT_explain_option ('b');
        fprintf (stderr, "        -C<pen attributes>\n");
        fprintf (stderr, "          offset focal mechanisms to the latitude and longitude specified in the last two columns of the input file before label.\n");
        fprintf (stderr, "          Default pen attributes is default pen.\n");
        fprintf (stderr, "          A line is plotted between both positions.\n");
        fprintf (stderr, "          A small circle is plotted at the initial location. Add P<pointsize value> to change the size of the circle.\n");
        fprintf (stderr, "        -Ddepmin/depmax Plot events between depmin and depmax deep.\n");
        fprintf (stderr, "        -E Set color used for extensive parts. [default is white]\n");
        fprintf (stderr, "        -G Set color used for compressive parts. [default is black]\n");
        fprintf (stderr, "           <r/g/b> (each 0-255) for color or <gray> (0-255) for gray-shade [0].\n");
        GMT_explain_option ('H');
        GMT_explain_option ('K');
        fprintf (stderr, "        -L draw line or symbol outline using the default pen (see -W) or sets pen attribute for outline.\n");
        fprintf (stderr, "        -M Same size for any magnitude. Size is given with -S.\n");
        fprintf (stderr, "        -N Do Not skip/clip symbols that fall outside map border [Default will ignore those outside]\n");
        GMT_explain_option ('O');
        GMT_explain_option ('P');
        fprintf (stderr, "        -r draw a box around text.\n");
        fprintf (stderr, "        -S Select format type and symbol size (in %s).\n", GMT_unit_names[gmtdefs.measure_unit]);
        fprintf (stderr, "           Choose format between\n");
        fprintf (stderr, "         (c) Focal mechanisms in Harvard CMT convention\n");
        fprintf (stderr, "             X, Y, depth, strike1, dip1, rake1, strike2, dip2, rake2, moment, newX, newY, event_title\n");
        fprintf (stderr, "             with moment in 2 columns : mantiss and exponent corresponding to seismic moment in dynes-cm\n");
        fprintf (stderr, "         (a) Focal mechanism in Aki & Richard's convention:\n");
        fprintf (stderr, "             X, Y, depth, strike, dip, rake, mag, newX, newY, event_title\n");
        fprintf (stderr, "         (p) Focal mechanism defined with\n");
        fprintf (stderr, "             X, Y, depth, strike1, dip1, strike2, fault, mag, newX, newY, event_title\n");
        fprintf (stderr, "             fault = -1/+1 for a normal/inverse fault\n");
        fprintf (stderr, "         (m) Sesmic moment tensor (Harvard CMT, with zero trace)\n");
        fprintf (stderr, "             X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp, newX, newY, event_title\n");
        fprintf (stderr, "         (z) Anisotropic part of seismic moment tensor (Harvard CMT, with zero trace)\n");
        fprintf (stderr, "             X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp, event_title\n");
        fprintf (stderr, "         (d) Best double couple defined from seismic moment tensor (Harvard CMT, with zero trace)\n");
        fprintf (stderr, "             X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp, newX, newY, event_title\n");
        fprintf (stderr, "         (x) Principal axis\n");
        fprintf (stderr, "             X, Y, depth, T_value, T_azim, T_plunge, N_value, N_azim, N_plunge\n");
        fprintf (stderr, "             P_value, P_azim, P_plunge, exp, newX, newY, event_title\n");
        fprintf (stderr, "         (t) Zero trace moment tensor defined from principal axis\n");
        fprintf (stderr, "             X, Y, depth, T_value, T_azim, T_plunge, N_value, N_azim, N_plunge\n");
        fprintf (stderr, "             P_value, P_azim, P_plunge, exp, newX, newY, event_title\n");
        fprintf (stderr, "         (y) Best double couple defined from principal axis\n");
        fprintf (stderr, "             X, Y, depth, T_value, T_azim, T_plunge, N_value, N_azim, N_plunge\n");
        fprintf (stderr, "             P_value, P_azim, P_plunge, exp, newX, newY, event_title\n");
        fprintf (stderr, "         Use -o option for old (psvelomeca) format (not depth in third column)\n");
        fprintf (stderr, "         Optionally add /fontsize[/offset][u]\n");
        fprintf (stderr, "      Default values are /%g/%f\n", default_fontsize, default_offset);
        fprintf (stderr, "      fontsize < 0 : no label written;\n");
        fprintf (stderr, "      offset is from the limit of the beach ball.\n");
        fprintf (stderr, "      By default label is above the beach ball. Add u to plot it under.\n");
        fprintf (stderr, "        -Tn[/<pen>] draw nodal planes and circumference only to provide a transparent beach ball using the default pen (see -W) or sets pen attribute. \n");
        fprintf (stderr, "         n = 1 the only first nodal plane is plotted\n");
        fprintf (stderr, "         n = 2 the only second nodal plane is plotted\n");
        fprintf (stderr, "         n = 0 both nodal planes are plotted.\n");
        fprintf (stderr, "         If moment tensor is required, nodal planes overlay moment tensor.\n");
        fprintf (stderr, "        -z overlays zero trace moment tensor.\n");
        GMT_explain_option ('U');
        GMT_explain_option ('V');
        fprintf (stderr, "        -W sets default pen attributes [width = %gp, color = (%d/%d/%d), texture = solid line].\n", 
            pen.width, pen.rgb[0], pen.rgb[1], pen.rgb[2]);
        fprintf (stderr, "        -Z Use cpt-file to assign colors based on depth-value in 3rd column\n");
        fprintf (stderr, "        -a plots axis. Default symbols are circles.\n");
        fprintf (stderr, "        -p draws P_symbol outline using the default pen (see -W) or sets pen attribute for outline.\n");
        fprintf (stderr, "        -t draws T_symbol outline using the default pen (see -W) or sets pen attribute for outline.\n");
        fprintf (stderr, "        -g Sets color used for P_symbol. [default is compressive parts color]\n");
        fprintf (stderr, "        -e Sets color used for T_symbol. [default is extensive parts color]\n");
        fprintf (stderr, "        -o Use psvelomeca format (Without depth in third column)\n");
        GMT_explain_option ('X');
        GMT_explain_option ('c');
        GMT_explain_option (':');
        GMT_explain_option ('.');
        exit (EXIT_FAILURE);
    }

    if (!def_lpen) lpen = pen;
    if (!def_cpen) cpen = pen;
    if (!def_ppen) ppen = pen;
    if (!def_tpen) tpen = pen;
    if (!def_tr_pen) tr_pen = pen;
    if (!def_tz_pen) tz_pen = pen;
    if (!def_pfill) pfill = fill;
    if (!def_tfill) tfill = efill;

    if (get_rgb) GMT_read_cpt (cpt);

    if (n_files > 0)
        nofile = FALSE;
    else
        n_files = 1;
    n_args = (argc > 1) ? argc : 2;
    
    GMT_err_fail (GMT_map_setup (west, east, south, north), "");

    GMT_plotinit (argc, argv);
    
    GMT_setpen (&pen);
    ps_setfont (gmtdefs.annot_font[0]);
    
    if (skip_if_outside) GMT_map_clip_on (GMT_no_rgb, 3);
    
    ix = (gmtdefs.xy_toggle[0]);
    iy = 1 - ix;

    done = FALSE;

    for (fno = 1; !done && fno < n_args; fno++) {    /* Loop over all input files */
        if (!nofile && argv[fno][0] == '-') continue;
        if (nofile) {
            fp = GMT_stdin;
            done = TRUE;
        }
        else if ((fp = GMT_fopen (argv[fno], "r")) == NULL) {
            fprintf (stderr, "psmeca: Cannot open file %s\n", argv[fno]);
            continue;
        }

        if (!nofile && gmtdefs.verbose) {
            fprintf (stderr, "psmeca: Working on file %s\n", argv[fno]);
            sprintf (line, "File: %s", argv[fno]);
            ps_comment (line);
        }
        if (GMT_io.io_header[GMT_IN]) 
            for (i = 0; i < GMT_io.n_header_recs; i++) not_used = GMT_fgets (line, BUFSIZ, fp);
        

        while (GMT_fgets (line, BUFSIZ, fp)) {
            n_rec++;
            if (read_cmt) {
                sscanf (line, "%s %s %s %s %s %s %s %s %s %s %s %s %[^\n]\n",
                    col[0], col[1], col[2], col[3], col[4], col[5], col[6], 
                    col[7], col[8], col[9], col[10], col[11], string); 
	    	last = 11;
            }
            else if (read_aki) {
                sscanf (line, "%s %s %s %s %s %s %s %s %[^\n]\n",
                    col[0], col[1], col[2], col[3], col[4], col[5], col[6], 
                    col[7], string);
	    	last = 7;
            }
            else if (read_planes) {
                sscanf (line, "%s %s %s %s %s %s %s %s %s %[^\n]\n",
                    col[0], col[1], col[2], col[3], col[4], col[5], col[6], 
                    col[7], col[8], string);
	    	last = 8;
            }
            else if (read_axis) {
                sscanf (line, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %[^\n]\n", 
                    col[0], col[1], col[2], col[3], col[4], col[5], col[6], col[7], 
                    col[8], col[9], col[10], col[11], col[12], col[13], string);
	    	last = 13;
            }
            else if (read_tensor) {
                sscanf (line, "%s %s %s %s %s %s %s %s %s %s %s %[^\n]\n", 
                    col[0], col[1], col[2], col[3], col[4], col[5], col[6], col[7], 
                    col[8], col[9], col[10], string);
	    	last = 10;
            }

	    /* Immediately skip locations outside of the map area */

            xy[ix] = atof (col[0]);
            xy[iy] = atof (col[1]);
            
            if (skip_if_outside) {
                GMT_map_outside (xy[0], xy[1]);
                if (GMT_abs (GMT_x_status_new) > 1 || GMT_abs (GMT_y_status_new) > 1) continue;
            }

	    /* In new (psmeca) input format, third column is depth.
	       Skip record when depth is out of range. Also read an extra column. */

            if (new) {
                depth = atof (col[2]);
		if (depth < depmin || depth > depmax) continue;
                if (get_rgb) GMT_get_rgb_from_z (depth, fill.rgb);
                sscanf (string, "%s %[^\n]\n", col[last+1], event_title);
            }
            else 
                strcpy(event_title, string);
            if (strlen(event_title) <= 0) sprintf(event_title,"\n");

	    /* Gather and transform the input records, depending on type */

            if (read_cmt) {
                meca.NP1.str = atof(col[2+new]);
                meca.NP1.dip = atof(col[3+new]);
                meca.NP1.rake = atof(col[4+new]);
                meca.NP2.str = atof(col[5+new]);
                meca.NP2.dip = atof(col[6+new]);
                meca.NP2.rake = atof(col[7+new]);
                meca.moment.mant = atof(col[8+new]);
                meca.moment.exponent = atoi(col[9+new]);
                if (meca.moment.exponent == 0) meca.magms = atof(col[8+new]);
            }
	    
	    else if (read_aki) {
                meca.NP1.str = atof(col[2+new]);
                meca.NP1.dip = atof(col[3+new]);
                meca.NP1.rake = atof(col[4+new]);
                meca.magms = atof(col[5+new]);    

                meca.moment.exponent = 0;
                define_second_plane(meca.NP1, &meca.NP2);
            }
	    
	    else if (read_planes) {
                meca.NP1.str = atof(col[2+new]);
                meca.NP1.dip = atof(col[3+new]);
                meca.NP2.str = atof(col[4+new]);
                fault = atof(col[5+new]);
                meca.magms = atof(col[6+new]);

                meca.moment.exponent = 0;
                meca.NP2.dip = computed_dip2(meca.NP1.str, meca.NP1.dip, meca.NP2.str);
                if (meca.NP2.dip == 1000.) {
                    not_defined = TRUE;
                    transparence_old = transparence;
                    n_plane_old = n_plane;
                    transparence = TRUE;
                    n_plane = 1;
                    meca.NP1.rake = 1000.;
                    if (gmtdefs.verbose)
                        fprintf(stderr, "WARNING : second plane is not defined for event %s only first plane is plotted.\n", line);
                }
                else
                    meca.NP1.rake = computed_rake2(meca.NP2.str, meca.NP2.dip, meca.NP1.str, meca.NP1.dip, fault);
                meca.NP2.rake = computed_rake2(meca.NP1.str, meca.NP1.dip, meca.NP2.str, meca.NP2.dip, fault);
            }
	    
	    else if (read_axis) {
                T.val = atof(col[2+new]);
                T.str = atof(col[3+new]);
                T.dip = atof(col[4+new]);
                T.e = atoi(col[11+new]);

                N.val = atof(col[5+new]);
                N.str = atof(col[6+new]);
                N.dip = atof(col[7+new]);
                N.e = atoi(col[11+new]);

                P.val = atof(col[8+new]);
                P.str = atof(col[9+new]);
                P.dip = atof(col[10+new]);
                P.e = atoi(col[11+new]);
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

                if (transparence || plot_dc) axe2dc(T, P, &meca.NP1, &meca.NP2);
            }

	    else if (read_tensor) {
                for (i=2+new, n=0; i<8+new; i++, n++) mt.f[n] = atof(col[i]);
		/* if (mt.f[0] + mt.f[1] + mt.f[2] < GMT_SMALL)		What is SMALL here??
			fprintf (stderr, "psmeca WARNING: Found mechanism with zero trace != 0\n"); */
                mt.expo = atoi(col[i]);
/*
F. A. Dahlen and Jeroen Tromp, Theoretical Seismology, Princeton, 1998, p.167.
Definition of scalar moment.
*/
                meca.moment.mant = sqrt(squared(mt.f[0]) + squared(mt.f[1]) + squared(mt.f[2]) + 2.*(squared(mt.f[3]) + squared(mt.f[4]) + squared(mt.f[5]))) / M_SQRT2;
                meca.moment.exponent = mt.expo;
                meca.magms = 0.;

/* normalization by M0 */
                for(i=0;i<=5;i++) mt.f[i] /= meca.moment.mant;

                GMT_momten2axe(mt, &T, &N, &P);

                if (transparence || plot_dc) axe2dc(T, P, &meca.NP1, &meca.NP2);
            }

	    /* Common to all input types ... */

            GMT_geo_to_xy (xy[0], xy[1], &plot_x, &plot_y);

	    /* If option -C is used, read the new position */

            if (read_change_position) {
                xynew[ix] = atof(col[last-1+new]);
                xynew[iy] = atof(col[last+new]);
                if (fabs(xynew[ix]) > EPSIL || fabs(xynew[iy]) > EPSIL) {
                    GMT_setpen(&cpen);
                    GMT_geo_to_xy(xynew[0], xynew[1], &plot_xnew, &plot_ynew);
                    ps_circle(plot_x, plot_y, pointsize, fill.rgb, 1);
                    ps_plot(plot_x, plot_y, PSL_PEN_MOVE);
                    ps_plot(plot_xnew, plot_ynew, PSL_PEN_DRAW_AND_STROKE);
                    plot_x = plot_xnew;
                    plot_y = plot_ynew;
	    	}
            }

            if (one_size) {
                meca.moment.mant = 4.;
                meca.moment.exponent = 23;
            }

            moment.mant = meca.moment.mant;
            moment.exponent = meca.moment.exponent;
            size = (computed_mw(moment, meca.magms) / 5.) * scale;

            get_trans(xy[0], xy[1], &t11, &t12, &t21, &t22);
            delaz = atan2d(t12,t11);

            if ((read_axis || read_tensor) && !plot_dc) {

                T.str = zero_360(T.str + delaz);
                N.str = zero_360(N.str + delaz);
                P.str = zero_360(P.str + delaz);

                GMT_setpen(&lpen);
                if (fabs(N.val) < EPSIL && fabs(T.val + P.val) < EPSIL) {
                    axe2dc(T, P, &meca.NP1, &meca.NP2);
                    ps_mechanism(plot_x,plot_y,meca,size,fill.rgb,efill.rgb,outline);
                }
                else
                    ps_tensor(plot_x,plot_y,size,T,N,P,fill.rgb,efill.rgb,outline,plot_zerotrace);
            }

            if (tr_zerotrace) {
                GMT_setpen(&tz_pen);
                ps_tensor(plot_x,plot_y,size,T,N,P,nofill.rgb,nofill.rgb,tr_zerotrace,tr_zerotrace);
            }

            if (transparence) {
                meca.NP1.str = zero_360(meca.NP1.str + delaz);
                meca.NP2.str = zero_360(meca.NP2.str + delaz);
                GMT_setpen(&tr_pen);
                switch (n_plane) {
                    case 0 : 
                        ps_meca(plot_x, plot_y, meca, size);
                        break;
                    default :
                        ps_plan(plot_x, plot_y, meca, size, n_plane);
                        break;
                }
                if (not_defined) {
                    not_defined = FALSE;
                    transparence = transparence_old;
                    n_plane = n_plane_old;
                }
            } else if (read_aki || read_cmt || read_planes || plot_dc) {
                meca.NP1.str = zero_360(meca.NP1.str + delaz);
                meca.NP2.str = zero_360(meca.NP2.str + delaz);
                GMT_setpen(&lpen);
                ps_mechanism(plot_x,plot_y,meca,size,fill.rgb,efill.rgb,outline);
            }                                        

            if (!no_label) {
                GMT_setpen(&pen);
                switch (justify) {
                    case 2 :
                        if (draw_box)
                            ps_rect(plot_x - size * 0.5, plot_y + size * 0.5 + offset + (fontsize / 72.0), plot_x + size * 0.5, plot_y + size * 0.5 + offset, bfill.rgb, FALSE);
                        ps_text(plot_x, plot_y + size * 0.5 + offset, fontsize, event_title, angle, justify, form);
                        break;
                    case 10 :
                        if (draw_box)
                            ps_rect(plot_x - size * 0.5, plot_y - size * 0.5 - offset - (fontsize / 72.0), plot_x + size * 0.5, plot_y - size * 0.5 - offset, bfill.rgb, FALSE);
                        ps_text(plot_x, plot_y - size * 0.5 - offset, fontsize, event_title, angle, justify, form);
                        break;
                }
            }

            if (plot_axis) {
                if (read_axis || read_tensor)
                    axis2xy(plot_x, plot_y, size, P.str, P.dip, T.str, T.dip, &P_x, &P_y, &T_x, &T_y);
                else
                    ps_pt_axis(plot_x, plot_y, meca, size, &P.str, &P.dip, &T.str, &T.dip, &P_x, &P_y, &T_x, &T_y);
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
        if (fp != stdin) GMT_fclose (fp);
    }
    
    if (gmtdefs.verbose) 
        fprintf (stderr, "psmeca: Number of records read: %li\n", n_rec);


    if (skip_if_outside) GMT_map_clip_off ();
    
    ps_setpaint (gmtdefs.basemap_frame_rgb);

    if (pen.texture) ps_setdash (CNULL, 0);
    
    GMT_map_basemap ();

    GMT_plotend ();
    
    GMT_end (argc, argv);

    exit (EXIT_SUCCESS);
}
