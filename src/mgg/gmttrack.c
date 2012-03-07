/*--------------------------------------------------------------------
 *	$Id: gmttrack.c,v 1.40 2011/07/11 19:22:06 guru Exp $
 *
 *    Copyright (c) 1991-2011 by P. Wessel and W. H. F. Smith
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * gmttrack will read *.gmt-files and and write PostScript code
 * that will create a navigation plot on a basemap using the
 * projection specified by the user. WESN must be specified on the 
 * command line along with other options for scaling, annotation etc.
 *
 * Author:	Paul Wessel
 * Date:	10-MAR-1988
 * Revised:	22-NOV-1988	Color Version
 *		13-JUN-1991	v2.0
 *		10-MAY-1995	v3.0b
 *		16-MAR-1999	v3.2, no longer -F
 *		30-JUN-1999	v3.3.1, Forgot to place INCH in ps_plotinit
 *				No mention of -W in explanation.
 *		16-FEB-2000	v3.3.4, Only write to PS file via pslib
 *		24-JUL-2001	v3.5, Added -N to not clip leg names
 *
 */
 
#include "gmt.h"
#include "pslib.h"
#include "gmt_mgg.h"

#define DST 0.25
#define ANSIZE 0.125
#define MPRDEG 111194.9

struct ANOT {
	int annot_int_dist;
	int annot_int_time;
	int tick_int_dist;
	int tick_int_time;
};

struct LEG_ANOT {	/* Structure used to annotate legs after clipping is terminated */
	double x, y;
	int rec;
	char text[16];
};

int main (int argc, char **argv)
{
	double west = 0.0, east = 0.0, south = 0.0, north = 0.0, lat, lon, angle;
	double last_lon = 0.0, last_lat = 0.0, annotsize = ANSIZE, factor = 0.001, dlon;
	double dx, dy, x, y, x0, y0, dist = 0.0, start_dist = 0.0, stop_dist = 1.0E100; 
	
	int n_records, leg_year;
	
	int i, iw, ie, is, in, annot_dist = 0, tick_dist = 0, annot_time = 0, tick_time = 0;
	int nfiles = 0, rec, pen_type = 0, k, white[3], black[3];
	int error = FALSE, first, last = FALSE, annot_name = FALSE, both = FALSE, no_clip = FALSE;
	int  annot_tick = FALSE, draw_tick = FALSE, start_time = 0, stop_time = 2000000000;
	int this_julian, last_julian, year, month, day, hour, minute, second, n_alloc_c = GMT_SMALL_CHUNK;
	int mon1, day1 = 0, year1, hour1, min1, mon2, day2 = 0, year2, hour2, min2, n_id = 0;
        
	char agency[10], gmtfile[BUFSIZ], label[50], comment[BUFSIZ];
	
	struct LEG_ANOT *cruise_id = VNULL;
	struct GMT_PEN pen;
	struct ANOT info;
	struct GMTMGG_REC *record = NULL;
	int *distance = NULL;
	FILE *fp = NULL;
	
	struct GMTMGG_TIME *gmt = NULL;
	
	double heading (int rec, int n_records, struct GMTMGG_REC *record);
	int get_annotinfo (char *args, struct ANOT *info);
	int gmttrack_outside (int lon, int lat, int w, int e, int s, int n);
	void annot_legname (double x, double y, int rec, int nrecs, char *text, double size, int w, int e, int s, int n, struct GMTMGG_REC *record, struct GMT_PEN *pen);

	GMT_init_pen (&pen, GMT_PENWIDTH);
	for (i = 0; i < 3; i++) black[i] = 0, white[i] = 255;
	
	argc = (int)GMT_begin (argc, argv);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {
				case 'B':
				case 'P':
				case 'O':
				case 'K':
				case 'R':
				case 'J':
				case 'U':
				case 'V':
				case 'X':
				case 'x':
				case 'Y':
				case 'y':
				case 'c':
				case '\0':
					error += (int)GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;
				case 'W':
					GMT_getpen (&argv[i][2], &pen);
					break;
				case 'A':
					annot_name = TRUE;
					if (argv[i][2]) annotsize = atof (&argv[i][2]) * GMT_u2u[GMT_PT][GMT_INCH];
					break;
				case 'M':
					error += get_annotinfo (&argv[i][2], &info);
					break;
				case 'N':
					no_clip = TRUE;
					break;
				case 'D':               /* Assign start/stop times for sub-section */
					if (argv[i][2] == 'a') {        /* Start date */
						sscanf(&argv[i][3], "%d/%d/%d/%d:%d", &mon1, &day1, &year1, &hour1, &min1);
					}
					else if (argv[i][2] == 'b')      {      /* Stop date */
						sscanf(&argv[i][3], "%d/%d/%d/%d:%d", &mon2, &day2, &year2, &hour2, &min2);
					}
					else
					error = TRUE;
					break;
					
				case 'S':               /* Assign start/stop position for sub-section */
					if (argv[i][2] == 'a')  /* Start position */
						start_dist = atof(&argv[i][3]) * 1000.0;	/* km -> m */
					else if (argv[i][2] == 'b')     /* Stop position */
						stop_dist = atof(&argv[i][3]) * 1000.0;		/* km -> m */
					else
						error = TRUE;
					break;
			}
		}
		else
			nfiles++;
	}
	if (start_dist > stop_dist || start_time > stop_time) error = TRUE;
	if ((day1 > 0 && start_dist > 0.) || (day2 > 0 && stop_dist < 1.0e100)) error = TRUE;
	if (nfiles == 0) error = TRUE;
	if (GMT_check_rgb (pen.rgb)) error = TRUE;
		 
	if (error) {
		fprintf(stderr,"usage: gmttrack leg(s) %s %s [-A[size]] [%s]\n", GMT_Rgeo_OPT, GMT_J_OPT, GMT_B_OPT);
		fprintf(stderr,"	[-Da<startdate>] [-Db<stopdate>] [-K] [-M<trackticks>] [-N] [-O] [-P] [-Sa<startdist>]\n");
		fprintf(stderr,"	[-Sb<stopdist>] [%s] [-V] [-W<pen>] [%s] [%s] [%s]\n\n", GMT_U_OPT, GMT_X_OPT, GMT_Y_OPT, GMT_c_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "	Give one or more names of .gmt files\n");
		GMT_explain_option ('J');
		GMT_explain_option ('R');
		fprintf(stderr,"	OPTIONS:\n\n");
		fprintf (stderr, "	-A will Annotate legs when they enter the grid. <size> is optional character size in points [9]\n");
		GMT_explain_option ('B');
		fprintf(stderr,"	-Da<date> only plots from date (given as mm/dd/yyyy/hh:mm)\n");
		fprintf(stderr,"	-Db<date> only plots up to date (given as mm/dd/yyyy/hh:mm)\n");
		GMT_explain_option ('K');
		fprintf (stderr, "	-M to put time/distance Marks on the track. E.g. a500ka24ht6h means (a)nnotate\n");
		fprintf (stderr, "	   every 500 km (k) and 24 h(ours), with (t)ickmarks every 500 km and 6 (h)ours\n");
		fprintf (stderr, "	-N Do Not clip leg name annotation that fall outside map border [Default will clip]\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf(stderr,"        -Sa<dist> only plots from dist (in km)\n");
		fprintf(stderr,"        -Sb<dist> only plots up to dist (in km)\n");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		fprintf (stderr, "	-W sets track pen attributes [width = %gp, color = (%d/%d/%d), texture = solid line].\n", 
			pen.width, pen.rgb[0], pen.rgb[1], pen.rgb[2]);
		GMT_explain_option ('X');
		GMT_explain_option ('c');
		GMT_explain_option ('.');

		exit (EXIT_FAILURE);
	}
	
	gmtmggpath_init();

	if (east < west) west -= 360.0;
		
	GMT_err_fail (GMT_map_setup (west, east, south, north), "");
	
	GMT_plotinit (argc, argv);
	
	ps_setpaint (gmtdefs.basemap_frame_rgb);
	
	GMT_setpen (&pen);
	distance = (int *) GMT_memory (VNULL, (size_t)1, sizeof (int), "gmttrack");
	record = (struct GMTMGG_REC *) GMT_memory (VNULL, (size_t)1, sizeof (struct GMTMGG_REC), "gmttrack");

	ps_setfont (7);	/* Times Italic */
	GMT_map_clip_on (GMT_no_rgb, 3);
	iw = irint (west * 1.0e6);	ie = irint (east * 1.0e6);	is = irint (south * 1.0e6);	in = irint (north * 1.0e6);
	both = (info.annot_int_time && info.annot_int_dist);
	ps_setpaint (pen.rgb);
	
	if (no_clip) {
		cruise_id = (struct LEG_ANOT *) GMT_memory (VNULL, (size_t)n_alloc_c, sizeof (struct LEG_ANOT), "gmttrack");
	}
	
	for (i = 1; i < argc; i++) {	/* Loop over arguments and open/plot the gmtfiles */
		if (argv[i][0] == '-') continue;
  		if (gmtmggpath_func(gmtfile,argv[i])) {
   			fprintf(stderr,"gmttrack : Cannot find leg %s\n", argv[i]);
     			continue;
  		}
		if ((fp = fopen(gmtfile, "rb")) == NULL) {
			fprintf(stderr, "gmttrack: Could not open file %s\n", gmtfile);
			continue;
		}
		
		/* Read first record of file containing start-year, n_records and info */
		
		if (fread((void *)(&leg_year), (size_t)4, (size_t)1, fp) != 1) {
			fprintf(stderr,"gmttrack: Error while reading first year\n");
			exit (EXIT_FAILURE);
		}
		if (fread((void *)(&n_records), (size_t)4, (size_t)1, fp) != 1) {
			fprintf(stderr,"gmttrack: Error while reading no of records\n");
			exit (EXIT_FAILURE);
		}
		if (fread((void *)agency, (size_t)10, (size_t)1, fp) != 1) {
			fprintf(stderr,"gmttrack: Error while reading info-header\n");
			exit (EXIT_FAILURE);
		}
		
		distance = (int *) GMT_memory ((void *)distance, (size_t)n_records, sizeof (int), "gmttrack");
		record = (struct GMTMGG_REC *) GMT_memory ((void *)record, (size_t)n_records, sizeof (struct GMTMGG_REC), "gmttrack");
		
		gmt = gmtmgg_init(leg_year);	/* Initialize gmt_structure */
		last_julian = -1;
                /* Decode date to time in sec if needed */

		if (day1 > 0) gmtmgg_time (&start_time, year1, mon1, day1, hour1, min1, 0, gmt);
		if (day2 > 0) gmtmgg_time (&stop_time, year2, mon2, day2, hour2, min2, 0, gmt);

		/* Read cruise data */
		
		for (rec = k = 0; rec < n_records; rec++) {
			if (fread((void *)(&record[k]), (size_t)18, (size_t)1, fp) != 1) {
				fprintf(stderr,"gmttrack: Error reading data record no %d\n",rec);
				exit (EXIT_FAILURE);
			}
			lon = record[k].lon * MDEG2DEG;
			lat = record[k].lat * MDEG2DEG;
			if (rec == 0) {
				dist = 0;
				last_lon = lon;
				last_lat = lat;
			}
			else {
				dlon = fabs (lon - last_lon);
				if (dlon > 180.0) dlon = 360.0 - dlon;
				dx = dlon * cosd (0.5 * (lat+last_lat));
				dy = (lat - last_lat);
				last_lon = lon;
				last_lat = lat;
				dist += hypot (dx, dy) * MPRDEG;	/* Dist in meters */
			}
			distance[k] = irint (dist);
			if (dist < start_dist || dist > stop_dist) continue;
			if (record[k].time < start_time || record[k].time > stop_time) continue;
			k++;
		}
		fclose (fp);
		n_records = k;
		
		if (gmtdefs.verbose) fprintf (stderr, "gmttrack: Plotting %s\n", agency);
		sprintf (comment, "Tracking %s", agency);
		ps_comment (comment);
		
		/* First draw the track line, clip segments outside the area */
		
		first = TRUE;
		last = FALSE;
		for (rec = 0; rec < n_records; rec++) {
			if (gmttrack_outside (record[rec].lon, record[rec].lat, iw, ie, is, in)) {
				if (last && rec < (n_records-1)) {
					lon = record[rec].lon * 1.0e-6;	lat = record[rec].lat * 1.0e-6;
					GMT_geo_to_xy (lon, lat, &x, &y);
					ps_plot (x, y, 2);
					last = FALSE;
				}
				first = TRUE;
				continue;
			}
			lon = record[rec].lon * 1.0e-6;	lat = record[rec].lat * 1.0e-6;
			GMT_geo_to_xy (lon, lat, &x, &y);
			if (first) {
				if (rec > 0) {
					lon = record[rec].lon * 1.0e-6;	lat = record[rec].lat * 1.0e-6;
					GMT_geo_to_xy (lon, lat, &x0, &y0);
					ps_plot (x0, y0, 3);
				}
				else
					ps_plot (x, y, 3);
				if (annot_name) {
					if (no_clip) {	/* Keep these in a list to plot after clipping is turned off */
						cruise_id[n_id].x = x;
						cruise_id[n_id].y = y;
						cruise_id[n_id].rec = rec;
						strcpy (cruise_id[n_id].text, argv[i]);
						n_id++;
						if (n_id == n_alloc_c) {
							n_alloc_c <<= 1;
							cruise_id = (struct LEG_ANOT *) GMT_memory ((void *)cruise_id, (size_t)n_alloc_c, sizeof (struct LEG_ANOT), "gmttrack");
						}
					}
					else
						annot_legname (x, y, rec, n_records, argv[i], GMT_u2u[GMT_INCH][GMT_PT] * 1.25 * annotsize, iw ,ie, is, in, record, &pen);
				}
				first = FALSE;
				if (info.annot_int_dist > 0) annot_dist = (distance[rec] / info.annot_int_dist + 1) * info.annot_int_dist;
				if (info.tick_int_dist > 0) tick_dist = (distance[rec] / info.tick_int_dist + 1) * info.tick_int_dist;
				if (info.annot_int_time > 0) annot_time = (record[rec].time / info.annot_int_time + 1) * info.annot_int_time;
				if (info.tick_int_time > 0) tick_time = (record[rec].time / info.tick_int_time + 1) * info.tick_int_time;
			}
			ps_plot (x, y, 2);
			last = TRUE;
			
			/* See if we need to annotate/tick the trackline for time/km marks */
			
			if (info.annot_int_time && (record[rec].time >= annot_time)) {
				annot_time += info.annot_int_time;
				annot_tick = 1;
			}
			if (info.annot_int_dist && (distance[rec] >= annot_dist)) {
				annot_dist += info.annot_int_dist;
				annot_tick += 2;
			}
			if (info.tick_int_time && (record[rec].time >= tick_time)) {
				tick_time += info.tick_int_time;
				draw_tick = 1;
			}
			if (info.tick_int_dist && (distance[rec] >= tick_dist)) {
				tick_dist += info.tick_int_dist;
				draw_tick += 2;
			}
			if (annot_tick) {
				angle = heading (rec, n_records, record);
				if (angle < 0.0)
					angle += 90.0;
				else
					angle -= 90.0;
				if (annot_tick & 1) {	/* Time mark */
					this_julian = gmtmgg_date (annot_time - info.annot_int_time,&year,&month,&day,&hour,&minute,&second,gmt);
					ps_command ("S");
					if (this_julian != last_julian) {
						sprintf (label, " %.2d:%.2d %d/%d", hour, minute, month, day);
						ps_circle (x, y, 0.2*annotsize, black, TRUE);
						ps_setfont (7);
					}
					else {
						sprintf (label, " %.2d:%.2d", hour, minute);
						ps_circle (x, y, 0.2*annotsize, white, TRUE);
					}
					ps_setpaint (gmtdefs.basemap_frame_rgb);
					ps_text (x, y, GMT_u2u[GMT_INCH][GMT_PT] * annotsize, label, angle, 5, 0);
					if (this_julian != last_julian) ps_setfont (6);
					last_julian = this_julian;
				}
				if (annot_tick & 2) {	/* Distance mark */
					ps_command ("S");
					sprintf (label, "%d km  ", (int)((annot_dist - info.annot_int_dist) * factor));
					ps_square (x, y, 0.4*annotsize, gmtdefs.basemap_frame_rgb, FALSE);
					ps_setpaint (gmtdefs.basemap_frame_rgb);
					ps_text (x, y, GMT_u2u[GMT_INCH][GMT_PT] * annotsize, label, angle, 7, 0);
				}
			}
			if (both && !(annot_tick & 1) && (draw_tick & 1)) {
				ps_command ("S");
				ps_circle (x, y, 0.2*annotsize, gmtdefs.basemap_frame_rgb, FALSE);
			}
			if (both && !(annot_tick & 2) && (draw_tick & 2)) {
				ps_command ("S");
				ps_square (x, y, 0.4*annotsize, gmtdefs.basemap_frame_rgb, FALSE);
			}
			if (draw_tick) {
				ps_command ("S");
				ps_setpaint (gmtdefs.basemap_frame_rgb);
				ps_cross (x, y, 0.4*annotsize);
			}
			if (annot_tick || draw_tick) {
				ps_setpaint (pen.rgb);
				ps_plot (x, y, 3);
				annot_tick = draw_tick = FALSE;
			}
		}
	}
	GMT_map_clip_off();
	if (pen_type > 0) ps_setdash (NULL, 0);

	GMT_map_basemap ();

	if (annot_name && no_clip) {	/* Plot leg names after clipping is terminated ( see -N) */
		int id;
		double size;
		size = GMT_u2u[GMT_INCH][GMT_PT] * 1.25 * annotsize;
		for (id = 0; id < n_id; id++) annot_legname (cruise_id[id].x, cruise_id[id].y, cruise_id[id].rec, n_records, cruise_id[id].text, size, iw ,ie, is, in, record, &pen);
		GMT_free ((void *)cruise_id);
	}

	GMT_plotend ();
	GMT_free ((void *)record);
	GMT_free ((void *)distance);
	
	gmtmgg_end ();
	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}


double heading (int rec, int n_records, struct GMTMGG_REC *record)
{
	int i1, i2;
	double angle, lon, lat, x1, x2, y1, y2;
	
	i1 = rec - 25;
	if (i1 < 0) i1 = 0;
	i2 = i1 + 25;
	if (i2 > (n_records-1)) i2 = n_records - 1;
	lon = record[i2].lon * 1.0e-6;	lat = record[i2].lat * 1.0e-6;
	GMT_geo_to_xy (lon, lat, &x2, &y2);
	lon = record[i1].lon * 1.0e-6;	lat = record[i1].lat * 1.0e-6;
	GMT_geo_to_xy (lon, lat, &x1, &y1);
	angle = atan2d (y2 - y1, x2 - x1);
	if (angle > 90.)
		angle -= 180;
	else if (angle < -90.0)
		angle += 180.0;
	return (angle);
}

int get_annotinfo (char *args, struct ANOT *info)
{
	int i1, i2, error = FALSE;
	int flag1, flag2, type;
	double value;
	
	info->annot_int_dist = info->tick_int_dist = 0;
	info->annot_int_time = info->tick_int_time = 0;

	i1 = 0;
	while (args[i1]) {
		flag1 = 'a';
		if (isalpha ((int)args[i1])) {
			flag1 = args[i1];
			if (isupper(flag1)) flag1 = tolower (flag1);
			i1++;
		}
		i2 = i1;
		while (args[i2] && strchr ("athkm", (int)args[i2]) == NULL) i2++;
		value = atof (&args[i1]);
		flag2 = args[i2];
		if (isupper(flag2)) flag2 = tolower (flag2);
		if (flag2 == 'h') {		/* Hours */
			value *= 3600;
			type = 't';
		}
		else if (flag2 == 'k') {	/* kilometers */
			value *= 1000;
			type = 'd';
		}
		else if (flag2 == 'm') {	/* Minutes */
			value *= 60;
			type = 't';
		}
		else				/* Default is seconds */
			type = 't';
		i2++;
		if (flag1 == 'a') {	/* Annotation interval */
			if (type == 'd')	/* Distance */
				info->annot_int_dist = (int)value;
			else
				info->annot_int_time = (int)value;
		}
		else {			/* Tickmark interval */
			if (type == 'd')	/* Distance */
				info->tick_int_dist = (int)value;
			else
				info->tick_int_time = (int)value;
		}
		i1 = i2;
	}
	if (info->annot_int_dist <= 0 && info->tick_int_dist <= 0 && info->annot_int_time <= 0 && info->tick_int_time <= 0)
		error = TRUE;
	if (info->annot_int_dist <= 0)
		info->annot_int_dist = info->tick_int_dist;
	else if (info->tick_int_dist <= 0)
		info->tick_int_dist = info->annot_int_dist;
	if (info->annot_int_time <= 0)
		info->annot_int_time = info->tick_int_time;
	else if (info->tick_int_time <= 0)
		info->tick_int_time = info->annot_int_time;
	return (error);
}

int gmttrack_outside (int lon, int lat, int w, int e, int s, int n)
{
	if (lat < s || lat > n) return (TRUE);
	lon -= 360000000;
	while (lon < w) lon += 360000000;
	if (lon > e) return (TRUE);
	return (FALSE);
}

void annot_legname (double x, double y, int rec, int nrecs, char *text, double size, int w, int e, int s, int n, struct GMTMGG_REC *record, struct GMT_PEN *pen)
{
	int just;
	double angle;
	
	if (rec == 0) {
		angle = 0.0;
		just = 6;
	}
	else {
		angle = heading (rec, nrecs, record);
		if (record[rec-1].lat < s)
			just = (angle >= 0.0) ? 1 : 3;
		else if (record[rec-1].lat > n)
			just = (angle >= 0.0) ? 11 : 9;
		else if (record[rec-1].lon < w)
			just = (angle >= 0.0) ? 9 : 1;
		else
			just = (angle >= 0.0) ? 3 : 11;
	}
	ps_setfont (5);
	ps_setpaint (gmtdefs.basemap_frame_rgb);
	ps_text (x, y, size, text, angle, just, 0);
	ps_setpaint (pen->rgb);
	ps_plot (x, y, 3);
	ps_setfont (7);
}
