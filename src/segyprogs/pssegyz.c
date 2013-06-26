/*--------------------------------------------------------------------
 *	$Id: pssegyz.c 9923 2012-12-18 20:45:53Z pwessel $
 *
 *    Copyright (c) 1999-2013 by T. Henstock
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/* pssegyz program to plot segy files in 3d in postscript with variable trace spacing option
 * uses routines from the GMT pslib (the postscript imagemask for plotting a
 * 1-bit depth bitmap which will not obliterate stuff underneath!
 *
 * Author:	Tim Henstock (then@noc.soton.ac.uk)
 * Date:	17-Nov-97
 * Version:	1.0
 *
 * heavily modified from pssegy version 1.2
 *
 * Bug fixes:	1.1, 11-6-96 remove undesirable normalization introduced by -Z option. Add -U option for reduction.
 *
 * add option for colored bitmap (modified psimagemask as well..., old version will segfault probably)
 * 
 * enhancements: 1.2 , 1/7/99 check number of samples trace by trace to cope with SEGY with variable trace length
 * NB that -T option from pssegy is _not_ transferred
 *
 * 		2.0, 6/7/99 update for GMT v 3.3.1
 *
 *              2.1 10/4/2001 fix unreduce bug, modify to byte-swap integers in the headers, make 64-bit safe
 *
 *              2.2 25/2/2002 fix bug with reduced data plotting, improve error checking on parameters
 *
 * This program is free software and may be copied, modified or redistributed
 * under the terms of the GNU public license, see http://www.gnu.org
 */

 
#include "gmt.h"
#include "pslib.h"
#include "segy_io.h"
 
/* internal function calls */
	double rms( float *data, int nsamp);
	void wig_bmap(double x0, double y0, float data0, float data1, double z0, double z1, double dev_x, double dev_y); 
	void shade_bmap( double x0, double y0, float data0, float data1, double z0, double z1, int negative, double dev_x, double dev_y);
	void plot_trace(float *data, double dz, double x0, double y0, int n_samp, int do_fill, int negative, int plot_wig, float toffset, double dev_x, double dev_y);
	void shade_tri (double apex_x, double apex_y, double edge_y, double slope, double slope0);
	void shade_quad (double x0, double y0, double x1, double y_edge, double slope1, double slope0);
	int paint(int ix, int iy);


unsigned char bmask[8]={128, 64, 32, 16, 8, 4, 2, 1};
unsigned char *bitmap;
int bm_nx, bm_ny;


int main (int argc, char **argv)
{
	GMT_LONG error = FALSE;
	int xplot_cdp = FALSE, xplot_offset = FALSE, fix_x = FALSE, byte_x = 0; /* parameters for location plot */
	int yplot_cdp = FALSE, yplot_offset = FALSE, fix_y = FALSE, byte_y = 0;
	int doclip = FALSE;
	int normalize = FALSE, do_fill = FALSE, negative = FALSE, plot_wig = FALSE;
	int no_zero = FALSE, polarity=1;
#ifdef WORDS_BIGENDIAN
	int swap_bytes = FALSE;
#else
	int swap_bytes = TRUE;
#endif

	int i, nm;
	int ix, iz, n_traces=10000, n_samp=0, n_sampr=0, shade[3]={0,0,0}, trans[3]={-1,-1,-1};

	int check;

	double w, e, s, n, dx = 1.0, dz = 0.0; /* dx, dz are trace and sample interval */
	double xlen, ylen, xpix, ypix;
	double x0 = 0.0 , y0 = 0.0;


	float bias = 0.0, scale = 1.0, clip = 0.0, dev_x = 0.0, dev_y = 0.0;
	double redvel=0.0;
	float toffset=0.0;

	char input[512];
	char *text = NULL;
	char *head = NULL;
	int head2;


	char reelhead[3200];
	float *data = NULL;
	SEGYHEAD *header = NULL;
	SEGYREEL binhead;

	FILE *fpi = NULL;


	input[0] = 0;
	w = e = s = n = 0.0;

	argc = (int)GMT_begin (argc, argv);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'V':
				case 'R':
				case 'J':
				case 'O':
				case 'K':
				case 'P':
				case '\0':
					error += GMT_parse_common_options (argv[i], &w, &e, &s, &n);
					break;
                                case 'E':
					error += GMT_get_proj3D (&argv[i][2], &z_project.view_azimuth, &z_project.view_elevation);
					break;
				/* parameters for wiggle mode */
				case 'F':
					do_fill = TRUE;
					if (GMT_getrgb (&argv[i][2], shade)) {
						++error;
						GMT_rgb_syntax ('F', " ");
					}
					break;
				case 'I':
					negative = TRUE;
					break;
				case 'W':
					plot_wig = TRUE;
					break;
				/* trace norm., clip and bias */
				case 'N':
					normalize = TRUE;
					break;
				case 'C':
					doclip = TRUE;
					clip = (float) atof (&argv[i][2]);
 					break;
				case 'B':
					bias = (float) atof (&argv[i][2]);
					break;
				case 'Z':
					no_zero = TRUE;
					break;
				/* variable spacing */
				case 'S':
				  if ((text = strstr(&argv[i][2], "/")) != NULL) {
				    	text = strtok(&argv[i][2], "/");
				    	if ((text[0] >='0' && text[0]<='9') || text[0] == '-'){
				      		fix_x = TRUE;
				      	x0 = (double) atof(text);
				    	}
				    	else{
						switch(text[0]){
							case 'o':
								xplot_offset = TRUE;
								break;
							case 'c':
								xplot_cdp = TRUE;
								break;
							case 'b':
					  			byte_x = atoi(&text[1]);
					  			break;
						}
				    	} /* now parameters for y */
				  	text = strtok(CNULL, "/");
				  	if (text != NULL){
				    		if ((text[0] >='0' && text[0]<='9') || text[0] == '-'){
				      			fix_y = TRUE;
				      			y0 = (double) atof(text);
				    		}
				    		else{
							switch(text[0]){
								case 'o':
									yplot_offset = TRUE;
									break;
								case 'c':
									yplot_cdp = TRUE;
									break;
								case 'b':
							  	byte_y = atoi(&text[1]);
							  	break;
							}
				    		}
				  	}
				  	else{
				    		fprintf(stderr,"%s:  Must specify parameters for x and y \n", GMT_program);
				    		++error;
				  	}
				}
				else{
					fprintf(stderr,"%s:  Must specify parameters for x and y \n", GMT_program);
					++error;
			 	}
				/*error += (!xplot_cdp && !xplot_offset && !byte_x && !fix_x);
				error += (!yplot_cdp && !yplot_offset && !byte_x && !fix_y);*/
				break;
				/* trace scaling */
				case 'D':
				  	if ((text = strstr(&argv[i][2], "/")) != NULL) {
				   	 	text = strtok(&argv[i][2], "/");
						dev_x = (float) atof (text);
				  		text = strtok(CNULL, "/");
						dev_y = (float) atof (text);
					}
					else
						dev_x = (float) atof (&argv[i][2]);
					/*error += ((dev_x < 0.0) && (dev_y < 0.0));*/
					error += (!dev_x && !dev_y);
			/*fprintf (stderr, "dev_x %f \t dev_y %f \n",dev_x, dev_y);*/
					break;
				/* over-rides of header info */
				case 'X': /* -X and -Y can be changed in gmt routines to lower case...*/
				case 'x':
					dx = atof (&argv[i][2]);
					break;
				case 'Y':
				case 'y':
					dz = atof (&argv[i][2]);
					break;
				case 'L':
					n_sampr = atoi (&argv[i][2]);
					break;
				case 'M':
					n_traces = atoi (&argv[i][2]);
					break;
				/* reduction velocity application */
				case 'U':
					redvel = atof (&argv[i][2]);
					break;
				case 'A':
					swap_bytes = !swap_bytes;
					break;
				default:
					error = TRUE;
					break;
			}
			}
		else if ((fpi = fopen (argv[i], "rb")) == NULL) {
			fprintf (stderr, "%s: Cannot find segy file %s\n",  GMT_program, argv[i]);
			exit (EXIT_FAILURE);
		}
	}


	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "pssegyz %s - Plot a segy file in PostScript\n\n", GMT_VERSION);
		fprintf (stderr, "usage: pssegyz [<segyfile>] %s %s -D<dev>\n", GMT_Jx_OPT, GMT_Rx_OPT);
		fprintf (stderr, "	[%s] [-C<clip>] [-B<bias>] [-N] [-Z]\n", GMT_E_OPT);
		fprintf (stderr, "	[-F<gray>|<r/g/b>] [-I] [-W] [-S<x>/<y>]\n");
		fprintf (stderr, "	[-X<dx>] [-Y<dz>] [-L<nsamp>] [-M<ntraces>] \n");
		fprintf (stderr, "	[-U<redvel>] [-A] [-O] [-K] [-P]\n");




		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, 
			"\n\t-Jx for projection.  Scale in INCH/units.  Specify one:\n\t -Jx<x-scale>              Linear projection\n\t-Jx<x-scale>l             Log10 projection\n\t  -Jx<x-scale>p<power>      x^power projection\n\tUse / to specify separate x/y scaling.\n\t If -JX is used then give axes lengths rather than scales\n\t regular map projections may not work!\n");
		GMT_explain_option ('R');
		fprintf (stderr, "	NB units for y are s or km\n");
		fprintf (stderr, "	-D<dev> to give deviation in X units of plot for 1.0 on scaled trace.\n");
		fprintf (stderr, "	<dev> is single number (applied equally in X and Y directions or <devX>/<devY>.\n");
		fprintf (stderr, "	 IEEE SEGY file [or standard input] \n\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_explain_option ('E');
		GMT_explain_option ('V');
		fprintf (stderr, "	-C<clip> to clip scaled trace excursions at <clip>, applied after bias\n");
		fprintf (stderr, "	-B<bias> to bias scaled traces (-B-0.1 subtracts 0.1 from values)\n");
		fprintf (stderr,"	-N to trace normalize the plot\n");
		fprintf (stderr,"		order of operations: [normalize][bias][clip](deviation)\n");
		fprintf (stderr,"	-Z to suppress plotting traces whose rms amplitude is 0 \n");
		fprintf (stderr,"	-F<gray>|<r/g/b> to fill variable area with shade\n");
		fprintf (stderr,"		only single color for the bitmap though\n");
		fprintf (stderr,"	-I to fill negative rather than positive excursions\n");
		fprintf (stderr,"	-W to plot wiggle trace\n");
		fprintf (stderr,"	must specify either -W or -F\n");
		fprintf (stderr,"	-X<mult> multiply trace locations by <mult>\n");
		fprintf (stderr,"	-Y<dz> to override sample interval\n");
		fprintf (stderr,"	-S<x/y> to set variable spacing\n");
		fprintf (stderr,"		x,y are (number) for fixed location, c for cdp, o for offset, b<n> for long int at byte n\n");
		fprintf (stderr,"	-L<nsamp> to override number of samples\n");
		fprintf (stderr,"	-M<ntraces> to fix number of traces. Default reads all traces.\n\t\t-M0 will read number in binary header, -Mn will attempt to read only n traces.\n");
		fprintf (stderr,"	-U<redvel> to apply reduction velocity (-ve removes reduction alreadz present)\n");
		fprintf (stderr,"	-A flips the default byte-swap state (default assumes data have a bigendian byte-order)\n");
		GMT_explain_option ('O');
		GMT_explain_option ('K');
		GMT_explain_option ('P');
		exit (EXIT_FAILURE);
	}

	if (negative && !do_fill){ /* negative with no fill */
		error++;
		fprintf(stderr,"%s: SYNTAX ERROR: Must specify -F with -I\n", GMT_program);
	}
	if (!do_fill && !plot_wig){ /* no plotting specified */
 		error++;
 		fprintf(stderr,"%s: SYNTAX ERROR: Must specify -F or -W\n", GMT_program);
 	}
 	if (dev_x < 0.0 || dev_y < 0.0){
 		error++;
 		fprintf(stderr,"%s: SYNTAX ERROR: Must specify a positive deviation\n",GMT_program);
 	}
 	if (!GMT_IS_LINEAR){
 		fprintf(stderr,"%s: WARNING: you asked for a non-rectangular projection. \n It will probably still work, but be prepared for problems\n",GMT_program);
 	}
  	if (!project_info.z_bottom && !project_info.z_top){ /* no z range in the -R parameter */
 		error++;
 		fprintf(stderr,"%s: ERROR: must specify z range in -R option.",GMT_program);
 	}
  	if (!project_info.region_supplied){ /* no -R parameter */
 		error++;
 		fprintf(stderr,"%s: ERROR: must specify -R option.",GMT_program);
 	}
	if (z_project.view_azimuth > 360.0 || z_project.view_elevation <= 0.0 || z_project.view_elevation > 90.0) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Enter azimuth in 0-360 range, elevation in 0-90 range\n", GMT_program);
                error++;
        }

 	if ((xplot_cdp && xplot_offset) || (yplot_cdp && yplot_offset) || (xplot_cdp && byte_x) || (yplot_cdp && byte_y) || (byte_x && xplot_offset) || (byte_y && yplot_offset)){
 		fprintf(stderr,"%s: SYNTAX ERROR: Can't specify more than one trace location key\n",GMT_program);
 		error++;
 	}
 
 	if (error) exit (EXIT_FAILURE);

	if (fpi == NULL) fpi = stdin;

/* set up map projection and PS plotting */
	GMT_err_fail (GMT_map_setup (w, e, s, n), "");
	GMT_plotinit (argc, argv);

        /*if (project_info.three_D) ps_transrotate (-z_project.xmin, -z_project.ymin, 0.0);*/


/* define area for plotting and size of array for bitmap */
	xlen = z_project.xmax-z_project.xmin;
	xpix = xlen*gmtdefs.dpi; /* pixels in x direction */
	/*xpix /= 8.0;
	bm_nx = 1 +(int) xpix;*/
	bm_nx = (int) ceil (xpix/8.0); /* store 8 pixels per byte in x direction but must have
				whole number of bytes per scan */
	ylen = z_project.ymax-z_project.ymin;
	ypix = ylen*gmtdefs.dpi; /* pixels in y direction */
	bm_ny = (int) ypix;
	nm = bm_nx*bm_ny;


/* read in reel headers from segy file */
	if ((check = get_segy_reelhd (fpi, reelhead)) != TRUE) exit(1);
	if ((check = get_segy_binhd (fpi, &binhead)) != TRUE) exit(1);

	if(swap_bytes){
/* this is a little-endian system, and we need to byte-swap ints in the reel header - we only
use a few of these*/
		if (gmtdefs.verbose) fprintf(stderr, "%s: swapping bytes for ints in the headers\n", GMT_program);
		binhead.num_traces = GMT_swab2(binhead.num_traces);
		binhead.nsamp = GMT_swab2(binhead.nsamp);
		binhead.dsfc = GMT_swab2(binhead.dsfc);
		binhead.sr = GMT_swab2(binhead.sr);
	}

/* set parameters from the reel headers */
	if (!n_traces)
		n_traces = binhead.num_traces;

	if (gmtdefs.verbose) fprintf(stderr, "%s: Number of traces in header is %d\n", GMT_program, n_traces);


	if (!n_sampr){/* number of samples not overridden*/
		n_sampr = binhead.nsamp;
		fprintf(stderr,"%s: Number of samples per trace is %d\n", GMT_program, n_sampr);
	}
	else if ((n_sampr != binhead.nsamp) && (binhead.nsamp))
		fprintf(stderr,"%s: warning nsampr input %d, nsampr in header %d\n", GMT_program, n_sampr, binhead.nsamp);

	if (!n_sampr){ /* no number of samples still - a problem! */
		fprintf(stderr, "%s: Error, number of samples per trace unknown\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if(gmtdefs.verbose) 
		fprintf(stderr, "%s: Number of samples is %dl\n", GMT_program, n_samp);

	if(binhead.dsfc != 5) fprintf(stderr, "%s: WARNING data not in IEEE format\n", GMT_program);

	if (!dz){
		dz = binhead.sr; /* sample interval of data (microseconds) */
		dz /= 1000000.0;
		fprintf(stderr,"%s: Sample interval is %f s\n", GMT_program, dz);
	}
	else if ((dz != binhead.sr) && (binhead.sr)) /* value in header overridden by input */
		fprintf(stderr, "%s: Warning dz input %f, dz in header %f\n", GMT_program, dz, (float)binhead.sr);

	if (!dz){ /* still no sample interval at this point is a problem! */
		fprintf(stderr, "%s: Error, no sample interval in reel header\n", GMT_program);
		exit (EXIT_FAILURE);
	}


	bitmap = (unsigned char *) GMT_memory (NULL, (size_t)nm, sizeof (char), "pssegyz");

	ix=0;
	while ((ix<n_traces) && (header = get_segy_header(fpi))){   /* read traces one by one */
	/* check true location header for x */
		if (xplot_offset){ /* plot traces by offset, cdp, or input order */
			int32_t offset = ((swap_bytes)? GMT_swab4(header->sourceToRecDist): header->sourceToRecDist);
			x0 = (double) offset;
		}
		else if (xplot_cdp){
			int32_t cdpval = ((swap_bytes)? GMT_swab4(header->cdpEns): header->cdpEns);
			x0 = (double) cdpval;
		}
		else if (byte_x){ /* ugly code - want to get value starting at byte_x of header into a double... */
		   	head = (char *) header;
			memcpy(&head2, &head[byte_x], 4); /* edited to fix bug where 8bytes were copied from head.
                                                Caused by casting to a long directly from char array*/ 
			x0 = (double) ((swap_bytes)? GMT_swab4(head2): head2);
		}
		else if (fix_x)
			x0 /= dx;
		else
			x0 = (1.0 + (double) ix); /* default x to input trace number */

	/* now do same for y */
		if (yplot_offset){ /* plot traces by offset, cdp, or input order */
			int32_t offset = ((swap_bytes)? GMT_swab4(header->sourceToRecDist): header->sourceToRecDist);
			y0 = (double) offset;
		}
		else if (yplot_cdp){
			int32_t cdpval = ((swap_bytes)? GMT_swab4(header->cdpEns): header->cdpEns);
			y0 = (double) cdpval;
		}
		else if (byte_y){
			head =  (char *) header;
			memcpy(&head2, &head[byte_y], 4); /* edited to fix bug where 8bytes were copied from head.
                                                Caused by casting to a long directly from char array*/ 
			y0 = (double) ((swap_bytes)? GMT_swab4(head2): head2);
		}
		else if (fix_y)
			y0 /= dx;
		else
			y0 = s / dx; /* default y to s edge of projection */

		x0 *= dx;
		y0 *= dx; /* scale x and y by the input dx scalar */

		if (swap_bytes){
/* need to permanently byte-swap some things in the trace header 
do this after getting the location of where traces are plotted in case the general byte_x case
overlaps a defined header in a strange way */
			 header->sourceToRecDist=GMT_swab4(header->sourceToRecDist);
			 header->sampleLength=GMT_swab2(header->sampleLength);
			 header->num_samps=GMT_swab4(header->num_samps);
		}

		if (gmtdefs.verbose) 
			fprintf(stderr, "%s: trace %d at x=%f, y=%f \n", GMT_program, ix+1, x0, y0);

		if (redvel){
			toffset = (float) -(fabs((double)(header->sourceToRecDist))/redvel);
			if (gmtdefs.verbose)
				fprintf(stderr, "%s: time shifted by %f\n", GMT_program, toffset);
		}

		data = (float *) get_segy_data(fpi, header); /* read a trace */
	/* get number of samples in _this_ trace (e.g. OMEGA has strange ideas about SEGY standard)
                or set to number in reel header */
                if ( !(n_samp = samp_rd(header)) ) n_samp = n_sampr;

		if(swap_bytes){ /* need to swap the order of the bytes in the data even though assuming IEEE format */
			int *intdata = (int *) data;
			for (iz=0; iz<n_samp; iz++){
				intdata[iz]=GMT_swab4(intdata[iz]);
			}
		}

		if(normalize || no_zero){
			scale= (float) rms(data, n_samp);
			if (gmtdefs.verbose) 
				fprintf(stderr, "%s: \t\t rms value is %f\n", GMT_program, scale);
		}
		for (iz=0; iz<n_samp; iz++){ /* scale bias and clip each sample in the trace */
			if (normalize) data[iz] /= scale;
			data[iz] += bias; 
			if(doclip && (fabs(data[iz]) > clip)) data[iz] = (float)(clip*data[iz]/fabs(data[iz])); /* apply bias and then clip */
		}

		if (!no_zero || scale)
			plot_trace (data, dz, x0, y0, n_samp, do_fill, negative, plot_wig, toffset, dev_x, dev_y);
		free (data);
		free (header);
		ix++;
	}

	/* map_clip_on (-1, -1, -1, 3); */
	/* set a clip at the map boundary since the image space overlaps a little */
	ps_bitimage (0.0,0.0,xlen, ylen, bitmap, 8*bm_nx, bm_ny, polarity, shade, trans);

	/* map_clip_off ();*/

	if (fpi != stdin) fclose (fpi);

	/*ps_rotatetrans (z_project.xmin, z_project.ymin, 0.0);*/
	GMT_plotend ();

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}


double rms(float *data, int n_samp)
{/* function to return rms amplitude of n_samp values from the array data */
	int ix;
	double sumsq=0.0;

	for (ix=0; ix<n_samp; ix++){
		sumsq += ((double) data[ix])*((double) data[ix]);
	}
	sumsq /= ((double) n_samp);
	sumsq = sqrt (sumsq);
	return (sumsq);
}

void plot_trace(float *data, double dz, double x0, double y0, int n_samp, int do_fill, int negative, int plot_wig, float toffset, double dev_x, double dev_y) 
	/* shell function to loop over all samples in the current trace, determine plot options
	 * and call the appropriate bitmap routine */ 
{
int iz;
int paint_wiggle;
float z0 = (float)project_info.z_bottom, z1;

	for(iz=1; iz<n_samp; iz++){ 	/* loop over samples on trace - refer to pairs iz-1, iz */
		z1 = (float )(dz * (float) iz + toffset);
		if(z1 >= project_info.z_bottom && z1 <= project_info.z_top){ /* check within z bounds specified */
	/*		fprintf(stderr,"x0, %f\t y0, %f\t,z1, %f\t data[iz], %f\t iz, %d\n",x0,y0,z1,data[iz],iz); */
			if (plot_wig) /* plotting wiggle */
				wig_bmap (x0, y0, data[iz-1],data[iz], z0, z1, dev_x, dev_y);
			if (do_fill){ /* plotting VA -- check data points first */
				paint_wiggle = ( (!negative && ((data[iz-1]>=0.)||(data[iz]>=0.))) || (negative && ((data[iz-1]<=0.0)||(data[iz]<=0.0))) );
				if (paint_wiggle)
						shade_bmap (x0, y0, data[iz-1], data[iz], z0, z1, negative, dev_x, dev_y);
			}
			z0=z1;
		}
	}
}


void wig_bmap(double x0, double y0, float data0, float data1, double z0, double z1, double dev_x, double dev_y) /* apply current sample with all options to bitmap */ 
{
double xp0, xp1, yp0, yp1, slope;
int px0, px1, py0, py1, ix, iy;


	GMT_geoz_to_xy (x0+(double)data0*dev_x, y0+(double)data0*dev_y, z0, &xp0, &yp0); /* returns 2 ends of line segment in plot coords */
	GMT_geoz_to_xy (x0+(double)data1*dev_x, y0+(double)data1*dev_y, z1, &xp1, &yp1);
	slope = (yp1-yp0)/(xp1-xp0);

	px0 = (int) ((xp0-z_project.xmin)*gmtdefs.dpi);
	px1 = (int) ((xp1-z_project.xmin)*gmtdefs.dpi);
	py0 = (int) ((yp0-z_project.ymin)*gmtdefs.dpi);
	py1 = (int) ((yp1-z_project.ymin)*gmtdefs.dpi);

/* now have the pixel locations for the two samples - join with a line..... */
	if (fabs(slope) <= 1.0){ /* more pixels needed in x direction */
		if (px0<px1){
			for (ix=px0; ix<=px1; ix++){
				iy = py0 + (int) (slope * (float) (ix - px0));
				paint(ix, iy);
			}
		}
		else{
			for (ix=px1; ix<=px0; ix++){
				iy = py0 + (int) (slope * (float) (ix - px0));
				paint(ix, iy);
			}

		}
	}
	else{ /* more pixels needed in y direction */
		if (py0<py1){
			for (iy=py0; iy<=py1; iy++){
				ix = px0 + (int) ( ((float) (iy-py0)) /slope);
				paint(ix, iy);
			}
		}
		else{
			for (iy=py1; iy<=py0; iy++){
				ix = px0 + (int) ( ((float) (iy-py0)) /slope);
				paint(ix, iy);
			}
		}
	}
}


void shade_bmap( double x0, double y0, float data0, float data1, double z0, double z1, int negative, double dev_x, double dev_y) /* apply current samples with all options to bitmap */ 
{
#define NPTS 4 /* 4 points for the general case here */
double xp[NPTS], yp[NPTS], temp, interp, slope01, slope02, slope12, slope13, slope23, slope03;
double slope0, slope1, slope2, slope3;
int ix, iy;

	if (data0 == 0.0 && data1 == 0.0)return; /* probably shouldn't strictly, but pathological enough I dont really want to deal with it! */

interp = 0.0;
	if ((data0*data1)<0.0){ 
/* points to plot are on different sides of zero - interpolate to find out where zero is */
		interp=z0+data0*((z0-z1)/(data1-data0));
		if(((data0<0.0) && negative) || ((data0>0.0)&& !negative)) { 
			/* plot from top to zero */
			z1=interp;
			data1=0.0;
		}
		else {
			z0=interp;
			data0=0.0;
		}
	}


	GMT_geoz_to_xy (x0+(double)data0*dev_x, y0+(double)data0*dev_y, z0, &xp[0], &yp[0]); /* returns 2 ends of line segment in plot coords */
	GMT_geoz_to_xy (x0+(double)data1*dev_x, y0+(double)data1*dev_y, z1, &xp[1], &yp[1]);
	GMT_geoz_to_xy (x0, y0, z0, &xp[2], &yp[2]); /* to get position of zero at each point*/
	GMT_geoz_to_xy (x0, y0, z1, &xp[3], &yp[3]); /* to get position of zero at each point*/

	/* now have four corner coordinates - need to sort them */
	for (ix=0; ix<NPTS-1; ix++)
		for (iy=ix+1; iy<NPTS; iy++)
			if(yp[ix]>yp[iy]){
				temp = yp[iy];
				yp[iy]=yp[ix];
				yp[ix]=temp;
				temp = xp[iy];
				xp[iy]=xp[ix];
				xp[ix]=temp;
			}


	/* have to fill the quadrilateral defined by 4 points (now ordered, but care with degenerate cases)*/

	slope01 = (xp[1]-xp[0])/(yp[1]-yp[0]);
	slope02 = (xp[2]-xp[0])/(yp[2]-yp[0]);
	slope12 = (xp[2]-xp[1])/(yp[2]-yp[1]);
	slope13 = (xp[3]-xp[1])/(yp[3]-yp[1]);
	slope23 = (xp[3]-xp[2])/(yp[3]-yp[2]);
	slope03 = (xp[3]-xp[0])/(yp[3]-yp[0]);
	if ((yp[0]!=yp[1]) && (yp[2]!=yp[3])){ /* simple case: tri-quad-tri */
		shade_tri(xp[0], yp[0], yp[1], slope01, slope02);
		shade_quad(xp[1], yp[1],xp[0]+slope02*(yp[1]-yp[0]), yp[2], slope02, slope13);
		shade_tri(xp[3], yp[3], yp[2], slope13, slope23);
	}
	if ((yp[0]==yp[1]) && (yp[2]!=yp[3])){
		if (xp[0]==xp[1]){ /* two triangles based on yp[1],yp[2]. yp[3] */
			shade_tri(xp[1], yp[1], yp[2], slope12, slope13);
			shade_tri(xp[3], yp[3], yp[2], slope23, slope13);
		}else{ /* quad based on first 3 points, then tri */
			slope0 = (((xp[0]<xp[1]) && (xp[3]<xp[2])) || ((xp[0]>xp[1])&&(xp[3]>xp[2])))*slope03 + (((xp[0]<xp[1])&&(xp[2]<xp[3])) || ((xp[0]>xp[1])&&(xp[2]>xp[3])))*slope02;
			slope1 = (((xp[1]<xp[0]) && (xp[3]<xp[2])) || ((xp[1]>xp[0]) && (xp[3]>xp[2])))*slope13 + (((xp[1]<xp[0])&&(xp[2]<xp[3])) || ((xp[1]>xp[0])&&(xp[2]>xp[3])))*slope12;
			slope3 = (((xp[1]<xp[0]) && (xp[3]<xp[2])) || ((xp[1]>xp[0]) && (xp[3]>xp[2])))*slope13 + (((xp[0]<xp[1])&&(xp[3]<xp[2])) || ((xp[0]>xp[1])&&(xp[3]>xp[2])))*slope03;
			shade_quad(xp[0], yp[0], xp[1], yp[2], slope0, slope1);
			shade_tri(xp[3], yp[3], yp[2], slope23, slope3);
		}
	}
	if ((yp[0]!=yp[1]) && (yp[2]==yp[3])){ 
		if(xp[2]==xp[3]){/* two triangles based on yp[0],yp[1]. yp[2] */
		shade_tri(xp[0], yp[0], yp[1], slope01, slope02);
		shade_tri(xp[2], yp[2], yp[1], slope12, slope02);
		}else{ /* triangle based on yp[0], yp[1], then quad based on last 3 points */
			slope0 = (((xp[0]<xp[1]) && (xp[3]<xp[2])) || ((xp[0]>xp[1]) && (xp[3]>xp[2])))*slope03 + (((xp[0]<xp[1])&&(xp[2]<xp[3])) || ((xp[0]>xp[1])&&(xp[2]>xp[3])))*slope02;
			shade_tri(xp[0], yp[0], yp[1], slope01, slope0);
			slope2 = (((xp[0]<xp[1]) && (xp[2]<xp[3])) || ((xp[0]>xp[1]) && (xp[2]>xp[3])))*slope02 + (((xp[0]<xp[1]) && (xp[3]<xp[2])) || ((xp[0]>xp[1]) && (xp[3]>xp[2])))*slope12;
			slope3 = (((xp[0]<xp[1]) && (xp[3]<xp[2])) || ((xp[0]>xp[1]) && (xp[3]>xp[2])))*slope03 + (((xp[0]<xp[1]) && (xp[2]<xp[3])) || 
((xp[0]>xp[1]) && (xp[2]>xp[3])))*slope13;
			shade_quad(xp[2], yp[2], xp[3], yp[1], slope2, slope3);
		}
	}

}

void shade_quad (double x0, double y0, double x1, double y_edge, double slope1, double slope0)
/* shade a quadrilateral with two sides parallel to x axis, one side at y=y0 with ends at x0 and x1,
with lines with gradients slope0 and slope1 respectively */
{
	int pedge_y, py0, iy, ix1, ix2, ix;

	if (y0 == y_edge) return;

	pedge_y = irint((y_edge-z_project.ymin)*gmtdefs.dpi);
	py0 = irint((y0-z_project.ymin)*gmtdefs.dpi);
	if (y0<y_edge){
		for (iy = py0; iy<pedge_y; iy++){
			ix1 = irint ((x0-z_project.xmin+(((double)iy /gmtdefs.dpi)+z_project.ymin - y0)*slope0)*gmtdefs.dpi);
			ix2 = irint ((x1-z_project.xmin+(((double)iy /gmtdefs.dpi)+z_project.ymin - y0)*slope1)*gmtdefs.dpi);
			if (ix1<ix2){
				for (ix=ix1; ix<ix2; ix++)
					paint(ix,iy);
			}else{
				for (ix=ix2; ix<ix1; ix++)
					paint(ix,iy);
			}
		}
	}else{
		for (iy = pedge_y; iy<py0; iy++){
			ix1 = irint ((x0-z_project.xmin+(((double)iy /gmtdefs.dpi)+z_project.ymin - y0)*slope0)*gmtdefs.dpi);
			ix2 = irint ((x1-z_project.xmin+(((double)iy /gmtdefs.dpi)+z_project.ymin - y0)*slope1)*gmtdefs.dpi);
			if (ix1<ix2){
				for (ix=ix1; ix<ix2; ix++)
					paint(ix,iy);
			}else{
				for (ix=ix2; ix<ix1; ix++)
					paint(ix,iy);
			}
		}
	}

}

void shade_tri (double apex_x, double apex_y, double edge_y, double slope, double slope0)
/* shade a triangle specified by apex coordinates, y coordinate of an edge (parallel to x-axis)
and slopes of the two other sides */
{
	int papex_y, pedge_y, iy, ix, x1, x2;

/* fprintf(stderr,"in shade_tri apex_x %f apex_y %f edge_y %f slope %f slope0 %f\n",apex_x, apex_y, edge_y, slope, slope0); */

	if (apex_y == edge_y) return;

	papex_y = irint((apex_y-z_project.ymin)*gmtdefs.dpi); /* location in pixels in y of apex and edge */
	pedge_y = irint((edge_y-z_project.ymin)*gmtdefs.dpi);
	if (apex_y < edge_y){
		for (iy = papex_y; iy<pedge_y; iy++){
			x1 = irint ((apex_x-z_project.xmin+(((double)iy /gmtdefs.dpi)+z_project.ymin - apex_y)*slope)*gmtdefs.dpi);
			x2 = irint ((apex_x-z_project.xmin+(((double)iy /gmtdefs.dpi)+z_project.ymin - apex_y)*slope0)*gmtdefs.dpi);
/*	fprintf(stderr,"apex_y<edge_y iy %d x1 %d x2 %d\n",iy,x1,x2);*/
				/* locations in pixels of two x positions for this scan line */
			if (x1<x2){
				for (ix=x1; ix<x2; ix++)
					paint(ix,iy);
			}else{
				for (ix=x2; ix<x1; ix++)
					paint(ix,iy);
			}
		}
	}else{
		for (iy = pedge_y; iy<papex_y; iy++){
			x1 = irint ((apex_x-z_project.xmin+(((double)iy /gmtdefs.dpi)+z_project.ymin - apex_y)*slope)*gmtdefs.dpi);
			x2 = irint ((apex_x-z_project.xmin+(((double)iy /gmtdefs.dpi)+z_project.ymin - apex_y)*slope0)*gmtdefs.dpi);
/*	fprintf(stderr,"apex_y>edge_y iy %d x1 %d x2 %d\n",iy,x1,x2); */
			if (x1<x2){
				for (ix=x1; ix<x2; ix++)
					paint(ix,iy);
			}
			else{
				for (ix=x2; ix<x1; ix++)
					paint(ix,iy);
			}
		}
	}
}

int paint(int ix, int iy)
/* ix iy is pixel to paint */
{
int byte, quot, rem;

/*fprintf(stderr,"painting ix %d iy %d\n",ix,iy);*/
	quot = ix/8;
	rem = ix - quot*8;

	if ((quot >= bm_nx-1) || (iy >= bm_ny-1) || (ix < 0) || (iy < 0))
			return (-1); /* outside bounds of plot array */

	byte = (bm_ny-iy-1)*bm_nx + quot; /* find byte to paint - flip vertical! */
	bitmap[byte] = bitmap[byte] | bmask[rem];
	return(0);
}



