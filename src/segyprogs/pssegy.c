/*--------------------------------------------------------------------
 *	$Id: pssegy.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *    Copyright (c) 1999-2014 by T. Henstock
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/* pssegy program to plot segy files in postscript with variable trace spacing option
 * uses the GMT pslib imagemask routines to plot a
 * 1-bit depth bitmap which will not obliterate stuff underneath!
 *
 * Author:	Tim Henstock (then@noc.soton.ac.uk)
 * Date:	1-July-1996
 * Version:	1.0
 *
 * Bug fixes:	1.1, 11-6-96 remove undesirable normalization introduced by -Z option. Add -U option for reduction.
 *
 * enhancements: 1.2, 7-20-98 add option to match location of traces from file
 *
 *               1.3, 1/7/99 check number of samples trace by trace to cope with SEGY with variable trace length
 *
 *               2.0, 5/7/99 update for GMT 3.3.1
 *
 *               2.1 10/4/2001 fix unreduce bug, modify to byte-swap if necessary, make 64-bit safe
 *
 *                   10/7/2009 fix bug when reading trace locations from arbitrary header locations,
 *                             8 bytes copied, should be 4 bytes 
 * 
 * This program is free software and may be copied or redistributed under the terms
 * of the GNU public license.
 */
 
#include "gmt.h"
#include "pslib.h"
#include "segy_io.h"

/* internal function calls */
	double rms( float *data, int nsamp);
	void wig_bmap(double x0, float data0, float data1, double y0, double y1);
	void shade_bmap(double x0, float data0, float data1, double y0, double y1, int negative);
	int paint( int ix, int iy);
	void plot_trace(float *data, double dy, double x0, int n_samp, int do_fill, int negative, int plot_wig, float toffset); 


unsigned char bmask[8]={128, 64, 32, 16, 8, 4, 2, 1};
unsigned char *bitmap;
int bm_nx, bm_ny;


int main (int argc, char **argv)
{
	GMT_LONG error = FALSE;
	int plot_cdp = FALSE, plot_offset = FALSE, byte_x = 0, doclip = FALSE;
	int normalize = FALSE, do_fill = FALSE, negative = FALSE, plot_wig = FALSE;
	int no_zero = FALSE, trace_file = FALSE;
#ifdef WORDS_BIGENDIAN
	int swap_bytes = FALSE;
#else
	int swap_bytes = TRUE;
#endif
	int i, nm;
	int ix, iy, n_traces=10000, n_samp=0, n_sampr=0, shade[3]={0,0,0}, trans[3]={-1,-1,-1}, n_tracelist=0;
	int check, polarity=1;
	int plot_it = FALSE;

	float bias = 0.0, scale = 1.0, clip = 0.0, deviation = 0.0;
	float toffset=0.0;
	double w, e, s, n, dx = 1.0, dy = 0.0; /* dx, dy are trace and sample interval */
	double xlen, ylen, xpix, ypix;
	double x0, test, err_head=0.0;
	double redvel=0.0;
	double *tracelist = VNULL;

	char input[512]="";
	char reelhead[3200];
	float *data = NULL;
	SEGYHEAD *header = NULL;
	char *head = NULL;
	int head2;
	SEGYREEL binhead;

	FILE *fpi = NULL, *fpt = NULL;


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
				/* parameters for wiggle mode */
				case 'F':
					do_fill = TRUE;
					if (GMT_getrgb (&argv[i][2], shade)) {
					  error++;
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
					switch(argv[i][2]){
						case 'o':
							plot_offset = TRUE;
							break;
						case 'c':
							plot_cdp = TRUE;
							break;
						case 'b':
							byte_x = atoi (&argv[i][3]);
							break;
					}
					break;
				/* trace scaling */
				case 'D':
					deviation = (float) atof (&argv[i][2]);
					break;
				/* over-rides of header info */
				case 'X': /* -X and -Y can be changed in gmt routines to lower case...*/
				case 'x':
					dx = atof (&argv[i][2]);
					break;
				case 'Y':
				case 'y':
					dy = atof (&argv[i][2]);
					fprintf(stderr," %s: Overriding sample interval dy = %f\n", GMT_program, dy);
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
				/* plot traces only at listed locations */
				case 'T':
					strcpy (input, &argv[i][2]);
					if ((trace_file = (int)strlen(input))){
						if ((fpt = fopen (input, "r")) == NULL){
							fprintf(stderr, "%s: Cannot find trace list file %s\n", GMT_program, input);
							error++;
						}
					}
					else{
						error++;
						fprintf(stderr, "%s: must give trace list filename\n", GMT_program);
					}
					break;
				case 'E':
					err_head = atof (&argv[i][2]);
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
			fprintf (stderr, "%s: Cannot find segy file %s\n", GMT_program,argv[i]);
			exit (EXIT_FAILURE);
		}
	}



	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "pssegy %s - Plot a segy file in PostScript\n\n", GMT_VERSION);
		fprintf (stderr, "usage: pssegy [<segyfile>] %s %s -D<dev> \n", GMT_Jx_OPT, GMT_Rx_OPT);
		fprintf (stderr, "	[-C<clip>] [-B<bias>] [-N] [-Z]\n");
		fprintf (stderr, "	[-F<gray>|<r/g/b>] [-I] [-W] [-S<header>]\n");
		fprintf (stderr, "	[-X<dx>] [-Y<dy>] [-L<nsamp>] [-M<ntraces>] \n");
		fprintf (stderr, "	[-U<redvel>] [-T<tracefile>] [-E<slop>] [-A] [-O] [-K] [-P]\n");

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, 
			"\n\t-Jx for projection.  Scale in INCH/units.  Specify one:\n\t -Jx<x-scale>              Linear projection\n\t-Jx<x-scale>l             Log10 projection\n\t  -Jx<x-scale>p<power>      x^power projection\n\tUse / to specify separate x/y scaling.\n\t If -JX is used then give axes lengths rather than scales\n\t regular map projections are not allowed!\n");
		GMT_explain_option ('R');
		fprintf (stderr, "	NB units for y are s or km\n");
		fprintf (stderr, "	-D<dev> to give deviation in X units of plot for 1.0 on scaled trace\n");
		fprintf (stderr, "	 IEEE SEGY file [or standard input] \n\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_explain_option ('V');
		fprintf (stderr, "	-C<clip> to clip scaled trace excursions at <clip>, applied after bias\n");
		fprintf (stderr, "	-B<bias> to bias scaled traces (-B-0.1 subtracts 0.1 from values)\n");
		fprintf (stderr,"	-N to trace normalize the plot\n");
		fprintf (stderr,"		order of operations: [normalize][bias][clip](deviation)\n");
		fprintf (stderr,"	-Z to suppress plotting traces whose rms amplitude is 0 \n");
		fprintf (stderr,"	-F<gray> to fill variable area with shade <gray>\n");
		fprintf (stderr,"	-Fr/g/b to fill variable area with color\n");
		fprintf (stderr,"		only a single color for the bitmap though!\n");
		fprintf (stderr,"	-I to fill negative rather than positive excursions\n");
		fprintf (stderr,"	-W to plot wiggle trace\n");
		fprintf (stderr,"	must specify either -W or -F\n");
		fprintf (stderr,"	-X<mult> multiply trace locations by <mult>\n");
		fprintf (stderr,"	-Y<dy> to override sample interval\n");
		fprintf (stderr,"	-S<header> to set variable spacing\n");
		fprintf (stderr,"		<header> is c for cdp or o for offset\n");
		fprintf (stderr,"	-L<nsamp> to override number of samples\n");
		fprintf (stderr,"	-M<ntraces> to fix number of traces. Default reads all traces.\n\t\t-M0 will read number in binary header, -Mn will attempt to read only n traces.\n");
		fprintf (stderr,"	-U<redvel> to apply reduction velocity (-ve removes reduction already present)\n");
		fprintf (stderr,"	-T<filename> to look in filename for a list of locations to select traces\n");
		fprintf (stderr,"		(same units as header * X, ie values printed by previous -V run)\n");
		fprintf (stderr,"	-E<error> slop to allow for -T. recommended in case of arithmetic errors!\n");
		fprintf (stderr,"	-A flips the default byte-swap state (default assumes data have a bigendian byte-order)\n");
		GMT_explain_option ('O');
		GMT_explain_option ('K');
		GMT_explain_option ('P');
		exit (EXIT_FAILURE);
	}

	if (err_head < 0.0){
		error++;
		fprintf(stderr,"%s: SYNTAX ERROR: slop cannot be negative\n", GMT_program);
	}
	if (negative && !do_fill){ /* negative with no fill */
		error++;
		fprintf(stderr,"%s: SYNTAX ERROR: Must specify -F with -I\n", GMT_program); 
	}
	if (!do_fill && !plot_wig){ /* no plotting specified */
		error++;
		fprintf(stderr,"%s: SYNTAX ERROR: Must specify -F or -W\n", GMT_program);
	}
	if (deviation <= 0.0){
		error++;
		fprintf(stderr,"%s: SYNTAX ERROR: Must specify a positive deviation\n",GMT_program);
	}
	if (!GMT_IS_LINEAR){
		fprintf(stderr,"%s: WARNING: you asked for a non-rectangular projection. \n It will probably still work, but be prepared for problems\n",GMT_program);
	}

	if (plot_cdp && plot_offset){
		fprintf(stderr,"%s: SYNTAX ERROR: Can't specify more than one trace location key\n",GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (fpi == NULL) fpi = stdin;

	if (trace_file){ /* must read in file of desired trace locations */
		tracelist = (double *) GMT_memory (CNULL, (size_t)GMT_CHUNK, sizeof(double), "pssegy");
		n_tracelist = GMT_CHUNK;
		ix=0;
		while ((fscanf(fpt, "%lf", &test)) != EOF){
			tracelist[ix] = test;
			ix++;
			if(ix == n_tracelist){ /* need more memory in array */
				n_tracelist += GMT_CHUNK;
				tracelist = (double *) GMT_memory ((char *)tracelist, (size_t)n_tracelist, sizeof(double), "pssegy");
			}
		}
		n_tracelist = ix;
		if(gmtdefs.verbose) fprintf(stderr, "%s: read in %d trace locations\n", GMT_program, n_tracelist);
	}

/* set up map projection and PS plotting */
	GMT_err_fail (GMT_map_setup (w, e, s, n), "");
	GMT_plotinit (argc, argv);

/* define area for plotting and size of array for bitmap */
	xlen = project_info.xmax-project_info.xmin;
	xpix = xlen*gmtdefs.dpi; /* pixels in x direction */
	/*xpix /= 8.0;
	bm_nx = 1 +(int) xpix;*/
	bm_nx = (int) ceil (xpix/8.0); /* store 8 pixels per byte in x direction but must have
				whole number of bytes per scan */
	ylen = project_info.ymax-project_info.ymin;
	ypix = ylen*gmtdefs.dpi; /* pixels in y direction */
	bm_ny = (int) ypix;
	nm = bm_nx*bm_ny;


/* read in reel headers from segy file */
	if ((check = get_segy_reelhd (fpi, reelhead)) != TRUE) exit(1);
	if ((check = get_segy_binhd (fpi, &binhead)) != TRUE) exit(1);

	if(swap_bytes){
/* this is a little-endian system, and we need to byte-swap ints in the reel header - we only
use a few of these*/
		if (gmtdefs.verbose) fprintf(stderr, "%s: swapping bytes for ints in the headers\n",GMT_program);
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
		fprintf(stderr,"%s: warning nsampr input %d, nsampr in header %d\n", GMT_program, n_sampr,  binhead.nsamp);

	if (!n_sampr){ /* no number of samples still - a problem! */
		fprintf(stderr, "%s: Error, number of samples per trace unknown\n", GMT_program);
		exit(EXIT_FAILURE);
	}

	if(gmtdefs.verbose) 
		fprintf(stderr, "%s: Number of samples for reel is %d\n", GMT_program, n_sampr);

	if(binhead.dsfc != 5) fprintf(stderr, "pssegy: WARNING data not in IEEE format\n");

	if (!dy){
		dy = (double) binhead.sr; /* sample interval of data (microseconds) */
		dy /= 1000000.0;
		fprintf(stderr,"%s: Sample interval is %f s\n", GMT_program, dy);
	}
	else if ((dy != binhead.sr) && (binhead.sr)) /* value in header overridden by input */
		fprintf(stderr, "%s: Warning dy input %f, dy in header %f\n", GMT_program, dy, (float)binhead.sr);

	if (!dy){ /* still no sample interval at this point is a problem! */
		fprintf(stderr, "%s: Error, no sample interval in reel header\n", GMT_program);
		exit(EXIT_FAILURE);
	}


	bitmap = (unsigned char *) GMT_memory (NULL, (size_t)nm, sizeof (unsigned char), "pssegy");

	ix=0;
	while ((ix<n_traces) && (header = get_segy_header(fpi))){   /* read traces one by one */

		if (plot_offset){ /* plot traces by offset, cdp, or input order */
			int32_t offset = ((swap_bytes)? GMT_swab4(header->sourceToRecDist): header->sourceToRecDist);
			x0 = (double) offset;
		}
		else if (plot_cdp){
			int32_t cdpval = ((swap_bytes)? GMT_swab4(header->cdpEns): header->cdpEns);
			x0 = (double) cdpval;
		}
		else if (byte_x){ /* ugly code - want to get value starting at byte_x of header into a double... */
			head = (char *) header;
			memcpy(&head2, &head[byte_x], 4); /* edited to fix bug where 8bytes were copied from head.
												Caused by casting to a long directly from char array*/ 
			x0 = (double) ((swap_bytes)? GMT_swab4(head2): head2);
		}
		else
			x0 = (1.0 + (double) ix);

		x0 *= dx;

		if (swap_bytes){
/* need to permanently byte-swap some things in the trace header 
do this after getting the location of where traces are plotted in case the general byte_x case
overlaps a defined header in a strange way */
			 header->sourceToRecDist=GMT_swab4(header->sourceToRecDist);
			 header->sampleLength=GMT_swab2(header->sampleLength);
			 header->num_samps=GMT_swab4(header->num_samps);
		}



/* now check that on list to plot if list exists */
		if (n_tracelist){
			plot_it = FALSE;
			for (i=0; i< n_tracelist; i++){
				if(fabs(x0-tracelist[i])<=err_head){
					plot_it = TRUE;
				}
			}
		}

		if (redvel){
			toffset = (float) -(fabs((double)(header->sourceToRecDist))/redvel);
			if (gmtdefs.verbose)
			fprintf(stderr, "pssegy: time shifted by %f\n",toffset);
		}

		data = (float *) get_segy_data(fpi, header); /* read a trace */
		/* get number of samples in _this_ trace (e.g. OMEGA has strange ideas about SEGY standard)
		or set to number in reel header */
		if ( !(n_samp = samp_rd(header)) ) n_samp = n_sampr;

		if(swap_bytes){ /* need to swap the order of the bytes in the data even though assuming IEEE format */
			int *intdata = (int *) data;
			for (iy=0; iy<n_samp; iy++){
				intdata[iy]=GMT_swab4(intdata[iy]);
			}
		}

		if(normalize || no_zero){
			scale=(float) rms(data, n_samp);
			if (gmtdefs.verbose) 
				fprintf(stderr, "pssegy: \t\t rms value is %f\n",scale);
		}
		for (iy=0; iy<n_samp; iy++){ /* scale bias and clip each sample in the trace */
			if (normalize) data[iy] /= scale;
			data[iy] += bias; 
			if(doclip && (fabs(data[iy]) > clip)) data[iy] = (float)(clip*data[iy]/fabs(data[iy])); /* apply bias and then clip */
			data[iy] *= deviation;
		}

		if ((!no_zero || scale) && (plot_it || !n_tracelist)){
			if (gmtdefs.verbose) 
				fprintf(stderr, "pssegy: trace %d plotting at %f \n", ix+1, x0);
			plot_trace (data, dy, x0, n_samp, do_fill, negative, plot_wig, toffset);
			}
		free (data);
		free (header);
		ix++;
	}

	GMT_map_clip_on (GMT_no_rgb, 3); /* set a clip at the map boundary since the image space overlaps a little */
	ps_bitimage (0.0,0.0,xlen, ylen, bitmap, 8*bm_nx, bm_ny, polarity, shade, trans);
/* have to multiply by 8 since pslib version of ps_imagemask is based on a _pixel_ count, whereas pssegy uses _byte_ count internally */
	GMT_map_clip_off ();

	if (fpi != stdin) fclose (fpi);

	GMT_plotend ();

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

double rms (float *data, int n_samp)
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

void plot_trace(float *data, double dy, double x0, int n_samp, int do_fill, int negative, int plot_wig, float toffset) 
	/* shell function to loop over all samples in the current trace, determine plot options
	 * and call the appropriate bitmap routine */ 
{
int iy;
int paint_wiggle;
double y0 = 0.0, y1;

	for(iy=1; iy<n_samp; iy++){ 	/* loop over samples on trace - refer to pairs iy-1, iy */
		y1 = dy * (float) iy + toffset;
		if (plot_wig) /* plotting wiggle */
			wig_bmap (x0, data[iy-1],data[iy], y0, y1);
		if (do_fill){ /* plotting VA -- check data points first */
			paint_wiggle = ( (!negative && ((data[iy-1]>=0.)||(data[iy]>=0.))) || (negative && ((data[iy-1]<=0.0)||(data[iy]<=0.0))) );
			if (paint_wiggle)
					shade_bmap (x0, data[iy-1], data[iy], y0, y1, negative);
		}
		y0=y1;
	}
}

void wig_bmap(double x0, float data0, float data1, double y0, double y1) /* apply current sample with all options to bitmap */
{
double xp0, xp1, yp0, yp1, slope;
int px0, px1, py0, py1, ix, iy;

	GMT_geo_to_xy (x0+ (double)data0, y0, &xp0, &yp0); /* returns 2 ends of line segment in plot coords */
	GMT_geo_to_xy (x0+ (double)data1, y1, &xp1, &yp1);
	slope = (yp1-yp0)/(xp1-xp0);

	px0 = (int) (xp0*gmtdefs.dpi);
	px1 = (int) (xp1*gmtdefs.dpi);
	py0 = (int) (yp0*gmtdefs.dpi);
	py1 = (int) (yp1*gmtdefs.dpi);

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


void shade_bmap(double x0, float data0, float data1, double y0, double y1, int negative) /* apply current samples with all options to bitmap */ 
{
double xp0, xp00, xp1, yp0, yp1, interp, slope;
int px0, px00, py0, py1, ixx, ix, iy;

	if ((data0*data1)<0.0){ 
/* points to plot are on different sides of zero - interpolate to find out where zero is */
		interp=y0+(double)data0*((y0-y1)/(double)(data1-data0));
		if(((data0<0.0) && negative) || ((data0>0.0)&& !negative)) { 
			/* plot from top to zero */
			y1=interp;
			data1=0.0;
		}
		else {
			y0=interp;
			data0=0.0;
		}
	}


	GMT_geo_to_xy (x0+(double)data0, y0, &xp0, &yp0); /* returns 2 ends of line segment in plot coords */
	GMT_geo_to_xy (x0+(double)data1, y1, &xp1, &yp1);
	GMT_geo_to_xy (x0, y0, &xp00, &yp0); /* to get position of zero */

	slope = (yp1-yp0)/(xp1-xp0);

	px0 = (int) (0.49+xp0*gmtdefs.dpi);
	px00 = (int) (0.49+xp00*gmtdefs.dpi);
	py0 = (int) (0.49+yp0*gmtdefs.dpi);
	py1 = (int) (0.49+yp1*gmtdefs.dpi);


/*  can rasterize simply by looping over values of y */
	if (py0<py1){
		for (iy=py0; iy<=py1; iy++){
			ixx = px0 + (int) ( ((float) (iy-py0)) /slope);
			if (ixx<px00){
				for (ix=ixx; ix<=px00; ix++)
					paint(ix, iy);
			}
			else{
				for (ix=px00; ix<=ixx; ix++)
					paint(ix, iy);
			}
		}
	}
	else{
		for (iy=py1; iy<=py0; iy++){
			ixx = px0 + (int) ( ((float) (iy-py0)) /slope);
			if (ixx<px00){
				for (ix=ixx; ix<=px00; ix++)
					paint(ix, iy);
			}
			else{
				for (ix=px00; ix<=ixx; ix++)
					paint(ix, iy);
			}
		}
	}


}

int paint(int ix, int iy)	/* pixel to paint */
{
int byte, quot, rem;


	quot = ix/8;
	rem = ix - quot*8;

	if ((quot >= bm_nx-1) || (iy >= bm_ny-1) || (ix < 0) || (iy < 0))
			return (-1); /* outside bounds of plot array */

	byte = (bm_ny-iy-1)*bm_nx + quot; /* find byte to paint - flip vertical! */
	bitmap[byte] = bitmap[byte] | bmask[rem];
	return (0);
}


