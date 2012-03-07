/*--------------------------------------------------------------------
 *	$Id: ps2raster.c,v 1.88 2011/07/08 21:27:06 guru Exp $
 *
 *	Copyright (c) 2006-2011 by J. Luis
 *	See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 *
 * ps2raster converts one or several PostScript file(s) to other formats using GhostScript.
 * It works by modifying the page size in order that the image will have a size
 * which is specified by the BoundingBox.
 * As an option, a tight BoundingBox may be computed.
 * ps2raster uses the ideas of the EPS2XXX.m from Primoz Cermelj published in MatLab Central
 * and of psbbox.sh of Remko Scharroo.
 *
 *
 *--------------------------------------------------------------------*/
/*
 * Authors:	Joaquim Luis and Remko Scharroo
 * Created:	15-FEB-2005
 * Modified:	15-NOV-2005 (It used to fail too many times in detecting if is GMT_ps)
 *		01-DEC-2005 Changed the generated script name to not overlap with -N
 *		27-DEC-2005 Integrated into GMT's misc supplement
 *		07-SEP-2006 Added EPS format; no more GMT specific; high-quality PDF output
 *		 1-APR-2007 Moved to GMT src directory
 *		 6-JUN-2007 Added -P option to force Portrait mode
 *		14-SEP-2007 Reduced bbox rendering from the default 4000 to 720 DPI. This
 *			    produces about 50% speedup for the -A option.
 *		 4-NOV-2007 Added -D to select alternative directory for output. Also added -V
 *			    and modifier -Au to remove GMT time-stamps.
 *		 6-NOV-2007 Can now process crappy Illu*er EPS files with its mix of \r and \n
 *		16-JUL-2008 Can now process %%HiResBoundingBox, if present
 *		14-NOV-2008 Added transparent PNG as another output format
 *		17-NOV-2008 New options to control text and graphics antialiasing
 *		28-DEC-2008 Added -F to select a different output file name
 *		17-JAN-2009 Added -W option to output world files
 *		05-MAR-2009 Added -Wg option to output a simple KML file
 *		28-MAR-2009 PW: Extended the -Wg option to select different altitude modes and to
 *					modify the document and overlay names
 */

#define GMT_WITH_NO_PS
#include "gmt.h"

#ifdef WIN32	/* Special for Windows */
#include <process.h>
#define getpid _getpid
#endif

#define N_GS_DEVICES		12	/* Number of supported GS output devices */
#define GS_DEV_EPS		0
#define GS_DEV_PDF		1
#define GS_DEV_JPG		2
#define GS_DEV_PNG		3
#define GS_DEV_PPM		4
#define GS_DEV_TIF		5
#define GS_DEV_BMP		6
#define GS_DEV_TPNG		7	/* PNG with transparency */
#define GS_DEV_JPGG		8	/* These are grayscale versions */
#define GS_DEV_PNGG		9
#define GS_DEV_TIFG		10
#define GS_DEV_BMPG		11

#define KML_GROUND_ABS		0
#define KML_GROUND_REL		1
#define KML_ABS				2
#define KML_SEAFLOOR_REL	3
#define KML_SEAFLOOR_ABS	4

struct PS2RASTER_CTRL {
	struct A {	/* -A[u][-] [Adjust boundingbox] */
		GMT_LONG active;
		GMT_LONG strip;	/* Remove the -U time-stamp */
		GMT_LONG reset;	/* The -A- turns -A off, overriding any automode in effect */
	} A;
	struct C {	/* -C<option> */
		GMT_LONG active;
	} C;
	struct D {	/* -D<dir> */
		GMT_LONG active;
		char *dir;
	} D;
	struct E {	/* -E<resolution> */
		GMT_LONG active;
		GMT_LONG dpi;
	} E;
	struct F {	/* -F<out_name> */
		GMT_LONG active;
		char *file;
	} F;
	struct G {	/* -G<GSpath> */
		GMT_LONG active;
		char *file;
	} G;
	struct L {	/* -L<listfile> */
		GMT_LONG active;
		char *file;
	} L;
	struct N {	/* -N */
		GMT_LONG active;
	} N;
	struct P2 {	/* -P */
		GMT_LONG active;
	} P;
	struct Q {	/* -Q[g|t]<bits> */
		GMT_LONG active;
		GMT_LONG on[2];	/* [0] for graphics, [1] for text antialiasing */
		GMT_LONG bits[2];
	} Q;
	struct S {	/* -S */
		GMT_LONG active;
	} S;
	struct T {	/* -T */
		GMT_LONG active;
		GMT_LONG eps;	/* TRUE if we want to make EPS (possibly in addition to another format) */
		GMT_LONG device;
	} T;
	struct V2 {	/* -V */
		GMT_LONG active;
	} V;
	struct W {	/* -W -- for world file production */
		GMT_LONG active;
		GMT_LONG warp;
		GMT_LONG kml;
		GMT_LONG mode;	/* 0 = clamp at ground, 1 is relative to ground, 2 is absolute 3 is relative to seafloor, 4 is clamp at seafloor */
		GMT_LONG min_lod, max_lod;	/* minLodPixels and maxLodPixels settings */
		GMT_LONG min_fade, max_fade;	/* minFadeExtent and maxFadeExtent settings */
		char *doctitle;		/* Name of KML document */
		char *overlayname;	/* Name of the image overlay */
		char *URL;		/* URL of remote site */
		double altitude;
	} W;
};

char * fgets2 (char *s, int n, FILE *stream)
{	/* Like fgets but breaks for either \r or \n */
	int c, k = 0, done = FALSE;
	n--;	/* Max chars to read */
	while (!done && (c = fgetc (stream)) != EOF) {
		if (c == '\r') c = '\n';	/* Replace ugly \r with nice \n */
		s[k++] = c;			/* Add to the growing buffer */
		if (c == '\n' || k == n) s[k] = '\0', done = TRUE;	/* Get \n or ran out of space */
	}
	return ((k == 0) ? NULL : s);
}

int main (int argc, char **argv) {
	GMT_LONG error = FALSE, grayscale = FALSE, found_proj = FALSE;

	GMT_LONG i, j, k, len, r, n_files = 0, pos_file, pos_ext, pix_w = 0, pix_h = 0, mode;
	GMT_LONG got_BB, got_HRBB, got_BBatend, file_has_HRBB, got_end, landscape, setup, i_unused = 0;
	double xt, yt, w, h, x0 = 0.0, x1 = 612.0, y0 = 0.0, y1 = 828.0;
	double	west = 0.0, east = 0.0, south = 0.0, north = 0.0;
	double	west_prj, south_prj, east_prj, north_prj;

	size_t n_alloc = GMT_SMALL_CHUNK;

	char ps_file[BUFSIZ];
	char **ps_names = NULL, *no_U_file = NULL, *clean_PS_file = NULL, *tmp_file = NULL, out_file[BUFSIZ], *BB_file = NULL;
	char line[BUFSIZ], c[BUFSIZ], *p = NULL, c1[20], c2[20], c3[20], c4[20], cmd[BUFSIZ], proj4_name[20], *quiet = NULL;
	char *gs_params = NULL, *gs_BB = NULL, *anti = NULL, text[BUFSIZ], gs_extra[BUFSIZ], *proj4_cmd = NULL;
	char *device[N_GS_DEVICES] = {"", "pdfwrite", "jpeg", "png16m", "ppmraw", "tiff24nc", "bmp16m", "pngalpha", "jpeggray", "pnggray", "tiffgray", "bmpgray"};
	char *ext[N_GS_DEVICES] = {".eps", ".pdf", ".jpg", ".png", ".ppm", ".tif", ".bmp", ".png", ".jpg", ".png", ".tif", ".bmp"};
	char *RefLevel[5] = {"clampToGround", "relativeToGround", "absolute", "relativeToSeaFloor", "clampToSeaFloor"};
#ifdef WIN32
	char at_sign[2] = "@";
#else
	char at_sign[2] = "";
#endif

#ifdef USE_GDAL
	struct GDALREAD_CTRL *to_gdalread = NULL;
	struct GD_CTRL *from_gdalread = NULL;
#endif

	FILE *fp = NULL, *fpo = NULL, *fpb = NULL, *fpl = NULL, *fp2 = NULL, *fpw = NULL;

	struct PS2RASTER_CTRL *Ctrl = NULL;

	void *New_ps2raster_Ctrl (), Free_ps2raster_Ctrl (struct PS2RASTER_CTRL *C);
	GMT_LONG parse_GE_settings (char *arg, struct PS2RASTER_CTRL *C);
	
	Ctrl = (struct PS2RASTER_CTRL *)New_ps2raster_Ctrl ();	/* Allocate and initialize a new control structure */

	/* Parameters for all the formats available */

	gs_params = "-q -dSAFER -dNOPAUSE -dBATCH -dUseFlateCompression=true -dPDFSETTINGS=/prepress -dEmbedAllFonts=true -dSubsetFonts=true -dMonoImageFilter=/FlateEncode -dAutoFilterGrayImages=false -dGrayImageFilter=/FlateEncode -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode";
	gs_BB = "-q -dSAFER -dNOPAUSE -dBATCH -sDEVICE=bbox"; /* -r defaults to 4000, see http://pages.cs.wisc.edu/~ghost/doc/cvs/Devices.htm#Test */
	memset ((void *)gs_extra, 0, BUFSIZ);
	
	/* Check and interpret the command line arguments */

	GMT_program = argv[0];
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {
				case '\0':
					GMT_give_synopsis_and_exit = TRUE;
					break;

				/* Supplemental parameters */
				case 'A':	/* Adjust BoundingBox */
					Ctrl->A.active = TRUE;
					if (argv[i][2] == 'u') Ctrl->A.strip = TRUE;
					if (argv[i][strlen(argv[i])-1] == '-') Ctrl->A.reset = TRUE;
					break;
				case 'C':	/* Append extra custom GS options */
					strcat (gs_extra, " ");		/* Append to list of extra GS options */
					strcat (gs_extra, &argv[i][2]);	/* Append to list of extra GS options */
					break;
				case 'D':	/* Change output directory */
					Ctrl->D.active = TRUE;
					Ctrl->D.dir = strdup (&argv[i][2]);
					break;
				case 'G':	/* Set GS path */
					Ctrl->G.active = TRUE;
					free ((void *)Ctrl->G.file);
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
				case 'E':	/* Set output dpi */
					Ctrl->E.active = TRUE;
					Ctrl->E.dpi = atoi(&argv[i][2]);
					break;
				case 'F':	/* Set explicitly the output file name */
					Ctrl->F.active = TRUE;
					Ctrl->F.file = strdup (&argv[i][2]);
					GMT_chop_ext (Ctrl->F.file);	/* Make sure file name has no extension */
					break;
				case 'L':	/* Give list of files to convert */
					Ctrl->L.active = TRUE;
					Ctrl->L.file = strdup (&argv[i][2]);
					break;
				case 'N':	/* Do NOT remove auxillary files used with GS */
					Ctrl->N.active = TRUE;
					break;
				case 'P':	/* Force Portrait mode */
					Ctrl->P.active = TRUE;
					break;
				case 'Q':	/* Anti-aliasing settings */
					Ctrl->Q.active = TRUE;
					if (argv[i][2] == 'g') {
						mode = 0;
						anti = "-dGraphicsAlphaBits=";
					}
					else if (argv[i][2] == 't') {
						mode = 1;
						anti = "-dTextAlphaBits=";
					}
					else {
						fprintf (stderr, "%s: GMT ERROR: Unrecognized option %s\n", GMT_program, argv[i]);
						error = TRUE;
						continue;
						
					}
					Ctrl->Q.on[mode] = TRUE;
					Ctrl->Q.bits[mode] = (argv[i][3]) ? atoi (&argv[i][3]) : 4;
					sprintf (text, " %s%ld", anti, Ctrl->Q.bits[mode]);
					strcat (gs_extra, text);	/* Append to list of extra GS options */
					break;
				case 'S':	/* Write the GS command to STDOUT */
					Ctrl->S.active = TRUE;
					break;
				case 'T':	/* Select output format (optionally also request EPS) */
					Ctrl->T.active = TRUE;
					grayscale = ((j = (int)strlen(argv[i])) > 3 && argv[i][j-1] == '-');
					for (j = 2; argv[i][j]; j++) {
						switch (argv[i][j]) {
							case 'e':	/* EPS */
								Ctrl->T.eps = TRUE;
								break;
							case 'f':	/* PDF */
								Ctrl->T.device = GS_DEV_PDF;
								break;
							case 'b':	/* BMP */
								Ctrl->T.device = (grayscale) ? GS_DEV_BMPG : GS_DEV_BMP;
								break;
							case 'j':	/* JPEG */
								Ctrl->T.device = (grayscale) ? GS_DEV_JPGG : GS_DEV_JPG;
								break;
							case 'g':	/* PNG */
								Ctrl->T.device = (grayscale) ? GS_DEV_PNGG : GS_DEV_PNG;
								break;
							case 'G':	/* PNG (transparent) */
								Ctrl->T.device = GS_DEV_TPNG;
								strcat (gs_extra, " -dMaxBitmap=100000000");	/* Append to list of extra GS options */
								break;
							case 'm':	/* PPM */
								Ctrl->T.device = GS_DEV_PPM;
								break;
							case 't':	/* TIFF */
								Ctrl->T.device = (grayscale) ? GS_DEV_TIFG : GS_DEV_TIF;
								break;
							case '-':	/* Just skip the trailing - for grayscale since it is handled separately */
								break;
							default:
								fprintf (stderr, "%s: GMT ERROR: Unrecognized option %s\n", GMT_program, argv[i]);
								error = TRUE;
								break;
						}
					}
					break;
				case 'V':	/* Verbose */
					Ctrl->V.active = TRUE;
					break;
				case 'W':	/* Save world file */
					error = parse_GE_settings (&argv[i][2], Ctrl);
					break;
				default:	/* Options not recognized */
					error = TRUE;
					break;
			}
		}
		else {
			strcpy (ps_file, argv[i]);
			n_files++;
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit || error) {
		fprintf (stderr,"ps2raster %s - Converts one or several [E]PS file(s) to raster formats using GhostScript.\n\n", GMT_VERSION);
		fprintf(stderr,"usage: ps2raster <psfile1> <psfile2> <...> [-A[u][-]] [-C<gs_command>] [-D<dir>]\n");
		fprintf(stderr,"       [-E<resolution>] [-F<out_name>] [-G<ghost_path>] [-L<listfile>] [-N] [-P]\n");
		fprintf(stderr,"       [-Q[g|t]1|2|4] [-S] [-Tb|e|f|g|G|j|m|t] [-V]\n");
		fprintf(stderr,"       [-W[+g][+k][+t<title>][+n<name>][+a<mode>[<alt]][+l<lodmin>/<lodmax>][+f<minfade>/<maxfade>][+u<URL>]]\n\n");
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr,"Works by modifying the page size in order that the resulting\n");
		fprintf (stderr,"image will have the size specified by the BoundingBox.\n");
		fprintf (stderr,"As an option, a tight BoundingBox may be computed.\n\n");
		fprintf (stderr,"	<psfile(s)> postscript file(s) to be converted.\n");
		fprintf (stderr,"\n\tOPTIONS:\n");
		fprintf (stderr,"\t-A Adjust the BoundingBox to the minimum required by the image contents\n");
		fprintf (stderr,"\t   Append u to strip out time-stamps (produced by GMT -U options)\n");
		fprintf (stderr,"\t   Append - to make sure -A is NOT activated by -W\n");
		fprintf (stderr,"\t-C Specify a single, custom option that will be passed on to GhostScript\n");
		fprintf (stderr,"\t   as is. Repeat to add several options [none].\n");
		fprintf (stderr,"\t-D Sets an alternative output directory (which must exist) [Default is same directory as PS files]\n");
		fprintf (stderr,"\t   Use -D. to place the output in the current directory.\n");
		fprintf (stderr,"\t-E Set raster resolution in dpi [default = 720 for PDF, 300 for others]\n");
		fprintf (stderr,"\t-F Force the output file name. By default output names are constructed\n");
		fprintf (stderr,"\t   using the input names as base, which are appended with an appropriate\n");
		fprintf (stderr,"\t   extension. Use this option to provide a different name, but WITHOUT\n");
		fprintf (stderr,"\t   extension. Extension is still determined automatically.\n");
		fprintf (stderr,"\t-G Full path to your ghostscript executable.\n");
		fprintf (stderr,"\t   NOTE: Under Unix systems this is generally not necessary.\n");
		fprintf (stderr,"\t   Under Windows, ghostscript is not added to the system's path.\n");
		fprintf (stderr,"\t   So either you do it yourself, or give the full path here.\n");
		fprintf (stderr,"\t   (e.g. -Gc:\\programs\\gs\\gs7.05\\bin\\gswin32c).\n");
		fprintf (stderr,"\t-L The <listfile> is an ASCII file with names of ps files to be converted\n");
		fprintf (stderr,"\t-N OBSOLETE. Use -S and/or -Te instead.\n");
		fprintf (stderr,"\t-P Force Portrait mode. All Landscape mode plots will be rotated back\n");
		fprintf (stderr,"\t   so that they show unrotated in Portrait mode.\n");
		fprintf (stderr,"\t   This is practical when converting to image formats or preparing\n");
		fprintf (stderr,"\t   EPS or PDF plots for inclusion in documents.\n");
		fprintf (stderr,"\t-Q Anti-aliasing setting for (g)raphics or (t)ext; append size (1,2,4) of sub-sampling box\n");
		fprintf (stderr,"\t   Default is no anti-aliasing, which is the same as specifying size 1.\n");
		fprintf (stderr,"\t-S Apart from executing it, also writes the ghostscript command to standard output.\n");
		fprintf (stderr,"\t-T Set output format [default is jpeg]\n");
		fprintf (stderr,"\t   b means BMP\n");
		fprintf (stderr,"\t   e means EPS\n");
		fprintf (stderr,"\t   f means PDF\n");
		fprintf (stderr,"\t   g means PNG\n");
		fprintf (stderr,"\t   G means PNG (with transparency)\n");
		fprintf (stderr,"\t   j means JPEG\n");
		fprintf (stderr,"\t   m means PPM\n");
		fprintf (stderr,"\t   t means TIF\n");
		fprintf (stderr,"\t   For b, g, j, t, append - to get a grayscale image [24-bit color].\n");
		fprintf (stderr,"\t   The EPS format can be combined with any of the other formats.\n");
		fprintf (stderr,"\t   For example, -Tef creates both an EPS and PDF file.\n");
		fprintf (stderr,"\t-V Provides progress report [default is silent] and shows the gdal_translate\n");
		fprintf (stderr,"\t   command, in case you want to use this program to create a geoTIFF file.\n");
		fprintf (stderr,"\t-W Write a ESRI type world file suitable to make (e.g.) .tif files be\n");
		fprintf (stderr,"\t   recognized as geotiff by softwares that know how to do it. Be aware,\n");
		fprintf (stderr,"\t   however, that different results are obtained depending on the image\n");
		fprintf (stderr,"\t   contents and if the -B option has been used or not. The trouble with\n");
		fprintf (stderr,"\t   -B is that it creates a frame and very likely its annotations and that\n");
		fprintf (stderr,"\t   introduces pixels outside the map data extent. As a consequence, the\n");
		fprintf (stderr,"\t   map extents estimation will be wrong. To avoid this problem, use the\n");
		fprintf (stderr,"\t   --BASEMAP_TYPE=inside option which plots all annotations related stuff\n");
		fprintf (stderr,"\t   inside the image and does not compromise the coordinate computations.\n");
		fprintf (stderr,"\t   The world file naming follows the convention of jamming a 'w' in the\n");
		fprintf (stderr,"\t   file extension. So, if the output is tif (-Tt) the world file is a .tfw,\n");
		fprintf (stderr,"\t   for jpeg a .jgw, and so on.\n");
		fprintf (stderr,"\t   Use -W+g to do a system call to gdal_translate and produce a true geoTIFF\n");
		fprintf (stderr,"\t   image right away. The output file will have the extension .tiff\n");
		fprintf (stderr,"\t   See the man page for other 'gotchas'. Automatically sets -A -P.\n");
		fprintf (stderr,"\t   Use -W+k to create a minimalist KML file that allows loading the image in\n");
		fprintf (stderr,"\t   Google Earth. Note that for this option the image must be in geographical\n");
		fprintf (stderr,"\t   coordinates. If not, a warning is issued but the KML file is created anyway.\n");
		fprintf (stderr,"\t   Several modifiers allow you to specify the content in the KML file:\n");
		fprintf (stderr,"\t   +t<doctitle> sets the document name [\"GMT KML Document\"]\n");
		fprintf (stderr,"\t   +n<layername> sets the name of this particular layer [\"GMT Image Overlay\"]\n");
		fprintf (stderr,"\t   +a<altmode>[<altitude>] sets the altitude mode of this layer, where\n");
		fprintf (stderr,"\t      <altmode> is one of 5 recognized by Google Earth:\n");
		fprintf (stderr,"\t      G clamped to the ground [Default]\n");
		fprintf (stderr,"\t      g Append altitude (in m) relative to ground\n");
		fprintf (stderr,"\t      A Append absolute altitude (in m)\n");
		fprintf (stderr,"\t      s Append altitude (in m) relative to seafloor\n");
		fprintf (stderr,"\t      S clamped to the seafloor\n");
		fprintf (stderr,"\t   +l<minLOD>/<maxLOD>] sets Level Of Detail when layer should be active [always active]\n");
		fprintf (stderr,"\t     Image goes inactive when there are fewer than minLOD pixels or more\n");
		fprintf (stderr,"\t     than maxLOD pixels visible.  -1 means never invisible.\n");
		fprintf (stderr,"\t   +f<minfade>/<maxfade>] sets distances over which we fade from opaque to transparent [no fading]\n");
		fprintf (stderr,"\t   +u<URL> prepands this URL to the name of the image referenced in the KML [local file]\n");

		exit (EXIT_FAILURE);
	}

	if (Ctrl->N.active) {
		fprintf (stderr, "%s: GMT WARNING: Option -N is obsolete. Will run with -S -Te.\n", GMT_program);
		Ctrl->S.active = Ctrl->T.eps = TRUE;
	}

	if (Ctrl->Q.on[0] && (Ctrl->Q.bits[0] < 1 || Ctrl->Q.bits[0] > 4)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR: Anti-aliasing for graphics requires sub-samplib box of 1,2, or 4\n", GMT_program);
		error = TRUE;
	}
	if (Ctrl->Q.on[1] && (Ctrl->Q.bits[1] < 1 || Ctrl->Q.bits[1] > 4)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR: Anti-aliasing for text requires sub-samplib box of 1,2, or 4\n", GMT_program);
		error = TRUE;
	}
	if (!Ctrl->T.active) Ctrl->T.device = GS_DEV_JPG;	/* Default output device if none is specified */

	if (error) exit (EXIT_FAILURE);

	if (n_files > 1 && Ctrl->L.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR: Cannot handle both a file list and multiple ps files in input\n", GMT_program);
		error = TRUE;
	}

	if (Ctrl->L.active && (fpl = fopen (Ctrl->L.file, "r")) == NULL) {
		fprintf (stderr, "%s: GMT ERROR: Cannot to open list file %s\n", GMT_program, Ctrl->L.file);
		error = TRUE;
	}

	if (Ctrl->F.active && (Ctrl->L.active || Ctrl->D.active)) {
		fprintf (stderr, "%s: GMT WARNING: Option -F and options -L OR -D are mutually exclusive. Ignoring option -F.\n", GMT_program);
		Ctrl->F.active = FALSE;
	}

	if (Ctrl->F.active && n_files > 1) {
		fprintf (stderr, "%s: GMT WARNING: Option -F is incompatible with multiple inputs. Ignoring option -F.\n", GMT_program);
		Ctrl->F.active = FALSE;
	}

	if (Ctrl->W.kml && !(Ctrl->T.device == GS_DEV_JPG || Ctrl->T.device == GS_DEV_JPGG || Ctrl->T.device == GS_DEV_TIF || 
				Ctrl->T.device == GS_DEV_TIFG || Ctrl->T.device == GS_DEV_PNG ||
				Ctrl->T.device == GS_DEV_TPNG || Ctrl->T.device == GS_DEV_PNGG) ) {
		fprintf (stderr, "%s: WARNERROR: As far as we know selected raster type is unsuported by GE\n", GMT_program);
	}

	if (Ctrl->W.active) {	/* Implies -P and -A (unless -A- is set ) */
		Ctrl->P.active = Ctrl->A.active = TRUE;
		if (Ctrl->A.reset) Ctrl->A.active = FALSE;
	}

	/* Use default DPI if not already set */
	if (Ctrl->E.dpi <= 0) Ctrl->E.dpi = (Ctrl->T.device == GS_DEV_PDF) ? 720 : 300;

	/* Multiple files in a file with their names */
	if (Ctrl->L.active) {
		ps_names = (char **) GMT_memory (VNULL, n_alloc, sizeof (char *), GMT_program);
		while (fgets2 (line, BUFSIZ, fpl) != NULL) {
			ps_names[n_files] = (char *) GMT_memory (VNULL, (size_t)BUFSIZ, sizeof (char), GMT_program);
			if (line[0] == '#' || line[0] == '\n') continue;
			sscanf (line, "%s", ps_names[n_files]);
			n_files++;
			if (n_files > (int)n_alloc) {
				n_alloc <<= 1;
				ps_names = GMT_memory ((void *)ps_names, n_alloc, sizeof (char *), GMT_program);
			}
		}
		fclose (fpl);
	}

	/* Multiple files given on command line */

	else if (n_files > 1) {
		ps_names = (char **) GMT_memory (VNULL, n_alloc, sizeof (char *), GMT_program);
		j = 0;
		for (k = 1; k < argc; k++) {
			if (argv[k][0] == '-') continue;
			ps_names[j] = (char *) GMT_memory (VNULL, (size_t)BUFSIZ, sizeof (char), GMT_program);
			ps_names[j] = argv[k];
			j++;
			if (n_files > (int)n_alloc) {
				n_alloc <<= 1;
				ps_names = GMT_memory ((void *)ps_names, n_alloc, sizeof (char *), GMT_program);
			}
		}
	}
	else {				/* Single file */
		ps_names = (char **) GMT_memory (VNULL, (size_t)1, sizeof (char *), GMT_program);
		ps_names[0] = ps_file;
	}

	/* Loop over all input files */

	for (k = 0; k < n_files; k++) {
		memset ((void *)out_file, 0, BUFSIZ);
		strcpy(ps_file,ps_names[k]);
		if ((fp = fopen (ps_file, "r")) == NULL) {
			fprintf (stderr, "%s: Cannot to open file %s\n", GMT_program, ps_file);
			continue;
		}

		if (Ctrl->V.active) fprintf (stderr, "%s: Processing %s:", GMT_program, ps_file);

		if (Ctrl->A.strip) {	/* Must strip off the GMT timestamp stuff */
			GMT_LONG dump = TRUE;
			if (Ctrl->V.active) fprintf (stderr, " Strip GMT time-stamp...");
			no_U_file = (char *) GMT_memory (VNULL, (size_t)BUFSIZ, sizeof(char), GMT_program);
			sprintf (no_U_file, "%s/ps2raster_%db.eps", Ctrl->D.dir, (int)getpid());
			if ((fp2 = fopen (no_U_file, "w+")) == NULL) {
				fprintf (stderr, "%s: Unable to create a temporary file\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			while (fgets2 (line, BUFSIZ, fp) != NULL) {
				if (!strncmp (line, "% Begin GMT time-stamp", (size_t)22)) dump = FALSE;
				if (dump) fprintf (fp2, "%s", line);
				if (!strncmp (line, "% End GMT time-stamp", (size_t)20)) dump = TRUE;
			}
			fclose (fp);	/* Close original PS file */
			rewind (fp2);	/* Rewind new file without timestamp */
			fp = fp2;	/* Set original file pointer to this file instead */
		}

		got_BB = got_HRBB = file_has_HRBB = got_BBatend = got_end = landscape = setup = FALSE;

		len = (GMT_LONG)strlen(ps_file);
		j = len - 1;
		pos_file = -1;
		pos_ext = -1;	/* In case file has no extension */
		for (i = 0; i < len; i++, j--) {
			if (pos_ext < 0 && ps_file[j] == '.') pos_ext = j;	/* Beginning of file extension */
			if (pos_file < 0 && (ps_file[j] == '/' || ps_file[j] == '\\')) pos_file = j + 1;	/* Beginning of file name */
		}
		if (pos_ext == -1) pos_ext = len - 1;	/* File has no extension */
		if (!Ctrl->D.active || pos_file == -1) pos_file = 0;	/* File either has no leading directory or we want to use it */

		/* Adjust to a tight BoundingBox if user requested so */

		if (Ctrl->A.active) {
			char *psfile_to_use;
			if (Ctrl->V.active) fprintf (stderr, " Find HiResBoundingBox ");
			BB_file = (char *) GMT_memory (VNULL, (size_t)BUFSIZ, sizeof(char), GMT_program);
			sprintf (BB_file, "%s/ps2raster_%ldc.bb", Ctrl->D.dir, (GMT_LONG)getpid());
			psfile_to_use = (Ctrl->A.strip) ? no_U_file : ((clean_PS_file) ? clean_PS_file : ps_file);
			sprintf (cmd, "%s%s %s %s %s 2> %s", at_sign, Ctrl->G.file, gs_BB, gs_extra, psfile_to_use, BB_file);
			i_unused = system (cmd);		/* Execute the command that computes the tight BB */
			if ((fpb = fopen (BB_file, "r")) == NULL) {
				fprintf (stderr, "%s: Unable to open file %s\n", GMT_program, BB_file);
				exit (EXIT_FAILURE);
			}
			while (fgets2 (line, BUFSIZ, fpb) != NULL && !got_BB) {	/* We only use the High resolution BB */
				if ((p = strstr (line,"%%HiResBoundingBox:"))) {
					sscanf (&p[19], "%s %s %s %s", c1, c2, c3, c4);
					x0 = atof (c1);		y0 = atof (c2);
					x1 = atof (c3);		y1 = atof (c4);
					if (x1 <= x0 || y1 <= y0) {
						fprintf (stderr, "%s: Unable to decode BoundingBox file %s\n", GMT_program, BB_file);
						fclose (fpb);
						remove (BB_file);	/* Remove the file */
						tmp_file = (char *) GMT_memory (VNULL, (size_t)BUFSIZ, sizeof(char), GMT_program);
						sprintf (tmp_file, "%s/", Ctrl->D.dir);
						strncat (tmp_file, &ps_file[pos_file], (size_t)(pos_ext - pos_file));
						strcat (tmp_file, ext[Ctrl->T.device]);
						sprintf (cmd, "%s%s %s %s -sDEVICE=%s -g1x1 -r%ld -sOutputFile=%s -f%s", 
							at_sign, Ctrl->G.file, gs_params, gs_extra, device[Ctrl->T.device],
							Ctrl->E.dpi, tmp_file, ps_file);
						i_unused = system (cmd);		/* Execute the GhostScript command */
						if (Ctrl->S.active) fprintf (stdout, "%s\n", cmd);

						GMT_free((void *) tmp_file);
						continue;
					}
					got_BB = got_HRBB = TRUE;
				}
			}
			fclose (fpb);
			remove (BB_file);	/* Remove the file with BB info */
			GMT_free ((void *)BB_file);
			if (got_BB && Ctrl->V.active) fprintf (stderr, "[%g %g %g %g]...", x0, y0, x1, y1);
		}

		/* Open temporary file to be processed by ghostscript. When -Te is used, tmp_file is for keeps */

		if (Ctrl->V.active && Ctrl->T.eps) fprintf (stderr, " Format EPS file...");
		tmp_file = (char *) GMT_memory (VNULL, (size_t)BUFSIZ, sizeof(char), GMT_program);
		if (Ctrl->T.eps) {
			if (Ctrl->D.active) sprintf (tmp_file, "%s/", Ctrl->D.dir);	/* Use specified output directory */
			if (!Ctrl->F.active)
				strncat (tmp_file, &ps_file[pos_file], (size_t)(pos_ext - pos_file));
			else
				strcat (tmp_file, Ctrl->F.file);
			strcat (tmp_file, ext[GS_DEV_EPS]);
			if ((fpo = fopen (tmp_file, "w")) == NULL) {
				fprintf (stderr, "%s: Unable to open file %s for writing\n", GMT_program, tmp_file);
				continue;
			}
		}
		else {
			sprintf (tmp_file, "%s/ps2raster_%ldd.eps", Ctrl->D.dir, (GMT_LONG)getpid());
			if ((fpo = fopen (tmp_file, "w+")) == NULL) {
				fprintf (stderr, "%s: Unable to create a temporary file\n", GMT_program);
				continue;
			}
		}

		/* Scan first 20 lines of input file for [HiRes]BoundingBox and Orientation statements.
		 * Since we prefer the HiResBB over BB we must continue to read until both are found or 20 lines have past */

		i = 0;
		while ((fgets2 (line, BUFSIZ, fp) != NULL) && i < 20 && !(got_BB && got_HRBB && got_end)) {
			i++;
			if (!line[0] || line[0] != '%') 
				{ /* Skip empty and non-comment lines */ }
			else if (!got_BB && (p = strstr (line, "%%BoundingBox:"))) {
				sscanf (&p[14], "%s %s %s %s",c1,c2,c3,c4);
				if (strncmp (c1, "(atend)", (size_t)7)) {	/* Got actual numbers */
					if (!got_HRBB) {	/* Only assign values if we havent seen the high-res version yet */
						x0 = atoi (c1);		y0 = atoi (c2);
						x1 = atoi (c3);		y1 = atoi (c4);
					}
					got_BB = TRUE;
				}
				else
					got_BBatend++;
			}
			else if ((p = strstr (line, "%%HiResBoundingBox:"))) {
				file_has_HRBB = TRUE;
				if (!got_HRBB) {
					sscanf (&p[19], "%s %s %s %s",c1,c2,c3,c4);
					if (strncmp (c1, "(atend)", (size_t)7)) {	/* Got actual numbers */
						x0 = atof (c1);		y0 = atof (c2);
						x1 = atof (c3);		y1 = atof (c4);
						got_HRBB = got_BB = TRUE;
					}
				}
			}
			else if ((p = strstr (line, "%%Orientation:"))) {
				sscanf (&p[14], "%s", c1);
				if (!strncmp (c1, "Landscape", (size_t)9)) landscape = TRUE;
			}
			else if ((p = strstr (line, "%%EndComments")))
				got_end = TRUE;
			if (got_BBatend == 1 && (got_end || i == 19)) {	/* Now is the time to look at the end of the file */
				got_BBatend++;			/* Avoid jumping more than once to the end */
				if (!fseek (fp, (long)-256, SEEK_END)) i = -30;
			}
		}

		/* Cannot proceed without knowing the BoundingBox */

		if (!got_BB) {
			fprintf (stderr, "%s: GMT FATAL ERROR: The file %s has no BoundingBox in the first 20 lines or last 256 bytes. Use -A option.\n", GMT_program, ps_file);
			continue;
		}

		/* Do the math on the BoundingBox and translation coordinates */

		if (Ctrl->P.active && landscape)
			xt = -x1, yt = -y0, w = y1-y0, h = x1-x0, r = -90;
		else
			xt = -x0, yt = -y0, w = x1-x0, h = y1-y0, r = 0;

		/* Rewind the input file and start copying and replacing */

		rewind (fp);
		while (fgets2 (line, BUFSIZ, fp) != NULL) {
			if (line[0] != '%') {	/* Copy any non-comment line, except one containing /PageSize in the Setup block */
				if (setup && strstr(line,"/PageSize") != NULL) continue;
				fprintf (fpo, "%s", line);
				continue;
			}
			else if (Ctrl->W.active && !found_proj) {
				if (!strncmp(&line[2], "PROJ", 4)) { /* Search for the PROJ tag in the ps file */
					char *ptmp = NULL, xx1[20], xx2[20], yy1[20], yy2[20];
					sscanf (&line[8], "%s %s %s %s %s %s %s %s %s",proj4_name,xx1,xx2,yy1,yy2,c1,c2,c3,c4);
					west_prj = atof (c1);		east_prj = atof (c2);
					south_prj = atof (c3);		north_prj = atof (c4);
					project_info.w = west = atof (xx1);		project_info.e = east = atof (xx2);
					if (project_info.w > 180.0 && project_info.e > 180.0) {
						project_info.w -= 360.0;
						project_info.e -= 360.0;
					}
					project_info.s = south = atof (yy1);	project_info.n = north = atof (yy2);
					found_proj = TRUE;
					if ((ptmp = strstr(&line[2], "+proj")) != NULL) {  /* Search for the +proj in the comment line */
						proj4_cmd = strdup(&line[(int)(ptmp - &line[0])]);
						GMT_chop (proj4_cmd);		/* Remove the new line char */
					}
					if (!strcmp(proj4_name,"latlong") || !strcmp(proj4_name,"xy") ||
						!strcmp(proj4_name,"eqc") ) {		/* Linear case, use original coords */
						west  = atof(xx1);		east  = atof(xx2);
						south = atof(yy1);		north = atof(yy2);
						/* One further test. +xy was found, but have we geog coords? Check that */
						if (!strcmp(proj4_name,"xy") && 
								(west >= -180) && ((east <= 360) && ((east - west) <= 360)) &&
								(south >= -90) && (north <= 90) ) {
							proj4_cmd = strdup("latlon");
							fprintf (stderr, "%s: WARNING: An unknown projection setting was found but since "
									"image coordinates seam to be geographical, a linear transformation "
									"will be used.\n", GMT_program);
						}
						else if (!strcmp(proj4_name,"xy") && Ctrl->W.warp) {	/* Do not operate on a twice unknown setting */
							fprintf (stderr, "%s: WARNERROR: You requested an automatic geotiff generation, but "
									"no recognized projection was used for the PS creation.\n", GMT_program); 
						}
					}
					else if (Ctrl->W.kml) {
						fprintf (stderr, "%s: WARNERROR: To GE images must be in geographical coords. Very likely "
									"this won't work as you wish inside GE.\n", GMT_program); 
					}
				}
			}
			sscanf (line, "%s",c);
			if (!strncmp(c, "%%BoundingBox:", (size_t)14)) {
				if (got_BB) fprintf (fpo, "%%%%BoundingBox: 0 0 %ld %ld\n", (GMT_LONG)ceil(w), (GMT_LONG)ceil(h));
				got_BB = FALSE;
				if (file_has_HRBB) continue;	/* High-res BB will be put elsewhere */
				if (got_HRBB) fprintf (fpo, "%%%%HiResBoundingBox: 0 0 %g %g\n", w, h);
				got_HRBB = FALSE;
				continue;
			}
			else if (!strncmp(c, "%%HiResBoundingBox:", (size_t)19)) {
				if (got_HRBB) fprintf (fpo, "%%%%HiResBoundingBox: 0 0 %g %g\n", w, h);
				got_HRBB = FALSE;
				continue;
			}
			else if (Ctrl->P.active && landscape && !strncmp(c, "%%Orientation:", (size_t)14)) {
				fprintf (fpo, "%%%%Orientation: Portrait\n");
				landscape = FALSE;
				continue;
			}
			else if (!strncmp(c, "%%BeginSetup", (size_t)12))
				setup=TRUE;
			else if (!strncmp(c, "%%EndSetup", (size_t)10))
				setup=FALSE;
			else if (!strncmp(c, "%%EndComments", (size_t)13)) {
				fprintf (fpo, "%s", line);
				if (r != 0) fprintf (fpo, "%ld rotate\n", r);
				if (!GMT_IS_ZERO(xt) || !GMT_IS_ZERO(yt)) fprintf (fpo, "%g %g translate\n", xt, yt);
				xt = yt = 0.0;
				r = 0;
				continue;
			}
#ifdef USE_GDAL
			else if (!strncmp(c, "%%PageTrailer", (size_t)13) && found_proj) {
				fgets2 (line, BUFSIZ, fp);
				fprintf (fpo, "%%%%PageTrailer\n");
				fprintf (fpo, "%s", line);

				/* Write a GeoPDF registration info */ 

				/* Allocate new control structures */
				to_gdalread = (struct GDALREAD_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GDALREAD_CTRL), "New_Gdalread_Ctrl");
				from_gdalread = (struct GD_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GD_CTRL), "New_Gd_Ctrl");
				to_gdalread->W.active = 1;
				from_gdalread->ProjectionRefPROJ4 = proj4_cmd;
				GMT_gdalread ( NULL, to_gdalread, from_gdalread);
				if (from_gdalread->ProjectionRefWKT != CNULL) {
					double x0, y0, x1;	/* Projected coordinates */
					double h0, v0, h1;	/* Correspnding point coordinates */
					double a, H, V;		/* a -> coeff of affine matrix; H,V -> origin shift in projected coords */
					double pX[4], pY[4], lptsX[4], lptsY[4];

					x0 = west_prj;		y0 = south_prj;		x1 = east_prj;		y1 = north_prj;
					h0 = v0 = 0;	/* because we used -A option so origin in points is at (0,0) */
					h1 = h0 + w;
					/* takes the projected coordinate system into the page coordinate system */
					a = (h1 - h0)/(x1 - x0);
					H = h0 - a*x0;
					V = v0 - a*y0;

					/* Use the above matrix to discover where the corners of the map will fall */
					pX[0] = a*x0 + H;	pX[1] = a*x0 + H;	pX[2] = a*x1 + H;	pX[3] = a*x1 + H;
					pY[0] = a*y0 + V;	pY[1] = a*y1 + V;	pY[2] = a*y1 + V;	pY[3] = a*y0 + V;

					/* Now compute the LPTS array */
					lptsX[0] = pX[0] / w;	lptsX[1] = pX[1] / w;	lptsX[2] = pX[2] / w;	lptsX[3] = pX[3] / w;
					lptsY[0] = pY[0] / h;	lptsY[1] = pY[1] / h;	lptsY[2] = pY[2] / h;	lptsY[3] = pY[3] / h;

					fprintf (fpo, "\n%% embed georegistation info\n");
					fprintf (fpo, "[ {ThisPage} <<\n");
					fprintf (fpo, "\t/VP [ <<\n");
					fprintf (fpo, "\t\t/Type /Viewport\n");
					fprintf (fpo, "\t\t/BBox[0 0 %.1f %.1f]\n", w, h);
					fprintf (fpo, "\t\t/Measure <<\n");
					fprintf (fpo, "\t\t\t/Type /Measure\n");
					fprintf (fpo, "\t\t\t/Subtype /GEO\n");
					fprintf (fpo, "\t\t\t/Bounds[0 0 0 1 1 1 1 0]\n");
					fprintf (fpo, "\t\t\t/GPTS[%f %f %f %f %f %f %f %f]\n", 
						south, west, north, west, north, east, south, east);
					fprintf (fpo, "\t\t\t/LPTS[%f %f %f %f %f %f %f %f]\n", 
						lptsX[0],lptsY[0], lptsX[1],lptsY[1], lptsX[2],lptsY[2], lptsX[3],lptsY[3]);
					fprintf (fpo, "\t\t\t/GCS <<\n");
					fprintf (fpo, "\t\t\t\t/Type /PROJCS\n");
					fprintf (fpo, "\t\t\t\t/WKT\n");
					fprintf (fpo, "\t\t\t\t(%s)\n", from_gdalread->ProjectionRefWKT);
					fprintf (fpo, "\t\t\t>>\n");
					fprintf (fpo, "\t\t>>\n");
					fprintf (fpo, "\t>>]\n");
					fprintf (fpo, ">> /PUT pdfmark\n\n");
				}
				GMT_free((void *)to_gdalread);
				GMT_free((void *)from_gdalread);
				continue;
			}
#endif
			fprintf (fpo, "%s", line);
		}

		fclose (fpo);
		fclose (fp);

		/* Build the converting ghostscript command and execute it */

		if (Ctrl->T.device != GS_DEV_EPS) {
			char tag[16];
			strcpy (tag, &ext[Ctrl->T.device][1]);
			GMT_str_toupper (tag);
			if (Ctrl->V.active) fprintf (stderr, " Convert to %s...", tag);

			if (!Ctrl->F.active) {
				if (Ctrl->D.active) sprintf (out_file, "%s/", Ctrl->D.dir);		/* Use specified output directory */
				strncat (out_file, &ps_file[pos_file], (size_t)(pos_ext - pos_file));
			}
			else
				strcpy (out_file, Ctrl->F.file);
			strcat (out_file, ext[Ctrl->T.device]);
			pix_w = (GMT_LONG)ceil (w * Ctrl->E.dpi / 72.0);
			pix_h = (GMT_LONG)ceil (h * Ctrl->E.dpi / 72.0);
			sprintf (cmd, "%s%s %s %s -sDEVICE=%s -g%ldx%ld -r%ld -sOutputFile=%s -f%s", 
				at_sign, Ctrl->G.file, gs_params, gs_extra, device[Ctrl->T.device],
				pix_w, pix_h, Ctrl->E.dpi, out_file, tmp_file);
			i_unused = system (cmd);		/* Execute the GhostScript command */
			if ((i_unused = access (out_file, R_OK)) != 0)
				fprintf(stderr, "\nPS2RASTER WARNING: file\n\t%s\nwas not created. Maybe you forgot to close the PS file (a -K in excess?)\n", out_file);
			if (Ctrl->S.active) fprintf (stdout, "%s\n", cmd);
		}
		if (Ctrl->V.active) fprintf (stderr, " Done.\n");

		if (!Ctrl->T.eps) remove (tmp_file);
		if (no_U_file) remove (no_U_file);
		if (clean_PS_file) remove (clean_PS_file);
		GMT_free ((void *)tmp_file);
		if (clean_PS_file) GMT_free ((void *)clean_PS_file);
		if (no_U_file) GMT_free ((void *)no_U_file);

		if ( Ctrl->W.active && found_proj && !Ctrl->W.kml ) {	/* Write a world file */
			double x_inc, y_inc;
			char *world_file, *wext, *s;

			x_inc = (east  - west)  / pix_w; 
			y_inc = (north - south) / pix_h; 
			if (Ctrl->V.active) fprintf(stderr, "width = %ld\theight = %ld\tX res = %f\tY res = %f\n", pix_w, pix_h, x_inc, y_inc);

			/* West and North of the world file contain the coordinates of the center of the pixel
			   but our current values are of the NW corner of the pixel (pixel registration). So 
			   we'll move halph pixel inward. */
			west  += x_inc / 2;
			north -= y_inc / 2;

			world_file = (char *) GMT_memory (VNULL, (size_t)BUFSIZ, sizeof(char), GMT_program);
			if (Ctrl->D.active) sprintf (world_file, "%s/", Ctrl->D.dir);	/* Use specified output directory */
			if (Ctrl->F.active) {		/* Must rip the raster file extension before adding the world one */
				for (i = (GMT_LONG)strlen(out_file) - 1; i > 0; i--) {
					if (out_file[i] == '.') { 	/* Beginning of file extension */
						pos_ext = i;
						break;
					}
				}
				out_file[pos_ext] = '\0';
				strcat(world_file, out_file);
			}
			else
				strncat (world_file, &ps_file[pos_file], (size_t)(pos_ext - pos_file));

			s = ext[Ctrl->T.device];
			wext = strdup(ext[Ctrl->T.device]);
			wext[1] = s[1];		wext[2] = s[3];		wext[3] = 'w';
			strcat (world_file, wext);

			if ((fpw = fopen (world_file, "w")) == NULL) {
				fprintf (stderr, "%s: Unable to open file %s for writing\n", GMT_program, world_file);
			}
			else {
				fprintf (fpw, "%.12f\n0.0\n0.0\n%.12f\n%.12f\n%.12f", x_inc, -y_inc, west, north);
				fclose(fpw);
				if (Ctrl->V.active) {
					fprintf(stderr, "Wrote world file %s\n", world_file);
					if (proj4_cmd) fprintf(stderr, "Proj4 definition: %s\n", proj4_cmd); 
				}
			
			}

			free ((void *)wext);	

			if (Ctrl->W.warp && proj4_cmd && proj4_cmd[1] == 'p') {	/* We got a usable Proj4 string. Run it (if gdal is around) */
				/* The true geotiff file will have the same base name plus a .tiff extension.
				   We will reuse the world_file variable because all it is need is to replace the extension */
				for (i = (GMT_LONG)strlen(world_file) - 1; i > 0; i--) {
					if (world_file[i] == '.') { 	/* Beginning of file extension */
						pos_ext = i;
						break;
					}
				}
				world_file[pos_ext] = '\0';
				strcat(world_file, ".tiff");

				if (!Ctrl->V.active)		/* Shut up the gdal_translate (low level) verbosity */
					quiet = " -quiet";
				else
					quiet = "";

#ifdef WIN32
				sprintf (cmd, "gdal_translate -a_srs \"%s\" -co COMPRESS=LZW -co TILED=YES %s %s %s", 
						proj4_cmd, quiet, out_file, world_file); 
#else
				sprintf (cmd, "gdal_translate -a_srs '%s' -co COMPRESS=LZW -co TILED=YES %s %s %s", 
						proj4_cmd, quiet, out_file, world_file); 
#endif
				free(proj4_cmd);
				i_unused = system (cmd);		/* Execute the gdal_translate command */
				if (Ctrl->V.active) fprintf(stderr, "\nThe gdal_translate command: \n%s\n", cmd);
			}
			else if (Ctrl->W.warp && !proj4_cmd)
				fprintf (stderr, "%s: Could not find the Proj4 command in the PS file. No conversion performed.\n", GMT_program);

			GMT_free ((void *)world_file);	
		}

		else if ( Ctrl->W.kml ) {	/* Write a basic kml file */
			char *kml_file;
			kml_file = (char *) GMT_memory (VNULL, (size_t)BUFSIZ, sizeof(char), GMT_program);
			if (Ctrl->D.active) sprintf (kml_file, "%s/", Ctrl->D.dir);	/* Use specified output directory */
			if (Ctrl->F.active) {		/* Must rip the raster file extension before adding the kml one */
				for (i = (GMT_LONG)strlen(out_file) - 1; i > 0; i--) {
					if (out_file[i] == '.') { 	/* Beginning of file extension */
						pos_ext = i;
						break;
					}
				}
				out_file[pos_ext] = '\0';
				strcat(kml_file, out_file);
				out_file[pos_ext] = '.';	/* Reset the extension */
			}
			else
				strncat (kml_file, &ps_file[pos_file], (size_t)(pos_ext - pos_file));

			strcat (kml_file, ".kml");

			if ((fpw = fopen (kml_file, "w")) == NULL) {
				fprintf (stderr, "%s: Unable to open file %s for writing\n", GMT_program, kml_file);
			}
			else {
				fprintf (fpw, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
				fprintf (fpw, "<kml xmlns=\"http://earth.google.com/kml/2.1\">\n");
				fprintf (fpw, "<Document>\n\t<name>%s</name>\n", Ctrl->W.doctitle);
				fprintf (fpw, "\t<GroundOverlay>\n\t\t<name>%s</name>\n", Ctrl->W.overlayname);
				fprintf (fpw, "\t\t<Icon>\n");
				fprintf (fpw, "\t\t\t<href>");
				if (Ctrl->W.URL) fprintf (fpw, "%s/", Ctrl->W.URL);
				fprintf (fpw, "%s</href>\n", out_file);
				fprintf (fpw, "\t\t</Icon>\n");
				fprintf (fpw, "\t\t<altitudeMode>%s</altitudeMode>\n", RefLevel[Ctrl->W.mode]);
				if (Ctrl->W.mode > KML_GROUND_ABS && Ctrl->W.mode < KML_SEAFLOOR_ABS)
					fprintf (fpw, "\t\t<altitude>%g</altitude>\n", Ctrl->W.altitude);
				fprintf (fpw, "\t\t<LatLonBox>\n");
				fprintf (fpw, "\t\t\t<north>%f</north>\n", project_info.n);
				fprintf (fpw, "\t\t\t<south>%f</south>\n", project_info.s);
				fprintf (fpw, "\t\t\t<east>%f</east>\n", project_info.e);
				fprintf (fpw, "\t\t\t<west>%f</west>\n", project_info.w);
				fprintf (fpw, "\t\t</LatLonBox>\n");
				fprintf (fpw, "\t\t<Region>\n");
				fprintf (fpw, "\t\t<LatLonAltBox>\n");
				fprintf (fpw, "\t\t\t<north>%f</north>\n", north);
				fprintf (fpw, "\t\t\t<south>%f</south>\n", south);
				fprintf (fpw, "\t\t\t<east>%f</east>\n", east);
				fprintf (fpw, "\t\t\t<west>%f</west>\n", west);
				fprintf (fpw, "\t\t</LatLonAltBox>\n");
				if (Ctrl->W.min_lod != Ctrl->W.max_lod) {	/* Control layer visibility */
 					fprintf (fpw, "\t\t<Lod>\n");
 					fprintf (fpw, "\t\t\t<minLodPixels>%ld</minLodPixels>\n", Ctrl->W.min_lod);
 					fprintf (fpw, "\t\t\t<maxLodPixels>%ld</maxLodPixels>\n", Ctrl->W.max_lod);
					if (Ctrl->W.min_fade) fprintf (fpw, "\t\t\t<minFadeExtent>%ld</minFadeExtent>\n", Ctrl->W.min_fade);
	 				if (Ctrl->W.max_fade) fprintf (fpw, "\t\t\t<maxFadeExtent>%ld</maxFadeExtent>\n", Ctrl->W.max_fade);
 					fprintf (fpw, "\t\t</Lod>\n");
				}
				fprintf (fpw, "\t\t</Region>\n");
				fprintf (fpw, "\t</GroundOverlay>\n");
				fprintf (fpw, "</Document>\n</kml>\n");
				fclose(fpw);
				if (Ctrl->V.active) {
					fprintf(stderr, "Wrote KML file %s\n", kml_file);
				}
			
			}
			GMT_free ((void *)kml_file);	
		}

		else if (Ctrl->W.active && !found_proj) {
			fprintf (stderr, "%s: Could not find the 'PROJ' tag in the PS file. No world file created.\n", GMT_program);
			fprintf (stderr, "This situation occurs when one of the two next cases is true:\n");
			fprintf (stderr, "1) the PS file was created with a pre-GMT v4.5.0 version\n");
			fprintf (stderr, "2) the PS file was not created by GMT\n");
		}
	}

	Free_ps2raster_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_free ((void *)ps_names);

	exit (EXIT_SUCCESS);
}

GMT_LONG parse_GE_settings (char *arg, struct PS2RASTER_CTRL *C)
{
	/* Syntax: -W[+g][+k][+t<doctitle>][+n<layername>][+a<altmode>][+f<fademin>/<fademax>][+l<lodmin>/<lodmax>][+u<URL>] */
	
	GMT_LONG error = FALSE;
	GMT_LONG pos = 0;
	char txt[BUFSIZ], p[BUFSIZ];
	
	C->W.active = TRUE;
	strcpy (txt, arg);
	while (!error && (GMT_strtok (txt, "+", &pos, p))) {
		switch (p[0]) {
			case 'a':	/* Altitude setting */
				switch (p[1]) {	/* Check which altitude mode we selected */
					case 'G':
						C->W.mode = KML_GROUND_ABS;
						break;
					case 'g':
						C->W.mode = KML_GROUND_REL;
						C->W.altitude = atof (&p[2]);
						break;
					case 'A':
						C->W.mode = KML_ABS;
						C->W.altitude = atof (&p[2]);
						break;
					case 's':
						C->W.mode = KML_SEAFLOOR_REL;
						C->W.altitude = atof (&p[2]);
						break;
					case 'S':
						C->W.mode = KML_SEAFLOOR_ABS;
						break;
					default:
						fprintf (stderr, "%s: GMT ERROR -W+a<mode>[par]: Unrecognized altitude mode %c\n", GMT_program, p[1]);
						error = TRUE;
						break;
				}
				break;
			case 'f':	/* Set fading options in KML */
				sscanf (&p[1], "%" GMT_LL "d/%" GMT_LL "d", &C->W.min_fade, &C->W.max_fade);
				break;
			case 'g':	/* Use gdal to make geotiff */
				C->W.warp = TRUE;
				break;
			case 'k':	/* Produce a KML file */
				C->W.kml = TRUE;
				break;
			case 'l':	/* Set KML level of detail for image */
				sscanf (&p[1], "%" GMT_LL "d/%" GMT_LL "d", &C->W.min_lod, &C->W.max_lod);
				break;
			case 'n':	/* Set KML document layer name */
				if (C->W.overlayname) free (C->W.overlayname);	/* Already set, free then reset */
				C->W.overlayname = strdup (&p[1]);
				break;
			case 't':	/* Set KML document title */
				if (C->W.doctitle) free (C->W.doctitle);	/* Already set, free then reset */
				C->W.doctitle = strdup (&p[1]);
				break;
			case 'u':	/* Specify a remote address for image */
				if (C->W.URL) free (C->W.URL);	/* Already set, free then reset */
				C->W.URL = strdup (&p[1]);
				break;
			default:
				fprintf (stderr, "%s: GMT ERROR -W+<opt>: Unrecognized option selection %c\n", GMT_program, p[1]);
				error = TRUE;
				break;
		}
	}
	return (error);
}

void *New_ps2raster_Ctrl () {	/* Allocate and initialize a new control structure */
	struct PS2RASTER_CTRL *C;

	C = (struct PS2RASTER_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PS2RASTER_CTRL), "New_ps2raster_Ctrl");

	/* Initialize values whose defaults are not 0/FALSE/NULL */
#ifdef WIN32
	C->G.file = strdup ("gswin32c");
#else
	C->G.file = strdup ("gs");
#endif
	C->D.dir = strdup (".");

	C->W.doctitle = strdup ("GMT KML Document");
	C->W.overlayname = strdup ("GMT Image Overlay");

	return ((void *)C);
}

void Free_ps2raster_Ctrl (struct PS2RASTER_CTRL *C) {	/* Deallocate control structure */
	if (C->D.dir) free ((void *)C->D.dir);
	if (C->F.file) free ((void *)C->F.file);
	if (C->G.file) free ((void *)C->G.file);
	if (C->L.file) free ((void *)C->L.file);
	free ((void *)C->W.doctitle);
	free ((void *)C->W.overlayname);
	if (C->W.URL) free ((void *)C->W.URL);
	GMT_free ((void *)C);
}
