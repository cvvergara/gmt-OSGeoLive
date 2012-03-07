/*--------------------------------------------------------------------
 *	$Id: gmt2bin.c,v 1.23 2011/07/11 19:22:05 guru Exp $
 *
 *    Copyright (c) 1991-2011 by P. Wessel and W. H. F. Smith
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * gmt2bin
 * program to read a leg list, open legid.gmt files, and determine what 1 x 1
 * degree bins of the earth are occupied by each leg and what data types exist
 * in each bin.  For each leg input, program writes a legid.bix file.
 * Note:  The program asks at the start if the legs in the leg list are inside
 * or outside legs; all legs in the list MUST be one or the other.  We flag the
 * bin_index file with this for each cruise in this session.
 *
 *	Input file :  a legid.gmt file, direct, recl = 18, the standard file
 *			of W. H. F. Smith and P. Wessel.  See makegmt.f
 *
 *	Output file:  a legid.bix file, as follows:
 *		first record: year1, total_gmt_flags (sum of 1(G), 2(M), 4(T)
 *		remaining recs: bindex gmt_flags 
 *
 *		year1 is the year of the first point in the file
 *		total_gmt_flags for existence of data over whole cruise.
 *
 *	New version 2.0	13-Jun-1991	P. Wessel
 *	Version 3.2, 12-APR-1999 (DOS compatible)
 */

#include "gmt.h"
#include "gmt_mgg.h"

#define MEM_SIZE 100000
#define MEM_INC 50000
#define MBINS 64800
#define DEG2MDEG   1000000
#define MIL_90   90000000
#define MIL_360 360000000

int main (int argc, char **argv)
{
	int year, n_records, bindex, error = FALSE, i, lastb, llat1, llat2, llon1, llon2, newbin, flag;
	size_t not_used = 0;
	double xedge, yedge, xx1, xx2, yy1, yy2, delx, dely, slope, xedgey, yedgex, rxedge, ryedge;
	char legid[16], agency[10], gmtfile[80], binfile[80];
	int inboxn[MBINS], inboxgmt[MBINS], totalgmt;
	char line[BUFSIZ], path[80];
	FILE *fpl = NULL, *fpo = NULL, *fp = NULL;
	struct GMTMGG_REC *record = NULL;

	argc = (int)GMT_begin (argc, argv);
	
	path[0] = 0;
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'V':
					gmtdefs.verbose = TRUE;
					break;
				case 'P':
					strcpy (path, &argv[i][2]);
					break;
				default:
					error = TRUE;
					fprintf (stderr, "GMT SYNTAX ERROR:  Unrecognized option -%c\n", argv[i][1]);
					break;
			}
		}
		else if ((fpl = fopen (argv[i], "r")) == NULL) {
			fprintf (stderr, "gmt2bin: Couldn't open file %s\n", argv[i]);
			exit (EXIT_FAILURE);
		}
	}
	
	if (argc == 1 || error) {
		fprintf (stderr, "usage: gmt2bin leglist [-Ppath -V]\n");
		fprintf (stderr, "	-P full path to directory with gmtfiles\n");
		fprintf (stderr, "	-V means verbose\n");
		exit (EXIT_FAILURE);
	}
	
	record = (struct GMTMGG_REC *) GMT_memory (VNULL, (size_t)1, sizeof (struct GMTMGG_REC), "gmt2bin");
		
	while (fgets (line, BUFSIZ, fpl)) {
		sscanf (line, "%s", legid);
		if (gmtdefs.verbose) fprintf (stderr, "Processing leg %s\n", legid);
		if (path[0])
			sprintf (gmtfile, "%s/%s.gmt", path, legid);
		else
			sprintf (gmtfile, "%s.gmt", legid);
		if ((fp = fopen (gmtfile, "rb")) == NULL) {
			fprintf (stderr, "gmt2bin: Cannot open file %s\n", gmtfile);
			continue;
		}
		not_used = fread ((void *)&year, (size_t)4, (size_t)1, fp);
		not_used = fread ((void *)&n_records, (size_t)4, (size_t)1, fp);
		not_used = fread ((void *)agency, (size_t)10, (size_t)1, fp);
		
		record = (struct GMTMGG_REC *) GMT_memory ((void *)record, (size_t)n_records, sizeof (struct GMTMGG_REC), "gmt2bin");
		
		i = 0;
		while (fread ((void *)&record[i], (size_t)18, (size_t)1, fp)) {
			
			while (record[i].lon < 0) record[i].lon += MIL_360;
			if (record[i].lon > MIL_360) continue;
			if (record[i].lat < -MIL_90 || record[i].lat > MIL_90) continue;
			record[i].lat += MIL_90;
			i++;
		}
		fclose (fp);
		n_records = i;
		
		totalgmt = 0;
		for (i = 0; i < MBINS; i++) inboxn[i] = inboxgmt[i] = 0;
		
		/*
		 * find occupied boxes, and also 'implied' occupancy for later xover checks.
		 * an implied box is one whose corner was cut by a ship going from one bin to
		 * another; we find these only when the deltat between points is < 15 minutes,
		 * so crossover interpolation applies.  This should arise only rarely.  The 15
		 * minute stipulation also means that in such cases there should be only one
		 * implied box; the following code tries to find only one.
		 */
		
		bindex = (record[0].lat / DEG2MDEG) * 360 + record[0].lon / DEG2MDEG;
		if (bindex >= MBINS) {
			fprintf (stderr, "gmt2bin: Index too big! (%d)\n", bindex);
			exit (EXIT_FAILURE);
		}
		lastb = bindex;
		inboxn[bindex] = 1;
		
		flag = 0;
		if (record[0].gmt[0] != GMTMGG_NODATA) flag |= 1;
		if (record[0].gmt[1] != GMTMGG_NODATA) flag |= 2;
		if (record[0].gmt[2] != GMTMGG_NODATA) flag |= 4;
		totalgmt |= flag;
		inboxgmt[bindex] |= flag;
		
		for (i = 1; i < n_records; i++) {
			bindex = (record[i].lat / DEG2MDEG) * 360 + record[i].lon / DEG2MDEG;
			if (bindex >= MBINS) {
				fprintf (stderr, "gmt2bin: Index too big! (%d)\n", bindex);
				exit (EXIT_FAILURE);
			}
			inboxn[bindex] = 1;
			flag = 0;
			if (record[i].gmt[0] != GMTMGG_NODATA) flag |= 1;
			if (record[i].gmt[1] != GMTMGG_NODATA) flag |= 2;
			if (record[i].gmt[2] != GMTMGG_NODATA) flag |= 4;
			totalgmt |= flag;
			inboxgmt[bindex] |= flag;
			
			if (lastb != bindex && (record[i].time - record[i-1].time) <= 900) {
				yy1 = record[i-1].lat * MDEG2DEG;
				yy2 = record[i].lat * MDEG2DEG;
 				xx1 = record[i-1].lon * MDEG2DEG;
				xx2 = record[i].lon * MDEG2DEG;
				llat1 = (int)yy1;
				llat2 = (int)yy2;
				llon1 = (int)xx1;
				llon2 = (int)xx2;
				if (llat1 != llat2 && llon1 != llon2) {
					dely = yy2 - yy1;
					yedge = (dely > 0.0) ? yy2 : yy1;
					delx = xx2 - xx1;
					if (delx > 180.0) {
						delx -= 360.0;
						xx1 += 360.0;
						xedge = xx1;
					}
					else if (delx < -180.0){
						delx += 360.0;
						xx2 += 360.0;
						xedge = xx2;
					}
					else
						xedge = (delx > 0.0) ? xx2 : xx1;
					slope = dely / delx;
					xedgey = slope * (xedge - xx1) + yy1;
					yedgex = ( (yedge - yy1) + slope * xx1 ) / slope;
					rxedge = hypot (xedge - xx1, xedgey - yy1);
					ryedge = hypot (yedgex - xx1, yedge - yy1);
					if (rxedge > ryedge)
						newbin = (record[i].lat / DEG2MDEG) * 360 + record[i-1].lon / DEG2MDEG;
					else
						newbin = (record[i-1].lat / DEG2MDEG) * 360 + record[i].lon / DEG2MDEG;
					if (newbin >= MBINS) {
						fprintf (stderr, "gmt2bin: Index too big! (%d)\n", bindex);
						exit (EXIT_FAILURE);
					}
					inboxn[newbin] = 1;
				}
			}
			lastb = bindex;
		}
	
		/*
		 * new feature:
		 */
	 
		if (totalgmt == 0)
			fprintf (stderr, "%s has empty data fields\n", gmtfile);
		else {
			sprintf (binfile, "%s.bix", legid);
			fpo = fopen (binfile, "w");
			fprintf (fpo, "%d %d\n", year, totalgmt);
			for (i = 0; i < MBINS; i++) if (inboxn[i]) fprintf (fpo, "%d %d\n", i, inboxgmt[i]);
			fclose (fpo);
		}
	}
	fclose (fpl);
	
	GMT_free ((void *)record);

	GMT_end (argc, argv);
	
	exit (EXIT_SUCCESS);
}
