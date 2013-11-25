Index: git/src/mgd77/mgd77header.c
===================================================================
--- git.orig/src/mgd77/mgd77header.c	2013-11-25 15:57:08.000000000 +0100
+++ git/src/mgd77/mgd77header.c	2013-11-25 17:19:37.000000000 +0100
@@ -759,7 +759,9 @@
 	not_used = fgets (line, BUFSIZ, F->fp);		/* Skip the column header  */
     
 	MGD77_header = (char *)GMT_memory (VNULL, (size_t)MGD77T_HEADER_LENGTH, sizeof (char), GMT_program);
-	// not_used = fgets (MGD77_header, BUFSIZ, F->fp);			/* Read the entire header record  */
+#if 0
+	not_used = fgets (MGD77_header, BUFSIZ, F->fp);			/* Read the entire header record  */
+#endif
     
 	for (i = 0; i < 2; i++) H->mgd77[i] = (struct MGD77_HEADER_PARAMS *) GMT_memory (VNULL, (size_t)1, sizeof (struct MGD77_HEADER_PARAMS), GMT_program);	/* Allocate parameter header */
     
