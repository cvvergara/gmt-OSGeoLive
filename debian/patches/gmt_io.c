Index: git/src/gmt_io.c
===================================================================
--- git.orig/src/gmt_io.c	2013-11-25 15:57:08.000000000 +0100
+++ git/src/gmt_io.c	2013-11-25 16:37:51.000000000 +0100
@@ -2772,8 +2772,9 @@
 
 		callen = strlen (s);
 		if (callen < 2) return (GMT_IS_NAN);	/* Maybe should be more than 2  */
-
-		//if ((p = strchr ( s, (int)('T'))) == NULL) {	/* This was too naive, being tricked by data like 12-OCT-20 (no trailing T, so OCT was it) */
+#if 0
+		if ((p = strchr ( s, (int)('T'))) == NULL) {	/* This was too naive, being tricked by data like 12-OCT-20 (no trailing T, so OCT was it) */
+#endif
 		if (s[0] == 'T') {	/* Got T<clock> presumably */
 			clocklen = callen - 1;
 			strncpy (clockstring, &s[1], (size_t)callen);
