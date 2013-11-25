Index: git/src/gmt_proj.c
===================================================================
--- git.orig/src/gmt_proj.c	2013-11-25 15:57:08.000000000 +0100
+++ git/src/gmt_proj.c	2013-11-25 16:45:37.000000000 +0100
@@ -1645,7 +1645,9 @@
 	project_info.g_outside = FALSE;
 
 	angle = M_PI - dlon;
-	// if (cosc < project_info.g_P_inverse) { /* over the horizon, but this was susceptible to minor roundoff. */
+#if 0
+	if (cosc < project_info.g_P_inverse) { /* over the horizon, but this was susceptible to minor roundoff. */
+#endif
 	if ((project_info.g_P_inverse - cosc) > GMT_CONV_LIMIT) { /* over the horizon */
 		project_info.g_outside = TRUE;
 
