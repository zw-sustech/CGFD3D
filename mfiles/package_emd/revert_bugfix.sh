#!/bin/sh

SUBDIR=EMDs/src/
MAINFILES='emdc.c emdc_fix.c cemdc.c cemdc_fix.c cemdc2.c cemdc2_fix.c'
for file in $MAINFILES
do
sed -i 's|#define _ALT_MEXERRMSGTXT_|/\*#define _ALT_MEXERRMSGTXT_\*/|' $SUBDIR$file
done
