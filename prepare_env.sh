#!/bin/bash
CITIES="http://download.geonames.org/export/dump/cities1000.zip"

set -e

mkdir cities

TMPFILE=`mktemp`
wget $CITIES -O $TMPFILE
unzip -d cities $TMPFILE && rm $TMPFILE

mkdir cache
