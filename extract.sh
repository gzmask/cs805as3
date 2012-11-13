sh compile.sh
bin/run
rawtopgm 512 512 output.raw > tmp.pgm
ppmtojpeg < tmp.pgm > result.jpg
exit 0
