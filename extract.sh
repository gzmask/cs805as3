noweb -t as2.nw
cp as2.nw latex/as2.tex
cd latex
texi2pdf as2.tex
rm as2.aux as2.log as2.out as2.tex
cd ..
sh compile.sh
bin/run
rawtopgm 512 512 output.raw > tmp.pgm
ppmtojpeg < tmp.pgm > result.jpg
exit 0
