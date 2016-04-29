#! /bin/bash
# ./quadrature_rule-3 9 20 && gnuplot quadrature_rule-3_draw.gnuplot && kpdf quadrature_rule-3_output.pdf 

./quadrature_rule-3 9 10 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qra.pdf
./quadrature_rule-3 9 20 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qrb.pdf
./quadrature_rule-3 9 30 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qrc.pdf
./quadrature_rule-3 9 40 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qrd.pdf
cp qra.pdf qraa.pdf
montage qr{a,b,c,d}.pdf -density 600x600 -geometry 1500x1500+5+5 qr_1.png

./quadrature_rule-3 9 50 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qra.pdf
./quadrature_rule-3 9 60 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qrb.pdf
./quadrature_rule-3 9 70 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qrc.pdf
./quadrature_rule-3 9 80 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qrd.pdf
cp qrb.pdf qrbb.pdf
montage qr{a,b,c,d}.pdf -density 600x600 -geometry 1500x1500+5+5 qr_2.png

./quadrature_rule-3 9 90 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qra.pdf
./quadrature_rule-3 9 100 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qrb.pdf
./quadrature_rule-3 9 110 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qrc.pdf
./quadrature_rule-3 9 120 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qrd.pdf
cp qrc.pdf qrcc.pdf
montage qr{a,b,c,d}.pdf -density 600x600 -geometry 1500x1500+5+5 qr_3.png


./quadrature_rule-3 9 130 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qra.pdf
./quadrature_rule-3 9 140 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qrb.pdf
./quadrature_rule-3 9 150 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qrc.pdf
./quadrature_rule-3 9 160 && gnuplot quadrature_rule-3_draw.gnuplot && cp quadrature_rule-3_output.pdf qrd.pdf
cp qra.pdf qrdd.pdf
montage qr{a,b,c,d}.pdf -density 600x600 -geometry 1500x1500+5+5 qr_4.png

montage qr{aa,bb,cc,dd}.pdf -density 600x600 -geometry 1500x1500+5+5 quadrature_rules.png
