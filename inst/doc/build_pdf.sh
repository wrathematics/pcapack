#!/bin/sh

rm *.aux *.bbl *.blg *.log *.out *.toc
pdflatex pbdPROF-guide.Rnw
pdflatex pbdPROF-guide.Rnw
rm *.aux *.bbl *.blg *.log *.out *.toc *.dvi
