TEXS = Sediment_model.tex

all: Sediment_model.pdf

#cover.pdf: cover.tex
#latex cover.tex; dvipdf cover.dvi

Sediment_model.pdf: $(TEXS)  Sediment_model.bib
	pdflatex Sediment_model; bibtex Sediment_model; pdflatex Sediment_model; pdflatex Sediment_model

#report.pdf: report.tex paper.pdf cover.pdf doc_pdf/Appendix.pdf
#pdflatex report.tex;


clean: 
	 \rm -f *.bbl *.dvi *.log *.aux *.blg *.out #Sediment_model.pdf  
spotless: 
	 \rm -f *.bbl *.dvi *.log *.aux *.blg *.out Sediment_model.pdf  

