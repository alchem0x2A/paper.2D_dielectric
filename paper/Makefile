BIB_FILE="ref.bib"
TEX_FILE="2D_dielectric.tex"
PDF_FILE="2D_dielectric.pdf"
all: compile clean

compile: 2D_dielectric.tex ref.bib
	latexmk -f -pdf -quiet -view=none -pdflatex='xelatex -interaction=nonstopmode --shell-escape' $(TEX_FILE)
	cp $(PDF_FILE) ../

clean:
	latexmk -c

clean-all:
	latexmk -C
