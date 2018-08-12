IMG_PATH=img
BIB_FILE=ref.bib
TEX_FILE=paper.tex
SI_FILE=SI.tex

all: convert_pdf compile clean

convert_pdf:
	python convert.py

compile: $(TEX_FILE) $(BIB_FILE)
	latexmk -f -pdf -quiet -view=none -pdflatex='xelatex -interaction=nonstopmode --shell-escape' $(TEX_FILE)
	cp $(PDF_FILE) ./build



clean:
	latexmk -c


clean-all:
	clean
	latexmk -C
	rm -f $(IMG_PATH)/*.pdf	#clean up pdf
