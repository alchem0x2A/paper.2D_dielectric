IMG_PATH=img
BIB_FILE=ref.bib
TEX_FILE=paper.tex
SI_FILE=SI.tex
GS_TAG= -sDEVICE=pdfwrite -dQUIET -sBATCH -dNOPAUSE



all: main SI compress transfer

main:
	latexmk -f -pdf -quiet -view=none -pdflatex='pdflatex -interactive=nonstopmode' $(TEX_FILE)

SI:
	latexmk -f -pdf -quiet -view=none -pdflatex='pdflatex -interactive=nonstopmode' $(SI_FILE)

compress:
	gs $(GS_TAG) -sOutputFile="paper-compak.pdf" paper.pdf
	gs $(GS_TAG) -sOutputFile="SI-compak.pdf" SI.pdf

transfer:
	cp $(TEX_FILE) $(SI_FILE) $(BIB_FILE) DOS_figs.tex raw_data.tex BS_figures.tex ./collab
	cp paper-compak.pdf ./collab/paper.pdf
	cp SI-compak.pdf ./collab/SI.pdf
	rsync -ahvz --exclude="*.png" --exclude="*converted-to.pdf" img collab/

Clean:
	latexmk -c


clean-all:
	latexmk -C
	rm -f $(IMG_PATH)/*.pdf	#clean up pdf
