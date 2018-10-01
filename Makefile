IMG_PATH=img
BIB_FILE=ref.bib
TEX_FILE=paper.tex
SI_FILE=SI.tex
GS_TAG= -sDEVICE=pdfwrite -dQUIET -sBATCH -dNOPAUSE



all: main SI

main: $(TEX_FILE) $(BIB_FILE)
	latexmk -f -pdf -quiet -view=none -pdflatex='xelatex -interaction=nonstopmode --shell-escape' $(TEX_FILE)

SI: $(SI_FILE) $(BIB_FILE) DOS_figs.tex raw_data.tex BS_figures.tex 
	latexmk -f -pdf -quiet -view=none -pdflatex='xelatex -interaction=nonstopmode --shell-escape' $(SI_FILE)

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
	clean
	latexmk -C
	rm -f $(IMG_PATH)/*.pdf	#clean up pdf
