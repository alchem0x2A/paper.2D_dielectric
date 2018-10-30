IMG_PATH=img
BIB_FILE=ref.bib
TEX_FILE=paper.tex
SI_FILE=SI.tex
GS_TAG= -sDEVICE=pdfwrite -dQUIET -sBATCH -dNOPAUSE
LATEXMK_TAG= -f -pdf -quiet -view=none -pdflatex='pdflatex -interaction=nonstopmode'
DIFF_TAG= --exclude-textcmd="section,subsection,figure" --config="PICTUREENV=(?:section|DIFnomarkup)[*]*" --graphics-markup=0 --disable-citation-markup



all: main SI diff

main:
	latexmk $(LATEXMK_TAG) $(TEX_FILE)

SI:
	latexmk $(LATEXMK_TAG) $(SI_FILE)

diff:
	cp ./paper.tex ./paper.tex.bak
	git checkout 45e0f99 ./paper.tex #old version
	cp paper.tex paper_old.tex
	cp ./SI.tex ./SI.tex.bak	
	git checkout 45e0f99 ./SI.tex
	cp SI.tex SI_old.tex
	git checkout HEAD ./paper.tex
	git checkout HEAD ./SI.tex
	cp ./paper.tex.bak paper.tex #restore to current mode
	cp ./SI.tex.bak SI.tex
	latexdiff $(DIFF_TAG) paper_old.tex paper.tex > paper_change.tex
	latexdiff $(DIFF_TAG) SI_old.tex SI.tex > SI_change.tex
	dos2unix paper_change.tex
	dos2unix SI_change.tex
	sed -e "s///" paper_change.tex > tmp.tex && mv tmp.tex paper_change.tex
	sed -e "s///" SI_change.tex > tmp.tex && mv tmp.tex SI_change.tex
	latexmk $(LATEXMK_TAG) paper_change.tex
	latexmk $(LATEXMK_TAG) SI_change.tex

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
