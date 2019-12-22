IMG_PATH=img
BIB_FILE=ref.bib
TEX_FILE=paper.tex
SI_FILE=SI.tex
GS_TAG= -sDEVICE=pdfwrite -dQUIET -sBATCH -dNOPAUSE
LATEXMK_TAG= -f -pdf -quiet -view=none -pdflatex='pdflatex -interaction=nonstopmode'
DIFF_TAG= --exclude-textcmd="section,subsection,figure,equation,subequation" --config="PICTUREENV=(?:section|DIFnomarkup)[*]*" --graphics-markup=0 --disable-citation-markup
PDFLATEX_TAG= -interaction=nonstopmode -draftmode
TAR_ROOT=submit


all: main SI

convert:
	python3 convert.py

main:
	pdflatex $(PDFLATEX_TAG) $(SI_FILE) #make SI.aux
	latexmk $(LATEXMK_TAG) $(TEX_FILE)

SI:
	pdflatex $(PDFLATEX_TAG) $(TEX_FILE) #make SI.aux	
	latexmk $(LATEXMK_TAG) $(SI_FILE)

diff:
	cp ./paper.tex ./paper.tex.bak
	git checkout 0ca1340 ./paper.tex #old version
	cp paper.tex paper_old.tex
	cp ./SI.tex ./SI.tex.bak	
	git checkout 0ca1340 ./SI.tex
	cp SI.tex SI_old.tex
	git checkout HEAD ./paper.tex
	git checkout HEAD ./SI.tex
	cp ./paper.tex.bak paper.tex #restore to current mode
	cp ./SI.tex.bak SI.tex
	latexdiff $(DIFF_TAG) paper_old.tex paper.tex > paper_change.tex
	latexdiff $(DIFF_TAG) SI_old.tex SI.tex > SI_change.tex
	sed -e "s///" paper_change.tex > tmp.tex && mv tmp.tex paper_change.tex
	sed -e "s///" SI_change.tex > tmp.tex && mv tmp.tex SI_change.tex
	echo "Now please manually run the compilation for changed LaTeX files!"
	latexmk $(LATEXMK_TAG) paper_change.tex
	latexmk $(LATEXMK_TAG) SI_change.tex

compress:
	gs $(GS_TAG) -sOutputFile="paper-compak.pdf" paper.pdf
	gs $(GS_TAG) -sOutputFile="SI-compak.pdf" SI.pdf

transfer:
	cp $(TEX_FILE) $(SI_FILE) $(BIB_FILE) DOS_figs.tex raw_data.tex BS_figures.tex ./collab
	cp paper*.pdf SI*.pdf ./collab
	rsync -ahvz --exclude="*.png" --exclude="*converted-to.pdf" img collab/

tar:
	echo "Archiving only for arXiv submission!"
	rm -rf $(TAR_ROOT)
	rm $(TAR_ROOT).zip
	mkdir $(TAR_ROOT)/
	cp $(TEX_FILE) $(SI_FILE) $(BIB_FILE) raw_data.tex paper.aux SI.aux $(TAR_ROOT)/
	cp 00README.XXX SI.bbl paper.bbl $(TAR_ROOT)/
	mv $(TAR_ROOT)/paper.tex $(TAR_ROOT)/main.tex
	mv $(TAR_ROOT)/paper.bbl $(TAR_ROOT)/main.bbl
	mv $(TAR_ROOT)/SI.tex $(TAR_ROOT)/suppl.tex
	mv $(TAR_ROOT)/SI.bbl $(TAR_ROOT)/suppl.bbl 
	mkdir $(TAR_ROOT)/img/
	cp img/*.pdf $(TAR_ROOT)/img/	
	zip -r $(TAR_ROOT).zip $(TAR_ROOT)/*

Clean:
	latexmk -c


clean-all:
	latexmk -C
	rm -f $(IMG_PATH)/*.pdf	#clean up pdf
