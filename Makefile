PRINCIPALE = main
PRINCIPALE_TEX = $(PRINCIPALE).tex
PRINCIPALE_PDF = $(PRINCIPALE).pdf

FILE_CLEAN = *.aux *.log *.out *.xdv *.toc *.fls *.fls *.fdb_latexmk *.synctex.gz *.synctex\(busy\) *.bbl *.bcf *.blg *.nav *.run.xml *.snm
#FILE_DISTCLEAN =  $(PRINCIPALE_PDF)

PRESE_NAME = Bucci_Cipriani_Pagella_Petruso_Puricelli_Venturini.pdf # nome standard di tutte le presentazioni
PRESE1_DIRECTORY = slides/presentation-2022-11-11

.PHONY: prese1 # distclean clean pdf

#pdf: $(PRINCIPALE_PDF)

prese1: $(PRESE1_DIRECTORY)/$(PRESE_NAME) # prese1 è il comando per richiedere il pdf

$(PRESE1_DIRECTORY)/$(PRESE_NAME): $(PRESE1_DIRECTORY)/*.tex # l'output in pdf della presentazione dipende dai tex nella cartella
	cd $(PRESE1_DIRECTORY) && latexmk -lualatex $(PRINCIPALE_TEX) && mv $(PRINCIPALE_PDF) $(PRESE_NAME) # se qualcosa è cambiato, ri-builda il documento pdf

# $(PRINCIPALE_PDF): $(PRINCIPALE_TEX) *.tex
# 	#git show -s --format=%H > commit_hash.part
# 	latexmk -lualatex $(PRINCIPALE_TEX) #-jobname=$(PRINCIPALE)
# 	mv $(PRINCIPALE_PDF) $(PRESE_NAME)
# 	#rm -f commit_hash.part

# clean:
# 	rm -f $(FILE_CLEAN) #commit_hash.part

# distclean: clean
# 	rm -f $(FILE_DISTCLEAN)