PRINCIPALE = main
PRINCIPALE_TEX = $(PRINCIPALE).tex
PRINCIPALE_PDF = $(PRINCIPALE).pdf

PRESE_NAME = Bucci_Cipriani_Pagella_Petruso_Puricelli_Venturini.pdf
PRESE1_DIRECTORY = slides/presentation-2022-11-11
PRESE2_DIRECTORY = slides/presentation-2023-01-09
PRESE3_DIRECTORY = slides/presentation-2023-02-14

REPORT_NAME = Stochastic_Block_Model_Prior_with_Ordering_Constraints_for_Gaussian_Graphical_Models.pdf
REPORT_DIRECTORY = report

.PHONY: prese1 prese2 prese3 report clean distclean pdf

pdf: prese1 prese2 prese3 report

# -----------------------------------------------------------------------------

prese1: $(PRESE1_DIRECTORY)/$(PRESE_NAME) # prese1 è il comando per richiedere il pdf

$(PRESE1_DIRECTORY)/$(PRESE_NAME): $(PRESE1_DIRECTORY)/*.tex # l'output in pdf della presentazione dipende dai tex nella cartella
	cd $(PRESE1_DIRECTORY)\
	&& latexmk -lualatex $(PRINCIPALE_TEX)\
	&& mv $(PRINCIPALE_PDF) $(PRESE_NAME) # se qualcosa è cambiato, ri-builda il documento pdf

# -----------------------------------------------------------------------------

prese2: $(PRESE2_DIRECTORY)/$(PRESE_NAME)

$(PRESE2_DIRECTORY)/$(PRESE_NAME): $(PRESE2_DIRECTORY)/*.tex
	cd $(PRESE2_DIRECTORY)\
	&& latexmk -lualatex $(PRINCIPALE_TEX)\
	&& mv $(PRINCIPALE_PDF) $(PRESE_NAME)

# -----------------------------------------------------------------------------

prese3: $(PRESE3_DIRECTORY)/$(PRESE_NAME)

$(PRESE3_DIRECTORY)/$(PRESE_NAME): $(PRESE3_DIRECTORY)/*.tex
	cd $(PRESE3_DIRECTORY)\
	&& latexmk -lualatex $(PRINCIPALE_TEX)\
	&& mv $(PRINCIPALE_PDF) $(PRESE_NAME)

# -----------------------------------------------------------------------------

report: $(REPORT_DIRECTORY)/$(REPORT_NAME)

$(REPORT_DIRECTORY)/$(REPORT_NAME): $(REPORT_DIRECTORY)/*.tex
	cd $(REPORT_DIRECTORY)\
	&& git show -s --format=%H > commit_hash.part\
	&& pdflatex -shell-escape $(PRINCIPALE_TEX)\
	&& rm -f commit_hash.part\
	&& mv $(PRINCIPALE_PDF) $(REPORT_NAME)

# -----------------------------------------------------------------------------

clean:
	find . \( -iname "*.aux" -o -iname "*.log" -o -iname "*.out" -o -iname "*.xdv" -o -iname "*.toc" -o -iname "*.fls" -o -iname "*.fls" -o -iname "*.fdb_latexmk" -o -iname "*.synctex.gz" -o -iname "*.synctex\(busy\)" -o -iname "*.bbl" -o -iname "*.bcf" -o -iname "*.blg" -o -iname "*.nav" -o -iname "*.run.xml" -o -iname "*.snm" -o -iname "*.dvi" -o -iname "*.ps" -o -iname "*commit_hash.part" \) -type f -delete

distclean: clean
	find . \( -iname "*$(PRESE_NAME)" -o -iname "*$(REPORT_NAME)" \) -type f -delete