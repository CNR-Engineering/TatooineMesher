
all: TUC2019_slides_Duron-et-al_TatooineMesher.pdf

%.pdf: %.tex
	cd $(dir $<) && xelatex $(notdir $<) && xelatex $(notdir $<)

clean:
	find -name '*.aux' -delete
	find -name '*.bbl' -delete
	find -name '*.blg' -delete
	find -name '*.log' -delete
	find -name '*.out' -delete
	find -name '*.toc' -delete
