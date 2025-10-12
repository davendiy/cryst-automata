

planar-latex-generate:
	python working/planar_srdegrees.py > latex/planar_groups/planar_groups.tex

planar-latex-compile:
	pdflatex -output-directory=latex/planar_groups latex/planar_groups/planar_groups.tex

planar-latex: planar-latex-generate planar-latex-compile
