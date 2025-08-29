

planar-latex:
	python working/planar_srdegrees.py > latex/planar_groups/planar_groups.tex
	pdflatex -output-directory=latex/planar_groups latex/planar_groups/planar_groups.tex
