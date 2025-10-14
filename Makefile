
.DEFAULT_GOAL := help

help:
	@echo "planar-latex-generate planar-latex-compile planar-latex article2025"

planar-latex-generate:
	python working/planar_srdegrees.py > latex/planar_groups/planar_groups.tex

planar-latex-compile:
	pdflatex -output-directory=latex/planar_groups latex/planar_groups/planar_groups.tex

article2025:
	pdflatex --output-directory=latex/computing-srdegrees latex/computing-srdegrees/article.tex

planar-latex: planar-latex-generate planar-latex-compile
