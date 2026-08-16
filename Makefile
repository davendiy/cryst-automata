
.DEFAULT_GOAL := help

help:
	@echo "planar-latex-generate planar-latex-compile planar-latex article2025"

planar-latex-generate:
	python working/planar_srdegrees.py > latex/planar_groups/planar_groups.tex

planar-latex-compile:
	pdflatex --output-directory=latex/planar_groups latex/planar_groups/planar_groups.tex

article2025:
	pdflatex --output-directory=latex/computing-srdegrees latex/computing-srdegrees/article.tex

article2026:
	pdflatex --output-directory=latex/ssdegrees latex/ssdegrees/ssarticle.tex

a: article2025
a2: article2026

planar-latex: planar-latex-generate planar-latex-compile
