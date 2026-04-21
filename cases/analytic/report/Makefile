# Makefile

NAME = main

INTERMEDIATES = $(NAME).log \
				$(NAME).aux \
				$(NAME).bbl \
				$(NAME).bcf \
				$(NAME).blg \
				$(NAME).out \
				$(NAME).run.xml \
				$(NAME).spl \
				$(NAME).acn \
				$(NAME).acr \
				$(NAME).alg \
				$(NAME).glg \
				$(NAME).glo \
				$(NAME).gls \
				$(NAME).ist

TABLES = ../summary_onegroup.tex ../summary_multigroup.tex

FIGS = 

LATEX = pdflatex
GLOSSARY = makeglossaries
PYTHON = python3

LATEXOPT = -halt-on-error -file-line-error -shell-escape

.PHONY : default all clean cleanpdf

default : all

all : $(NAME).pdf

$(NAME).pdf : $(NAME).tex $(TABLES)
	$(LATEX) $(LATEXOPT) -draftmode $(NAME)
	$(LATEX) $(LATEXOPT) -draftmode $(NAME)
	$(LATEX) $(LATEXOPT) $(NAME)
	@grep -i --color=always 'warn' $(NAME).log || true
	@grep -i --color=always 'overfull' $(NAME).log || true
	@grep -i --color=always 'underfull' $(NAME).log || true

clean :
	$(RM) $(INTERMEDIATES)

cleanpdf :
	$(RM) $(NAME).pdf
