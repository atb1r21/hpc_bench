Joseph Prentice, September 2020

=============
BASIC OUTLINE
=============

Documentation is stored in .rst format (reStructuredText). Much of the
documentation has been converted from .tex (and some from .docx)
format into .rst format using pandoc (https://pandoc.org), followed by
a certain amount of tweaking of the resulting .rst files.

The resulting set of .rst files are then combined and converted into a
variety of formats using the Sphinx Python3 package
(https://www.sphinx-doc.org/en/master/). This can generate a combined
HTML, TeX/PDF, man pages, plain tex, etc. To generate a new version of
the documentation you will need Sphinx installed (v1.6.7 was used to
initially generate the documentation). In this directory, you then
simply run

make latexpdf

if you want the TeX/PDF documentation. The resulting .tex and .pdf
files will be found in the _build/latex directory. For a new version
of the HTML documentation, you simply run

make html

and the resulting files will be found in _build/html. (You might need
to run this twice to get all the references correct.) index.html is
the core file that links to all the others.

============
CONTRIBUTING
============

Editing existing pages should be fairly simple -- the .rst format is
well-documented elsewhere. Equations can be written using standard TeX
(no fancy packages).

If you want to contribute a new page, you can write it natively in
.rst (preferred), or you can convert a .tex file using pandoc. A
example command for doing this would be

pandoc EMFT_in_ONETEP.tex -f latex -t rst -s -o EMFT_in_ONETEP.rst

The resulting .rst file will likely need some tweaking -- common
issues are:

1. Unknown LaTeX commands (likely caused by using fancy packages)

2. Intendation immediately after an equation (should be no indentation
for text immediately after an equation, otherwise the code thinks it
is part of the equation)

3. Labels and referencing of equations (see existing .rst files for
examples)

4. Citations of papers (see existing .rst files for examples)

You will likely want to place your contribution within a subsection of
the documentation. Each subsection has its own index_*.rst file, that
is referred to from the main index.rst file. You can either place your
contribution within one of the existing subsections, in which case you
simply need to add it to the appropriate index_*.rst file, or you can
create your own subsection. To do this, you will need to create a new
index_*.rst file (likely copying one of the existing ones), and place
the name of your .rst documentation file within it. You will also need
to add the new index_*.rst file to the main index.rst file.

=======================
CHANGING OVERALL OUTPUT
=======================

If you want to adjust overall aspects of the documentation
(e.g. version number, HTML theme, LaTeX options), you can do this by
editing the conf.py file.
