# Tropical Elliptic Curves

This repository contains code for dealing with the output of mptopcom [[3]](#3)
in polymake [[2]](#2). Furthermore the extracted core data can be found here,
as well as scripts for converting this into human readable form.

## Generating latex files
For this switch to the root of this repository and run the following command:
```
polymake --script scripts/to_latex.pl
```
(This script has hard coded paths, so it will not work from elsewhere).

Afterwards switch into the latex directory and run the following:
```
cd latex/
pdflatex main.tex
```
Note that this will take a while, since this document contains 4009 TikZ pictures.

## References
<a id="2">[2]</a>
[polymake](https://polymake.org/) - open source software for research in
polyhedral geometry.

<a id="3">[3]</a>
[mptopcom](https://polymake.org/mptopcom) - computing triangulations of point configurations in parallel.

