# Tropical Elliptic Curves

This repository contains the atlas of tropical elliptic curves `atlas.pdf` as
described in our article [??], as well as code and links to datasets we used to
generate the atlas.

## Datasets
In order to arrive at the atlas we first used mptopcom [[2]](#2) to compute at
the regular full triangulations of the Cayley polytope $C(2\Delta_3, 2\Delta_3)$.

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
<a id="1">[1]</a>
[polymake](https://polymake.org/) - open source software for research in
polyhedral geometry.

<a id="2">[2]</a>
[mptopcom](https://polymake.org/mptopcom) - computing triangulations of point configurations in parallel.

<a id="3">[3]</a>
[Regular full triangulations of the Cayley polytope C(2\Delta_3, 2\Delta_3)](https://doi.org/10.5281/zenodo.15682603) Dataset on Zenodo.

<a id="4">[4]</a>
[Regular unimodular triangulations of the Cayley polytope C(2\Delta_3, 2\Delta_3)](https://doi.org/10.5281/zenodo.12820155) Dataset on Zenodo.
