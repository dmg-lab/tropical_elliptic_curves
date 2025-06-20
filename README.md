# Tropical Elliptic Curves

This repository contains the atlas of tropical elliptic curves `atlas.pdf` as
described in our article [??], as well as code and links to datasets we used to
generate the atlas.

## Datasets
In order to arrive at the atlas we first used mptopcom [[2]](#2) to compute at
the regular full triangulations of the Cayley polytope $C(2\Delta_3,
2\Delta_3)$. From this dataset we filtered the (regular) unimodular
triangulations into dataset [[4]](#4). 

On these triangulations we ran the polymake [[1]](#1) script
`cycle_length_iso.pl` to arrive at the isomorphism classes of graphs that
appear among these tropical elliptic curves. The triangulations are numbered
and we refer to these numbers as IDs. The resulting data is contained four
polymake data files in the folder `data` of this repository:
- `regular_unimodular.representative_graphs`: For every ID of a representative
  of a class of tropical elliptic curves up to graph isomorphism, map this ID
  to the corresponding graph. The polymake datatype is `Map<Int,
  Set<Array<Int>>>`.
- `regular_unimodular.representative_triangulations`: For an ID as above, map
  it to the triangulation realizing this graph. The polymake datatype is
  `Map<Int, Array<Set<Int>>>`.
- `regular_unimodular.cl2ids`: The "cl" stands for "cycle length". Map every
  cycle length to the set of all IDs realizing a tropical elliptic curve with
  that cycle length. The polymake datatype is `Map<Int, Set<Int>>`.
- `regular_unimodular.classes_histogram`: For every ID as in the first two
  entries, map this ID to the number of regular unimodular triangulations with
  the same graph. The polymake datatype is `Map<Int, Int>`

These can then be turned into a large latex document as described below.

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
