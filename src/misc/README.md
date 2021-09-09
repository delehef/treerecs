###### Treerecs – Copyright © by INRIA – All rights reserved – Nicolas Comte, 2017
# Treerecs - GeneTreeEditor
Treerecs - GeneTreeEditor is a small program which provides edition of gene names 
using species names and format conversion (Newick, Nhx or PhyloXML) . 
According to a mapping (using gene names and a separator, a smap or trees only),
GeneTreeEditor can modify gene names, by changing position of the species name, 
change the separator character or simply add a species name. 
GeneTreeEditor can creates an smap of the mapping.

# Usage
## Syntax
*    `./genetreeEditor -h` or `--help`
*    `./genetreeEditor -V` or `--version`
*    `./genetreeEditor -g GENETREE_FILE`[` -s SPECIESTREE_FILE `[`-S SMAP_FILE`]` `]

## Options and details
* `-h`, `--help`                print help, then exit
* `-V`, `--version`             print version number, then exit
* `-v`, `--verbose`             verbose
* `-g`, `--genetree`            gene tree(s) filename in a [Newick](https://en.wikipedia.org/wiki/Newick_format), [NHX](https://sites.google.com/site/cmzmasek/home/software/forester/nhx) or [PhyloXML](http://phyloxml.org/) format
* `-s`, `--speciestree`         species tree filename in a [Newick](https://en.wikipedia.org/wiki/Newick_format), [NHX](https://sites.google.com/site/cmzmasek/home/software/forester/nhx) or [PhyloXML](http://phyloxml.org/) format
* `-S`, `--smap`                map filename. First column as the gene names and the second, species names
* `-N`, `--tree-index`          choose one gene tree in the file defined by its index
* `-o`, `--output`              output filename, default "output.txt"
* `-c`, `--sep`                 separator character in gene names to find the species name (default to '_')
* `-p`, `--prefix`              position of the species_name (before the gene_name: yes (y), after no (n), default yes)
* `-m`, `--create-map`          save map used to reconcile trees in file
* `--new-separator`             change separator between gene name and species name
* `--new-position`              change position of the species name (give prefix or postfix)
* `--case-sensitive`            use a case-sensitive mapping
* `--translate`                 translate tree in a specific format (newick, nhx or phyloxml)

# Treerecs - ALEevaluate_undated
Treerecs - ALEevaluate_undated is an executable of the ALE ([Amalgamated likelihood estimation](https://github.com/ssolo/ALE), by Szollosi GJ et al)
ALEevaluate_undated program. This version has other mapping options and a different implementation in Bio++ use.

# Usage
## Syntax
*   `ALEevaluate_undated species_tree_file gene_tree_file`
*   `ALEevaluate_undated species_tree_file gene_tree_file [separators=gene_name_separator [position=(after, before)]] [smap=leaves_map_file] O_R=OriginationAtRoot delta=DuplicationRate tau=TransferRate lambda=LossRate beta=weight_of_sequence_evidence outputFiles=n`

## Options and details
* `species_tree_file`           species tree filename in a [Newick](https://en.wikipedia.org/wiki/Newick_format) or [NHX](https://sites.google.com/site/cmzmasek/home/software/forester/nhx) format
* `gene_tree_file`              gene tree(s) filename in a [Newick](https://en.wikipedia.org/wiki/Newick_format) or [NHX](https://sites.google.com/site/cmzmasek/home/software/forester/nhx) format
* `smap`                        map filename. First column as the gene names and the second, species names
* `separators`                  gene name separator
* `position`                    species name position "after" or "before" gene name
* `O_R`                         origination at root
* `delta`                       gene duplication rate
* `lambda`                      gene loss rate
* `tau`                         transfer rate
* `beta`                        weight of sequence evidence

More details in ([ALE github page](https://github.com/ssolo/ALE), by Szollosi GJ et al)
