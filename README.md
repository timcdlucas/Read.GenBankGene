# Read.GenBankGene

A hacked function that lets you collect GenBank data while specifying which gene you want. 

+ Give the `gene` argument a string and it takes only the genetic data for that gene.
+ Use `CDS=TRUE` if you only want the coding region of the gene.


+ Definitely not robust and I haven't tested it much. I used it for one analysis.
+ Uses a pretty hacky way to access the genes
+ Can't accept multiple gene names.

The function is hacked from `read.GenBank` [here ](http://svitsrv25.epfl.ch/R-doc/library/ape/html/read.GenBank.html) from the [ape package](http://cran.r-project.org/web/packages/ape/). I've marked in the code which contributions are mine.


Feel free to use my bits of the code. Please check the license for `ape` to use those bits of the code.
