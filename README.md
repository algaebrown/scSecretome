# scSecretome

Single Cell analysis on the secretome aspect of haematopoietic stem cells
[Presentation Here](https://drive.google.com/file/d/1WaI8buVv79Nqwmrgr63iH0zvlX6cMGJn/view?usp=sharing)

1. `sc_notebook`: exploratory analysis on single cell data on stromal/haematopoietic cells in human/mice bone marrow
    - 5 datasets are available: `pellin` contains human haematopoietic stem cell. `wolock`, `tikhonova`, `barynowano` contain mouse stroma. Each dataset are enriched with different cell type, since different FAC enrichment strategy used.
    - `bbknn` is a method to integrate single cell data of different source.
    - mapping mouse to human homolog are based on [MGI database](http://www.informatics.jax.org/homology.shtml)
2. `secM_notebook`: secretory machinery changes in different cell population
    - annotation of secM machinery are from the Lewis lab, Feizi's list
3. `secP_notebook`: secretory protein/membrane proteins in different cell popultaion
    - annotation of secretory protein relies on [Human Protein Atlas(HPA)](https://www.proteinatlas.org/)
4. `cell2cell_notebook`: cell-cell interaction between cell types
    - based on [cell2cell](https://github.com/earmingol/cell2cell) package
    - reactome enrichment analysis by [ReactomePA](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html) in R
    - GO enrichment in [GOenrich](https://github.com/jdrudolph/goenrich)

scripts are in `scSecretome`
