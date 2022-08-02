# Stability simulation
Simulation script from [Cadiou et al., 2021](https://pubmed.ncbi.nlm.nih.gov/33652445/).
Instructions on how to use Rscripts are available on the Rscripts themselves.
Seeds need to be specified in order to make results reproducible.

# General information
Scripts were developed to perform a two-layer Monte-Carlo simulation to assess stability and performance in terms of sensitivity and false detection proportion of various variable selection algorithms used for inference purpose in epidemiology.

For each iteration, a realistic exposome is generated using real exposome data. An continuous Gaussian outcome is then generated assuming it is influenced to varying extents by a linear combination of factors drawn from this simulated set.

In [(Cadiou et al., 2021)](https://journals.lww.com/epidem/Abstract/2021/05000/Instability_of_Variable_selection_Algorithms_Used.13.aspx), real data from the Helix project [(Vrihjeid et al., 2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4048258/) were used. These data are not publicly available and thus not provided here. 

For more details, please refer to [(Cadiou et al., 2021)](https://journals.lww.com/epidem/Abstract/2021/05000/Instability_of_Variable_selection_Algorithms_Used.13.aspx).
