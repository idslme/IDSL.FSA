\name{UFSA_formula_vector_generator}
\alias{UFSA_formula_vector_generator}
\title{
Molecular Formula Vector Generator
}
\description{
This function convert a molecular formulas into a numerical vector
}
\usage{
UFSA_formula_vector_generator(molecular_formula, Elements, LElements = length(Elements),
allowedRedundantElements = FALSE)
}
\arguments{
  \item{molecular_formula}{
molecular formula
}
  \item{Elements}{
a string vector of elements. This value must be driven from the `element_sorter` function.
}
  \item{LElements}{
number of elements. To speed up loop calculations, consider calculating the number of elements outside of the loop.
}
  \item{allowedRedundantElements}{
`TRUE` should be used to deconvolute molecular formulas with redundant elements (e.g. CO2CH3O), and `FALSE` should be used to skip such complex molecular formulas.(default value)
}
}
\value{
a numerical vector for the molecular formula. This function returns a vector of -Inf values when the molecular formula has elements not listed in the `Elements` string vector.
}
\examples{
molecular_formula <- "C12H2Br5Cl3O"
Elements <- UFSA_element_sorter()
mol_vec <- UFSA_formula_vector_generator(molecular_formula, Elements)
##
regenerated_molecular_formula <- UFSA_hill_molecular_formula_printer(Elements, mol_vec)
}
