# RecombGenesSex

This repository contains the data and code associated with the paper titled "Recombination enables higher numbers of recessive genes, contributing to the emergence of sexual mating in complex organisms."

## Abstract
The paper explores the drift-barrier hypothesis and the role of recombination in maintaining gene pool stability. It investigates how random genetic drift constrains phenotype refinement under natural selection and highlights the influence of factors like effective population size, mutation rate, and recessive gene count on genomic architecture.

## Reproducing the Findings
You can reproduce the findings presented in the paper using the provided code. The simulations conducted using an adaptive dynamics model with constant overall birth rate (equivalent to Wright-Fisher model) reveal transitions in mutation burden and prevalence, linked to the extinction of wild-type haplotypes. The study confirms the drift-barrier hypothesis and demonstrates the essential role of recombination in preventing population collapse and stabilizing the gene pool in complex organisms.

## Contents
- `code/`: Contains the code for simulations and analysis.
- `data/`: Contains the data to create the figures in the paper.
- `figures/`: Contains the files to generate the figures in the paper.

## Usage
1. Clone this repository.
2. Navigate to the `code/` directory.
3. Include the file `runandsafe_IBM.jl` into a Julia script.
4. Execute the function `dosimulation(N,dni,filename,abs_path)`
5. The function arguments are as follows

| variable   | description                                                                                         
|
|------------|---------------------------------------------------------------------------------------------------------|
| N          | number of genes                                                                                         
|
| dni        | genome-wide mutation rate                                                                               |
| filename   | name of the output file                                                                                 |
| abs_path   | path where the output should be saved                                                                   |
| K          | total population size                                                                                   |
| tend       | length of simulation (in generations)                                                                   |
| b          | individual birth rate                                                                                   |
| d          | individual death rate                                                                                   |
| c          | competition between individual (c=nothing results in uniform competition)                               |
| rec        | Recombination rate between 0 (no recombinaton) and 1 (full recombination)                               |
| birthrates | Choose between "allbirthrates!" (const. population size) and "truerates!" (fluctuation population size) |
| nruns      | Number of independent runs of the simulation                                                            |

## Citation
If you use the code or data from this repository in your research, please cite the paper:
[Insert citation information here]

## Contact
For any questions or inquiries, please contact Luis A. La Rocca at luis.larocca@uni-bonn.de.
