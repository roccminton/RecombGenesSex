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

| variable                      | description                                                                                             |
|-------------------------------|---------------------------------------------------------------------------------------------------------|
| N                             | number of genes                                                                                         |
| dni                           | genome-wide mutation rate                                                                               |
| filename                      | name of the output file                                                                                 |
| abs_path                      | path where the output should be saved                                                                   |
| K = 10_000                    | total population size                                                                                   |
| tend = 100_000                | length of simulation (in generations)                                                                   |
| b = 1.0                       | individual birth rate                                                                                   |
| d = 0.9                       | individual death rate                                                                                   |
| c = nothing                   | competition between individual (c=nothing results in uniform competition)                               |
| rec = 0                       | Recombination rate between 0 (no recombinaton) and 1 (full recombination)                               |
| birthrates = "allbirthrates!" | Choose between "allbirthrates!" (const. population size) and "truerates!" (fluctuation population size) |
| nruns = 3                     | Number of independent runs of the simulation                                                            |

Only N,dni,filename and abs_path are mandatory variables. The others are optional and have a default value as shown.
The output of the simulations is safed in a seperate file for each run and an overview plot over the summery statistics (mutation burden and prevalence) is saved in one pdf file for all runs together.

**Remark:** Depending on the parameter configurations and the length of the simulations, one run can take up to 24 hours. Especially after the increase in mutation burden, the simulations slow down significantly.

## Citation
If you use the code or data from this repository in your research, please cite the paper:
[Insert citation information here]

## Contact
For any questions or inquiries, please contact Luis A. La Rocca at luis.larocca@uni-bonn.de.
