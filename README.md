# Scripts for ELAI validation and application on VN Robusta

1. evaluate [genetic structure](./genetic_structure) of the African and Vietnamese groups
2. simulate [source populations](./validate_elai/simulate_source)
3. simulate [hybrids](./validate_elai/simulate_hybrids) for testing
4. [run elai](https://github.com/vibaotram/snakelai.git) with different sets of parameters and snps
5. [evaluate](./validate_elai/validate_elai.R) elai results
6. [run elai](https://github.com/vibaotram/snakelai.git) on Vietnamese dataset
7. [obtain final inference](./elai_TR/elai_TR_results.R) of the Vietnamese individuals by merging results of 2 SNP sets
