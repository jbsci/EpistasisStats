# Script to generate simulated loglikelihood ratio data for separation distance and epistasis

## General functionality

To produce a dataset similar to that used in this study, use the following flags (example for binding):

`./loglikelihood_separation.py --datfile <(cut -f7,9 -d"	" ../../data/processed/skempi_bind_processed.txt)`

This will write the resulting set of simulated lambdas to output.txt and the experimental likelihood to exp.txt

To generate the resulting figures, please see the figure generation script.
