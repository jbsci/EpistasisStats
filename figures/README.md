# Plot generation scripts

Scripts to generate the plots in the manuscript are provided in the `$repositoryroot/figures/scripts` directory. 

## Scripts

1. `basic_epistasis_plot.py`: Used for Figure 1, the scatterplots which show the shape and direction of epistasis present in the dataset(s)
2. `sepdist_compare.py`: Used for Figure 2, generates the comparison of separation distance vs. observed epistasis
3. `loglikealtdiag.py`: Used for Figure 3, generates the diagram to compare the alternative and null models for the likelihood ratio analysis
4. `loglikeplotter.py`: Used for Figure 3, generates the analysis figures for where the experimental data falls on the simulated likelihood ratio distribution
5. `feature_cat_compare.py`: Used for Figure 4, shows the breakdown of categroical features. Also used for the supplemental figures for all cateogrical breakdowns and epistasis presence/breakdown barplots.
6. `rsqpare.py`: Used for supplemental figure for the full validation results

## Flags
Here are the specific flags to generate the figures in the manuscript. Assumes repository pathing.


### Figure 1

To generate Figure 1A: 

`./scripts/basic_epistasis_plot.py -i ../data/processed/skempi_bind_processed.txt -o ./Fig1/Fig1A -ptype bind`

To generate Figure 1B:
`./scripts/basic_epistasis_plot.py -i ../data/processed/protherm4_fold_processed.txt -o ./Fig1/Fig1B -ptype fold`

To combined and generate Figure 1:

`convert -density 600 ./Fig1/Fig1A.pdf ./Fig1/Fig1B.pdf +append +repage -quality 600 ./Fig1/Fig1.pdf`

### Figure 2

To generate Figure 2:

`./scripts/sepdist_compare.py --bind ../data/processed/skempi_bind_processed.txt --fold ../data/processed/protherm4_fold_processed.txt --fs 31 --out ./Fig2/Fig2`

### Figure 3

To generate Figure 3A:

`./scripts/loglikealtdiag.py -o ./Fig3/Fig3A`

To generate Figure 3B:

`./scripts/loglikeplotter.py -sim ../data/loglikedata/bind/sim_1000.txt -exp ../data/loglikedata/bind/exp_1000.txt -ptype bind -o ./Fig3/Fig3B`

To generate Figure 3C:

`./scripts/loglikeplotter.py -sim ../data/loglikedata/fold/sim_1000.txt -exp ../data/loglikedata/fold/exp_1000.txt -ptype fold -o ./Fig3/Fig3C`

To combine to make Figure 3:

`convert -density 600 ./Fig3/Fig3A.pdf ./Fig3/Fig3B.pdf ./Fig3/Fig3C.pdf +append +repage -quality 600 ./Fig3/Fig3.pdf`

### Figure 4

To generate Figure 4A:

`./scripts/feature_cat_compare.py --feature charge_ab2 --system binding --simplename "Charge" --output ./Fig4/Fig4A  --fs 28 --figsize "8x8" --data ../data/processed/skempi_bind_processed.txt --titleside "right" --titleshift 0.035`

To generate Figure 4B:

`./scripts/feature_cat_compare.py --feature cplx_type --system binding --simplename "Complex Type" --output ./Fig4/Fig4B --rot 16 --fs 28 --figsize "8.6x8.82" --data ../data/processed/skempi_bind_processed.txt --titleside "right" --stripnum --titleshift 0.065`

To combine to make Figure 4:

`convert -density 600 ./Fig4/Fig4A.pdf ./Fig4/Fig4B.pdf +repage -quality 600 +append +repage ./Fig4/Fig4.pdf`


### Supplemental Figures


To generate Figure S1A and S1B:

`./scripts/feature_cat_compare.py --data ../data/processed/skempi_bind_processed.txt --presence_bar --output ./Supplement/S1A --system binding`
`./scripts/feature_cat_compare.py --data ../data/processed/protherm4_fold_processed.txt --presence_bar --output ./Supplement/S1B --system folding`

To combined to Figure S1:

`convert -density 600 ./Supplement/S1A.pdf ./Supplement/S1B.pdf +append +repage -quality 600 ./Supplement/S1.pdf`

To generate Figure S2A and S2B:

S2A components:

`./scripts/rsqparse.py -dataroot ../data/validation_results_parsed/binding_runs/ -outdir ./Supplement/ -system bind -fs 24`

S2B components:

`./scripts/rsqparse.py -dataroot ../data/validation_results_parsed/folding_runs/ -outdir ./Supplement/ -system fold -fs 24`

Combining them together:

`convert -density 600 ./Supplement/leave_10per_out_bind.pdf ./Supplement/leave_10per_out_bind_errb.pdf -append +repage -quality 600 ./Supplement/S2A.pdf`
`convert -density 600 ./Supplement/leave_10per_out_fold.pdf ./Supplement/leave_10per_out_fold_errb.pdf -append +repage -quality 600 ./Supplement/S2B.pdf`

`convert -density 600 ./Supplement/S2A.pdf ./Supplement/S2B.pdf +append +repage -quality 600 ./Supplement/S2.pdf`


To generate Figure S3:


S3A

`./scripts/feature_cat_compare.py --feature intside_int --system binding --simplename "Interaction Side" --output ./Supplement/S3A  --fs 28 --figsize "8x8" --data ../data/processed/skempi_bind_processed.txt --titleside "left" --titleshift -0.035`

S3B

`./scripts/feature_cat_compare.py --feature ss_ab3 --system binding --simplename "Secondary Structure" --output ./Supplement/S3B  --fs 28 --figsize "8x8" --data ../data/processed/skempi_bind_processed.txt --titleside "left" --titleshift -0.035`

S3C

`./scripts/feature_cat_compare.py --feature charge_ab2 --system folding --simplename "Charge" --output ./Supplement/S3C  --fs 28 --figsize "8x8" --data ../data/processed/protherm4_fold_processed.txt --titleside "right" --titleshift 0.08`

S3D

`./scripts/feature_cat_compare.py --feature hp_ab3 --system folding --simplename "Hydrophobicity" --output ./Supplement/S3D  --fs 28 --figsize "8x8" --data ../data/processed/protherm4_fold_processed.txt --titleside "right" --titleshift 0.08`

To combine:

S3top:

`convert -density 600 ./Supplement/S3A.pdf ./Supplement/S3B.pdf +append -quality 600 +repage ./Supplement/S3top.pdf`

S3bottom:

`convert -density 600 ./Supplement/S3C.pdf ./Supplement/S3D.pdf +append -quality 600 +repage ./Supplement/S3bottom.pdf`

S3:

`convert -density 600 ./Supplement/S3top.pdf ./Supplement/S3bottom.pdf -append -quality 600 +repage ./Supplement/S3.pdf`


Note: A,B,C,D,...etc designations added via postprocessing with inkscape. 
