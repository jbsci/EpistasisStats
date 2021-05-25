# Model selection and processing

## Dependencies:

1. R
2. Mass
3. gtools
4. stringi
5. comprehenr
6. rlist
7. StepAICc (provided in stepaicc direcotry, sourced from: https://github.com/biometry/APES/blob/master/Data/Dormann2013/stepAICc.r, unmodified, license: https://github.com/biometry/APES/blob/master/LICENSE.md)

## Usage:

1. Install dependencies 
2. For both binding and folding, set the path to `${reporoot}/data/processed/` (absolute, not relative).
3. Run all, it will take some time but generate two files: best_model_{bind,fold}.txt and delta_rsq_best_{bind,fold}.txt
