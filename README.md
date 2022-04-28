# Online FDR in RNASeq Simulation Study Code

## Description

Code for Figures 1 and 2 will be shared upon request since it involves proprietary Merck data. However, we are happy to provide the simulation code in this repo, which is the core of the paper anyway.

* run_change_pi_sim.R simulates the RNAseq data, performs differential expression, and applies FDR control methods
* run_change_pi_sim_reorder.R simulates the RNAseq data, performs differential expression, and applies FDR control methods for the reordered setting
* change_pi_sim_results.Rmd generates Figure 3 and Figure 4 using the data generated from run_change_pi_sim_...

## Getting Started

### Dependencies

* The scripts ("run_change_pi_sim.R" and "run_change_pi_sim_reorder.R") were designed to run on a Sun Grid Engine (SGE) task scheduler
* See the renv.lock file or use `renv::restore()` to use the exact package versions we used in our scripts

### Executing program

* If using SGE, you can the following command the submit the job:

```
qsub -t 1-100 -N change_pi_sim ~/runr.sh run_change_pi_sim.R
```

* -t is the range of simulation runs you want to do (e.g. 1-100 is 100 simulation runs)
* -N is the name of the job
* runr.sh is the shell script that initializes the R script
* Note the Batch procedures are quite computationally heavy

## Authors

Lathan Liou
[@lathanliou](https://twitter.com/lathanliou)
