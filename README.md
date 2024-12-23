## Larix demographics

### Collaboration with Dan Turck & Dave Tank

<br>

This repository contains my code for fitting demographic models and estimating model parameters for Dan's Larix data.


### Create SFS

FastSimCoal takes the site frequency spectrum (SFS) as input. To generate the SFS, I used Isaac Overcast's `easySFS.py` script from here: [https://github.com/isaacovercast/easySFS](https://github.com/isaacovercast/easySFS ). To determine how many alleles to project down (downsample to) to account for missing data, I ran the first line to preview the sizes of datasets at various projections. The second lines generates the SFS is the projection I used:


```
	easySFS.py -i 4LAR_no_outgroups.recode_dp07mis07.recode.vcf -p LAR_Pop_file_without_outgroups_3_populations.txt --preview

	easySFS.py -i 4LAR_no_outgroups.recode_dp07mis07.recode.vcf -p LAR_Pop_file_without_outgroups_3_populations.txt --proj=8,14,22 --prefix LarK3 -o /Users/harrington/Active_Research/Larix_data_and_outs/SFS
```

The file output from this is `LarK3_MSFS.obs`.

### Divergence time

Divergence time at the root is fixed to either 20k years or 200k years. Assuming a generation time of 180 years, these times become 111 or 1,111 generations. We fit all models under each of these assumed divergence times. 



### Models


1. Make the models. Making the models and ensuring that they are specified correctly is the most tedious part of this. Models are contained in the `Models` directory. Each contains a `.tpl` and `.est` file that specify the coalescent model for FastSimCoal. We use fastsimcoal26 v2.6.0.3.  

To test that each model is correctly specified, I execute a very short test run and visualize the output `.par` file using the `ParFileInterpreter-v6.3.1.r` script. This was previously available on the site hosting FSC2, but has since been replaced by a new version. It is included here in this repo.


* Example of checking a model with a short run

A very short FSC run:
`fsc26 -t *tpl -e *est -n 1000 -L 40 -m -M -u -x  -q`

Check that it looks right - don't worry about parameter estimates, it didn't run long enough to optimize (note this uses a path to the script that exists on my computer):
`Rscript /Applications/fsc26_mac64/ParFileInterpreter-v6.3.1.r *par`


2. Get the linux executable for FSC: http://cmpg.unibe.ch/software/fastsimcoal2/ (for running this on a Linux cluster, as I did).

3. Put the `Models` directory (without any name changes or scripts below won't work) onto the cluster. The `Models` directory must contain a directory for each model to be estimated (each containing a `.tpl` and `.est` file) - the current `Models` directory does this - as well as the sfs file (i.e., a single SFS file in the Models directory outside dirs for each model). 

## EDIT THIS: You will need to add the sfs `.obs` file from Dryad if running this.

4. To run multiple replicates of FSC parameter estimation, use the `Prep_FSC_reps.sh` script. Details of how to use this script are in the comments inside it. Briefly, run it from inside the Models folder made in step 6, supplying an argument specifying if you are using either multi-dimensional SFS (argument MSFS) or two-dimensional SFS (argument jointMAF). This will make 50 replicates in a directory "Reps" in the parent directory of Models, each containing the necessary tpl, est, and obs file. e.g., 

`cd /project/inbreh/turck_fsc/Models`

`Prep_FSC_reps.sh MSFS`

5. Submit all of the jobs to the cluster as a big job array on a SLURM cluster using `FSC_Larix_Mods.slurm` (or similar named slurm scripts for rate/stacks datasets).


6. For each of the models find the single run with the best likelihood. The script `Get_best_FSCacross_mods.sh` wraps the `fsc-selectbestrun.sh` script from here: [https://raw.githubusercontent.com/speciationgenomics/scripts/master/fsc-selectbestrun.sh](https://raw.githubusercontent.com/speciationgenomics/scripts/master/fsc-selectbestrun.sh) (and copied here with some slight modifications) across each model. Comments in this script further describe its usage--run it from within the `Reps` directory. 

```
cd /project/inbreh/turck_fsc/Reps
Get_best_FSCacross_mods.sh
```

To get AIC from the bestlhoods files, go inside `best_L_allMods` directory, created by last script, and run the script `Get_AIC_across_mods.R`(!!!NOTE: This AIC calculation will only work if you output 1 and exactly 1 parameter to this file for each estimated parameter otherwise the AIC values will be incorrect -- i.e., if you have a complex parameter, and you output both the estimated parameter and the complex parameter that is a transformation of the estimated parameter, you'll need a different solution here).

```
cd /project/inbreh/turck_fsc/Reps/best_L_allMods
module load gcc/12.2.0 r/4.2.2
Rscript ../../Get_AIC_across_mods.R
```

7. Convert the parameter estimates into useful units (e.g., years, number of migrants per generation, etc.) -- I still don't have a fully automated solution for this, so has to be done on an ad hoc basis depending on what parameters are included, etc. I've been doing this locally, downloading files in `/project/inbreh/turck_fsc/Reps/best_L_allMods` and for the best fit model determined by AIC, go into that directory and download the files from the `bestrun` directory. `Par_conv_FSC_Larix.R` will do these conversions from these files - see comments in the script for more info. This script will also start some prep for parametric bootstrap estimation of confidence intervals around parameter estimates for the best fit model.

8. Prep parametric bootstrap replicates. `Par_conv_FSC_Larix.R` from previous step generates input for bootstraps. Copy that directory to the cluster. Go into the directory that contains each of the 100 bootstrap reps to be executed and run `Prep_FSC_reps_of_bootreps.sh` to generate 50 fsc replicates of each bootstrap replicate.

```
cd /project/inbreh/turck_fsc/boot_inputs/31k_AncMCoInK3_maxL
../../Prep_FSC_reps_of_bootreps.sh
```

This is because we want to run 50 independent FSC searches for the best likelihood for each of the 100 bootstrap replicates. The replicates are created in a new directory in the parent of the current directory. Use `FSC_Larix_boot.slurm `to run FSC on these reps.

9. Summarize bootstrap replicates. For each bootstrap replicate, need to get the best of the runs, pull the bestlhoods file out, and then summarize these all together. `Get_best_FSCacross_boots.sh` will get the best run within each bootstrap rep - this should be run from within the `Reps` directory for the bootstraps. Run the R script `Get_Larix_pars_across_bootreps.R` on files in the `best_L_allMods` directory created by the previous script.



## NOTE: I think there are other times that aren't being converted - times for when ancient migration ends, possibly others - need to find these 


Population designations:

0 - Cascade
1 - N. Rockies
2 - S Rockies


Hunting down failed runs:
`grep CANCELLED *Larix_7145250*  | cut -d "_" -f 5 | cut -d "." -f 1 | sort > failed2_fsc.txt`

`grep 'Bad parameters'  *.err` This was only one, and also contains "CANCELLED", don't need to get this separately.



- edit sample sizes
- check that all params exist match up in the est file

# just notes to me:

models made (Xs have been checked):
numbers correspond to numbers in Dan's diagram

20k div time:
1. 20_noMK3 - X
2. 20_AncMK3 - X
3. 20_MallAsymK3 - X
4. 20_AncMintK3 - X
5. 20_SecK3 - X
6. 20_SecIntK3 - X
7. 20_AncMCoInK3 - X

200k div time:
1. 200_noMK3 - X
2. 200_AncMK3 - X
3. 200_MallAsymK3 - X
4. 200_AncMintK3 - X
5. 200_SecK3 - X
6. 200_SecIntK3 - X
7. 200_AncMCoInK3 - X
16. 200_secLGMintK3 - X

200k div time with population resize (expecting expansion)added in:
200_AncMCoInExpK3 - 200_AncMCoInK3 with population resize parameter added in - X
200_MallAsymExpK3 - 200_MallAsymK3 with population resize parameter added in - X
200_noMExpK3 - 200_noMK3 with population resize parameter added in - X
200_SecExpK3 - 200_SecK3 with population resize parameter added in - X




## Running analyses using a rate instead of fixed times:

Mutation rate from [Torre et al. 2017](https://academic.oup.com/mbe/article/34/6/1363/3053316)  - 1.19939E–09 per year -> 2.158902e-07 per generation


Rate-based div time with population resize (expecting expansion)added in:
Rate_AncMCoInExpK3 - 200_AncMCoInK3 with population resize parameter added in - X
Rate_MallAsymExpK3 - 200_MallAsymK3 with population resize parameter added in - X
Rate_noMExpK3 - 200_noMK3 with population resize parameter added in - X
Rate_SecExpK3 - 200_SecK3 with population resize parameter added in - X




## May need to update some stuff above here with better description, but I think everything is basically documented - gotta be better at documenting **everything**


## Running FSC with a rate but with output from Stacks

- use `04LAR_no_outgroups.recode_dp07mis07.recode.vcf`  in `~/Active_Research/Larix_data_and_outs/Actually use these/Stacks/Stacks_runs_for_fast_simcoal`


EasySFS, as above:

```
	easySFS.py -i 04StksLAR_no_outgroups.recode_dp07mis07.recode.vcf -p 4Stks_Pop_file_no_outgroups_3pops.txt --preview

	easySFS.py -i 04StksLAR_no_outgroups.recode_dp07mis07.recode.vcf -p 4Stks_Pop_file_no_outgroups_3pops.txt --proj=8,12,19 --prefix LarStkK3 -o /Users/harrington/Active_Research/Larix_data_and_outs/SFS_stacks
```

The file output from this is `LarStkK3_MSFS.obs`.

Then, just follow the "Models" steps above. These are the same models as the rate ones for ipyrad, just with different SFS input, and need to change the number of individuals per pop.


## Running FSC with output from stacks, but calibrated to 31k years at the root

- Use same SFS as for Stacks with a rate
- calibrate root to 172 generations (31k years / 180 years/generation)

Then, just follow the "Models" steps above. These are the same models as the rate ones for ipyrad, just with different SFS input, and need to change the number of individuals per pop.

- These seem best

- Done with running all the models

For bootstraps, 4 replicates returned negative TRESIZE parameters, just exlucded these when summarizing bootstraps



## Testing no monomorphics

The `-0` option does not take monomorphic sites into account when estimating parameters.

Follow the steps above for setting up and running models, except that the FSC command will have `-0` added to the call in the slurm script.

Note that the documentation says that it doens't use a rate in this mode and requires a fixed parameter, so for the rate analysis, this might just ignore the -0. 


This did not seem better, estimates seemed weird, if anything.


<br>
<br>

### stairwayplot2

* Directory `stairwayplot2`

`Larix_stairwayplot2.ipynb` contains Jupyter notebook for estimation of population size through time for populations of *Larix* over multiple different datasets.


### Testing out gadma

For sequence length trying this:

Original seq length from ipyrad stats file was:

```
snps matrix size: (33, 3732122)
sequence matrix size: (33, 110256726)
```


and then vcf filtering log contains:


```
After filtering, kept 35749 out of a possible 3607614 Sites
```

So there's some other loss of sites I don't know about, but the starting numbers for SNP matrix size and possible sites in the vcf filtering log are similar enough. 

So `35749 / 3607614 = 0.0099` so for sequence length, let's multiply by that


`110256726 * 0.0099 = 1091542` for the length



- Initial GADMA run

`gadma_Larix_K3_params.txt` defines parameters, specifies only a single event/set of parameters between the splits and since the most recent split

`gadma_Larix.slurm` runs that paramter set






Plotting of final model failed initially, so had to make an edit in file `/home/sharrin2/.conda/envs/gadma/lib/python3.11/site-packages/moments/ModelPlot.py`

```
ax.grid(b=True, which="major", axis="y", color=gridline_color)
```

to 

```
ax.grid(visible=True, which="major", axis="y", color=gridline_color)
```



and:

```
    fig_kwargs = {
        "figsize": (9.6, 5.4),
        "dpi": 200,
        "facecolor": fig_bg_color,
        "edgecolor": fig_bg_color,
    }

```


to 

```
    fig_kwargs = {
        "dpi": 200,
        "facecolor": fig_bg_color,
        "edgecolor": fig_bg_color,
    }

```


then it works to plot models.




- Initial GADMA run

`gadma_Larix_K3_2event_params.txt` defines parameters, specifies 2 events/sets of parameters between the splits and since the most recent split to allow for things like secondary contact or ancient migration that has ceased in either time period.

`gadma_Larix_2event.slurm` runs that paramter set









