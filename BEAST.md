# What I've learned about BEAST

First of all, the [BEAST book](https://www.beast2.org/book/) and the [BEAST2 website](https://www.beast2.org/) are very helpful, so I'd defer to those for any questions that may come up. Google Group responses are also often very helpful (found by Googling). And feel free to ask me questions, but I'm by no means an expert!

## Overview
1. [Short background on Bayesian statistics](#short-background-on-bayesian-statistics)
2. [How do you set up a BEAST run?](#how-do-you-set-up-a-beast-run)
3. [How do you know what parameters and priors to use?](#how-do-you-know-what-parameters-and-priors-to-use)
4. [Running BEAST](#running-beast)
5. [Analyzing BEAST output](#analyzing-beast-output)


## Short background on Bayesian statistics

### Frequentist vs. Bayesian statistics

There are two main types of statistics: frequentist and Bayesian. Frequentist statistics uses only the data to generate a point estimate of what you're interested in. For example, maximum likelihood (ML) methods, used by RAxML, are a frequentist approach to statistics-RAxML provides you with a single tree, the most likely tree (maximum likelihood tree). Bayesian statistical methods, on the other hand, use prior knowledge and data to generate a distribution of possible outcomes. BEAST, for instance, provides you with a collection (distribution) of trees.

Takeaway points:
* Frequentist statistics: only data, point estimate
* Bayesian statistics: data and prior knowledge, distribution of estimates

### How does Bayesian statistics take into account prior knowledge?

In Bayesian statistics, each unknown parameter of interest has a prior distribution. In a nutshell, you use the data you have to update the prior distribution, thus generating a posterior distribution that takes into account both the prior distribution and the data.

BEAST uses a method called Markov Chain Monte Carlo (MCMC) to sample from the posterior distribution.


## How do you set up a BEAST run?

The way BEAST is set up, the BEAUTi GUI is used to generate an xml file that can then be run using BEAST. Unfortunately, there is no command line version of BEAUTi, and no easy way to generate BEAST xml files from the command line. However, once the xml files are created, BEAST itself can be run from the command line.

### Pre-BEAUTi

BEAUTi takes a nexus file of the sequences you want to run in BEAST.

The way I set up my sequence names is isolateName_location_date to make it easy to parse them out in BEAUTi. The location is to be able to do ancestral reconstruction and the date is to be able to use temporal information and get a dated tree.

You'll probably have to change the dates to decimal format before loading the nexus file into BEAUTI. Here's an R function I made to do that:

```/nfs/esnitkin/Zena/lib/dateToDecimalYear.R```

If you're adding discrete traits and/or tip dates, here's a good tutorial:

[Ancestral reconstruction/Discrete phylogeography (Version 2.4.1)](https://github.com/BEAST2-Dev/beast-classic/releases/download/v1.3.0/ARv2.4.1.pdf)

### Overview of BEAUTi

In BEAUTi, I set the substitution, clock, and population models, add the dates, add the discrete trait, and set the chain length and amount to log. I change the amount to log in proportion to the amount I change the chain length. It's pretty easy to make a bunch of different ones once you have one done, because you can just change the model you want and then save it as a different name. (Although it'd be nice to make a script to do it at some point!)

I then copy all of the BEAUTi files over to flux to run BEAST.

### What input do you need to generate the BEAST xml file using BEAUTi?

A nexus file of the sequences you want to use. The easiest way to have dated tips or add a discrete trait is to be smart in the way you name the sequences. For instance, you could name them seqName_discreteTrait_dateInDecimalFormat. This makes creating the xml file using BEAUTi *a lot* easier. 

### Should I use an outgroup when using BEAST?

The short answer is probably not. Outgroups often introduce long branches into the tree, which may cause more difficulties in estimating various parameters (see BEAST book p98-99 for more details). 

### What is the maximum number of sequences I can use?

Good question. I'm not entirely sure, but it definitely depends on what types of sequences you have. I'm having trouble using ~400 sequences on a data set with little variation (low ESS - effective sample size, see below). Other people have used this many sequences (or more) with success. 

### If I'm using discrete traits, how many can I use?

Another good question that I don't quite know the answer to. If there are a lot, you'll probably have too many parameters and the BEAST run won't even start. If your discrete trait is location, you may consider combining locations in a similar region together.

## How do you know what parameters and priors to use?

You test different ones out. Start with the simplest (hky nucelotide substitution model, strict molecular clock, coalescent constant population size). If that doesn't work (low ESS), then you have to troubleshoot from there. In theory. I usually start a bunch of different combinations at once so I don't have to wait around after I analyze just one model. It kind of depends on how much computational time you have/want to use. After the overview are some things you can test.

### Overview of choosing the correct models for the data

In a nutshell, I try out different substitution, clock, and population models (but don't change many of the other parameters, at least at first). The important ones to look at are:

- exponential population growth: if zero is not in the 95% HPD, then you can reject a constant population growth model
- lognormal molecular clock: if the distribution of the coefficient of variation for the rate is not near zero, then a strict molecular clock can be rejected (in Tracer, look at the distribution of rate.fileName.coefficientOfVariation)
- could be important: if a lot of the nucleotide substitution rates go to zero in the GTR model, it might be too complex for the data

I also do model comparison by doing path sampling. More info coming soon.

### Nucleotide Substitution Model (GTR vs. HKY)
The prior on the nucleotide substitution rates is Jeffrey's prior (1/x), which is invalid when x &rarr; 0. So if one of the nucleotide substitution rates for the GTR model goes to zero, then you can either:
* Use the HKY model of nucleotide substitution (because the GTR model is probably too parameter-rich for the data)
* Change Jeffrey's prior to a uniform distribution or a lognormal distribution (note: I haven't tried this myself)
When Jeffrey's prior is invalid for one of the substitution rates, the prior and posterior distributions in the trace file have low ESSs. 

### Clock Model (Strict vs. Lognormal)

Look at the coefficient of variation for the relaxed uncorrelated lognormal (UCLN) molecular clock. If it's >> 0.1, then it's not clocklike and a strict molecular clock should not be used. If, on the other hand, the 95% HPD of the coefficient of variation is close to zero, then the substitution rate can be considered clock-like and a strict molecular clock can be used. Also, if the distribution of the coefficient of variation bumps up against zero, then a strict molecular clock cannot be rejected.

When using the UCLN clock model, I usually change the prior on the mean from a Uniform prior to a Gamma(0.001,1000) prior. The Uniform is an improper prior and can cause issues if you want to do Path Sampling analysis for model selection.

### Population Growth Model (Coalescent Constant vs. Exponential vs. Bayesian Skyline)

To decide whether a coalescent constant model can be used, look at the growth rate of the coalescent exponential population growth model in Tracer. If zero is not in the 95% HPD, then a constant population size can be rejected and the coalescent constant population model should not be used. Alternatively, you can run an extended bayesian population growth model and look at the sum of the indicators of all trees. Again, if zero is not in the 95% HPD, then a constant population size can be rejected. I'm still trying to determine the best way to determine whether to use exponential or bayesian skyline. I think the answer is you have to do model comparison. More on that once I've learned more.

## Running BEAST

I've created some scripts to run BEAST on flux once you have the xml files. Currently, they're located in the directory `/nfs/esnitkin/Zena/lib`.

From the directory where I want to generate the BEAST directories (and where the BEAUTi xml files are), I create an input file using (for example) this command:

```/nfs/esnitkin/Zena/lib/generate_input_file_beast_pbs.sh '*.xml' 3``` 

This example would generate an input_beast.txt file with all of the xml files in the directory listed 3 times.

Then, I run this to start the BEAST PBS jobs:
```perl /nfs/esnitkin/Zena/lib/beast_pbs.pl input_beast.txt```

This will create a directory for each run. If you want to see how far a run has gotten, or do preliminary analyses (such as looking at the log file in Tracer), you can find the partial log file (and trees file) in the respective directory. 


## Analyzing BEAST output

After the runs have finished, the first thing I do is copy the log files over to my computer and take a look at them in Tracer. If you load a set of files that have the exact same models (i.e. if you did a run in triplicate), then you can easily look at the individual traces and ESSs as well as the combined traces and ESSs. It's also useful to load different ones to compare probability distributions for the likelihood (more positive is better), etc.

If the trace files look good, then you can combine the logs and/or the trees using this script:

```/nfs/esnitkin/Zena/lib/combine_logs_get_consensus_tree.sh```

You can use `-h` to get more info. Example:

```/nfs/esnitkin/Zena/lib/combine_logs_get_consensus_tree.sh -ltLT '*gtr_logn_exp' -o penn-rm_gtr_logn_exp_100mil```

This would generate a combined log file, a combined trees file, and an mcc tree file (from the location trees). You might have to resample the trees if you run out of memory (use `-r`). (It might be useful to figure out how to use more memory at some point, but I haven't gotten around to that.)


### I have low ESSs. What do I do?

Look [here](https://www.beast2.org/how-do-i-increase-the-ess-of-a-parameter/) and [here](https://www.beast2.org/increasing-esss/). 

If the ESSs for a small number of  parameters are low, then it might not matter if those aren't what you're interested in. If one of those is what you're interested in, then I'd also try Googling about low ESSs for that specific parameter.

