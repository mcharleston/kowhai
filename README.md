# Kowhai

C++ program to both model coevolution at species level, and (later) do cophylogenetic analysis.

*Kowhai is named after the iconic Aotearoa / New Zealand tree, and rhymes with "co-phy"(logenetics).*

## Introduction 

This C++ program implements my stochastic simulation model for generating species tree - gene tree complexes with *codivergence* (cospeciation), *duplication* including **segmental** duplications, *extinction* and *lineage sorting* which together appear as **loss* events.

It was written and is maintained by Michael A. Charleston.


## Usage

`kowhai` must be compiled using a C++11 or later compliant compiler, such as gcc. 
The compiled binary may be put anywhere in your directory structure, such as in users' home directory or `~/bin/`.

It should be run from command-line / terminal using `kowhai` in such a way that the shell can find the executable (e.g. as `~/bin/kowhai` or `kowhai` if `~/bin` is in your $PATH variable), and arguments as below:

### Input arguments

`kowhai`

*  `-h` or `--help` to print this help message
*  `--sim` [options] 
	* `-nH <int>`
	to set the number of LEAVES/tips in a SINGLE Host/species tree, generated under a Yule model 
*  `-nP <int>` to set the number of Parasite/gene TREES in each replicate
*  `-nR <int>` to set the number of simulation replicates to do (default value 1)
*	`-pC <float>` to set the probability of codivergence at each host node
*	`-pJ <float>` to set the conditional probability of joint duplication
*	`-rB <float>` to set the birth / duplication rate in the dependent phylogenies (default value: )
*	`-rHS <float>` to set the host switch rate in the dependent phylogenies (default value )
*	`-rX <float>` to set the death rate in the dependent phylogenies (default value )
*	`--for-multrec` to provide output suitable for segdup input (default: OFF)
*	`--for-segdup` to provide output suitable for segdup input (default: OFF).
*	`--host-sets-rate` to set the rates of the dependent phylogenies as determined by the HOST lineage (default: OFF)
*	`--verbose` set verbosity to output more stuff.
