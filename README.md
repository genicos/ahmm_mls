# Ancestry_HMM_MLS

### Overview

Ancestry_HMM_MLS (Ancestry HMM Multi Locus Selection), is a program to infer the strength and presence of multiple loci under selection in admixed populations. Ancestry_HMM_MLS is an extension to Ancestry_HMM, and more information on the format of the input files can be found at https://github.com/russcd/Ancestry_HMM.


#### Basic usage
> ahmm_mls -i [input_file] –s [sample_file] -m [model_file]

For a list of options and arguments, see the built in help message:

> ahmm_mls --help

This will print the following:

> ahmm_mls usage:
> 
>	required:
>		-i [string]		input file
>		-s [string]		sample id and ploidy file
>		-m [string]		model file
>	optional:
>		--help			print this help statement
>		-g			samples are specified with genotypes rather than read counts
>		-c [int]		number of cores used
>		-f			output relative fitnesses rather than selection coefficients
>		-vo			verbose stderr output
>		-R [float]		specify morgan distance from selected site after which they are ignored
>		-t1 [float]		threshold of lnL ratio range in simplex for first stage of optimization
>		-t2 [float]		threshold of lnL ratio range in simplex for second stage of optimization
>		-k [int]		number of skipped regions between adjacent sampled sites
					for each one calculated.

#### Input File and Sample File Format

AHMM_MLS uses the same file format for genotype data as Ancestry_HMM. See https://github.com/russcd/Ancestry_HMM for further details.

#### Model File Format

Each line in the model file corresponds to a model being fit to the genotype data, which may or may not include optimization.

Each line must have the following 5 required columns

1. Site coordinates used in this line (0 for bp, 1 for Morgans)
2. Name of model
3. Local Ancestry decoding method

This is specified as a string containing the characters 'd', 'm', or 's'. If 'd' is present, then the average local ancestry of the samples is computed. If 'm' is present, the expected local ancestry from the model is computed. If 's' is present, then the forward-packward posterior probabilities are decoded for each sample. If neither of these characters are present, no local ancestry decoding takes place (we recommed simply placing an 'n', for no decoding).

4. Admixture fraction (proportion of admixed population from parental population 0)
5. Time since admixture (generations)

After the required columns, you may optionally define selected sites to be included in the model. A Site is defined by its location, its dominance coefficient, and its selection coefficient.
For a fixed site with no optimization, the following 6 columns should be added.

> l \<location of site\> h \<dominance coefficient\> s \<selection coefficient\>

To search for the location of the site, place the location within parenthesis like so:

> l (0.134) h 0.5 s 0.01

To specify a range for a bounded search for the location, add two columns after the location, the first being the lower bound, the second being the upper bound, like so:

> l (0.134) 0.133 0.135 h 0.5 s 0.01

To search for the dominance coefficient of the site, replace the dominance coefficient with a pair of parenthesis, or completely ommit the 'h' and dominance coefficient columns, as in the following two examples:

> l 0.134 h () s 0.01
> l 0.134 s 0.01

To search for the selection coefficient, replace the selection coefficient with a pair of parenthesis, or ommmit the 's' and selection coefficient columns, as in the following two examples:

> l 0.134 h 0.5 s ()
> l 0.134 h 0.5

To restrict the selection coefficient to positive or negative while searching, place a '+' or '-' after the parenthesis, like so:

> l 0.134 h 0.5 s ()+
> l 0.134 h 0.5 s ()-


Here are a few examples of lines that may be in a model file.

> 1 neut_model s 0.3 620

This line would fit a neutral model in with an admixture fraction of 0.3, and a time since admixture of 620 generations. It would then decode the posterior-probabilities for each sample.

> 0 la_compare dm 0.3 620 l 2048541 h 0.5 s 0.01 l 8836891 h 0 s 0.02

This line would fit a model with two selected sites, defined by their bp position. It would then compute the expected local ancestry from the model, and the average local ancestry from the forward backward decoded posterior probabilities accross each sample.

> 1 fit_model n 0.3 620 l (0.134) h 0.5 s ()+ l 0.25

This line would optimize a two site model, with the sites defined in morgans. The first site would not have its position fixed, and would also be optimized for its selection coefficient, while being restricted to only having a positive selection coefficient, and additive selection. The second site has its position fixed, but both the dominance and selection coefficients are optimized. 
