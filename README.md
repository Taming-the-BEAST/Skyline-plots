
# Background

Population dynamics influence the shape of the tree and consequently, the shape of the tree contains some information about past population dynamics. The so-called Skyline methods allow to extract this information from phylogenetic trees in a non-parametric manner. It is non-parametric since there is no underlying system of differential equations governing the inference of these dynamics. In this tutorial we will look at two different methods to infer these dynamics from sequence data. The first one is the Bayesian Coalescent Skyline plot {% cite Drummond2005 --file Skylines/master_refs %}, which is based on the coalescent model, and the second one is the Birth-Death skyline {% cite Stadler2013 --file Skylines/master_refs %} plot based on the birth-death model. The conceptual difference between coalescent and birth-death approaches lies in the direction of the flow of time. In the coalescent, the time is modeled to go backwards, from present to past, while in the birth-death approach it is modeled to go forwards. Two other fundamental differences are the parameters that are inferred and the way sampling is treated. 

----

# Programs used in this Exercise


### BEAST2 - Bayesian Evolutionary Analysis Sampling Trees 2

BEAST2 is a free software package for Bayesian evolutionary analysis of molecular sequences using MCMC and strictly oriented toward inference using rooted, time-measured phylogenetic trees {% cite Bouckaert2014 --file Skylines/master_refs %}. This tutorial uses the BEAST2 version 2.4.2.


### BEAUti - Bayesian Evolutionary Analysis Utility

BEAUti is a graphical user interface tool for generating BEAST2 XML configuration files.


### Tracer

[Tracer](http://tree.bio.ed.ac.uk/software/tracer) is used to summarise the posterior estimates of the various parameters sampled by the Markov chain. This program can be used for visual inspection and assessment of convergence. It helps to quickly view median estimates 95% highest posterior density intervals of the parameters, and calculates the effective sample sizes (ESS) of parameters. It also helps to visualise potential parameter correlations.


### R

We will be using [R](\href{https://www.r-project.org) to analyze the output of the birth-death skyline plot. We will start the provided R script from the terminal, hence there is no need for applications like RStudio, which provides a graphical user interface for R. If you prefer using RStudio feel free to do so.

----


# Practical: Bayesian and birth-death skyline plot

In this tutorial we will estimate the dynamics of the Egyptian Hepatitis C epidemic from genetic sequence data collected in 1993.

The aim of this tutorial is to:
- Learn how to infer population dynamics;
- Get to know how to choose the set-up of a skyline analysis;
- Get to know the advantages and disadvantages of the Bayesian Coalescent Skyline and the Birth-Death Skyline.


## The Data
The dataset consists of an alignment of 63 Hepatitis C sequences sampled in 1993 in Egypt {% cite Ray2000 --file Skylines/master_refs %}. This dataset has been used previously to test the performance of skyline methods {% cite Pybus2003,Drummond2005,Stadler2013 --file Skylines/master_refs %}.

With an estimated 15-25%, Egypt has the highest Hepatits C prevalence in the world. In the mid 20^(th) century, the prevalence of Hepatitis C increased drastically (see Figure 1](#fig:prevalence) for estimates). We will try to infer this increase from sequence data. 

\begin{figure}[h!]
\centering
\fbox{\includegraphics[scale=0.5]{figures/Estimated_number_hcv.png}}
\caption{\small The estimated number of Hepatitis C cases in Egypt~\citep{Pybus2003}.}
\label{fig:prevalence}
\end{figure}

<figure>
	<a id="fig:prevalence"></a>
	<img src="figures/Estimated_number_hcv.png" alt="">
	<figcaption>Figure 1: The estimated number of Hepatitis C cases in Egypt {% citep Pybus2003 --file Skylines/master_refs.bib %}.</figcaption>
</figure>


## Creating the Analysis File with BEAUti

We will use BEAUti to generate the configuration file for BEAST2 from the sequence alignment.

### Install BEAST 2 Plug-Ins

While the Bayesian Coalescent Skyline plot is integrated in the core of BEAST2, we need to install the BDSKY package, which contains the birth-death skyline plot functionality. Installation of packages is done using the package manager, which is integrated into BEAUti. Open the package manager with `File > Manage Packages` in BEAUti. Select the BDSKY package and install it using the `Install/Upgrade` button ([Figure 2](fig:install)).

<figure>
	<a id="fig:install"></a>
	<img src="figures/install_bdsky.png" alt="">
	<figcaption>Figure 2: Install the package BDSKY which contains the birth-death skyline functionality.</figcaption>
</figure>


### Setting up the analysis with Bayesian Coalescent Skyline

To import the aligned sequences into BEAUti, use `File > Import Alignment` to select the `*.nexus` file.

<figure>
	<a id="fig:import"></a>
	<img src="figures/import_alignment.png" alt="">
	<figcaption>Figure 3: Import the alignment.</figcaption>
</figure>

BEAUti will recognize the sequences from the `*.nexus` file as nucleotide data. It will do so for sequence files with the character set of ** A | C | G | T | N **, where ** N ** indicates an unknown nucleotide. As soon as other non-gap characters are included (e.g. using **R** or **Y** to indicate purines and pyramidines) BEAUti will not recognize the data as nucleotides anymore, unless the type of data is specified in the `*.nexus` file.

After we have loaded the sequences into BEAUti, we have to specify the evolutionary model. We will be using the very general GTR model ([Figure 4](fig:model)), which estimates transition probabilities between individual nucleotides separately, meaning that transition probabilities between e.g. **A** and **T** will be inferred separately to the ones between **A** and **C**. Additionally, we should allow for rate heterogeneity among sites. We can do this by changing the Gamma Category Count to 4 (normally between 4 and 6).

<figure>
	<a id="fig:model"></a>
	<img src="figures/choose_gtr.png" alt="">
	<figcaption>Figure 4: Set GTR as a site model. Also use a Gamma Category Count of 4.</figcaption>
</figure>

As we use sequences that were sampled at the same point in time, we need to fix the clock rate (for more information on this please refer to the tutorial on molecular clocks). We will use an estimate inferred in {% cite Pybus2001 --file Skylines/master_refs %} to fix the clock rate. In this case all the samples were contemporaneous (at the same time) and the clock rate works as a mapping of the estimated tree branch lengths into calendar time.

We will keep the strict clock model and will set `Clock.rate` to 0.00079.

<figure>
	<a id="fig:clockrate"></a>
	<img src="figures/set_clockrate.png" alt="">
	<figcaption>Figure 5: Set the clock rate to 0.00079.</figcaption>
</figure>

Next, we need to go the `Priors` tab and set the Bayesian Coalescent Skyline as a tree prior (Figure~\ref{fig:coalescent}).

\begin{figure}[h!]
\centering
\fbox{\includegraphics[width=\textwidth]{figures/choose_coalescentSkyline.png}}
\caption{\small Choose the Coalescent Bayesian Skyline as a population prior.}
\label{fig:coalescent}
\end{figure}

%\fixme{It has not been explained previously. Or do you mean that it should have been explained in class? Explain it once again here quickly anyway.}
The Bayesian Coalescent Skyline works by dividing the time between the present and the root of the tree into intervals, thus the number of these intervals has to be defined. Each interval will have a different effective population size. 
The Bayesian Coalescent Skyline will estimate the number of coalescent events within each interval (which is captured in the Group Size parameter) as well as the effective population size for that interval. The number of intervals is equal to the dimension specified. If we have $d$ intervals, the effective population size is allowed to change $d-1$ times. To specify the number of dimensions, we need to first go to the initialization panel. This is by default not visible \texttt{View > Show Initialization Panel}.

%\begin{figure}[h!]
%\centering
%\fbox{\includegraphics[width=0.2\textwidth]{figures/goto_initialization.png}}
%\caption{\small go to the initialization panel.}
%\label{fig:initialization}
%\end{figure}

%\fixme{This procedure needs to be explained better (what do you expand and what do you write where) as we didn't do this before. Also, 5 is the default value, you need to say it.}

For this analysis we will set the number of dimensions to $4$ (the default value is 5). Keep in mind that one has to change the dimension of \textbf{bPopSizes} as well as \textbf{bGroupSizes}. The dimension of both parameters has to be the same (Figure~\ref{fig:dimensions}).

\begin{figure}[h!]
\centering
\fbox{\includegraphics[width=\textwidth]{figures/set_dimension.png}}
\caption{\small Set the dimension of the two parameters, bPopSizes and bGroupSizes, to $4$.}
\label{fig:dimensions}
\end{figure}


Choosing the dimension for the Bayesian Coalescent Skyline can be rather arbitrary. If the dimension is chosen too low, not all population changes are captured, if it is chosen too large, there might be too little information in an interval to support an estimate of a population size. There are implementations in BEAST of the coalescent skyline that either sample dimensions (Extended Bayesian Skyline~\citep{Heled2008}) or do not require dimensions to be specified (Skyride~\citep{Minin2008}).

We can leave the priors as they are and save the settings to *.xml

%The MCMC will now run for a couple of minutes, depending on the speed of the computer. In the meantime, we can open tracer to have a look at the *.log file.








----

This tutorial was written by Author Name for [Taming the BEAST](https://taming-the-beast.github.io) and is licensed under a [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0/). 

----

# Relevant References

{% bibliography --cited --file Skylines/master_refs %}







