---
author: Nicola F. Müller,Louis du Plessis
level: Beginner
title: Skyline plots
subtitle: Inference of past population dynamics using Bayesian Coalescent Skyline and Birth-Death Skyline plots.
beastversion: 2.5.x
tracerversion: 1.7.x
---


# Background

Population dynamics influence the shape of the tree and consequently, the shape of the tree contains some information about past population dynamics. The so-called Skyline methods allow to extract this information from phylogenetic trees in a non-parametric manner. It is non-parametric since there is no underlying system of differential equations governing the inference of these dynamics. In this tutorial we will look at two different methods to infer these dynamics from sequence data. The first one is the Coalescent Bayesian Skyline plot {% cite Drummond2005 --file Skyline-plots/master-refs %}, which is based on the coalescent model, and the second one is the Birth-Death skyline {% cite Stadler2013 --file Skyline-plots/master-refs %} plot based on the birth-death model. The conceptual difference between coalescent and birth-death approaches lies in the direction of the flow of time. In the coalescent, the time is modeled to go backwards, from present to past, while in the birth-death approach it is modeled to go forwards. Two other fundamental differences are the parameters that are inferred and the way sampling is treated. 

----

# Programs used in this Exercise


### BEAST2 - Bayesian Evolutionary Analysis Sampling Trees 2

BEAST2 ([http://www.beast2.org](http://www.beast2.org)) is a free software package for Bayesian evolutionary analysis of molecular sequences using MCMC and strictly oriented toward inference using rooted, time-measured phylogenetic trees. This tutorial is written for BEAST v{{ page.beastversion }} {% cite Bouckaert2014 --file Skyline-plots/master-refs %}.


### BEAUti - Bayesian Evolutionary Analysis Utility

BEAUti2 is a graphical user interface tool for generating BEAST2 XML configuration files.

Both BEAST2 and BEAUti2 are Java programs, which means that the exact same code runs on all platforms. For us it simply means that the interface will be the same on all platforms. The screenshots used in this tutorial are taken on a Mac OS X computer; however, both programs will have the same layout and functionality on both Windows and Linux. BEAUti2 is provided as a part of the BEAST2 package so you do not need to install it separately.


### Tracer

Tracer ([http://beast.community/tracer](http://beast.community/tracer)) is used to summarise the posterior estimates of the various parameters sampled by the Markov Chain. This program can be used for visual inspection and to assess convergence. It helps to quickly view median estimates and 95% highest posterior density intervals of the parameters, and calculates the effective sample sizes (ESS) of parameters. It can also be used to investigate potential parameter correlations. We will be using Tracer v{{ page.tracerversion }}.


### R

We will be using [R](\href{https://www.r-project.org) to analyze the output of the birth-death skyline plot. We will start the provided R-script from the terminal, hence there is no need for applications like RStudio, which provides a graphical user interface for R. If you prefer using RStudio feel free to do so.

----


# Practical: Bayesian and Birth-death Skyline Plots

In this tutorial we will estimate the dynamics of the Egyptian Hepatitis C epidemic from genetic sequence data collected in 1993.

The aim of this tutorial is to:

- Learn how to infer population dynamics;
- Get to know how to choose the set-up of a skyline analysis;
- Get to know the advantages and disadvantages of the Coalescent Bayesian Skyline Plot and the Birth-Death Skyline.


## The Data
The dataset consists of an alignment of 63 Hepatitis C sequences sampled in 1993 in Egypt {% cite Ray2000 --file Skyline-plots/master-refs %}. This dataset has been used previously to test the performance of skyline methods {% cite Drummond2005 --file Skyline-plots/master-refs %} {% cite Stadler2013 --file Skyline-plots/master-refs %}.

With an estimated 15-25%, Egypt has the highest Hepatits C prevalence in the world. In the mid 20^(th) century, the prevalence of Hepatitis C increased drastically (see [Figure 1](#fig:prevalence) for estimates). We will try to infer this increase from sequence data. 

The alignment file can be downloaded from the Taming the BEAST website at [https://taming-the-beast.org/tutorials/Skyline-plots/](https://taming-the-beast.org/tutorials/Skyline-plots/) by downloading the file `hcv.nexus` from the left-hand panel, under the heading **Date**.

<figure>
	<a id="fig:prevalence"></a>
	<img style="width:50%;" src="figures/Estimated_number_hcv.png" alt="">
	<figcaption>Figure 1: The growth of the effective population size of the Hepatitis C epidemic in Egypt {% cite Pybus2003 --file Skyline-plots/master-refs.bib %}.</figcaption>
</figure>
<br>


## Creating the Analysis Files with BEAUti

We will use BEAUti to generate the configuration file for BEAST2 from the sequence alignment.

### Install BEAST 2 packages

While the coalescent-based Bayesian Skyline Plot is integrated in the BEAST2 core, we need to install the BDSKY package, which contains the birth-death skyline model. Installation of packages is done using the package manager, which is integrated into BEAUti. 


> Open the **BEAST2 Package Manager** by navigating to **File > Manage Packages**. 
> 
> Install the **BDSKY** package by selecting it and clicking the **Install/Upgrade** button ([Figure 2](#fig:install)).
>

After the installation of a package, the program is on your computer, but BEAUti is unable to load the template files for the newly installed model unless it is restarted. So, let's restart BEAUti to make sure we have the **BDSKY** model at hand.


> Close the **BEAST2 Package Manager** and **_restart_** BEAUti to fully load the **BDSKY** package.
> 

<figure>
	<a id="fig:install"></a>
	<img style="width:75%;" src="figures/install_bdsky.png" alt="">
	<figcaption>Figure 2: Install the BDSKY package which contains the birth-death skyline model.</figcaption>
</figure>
<br>

----

### Setting up the Coalescent Bayesian Skyline analysis

To start we have to import the alignment into BEAUti.

> In the **Partitions** panel, import the nexus file with the alignment by navigating to **File > Import Alignment** in the menu and then finding the `hcv.nexus` file on your computer **or** simply drag and drop the file into the **BEAUti** window.
> 

BEAUti will recognize the sequences from the `*.nexus` file as nucleotide data. It will do so for sequence files with the character set of **A | C | G | T | N**, where **N** indicates an unknown nucleotide. As soon as other non-gap characters are included (e.g. using **R** or **Y** to indicate purines and pyramidines) BEAUti will not recognize the data as nucleotides anymore (unless the type of data is specified in the `*.nexus` file) and open a dialogue box to confirm the data type.

The sequences were all sampled in 1993 so we are dealing with a homochronous alignment and do not need to specify tip dates. 

> Skip the **Tip Dates** panel and navigate to the **Site Model** panel.
>

The next step is to specify the model of nucleotide evolution (the site model). We will be using the GTR model, which is the most general reversible model and estimates transition probabilities between individual nucleotides separately. That means that the transition probabilities between e.g. **A** and **T** will be inferred separately to the ones between **A** and **C**, however transition probabilities from **A** to **C** will be the same as **C** to **A** etc. Additionally, we allow for rate heterogeneity among sites. We do this by changing the Gamma Category Count to 4 (normally between 4 and 6).

> Change the **Gamma Category Count** to 4, make sure that the estimate box next to the **Shape** parameter of the Gamma distribution is ticked and set **Subst Model** to **GTR**. Make sure that the estimate box is ticked for all but one of the 6 rates (there should be 5 ticked boxes) and that **Frequencies** are estimated ([Figure 3](#fig:model)).
>

<figure>
	<a id="fig:model"></a>
	<img src="figures/choose_gtr.png" alt="">
	<figcaption>Figure 3: Set a GTR site model with a Gamma Category Count of 4.</figcaption>
</figure>
<br>

> **Topic for discussion:** Why are only 5 of the 6 rates of the GTR model estimated? 
>

Because our sequences are homochronous there is no information in our dataset to estimate the clock rate (for more information on this refer to the [prior-selection](../Prior-selection/) tutorial) and we have to use external information to calibrate the clock. We will use an estimate inferred in {% cite Pybus2001 --file Skyline-plots/master-refs %} to fix the clock rate. In this case all the samples were contemporaneous (at the same time) and the clock rate is simply a scaling of the estimated tree branch lengths (in substitutions/site) into calendar time.

> Navigate to the **Clock Model** panel.
> 
> Leave the clock model as a **Strict Clock** and set **Clock.rate** to 0.00079 s/s/y ([Figure 4](#fig:clockmodel)). 
> (Note that BEAUti is smart enough to know that the clock rate cannot be estimated on this dataset and grays out the estimate checkbox).
>

<figure>
	<a id="fig:clockmodel"></a>
	<img src="figures/set_clockrate.png" alt="">
	<figcaption>Figure 4: Set the clock rate to 0.00079 s/s/y.</figcaption>
</figure>
<br>

Now we are ready to set up the Bayesian Skyline Plot as a tree-prior. 

> Navigate to the **Priors** panel and select **Coalescent Bayesian Skyline** as the tree prior ([Figure 5](#fig:coalescent)).
>

<figure>
	<a id="fig:coalescent"></a>
	<img src="figures/choose_bsp.png" alt="">
	<figcaption>Figure 5: Choose the Coalescent Bayesian Skyline as a tree prior.</figcaption>
</figure>
<br>

The Coalescent Bayesian Skyline divides the time between the present and the root of the tree (the tMRCA) into segments, and estimates a different effective population size ({% eqinline N_e %}) for each segment. The endpoints of segments are tied to the branching times (also called coalescent events) in the tree ([Figure 6](#fig:coal_events)), and the size of segments is measured in the number of coalescent events in each segment. The Coalescent Bayesian Skyline groups coalescent events into segments and jointly estimates the {% eqinline N_e %} (**bPopSizes** parameter in BEAST) and the size (**bGroupSizes** parameter in BEAST) of each segment. To set the number of segments we have to change the dimension of **bPopSizes** and **bGroupSizes** (note that the dimension of both parameters has to be the same).


<figure>
	<a id="fig:coal_events"></a>
	<img style="width:50%;" src="figures/coalescent_intervals.png" alt="">
	<figcaption>Figure 6: Example tree where the red dotted lines show the time-points of coalescent events.</figcaption>
</figure>
<br>


> To change the number of segments we have to navigate to the **Initialialization** panel, which is by default not visible. Navigate to **View > Show Initialization Panel** to make it visible and navigate to it ([Figure 7](#fig:init)).
>
> Set the dimension of **bPopSizes** and **bGroupSizes** to 4 (the default value is 5) ([Figure 8](#fig:dimensions)). 
>

<figure>
	<a id="fig:init"></a>
	<img style="width:25%;" src="figures/goto_initialization.png" alt="">
	<figcaption>Figure 7: Show the initialization panel.</figcaption>
</figure>
<br>

<figure>
	<a id="fig:dimensions"></a>
	<img src="figures/set_dimension.png" alt="">
	<figcaption>Figure 8: Set the dimension of bPopSizes and bGroupSizes to 4.</figcaption>
</figure>
<br>

This sets the number of segments equal to 4 (the parameter dimension), which means {% eqinline N_e %} will be allowed to change 3 times between the tMRCA and the present (if we have {% eqinline d %}  segments, {% eqinline N_e %} is allowed to change {% eqinline d-1 %} times. 

We can leave the rest of the priors as they are and save the XML file. We want to shorten the chain length and decrease the sampling frequency so the analysis completes in a reasonable time and the output files stay small. (Keep in mind that it will be necessary to run a longer chain for parameters to mix properly).

> Navigate to the **MCMC** panel. 
>
> Change the **Chain Length** from 10'000'000 to 3'000'000.
> 
> Click on the arrow next to the **tracelog** and change the **File Name** to `$(filebase).log` and set the 
> **Log Every** to 3'000.
> 
> Click on the arrow next to the **treelog** and change the **File Name** to `$(filebase).$(tree).log` and set the 
> **Log Every** to 3'000.
> 
> Leave all other settings at their default values and save the file as `hcv_coal.xml`.
>

When we run the analysis `$(filebase)` in the name of the `*.log` and `*.trees` files will be replaced by the name of the XML file. This is a good idea, since it makes it easy to keep track of which XML files produced which output files.

Now we are ready to run the analysis.

> Start **BEAST2** and choose the file `hcv_coal.xml`. 
> 
> If you have **BEAGLE** installed tick the box to **Use BEAGLE library if available**, which will make the run faster.
>
> Hit **Run** to start the analysis.
>

The analysis will take 5-10 minutes to run. Read through the next section while waiting for your results or start preparing the XML file for the [birth-death skyline](#sec:bdsky) analysis.



### Effective Population Size


---

This makes the Bayesian Skyline flexible enough to model very complicated {% eqinline N_e %} dynamics, provided that enough segments are specified. 

What the Bayesian Coalescent Skyline is doing is dividing the full tree into dimension {% eqinline d %} segments between coalescent events. It then estimates based on what happens in those segments, what the effective population size is 

(the actual estimation additionally involves smoothing priors which restrict the difference in effective population sizes between two segments).

The length of an segment is therefore not a fixed value but dependent on where the coalescent events are ([Figure 6](#fig:coal_principle), compare later to [Figure 18](#fig:bdsky_principle) for the birth-death skyline), as well as the number of events contained within an segment (the `GroupSize` parameter). 

segments are grouped together, because as [Figure 6](#fig:coal_principle) shows, some of the segments can be very small, which may lead to erratic estimates. Grouping segments together leads to smoother estimates.

---

The effective population size ({% eqinline N_e %}) is the inverse of the rate of coalescence {% eqinline \lambda %}. The larger {% eqinline N_e %} is the less likely lineages are to coalesce. Thus, intervals in a sampled tree with many branching events often coincides with periods when the population size was small. Similarly, few branching events occur during periods of large population size. (Note that these results are conditioned on building a tree by sampling a small fraction of the population). 

{% eq 

\lambda = \frac{1}{N_e}

%}

For an SIR model (**S**usceptible, **I**nfected and **R**ecovered), {% eqinline N_e %} is proportional to the overall population size {% eqinline N %} and the number of infected {% eqinline I %} and inversely proportional to the transmission rate {% eqinline \theta %}. 

{% eq 
	N_e = \frac{I}{\theta} \frac{N}{S}
%}


Estimates of {% eqinline N_e %} therefore do not directly tell us something about the number of infected, nor the transmission rate. However, changes in {% eqinline N_e %} can be informative about changes in the transmission rate or the number of infected (if they do not cancel out).




### Exploring the results of the Coalescent Bayesian Skyline analysis

For the reconstruction of the population dynamics, we need two files, the `*.log` file and the `*.trees` file. The log file contains the information about the group sizes and population sizes of each segment, while the trees file is needed for the times of the coalescent events. 

> Load the logfile into **Tracer** to check mixing and parameter estimates ([Figure 9](#fig:tracer_bsp)).  
>

<figure>
	<a id="fig:tracer_bsp"></a>
	<img src="figures/bsp_tracer.png" alt="">
	<figcaption>Figure 9: Loading the log file into Tracer.</figcaption>
</figure>
<br>

Because we shortened the chain most parameters have very low ESS values. If you like, you can compare your results with the example results we obtained with identical settings and a chain of 30,000,000 (`hcv_coal_30M.log`).  

> Navigate to **Analysis > Bayesian Skyline Reconstruction**. From there open the `*.trees` file. To get the correct dates in the analysis we should specify the **Age of the youngest tip**. In our case it is 1993, the year where all the samples were taken. If the samples were taken through time (heterochronous data), the age of the youngest tip is the time when the most recent sample was taken. 
> 
> Press **OK** to reconstruct the past population dynamics ([Figure 10](#fig:trees)).
>

<figure>
	<a id="fig:trees"></a>
	<img style="width:75%;" src="figures/open_trees.png" alt="">
	<figcaption>Figure 10: Reconstructing the Bayesian Skyline plot in Tracer.</figcaption>
</figure>
<br>

The output will have the years on the x-axis and the effective population size on the y-axis. By default, the y-axis is on a log-scale. If everything worked as it is supposed to work you will see a sharp increase in the effective population size in the mid 20^(th) century, similar to what is seen on [Figure 11](#fig:skyline).

(Note that the reconstruction will only work if the `*.log` and `*.trees` files contain the same number of states and both files were logged at the same frequency).

<figure>
	<a id="fig:skyline"></a>
	<img style="width:75%;" src="figures/skyline_analysis.png" alt="">
	<figcaption>Figure 11: Coalescent Bayesian Skyline analysis output. The black line is the median estimate of the estimated effective population size (can be changed to the mean estimate). The two blue lines are the upper and lower bounds of the 95% HPD interval. The x-axis is the time in years and the y-axis is on a log-scale.</figcaption>
</figure>
<br>

There are two ways to save the analysis, it can either be saved as a `*.pdf` for display purposes or as a tab delimited file. 

> Navigate to **File > Export Data Table**. 
> 
> Enter the filename as `hcv_coal.tsv` and save the file.
>

The exported file will have five rows, the time, the mean, median, lower and upper boundary of the 95% HPD interval of the estimates, which you can use to plot the data with other software (R, Matlab, etc).




### Choosing the Dimension

Choosing the dimension for the Bayesian Skyline can be rather arbitrary. If the dimension is chosen too low, not all population size changes are captured (the resulting skyline plot is too smooth), but if it is chosen too large, there may be too little information in an segment to support a robust estimate of the population size (the resulting skyline plot is too noisy). There are extensions in BEAST of the coalescent skyline that either sample the number of segments (Extended Bayesian Skyline {% cite Heled2008 --file Skyline-plots/master-refs %}) or do not require the number of segments to be specified (Skyride {% cite Minin2008 --file Skyline-plots/master-refs %}).



If we compare the estimates of the population dynamics using different dimensions, we see that most of the dynamics are already captured with having only 2 dimensions, as shown in [Figure 9](#fig:comparison). Adding more dimensions only changes the inferred effective population size before 1900. Note that adding more dimensions adds a slight dip before the increase in the effective population size (around 1900). When comparing to the HPD intervals ([Figure 8](#fig:skyline)) we see that this dip is not significant and may not be indicative of a real decrease in the effective population size before the subsequent increase.


<figure>
	<a id="fig:comparison"></a>
	<img style="width:100%;" src="figures/comparison_dimension.png" alt="">
	<figcaption>Figure 9: Estimated mean effective population sizes using different dimensions.</figcaption>
</figure>
<br>

The choice of the number of dimensions can also have a direct effect on how fast the MCMC converges ([Figure 10](#fig:ess)). The slower convergence with increasing dimension can be caused by e.g. less information in intervals. To some extent it is simply caused by the need to estimate more parameters though.

<figure>
	<a id="fig:ess"></a>
	<img style="width:50%;" src="figures/ess_vs_dim_coal.png" alt="">
	<figcaption>Figure 10: The ESS value of the posterior after running an MCMC chain with {% eqinline 10^7 %} samples, logged every {% eqinline 10^3 %} steps and a burnin of 10% for using different dimensions of the Coalescent Bayesian Skyline.</figcaption>
</figure>
<br>

----

<a id="sec:bdsky"></a>

### Setting up the Birth-death Skyline analysis 


In the first analysis, we used the coalescent approach to estimate population dynamics. We now want to repeat the analysis using the birth-death skyline model. We will use the same model setup as in the previous analysis and only change the tree prior. 

> Restart **BEAUti**, load `hcv.nexus` as before and set up the same site and clock model as in the Coalescent Bayesian Skyline analysis.
>

We will need to set the prior to **Birth Death Skyline Contemporary**, since the sequences were all sampled at the same point in time. For samples that were taken through time, we would use **Birth Death Skyline Serial**. As with the Coalescent Bayesian Skyline, we need to pick the number of dimensions. Here we set the dimensions for {% eqinline R_e %}, the effective reproduction number, which denotes the average number of secondary infections caused by an infected person at a given time during the epidemic, i.e. an {% eqinline R_e %} of 2 would mean that every infected person causes two new infections on average. In other words, an {% eqinline R_e %} above 1 means that the number of cases are increasing, therefore the disease will cause an exponentially growing epidemic, and an {% eqinline R_e %} below 1 means that the epidemic will die out. 


> Navigate to the **Priors** panel and select **Birth Death Skyline Contemporary** as the tree prior ([Figure 11](#fig:bdsky)).
>
> Then, click on the button where it says **initial = [2.0] [0.0, Infinity]** next to **reproductiveNumber**. A pop-up window will open which allows us to change the dimension of the parameter ([Figure 9](#fig:dimensions_bdsky)). In this case we will keep the default dimension of 10. 
> 
> Press **OK** to close the pop-up window.
>

<figure>
	<a id="fig:bdsky"></a>
	<img src="figures/choose_bdsky.png" alt="">
	<figcaption>Figure 11: Setting the prior on the tree to the birth-death skyline.</figcaption>
</figure>
<br>

<figure>
	<a id="fig:dimensions_bdsky"></a>
	<img style="width:75%;" src="figures/choose_dimension_bdsky.png" alt="">
	<figcaption>Figure 12: Setting the dimension of the reproductiveNumber parameter.</figcaption>
</figure>
<br>

As with setting the dimension of the Coalescent Bayesian Skyline the dimension of {% eqinline R_e %} can also be set in the initialization panel. Choosing this dimension can again be arbitrary and may require the testing of a few different values. Too few intervals and not all rate shifts are captured. Too many intervals and the intervals may not contain enough information to infer parameters.

The **Birth Death Skyline Contemporary** model has 4 parameters, **reproductiveNumber** ({% eqinline R_e %}), **becomeUninfectiousRate** (the rate at which infected patients become uninfectious, {% eqinline \delta %}, through recovery, death or isolation), **rho** (the proportion of lineages sampled in the present, {% eqinline \rho %}) and **origin** (the time at which the index case became infected, which is always earlier than the tMRCA of the tree). We may know some of these parameters from literature or be able to estimate them from external sources. For example, the average time that patients are able to transmit a disease is informative about the **becomeUninfectiousRate**. This prior knowledge we can incorporate in our analysis by setting appropriate priors for these parameters.

 

We will use a lognormal prior for {% eqinline R_e %}. This is a good prior distribution to use for rates since it is always positive (a rate cannot be negative) and has a long tail defined over all positive numbers. The long tail allows arbitrarily high estimates of {% eqinline R_e %}, but does not place much weight on very high rates. This agrees with our prior knowledge about the {% eqinline R_e %} (most diseases have an {% eqinline R_e %} between 1.2 and 5. Measles is one of the most infectious diseases we know about and has {% eqinline R_e \approx 18 %}). If an epidemic is neither growing or declining, it has an {% eqinline R_e %} of 1, which we will use as a null hypothesis, by setting a prior on {% eqinline R_e %} centered around 1 (we assume that if there isn't a strong signal in an interval for an epidemic to grow or decline that {% eqinline R_e = 1 %}, i.e. the epidemic stays constant). 

> Select a **Log Normal** distribution for the **reproductiveNumber** prior.
>
> Click on the arrow to the left of **reproductiveNumber** to open all the options for {% eqinline R_e %} settings 
>
> Set **M** to 0, which results in a median of 1. We set the variance to 1.25, which places most weight below 7.82 (95% quantile). ([Figure 14](#fig:r0prior)).
>

<figure>
	<a id="fig:r0prior"></a>
	<img src="figures/bdsky_prior_r0.png" alt="">
	<figcaption>Figure 14: Setting the {% eqinline R_e %} prior.</figcaption>
</figure>
<br>

For the becoming uninfectious rate we will again use a log normal prior. The inverse of the becoming noninfectious rate is the average infectious period. In some patients an HCV infection only lasts a few weeks, while in others it is a chronic infection lasting for many years. Setting {% eqinline M=0 %} and {% eqinline S=1.25 %} results in the same prior we used for the {% eqinline R_e %}.  In terms of the becoming noninfectious rate, this translates to the 95% quantiles for the infectious period falling between 0.128 years (46.67 days) and 11.59 years, with a median of 1 year. We will see later that there is a strong signal in the data for a longer becoming noninfectious period. 

> Set the same prior for **becomeUninfectiousRate** as for **reproductiveNumber** (Log Normal, with M=0.0, S=1.25) ([Figure 15](#fig:bURprior))
>

<figure>
	<a id="fig:bURprior"></a>
	<img src="figures/bdsky_prior_uninf.png" alt="">
	<figcaption>Figure 15: Setting the becoming uninfectious rate prior.</figcaption>
</figure>
<br>

The sampling proportion, {% eqinline \rho %}, represents the proportion of HCV cases in Egypt in 1993 that are included in the analysis. In 1993 Egypt had a population of roughly 60 million people, and with a prevalence of at least 15% this translates into millions of cases, while we only have 63 sequences. 

We will use a beta distribution for the prior on {% eqinline \rho %}. Beta distributions are a very flexible class of distributions that are only defined between 0 and 1, making them ideal to use for proportions. 

> Select a **Beta** distribution for the **rho** prior.
>
> Click on the arrow to the left of **rho** to open all the options for the prior settings.
> 
> Alpha to 1 and Beta to 9999, reflecting our prior knowledge that our dataset represents only a miniscule fraction of cases ([Figure 17](#fig:rhoprior)).
>

<figure>
	<a id="fig:rhoprior"></a>
	<img src="figures/bdsky_prior_rho.png" alt="">
	<figcaption>Figure 17: Setting the prior on {% eqinline \rho %}.</figcaption>
</figure>
<br>

Finally, we need to set a prior for the origin of the epidemic. We will once again use a log normal distribution for this parameter.  Note that the origin also has to be positive and needs to be bigger than the MRCA of the tree. We know that HCV has been circulating in Egypt for at least a hundred years, so we set a prior with a median value greater than 100. 

> Set a **Log Normal** prior for **origin** with **M = 5** and **S = 0.5** ([Figure 16](#fig:oriprior)), resulting in a median prior estimate for the origin of 148 years.
>

<figure>
	<a id="fig:oriprior"></a>
	<img src="figures/bdsky_prior_ori.png" alt="">
	<figcaption>Figure 16: Setting the prior on the origin of the epidemic.</figcaption>
</figure>
<br>

The rest of the priors pertain to the site model parameters and we can leave them as they are. 

> Navigate to the **MCMC** panel. 
>
> Change the **Chain Length** from 10'000'000 to 3'000'000.
> 
> Click on the arrow next to the **tracelog** and change the **File Name** to `$(filebase).log` and set the 
> **Log Every** to 3'000.
> 
> Click on the arrow next to the **treelog** and change the **File Name** to `$(filebase).$(tree).log` and set the 
> **Log Every** to 3'000.
> 
> Leave all other settings at their default values and save the file as `hcv_bdsky.xml`.
>

Now we are ready to run the analysis.

> Start **BEAST2** and choose the file `hcv_bdsky.xml`. 
> 
> If you have **BEAGLE** installed tick the box to **Use BEAGLE library if available**, which will make the run faster.
>
> Hit **Run** to start the analysis.
>




### The parameterization of the Birth-Death Model

The Birth-Death model is parameterized very differently from the coalescent model. While the coalescent uses the effective population size, which as the name already tells us is defined on a population level, the birth-death model uses per lineage rates. The transmission rate {% eqinline \lambda %} tells us at which rate infected individuals infect susceptibles. This rate is also referred to as the birth rate. The sampling rate {% eqinline \psi %} and the sampling probability {% eqinline \rho %} describe how likely it is for an infected individual to be sampled and therefore how likely they are to appear in the tree as tips. The becoming noninfectious rate {% eqinline \delta %} is the sum of the death rate {% eqinline \mu %} and the sampling rate {% eqinline \psi %}.

{% eq
\delta = \psi + \mu
%}

The death rate {% eqinline \mu %} is the rate at which lineages disappear (go extinct) from a population without being sampled. You can also see from the above equation that we assume that a sampled lineage cannot transmit anymore. The consequence for the phylogeny is that a sampled lineage cannot be an ancestor of any lineage. This assumption can be relaxed, but we will not do so during this tutorial.

The {% eqinline R_e %} we estimate is then defined as follows:

{% eq
R_{0} = \frac{\lambda}{\psi + \mu} = \frac{\lambda}{\delta}
%}

    |                                                          |         
--------------------------------------------------------------:|:----------------------
if {% eqinline \lambda > \delta %} then {% eqinline R_e > 1 %} | epidemic grows
if {% eqinline \lambda = \delta %} then {% eqinline R_e = 1 %} | epidemic stays constant
if {% eqinline \lambda < \delta %} then {% eqinline R_e < 1 %} | epidemic declines

The birth-death skyline allows these rates to change over time. This is done by dividing the time from the origin (of the epidemic, which is not necessarily the same as the root of the tree) to the most recent sample into dimension {% eqinline d %} equally spaced intervals (see [Figure 18](#fig:bdsky_principle)). The rates are then allowed to change between two intervals. Within an interval rates are constant. In principle, all rates in all intervals could be different. Since the transmission rate and the becoming noninfectious rate are highly correlated, this is not always practical. Often we assume the becoming noninfectious rate to be the same in all intervals while changing the transmission rate {% eqinline \lambda %}. 

<figure>
	<a id="fig:bdsky_principle"></a>
	<img style="width:50%;" src="figures/bdsky_intervals5.png" alt="">
	<figcaption>Figure 18: Example tree where the red dotted lines are an example of where rates could be allowed to change on the tree. The branch at the root (compare Figure 6) is indicating the origin of the epidemic, which is also estimated in the BDSKY.</figcaption>
</figure>
<br>


There are some clear differences between the birth-death and the coalescent skyline. First, the way intervals are defined. In the coalescent skyline, intervals are always between coalescent events, while this restriction does not exist for the birth-death skyline. Second, the birth-death skyline does not infer changes in population sizes.

The coalescent on the other hand does infer the effective population size, which is a parameter proportional to absolute sizes. The third difference is the inference of the origin of an epidemic by the birth death model, which is not done by the coalescent.


### Visualizing the Birth-Death Skyline Output

There is no equivalent visualization of the skyline plot of a Birth-death Skyline (BDSKY) analysis in Tracer as there is for the Coalescent Bayesian Skyline. But because BDSKY separates the full tree into equally spaced intervals, we can already get an idea of the inference just by looking at the inferred {% eqinline R_e %} values (see [Figure 19](#fig:bdsky_dynamics)). This gives us a good idea of the trend, but it is not completely accurate. Since we are also estimating the origin parameter, the interval times are slightly different in each posterior sample. Thus, the different intervals overlap slightly. The advantage of this is that we get a smooth estimate through time. The disadvantage is that we need to do some extra post-processing to plot the smooth skyline.

As with the Coalescent Bayesian Skyline, because we shortened the chain most parameters have very low ESS values. If you like, you can compare your results with the example results we obtained with identical settings and a chain of 30,000,000 (`hcv_bdsky_30M.log`).  

<figure>
	<a id="fig:bdsky_dynamics"></a>
	<img src="figures/bdsky_tracer.png" alt="">
	<figcaption>Figure 19: Estimated population dynamics by BDSKY in Tracer.</figcaption>
</figure>
<br>


We will use the R-package `bdskytools` to plot the output of the BDSKY analysis. The package is still in development and currently not available over CRAN. Thus, we have to install the package directly over GitHub (note that you only have to install the package once): 

```{R}
install.packages("devtools")
library(devtools)

devtools::install_github("laduplessis/bdskytools")
```

Once the package is installed we have to load the package into our R workspace before we can use the functions in the package.
To plot the results, we need to first tell R where to find the `*.log` file of our run and then load it into R (discarding 10% of samples as burn-in). If you are using RStudio, you can change the working directory to the directory where you stored your log files, which makes it easier to load the files in R.

```{R}
library(bdskytools)

# Navigate to Session > Set Working Directory > Choose Directory (on RStudio)
# or change fname to the full path to the log file
fname <- "hcv_bdsky.log"   
lf    <- readLogfile(fname, burnin=0.1)
```

Next, we can extract the HPDs of {% eqinline R_e %} and the becoming uninfectious rate: 

```{R}
Re_sky    <- getSkylineSubset(lf, "reproductiveNumber")
Re_hpd    <- getMatrixHPD(Re_sky)
delta_hpd <- getHPD(lf$becomeUninfectiousRate)
```


Next we plot the raw HPD intervals of {% eqinline R_e %}. This is equivalent to the output in Tracer. 

```{R}
plotSkyline(1:10, Re_hpd, type='step', ylab="R")
```

<figure>
	<a id="fig:bdsky_hpds"></a>
	<img style="width:80%;" src="figures/bdsky_hpds.png" alt="">
	<figcaption>Figure 20: The HPDs of {% eqinline R_e %} (equivalent to the previous figure from Tracer).</figcaption>
</figure>
<br>

In order to plot the smooth skyline we have to marginalise our {% eqinline R_e %} estimates on a regular timegrid and calculate the HPD at each gridpoint. It is usually a good idea to use a grid with more cells than the dimension of {% eqinline R_e %}. To do this we first calculate the marginal posterior at every time of interest using the function `gridSkyline` and then calculate the HPD for each of the finer time intervals. The times to grid the skyline on (`timegrid`), refers to years in the past. 

```R
timegrid       <- seq(0,400,length.out=101)
Re_gridded     <- gridSkyline(Re_sky,    lf$origin, timegrid)
Re_gridded_hpd <- getMatrixHPD(Re_gridded)
```

Now we are ready to plot the smooth skyline (remember that the sequences were sampled in 1993):

```R
times     <- 1993-timegrid
plotSkyline(times, Re_gridded_hpd, type='smooth', xlab="Time", ylab="R")
```

<figure>
	<a id="fig:bdsky_smooth"></a>
	<img style="width:80%;" src="figures/bdsky_smooth.png" alt="">
	<figcaption>Figure 21: The smooth {% eqinline R_e %} skyline.</figcaption>
</figure>
<br>

We can plot the gridded {% eqinline R_e %} skyline (not its HPDs) for a few of the MCMC samples to see what it really looks like as the Markov chain samples parameters. Note that the intervals overlap between different posterior samples. This is because the origin is different in each of the plotted samples. As we add more samples to the plot we start to see the smooth skyline appear. 

```R
plotSkyline(times, Re_gridded, type='steplines', traces=1, col=pal.dark(cblue,1),ylims=c(0,5), 
            xlab="Time", ylab="R", main="1 random sample")
plotSkyline(times, Re_gridded, type='steplines', traces=10, col=pal.dark(cblue,0.5),ylims=c(0,5), 
            xlab="Time", ylab="R", main="10 random samples")
plotSkyline(times, Re_gridded, type='steplines', traces=100, col=pal.dark(cblue,0.5),ylims=c(0,5), 
            xlab="Time", ylab="R", main="100 random samples")
plotSkyline(times, Re_gridded, type='steplines', traces=1000, col=pal.dark(cblue,0.1),ylims=c(0,5), 
            xlab="Time", ylab="R", main="1000 random samples")
```

<figure>
	<a id="fig:bdsky_traces"></a>
	<img style="width:80%;" src="figures/bdsky_traces.png" alt="">
	<figcaption>Figure 22: Increasing the number of traces plotted from 1 to 10, to 100, to 1000.</figcaption>
</figure>
<br>

Finally, we can plot both the {% eqinline R_e %} and {% eqinline \delta %} (the becoming uninfectious rate) on a single set of axes. Since we left the dimension of the becoming uninfectious rate at 1, it is constant through time. (Normally we would not plot constant parameters over a time period)! The output should be similar to [Figure 23](#fig:bdsky_out).

```R
par(mar=c(5,4,4,4)+0.1)

plotSkylinePretty(range(times), as.matrix(delta_hpd), type='step', axispadding=0.0, 
                  col=pal.dark(cblue), fill=pal.dark(cblue, 0.5), col.axis=pal.dark(cblue), 
                  ylab=expression(delta), side=4, yline=2, ylims=c(0,1), xaxis=FALSE)

plotSkylinePretty(times, Re_gridded_hpd, type='smooth', axispadding=0.0, 
                  col=pal.dark(corange), fill=pal.dark(corange, 0.5), col.axis=pal.dark(corange), 
                  xlab="Time", ylab=expression("R"[e]), side=2, yline=2.5, xline=2, xgrid=TRUE, 
                  ygrid=TRUE, gridcol=pal.dark(cgray), ylims=c(0,3), new=TRUE, add=TRUE)
```

<figure>
	<a id="fig:bdsky_out"></a>
	<img style="width:80%;" src="figures/bdsky_output_pretty.png" alt="">
	<figcaption>Figure 23: Estimates of the inferred {% eqinline R_e %} (orange) over time and the estimate of the becoming uninfectious rate (blue), for which we only used one value.</figcaption>
</figure>
<br>

An R-script with the above commands (and a few extras) is in the `scripts/` directory (`Skyline_Example.R`). 
If the bdskytools package cannot be installed from GitHub the relevant scripts are also provided in the `scripts/` directory.


----

## Some considerations for using skyline plots

Both the coalescent and the birth-death skylines assume that the population is well-mixed. That is, they assume that there is no significant population structure and that the sequences are a random sample from the population. However, if there is population structure, for instance sequences were sampled from two different villages and there is much more contact within than between villages, then the results will be biased {% cite Heller2013 --file Skyline-plots/master-refs %}. Instead a structured model should then be used to account for these biases.


----

# Useful Links

- [Bayesian Evolutionary Analysis with BEAST 2](http://www.beast2.org/book.html) {% cite BEAST2book2014 --file Skyline-plots/master-refs %}
- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)
- BEAST 1 website and documentation: [http://beast.bio.ed.ac.uk](http://beast.bio.ed.ac.uk)
- Join the BEAST user discussion: [http://groups.google.com/group/beast-users](http://groups.google.com/group/beast-users) 

----

# Relevant References

{% bibliography --cited --file Skyline-plots/master-refs %}

