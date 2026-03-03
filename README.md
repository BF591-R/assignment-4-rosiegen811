# Assignment Data Science Basics

## Problem Statement

High dimensional data is typically time and resource-intensive to analyze,
interpret and visualize. Gene expression data (including microarray data) often
consists of measurements for tens of thousands of known genes for every sample.
Learning methods to inspect, reduce and display high dimensional data is
necessary for many machine learning and bioinformatics problems. It is also a
critical exploratory data analysis tool and one of the most widely used
methodologies for this purpose. 

## Required Readings 

Please read sections 7.4 and 8.5.6-8.5.7 in the textbook.

## Learning Objectives
1. Understand the basic principles of Principal Component Analysis (PCA) and 
   hierarchical clustering / heatmaps
2. Perform PCA in R and evaluate the basic outputs of PCA
3. Generate basic clustered heatmaps

## Skill List
1. R: `scale()`, transpose `t()`, `prcomp()`, `heatmap()`, ggplot2

## Background on Microarrays

A microarray consists of thousands of specifically designed single stranded DNA 
sequences affixed or bound to a solid surface (glass, nylon, etc.). Extracted RNA
 or DNA from samples of interest are labeled with fluorescent dyes, and hybridized 
with the bound probes by flowing across the surface. Non-hybridized molecules are 
washed away and b a laser will excite the attached dye to produce light to be detected
by a scanner and converted to a digital image. Finally, image processing will 
transform each spot (bound probe) to a numerical measurement that can be used to 
infer expression levels for genes.

## Background on Principal Component Analysis

Principal component analysis is an exploratory data analysis technique commonly 
used to reduce the dimensionality of data while trying to simultaneously minimize 
information loss. To understand the inner workings of PCA will require an in-depth
review of linear algebra and is beyond the scope of this specific assignment. 
We will provide several links and references if you wish to do so on your own. 
At its core, PCA creates new uncorrelated linear combinations of the original 
variables that successively maximize variance. Thus, the first principle component 
(PC) represents the direction of the data that explains a maximal amount of variance.
The second PC represents the direction that captures the second most variance and 
so on and so forth. As you can infer, if a small number of PCs capture a majority 
of variance contained within a dataset, they can be used as a lower dimensional
representation of the data.

## Marisa et al. Gene Expression Classification of Colon Cancer into Molecular Subtypes: Characterization, Validation, and Prognostic Value. PLoS Medicine, May 2013. PMID: 23700391

The example intensity data was taken from the listed publication. In it, the
authors proposed the use of gene expression profiling (through microarray
technology) to generate a robust and reproducible classifier to identify
subtypes of colorectal cancer samples. They identified six subtypes that they
demonstrated to be significantly associated with distinct molecular pathways and
clinical pathologies. We have provided you a sample expression matrix that was
normalized and subjected to batch correction. We will use the
`data/example_intensity_data.csv` as our data matrix going forward.

## Scaling data using R `scale()`

Oftentimes, it is necessary to analyze or plot data that contains variables that
differ by multiple orders of magnitude. In bioinformatics, this is an extremely
common situation as different genes may be expressed at wildly different levels.
For instance, in typical HTS experiments and as you will see later on in the
course, you may detect genes with counts ranging from the single digits to those
with counts in the tens of thousands. Heatmaps of gene expression data are
commonly transformed by standard scaling in order to improve the visualization
while retaining the underlying pattern of expression between samples.

One such way to standardize data is the `scale()` function in R. Scale
essentially takes a vector of values and determines the mean and standard
deviation of the entire vector. Scale() will then subtract the mean from each
element in the vector and divide each by the standard deviation. Z-scores thus
correspond to the number of standard deviations by which a data point is above
or below the mean of all the observations (i.e. a Z-score of 0 represents a
value equal to the mean and a Z-score of 1 represents a value one standard
deviation greater than the mean value). This has the effect of representing all
of your data points on the same "scale" while preserving the pattern or profile
of differences inherent between values across different observations.

The [transpose function `t()`](https://bu-bioinfo.github.io/r_for_biological_sciences/prog-basics.html#matrices)
is a function covered in the textbook for this course. As covered earlier, it
converts a $m \times n$ matrix to $n \times m$. The built-in R function
`scale()` operates on a column-wise basis. 

## Proportion of variance explained

As mentioned earlier, principal components represent a successively maximal
amount of variance in your dataset. It is often helpful to visualize the
percentage of total variance explained by each principal component.

## Plotting and visualization of PCA

In typical HTS experiments, PCA is a common tool used to analyze the similarity
between samples and determine if the experimental variable of interest
(genotype, knockout, etc.) show any patterning or clustering in the PC biplot.
The biplot is one of the most common visualizations and is created by plotting
the first two principal components (which may not always represent a majority of
variance in your dataset) against each other. Typically, one will examine the
pattern of clustering in this plot to see if samples are grouped together in any
meaningful pattern according to important experimental variables. It is
important to remember that these components are linear combinations of the
original features and are not easily interpretable. The value or score of each
sample in terms of its principal components may be found in the pca results
object accessed by $x.

## Hierarchical Clustering and Heatmaps

Now that we've discussed the basics of PCA, we will move on to discuss the use
of hierarchical clustering and heatmaps. Heatmaps are a common plot used to
quickly visualize the expression and pattern of many relevant observations
across samples. For HTS data, they essentially depict the expression of genes
(rows) for each sample (columns) as colors that denote the magnitude of
expression to enable quick identification of patterns and changes between
samples/experimental groups. Prior to constructing a heatmap, it is common to
use some form of clustering algorithm to learn potential groups and patterns in
the expression data.

One such method is known as hierarchical clustering (specifically focusing on
agglomerative clustering), a method which attempts to cluster data into
hierarchies where single observations begin as separate points and are
successively merged into larger clusters until all points have been grouped.
Typically, agglomerative clustering will begin with a dissimilarity matrix that
defines how "close" all observations are by calculating some "distance" metric
for every pair of observations. There are multiple "distance" metrics that can
be chosen, and one of the most commonly used is simple euclidean distance. A
linkage function will then use these distance metrics to define how clusters are
formed based on these distance values. For example, complete linkage will define
clusters by the maximum value of all pairwise distances between two different
clusters.   

We will be using the base R function `heatmap()` which by default
uses euclidean distance and the `hclust()` function to produce a clustered
heatmap. The `heatmap()` function will also scale your data row-wise to aid in
visualization.

To begin, we have provided you with a representative output of differentially
expressed probes (`data/differential_expression_results.csv`) between the C3 and C4
subtypes. 

## References

<small>
https://www.genome.gov/about-genomics/fact-sheets/DNA-Microarray-Technology

Govindarajan, R., Duraiyan, J., Kaliyappan, K. & Palanisamy, M. Microarray and its applications. Journal of Pharmacy & Bioallied Sciences 4, S310 (2012).

Trevino, V., Falciani, F. & Barrera-Saldaña, H. A. DNA microarrays: A powerful genomic tool for biomedical and clinical research. Molecular Medicine 13, 527–541 (2007).

Shlens, J. A Tutorial on Principal Component Analysis.(https://arxiv.org/abs/1404.1100)

Lever, J., Krzywinski, M. & Altman, N. Points of Significance: Principal component analysis. Nature Methods 14, 641–642 (2017).

Hastie, T., Hastie, T., Tibshirani, R., & Friedman, J. H. (2001). The elements of statistical learning: Data mining, inference, and prediction. New York: Springer.https://hastie.su.domains/ElemStatLearn/printings/ESLII_print12_toc.pdf

</small>
