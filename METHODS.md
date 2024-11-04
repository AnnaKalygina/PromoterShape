# Methods
### Raw data collection
The following species are chosen for the promoter analysis: *Homo sapiens* (29598 promoters), *Drosophila melanogaster* (16972 promoters), *Caenorhabditis elegans* (7120 promoters), *Plasmodium 
falciparum* (5597 promoters), *Arabidopsis thaliana* (22703 promoters), and *Saccharomyces cerevisiae* (5117 promoters).

The promoter sequences are downloaded from the [Eukaryotic Promoter 
Database (EPD)](https://epd.expasy.org/epd/). The selection was restricted 
to the EPDnew IDs. The sequence of each promoter includes -2000 base pairs 
upstream and 2000 base pairs downstream of a reported TSS. 

### Controls

There are two types of controls I am using in this analysis:
- random shuffling 
- HMM generation

**Random shuffling**

To generate controls for statistical analysis, all DNA sequences were randomly shuffled, preserving the nucleotide composition while eliminating any positional information. These shuffled sequences serve as a control set for comparison with the real promoter sequences.

To shuffle the sequences, the [shuffle_fasta.py](/py_scripts/shuffle_fasta.py) script was used. This script parses FASTA files and shuffles the nucleotides within each sequence independently, resulting in the same nucleotide composition but with the bases randomly ordered.

**HMM generation**

Another method to create control sequences is by using HMMs (Hidden Markov Models) to generate random sequences based on a profile built from the input sequences. This approach conserves the nucleotide composition based on positional probabilities derived from the input sequences.

The following commands are used to first build an HMM profile from a multiple-sequence alignment in a raw FASTA file and then generate random sequences based on this profile:

``` bash
hmmbuild --dna control_promoter_sequences/athaliana_200_hmm_profile.txt raw_promoter_sequences/athaliana_200.fa 
hmmemit -o control_promoter_sequences/athaliana_200_hmm.fa -N 22703 \ 
control_promoter_sequences/athaliana_200_hmm_profile.txt 
```
- `hmmbuild`: creates an HMM profile from the multiple-sequence alignment in the raw_promoter_sequences/athaliana_200.fa file.
- `hmmemit`: generates random sequences based on the HMM profile. The -N 22703 option specifies the number of sequences to generate, matching the original number of sequences.

P.S: A problem was encountered at this step: the `hmmemit` produces non-aligned sequences of varying length, even though the hmm_profile specifies the required sequence to be 401 nucleotides long. For example, there are 38 (out of 29k) sequences in __hsapiens__ that are longer than the requested length. Therefore a higher total number of promoters was generated through `hmmemit` and then the relevant number of sequences with a consistent length of 401 nucleotide was taken.

**Skylign comparison**
The nucleotide frequencies (skylign) of all tree alignments: raw promoter sequences, hmm profile and shuffled sequences were generated and compared:
1. Raw promoter sequence 
![image](skylign_logos/scerevisiae_raw_skylign_logo.png)
2. HMM sequences
![image](skylign_logos/scerevisiae_hmm_skylign_logo.png)
3. Shuffled sequences
![image](skylign_logos/scerevisiae_shuffled_skylign_logo.png)

### Promoter shape prediciton
The shape of the promoters from the six species was analyzed to investigate conservation. Promoter shape prediction was focused on the 400 bp region surrounding the TSS. The  [deepDNAshape package](https://github.com/JinsenLi/deepDNAshape/blob/main/README.md) was used to predict DNA shape properties for each promoter with the following script:

``` bash
species=( athaliana celegans dmelanogaster hsapiens pfalciparum scervisiae )
properties=( MGW Shear Stretch Stagger Buckle ProT Opening Shift Slide Rise Tilt Roll HelT )

input_dir="/local/home/quee4387/raw_promoter_sequences/"
output_dir="/local/home/quee4387/dna_shape/"
for spec in "${species[@]}"
do
    for prop in "${properties[@]}"
    do
        echo "Predicting $prop for $spec"
        python /local/home/quee4387/.conda/envs/myenv/bin/deepDNAshape \
        --file ${input_dir}${spec}_2000.fa \
        --feature $prop \
        --output ${output_dir}${spec}_${prop}_2000.txt
    done
done
```
In this analysis, the following DNA shape features were predicted:
- MGW (hufflingminor groove width)
- Shear
- Stretch
- Stagger
- Buckle
- ProT (propeller twist)
- Opening
- Shift
- Slide
- Rise
- Tilt
- Roll
- HelT (hellical twist)

This script generates 78 output `.txt` files, one for each combination of species and DNA shape property.

### Statistical analysis
#### Collapsing predicted shape matrices 
The predicted shape matrices contain $n$ rows and $m$ columns, where $n$ - number of promoters fed in the model, $m$ - length of promoter sequences. The dataset could be collapsed into $1 \cdot m$ matrix by averaging parameter prediction across all promoters. The parameters could be averaged by simply taking the mean across all promoter sequences per position, or by taking a z-score across all promoter sequences per position.

The data for each species and each property is stored in a [long table format](/combined_dna_shape_table.csv) with the following columns: species, property, source (raw, hmm, shuffled), position and value. This data is later used for plotting comparison graphs. The [combine_shapes.py](/py_scripts/combine_shapes.py) was used for combining predicitons into single dataframe.

### Comparing raw promoter sequences with the controls
While shuffled sequences didn't show any descrete shape pattern, the hmm emited sequnces followed a pattern almost identical to the one displayed by the real raw promoters.

![image](img/raw_vs_control_athaliana_MGW.png)
The following [code](/r_scripts/raw_hmm_shuffled_comparison_plot.R) was used for plotting the graph above.

In conclusion, the results suggest that nucleotide composition influences the shape of a promoter sequence. So the shape features are the product of nucleotide composiiton, however, in some cases the shapes of raw promoters is more exaggerated than of the hmm sequences.
Next question is: when comparing different species do we see that the shape is more conserved than the nucleotide compositon among species?

### Plotting shape conservation within species
To visualize shape conservation, the .txt files were randomly sampled, resulting in matrices of size 1000 × 400, where 1000 represents the number of promoters and 400 represents the length of each promoter sequence. These matrices were used to plot the predicted shape features for each species and property, providing insights into promoter shape conservation across the selected species.
The following command was used: 

``` bash 
shuf -n 1000 <athaliana_MGW_200.txt> > <athaliana_MGW_200_sample_1k.txt>
```
The following R script is used to plot the figures:

``` R
library(ggplot2)
library(reshape2)
library(gridExtra)

properties <- c("MGW", "Buckle","Opening", "Tilt")


input_dir <- "~/Downloads/athaliana_"
plot_list <- list()

for (prop in properties) {

  file_path <- paste0(input_dir, prop, "_200.txt")
  data <- read.table(file_path, header = FALSE)
  average_data <- colMeans(data)
  positions <- seq(-200, 200, length.out = ncol(data))
  average_df <- data.frame(Position = positions, Mean_Value = average_data)

  p <- ggplot(average_df, aes(x = Position, y = Mean_Value)) +
    geom_line(color = "blue") +
    theme_minimal() +
    labs(x = "Position", y = "Average Feature Value", 
         title = paste("Average DNA Shape Prediction for", prop)) +
    scale_x_continuous(breaks = seq(-200, 200, 50))  # Adjust x-axis ticks

  plot_list[[prop]] <- p
}

grid.arrange(grobs = plot_list, ncol = 1)
```
The following figures are obtained:
![image](https://github.com/user-attachments/assets/fe443fd8-0c46-4694-87a5-aa164227bcb3)

### Data storage and accessing 


# Statistics
### Within species conservation 
I hypothesize that the shapes of the promoters within species are conserved. In other words, I want to show that the values of a given shape parameter across different promoters at a given position are consistent and are not significantly different. For this purpose, I perform tests to assess low variability across promoters at each position. The shape of promoter sequences is also compared to control sequences, with the expectation that real promoters have lower variability than the randomly shuffled or HMM-generated control sequences. 

The following tests are used:
- **Coefficient of variation**:
$$CV = \dfrac{standard deviation}{mean}$$
The coefficient of variation is calculated for each base pair position across promoters. A low CV indicates that the shape parameter at a specific position is consistent across different promoters, suggesting conservation. This test is performed on both the real promoter sequences and the control sequences (randomly shuffled or generated by HMMs).

  - Real promoters: if the CV is consistently low across all positions, it supports the hypothesis that the shape is conserved at specific regions, especially around the transcription start site (TSS).
  - Control sequences: higher CV values for shuffled or HMM-generated sequences would indicate that the real promoters are more conserved in their shape, confirming the hypothesis of conservation.

- **One-way ANOVA**
One-way ANOVA is used to test whether the mean shape parameter values at a given base pair position are significantly different across all promoters. The null hypothesis is that the mean values at each position are the same across promoters, which would support conservation.

  - Null Hypothesis (H0): The means of the shape parameter values across promoters at the same position are equal.
  Alternative Hypothesis (H1): The means are significantly different.

A non-significant ANOVA result (high p-value) suggests that the shape parameter values at the same position are not significantly different, supporting the conservation hypothesis. For this analysis, ANOVA is run across all base pair positions.


### Across species conservation
**Is shape more conserved than the sequence?** 


I am aiming to demonstrate that promoter shapes are more conserved than promoter nucleotide composition across different eukaryotic species.
Through comparative statistical analyses, I will quantify the extent of conservation in shape features and nucleotide composition to establish whether promoter shapes exhibit higher consistency across species.

Null Hypothesis (H0): There is no significant difference in conservation between nucleotide composition and shape features across species.

Alternative Hypothesis (H1): Shape features are more conserved than nucleotide composition across species.

#### Comparing nucleotide composition across species

#### Standartization
When working with the shapes of different promoter sequences across many different species, we want the values to be normalized.
Promoters within and across species have high levels of variation, but we want to be able to assess the base pair difference while keeping the promoter-wise variation negligible. For this purpose the values across base pairs in each promoter are z-scored and then the average across promoters is calculated. With this approach we preserve the base pair -wise variation while centering the values around $0$ for more more efficient cross-species comparison. 

The normalisation yields the same shape but on a different scale:
![image](img/normalised_comparison.png) 

### Measuring conservation
https://mathoverflow.net/questions/140813/what-is-a-good-algorithm-to-measure-similarity-between-two-dynamic-graphs

Here is a species-wise comparison:
![image](img/comparison_buckle_dm_hs.png)

### Kullback–Leibler divergence
The KL divergence is used to calculate how the shape probability distributions and nucleotide composition distributions are different across species (T.D. Pham & J.Zuegg, 2004).

The apprdoach used derives a probabilistic distance between two sequences using a symmetrized version of the KLD, which directly compares two Markov
models built for the two corresponding biological sequences. It measures the similarity between any two sequences by computing the log-likelihood difference between two Markov models with the same observation data.

- SimMM (SIMilarity of Markov Models)
- CLUSTAL W