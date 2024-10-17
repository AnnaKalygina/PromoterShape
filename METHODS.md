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

To shuffle the sequences, the [shuffle_fasta.py](shuffle_fasta.py) script was used. This script parses FASTA files and shuffles the nucleotides within each sequence independently, resulting in the same nucleotide composition but with the bases randomly ordered.

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

# Statistics

### Standartization
When working with the shapes of different promoter sequences, we want the values to be normalized.
Promoters within and across species have high levels of variation, but we want to be able to assess the base pair difference while keeping the promoter-wise variation negligible. For this purpose the values across base pairs in each promoter are z-scored and then the average across promoters is calculated. With this approach we preserve the base pair -wise variation while centering the values around $0$ for more more efficient cross-species comparison. 

The normalisation yields the same shape but on a different scale:
![image](img/normalised_comparison.png) 

### Measuring conservation
https://mathoverflow.net/questions/140813/what-is-a-good-algorithm-to-measure-similarity-between-two-dynamic-graphs
