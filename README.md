# EstimateChromosome

This repository was created for Cranfield BIX Group Project, Group 2.
It contains the R script that performs chromosome estimation on the server. 
Two versions of the code exist on the server, EstimateChromosomes_v2.R which reads a fasta file, and Delta_EstimateChromosomes.R which reads a delta file.
Both are run via command line.

In both scripts, the code extracts the lengths of each sequence from its header and is clustered by kmeans clustering.
The cluster with the largest y-value center is deemed to be the chromosomes. These are sorted by ID and exported to a CSV file.

Jia Ying Chua
April 2023
