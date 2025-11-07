# Aligned_data_analysis_SynapseClustering
Data analysis for 'Eye-specific differences in active zone addition during synaptic competition in the developing visual system'
Zhang et al. eLife 2023;12:RP91431. DOI: https://doi.org/10.7554/eLife.91431

## Cluster Volume: 
Extract volume information from complex or simple clusters for visualization and statistical analysis. 
## Dist_Matrix: 
An old method for data inspection. From each individual synaptic clsuter, inspect 20 cluster nearby and list their ordered distances. 

At the index of 2-3, significant differences in distance can be found for samples of various ages/genotypes. It verifies the local synapse arrangement inlcude 2-3 synaptic clusters, which are typically complex synapses. 
## RDF_Function: 
Calculate the radius distribution function from each individual synapses. This indicates whether complex synapse is a common local feature in this sample and offers a (rather complex) method to reveal typicl size of a complex synapse. 
## Eye-specific shell: 
Shell analysis of complex synapses similar to [previous publication](https://www.cell.com/cell-reports/fulltext/S2211-1247(23)00096-7). 

The analysis shows no significant results, thouhgh. 
## Large_scale_local_density: 
A large-scale correlation analysis of spatial distribution of synapses. It creates Z-projection of the image stack using customized parameters, montage the projection, and calcualte local denstiy self-correlation. A comparison of original and randomized synapse distritution is included. 

In the SCN, this method reveals signification results (indicated by the correlation coefficient R) between WT and OPN4 Cre/Cre mice, indicating a local synapse arrangement at the scale of um. But it does not show signification differences in the dLGN dataset. 
## Eye_specific
Core analysis
* Dist_matrix: similar to 'Dist_matrix' described above.
* Non-ret: Complex synapse analysis for non-retinal input.
* Complex_Syn_Density.m and Complex_Syn_Batch:
Calculate the size of complex synapses.
* Cluster_Around_Complex*.m
Inspect eye-specific simple clusters that are close to a complex synapses from the same or different eye. The inspection includes ratio/number of simple synapses and their properties depending on their relation to the closest complex synapse.
* Simp_near_comp_comp_*.m
Similar as above, but inspecting the properties of complex synapses.
