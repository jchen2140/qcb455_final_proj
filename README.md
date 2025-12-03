The first part of the code includes the section where I create the UMAP plots based on the bias-corrected data.
I do this through taking the average co-essentiality scores from the genes in the tissues of interest and then using these values
for the coordinates and for the color to match the original paper (although I changed the colors to red blue instead of yellow-blue
for easier visualization). I also went through all 20 cancer types shown in Jessie's figure, analyzed the top 3 most co-expressed genes,
and clustered these together to gain a broader biological understanding of the potential genetic mechanisms at play. This led to
figure 7a. I manually annotated the genes with the broad functions based on GO terms and the pathways in which they are involved in 
Supplementary Table 2. From here, I plotted the clusters as a bar chart to show the relative size of each type of biological pathway
for highly co-expressed genes in cancer.

My code document is heavily commented with more details about how I do each of the steps and what I am doing at each step.
