The first part of the code includes the section where I create the UMAP plots based on the bias-corrected data.
I do this through taking the average co-essentiality scores from the genes in the tissues of interest and then using these values
for the coordinates and for the color to match the original paper (although I changed the colors to red blue instead of yellow-blue
for easier visualization). I also went through all 20 cancer types shown in Jessie's figure, analyzed the top 3 most co-expressed genes,
and clustered these together to gain a broader biological understanding of the potential genetic mechanisms at play. This led to
figure 7a. I manually annotated the genes with the broad functions based on GO terms and the pathways in which they are involved in 
Supplementary Table 2. From here, I plotted the clusters as a bar chart to show the relative size of each type of biological pathway
for highly co-expressed genes in cancer.

My code document is heavily commented with more details about how I do each of the steps and what I am doing at each step.

I have a .gitignore file in my branch that includes the file with the data I am analyzing. It was too large to upload into code,
so I used the .gitignore feature so that it wouldn't create issues with the overall repository. In this .gitignore file, I have a
folder called data and another file called omnipathr-log/omnipathr-20251130-1530.log. Data contained the files "bias_corrected_data.csv"
and "DepMap-2018q3-celllines.csv". The former contains a matrix of CERES scores for each of the genes from every cell line,
after the bias-correction done in Yubi's step of the data processing. The latter file contains information about all of the cell
lines used in the study. While I didn't directly use this file in my code, I read through the file carefully to get a stronger
understanding about what these cell lines represented. Finally, the file omnipathr-log/omnipathr-20251130-1530.log was downloaded
to my VSCode locally during my coding when I was originally planning to replicate figure 2a rather than extend the analysis of the
genes. I ultimately did not use this file, but I received several errors when I tried to delete it either locally or globally from
VSCode or GitHub. I moved it to the .gitignore file so that it would not get in the way of the rest of my branch. However, it was
not used and is not a component of this QCB455 final project.
