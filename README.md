# Semi-Manual-MD-ECD-De-Novo

Anthropic’s Claude Sonnet 4.6 and Anaconda’s Jupyter Notebook were used to create Python 3 code for the graphing of deconvoluted top-down ECD data into averagine scaled mass defect plots (averagine scaled mass defect vs monoisotopic mass). The code then allowed the user to visually assign c- and z-ions from the plot produced. After the assignment, the code calculated the mass differences between all the possible c-ion fragments. It then identified amino acid residues that match the observed mass differences, and outputted the residue with the condition that the mass difference falls within a 10 ppm error range. This corresponds to an absolute mass error of no more than 2 mDa. This process was then repeated by the code for the assigned z-ions.

If using Anaconda’s Jupyter Notebook and figures are not interative run line 2 of code in a separate cell and fully restart the application.

This code was developped for a masters thesis at the University of Edinburgh.
