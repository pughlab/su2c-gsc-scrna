## Gradient of developmental and injury-response transcriptional states defines roots of glioblastoma heterogeneity

*Laura M. Richards, Owen Whitley, Graham MacLeod, Florence M.G. Cavalli, Fiona J. Coutinho , Michelle Kushida, Nataliia Svergun, Danielle Croucher, Kenny Yu, Naghmeh Rastegar, Moloud Ahmadi, Danielle A. Bozek, Julia E. Jaramillo, Naijin Li, Erika Luis, Nicole I. Park, Julian Spears, Michael D. Cusimano, Benjamin Haibe-Kains, H. Artee Luchman, Samuel Weiss, Stephane Angers, Peter B. Dirks, Gary D. Bader, Trevor J. Pugh*

Glioblastomas (GBM) harbour diverse populations of cells, including a rare subpopulation of glioblastoma stem cells (GSCs) that drive tumorigenesis. To characterize functional diversity within the tumor-initiating fraction of GBM, we performed single-cell RNA-sequencing on >69,000 GSCs cultured from the tumors of 26 patients. We observed a high degree of inter- and intra-GSC transcriptional heterogeneity that could not be fully explained by somatic alterations at the DNA level. Instead, we found GSCs mapped along a  transcriptional gradient spanning two cellular states reminiscent of normal neural developmental and inflammatory wound response processes. Genome-wide CRISPR/Cas9 drop-out screens recapitulated this observation, as essential genes within each GSC comprised a mixture of these states depending on the position along this axis as measured by gene expression profiling. Further single-cell RNA-sequencing of patient tumors found the majority of cancer cells organize along an astrocyte maturation gradient orthogonal to the GSC gradient yet retained expression of founder transcriptional programs found in GSCs. Our work supports a model whereby GBM heterogeneity is rooted in a fundamental GSC-based neural tissue regeneration and wound response transcriptional continuum.

***



**Cohort Overview:**  
XX cells from 29 BTSC culturers, derived from XX patients.     
XX cells from 24 primary GBM samples from 7 patients (some tumours have multiple regions/sample)  
XX/XX cells from primary tumours are malignant, XXX are normal brain, and XXX are of immune origin.   
X samples have matched BTSC culture and primary GBM.   
  
**Manuscript:**   
https://docs.google.com/document/d/1qLlc83DX5s8qPE3pJI2Rfcaq0i0Zkn0_2sxnk3msaGY/edit?usp=sharing    

---
##  1.0 Scripts

runCellRanger.sh  
> Run CellRanger   
  
runBamTagHistogram.sh    
>     

Dropbead.R    
QC.R    
SeuratPreprocessing.R    

##  2.0 Figures

Figures for this manuscript can be accessed on google drive:
https://drive.google.com/open?id=1o8pcIJsZ3Znnoy2uowlab0VcohdVvrVp


##  3.0 Data

To easily access and query single cell data related to this study, we have uploaded data to the Broad Single Cell Portal.   
Available data includes:
> Expression Matrix for Global BTSC data set  
> Meta-data for Global BTSC data set  
> Expression Matric for Global Tumour data set  
> Meta-data for Global BTSC data set  
> Gene signatures used in study  


For upload to EGA, we will be depositing fastqs for all scRNA-seq, bulk RNA and WGS.


