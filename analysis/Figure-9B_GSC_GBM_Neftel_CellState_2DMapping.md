
---
# Project GSCs and Tumour cells onto cell state map
---

No G800_L

"/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/Neftel_GBM_Mapping"

---

## 1.0 Score all GSCs and Tumour cells with sigs

---

Use addmodule score -- already done

---
## 2) Define y-axis (D): classify cells as OPC/NPC or AC/MES
---

3. Average NPC scores and MES scores to combine groups (now n=4)
4. D = max(SCopc,SCnpc) - max(SCac,SCmes)
5. If D is negative , then AC-MES, If D is positive (D>0) then OPC,NPC


```R
options(repos='http://cran.rstudio.com/')
library(ggplot2)
library(gridExtra)
#install.packages("viridis")
library(viridis)
library(MASS)
library(ggpubr)
```

    Loading required package: viridisLite
    



```R
setwd("~/Desktop/H4H/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/Neftel_GBM_Mapping")
Neftel <- readRDS("GSC_GBM_scRNA_AddModuleScore_Neftel_meta.rds")
```


```R
### Average NPC1/2 and MES1/2 scores for each cell
Neftel$Neftel_NPC <- apply(cbind(Neftel$Neftel_NPC1, Neftel$Neftel_NPC2), 1, mean)
Neftel$Neftel_MES <- apply(cbind(Neftel$Neftel_MES1, Neftel$Neftel_MES2), 1, mean)
head(Neftel)
```


<table>
<caption>A data.frame: 6 × 23</caption>
<thead>
	<tr><th></th><th scope=col>CellType</th><th scope=col>PatientID</th><th scope=col>SampleID</th><th scope=col>Neftel_MES2</th><th scope=col>Neftel_MES1</th><th scope=col>Neftel_AC</th><th scope=col>Neftel_OPC</th><th scope=col>Neftel_NPC1</th><th scope=col>Neftel_NPC2</th><th scope=col>Neftel_G1.S</th><th scope=col>⋯</th><th scope=col>Developmental_GSC_AUC</th><th scope=col>InjuryResponse_GSC_AUC</th><th scope=col>Developmental_GSC_AUC_z</th><th scope=col>InjuryResponse_GSC_AUC_z</th><th scope=col>Developmental_GSC_AUC_norm</th><th scope=col>InjuryResponse_GSC_AUC_norm</th><th scope=col>Dev_IR_Diff</th><th scope=col>Orig.ID</th><th scope=col>Neftel_NPC</th><th scope=col>Neftel_MES</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl[,1]&gt;</th><th scope=col>&lt;dbl[,1]&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCACGGACAA</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td>-0.01913800</td><td>-0.071709632</td><td> 0.06980817</td><td>-0.03199840</td><td>-0.06716432</td><td> 0.095209179</td><td> 0.23451153</td><td>⋯</td><td>0.1833756</td><td>0.1691589</td><td>1.5502934</td><td>-0.7927835</td><td>0.5502658</td><td>0.2204938</td><td>0.32977199</td><td>BT127_L</td><td> 0.0140224276</td><td>-0.04542382</td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCATCCGGGT</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td>-0.08061522</td><td> 0.033282546</td><td>-0.10055778</td><td>-0.12732257</td><td> 0.07976724</td><td>-0.005820044</td><td> 0.16263892</td><td>⋯</td><td>0.1225561</td><td>0.1561081</td><td>0.1716288</td><td>-1.0550278</td><td>0.3172655</td><td>0.1786540</td><td>0.13861148</td><td>BT127_L</td><td> 0.0369735970</td><td>-0.02366634</td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGGTACAGTTC</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td> 0.11935181</td><td>-0.007389301</td><td>-0.09355375</td><td> 0.13067802</td><td> 0.31634315</td><td> 0.165309543</td><td>-0.08324292</td><td>⋯</td><td>0.1403795</td><td>0.1532729</td><td>0.5756518</td><td>-1.1119984</td><td>0.3855471</td><td>0.1695646</td><td>0.21598251</td><td>BT127_L</td><td> 0.2408263459</td><td> 0.05598125</td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGTCTACGAGT</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td> 0.16918706</td><td> 0.091923570</td><td> 0.02553869</td><td> 0.09741310</td><td> 0.49856022</td><td> 0.477332226</td><td>-0.34891006</td><td>⋯</td><td>0.1531699</td><td>0.1659836</td><td>0.8655863</td><td>-0.8565874</td><td>0.4345473</td><td>0.2103142</td><td>0.22423312</td><td>BT127_L</td><td> 0.4879462244</td><td> 0.13055531</td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGAGTGGTAAT</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td> 0.14816768</td><td> 0.204905229</td><td>-0.02794300</td><td>-0.03860135</td><td> 0.05134542</td><td> 0.006812182</td><td> 0.62688342</td><td>⋯</td><td>0.1330963</td><td>0.1653948</td><td>0.4105539</td><td>-0.8684189</td><td>0.3576448</td><td>0.2084265</td><td>0.14921831</td><td>BT127_L</td><td> 0.0290787990</td><td> 0.17653645</td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGCAGGACGTA</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td> 0.37944621</td><td> 0.255718588</td><td> 0.46389484</td><td> 0.28277900</td><td> 0.04921320</td><td>-0.050914580</td><td>-0.22336134</td><td>⋯</td><td>0.1597614</td><td>0.2143644</td><td>1.0150026</td><td> 0.1155828</td><td>0.4597993</td><td>0.3654192</td><td>0.09438006</td><td>BT127_L</td><td>-0.0008506884</td><td> 0.31758240</td></tr>
</tbody>
</table>




```R
### Define D (aka y-axis)

#take max of these two scores
a <- apply(cbind(Neftel$Neftel_OPC, Neftel$Neftel_NPC), 1, max)
b <- apply(cbind(Neftel$Neftel_AC, Neftel$Neftel_MES), 1, max)

```


```R
Neftel$Y.axis <-  a - b
Neftel$CellClass <- ifelse(Neftel$Y.axis > 0, "OPC.NPC", "AC.MES")
```


```R
### classify cells by their max score

cc <- cbind(Neftel$Neftel_AC, Neftel$Neftel_MES, Neftel$Neftel_NPC, Neftel$Neftel_OPC)
colnames(cc) <- c("AC", "MES", "NPC", "OPC")
#head(cc)

Neftel$MaxClass <- apply(cc,1,function(x) which(x==max(x)))
Neftel$MaxClass <- gsub(1, "AC", Neftel$MaxClass)
Neftel$MaxClass <- gsub(2, "MES", Neftel$MaxClass)
Neftel$MaxClass <- gsub(3, "NPC", Neftel$MaxClass)
Neftel$MaxClass <- gsub(4, "OPC", Neftel$MaxClass)
    
Neftel$ClassSign_X <- Neftel$MaxClass
Neftel$ClassSign_X <- gsub("AC", -1, Neftel$ClassSign_X )
Neftel$ClassSign_X <- gsub("MES", 1,Neftel$ClassSign_X )
Neftel$ClassSign_X <- gsub("OPC", -1, Neftel$ClassSign_X )
Neftel$ClassSign_X <- gsub("NPC", 1, Neftel$ClassSign_X )
head(Neftel)
    
```


<table>
<caption>A data.frame: 6 × 27</caption>
<thead>
	<tr><th></th><th scope=col>CellType</th><th scope=col>PatientID</th><th scope=col>SampleID</th><th scope=col>Neftel_MES2</th><th scope=col>Neftel_MES1</th><th scope=col>Neftel_AC</th><th scope=col>Neftel_OPC</th><th scope=col>Neftel_NPC1</th><th scope=col>Neftel_NPC2</th><th scope=col>Neftel_G1.S</th><th scope=col>⋯</th><th scope=col>Developmental_GSC_AUC_norm</th><th scope=col>InjuryResponse_GSC_AUC_norm</th><th scope=col>Dev_IR_Diff</th><th scope=col>Orig.ID</th><th scope=col>Neftel_NPC</th><th scope=col>Neftel_MES</th><th scope=col>Y.axis</th><th scope=col>CellClass</th><th scope=col>MaxClass</th><th scope=col>ClassSign_X</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCACGGACAA</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td>-0.01913800</td><td>-0.071709632</td><td> 0.06980817</td><td>-0.03199840</td><td>-0.06716432</td><td> 0.095209179</td><td> 0.23451153</td><td>⋯</td><td>0.5502658</td><td>0.2204938</td><td>0.32977199</td><td>BT127_L</td><td> 0.0140224276</td><td>-0.04542382</td><td>-0.05578574</td><td>AC.MES </td><td>AC </td><td>-1</td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCATCCGGGT</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td>-0.08061522</td><td> 0.033282546</td><td>-0.10055778</td><td>-0.12732257</td><td> 0.07976724</td><td>-0.005820044</td><td> 0.16263892</td><td>⋯</td><td>0.3172655</td><td>0.1786540</td><td>0.13861148</td><td>BT127_L</td><td> 0.0369735970</td><td>-0.02366634</td><td> 0.06063994</td><td>OPC.NPC</td><td>NPC</td><td>1 </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGGTACAGTTC</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td> 0.11935181</td><td>-0.007389301</td><td>-0.09355375</td><td> 0.13067802</td><td> 0.31634315</td><td> 0.165309543</td><td>-0.08324292</td><td>⋯</td><td>0.3855471</td><td>0.1695646</td><td>0.21598251</td><td>BT127_L</td><td> 0.2408263459</td><td> 0.05598125</td><td> 0.18484509</td><td>OPC.NPC</td><td>NPC</td><td>1 </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGTCTACGAGT</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td> 0.16918706</td><td> 0.091923570</td><td> 0.02553869</td><td> 0.09741310</td><td> 0.49856022</td><td> 0.477332226</td><td>-0.34891006</td><td>⋯</td><td>0.4345473</td><td>0.2103142</td><td>0.22423312</td><td>BT127_L</td><td> 0.4879462244</td><td> 0.13055531</td><td> 0.35739091</td><td>OPC.NPC</td><td>NPC</td><td>1 </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGAGTGGTAAT</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td> 0.14816768</td><td> 0.204905229</td><td>-0.02794300</td><td>-0.03860135</td><td> 0.05134542</td><td> 0.006812182</td><td> 0.62688342</td><td>⋯</td><td>0.3576448</td><td>0.2084265</td><td>0.14921831</td><td>BT127_L</td><td> 0.0290787990</td><td> 0.17653645</td><td>-0.14745765</td><td>AC.MES </td><td>MES</td><td>1 </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGCAGGACGTA</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td> 0.37944621</td><td> 0.255718588</td><td> 0.46389484</td><td> 0.28277900</td><td> 0.04921320</td><td>-0.050914580</td><td>-0.22336134</td><td>⋯</td><td>0.4597993</td><td>0.3654192</td><td>0.09438006</td><td>BT127_L</td><td>-0.0008506884</td><td> 0.31758240</td><td>-0.18111584</td><td>AC.MES </td><td>AC </td><td>-1</td></tr>
</tbody>
</table>



##### 3) Define the x-axis

1. For AC-MES cells (D < 0), the x axis was defined as log2(jSCac – SCmesj)
2. For OPC-NPC cells (D>0), the x-axis was defined as log2(jSCopc – SCnpcj+1)



```R
Neftel$OPC.NPC_x <- log2(abs(Neftel$Neftel_OPC - Neftel$Neftel_NPC)+1)
Neftel$AC.MES_x <- log2(abs(Neftel$Neftel_AC - Neftel$Neftel_MES)+1)
#head(Neftel[10:ncol(Neftel)])

Neftel$X.axis <- NA
Neftel$X.axis[grep("OPC.NPC", Neftel$CellClass)] <- Neftel$OPC.NPC_x[grep("OPC.NPC", Neftel$CellClass)]
Neftel$X.axis[grep("AC.MES", Neftel$CellClass)] <- Neftel$AC.MES_x[grep("AC.MES", Neftel$CellClass)]
#head(Neftel[10:ncol(Neftel)])

#now change the sign based on maxclass for cell

#if AC is the max value, make postive; if MES is the max, make negative
#OPC = +, NPC = -

Neftel$X.axis_Class <- Neftel$X.axis * as.numeric(Neftel$ClassSign_X)
head(Neftel)

```


<table>
<caption>A data.frame: 6 × 31</caption>
<thead>
	<tr><th></th><th scope=col>CellType</th><th scope=col>PatientID</th><th scope=col>SampleID</th><th scope=col>Neftel_MES2</th><th scope=col>Neftel_MES1</th><th scope=col>Neftel_AC</th><th scope=col>Neftel_OPC</th><th scope=col>Neftel_NPC1</th><th scope=col>Neftel_NPC2</th><th scope=col>Neftel_G1.S</th><th scope=col>⋯</th><th scope=col>Neftel_NPC</th><th scope=col>Neftel_MES</th><th scope=col>Y.axis</th><th scope=col>CellClass</th><th scope=col>MaxClass</th><th scope=col>ClassSign_X</th><th scope=col>OPC.NPC_x</th><th scope=col>AC.MES_x</th><th scope=col>X.axis</th><th scope=col>X.axis_Class</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCACGGACAA</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td>-0.01913800</td><td>-0.071709632</td><td> 0.06980817</td><td>-0.03199840</td><td>-0.06716432</td><td> 0.095209179</td><td> 0.23451153</td><td>⋯</td><td> 0.0140224276</td><td>-0.04542382</td><td>-0.05578574</td><td>AC.MES </td><td>AC </td><td>-1</td><td>0.06491158</td><td>0.1573438</td><td>0.1573438</td><td>-0.1573438</td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCATCCGGGT</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td>-0.08061522</td><td> 0.033282546</td><td>-0.10055778</td><td>-0.12732257</td><td> 0.07976724</td><td>-0.005820044</td><td> 0.16263892</td><td>⋯</td><td> 0.0369735970</td><td>-0.02366634</td><td> 0.06063994</td><td>OPC.NPC</td><td>NPC</td><td>1 </td><td>0.21945809</td><td>0.1068728</td><td>0.2194581</td><td> 0.2194581</td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGGTACAGTTC</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td> 0.11935181</td><td>-0.007389301</td><td>-0.09355375</td><td> 0.13067802</td><td> 0.31634315</td><td> 0.165309543</td><td>-0.08324292</td><td>⋯</td><td> 0.2408263459</td><td> 0.05598125</td><td> 0.18484509</td><td>OPC.NPC</td><td>NPC</td><td>1 </td><td>0.15075244</td><td>0.2010504</td><td>0.1507524</td><td> 0.1507524</td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGTCTACGAGT</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td> 0.16918706</td><td> 0.091923570</td><td> 0.02553869</td><td> 0.09741310</td><td> 0.49856022</td><td> 0.477332226</td><td>-0.34891006</td><td>⋯</td><td> 0.4879462244</td><td> 0.13055531</td><td> 0.35739091</td><td>OPC.NPC</td><td>NPC</td><td>1 </td><td>0.47563811</td><td>0.1440681</td><td>0.4756381</td><td> 0.4756381</td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGAGTGGTAAT</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td> 0.14816768</td><td> 0.204905229</td><td>-0.02794300</td><td>-0.03860135</td><td> 0.05134542</td><td> 0.006812182</td><td> 0.62688342</td><td>⋯</td><td> 0.0290787990</td><td> 0.17653645</td><td>-0.14745765</td><td>AC.MES </td><td>MES</td><td>1 </td><td>0.09447952</td><td>0.2684098</td><td>0.2684098</td><td> 0.2684098</td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGCAGGACGTA</th><td>GSC</td><td>BT127</td><td>BT127_L</td><td> 0.37944621</td><td> 0.255718588</td><td> 0.46389484</td><td> 0.28277900</td><td> 0.04921320</td><td>-0.050914580</td><td>-0.22336134</td><td>⋯</td><td>-0.0008506884</td><td> 0.31758240</td><td>-0.18111584</td><td>AC.MES </td><td>AC </td><td>-1</td><td>0.36022906</td><td>0.1970003</td><td>0.1970003</td><td>-0.1970003</td></tr>
</tbody>
</table>




```R
saveRDS(Neftel, file = "Neftel_metadata.rds")
```

---
## 3.0 Plot results
---



```R
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


```


```R
colnames(Neftel)
Neftel$Diff_z <- scale(Neftel$Dev_IR_Diff)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CellType'</li><li>'PatientID'</li><li>'SampleID'</li><li>'Neftel_MES2'</li><li>'Neftel_MES1'</li><li>'Neftel_AC'</li><li>'Neftel_OPC'</li><li>'Neftel_NPC1'</li><li>'Neftel_NPC2'</li><li>'Neftel_G1.S'</li><li>'Neftel_G2.M'</li><li>'Developmental_GSC'</li><li>'InjuryResponse_GSC'</li><li>'Developmental_GSC_AUC'</li><li>'InjuryResponse_GSC_AUC'</li><li>'Developmental_GSC_AUC_z'</li><li>'InjuryResponse_GSC_AUC_z'</li><li>'Developmental_GSC_AUC_norm'</li><li>'InjuryResponse_GSC_AUC_norm'</li><li>'Dev_IR_Diff'</li><li>'Neftel_NPC'</li><li>'Neftel_MES'</li><li>'Y.axis'</li><li>'CellClass'</li><li>'MaxClass'</li><li>'ClassSign_X'</li><li>'OPC.NPC_x'</li><li>'AC.MES_x'</li><li>'X.axis'</li><li>'X.axis_Class'</li></ol>




```R
#### Cell state map with density TUMOUR CELLS

subset <- Neftel[Neftel$CellType == "Tumor", ]
subset$density <- get_density(subset$X.axis_Class, subset$Y.axis, n = 100)

b <- ggplot(subset, aes(x = X.axis_Class, y = Y.axis, color = CellType)) +
geom_point(aes(X.axis_Class, Y.axis, color = density), alpha = 1, size = 0.7) +  scale_color_viridis() + 
#xlab("<-- AC-like      MES-like -->") +
#ylab("<-- AC/MES-like      NPC/OPC-like -->") + 
geom_vline(xintercept = 0, lty = 2, col = "black") + 
geom_hline(yintercept = 0, lty = 2, col = "black") + 
#ggtitle("<-- OPC-like      NPC-like -->")  +
 theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))



#### Cell state map with Dev-IR TUMOUR CELLS
limit <- max(abs(subset$Dev_IR_Diff)) * c(-1, 1)
c <- ggplot(subset, aes(x = X.axis_Class, y = Y.axis, color = Dev_IR_Diff)) +
geom_point(alpha = 1, size = 0.7)  + 
#xlab("<-- AC-like      MES-like -->") +
#ylab("<-- AC/MES-like      NPC/OPC-like -->") + 
geom_vline(xintercept = 0, lty = 2, col = "black") + 
geom_hline(yintercept = 0, lty = 2, col = "black") + 
#ggtitle("<-- OPC-like      NPC-like -->")  +
 theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) +
scale_color_gradientn(colors = c("blue", "white", "red"), limit = limit)


print(b)
print(c)


```


![png](output_15_0.png)



![png](output_15_1.png)



```R
#### Cell state map with density GSCs

subset <- Neftel[Neftel$CellType == "GSC", ]
subset$density <- get_density(subset$X.axis_Class, 
                              subset$Y.axis, 
                              n = 100)

a <- ggplot(subset, aes(x = X.axis_Class, y = Y.axis, color = CellType)) +
geom_point(aes(X.axis_Class, Y.axis, color = density), alpha = 1, size = 0.7) +  scale_color_viridis() + 
#xlab("<-- AC-like      MES-like -->") +
#ylab("<-- AC/MES-like      NPC/OPC-like -->") + 
geom_vline(xintercept = 0, lty = 2, col = "black") + 
geom_hline(yintercept = 0, lty = 2, col = "black") + 
#ggtitle("<-- OPC-like      NPC-like -->")  +
 theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) 


print(a)


#### Cell state map with Dev-IR TUMOUR CELLS
limit <- max(abs(subset$Dev_IR_Diff)) * c(-1, 1)
d <- ggplot(subset, aes(x = X.axis_Class, y = Y.axis, color = Dev_IR_Diff)) +
geom_point(alpha = 1, size = 0.7)  + 
#xlab("<-- AC-like      MES-like -->") +
#ylab("<-- AC/MES-like      NPC/OPC-like -->") + 
geom_vline(xintercept = 0, lty = 2, col = "black") + 
#geom_hline(yintercept = 0, lty = 2, col = "black") + 
#ggtitle("<-- OPC-like      NPC-like -->")  +
 theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) +
scale_color_gradientn(colors = c("blue", "white", "red"), limit = limit)

print(d)
```


![png](output_16_0.png)



![png](output_16_1.png)



```R
pdf("~/Desktop/GSC_GBM_NeftelMapping.pdf", width = 10, height = 8)
ggarrange(a,b,d,c)
dev.off()
```


<strong>pdf:</strong> 2



```R
pdf("~/Desktop/GSC_GBM_NeftelMapping_noleg.pdf", width = 10, height = 8)
ggarrange(a,b,d,c, legend = "none")
dev.off()
```


<strong>pdf:</strong> 2



```R
pdf("~/Desktop/GSC_GBM_NeftelMapping_noleg_nolabels.pdf", width = 10, height = 8)
ggarrange(a,b,d,c, legend = "none")
dev.off()
```


<strong>pdf:</strong> 2

samples <- unique(Neftel$SampleID)

for (i in 1:length(samples)){
    
    
   print(samples[i])
    
    d <- ggplot(Neftel[Neftel$SampleID == samples[i], ], aes(x = X.axis_Class, y = Y.axis, color = CellType)) +
    geom_point(alpha = 0.1, size = 0.7) + scale_color_manual(values=c("black")) + xlab("<-- AC-like      MES-like -->") +
    ylab("<-- AC/MES-like      NPC/OPC-like -->") + geom_vline(xintercept = 0, lty = 2, col = "darkgrey") + 
    geom_hline(yintercept = 0, lty = 2, col = "darkgrey") + ggtitle("<-- OPC-like      NPC-like -->") +
    geom_rug(alpha = 0.1) + theme_classic() + xlim(c(-2,2)) + ylim(c(-2,2))
    print(d)
    
    
}

dev.off()
---
### Plot proportion cell types per sample
---


```R
sample.order <- as.character(read.table("~/Downloads/increasingproneural.txt", sep = "\t")[,1])
sample.order <- sample.order[!sample.order == ""]
sample.order <- unique(gsub('.{3}$', '', sample.order))
sample.order <- c(sample.order, unique(as.character(meta$Orig.ID[grep("_T", meta$Orig.ID)])))
sample.order <- sample.order[!sample.order == "G800_L"]
sample.order
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'G945-I_L'</li><li>'G837_L'</li><li>'G945-K_L'</li><li>'G566_L'</li><li>'G729_L'</li><li>'G945-J_L'</li><li>'G583_L'</li><li>'G564_L'</li><li>'G797_L'</li><li>'G549_L'</li><li>'G851_L'</li><li>'G876_L'</li><li>'G799_L'</li><li>'BT73_L'</li><li>'G885_L'</li><li>'G895_L'</li><li>'G523_L'</li><li>'G946-J_L'</li><li>'BT147_L'</li><li>'BT67_L'</li><li>'G620_L'</li><li>'G637_L'</li><li>'BT127_L'</li><li>'BT89_L'</li><li>'BT94_L'</li><li>'BT84_L'</li><li>'BT48_L'</li><li>'G946-K_L'</li><li>'G1003-A_T'</li><li>'G1003-B_T'</li><li>'G1003-C_T'</li><li>'G1003-D_T'</li><li>'G620_T'</li><li>'G910-A_T'</li><li>'G910-B_T'</li><li>'G910-C_T'</li><li>'G910-D_T'</li><li>'G910-E_T'</li><li>'G945-I_T'</li><li>'G945-J_T'</li><li>'G945-K_T'</li><li>'G946-I_T'</li><li>'G946-J_T'</li><li>'G946-K_T'</li><li>'G967-A_T'</li><li>'G967-B_T'</li><li>'G967-C_T'</li><li>'G967-D_T'</li><li>'G983-A_T'</li><li>'G983-B_T'</li><li>'G983-C_T'</li></ol>


sample.order <- c("G945-I_L",
                  "G837_L",
                  "G945-K_L",
                  "G566_L",
                  "G729_L",
                  "G945-J_L",
                  "G583_L",
                  "G564_L",
                  "G797_L",
                  "G549_L",
                  "G851_L",
                  "G876_L",
                  "G799_L",
                  "BT73_L",
                  "G885_L",
                  "G895_L",
                  "G523_L",
                  "G946-J_L",
                  "BT147_L",
                  "BT67_L",
                  "G620_L",
                  "G637_L",
                  "BT127_L",
                  "BT89_L",
                  "BT94_L",
                  "BT84_L",
                  "BT48_L",
                  "G1003_T",
                  "G620_T",
                  "G910_T",
                  "G945_T",
                  "G946_T",
                  "G967_T",
                  "G983_T"
                 )

```R
meta <- readRDS("GSC_GBM_scRNA_AddModuleScore_Neftel_meta.rds")
meta <- cbind(Neftel, meta$Orig.ID)
colnames(meta)[32] <- "Orig.ID"
```


```R
counts <-  table(meta$Orig.ID, meta$MaxClass)
prop <- prop.table(counts, margin = 1)
prop <- data.matrix(prop)

prop <- prop[sample.order, ]
head(prop)
```


              
                         AC          MES          NPC          OPC
      G945-I_L 0.0000000000 0.9904379422 0.0080321285 0.0015299292
      G837_L   0.0000000000 0.9398335822 0.0595263495 0.0006400683
      G945-K_L 0.0013583265 0.9163270850 0.0684596577 0.0138549307
      G566_L   0.0042016807 0.9957983193 0.0000000000 0.0000000000
      G729_L   0.0076362374 0.9909753558 0.0010413051 0.0003471017
      G945-J_L 0.0038543072 0.9668529582 0.0148390827 0.0144536520



```R
library(RColorBrewer)
```


```R
pdf("~/Desktop/Proportions_Suva_noG800_Sep2020.pdf", width = 35, height = 8)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow=c(2,2))

barplot(t(prop),
        las = 2,
        col = brewer.pal(n = 4, name = "Dark2"),
        main = "\n <-----Immes        Proneural------>",
        #ylab = "% cells",
        cex.names = 2,
        cex.axis = 2,
        #ylab = ""
        
       )
legend(65,1, c(colnames(prop)), 
       fill = brewer.pal(n = 4, name = "Dark2")
      )
dev.off()
```


<strong>pdf:</strong> 2



```R

```


```R



########


##plot BTSCs
df <- hybrid[hybrid$SampleType == "GBM", ]

#order BTSCs with increasing proneural score
#sample.order <- as.character(read.table("~/Downloads/increasingproneural.txt", sep = "\t")[,1])
#sample.order <- sample.order[!sample.order == ""]
#sample.order <- unique(gsub('.{3}$', '', sample.order))
#sample.order



counts <-  table(df$SampleID, df$FirstClass)
prop <- prop.table(counts, margin = 1)
prop <- data.matrix(prop)
rownames(prop) <- gsub("_GBM", "_T", rownames(prop))

#order <- sample.order[sample.order%in% rownames(prop)]
#prop <- prop[order, ]
head(prop)


#par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

barplot(t(prop),
        las = 2,
        col = c("gold", "Pink", "lightblue", "lightgreen"),
        main = "Tumours (not sorted in a smart way....)",
        ylab = "% cells"
       )
legend(30,1, c(colnames(prop)), 
       fill = c("gold", "Pink", "lightblue", "lightgreen")
      )

dev.off()
```


            
                       AC          MES          NPC          OPC
      G837_L 0.0000000000 0.9411137188 0.0586729251 0.0002133561
      G566_L 0.0042016807 0.9957983193 0.0000000000 0.0000000000
      G729_L 0.0072891357 0.9916695592 0.0006942034 0.0003471017
      G583_L 0.0125642490 0.7978298115 0.1444888635 0.0451170760
      G564_L 0.0536822228 0.9209603453 0.0240086323 0.0013487996
      G797_L 0.1094827586 0.8905172414 0.0000000000 0.0000000000



               
                         AC         MES         NPC         OPC
      G1003-A_T 0.864158830 0.127830024 0.002438175 0.005572971
      G1003-B_T 0.627906977 0.369509044 0.001291990 0.001291990
      G1003-C_T 0.571296296 0.428703704 0.000000000 0.000000000
      G1003-D_T 0.767783657 0.212228101 0.005291005 0.014697237
      G620_T    0.499145299 0.052991453 0.131623932 0.316239316
      G910-A_T  0.437956204 0.248175182 0.299270073 0.014598540



<strong>pdf:</strong> 2



```R
##plot BTSCs
df <- hybrid[hybrid$SampleType == "GBM", ]

#order BTSCs with increasing proneural score
#sample.order <- as.character(read.table("~/Downloads/increasingproneural.txt", sep = "\t")[,1])
#sample.order <- sample.order[!sample.order == ""]
#sample.order <- unique(gsub('.{3}$', '', sample.order))
#sample.order



counts <-  table(df$SampleID, df$FirstClass)
prop <- prop.table(counts, margin = 1)
prop <- data.matrix(prop)
rownames(prop) <- gsub("_GBM", "_T", rownames(prop))

#order <- sample.order[sample.order%in% rownames(prop)]
#prop <- prop[order, ]
head(prop)


par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

barplot(t(prop),
        las = 2,
        col = c("gold", "Pink", "lightblue", "lightgreen"),
        main = "Tumours (alphabetical)",
        ylab = "% cells"
       )
legend(30,1, c(colnames(prop)), 
       fill = c("gold", "Pink", "lightblue", "lightgreen")
      )
```


               
                         AC         MES         NPC         OPC
      G1003-A_T 0.864158830 0.127830024 0.002438175 0.005572971
      G1003-B_T 0.627906977 0.369509044 0.001291990 0.001291990
      G1003-C_T 0.571296296 0.428703704 0.000000000 0.000000000
      G1003-D_T 0.767783657 0.212228101 0.005291005 0.014697237
      G620_T    0.499145299 0.052991453 0.131623932 0.316239316
      G910-A_T  0.437956204 0.248175182 0.299270073 0.014598540



![png](output_30_1.png)


---
### Indexed Ridge Plot of Cell compositions. 
---


```R
library(ggridges)
```

    Warning message:
    “package ‘ggridges’ was built under R version 3.4.4”
    Attaching package: ‘ggridges’
    
    The following object is masked from ‘package:ggplot2’:
    
        scale_discrete_manual
    



```R
df  <- as.data.frame.matrix(prop) 
df$Sample <- rownames(df)
df$Index <- 1:nrow(df)
head(df)
```


<table>
<thead><tr><th></th><th scope=col>AC</th><th scope=col>MES</th><th scope=col>NPC</th><th scope=col>OPC</th><th scope=col>Sample</th><th scope=col>Index</th></tr></thead>
<tbody>
	<tr><th scope=row>G945-I_L</th><td>0.000000000 </td><td>0.9906292   </td><td>0.0084146108</td><td>0.0009562058</td><td>G945-I_L    </td><td>1           </td></tr>
	<tr><th scope=row>G837_L</th><td>0.000000000 </td><td>0.9417538   </td><td>0.0578195007</td><td>0.0004267122</td><td>G837_L      </td><td>2           </td></tr>
	<tr><th scope=row>G945-K_L</th><td>0.001358327 </td><td>0.9231187   </td><td>0.0687313230</td><td>0.0067916327</td><td>G945-K_L    </td><td>3           </td></tr>
	<tr><th scope=row>G566_L</th><td>0.004201681 </td><td>0.9957983   </td><td>0.0000000000</td><td>0.0000000000</td><td>G566_L      </td><td>4           </td></tr>
	<tr><th scope=row>G729_L</th><td>0.007636237 </td><td>0.9913225   </td><td>0.0006942034</td><td>0.0003471017</td><td>G729_L      </td><td>5           </td></tr>
	<tr><th scope=row>G945-J_L</th><td>0.004047023 </td><td>0.9695510   </td><td>0.0156099441</td><td>0.0107920601</td><td>G945-J_L    </td><td>6           </td></tr>
</tbody>
</table>




```R
ggplot(df, aes(x = AC, y = Index)) + geom_density_ridges()

```

    Picking joint bandwidth of 0.124
    ERROR while rich displaying an object: Error: geom_density_ridges requires the following missing aesthetics: y
    
    Traceback:
    1. FUN(X[[i]], ...)
    2. tryCatch(withCallingHandlers({
     .     rpr <- mime2repr[[mime]](obj)
     .     if (is.null(rpr)) 
     .         return(NULL)
     .     prepare_content(is.raw(rpr), rpr)
     . }, error = error_handler), error = outer_handler)
    3. tryCatchList(expr, classes, parentenv, handlers)
    4. tryCatchOne(expr, names, parentenv, handlers[[1L]])
    5. doTryCatch(return(expr), name, parentenv, handler)
    6. withCallingHandlers({
     .     rpr <- mime2repr[[mime]](obj)
     .     if (is.null(rpr)) 
     .         return(NULL)
     .     prepare_content(is.raw(rpr), rpr)
     . }, error = error_handler)
    7. mime2repr[[mime]](obj)
    8. repr_text.default(obj)
    9. paste(capture.output(print(obj)), collapse = "\n")
    10. capture.output(print(obj))
    11. evalVis(expr)
    12. withVisible(eval(expr, pf))
    13. eval(expr, pf)
    14. eval(expr, pf)
    15. print(obj)
    16. print.ggplot(obj)
    17. ggplot_build(x)
    18. ggplot_build.ggplot(x)
    19. by_layer(function(l, d) l$compute_geom_1(d))
    20. f(l = layers[[i]], d = data[[i]])
    21. l$compute_geom_1(d)
    22. f(..., self = self)
    23. check_required_aesthetics(self$geom$required_aes, c(names(data), 
      .     names(self$aes_params)), snake_class(self$geom))
    24. stop(name, " requires the following missing aesthetics: ", paste(missing_aes, 
      .     collapse = ", "), call. = FALSE)





![png](output_34_2.png)



```R
ggplot(df, aes(x = AC, y = Sample)) +
  geom_density_ridges(aes(fill = Sample)) 
```

    Picking joint bandwidth of NaN





![png](output_35_2.png)



```R
c("gold", "Pink", "lightblue", "lightgreen")
```


<table>
<thead><tr><th></th><th scope=col>AC</th><th scope=col>MES</th><th scope=col>NPC</th><th scope=col>OPC</th><th scope=col>Sample</th><th scope=col>Index</th></tr></thead>
<tbody>
	<tr><th scope=row>G945-I_L</th><td>0.000000000 </td><td>0.9906292   </td><td>0.0084146108</td><td>0.0009562058</td><td>G945-I_L    </td><td>1           </td></tr>
	<tr><th scope=row>G837_L</th><td>0.000000000 </td><td>0.9417538   </td><td>0.0578195007</td><td>0.0004267122</td><td>G837_L      </td><td>2           </td></tr>
	<tr><th scope=row>G945-K_L</th><td>0.001358327 </td><td>0.9231187   </td><td>0.0687313230</td><td>0.0067916327</td><td>G945-K_L    </td><td>3           </td></tr>
	<tr><th scope=row>G566_L</th><td>0.004201681 </td><td>0.9957983   </td><td>0.0000000000</td><td>0.0000000000</td><td>G566_L      </td><td>4           </td></tr>
	<tr><th scope=row>G729_L</th><td>0.007636237 </td><td>0.9913225   </td><td>0.0006942034</td><td>0.0003471017</td><td>G729_L      </td><td>5           </td></tr>
	<tr><th scope=row>G945-J_L</th><td>0.004047023 </td><td>0.9695510   </td><td>0.0156099441</td><td>0.0107920601</td><td>G945-J_L    </td><td>6           </td></tr>
</tbody>
</table>




```R
pdf("~/Desktop/proprs_Cells.pdf", width = 20, height = 10)

par(mfrow=c(4,1))

barplot(df$AC,
        names = rownames(df),
        col = "gold",
        las = 2, ylim=c(0,1)
       )
barplot(df$MES,
        names = rownames(df),
        col = "Pink",
        las = 2, ylim=c(0,1)
       )

barplot(df$NPC,
        names = rownames(df),
        col = "lightblue",
       las = 2, ylim=c(0,1))

barplot(df$OPC,
        names = rownames(df),
        col = "lightgreen",
        las = 2, ylim=c(0,1)
       )

dev.off()
```


<strong>pdf:</strong> 2



```R
plot(density(df$AC))
plot(density(df$MES))
```


![png](output_38_0.png)



![png](output_38_1.png)



```R

```


```R

```


```R

```


```R

```


```R

```

---
### Calculate distance between points to the core of the map (0,0)
---

We want to show that GSCs localize to the centre of the map (more so than tumour cells)

-- Make a distance to the core metric


```R
Neftel$Core_dist <- NA
head(Neftel)
```


<table>
<thead><tr><th></th><th scope=col>nGene</th><th scope=col>nUMI</th><th scope=col>percent.mito</th><th scope=col>PatientID</th><th scope=col>SampleID</th><th scope=col>SampleType</th><th scope=col>Neftel_MES2</th><th scope=col>Neftel_MES1</th><th scope=col>Neftel_AC</th><th scope=col>Neftel_OPC</th><th scope=col>⋯</th><th scope=col>MES</th><th scope=col>Y.axis</th><th scope=col>CellClass</th><th scope=col>MaxClass</th><th scope=col>ClassSign_X</th><th scope=col>OPC.NPC_x</th><th scope=col>AC.MES_x</th><th scope=col>X.axis</th><th scope=col>X.axis_Class</th><th scope=col>Core_dist</th></tr></thead>
<tbody>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCACGGACAA</th><td> 640       </td><td>  875      </td><td>0.043428571</td><td>BT127      </td><td>BT127_L    </td><td>BTSC       </td><td>-0.02128193</td><td>-0.06049909</td><td> 0.07960276</td><td>-0.03909653</td><td>⋯          </td><td>-0.04089051</td><td>-0.05471932</td><td>AC.MES     </td><td>AC         </td><td>-1         </td><td>0.08947099 </td><td>0.1641340  </td><td>0.1641340  </td><td>-0.1641340 </td><td>NA         </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCATCCGGGT</th><td>1036       </td><td> 2408      </td><td>0.002076412</td><td>BT127      </td><td>BT127_L    </td><td>BTSC       </td><td>-0.07451009</td><td> 0.02618286</td><td>-0.09767411</td><td>-0.13332228</td><td>⋯          </td><td>-0.02416361</td><td> 0.06967643</td><td>OPC.NPC    </td><td>NPC        </td><td>1          </td><td>0.23736193 </td><td>0.1023363  </td><td>0.2373619  </td><td> 0.2373619 </td><td>NA         </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGGTACAGTTC</th><td>3240       </td><td>10058      </td><td>0.078047326</td><td>BT127      </td><td>BT127_L    </td><td>BTSC       </td><td> 0.12474159</td><td> 0.01255947</td><td>-0.09031973</td><td> 0.10493229</td><td>⋯          </td><td> 0.06865053</td><td> 0.16629791</td><td>OPC.NPC    </td><td>NPC        </td><td>1          </td><td>0.17634340 </td><td>0.2128436  </td><td>0.1763434  </td><td> 0.1763434 </td><td>NA         </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGTCTACGAGT</th><td>3337       </td><td>10798      </td><td>0.061863308</td><td>BT127      </td><td>BT127_L    </td><td>BTSC       </td><td> 0.16997627</td><td> 0.10159148</td><td> 0.03591002</td><td> 0.07576596</td><td>⋯          </td><td> 0.13578387</td><td> 0.34807302</td><td>OPC.NPC    </td><td>NPC        </td><td>1          </td><td>0.49374051 </td><td>0.1373381  </td><td>0.4937405  </td><td> 0.4937405 </td><td>NA         </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGAGTGGTAAT</th><td>4140       </td><td>14601      </td><td>0.081501267</td><td>BT127      </td><td>BT127_L    </td><td>BTSC       </td><td> 0.15017229</td><td> 0.22260042</td><td>-0.03413590</td><td>-0.07153851</td><td>⋯          </td><td> 0.18638636</td><td>-0.16201549</td><td>AC.MES     </td><td>MES        </td><td>1          </td><td>0.13212852 </td><td>0.2874986  </td><td>0.2874986  </td><td> 0.2874986 </td><td>NA         </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGCAGGACGTA</th><td> 543       </td><td>  820      </td><td>0.108536585</td><td>BT127      </td><td>BT127_L    </td><td>BTSC       </td><td> 0.38546143</td><td> 0.24715598</td><td> 0.45306749</td><td> 0.27539241</td><td>⋯          </td><td> 0.31630870</td><td>-0.17767508</td><td>AC.MES     </td><td>AC         </td><td>-1         </td><td>0.35390705 </td><td>0.1849262  </td><td>0.1849262  </td><td>-0.1849262 </td><td>NA         </td></tr>
</tbody>
</table>




```R
for (i in 1:nrow(Neftel)){
    
    core <- c(0,0)
    cell <- c(Neftel$X.axis_Class[i], Neftel$Y.axis[i])

    x <- (core[1] - cell[1])^2
    y <- (core[2] - cell[2])^2

    Neftel$Core_dist[i] <- sqrt(x+y)
    
}


```


```R
head(Neftel)
```


<table>
<thead><tr><th></th><th scope=col>nGene</th><th scope=col>nUMI</th><th scope=col>percent.mito</th><th scope=col>PatientID</th><th scope=col>SampleID</th><th scope=col>SampleType</th><th scope=col>Neftel_MES2</th><th scope=col>Neftel_MES1</th><th scope=col>Neftel_AC</th><th scope=col>Neftel_OPC</th><th scope=col>⋯</th><th scope=col>MES</th><th scope=col>Y.axis</th><th scope=col>CellClass</th><th scope=col>MaxClass</th><th scope=col>ClassSign_X</th><th scope=col>OPC.NPC_x</th><th scope=col>AC.MES_x</th><th scope=col>X.axis</th><th scope=col>X.axis_Class</th><th scope=col>Core_dist</th></tr></thead>
<tbody>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCACGGACAA</th><td> 640       </td><td>  875      </td><td>0.043428571</td><td>BT127      </td><td>BT127_L    </td><td>BTSC       </td><td>-0.02128193</td><td>-0.06049909</td><td> 0.07960276</td><td>-0.03909653</td><td>⋯          </td><td>-0.04089051</td><td>-0.05471932</td><td>AC.MES     </td><td>AC         </td><td>-1         </td><td>0.08947099 </td><td>0.1641340  </td><td>0.1641340  </td><td>-0.1641340 </td><td>0.1730149  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCATCCGGGT</th><td>1036       </td><td> 2408      </td><td>0.002076412</td><td>BT127      </td><td>BT127_L    </td><td>BTSC       </td><td>-0.07451009</td><td> 0.02618286</td><td>-0.09767411</td><td>-0.13332228</td><td>⋯          </td><td>-0.02416361</td><td> 0.06967643</td><td>OPC.NPC    </td><td>NPC        </td><td>1          </td><td>0.23736193 </td><td>0.1023363  </td><td>0.2373619  </td><td> 0.2373619 </td><td>0.2473772  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGGTACAGTTC</th><td>3240       </td><td>10058      </td><td>0.078047326</td><td>BT127      </td><td>BT127_L    </td><td>BTSC       </td><td> 0.12474159</td><td> 0.01255947</td><td>-0.09031973</td><td> 0.10493229</td><td>⋯          </td><td> 0.06865053</td><td> 0.16629791</td><td>OPC.NPC    </td><td>NPC        </td><td>1          </td><td>0.17634340 </td><td>0.2128436  </td><td>0.1763434  </td><td> 0.1763434 </td><td>0.2423881  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGTCTACGAGT</th><td>3337       </td><td>10798      </td><td>0.061863308</td><td>BT127      </td><td>BT127_L    </td><td>BTSC       </td><td> 0.16997627</td><td> 0.10159148</td><td> 0.03591002</td><td> 0.07576596</td><td>⋯          </td><td> 0.13578387</td><td> 0.34807302</td><td>OPC.NPC    </td><td>NPC        </td><td>1          </td><td>0.49374051 </td><td>0.1373381  </td><td>0.4937405  </td><td> 0.4937405 </td><td>0.6040981  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGAGTGGTAAT</th><td>4140       </td><td>14601      </td><td>0.081501267</td><td>BT127      </td><td>BT127_L    </td><td>BTSC       </td><td> 0.15017229</td><td> 0.22260042</td><td>-0.03413590</td><td>-0.07153851</td><td>⋯          </td><td> 0.18638636</td><td>-0.16201549</td><td>AC.MES     </td><td>MES        </td><td>1          </td><td>0.13212852 </td><td>0.2874986  </td><td>0.2874986  </td><td> 0.2874986 </td><td>0.3300068  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGCAGGACGTA</th><td> 543       </td><td>  820      </td><td>0.108536585</td><td>BT127      </td><td>BT127_L    </td><td>BTSC       </td><td> 0.38546143</td><td> 0.24715598</td><td> 0.45306749</td><td> 0.27539241</td><td>⋯          </td><td> 0.31630870</td><td>-0.17767508</td><td>AC.MES     </td><td>AC         </td><td>-1         </td><td>0.35390705 </td><td>0.1849262  </td><td>0.1849262  </td><td>-0.1849262 </td><td>0.2564491  </td></tr>
</tbody>
</table>




```R
library(ggpubr)
```

    Warning message:
    “package ‘ggpubr’ was built under R version 3.4.4”Loading required package: magrittr



```R
p <- ggplot(Neftel, aes(x=SampleType, y=Core_dist, fill = SampleType)) + 
  geom_boxplot(outlier.colour=alpha("black", 0.8), outlier.size=0.1) + theme_minimal() +
 scale_fill_brewer(palette="RdBu")

pdf("~/Desktop/Distance_Core.pdf", width = 3, height = 3)
p
dev.off()
```




<strong>pdf:</strong> 2



```R
p <- ggplot(Neftel, aes(x=SampleType, y=Core_dist, fill = SampleType)) + geom_violin() + 
  geom_boxplot(outlier.colour=alpha("black", 0.8), outlier.size=0.1, width = 0.025) + theme_minimal() +
 scale_fill_brewer(palette="RdBu")
p
```




![png](output_50_1.png)



```R

```
