###############################################################################
### Create Universal Metadata for SU2C Project ###

### OVERVIEW: Organize metadata for all SU2C project samples
###           for 'raw' data, see their descriptions in associated README.txt
###           files in their directory location
print('BEGIN SCRIPT')
print(Sys.time())
### Load required packages + functions ----------------------------------------
require(readxl)
### Load in Data
# data_folder <- '../data_01_2018'
## setup directories
top_dir <- '~/projects/su2c_v2'
## raw data directories
raw_RNA_dir <- file.path(top_dir, 'data/raw/RNA')
raw_coh123_nonSU2C <- file.path(raw_RNA_dir, 
                                'RNA_florence_06_2018_SU2C_coh123_nonSU2C')
raw_coh123 <- file.path(raw_RNA_dir,
                        'RNA_florence_01_2018_SU2C_coh123')
raw_coh4 <- file.path(raw_RNA_dir,
                      'RNA_florence_07_2018_SU2C_cohort4')
## directory to output preprocessed RNA data
preprocess_metadata_dir <- file.path(top_dir, 'data/preprocessed/metadata')
if (!dir.exists(preprocess_metadata_dir))
{
  dir.create(preprocess_metadata_dir)
}
## sample metadata, cohorts 1-3
sample_meta_file <- 'M_info_cohort_1_2_3_nonSU2C_init.RData'
## sample metadata, 'cohort' 4
RNA_count_meta_cohort_4_file <- 'M_count_info_cohort4_GBM.RData'
## patient metadata annotations
# patient_info_file <- 'SU2C_sample_annotations_formatted_Oct2017.xlsx'
patient_info_file <- 'M_samples_info_cohort1_2_3_nonSU2C_91samples.RData'
load(file.path(raw_coh123_nonSU2C, sample_meta_file))
# M_samples <- readxl::read_xlsx(file.path(raw_coh123, patient_info_file))
load(file.path(raw_coh123, patient_info_file))
load(file.path(raw_coh4, RNA_count_meta_cohort_4_file))
## convert M_samples to dataframe
# M_samples_temp <- M_samples
M_samples <- as.data.frame(M_samples)
# if (!all(M_samples == M_samples_temp))
# {
#   stop('problem in conversion of tibble to dataframe')
# }
# rm(M_samples_temp)


dim(M_info_cohort_1_2_3_nonSU2C_init)
# [1] 189  12
dim(M_samples)
# [1] 91 12
# dim(M_count_cohort4_GBM)
# # [1] 26475     4
dim(M_info_cohort4_GBM)
# [1]  4 13

### Combine metadata for cohorts 1-3 + nonSU2C + cohort 4 ---------------------
common_metadata <- intersect(colnames(M_info_cohort_1_2_3_nonSU2C_init),
                             colnames(M_info_cohort4_GBM))
M_info_coh1234_nonSU2C <- rbind(M_info_cohort_1_2_3_nonSU2C_init[, common_metadata],
                                M_info_cohort4_GBM[, common_metadata])
M_info_coh1234_nonSU2C <- as.data.frame(M_info_coh1234_nonSU2C, stringsAsFactors = F)
# colnames(M_info_coh1234_nonSU2C)
# [1] "Sample"          "Lib"             "Lab"             "Patient"         "Type_init"      
# [6] "Type_bis"        "Patient_used_id" "SampleID_used"   "Tumor_type"      "Cohort"         
# [11] "Type"            "condition"     

### Assign patient age, age group, primary recurrent status to samples --------
## for patient metadata, patient id is unique identifier
# any(duplicated(M_samples$Patient_ID))
# [1] FALSE
# colnames(M_samples)
# [1] "Cohort"                                   "Patient_ID"                              
# [3] "Sex"                                      "Age"                                     
# [5] "Age_group"                                "Tumour"                                  
# [7] "Tumour_Type"                              "Primary_recurrent"                       
# [9] "Primary_for_recurrent"                    "Sphere_forming_frequency_percentage"     
# [11] "Population_doubling_time_hrs"             "Xenograft_formation_median_survival_days"
## subset only for samples corresponding to patient in patient metadata
rownames(M_samples) <- M_samples$Patient_ID
patient.ids.coh1234_nonSU2C <- M_info_coh1234_nonSU2C$Patient_used_id
head(patient.ids.coh1234_nonSU2C)
# [1] "G549" "G583" "PFA1" "PFA6" "G729" "PFA2"
all(patient.ids.coh1234_nonSU2C %in% M_samples$Patient_ID)
patient.ids.coh1234_nonSU2C[which(!patient.ids.coh1234_nonSU2C %in%
                                    M_samples$Patient_ID)]
# [1] "BT1998001" "BT2010114" "TL3Adult"  "TL1Adult"  "G938"      "G938"      "G938" 

M_info_coh1234_nonSU2C <- M_info_coh1234_nonSU2C[patient.ids.coh1234_nonSU2C 
                                                 %in% M_samples$Patient_ID,]
## assign patient age
M_info_coh1234_nonSU2C$Age <- M_samples[M_info_coh1234_nonSU2C$Patient_used_id,
                                        'Age']
## convert age to numeric type if string denotes an integer
## convert to NA otherwise. 
integer.age <- grep('^[0-9]*$', M_info_coh1234_nonSU2C$Age)
temp.age <- M_info_coh1234_nonSU2C$Age
M_info_coh1234_nonSU2C$Age <- as.numeric(M_info_coh1234_nonSU2C$Age)
## Check that this is done correctly
if (!all(M_info_coh1234_nonSU2C$Age[integer.age] == as.numeric(temp.age[integer.age]))) {
  stop('expect any age that was denoted by an integer value in a string was converted
       to proper value upon coercion to numeric type')
}
# M_info_coh1234_nonSU2C[-integer.age,c('SampleID_used', 'Age')]
# A61501_Taylor_1760_Primary          PFA6_T   NA
# A61503_Taylor_1705_Primary          PFA2_T 0.75
# A61521_Control_NS2               HF7450_NL   NA
# A61522_Taylor_1705_Primary          PFA2_T 0.75
# A61527_Taylor_1756_Primary          PFA5_T   NA
# A61530_Taylor_1705_Line             PFA2_L 0.75
# A61541_Control_NS3               HF5205_NL   NA
# A61545_Control_NS1               HF6562_NL   NA
# A61549_Taylor_1756_Line             PFA5_L   NA
# A61552_Taylor_1760_Line             PFA6_L   NA
# A61553_Taylor_1705_Primary          PFA2_T 0.75
# A67957_PFA3_tissue                  PFA3_T   NA
# A67970_Control_NS1               HF6562_NL   NA
# A67973_PFA5_normoxia       PFA5_L_normoxia   NA
# A67983_PFA3_line                    PFA3_L   NA
# A67988_Control_NS3               HF5205_NL   NA
# A67990_Control_NS2               HF7450_NL   NA
# A68002_PFA2_normoxia       PFA2_L_normoxia 0.75
# A68008_PFA2_hypoxia         PFA2_L_hypoxia 0.75
# A68009_PFA5_hypoxia         PFA5_L_hypoxia   NA
# A77699_HF7450_line               HF7450_NL   NA
# A77704_HF6_line                     HF6_NL   NA
# A77712_HF5205_line               HF5205_NL   NA
# A77722_HF5_line                     HF5_NL   NA
# A77723_HF6562_line               HF6562_NL   NA
# A77724_HF11_line                   HF11_NL   NA

## assign age group
M_info_coh1234_nonSU2C$Age_group <- M_samples[M_info_coh1234_nonSU2C$Patient_used_id,
                                              'Age_group']
# unique(M_info_coh1234_nonSU2C$Age_group)
# [1] "adult"     "pediatric" "fetal"  
## assign primary recurrent status + primary for recurrent
## primary for recurrent = primary tumor associated with recurrent tumor
M_info_coh1234_nonSU2C$Primary_recurrent <- M_samples[M_info_coh1234_nonSU2C$Patient_used_id,
                                                      'Primary_recurrent']
M_info_coh1234_nonSU2C$Primary_for_recurrent <- M_samples[M_info_coh1234_nonSU2C$Patient_used_id,
                                                      'Primary_for_recurrent']
## assign sex
M_info_coh1234_nonSU2C$Sex <- M_samples[M_info_coh1234_nonSU2C$Patient_used_id,
                                                          'Sex']
### Save Metadata -------------------------------------------------------------
sample_meta_coh1234_nonSU2C <- M_info_coh1234_nonSU2C
save(sample_meta_coh1234_nonSU2C, file = file.path(preprocess_metadata_dir,
                                                   'sample_meta_coh1234_nonSU2C.RData'))
print('END SCRIPT')
print(Sys.time())
sessionInfo()
