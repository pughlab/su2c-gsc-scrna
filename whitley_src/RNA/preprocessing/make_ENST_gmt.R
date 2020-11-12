###############################################################################
### Make GMT File For DNAm probes

## We shall make a gmt file mapping pathways to ENST ids for GRCh37

library(biomaRt)
# get grch37 assembly
listEnsemblArchives()
listMarts(host = 'http://grch37.ensembl.org')
ensemblGRCH37 <- useMart(host='grch37.ensembl.org', 
                         biomart='ENSEMBL_MART_ENSEMBL', 
                         dataset='hsapiens_gene_ensembl')
attr <- listAttributes(ensemblGRCH37)
attr$name[grep('^go|name_1006', attr$name)]
# get annotations mapping ENST ids to go ids
output_dir <- file.path('../../../data/preprocessed/gmt_files')
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
print('Retrieving following attributes')
attr_retrieve <- c('ensembl_transcript_id', 'go_id', 'go_linkage_type', 'name_1006', 'namespace_1003')
attr[attr$name %in% attr_retrieve, ]
print('submitting query to biomart')
anno_GRCH37 <- getBM(attributes = attr_retrieve, 
                     mart = ensemblGRCH37)

# function to make gmt file
make_gmt <- function(biomart_anno = anno_GRCH37,
                     output_file = file.path(output_dir, 'biomart_GRCh37_ENST_GO_no_iea.gmt')) {
  # take biomart annotations of ENST ids to GO terms and make gmt file
  
  # get rid of all GO annotations that are IEA, RCA, ND
  biomart_anno <- biomart_anno[which(!biomart_anno$go_linkage_type %in% c('IEA', 'RCA', 'ND')), ]
  # get rid of all empty GO annotations
  biomart_anno <- biomart_anno[which(!biomart_anno$go_id == ''),]
  
  # steps for constructing file
  # remove output file if already exists
  if (file.exists(output_file)) {
    system(paste('rm', output_file))
  }
  print(paste('writing gmt file', output_file))
  output_fileconn <- file(output_file)
  file_str <- ''
  # tally how many pathways we've kept
  num_pathways <- 0
  # loop through go ids and write file
  unique_go_ids <- unique(biomart_anno$go_id)
  n_go_ids <- length(unique_go_ids)
  pb <- txtProgressBar(min = 0, max = n_go_ids, style = 3)
  for (i in 1:n_go_ids) {
    go_id_i <- unique_go_ids[i]
    row_ind <- biomart_anno$go_id == go_id_i
    ensembl_transcript_ids_i <- unique(biomart_anno$ensembl_transcript_id[row_ind])
    if (length(ensembl_transcript_ids_i) > 0) {
      # add go term + probes to file
      go_type <- unique(biomart_anno$namespace_1003[biomart_anno$go_id == go_id_i])
      go_term <- unique(biomart_anno$name_1006[biomart_anno$go_id == go_id_i])
      line_write <- paste(paste(paste(toupper(go_term), 
                                      toupper(go_type), 
                                      go_id_i, sep = '%'),
                                go_term,
                                sep = '\t'), 
                          paste(ensembl_transcript_ids_i, collapse = '\t'), 
                          sep = '\t')
      line_write <- paste0(line_write, '\n')
      file_str <- paste0(file_str, line_write)
      num_pathways <- num_pathways + 1
    }
    setTxtProgressBar(pb, i)
  }
  writeLines(file_str, con = output_fileconn)
  print(paste('Finished writing gmt file with', num_pathways, 'pathways'))
  close(output_fileconn)
}

# make gmt file for all probes
make_gmt()

Sys.time()
sessionInfo()
