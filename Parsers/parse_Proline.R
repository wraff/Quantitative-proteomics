library(stringr)
library(readxl)

# Parameters set in launcher script (see "Lauchers" section)
# quantif_file = quantification file to parse
# output_file = output table
# intensities_type = abundance/raw_abundance
# sheet = sheet number if quantif_file is a xlsx file   # (wr9sep21) set to NULL or NA for automatic determination (see below)

## Load data
if(grepl(".xlsx",quantif_file)){
  ## (wr9sep21) you can automatically deduce the sheet number if the name 'Protein sets' is always there (no longer need to specify 'sheet' as argument)
  if(length(sheet) >0) if(any(is.na(sheet))) sheet <- NULL
  if(length(sheet) <1) sheets <- readxl::excel_sheets(quantif_file)
  sheet <- which(sheets == "Protein sets") 
  proteinGroupsInput <- read_excel(quantif_file, sheet, col_names = TRUE)
} else {
  proteinGroupsInput <- read.table(quantif_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
}


## Extract data
Id <- proteinGroupsInput[grepl("id$",names(proteinGroupsInput))][1]
colnames(Id) <- "Id"

Accession <- proteinGroupsInput[,which(colnames(proteinGroupsInput) %in% "accession")][1]    # (wr9sep21 change) 
colnames(Accession) <- "Accession"

Gene_name <- proteinGroupsInput[grepl("gene_name",names(proteinGroupsInput))]
if(length(Gene_name) >0) {
  colnames(Gene_name) <- "Gene_name"
}

## (wr9sep21)  initial Fasta-line from 'description'
chDesc <- any(colnames(proteinGroupsInput) %in% "description") 
if(any(chNpep)) { Description <- proteinGroupsInput[,which(colnames(proteinGroupsInput) %in% "description")][1]
} else stop("Could not find column entitled 'description' in ",quantif_file," !!" )

## (wr9sep21)  total number of peptides
chNpep <- any(colnames(proteinGroupsInput) %in% "#peptides") 
if(any(chNpep)) { Npeptides <- proteinGroupsInput[,which(chNpep)[1]] 
} else { Npeptides <- NULL; message(" .. Cound not find column '#peptides'")}


intensities <- proteinGroupsInput[grepl(paste0("^",intensities_type,".+"),names(proteinGroupsInput))]
colnames(intensities) <- paste0("Intensity",sub(intensities_type, "", colnames(intensities)))

identification_types <- proteinGroupsInput[grepl("psm_count",names(proteinGroupsInput))]
samples = sub("psm_count_", "", colnames(identification_types))
for(sample in samples){
  identification_type = identification_types[grepl(paste0("_",sample),names(identification_types))]
  colname = colnames(identification_type)[1]
  identification_type[[colname]][identification_type[[colname]] >0 & !(is.na(identification_type[[colname]]))] <- "By MS/MS"
  identification_type[[colname]][identification_type[[colname]]==0 & !(is.na(identification_type[[colname]]))] <- "By matching"
  identification_type[[colname]][is.na(identification_type[[colname]])] <- NA
  identification_types[grepl(paste0("_",sample),names(identification_types))] = identification_type[[colname]]
}
colnames(identification_types) <- paste0("Identification_type_",samples)

specific_peptides = proteinGroupsInput[grepl("specific_peptide_matches",names(proteinGroupsInput))]
if(length(specific_peptides)==0){
  specific_peptides = proteinGroupsInput[grepl("specific_spectral_count",names(proteinGroupsInput))]
}

if(length(specific_peptides) >0){
  colnames(specific_peptides) <- "Specific_peptides"
}

## Build dataframe
proteinGroupsOutput = cbind(Id,Accession,Description,Gene_name,identification_types,Npeptides,specific_peptides,intensities)  # (wr9sep21 update)
write.table(proteinGroupsOutput,output_file,row.names=FALSE,sep="\t")
