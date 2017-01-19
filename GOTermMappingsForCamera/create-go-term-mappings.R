# R versions of the MSigDB data sets are available here: 
# http://bioinf.wehi.edu.au/software/MSigDB/

# The data are provided as lists, where each element is a named gene set. Within
# the gene set are Entrez IDs for the genes. 

# For the GO terms, it's a little annoying to map from the titles given to the
# gene sets back to the GO ID and to the ontology (biological process, molecular
# function, or cellular component) it belongs to. 

# This script takes the GO gene sets from that webpage and splits them into 
# three separate gene sets, by ontology. It also creates a data frame that
# can be used to map back to GO ID and so on. 

library(GO.db)
library(dplyr)
library(stringr)

load("data/human_c5_v5p2.rdata")

# get all of the current GO terms and turn into data frame
terms <- as.list(GOTERM)
terms <- data.frame(Term = sapply(terms, function(x) x@Term), 
   Ontology = sapply(terms, function(x) x@Ontology), 
   GOId = sapply(terms, function(x) x@GOID))

# terms from camera results and from GO.db have different punctuation
# making matching complicated, so remove everything that's not 
# alphanumeric to match them up
database_terms <- terms %>% 
   mutate(Term = as.character(Term)) %>% 
   mutate(GOTermMatch = str_replace_all(Term, "[^[:alnum:]]", ""), 
      GOTermMatch = str_to_lower(GOTermMatch)) 

msigdb_terms <- data.frame(TermFromMSigDB = names(Hs.c5)) %>% 
   mutate(GOTermMatch = as.character(TermFromMSigDB)) %>% 
   mutate(GOTermMatch = str_replace(GOTermMatch, "^GO", ""),
      GOTermMatch = str_replace_all(GOTermMatch, "[^[:alnum:]]", ""), 
      GOTermMatch = str_to_lower(GOTermMatch)) 

# merge the term information from the GO database with the term names from 
# msigdb
matched <- merge(database_terms, msigdb_terms, 
   by = "GOTermMatch", all.y = TRUE) %>% 
   # a few terms may not match the database
   dplyr::mutate(Term = ifelse(is.na(Term), 
      as.character(TermFromMSigDB), as.character(Term)))

# Missing ----------------------------------------------------------------------
# How many of the terms in the Hs.c5 data are missing from the GO package?
nrow(filter(matched, str_detect(Term, "^GO_")))
filter(matched, str_detect(Term, "^GO_")) %>% View()
# I looked some of the missing terms up and they've all been replaced by or are
# synonyms for other terms.

# Duplicates -------------------------------------------------------------------
# How many of the terms end up as duplicates?
duplicate_matches <- filter(matched, duplicated(GOTermMatch)) %>% 
   .$GOTermMatch
sum(matched$GOTermMatch %in% duplicate_matches)
filter(matched, GOTermMatch %in% duplicate_matches)
# This happens when the term names differ only by non-alphanumeric characters, 
# so the name used to match matches more than one term. I looked them up. 
# Note that you have to specify human only

# length(Hs.c5[["GO_ALCOHOL_DEHYDROGENASE_NADP_ACTIVITY"]])
#> 16
# GO:0018455 has 1 gene
# GO:0008106 has 25 genes

# length(Hs.c5[["GO_NAD_BINDING"]])
# > 53
# GO:0070403 has 14 genes
# GO:0051287 has 56 genes

# length(Hs.c5[["GO_NADP_BINDING"]])
# > 44
# GO:0050661 has 50 genes
# GO:0070401 has 1 gene

# so remove GO:0018455, GO:0070403, and GO:0070401 
matched <- matched %>% 
   filter(!(GOId %in% c("GO:0018455", "GO:0070403", "GO:0070401")))

# Save files -------------------------------------------------------------------
Hs.c5.BP <- Hs.c5[names(Hs.c5) %in% filter(matched, Ontology == "BP")$TermFromMSigDB]
Hs.c5.MF <- Hs.c5[names(Hs.c5) %in% filter(matched, Ontology == "MF")$TermFromMSigDB]
Hs.c5.CC <- Hs.c5[names(Hs.c5) %in% filter(matched, Ontology == "CC")$TermFromMSigDB]

save(Hs.c5.BP, file = "GOTermMappingsForCamera/Hs.c5.BP.rdata", 
   compress = "bzip2", compression_level = 9)
save(Hs.c5.MF, file = "GOTermMappingsForCamera/Hs.c5.MF.rdata", 
   compress = "bzip2", compression_level = 9)
save(Hs.c5.CC, file = "GOTermMappingsForCamera/Hs.c5.CC.rdata", 
   compress = "bzip2", compression_level = 9)
write.csv(matched, file = "GOTermMappingsForCamera/GOTermMappings.csv", 
   row.names = FALSE)