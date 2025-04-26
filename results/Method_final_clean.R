
library(DspikeIn)
library(ggplot2)
library(phyloseq)
library(decontam)
library(ggrepel)

################## Prok correct 12 Jan 2025

## (Tab.1_Method_metadata_18March_2025_blanks_refined_noChelex.csv)


#Prokaryotes
setwd("/Users/mghotbi/Desktop/Publications/Second Chapter_96-well plate_1.5Timeseries /2A_method_HPC/2A_method")

OTU<-read.csv("ASV_table_for_correct_prok.csv",header=T, row.names = 1)
#TAX<-read.csv("gg2_taxonomy_for_silva_selected_codes_without_mit_chl_Unassigned_prok_correct_parsed.csv",header=T, row.names = 1)
TAX<-read.csv("gg2-taxonomy_silvaselectedcodes_No-mit-chlo-Unassigned_prok_correct_parsed_tidy_23March2025.csv",header=T, row.names = 1)
#TAX<-read.csv("Silva_taxonomy_No_mit_no_chloroplast_Unassigned_prok_correct_parsed31March_fromDavidfile.csv",header=T, row.names = 1)

#Meta <-read.csv("Method_metadata_with_all_blanks_4_equal_replicates_No_10X_noBeadclean_Final_12Jan2025.csv", header=T, row.names = 1)
Meta <-read.csv("Tab.1_Method_metadata_18March_2025_blanks_refined_noChelex.csv", header=T, row.names = 1)
#Meta <-read.csv("Method_metadata_for_rare_at5000_67samples.csv", header=T, row.names = 1)


TAX
Meta

#####################################################

#chloroplast 

setwd("/Users/mghotbi/Desktop/Publications/Second Chapter_96-well plate_1.5Timeseries /2A_method_HPC/method.gpt.counterparts-full-rep-chloroplast-table.qza.biom")
#Chloroplast

OTU<-read.csv("Phytoref_ASV_table_GPT_counterparts.csv",header=T, row.names = 1)
dim(OTU)
#TAX<-read.csv("TAX_Phytoref_GPT_counterparts.csv",header=T, row.names = 1)
TAX<-read.csv("TAX_phytoref_7cat_to_fool_22jan25.csv",header=T, row.names = 1)

TAX
#Meta <-read.csv("~/Desktop/Publications/Second Chapter_96-well plate_1.5Timeseries /2A_method_HPC/method.gpt.counterparts-full-rep.qza.biom/Mtadata_fullrep_GPT_counterpart_corrected_allmethods.csv",header=T, row.names = 1)
#Meta <-read.csv("~/Desktop/Publications/Second Chapter_96-well plate_1.5Timeseries /2A_method_HPC/method.gpt.counterparts-full-rep.qza.biom/Mtadata_fullrep_GPT_counterpart_corrected_allmethods.csv",header=T, row.names = 1)
#Meta <-read.csv("/Users/mghotbi/Desktop/Publications/Second Chapter_96-well plate_1.5Timeseries /2A_method_HPC/2A_method/Method_metadata_with_all_blanks_4_equal_replicates_No_10X_noBeadclean_Final_12Jan2025.csv",header=T, row.names = 1)
Meta <-read.csv("Tab.1_Method_metadata_18March_2025_blanks_refined_noChelex.csv", header=T, row.names = 1)
Meta 

#########################
#18S

setwd("/Users/mghotbi/Desktop/Publications/Second Chapter_96-well plate_1.5Timeseries /2A_method_HPC/2A_method/marjan_merged_euk_table.qza.biom")
#Chloroplast

#OTU<-read.csv("SK.learn.ASV_table.csv",header=T, row.names = 1)
OTU<-read.csv("~/Desktop/Publications/Second Chapter_96-well plate_1.5Timeseries /2A_method_HPC/2A_method/18S_ASV_table_No_metazoa.csv",header=T, row.names = 1)
dim(OTU)
#TAX<-read.csv("SK.learn_tax_table.csv",header=T, row.names = 1)
TAX<-read.csv("~/Desktop/Publications/Second Chapter_96-well plate_1.5Timeseries /2A_method_HPC/2A_method/18S_tax_table_parsed_Nometazoa.csv",header=T, row.names = 1)
#Meta <-read.csv("~/Desktop/Publications/Second Chapter_96-well plate_1.5Timeseries /2A_method_HPC/method.gpt.counterparts-full-rep.qza.biom/Mtadata_fullrep_GPT_counterpart_corrected_allmethods.csv",header=T, row.names = 1)
#Meta <-read.csv("/Users/mghotbi/Desktop/Publications/Second Chapter_96-well plate_1.5Timeseries /2A_method_HPC/2A_method/Method_metadata_with_all_blanks_4_equal_replicates_No_10X_noBeadclean_Final_12Jan2025.csv",header=T, row.names = 1)
Meta <-read.csv("Tab.1_Method_metadata_18March_2025_blanks_refined_noChelex.csv", header=T, row.names = 1)
unique(Meta$Method) 


##################################
#spike cell no detection 

## (Tab.1_Method_metadata_18March_2025_blanks_refined_noChelex.csv)


#Prokaryotes
setwd("/Users/mghotbi/Desktop/Publications/Second Chapter_96-well plate_1.5Timeseries /2A_method_HPC/2A_method")

OTU<-read.csv("ASV_Table_marjan_prok_merge_29March_2025_allsamples.csv",header=T, row.names = 1)
#TAX<-read.csv("gg2_taxonomy_for_silva_selected_codes_without_mit_chl_Unassigned_prok_correct_parsed.csv",header=T, row.names = 1)
TAX<-read.csv("gg2-taxonomy_silvaselectedcodes_No-mit-chlo-Unassigned_prok_correct_parsed_tidy_23March2025.csv",header=T, row.names = 1)
#Meta <-read.csv("Method_metadata_with_all_blanks_4_equal_replicates_No_10X_noBeadclean_Final_12Jan2025.csv", header=T, row.names = 1)
#Meta <-read.csv("Tab.1_Method_metadata_18March_2025_blanks_refined_noChelex.csv", header=T, row.names = 1)
Meta <-read.csv("Method_Metadata_for_extraction_detection_treshold_cell_number.csv", header=T, row.names = 1)
TAX
Meta



identical(sample_names_to_keep,rownames(Meta))

any(sample_names_to_keep %in% colnames(OTU))








Meta$SampleID <- rownames(Meta)
sample_names_to_keep <-Meta$SampleID
sample_names_to_keep


colnames(Meta)
rownames(Meta)
dim(Meta)
#[1] 211  12

filtered_OTU <- OTU[, colnames(OTU) %in% sample_names_to_keep]
dim (filtered_OTU)
OTU<-filtered_OTU

OTU
colnames(Meta)
rownames(Meta)

OTU <- otu_table(as.matrix(OTU), taxa_are_rows = TRUE)
rownames(OTU)
colnames(OTU)

TAX <- tax_table(as.matrix(TAX))
rownames(TAX)
colnames(TAX)

Meta <- sample_data(Meta)
dim(Meta)
sample_names(Meta)
head(Meta)

TAX
OTU
Meta


ps<- merge_phyloseq (OTU,TAX,Meta)
ps

Meta <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
Meta$LibrarySize <- sample_sums(ps)
Meta
Meta[,c(1,2:5,12)]
#write_csv(Meta, "18S_marjan_nobeadclean_noX10_12Jan25_Physeq_with_librarysize_equal.csv")
# now we want to do the decontam 


metadata <-sample_data(ps)
head(metadata)
dim(metadata)

#Inspect Library Sizes
#whether that sample was a true positive sample or a negative control:


library(ggplot2)
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_blank)) + geom_point()

Meta<-df
### there are two blanks which are not blank 
# I want to find them 


blank <-subset_samples(ps, Sample_blank == "Blank" )
blank 

sample_data(blank)

df2 <- as.data.frame(sample_data(blank)) # Put sample_data into a ggplot-friendly data.frame
df2$LibrarySize <- sample_sums(blank)
df2$SampleID <- rownames(df2)
df2
df2 <- df2[order(df2$LibrarySize),]
df2$Index <- seq(nrow(df2))
ggplot(data=df2, aes(x=Index, y=LibrarySize, color=Sample_blank)) + geom_point()+
  geom_text(aes(label = SampleID), hjust = 1.1, vjust = -0.5, size = 3) +  # Add labels
  theme_minimal() +
  labs(x = "Index", y = "Library Size", title = "Library Size by Sample")

write_csv(df, "Method_lib_size_standard.csv")

#These 2 are spiked bacteria added to blank samples , to test the efficiency of extraction of the spiked ones 
#Filterblank.26Oct22.spiked.32     Filterblank  blank        bulk        Blank conventional  Blank Blank       32306    11
#extractionblank.7Nov22.46     Extractionblank  blank        bulk        Blank conventional  Blank Blank       40509    12



library(decontam)

#Identify Contaminants - Prevalence

sample_data(ps)$is.neg <- sample_data(ps)$Sample_blank == "Blank"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.05)
contamdf.prev05

table(contamdf.prev05$contaminant)


contaminants <- rownames(contamdf.prev05[contamdf.prev05$contaminant == TRUE, ])
contaminants  # List of contaminants (e.g., ASV IDs or OTU IDs)


contaminants_prevelance <- contamdf.prev05[contamdf.prev05$contaminant == TRUE, ]
contaminants_prevelance
contamdf.prev05
write.csv(contaminants_prevelance, "contaminants_prevelance_0.05_rarefaction4000.csv")



@#Final conclusion based on inspection


#Spiked bacteria ASVs
# "1ae1bd31e0995f907062cc851d803bd7","39e27cbe03e2b71dbc8a5a1c8b8e7171","ac3a3e27aa0a24e16938ecd9de1c1060","1b468de6094a372d7fbc7d58b16f81c6"

#16S
###the 3 and all 13+9
#the ones from inspection and decontam and despiked 


###final 0.5 prev checked that I will use 

ASVs_to_remove <- c(
  "ac164649e33237f799c1ffee4f344d5d",
  "92fb114641e27921122496428e2ef0dd",
  "3a8a7409bdc33da14a01ec3e3313b9a7",
  "3fcfe36ad281b523c2a5b9793512f5ad",
  "863728e1cec6befd5ba02d15baef4c36",
  "7c081f19c943c8d11819cdce94984719",
  "1c16af05dfb4ec260717d24d0b9b8274",
  "9076f0dc1a57c84a55fd88e445511903",
  "1af47a3186d12f284691296661bc7310",
  "f2d099ef5556f02539ded15181a3d994",
  "fd55e4491e6b3cf2259827d9f78367d6",
  "01cc2516697cf1b59d3bdcdc881c233e",
  "637ddfd0682928e2e0c4e0001883425e",
  "1ae1bd31e0995f907062cc851d803bd7",
  "39e27cbe03e2b71dbc8a5a1c8b8e7171",
  "ac3a3e27aa0a24e16938ecd9de1c1060",
  "1b468de6094a372d7fbc7d58b16f81c6",
  "c728db0a39ffcd40e58dd84831880fde",
  "e1aee885aa820fc1bbc8eea6c27cdc3d",
  "ac70ecd361a0c67c771380f7aac46cc7",
  "79f37fee0660e917bd1debe546718bad",
  "4db5ac9fbbb89c5fc979475084c8c596",
  "645dd8b575c7e7c0952933bbe90b4bb9",
  "dd4ac7657b2665699cffd5a9b939a260",
  "ac9b52e609fb05e7ede04ae1e6a6ff2f",
  "e157cf6c27f73286aa20b0107e08e7d1"
)

##last 9 are standard com which remained undetected

###all 0.05 plus last ones 

ASVs_to_remove <- c(
  "ac164649e33237f799c1ffee4f344d5d",
  "92fb114641e27921122496428e2ef0dd",
  "863728e1cec6befd5ba02d15baef4c36",
  "9076f0dc1a57c84a55fd88e445511903",
  "1c16af05dfb4ec260717d24d0b9b8274",
  "52fbe94fdd249c3c74621ff9ca5a4487",
  "63f512c872fcd1c21b3b7ff7c1fc31fa",
  "81149eedf11c7f9f9ba8707082ee993a",
  "a83209dccceda9136aeabc58ede7df7b",
  "bae4b88b7cff117c98b2395ac7ab1038",
  "01cc2516697cf1b59d3bdcdc881c233e",
  "637ddfd0682928e2e0c4e0001883425e",
  "f2d099ef5556f02539ded15181a3d994",
  "fd55e4491e6b3cf2259827d9f78367d6",
  "0584d6ff6b1991a827c7ed32025d3d57",
  "05c159796be01e86e66b7ee7a9748b2a",
  "0adbc7e3ad3a8d3cac6ef2d3135ac267",
  "1f37e19c844e9bd8878b16056a44d9b4",
  "30aa14a99983077dcc7e593a6371440f",
  "3fcfe36ad281b523c2a5b9793512f5ad",
  "5f42e03097ae97d5659271915b275879",
  "654639f77838c145e2566123862ea66b",
  "7363aac88b5325ac49cb469e284e10dd",
  "7672f20fb33102ef6f580077260ba892",
  "87d12797d5311573897bdb415cff26be",
  "9110c4bc0570a8f5bf9658b76a52b848",
  "9426cd8b391e36a748e9014eb26ada1a",
  "c833a3b511abe63cbb783e2671d5e919",
  "f9e709725cf3c5fdde339d139c386406",
  "3a8a7409bdc33da14a01ec3e3313b9a7",
  "7c081f19c943c8d11819cdce94984719",
  "1af47a3186d12f284691296661bc7310",
  "1ae1bd31e0995f907062cc851d803bd7",
  "39e27cbe03e2b71dbc8a5a1c8b8e7171",
  "ac3a3e27aa0a24e16938ecd9de1c1060",
  "1b468de6094a372d7fbc7d58b16f81c6",
  "c728db0a39ffcd40e58dd84831880fde",
  "e1aee885aa820fc1bbc8eea6c27cdc3d",
  "ac70ecd361a0c67c771380f7aac46cc7",
  "79f37fee0660e917bd1debe546718bad",
  "4db5ac9fbbb89c5fc979475084c8c596",
  "645dd8b575c7e7c0952933bbe90b4bb9",
  "dd4ac7657b2665699cffd5a9b939a260",
  "ac9b52e609fb05e7ede04ae1e6a6ff2f",
  "e157cf6c27f73286aa20b0107e08e7d1"
)




###chloroplast contam 
#0



#FOR 18S 0.05 prevelance 

#bfb7f4c25584efd002c8192a00bac3a7 
#ca50076dd76b0097d37c4967b039586e

ASVs_to_remove <- c("bfb7f4c25584efd002c8192a00bac3a7",
                    "ca50076dd76b0097d37c4967b039586e"
)



#Filter out the contaminant ASVs
ps_filtered <- prune_taxa(!(taxa_names(ps) %in% ASVs_to_remove), ps)
# Check results
#ps_filtered <-ps



########################################################

##### Load your data

#To calculate Spearman's correlation between read counts from microvolume extractions (mv) and bulk extractions ()bulk
#in your column Sample_size, follow these steps:

#Ensure you have two groups of data: mv read counts and bulk read counts.
#Organize your data into a table with two columns:
#Column 1: MV read counts
#Column 2: Bulk read counts

#physeq <- ps
physeq <- ps_filtered
Meta <- as.data.frame(sample_data(physeq)) # Put sample_data into a ggplot-friendly data.frame
Meta$LibrarySize <- sample_sums(physeq)
Meta

#write_csv(Meta, "16S_marjan_Physeq_with_librarysize_after_Decontam_Despike2_12Jan25_equal_right.csv")
dim(Meta)

######################################
library(phyloseq)
physeq

physeq2 <- prune_samples(sample_sums(physeq) > 0, physeq)
physeq2

physeq3<-subset_samples(physeq2, Purification != "bead-cleaned")
physeq3
#137 -2 chelexX10
#134
#69, 18s
physeq4<-subset_samples(physeq3, Type != "Blank")
physeq4

physeq5<-subset_samples(physeq4, KFTno != "Blank")
physeq5

physeq6 <-subset_samples(physeq5, Type != "Spike")
physeq6

physeq6 <-subset_samples(physeq6, Method != "FTPIG-ch-X10")
physeq6

#85 prok, rare 79
#84 chloroplast , rare 76 
#63 18s, rare 59
physeq6@sam_data$Method

sample_depth<- sample_sums(physeq6)
min_depth<- min(sample_depth)
min_depth


sample_data(physeq6)

physeq_rare<-rarefy_even_depth(physeq6,rngseed = 02032020,sample.size = 4000) #prokaryotes
physeq_rare<-rarefy_even_depth(physeq6,rngseed = 02032020,sample.size = 2000) #prokaryotes
physeq_rare<-rarefy_even_depth(physeq6,rngseed = 02032020,sample.size = 100) #chloroplast
physeq_rare<-rarefy_even_depth(physeq6,rngseed = 02032020,sample.size = 100) #chloroplast

physeq_rare

saveRDS(physeq_rare, "physeq_rare_4000rared.rds")

otu.rare = otu_table(physeq6)
otu.rare = as.data.frame(t(otu.rare))
sample_names = rownames(otu.rare)

# we will use vegan rarecurve 
library(vegan)
otu.rarecurve = rarecurve(otu.rare, step = 100,cex=0., label = T)
otu.rarecurve
physeq6.rarefied = rarefy_even_depth(physeq6, rngseed=1, sample.size=0.9*min(sample_sums(physeq6)), replace=F)
physeq6.rarefied

rarecurve(as.data.frame(t(otu_table(physeq6))), step=1000, cex=0)


# Load required libraries
library(phyloseq)
library(ggplot2)

# Perform rarefaction and store results
rare_data <- rarecurve(as.data.frame(t(otu_table(physeq6))), step=1000, cex=0.2, plot=FALSE)


phy.vegan <- as(otu_table(physeq6), "matrix")
phy.vegan
vegan::rarecurve(phy.vegan, 
                 step=100, 
                 main = "Alpha Rarefaction Curve",
                 cex=0.0,
                 label = T)





#physeq_transformed <- transform_sample_counts(physeq_rare, sqrt)
physeq_transformed
physeq_rare<-physeq_transformed

otu_table(physeq_transformed)

  
Meta_rare<-sample_data(physeq_rare)
levels(Meta_rare$Method)
dim(Meta_rare)

  
                      
ord <- phyloseq::ordinate(physeq_rare, method = "PCoA", distance = "bray")

# Inspect the ordination results
ord

plot_ordination(physeq_rare, ord, type = "samples")

Extracted_volume<- as.factor(Meta$Extracted_Volume)


physeq_rare@sam_data$KFTno <- factor(physeq_rare@sam_data$KFTno, 
                                     levels = c(83,85,92,100,101,102,114))  # Sort in ascending order (or use decreasing = TRUE for descending)


pcoa_plot <- plot_ordination(physeq_rare, ord, type = "samples", color = "KFTno", shape = "Method") +
  scale_colour_manual(values = c(
    "100" = "red",  # Orange
    "101" = "#56B4E9",  # Blue
    "102" = "#009E73",  # Green
    "11"  = "#F0E442",  # Yellow
    "114" = "#0072B2",  # Dark Blue
    "31"  = "#D55E00",  # Vermillion
    "83"  = "#CC79A7",  # Reddish Purple
    "84"  = "#CC33CC",  # Magenta (replaced gray)
    "85"  = "pink",  # Magenta (replaced gray)
    "86"  = "#FDE725",  # Bright Yellow-Green
    "90"  = "orange",  # Soft Blue
    "92"  = "blue"  ,
    "95"  = "darkgreen"   
  )) +
  #geom_text(aes(label = Sample_code), vjust = -0.5, hjust = -0.5, size = 3,  max.overlaps = 10) +  # Add Sample_code as labels
  
  geom_point(size = 5, alpha = 0.7) +  # Increased alpha for better visibility
  theme_minimal() +
  labs(title = "PCoA of Bray-Curtis Distances",
       color = "Your Metadata Label") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))+
  #geom_text_repel(aes(label = SampleID), size = 3, max.overlaps = 20) +  # Add labels
  theme_bw()
pcoa_plot

unique(Meta$KFTno)



#sample_data(physeq5)$Extracted_volume <- as.numeric(Meta5$Extracted_volume) 
head(sample_data(physeq_rare))
s<-sample_data(physeq_rare)
s$Extracted_volume




##################


unique(physeq_rare@sam_data$Method) 

# Convert 'Extracted_Volume' to a factor to treat it as a discrete variable
physeq_rare@sam_data$Extracted_volume <- factor(physeq_rare@sam_data$Extracted_volume, 
                                                levels = c(100, 200, 400, 1000, 50000))


# Now plot the PCoA with 'Extracted_Volume' as a discrete size
pcoa_plot <- plot_ordination(physeq_rare, ord, type = "samples", color = "KFTno",  shape = "Method") +
  scale_shape_manual(values = c(
    "conventional" = 19,  # Circle
    "FT" = 15,            #Square
    "FTP" = 17,           # Plus
    "FTPIG" = 9,         # Diamond_Cross
    "FTPIG-X10" = 16,     #  not there
    "FTPIG-ch-X10" = 23   # Triangle
  )) +
  
  scale_colour_manual(values = c(
    "100" = "red",  # Orange
    "101" = "#56B4E9",  # Blue
    "102" = "#009E73",  # Green
    "11"  = "#F0E442",  # Yellow
    "114" = "#0072B2",  # Dark Blue
    "31"  = "#D55E00",  # Vermillion
    "83"  = "#CC33CC",  # Reddish Purple 
    "84"  = "#CC79A7",  # Magenta (replaced gray)
    "85"  = "#FFA07A",  # light salmon
    "86"  = "#FDE725",  # Bright Yellow-Green
    "90"  = "orange",  # Soft Blue
    "92"  = "blue",
    "95"  = "darkgreen"
  )) +
  geom_point(aes(size = Extracted_volume),  alpha = 0.5) +  # Map size to 'Extracted_Volume'
  scale_size_manual(
    name = "Extracted Volume (µL)",  # Legend title
    values = c("100" = 6, "200" = 7, "400" = 8, "1000" = 9, "50000" = 10),  # Assign custom sizes
    breaks = c("100", "200", "400", "1000", "50000"),  # Show specific values in the legend
    labels = c("100 µL", "200 µL", "400 µL", "1000 µL","50000 µL")  # Custom labels
  ) +
  #geom_text(aes(label = Sample_code), vjust = -0.5, hjust = -0.5, size = 3,  max.overlaps = 10) +  # Add Sample_code as labels
  theme_minimal() +
  stat_ellipse(aes(color = KFTno, group = KFTno), type = "t", linetype = "solid", size = 0.8, level = 0.95, alpha = 0.3) +
  theme_minimal() +
  labs(title = "PCoA of Bray-Curtis Distances",
       color = "KFT_sample_numbers") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20)) +
  
  theme_bw()

pcoa_plot





#write_csv(Metaps_after_equal, "Metaps_after_equal.csv")


#################################################################################################

##Barplots 

library(DspikeIn)
physeq_rare@sam_data
physeq_rare2 <- tidy_phyloseq_tse(physeq_rare) 
#physeq_rare2
tax_table(physeq_rare)


#physeq_rare<-readRDS("~/Downloads/Marjan_cleaned_phyloseq.rda" )
physeq_rare

###Order demonstrates the best
#But bars not ordered for abundance 
tax_table(physeq_rare)
color_palette$extended_palette
physeq_rare@sam_data$KFTno <- factor(physeq_rare@sam_data$KFTno, 
                                     levels = c(83,85,92,100,101,102,114))  # Sort in ascending Phylum (or use decreasing = TRUE for descending)

bp_ab <- taxa_barplot(physeq_rare,
                      target_glom = "Order",              
                      treatment_variable = "Code", 
                      abundance_type = "relative",        
                      x_angle = 0,                       
                      fill_variable = "Order",            
                      facet_variable = "KFTno",           
                      top_n_taxa = 30,                   
                      palette = color_palette$extended_palette)

# Modify the ggplot object directly from the returned list
bp_ab$barplot <- bp_ab$barplot + 
  theme(legend.position = "right",                    
        legend.box = "vertical",                    
        legend.spacing.x = unit(0.5, "cm")) +         
  guides(fill = guide_legend(nrow = 30, byrow = TRUE))  

# Print the adjusted plot
print(bp_ab$barplot)


bp_ab$barplot <- bp_ab$barplot + 
  facet_grid(~KFTno, scales = "free_x") + 
  theme(legend.position = "right",                    
        legend.box = "vertical",                    
        legend.spacing.x = unit(0.5, "cm"),         # Adjust legend spacing
        axis.text.x = element_text(angle = 90,      # Rotate x-axis text vertically
                                   vjust = 0.5,    # Adjust vertical alignment
                                   hjust = 1,      # Adjust horizontal alignment
                                   size = 14),     # Set font size for x-axis text
        axis.text.y = element_text(size = 14))      # Increase y-axis text size

# Print the adjusted plot
print(bp_ab$barplot)


# Barplots

```{r}

# Convert KFTno to factor with defined order
physeq_rare_4000@sam_data$KFTno <- factor(
  physeq_rare_4000@sam_data$KFTno,
  levels = c(83, 85, 92, 100, 101, 102, 114)
)

# Create relative abundance barplot grouped at Order level
bp_ab <- DspikeIn::taxa_barplot(
  physeq_rare_4000,
  target_glom        = "Order",            # Aggregate at Order level
  treatment_variable = "Code",             # Sample-level grouping for x-axis
  abundance_type     = "relative",         # Show relative abundance
  x_angle            = 0,                  # Horizontal x-axis labels
  fill_variable      = "Order",            # Color by Order
  facet_variable     = "KFTno",            # Facet by KFTno (sample origin)
  top_n_taxa         = 30,                 # Show top 30 most abundant Orders
  palette            = color_palette$mix_MG
)

# First round of plot customization
bp_ab$barplot <- bp_ab$barplot +
  theme_minimal(base_size = 14) +
  theme(
    legend.position   = "right",
    legend.box        = "vertical",
    legend.spacing.x  = unit(0.5, "cm")
  ) +
  guides(
    fill = guide_legend(nrow = 30, byrow = TRUE)
  ) +
  labs(
    title    = "Relative Abundance of Top 30 Orders per Sample (Faceted by KFTno)",
    y        = "Relative Abundance (%)",
    x        = "Sample Code",
    fill     = "Order"
  )

# Second round: fine-tune x-axis and facets
bp_ab$barplot <- bp_ab$barplot +
  facet_grid(~KFTno, scales = "free_x") +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 14, face = "bold")
  )

# Display final plot
print(bp_ab$barplot)
#############################################more than 7 ranks so 
########################### 

Order_to_filter <- "Embryophyceae"
physeq_rare <- subset_taxa(physeq_rare, Order != Order_to_filter)
tax_ph<-physeq_rare@tax_table
show(tax_ph)



library(DspikeIn)
library(phyloseq)
library(ggplot2)
library(microbiome)
#physeq_rare<-physeq5
tax_ph<-tax_table(physeq_rare)

color_palette$light_MG
physeq_rare@sam_data$KFTno <- factor(physeq_rare@sam_data$KFTno, 
                                     levels = c(83,85,92,100,101,102,114))  # Sort in ascending Phylum (or use decreasing = TRUE for descending)


# Collapse data to the desired rank (e.g., Phylum)
physeq_glom <- tax_glom(physeq_rare, "Genus")

ntaxa(physeq_rare) # 2110
nsamples(physeq_rare) #86
sum(otu_table(physeq_rare)) #86000


# Convert to a data frame for ggplot2
otu_data <- psmelt(physeq_glom)
dim(otu_data)
# Combine low-abundance taxa into "Other"
top_n_taxa <- 30
top_taxa <- names(sort(tapply(otu_data$Abundance, otu_data$Genus, sum), decreasing = TRUE)[1:top_n_taxa])
otu_data$Genus <- ifelse(otu_data$Genus %in% top_taxa, otu_data$Genus, "Other")

# Ensure 'Other' is on top of the stack
otu_data$Genus <- factor(otu_data$Genus, levels = c(setdiff(unique(otu_data$Genus), "Other"), "Other"))

# Create a stacked bar plot using ggplot2
bp_ab <- ggplot(otu_data, aes(x = Code, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~KFTno, scales = "free_x") +  # Group by GPTno while keeping Methods together
  scale_fill_manual(values = color_palette$extended_palette) +
  theme_bw() +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.spacing.x = unit(0.5, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # Rotate x-axis text vertically
  guides(fill = guide_legend(nrow = 55, byrow = TRUE))

bp_ab


###################################################


# Generate the barplot with adjusted legend layout
bp_ab <- taxa_barplot(physeq_rare,
                      target_glom = "Order",              
                      treatment_variable = "Sample_code", 
                      abundance_type = "relative",        
                      x_angle = 90,                       
                      fill_variable = "Order",            
                      facet_variable = "KFTno",           
                      top_n_taxa = 50,                   
                      palette = color_palette$extended_palette)

# Extract the ggplot object from the returned list and customize it
bp_ab$barplot <- bp_ab$barplot + 
  facet_grid(~KFTno, scales = "free_x", drop = TRUE) +
  theme(legend.position = "bottom",                  
        legend.box = "horizontal",                    
        legend.spacing.x = unit(0.5, "cm")) +         
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) 

# Print the modified plot
print(bp_ab$barplot)

physeq@sam_data$Sample_code <- factor(physeq@sam_data$Sample_code)

physeq@sam_data


################################################################

# samples need to be omitted 11,31 ,92
#corrected 114, 95


##########Betadispers 
physeq_rare
table(meta(physeq_rare)$Method)

# conventional           FT          FTP        FTPIG 
#     24                 8           22           18 

#phyto
# conventional           FT          FTP        FTPIG 
#      22                13           22           19 

#microbial Eukaryotes
#   conventional           FT          FTP        FTPIG 
#        15                5           13           16 

bray_dist = phyloseq::distance(physeq_rare, method="bray", weighted=T)

metadata<-sample_data(physeq_rare)
betadispmod <- betadisper(bray_dist, metadata$Method)
#betadispmod <- betadisper(bray_dist, metadata$Extracted_volume)
betadispmod
anova(betadispmod)

# Tukey's HSD

perm<-permutest(betadispmod, pairwise = TRUE, permutations = 9999, p.adjust = "fdr")
perm


anova(betadispmod)


HSD <- TukeyHSD(betadispmod)
HSD_results <- as.data.frame(HSD$group) 
HSD_results$Significance <- cut(HSD_results$`p adj`, 
                                breaks = c(0, 0.001, 0.01, 0.05, 1), 
                                labels = c("***", "**", "*", "ns"))

HSD
plot(HSD)
print(HSD_results)


####Betadispers visualizytion

# these give ugly base R plots

boxplot(betadispmod)

# the blow code block was used to explore the betdispmod object, to identify the correct data to manually curate a dataframe from for ggplot2
# ggvegan is a R package for plotting vegan objects in ggplot2, but it doesn't work for 'class betadisper' objects, as from the betadisper function above. 
centroids=betadispmod$centroids
#write.csv(centroids,file="centroids.csv")
#centroids=read.csv("centroids.csv",header = T,row.names = 1,sep=",")
betadispmod$group.distances
betadispmod$vectors
betadispmod$eig
betadispmod$distances

library(broom)

# Try to make a dataframe with the required information from the betadispermod class object to plot in ggplot2
scores <- vegan::scores(betadispmod) # these are the PCoA loading scores
str(scores)
betscores <- as.data.frame(vegan::scores(betadispmod, 1:4, display="sites")) # by site is by samples
betdist <- as.data.frame(betadispmod$distances) # these are the distances from centroid. The sample IDs are the row names, and the distances are in row 1 (betadispmod$distances). 
betdist
betdistT <- t(betdist) #Transpose it so that the distances are in column one - when you corrected for biases with "bias.adjust = TRUE" the formatting changed, so you no longer need to transpose
betscores$Distance <- betadispmod$distances # now take the distances and add them to the "betscores" dataframe. Make sure that the sampleIDs are in the same order
betcentroid <- as.data.frame(vegan::scores(betadispmod,display="centroid")) # this makes a small dataframe with the centroid values PER HISTORY GROUP
t(betcentroid)
# now add history group to the newaly produced dataframe
betscores$Method <- as.factor(physeq_rare@sam_data$Method)
betscores$Method <- as.factor(physeq_rare@sam_data$Method)


length(rownames(betscores))


centroids


labs(shape = "Extraction Method", color = "Sample (KFTno)")

betscores$Method <- factor(betscores$Method, levels = c("conventional", "FT", "FTP", "FTPIG"))
betscores$KFTno <- physeq_rare@sam_data$KFTno
betscores$Method <- physeq_rare@sam_data$Method
betscores$Extracted_volume <- physeq_rare@sam_data$Extracted_volume

dispersplot<-ggplot(betscores, aes(x = Method, y = Distance)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(aes(shape = Method, color = KFTno,  size = Extracted_volume)) +
  
  scale_shape_manual(values = c(
    "conventional" = 19,
    "FT" = 15,
    "FTP" = 17,
    "FTPIG" = 9
  )) +
  
  scale_colour_manual(values = c(
    "100" = "red", 
    "101" = "#56B4E9",  
    "102" = "#009E73",  
    "11"  = "#F0E442",  
    "114" = "#0072B2",  
    "31"  = "#D55E00",  
    "83"  = "#CC33CC",  
    "84"  = "#CC79A7",  
    "85"  = "#FFA07A",  
    "86"  = "#FDE725",  
    "90"  = "orange",  
    "92"  = "blue",
    "95"  = "darkgreen"
  )) +
  
  labs(shape = "Extraction Method", color = "Sample (KFTno)") +
  theme_bw() +
  ylab("Distance to Centroid") +
  ggtitle("Beta Dispersion") +
  theme(
    text = element_text(family = "Arial", size = 12),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    strip.background = element_rect(colour = "black", fill = "gray98")+
      scale_size_manual(
        name = "Extracted Volume (µL)",  # Legend title
        values = c("100" = 5.5, "200" = 6, "400" = 7, "1000" = 8, "50000" = 8),  # Assign custom sizes
        breaks = c("100", "200", "400", "1000", "50000"),  # Show specific values in the legend
        labels = c("100 µL", "200 µL", "400 µL", "1000 µL","50000 µL")  # Custom labels
      )
  )

dispersplot
ggsave("Betadispersisolationsource_5-6.pdf", height=5, width=6, device="pdf") 

#############

dispersplot <- ggplot(betscores, aes(x = Method, y = Distance)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(aes(shape = Method, color = KFTno, size = 4)) +
  
  scale_shape_manual(values = c(
    "conventional" = 19,
    "FT" = 15,
    "FTP" = 17,
    "FTPIG" = 9
  )) +
  
  scale_colour_manual(values = c(
    "100" = "red", 
    "101" = "#56B4E9",  
    "102" = "#009E73",  
    "11"  = "#F0E442",  
    "114" = "#0072B2",  
    "31"  = "#D55E00",  
    "83"  = "#CC33CC",  
    "84"  = "#CC79A7",  
    "85"  = "#FFA07A",  
    "86"  = "#FDE725",  
    "90"  = "orange",  
    "92"  = "blue",
    "95"  = "darkgreen"
  )) +
  
  labs(shape = "Extraction Method", color = "Sample (KFTno)") +
  theme_bw() +
  ylab("Distance to Centroid") +
  ggtitle("Beta Dispersion") +
  theme(
    text = element_text(family = "Arial", size = 12),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    strip.background = element_rect(colour = "black", fill = "gray98")
  ) 

dispersplot



###### now we run PERMANOVA to see community differences 
library(phyloseq)
library(ggplot2)
library(tidyr)
library(dplyr)
library(microbiome)
library(vegan)
library(pairwiseAdonis)

# I make 3 subsets to have balanced design 


subset_5_12Aug <- subset_samples(physeq_rare, KFTno %in% c(83,85) )
subset_5_12Aug


Extracted_volume <- factor(meta(subset_5_12Aug)$Extracted_volume)
Month <- factor(meta(subset_5_12Aug)$Month)
Temp <- factor(meta(subset_5_12Aug)$Temp)
Salinity <- factor(meta(subset_5_12Aug)$Salinity)
KFTno<- factor(meta(subset_5_12Aug)$KFTno)
Method<- factor(meta(subset_5_12Aug)$Method)
class(Temp) #factor othervise vector

Month
Temp
Salinity
KFTno
Extracted_volume
Method

################beta-dispers for subsets
##I decided betadispers doesnt neeed subsampling 

subset_5_12Aug <- subset_samples(physeq_rare, KFTno %in% c(83,85) )
subset_5_12Aug

subset_5_12Aug
bray_dist = phyloseq::distance(subset_5_12Aug, method="bray", weighted=T)

metadata<-sample_data(subset_5_12Aug)
betadispmod <- betadisper(bray_dist, metadata$Method)
betadispmod <- betadisper(bray_dist, metadata$Extracted_volume)
betadispmod
anova(betadispmod)

# Tukey's HSD
perm<-permutest(betadispmod, pairwise = TRUE, permutations = 9999, p.adjust = "fdr")
perm

#summary(anova_result)

HSD <- TukeyHSD(betadispmod)
HSD_results <- as.data.frame(HSD$group) 
HSD_results$Significance <- cut(HSD_results$`p adj`, 
                                breaks = c(0, 0.001, 0.01, 0.05, 1), 
                                labels = c("***", "**", "*", "ns"))

HSD
plot(HSD)
print(HSD_results)

table(meta(subset_5_12Aug)$Method)

#KFT
adonis_result <- adonis2(phyloseq::distance(subset_5_12Aug, method = "bray") ~  KFTno ,  data = meta(subset_5_12Aug), permutations = 9999)
adonis_result 



#Method
subset_5_12Aug
adonis_result <- adonis2(phyloseq::distance(subset_5_12Aug, method = "bray") ~  Method , strata = KFTno, data = meta(subset_5_12Aug), permutations = 9999)
adonis_result 

pairwise_result <- pairwise.adonis2(phyloseq::distance(subset_5_12Aug,  method = "bray") ~ Method  , 
                                    p.adjust= "Holm", 
                                    strata = meta(subset_5_12Aug)$KFTno,
                                    data = meta(subset_5_12Aug))


pairwise_result
meta(subset_5_12Aug)$KFTno

library(cluster)
library(pairwiseAdonis)
#Volume 
adonis_result <- adonis2(phyloseq::distance(subset_5_12Aug, method = "bray") ~ Extracted_volume , strata = KFTno, data = meta(subset_5_12Aug), permutations = 9999)
adonis_result 

pairwise_result <- pairwise.adonis2(phyloseq::distance(subset_5_12Aug,  method = "bray") ~ Extracted_volume , 
                                    p.adjust= "Holm", 
                                    strata = meta(subset_5_12Aug)$KFTno,
                                    data = meta(subset_5_12Aug))


pairwise_result





##########

subset_2_30Sep_15Nov <- subset_samples(physeq_rare, KFTno %in% c(92,100,114))
subset_2_30Sep_15Nov


Extracted_volume <- factor(meta(subset_2_30Sep_15Nov)$Extracted_volume)
Month <- factor(meta(subset_2_30Sep_15Nov)$Month)
Temp <- factor(meta(subset_2_30Sep_15Nov)$Temp)
Salinity <- factor(meta(subset_2_30Sep_15Nov)$Salinity)
KFTno<- factor(meta(subset_2_30Sep_15Nov)$KFTno)
class(Temp) #factor othervise vector

Month
Temp
Salinity
GPTno
Extracted_volume

###betadispers
subset_2_30Sep_15Nov
bray_dist = phyloseq::distance(subset_2_30Sep_15Nov, method="bray", weighted=T)

metadata<-sample_data(subset_2_30Sep_15Nov)
betadispmod <- betadisper(bray_dist, metadata$Method)
betadispmod <- betadisper(bray_dist, metadata$Extracted_volume)
betadispmod
anova(betadispmod)

# Tukey's HSD

perm<-permutest(betadispmod, pairwise = TRUE, permutations = 9999, p.adjust = "fdr")
perm

#anova_result <- aov(betadispmod)
#summary(anova_result)
anova(betadispmod)

HSD <- TukeyHSD(betadispmod)
HSD_results <- as.data.frame(HSD$group) 
HSD_results$Significance <- cut(HSD_results$`p adj`, 
                                breaks = c(0, 0.001, 0.01, 0.05, 1), 
                                labels = c("***", "**", "*", "ns"))

HSD
plot(HSD)
print(HSD_results)


############
adonis_result <- adonis2(phyloseq::distance(subset_2_30Sep_15Nov, method = "bray") ~  KFTno , data = meta(subset_2_30Sep_15Nov), permutations = 9999)
adonis_result 

pairwise_result <- pairwise.adonis2(phyloseq::distance(subset_2_30Sep_15Nov,  method = "bray") ~ KFTno , 
                                    p.adjust= "Holm", 
                                    data = meta(subset_2_30Sep_15Nov))


pairwise_result


adonis_result <- adonis2(phyloseq::distance(subset_2_30Sep_15Nov, method = "bray") ~  Method , strata = KFTno, data = meta(subset_2_30Sep_15Nov), permutations = 9999)
adonis_result 

pairwise_result <- pairwise.adonis2(phyloseq::distance(subset_2_30Sep_15Nov,  method = "bray") ~ Method  , 
                                    p.adjust= "Holm", 
                                    strata = meta(subset_2_30Sep_15Nov)$KFTno,
                                    data = meta(subset_2_30Sep_15Nov))


pairwise_result
meta(subset_2_30Sep_15Nov)$KFTno

adonis_result <- adonis2(phyloseq::distance(subset_2_30Sep_15Nov, method = "bray") ~ Extracted_volume , strata = KFTno, data = meta(subset_2_30Sep_15Nov), permutations = 9999)
adonis_result 

pairwise_result <- pairwise.adonis2(phyloseq::distance(subset_2_30Sep_15Nov,  method = "bray") ~ Extracted_volume , 
                                    p.adjust= "Holm", 
                                    strata = meta(subset_2_30Sep_15Nov)$KFTno,
                                    data = meta(subset_2_30Sep_15Nov))


pairwise_result


############################

subset_4_7Oct <- subset_samples(physeq_rare, KFTno %in% c(101,102))
subset_4_7Oct

Extracted_volume <- factor(meta(subset_4_7Oct)$Extracted_volume)
Month <- factor(meta(subset_4_7Oct)$Month)
Temp <- factor(meta(subset_4_7Oct)$Temp)
Salinity <- factor(meta(subset_4_7Oct)$Salinity)
KFTno<- factor(meta(subset_4_7Oct)$KFTno)
class(Temp) #factor othervise vector

Month
Temp
Salinity
GPTno
Extracted_volume


###betadispers
subset_4_7Oct
bray_dist = phyloseq::distance(subset_4_7Oct, method="bray", weighted=T)

metadata<-sample_data(subset_4_7Oct)
betadispmod <- betadisper(bray_dist, metadata$Method)
betadispmod <- betadisper(bray_dist, metadata$Extracted_volume)
betadispmod
anova(betadispmod)

# Tukey's HSD

perm<-permutest(betadispmod, pairwise = TRUE, permutations = 9999, p.adjust = "fdr")
perm

#anova_result <- aov(betadispmod)
#summary(anova_result)
anova(betadispmod)

HSD <- TukeyHSD(betadispmod)
HSD_results <- as.data.frame(HSD$group) 
HSD_results$Significance <- cut(HSD_results$`p adj`, 
                                breaks = c(0, 0.001, 0.01, 0.05, 1), 
                                labels = c("***", "**", "*", "ns"))

HSD
plot(HSD)
print(HSD_results)




adonis_result <- adonis2(phyloseq::distance(subset_4_7Oct, method = "bray") ~  KFTno , data = meta(subset_4_7Oct), permutations = 9999)
adonis_result 


adonis_result <- adonis2(phyloseq::distance(subset_4_7Oct, method = "bray") ~  Method , strata = KFTno, data = meta(subset_4_7Oct), permutations = 9999)
adonis_result 

pairwise_result <- pairwise.adonis2(phyloseq::distance(subset_4_7Oct,  method = "bray") ~ Method  , 
                                    p.adjust= "Holm", 
                                    strata = meta(subset_4_7Oct)$KFTno,
                                    data = meta(subset_4_7Oct))


pairwise_result



adonis_result <- adonis2(phyloseq::distance(subset_4_7Oct, method = "bray") ~ Extracted_volume , strata = meta(subset_4_7Oct)$KFTno , data = meta(subset_4_7Oct), permutations = 9999)
adonis_result 

pairwise_result <- pairwise.adonis2(phyloseq::distance(subset_4_7Oct,  method = "bray") ~ Extracted_volume , 
                                    p.adjust= "Holm", 
                                    strata = meta(subset_4_7Oct)$KFTno ,
                                    data = meta(subset_4_7Oct))


pairwise_result


nrow(meta(subset_4_7Oct)) # Check number of rows in metadata
length(KFTno)   
meta(subset_4_7Oct)$KFTno

#######################################################################

##All data together 
physeq_rare
Meta$Extracted_volume
Extracted_volume <- factor(meta(physeq_rare)$Extracted_volume)
Month <- factor(meta(physeq_rare)$Month)
Temp <- factor(meta(physeq_rare)$Temp)
Salinity <- factor(meta(physeq_rare)$Salinity)
KFTno<- factor(meta(physeq_rare)$KFTno)
class(Temp) #factor othervise vector

Month
Temp
Salinity
KFTno
Extracted_volume

table(meta(physeq_rare)$Method, useNA = "always")
table(sample_data(physeq_rare)$Method, useNA = "always")
table(sample_data(physeq_rare)$Extracted_volume, useNA = "always")
table(sample_data(physeq_rare)$Month, useNA = "always")
table(sample_data(physeq_rare)$Temp, useNA = "always")
table(sample_data(physeq_rare)$Salinity, useNA = "always")
#     conventional      FT          FTP         FTPIG         <NA> 
#          27           11           22           19            0                 

# FTP and FTPIG are common between all samples and can be tested for all but the rest need to be excluded and subsamples separately 
# Comparing Methods within Each Time Point:
#  If you want to test differences for each time point individually (e.g., are the extraction methods significantly different at each time point), you could split your dataset by time point and run individual PERMANOVA tests.




#KFTno
adonis_result <- adonis2(phyloseq::distance(physeq_rare, method = "bray") ~  KFTno ,  data = meta(physeq_rare), permutations = 9999)
adonis_result 

pairwise_result <- pairwise.adonis2(phyloseq::distance(physeq_rare,  method = "bray") ~ KFTno , 
                                    p.adjust= "Holm", 
                                    data = meta(physeq_rare))


pairwise_result


adonis_result <- adonis2(phyloseq::distance(physeq_rare, method = "bray") ~  Temp*Salinity*Month ,  data = meta(physeq_rare), permutations = 9999)
adonis_result 



#Method 
adonis_result <- adonis2(phyloseq::distance(physeq_rare, method = "bray") ~  Method , strata = KFTno, data = meta(physeq_rare), permutations = 9999)
adonis_result 



pairwise_result <- pairwise.adonis2(phyloseq::distance(physeq_rare,  method = "bray") ~ Method , 
                                    p.adjust= "Holm", 
                                    strata = KFTno,
                                    data = meta(physeq_rare))


pairwise_result





#Volume 
adonis_result <- adonis2(phyloseq::distance(physeq_rare, method = "bray") ~ Extracted_volume , strata = KFTno, data = meta(physeq_rare), permutations = 9999)
adonis_result 

pairwise_result <- pairwise.adonis2(phyloseq::distance(physeq_rare,  method = "bray") ~ Extracted_volume , 
                                    p.adjust= "Holm", 
                                    strata = KFTno,
                                    data = meta(physeq_rare))


pairwise_result

"Holm"
"bonferroni"

#####################################################################

#Dseq2 between conventional and FTPIG 
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2", force = TRUE)

library(DESeq2)
library(microbiome)



Extracted_volume <- factor(meta(Subset_Dseq)$Extracted_volume)
Month <- factor(meta(Subset_Dseq)$Month)
Temp <- factor(meta(Subset_Dseq)$Temp)
Salinity <- factor(meta(Subset_Dseq)$Salinity)
KFTno<- factor(meta(Subset_Dseq)$KFTno)


#Subset_1 <- subset_samples(physeq_rare, KFTno %in% c ( "101", "102") )

Subset_Dseq <- subset_samples(physeq_rare, Method %in% c ( "conventional", "FTPIG") )
meta(Subset_Dseq)
Subset_Dseq
levels(as.factor(sample_data(Subset_Dseq)$Method))

dds <- phyloseq_to_deseq2(Subset_Dseq, ~ Method)
dds

dds <- DESeq2::DESeq(dds)
dds

# ✅ Step 4: Apply Variance Stabilizing Transformation (VST)
vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
vsd



otu_vst <- otu_table(assay(vsd), taxa_are_rows = TRUE)
ps_vst <- phyloseq(otu_vst, tax_table(Subset_Dseq), sample_data(Subset_Dseq))
ps_vst
otu_vst


#Subset_Dseq_rel = transform_sample_counts(Subset_Dseq, function(x)100* x/sum(x) )
physeq_transformed <- transform_sample_counts(Subset_Dseq, sqrt)




results_DESeq2 <- perform_and_visualize_DA(
  obj = ps_vst  ,
  method = "DESeq2",
  group_var = "Method",
  contrast = c("FTPIG", "conventional" ),
  output_csv_path = "DA_DESeq2.csv",
  target_glom = "Genus",
  significance_level = 0.05,
)

results_DESeq2

### or 

# ✅ Step 5: Extract differential abundance results
#res <- DESeq2::results(dds, contrast = c("Method", "FTPIG", "conventional"), alpha = 0.05)
#res



write.csv(results_DESeq2$results, "results_DESeq2_Plagibacter_UBA_1268.csv")

head(results_DESeq2$results)  #  sig taxa
results_DESeq2$results
results_DESeq2_rel$plot
results_DESeq2$bar_plot
tax_table(results_DESeq2$obj_significant) 
tax_data <- as.data.frame(tax_table(results_DESeq2$obj_significant))
tax_data
le(Subset_Dseq)
################

Phy_deseq<- phyloseq_to_deseq2(Subset_Dseq, ~ Method)
library(edgeR)
sizeFactors(Phy_deseq) <- calcNormFactors(counts(Phy_deseq))
bac.diff = DESeq(Phy_deseq, test="Wald", fitType="parametric")
bac.diff
bac.diff1 <- bac.diff[which(bac.diff$padj < alpha), ]
bac.diff1 <- cbind(as(bac.diff1, "data.frame"), as(tax_table(bac.crap.even)[rownames(bac.diff1), ], "matrix"))
x = tapply(bac.sigtabFUN1$log2FoldChange, bac.sigtabFUN1$Phylum, function(x) max(x))
x = sort(x, TRUE)

######################################################################

#Shannon

# 1.We calculate Shannon Diversity:

diversity <- estimate_richness(physeq_rare, measures = "Shannon")
diversity

diversity_data <- cbind(diversity, meta(physeq_rare))
diversity_data

shapiro.test(diversity_data$Shannon)
#data:  diversity_data$Shannon
#W = 0.95897, p-value = 0.0046
#shannon diversity is not normally distributed so needs scaling
#W closer to one and pvalue greater than 0.05 mean normal distribution


# Chloroplast
#data:  diversity_data$Shannon
#W = 0.97699, p-value = 0.08609


#So we scale the data 
scale_diversity <- diversity_data %>% mutate_at(c("Shannon"), ~(scale(.) %>% as.vector)) # it only scales prok_shannon col
scale_diversity

dat_diversity<-as.data.frame(scale_diversity)
dat_diversity
dat_diversity$Shannon



#########################################################<<<<<<>>>>######## shapes are methods



bp_shannon <- ggplot(diversity_data, aes(x = KFTno, y = Shannon, fill = as.factor(KFTno))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.8)) +  
  geom_jitter(aes(color = as.factor(KFTno), shape = Method), size =  5, 
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) +  
  theme_minimal() +
  facet_grid(~KFTno, scales = "free_x") +  
  theme(
    legend.position = "right",                    
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  labs(x = "Sample (Code)", y = "Shannon Diversity", title = "Shannon Diversity Boxplots with Custom Colors & Shapes") +
  
  # Custom colors for GPTno (Boxplot & Jitter Points)
  scale_fill_manual(values = c(
    "100" = "red", 
    "101" = "#56B4E9",  
    "102" = "#009E73",  
    "11"  = "#F0E442",  
    "114" = "#0072B2",  
    "31"  = "#D55E00",  
    "83"  = "#CC33CC",  
    "84"  = "#CC79A7",  
    "85"  = "#FFA07A",  
    "86"  = "#FDE725",  
    "90"  = "orange",  
    "92"  = "blue",
    "95"  = "darkgreen")) +
  
  scale_colour_manual(values = c(
    "100" = "red", 
    "101" = "#56B4E9",  
    "102" = "#009E73",  
    "11"  = "#F0E442",  
    "114" = "#0072B2",  
    "31"  = "#D55E00",  
    "83"  = "#CC33CC",  
    "84"  = "#CC79A7",  
    "85"  = "#FFA07A",  
    "86"  = "#FDE725",  
    "90"  = "orange",  
    "92"  = "blue",
    "95"  = "darkgreen")) +
  
  # Custom shapes for Methods
  scale_shape_manual(values = c(
    "conventional" = 19,  # Circle
    "FT" = 15,            # Square
    "FTP" = 17,           # Triangle
    "FTPIG" = 9,         # Diamond_Cross
    "FTPIG-X10" = 16,     # Filled Circle
    "FTPIG-ch-X10" = 3   # Plus
  ))

# Print the plot
print(bp_shannon)

#############################This is the final shannon plot based on GPTno color_shape_black round rim

bp_shannon <- ggplot(diversity_data, aes(x = KFTno, y = Shannon, fill = as.factor(KFTno))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +  
  geom_jitter(aes(color = as.factor(KFTno), shape = Method, fill = as.factor(KFTno), size = Extracted_volume), 
              stroke = 1.5, alpha = 0.7,  # Controls outline thickness & transparency
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) +  
  theme_minimal() +
  facet_grid(~KFTno, scales = "free_x") +  
  theme(
    legend.position = "right",                    
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  labs(x = "Sample (Code)", y = "Shannon Diversity", title = "Shannon Diversity Boxplots with Extracted Volume Mapping") +
  
  # Custom colors for GPTno (Boxplot & Jitter Points)
  scale_fill_manual(values = c(
    "100" = "red", 
    "101" = "#56B4E9",  
    "102" = "#009E73",  
    "11"  = "#F0E442",  
    "114" = "#0072B2",  
    "31"  = "#D55E00",  
    "83"  = "#CC33CC",  
    "84"  = "#CC79A7",  
    "85"  = "#FFA07A",  
    "86"  = "#FDE725",  
    "90"  = "orange",  
    "92"  = "blue",
    "95"  = "darkgreen")) +
  
  scale_colour_manual(values = rep("black", length(unique(diversity_data$KFTno)))) +  # Black outlines
  
  # Custom shapes for Methods
  scale_shape_manual(values = c(
    "conventional" = 21,  # Circle
    "FT" = 22,            # Square
    "FTP" = 24,           # Triangle
    "FTPIG" = 23,         # Diamond
    "FTPIG-X10" = 25,     # Inverted Triangle
    "FTPIG-ch-X10" = 21   # Circle
  )) +
  #coord_cartesian(ylim = c(3.5, 5)) +
  # Custom sizes for Extracted Volume
  scale_size_manual(
    name = "Extracted Volume (µL)",  # Legend title
    values = c("100" = 3, "200" = 5, "400" = 7, "1000" = 9, "50000" = 12),  # Assign custom sizes
    breaks = c("100", "200", "400", "1000", "50000"),  # Show specific values in the legend
    labels = c("100 µL", "200 µL", "400 µL", "1000 µL","50000 µL")  # Custom labels
  )

# Print the plot
print(bp_shannon)



#############################################################################
dat_diversity #scaled diversity 


#Shannon

# 1.We calculate Shannon Diversity:

diversity <- estimate_richness(physeq_rare, measures = "Shannon")
diversity

diversity_data <- cbind(diversity, meta(physeq_rare))
diversity_data

shapiro.test(diversity_data$Shannon)
#data:  diversity_data$Shannon
#W = 0.95897, p-value = 0.0046
#shannon diversity is not normally distributed so needs scaling
#W closer to one and pvalue greater than 0.05 mean normal distribution


# Chloroplast
#data:  diversity_data$Shannon
#W = 0.97699, p-value = 0.08609


Extracted_volume <- factor(meta(dat_diversity)$Extracted_volume)
Month <- factor(meta(dat_diversity)$Month)
Temp <- factor(meta(dat_diversity)$Temp)
Salinity <- factor(meta(dat_diversity)$Salinity)
KFTno<- factor(meta(dat_diversity)$KFTno)

#So we scale the data 
scale_diversity <- diversity_data %>% mutate_at(c("Shannon"), ~(scale(.) %>% as.vector)) # it only scales prok_shannon col
scale_diversity

dat_diversity<-as.data.frame(scale_diversity)
dat_diversity
dat_diversity$Shannon





library(lme4)
library(lmerTest)
library(emmeans)

# Run a linear mixed-effects model
#dat_diversity. --> scaled
#diversity_data. -> Nnot scaled

?????!!!!! here

lmm1 <- lmer(dat_diversity$Shannon ~ KFTno + Method +  (1|Method) , data = dat_diversity)
lmm1
Anova(lmm1)
pairwise_comparisons1 <- emmeans(lmm1, pairwise ~ Method,  adjust= "holm")
pairwise_comparisons1


lmm2 <- lmer(dat_diversity$Shannon ~ KFTno  + Method + (1 | Method) + (1 | Extracted_volume) , data = dat_diversity)
lmm2
Anova(lmm2)
pairwise_comparisons2 <- emmeans(lmm2, pairwise ~ Method,  adjust= "holm")
pairwise_comparisons2



lmm3 <- lmer(dat_diversity$Shannon ~ Method  + (1|Temp:Salinity:Month) , data = dat_diversity)
lmm3
Anova(lmm3)
pairwise_comparisons3 <- emmeans(lmm3, pairwise ~ Method,  adjust= "holm")
pairwise_comparisons3

lmm6 <- lmer(dat_diversity$Shannon ~ KFTno + Method  + (1|Month:Temp) , data = dat_diversity)
lmm6
Anova(lmm6)
pairwise_comparisons6 <- emmeans(lmm6, pairwise ~ Method,  adjust= "holm")
pairwise_comparisons6

pairwise_comparisons6_5 <- emmeans(lmm6, pairwise ~ KFTno:Method,  adjust= "holm")
pairwise_comparisons6_5


lmm4 <- lmer(dat_diversity$Shannon ~ Method  + (1|KFTno) , data = dat_diversity)
lmm4
Anova(lmm4)
pairwise_comparisons4 <- emmeans(lmm4, pairwise ~ Method,  adjust= "holm")
pairwise_comparisons4



lmm5 <- lmer(dat_diversity$Shannon ~ Method  + (1|Temp:Salinity:Month) , data = dat_diversity)
lmm5
Anova(lmm5)
pairwise_comparisons5 <- emmeans(lmm5, pairwise ~ Method,  adjust= "holm")
pairwise_comparisons5



lmm7 <- lmer(dat_diversity$Shannon ~ KFTno + Method + (1| Method/Extracted_volume) , data = dat_diversity)
lmm7
Anova(lmm7)
pairwise_comparisons7 <- emmeans(lmm7, pairwise ~ Method,  adjust= "holm")
pairwise_comparisons7


lmm8 <- lmer(dat_diversity$Shannon ~ KFTno + Method + (1| Extracted_volume) , data = dat_diversity)
lmm8
Anova(lmm8)
pairwise_comparisons8 <- emmeans(lmm8, pairwise ~ Method,  adjust= "holm")
pairwise_comparisons8




lmm9 <- lmer(dat_diversity$Shannon ~ KFTno + Extracted_volume + (1| Method) , data = dat_diversity)
lmm9
Anova(lmm9)
pairwise_comparisons9 <- emmeans(lmm8, pairwise ~ Method,  adjust= "holm")
pairwise_comparisons9


AIC(lmm1, lmm2, lmm3, lmm4, lmm5,lmm6,lmm7,lmm8,lmm9) #--> #based on AIC lmm6 is the best
BIC(lmm1, lmm2, lmm3, lmm4, lmm5,lmm6,lmm7,lmm8)







#lmm <- lmer(dat_diversity$Shannon ~ Method + (1|KFTno:Temp), data = dat_diversity)
summary(lmm1)

# Check significance of fixed effects
anova(lmm)
Anova(lmm)

# Perform pairwise comparisons for Method
pairwise_comparisons <- emmeans(lmm1, pairwise ~ KFTno + Method,  adjust= "holm")
summary(pairwise_comparisons)

#or
pairs(emmeans(lmm, specs = "Method"), adjust = "Holm")



AIC(lmeModel, lmeModel0,lmeModel1)
BIC(lmeModel, lmeModel0,lmeModel1)




#pairs(emmeans(lmm, specs = "Method"), adjust = "Holm")
#contrast              estimate    SE   df t.ratio p.value
#conventional  - FT      -1.031 0.319 88.9  -3.237  0.0085
#conventional  - FTP     -1.156 0.243 89.5  -4.757  <.0001
#conventional  - FTPIG   -0.459 0.231 90.9  -1.988  0.1495
#FT - FTP                -0.125 0.335 88.8  -0.373  0.7099
#FT - FTPIG               0.572 0.317 89.5   1.804  0.1495
#FTP - FTPIG              0.697 0.243 88.9   2.873  0.0204

#Degrees-of-freedom method: kenward-roger 
#P value adjustment: holm method for 6 tests 


#lmm <- lmer(subset_all_methods_101_102$Shannon ~ Method + (1|GPTno:Temp), data = subset_all_methods_101_102)
#summary(lmm)


#################################################################################
#########Samples which had FT, FTP, and FTPIG 

# "101","102","92","114" but I need to exclude Chelex from method here 
#We need to exclude these from our subset 
#FTPIG-ch-X10,FTPIG-X10
#Also exclude those 100ul extra 101 to check consistency of the method 

subset_FT_FTP_FTPIG <- subset(diversity_data,
                              GPTno %in% c("101","102","92","114") &
                                !(Method %in% c("FTPIG-ch-X10","FTPIG-X10")))


#&                                                                              
#!(SampleID %in% c( "D1.SW1", "D1.SW2", "D1.FSW1","D1.FSW2")))                          


table((subset_FT_FTP_FTPIG)$Method, useNA="always")

#prok
# conventional           FT          FTP        FTPIG         <NA> 
#.   16                  6           13           14            0 

#euk
#    conventional       FT          FTP        FTPIG         <NA> 
#         11            7           13           19           0 




ggplot(subset_FT_FTP_FTPIG, aes(x = Method, y = Shannon)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, aes(color = GPTno)) +
  labs(title = "Shannon Diversity Index by Method and Timepoint",
       x = "Method",
       y = "Shannon Diversity Index") +
  theme_bw() +
  scale_fill_brewer(palette = "Set3") +  # Boxplot fill colors for Method
  scale_colour_manual(values = c(
    "100" = "red",  # Orange
    "101" = "#56B4E9",  # Blue
    "102" = "#009E73",  # Green
    "11"  = "#F0E442",  # Yellow
    "114" = "#0072B2",  # Dark Blue
    "31"  = "#D55E00",  # Vermillion
    "83"  = "#CC33CC",  # Reddish Purple 
    "84"  = "#CC79A7",  # Magenta (replaced gray)
    "85"  = "#FFA07A",  # light salmon
    "86"  = "#FDE725",  # Bright Yellow-Green
    "90"  = "orange",  # Soft Blue
    "92"  = "blue",
    "95"  = "darkgreen"
  ), name = "KFT_sample_numbers") +  # Add legend title
  theme(axis.text.x = element_text(angle = 0, hjust = 1))




library(lme4)
library(lmerTest)
library(emmeans)

# Run a linear mixed-effects model

lmm <- lmer(subset_FT_FTP_FTPIG$Shannon ~ Method + (1|GPTno:Temp), data = subset_FT_FTP_FTPIG)
summary(lmm)

# Check significance of fixed effects
anova(lmm)

# Perform pairwise comparisons for Method
pairwise_comparisons <- emmeans(lmm, pairwise ~ Method,  adjust= "holm")
summary(pairwise_comparisons)

#$emmeans
#Method       emmean    SE   df lower.CL upper.CL
#conventional   1.96 0.146 24.6     1.66     2.26
#FT             2.33 0.195 19.8     1.92     2.73
#FTP            2.27 0.131 20.6     2.00     2.55
#FTPIG          2.22 0.121 18.2     1.96     2.47

#Degrees-of-freedom method: kenward-roger 
#Confidence level used: 0.95 

#$contrasts
#contrast             estimate    SE   df t.ratio p.value
#conventional - FT     -0.3641 0.242 40.1  -1.507  0.7092
#conventional - FTP    -0.3107 0.195 42.0  -1.595  0.7092
#conventional - FTPIG  -0.2533 0.188 42.0  -1.345  0.7429
#FT - FTP               0.0534 0.239 37.1   0.224  1.0000
#FT - FTPIG             0.1109 0.229 39.9   0.484  1.0000
#FTP - FTPIG            0.0574 0.176 39.1   0.327  1.0000

#Degrees-of-freedom method: kenward-roger 
#P value adjustment: holm method for 6 tests 


#####################################################################

####Now we want to check all samples with FT, FTP

#We need to exclude these from our subset 
#FTPIG-ch-X10,FTPIG-X10
#Also exclude those 100ul extra 101 to check consistency of the method 

subset_FT_FTP <- subset(diversity_data,
                        GPTno %in% c("100","101","102","92","114","85","83") &
                          !(Method %in% c("FTPIG-ch-X10","FTPIG-X10","FTPIG")))


#&
#!(SampleID %in% c( "D1.SW1", "D1.SW2", "D1.FSW1","D1.FSW2")))



table((subset_FT_FTP)$Method, useNA="always")

#conventional       FT          FTP         <NA> 
#     20            13           18            0 



ggplot(subset_FT_FTP, aes(x = Method, y = Shannon)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, aes(color = GPTno)) +
  labs(title = "Shannon Diversity Index by Method and Timepoint",
       x = "Method",
       y = "Shannon Diversity Index") +
  theme_bw() +
  scale_fill_brewer(palette = "Set3") +  # Boxplot fill colors for Method
  scale_colour_manual(values = c(
    "100" = "red",  # Orange
    "101" = "#56B4E9",  # Blue
    "102" = "#009E73",  # Green
    "11"  = "#F0E442",  # Yellow
    "114" = "#0072B2",  # Dark Blue
    "31"  = "#D55E00",  # Vermillion
    "83"  = "#CC33CC",  # Reddish Purple 
    "84"  = "#CC79A7",  # Magenta (replaced gray)
    "85"  = "#FFA07A",  # light salmon
    "86"  = "#FDE725",  # Bright Yellow-Green
    "90"  = "orange",  # Soft Blue
    "92"  = "blue",
    "95"  = "darkgreen"
  ), name = "KFT_sample_numbers") +  # Add legend title
  theme(axis.text.x = element_text(angle = 0, hjust = 1))




library(lme4)
library(lmerTest)
library(emmeans)

# Run a linear mixed-effects model

lmm <- lmer(subset_FT_FTP$Shannon ~ Method + (1|GPTno:Temp), data = subset_FT_FTP)
summary(lmm)

# Check significance of fixed effects
anova(lmm)

# Perform pairwise comparisons for Method
pairwise_comparisons <- emmeans(lmm, pairwise ~ Method,  adjust= "holm")
summary(pairwise_comparisons)

#No significant differences 
######################





######################################

#This is the same code but with mean and error bar 

library(ggplot2)

ggplot(diversity_data, aes(x = Method, y = Shannon)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, aes(color = GPTno)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), alpha=0.1) +
  
  # Mean dots
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", position = position_dodge(0.8)) +
  
  # Error bars for standard error
  stat_summary(
    fun.data = function(x) {
      data.frame(
        ymin = mean(x) - sd(x)/sqrt(length(x)),  # Lower bound: mean - SE
        ymax = mean(x) + sd(x)/sqrt(length(x))   # Upper bound: mean + SE
      )
    },
    geom = "errorbar", width = 0.25, color = "black", position = position_dodge(0.8)) +
  
  labs(title = "Shannon Diversity Index by Method and Timepoint",
       x = "Method",
       y = "Shannon Diversity Index") +
  theme_bw() +
  scale_fill_brewer(palette = "Set3") +  # Boxplot fill colors for Method
  scale_colour_manual(values = c(
    "100" = "red",  # Orange
    "101" = "#56B4E9",  # Blue
    "102" = "#009E73",  # Green
    "11"  = "#F0E442",  # Yellow
    "114" = "#0072B2",  # Dark Blue
    "31"  = "#D55E00",  # Vermillion
    "83"  = "#CC33CC",  # Reddish Purple 
    "84"  = "#CC79A7",  # Magenta (replaced gray)
    "85"  = "#FFA07A",  # light salmon
    "86"  = "#FDE725",  # Bright Yellow-Green
    "90"  = "orange",  # Soft Blue
    "92"  = "blue",
    "95"  = "darkgreen"
  ), name = "KFT_sample_numbers") +  # Add legend title
  theme(axis.text.x = element_text(angle = 0, hjust = 1))




#####in case I want violin plot but I think box is better


library(ggplot2)

ggplot(diversity_data, aes(x = Method, y = Shannon)) +
  # Violin plot
  geom_violin(trim = FALSE, position = position_dodge(0.8), fill = "lightgray", alpha = 0.5) +
  
  # Mean dots
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", position = position_dodge(0.8)) +
  
  # Error bars for standard error
  stat_summary(
    fun.data = function(x) {
      data.frame(
        ymin = mean(x) - sd(x)/sqrt(length(x)),  # Lower bound: mean - SE
        ymax = mean(x) + sd(x)/sqrt(length(x))   # Upper bound: mean + SE
      )
    },
    geom = "errorbar", width = 0.25, color = "black", position = position_dodge(0.8)) +
  
  # Median line
  stat_summary(fun = median, geom = "crossbar", width = 0.5, size = 0.5, color = "black", position = position_dodge(0.8)) +
  
  # Jitter plot to show individual points
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, aes(color = GPTno)) +
  
  labs(title = "Shannon Diversity Index by Method and Timepoint",
       x = "Method",
       y = "Shannon Diversity Index") +
  theme_bw() +
  scale_fill_brewer(palette = "Set3") +  # Violin fill colors for Method
  scale_colour_manual(values = c(
    "100" = "red",  # Orange
    "101" = "#56B4E9",  # Blue
    "102" = "#009E73",  # Green
    "11"  = "#F0E442",  # Yellow
    "114" = "#0072B2",  # Dark Blue
    "31"  = "#D55E00",  # Vermillion
    "83"  = "#CC33CC",  # Reddish Purple 
    "84"  = "#CC79A7",  # Magenta (replaced gray)
    "85"  = "#FFA07A",  # Light salmon
    "86"  = "#FDE725",  # Bright Yellow-Green
    "90"  = "orange",  # Soft Blue
    "92"  = "blue",
    "95"  = "darkgreen"
  ), name = "KFT_sample_numbers") +  # Add legend title
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


##############################




head(Meta5)

physeq_rare_class<- tax_glom(physeq_rare, taxrank = "Genus" )
ps_class_df <-psmelt(physeq_rare_class)

sample_names(physeq_rare_class)<- ps_class_df$sample.code

ps_class_df$GPTno <- sample_data(physeq_rare_class)$GPTno[match(ps_class_df$Sample, sample_names(physeq_rare_class))]
ps_class_df$GPTno

class_abundance <- ps_class_df %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance))

# Get the top 30 classes
top_classes <- class_abundance$Genus[1:30]
top_classes 
# Filter the dataset to include only top 30 or others
ps_class_df$Genus<- ifelse(ps_class_df$Genus %in% top_classes, ps_class_df$Genus, "Others")


colors <- c(brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), brewer.pal(9, "Set2"))


ggplot(ps_class_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~GPTno, scales = "free_x") +  # Wrap by 'GPTno'
  scale_fill_viridis_d() +  # Use color-blind-friendly color scale
  theme_minimal() +
  labs(title = "Sample Composition at Class Level",
       x = "Sample",
       y = "Abundance",
       fill = "Class") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability



##asign colors 


ps_class_df
# Ensure all classes are assigned a color
top_classes <- unique(ps_class_df$Class)
top_classes

# Define a color palette that includes all of the top 30 classes
# Adding a broader set of distinguishable colors for all the classes.


# Create the plot with the updated color scheme
Meta
ps_class_df$SampleCode <- metadata$sample.code


ggplot(ps_class_df, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~GPTno) +
  scale_fill_manual(values = color_palette) +
  labs(title = "Sample Composition at Class Level (Top 30)", 
       fill = "Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) 






table(sample_data(physeq_rare)$Method)
sample_data(physeq_rare)$Method
unique(sample_data(physeq_rare)$Method)
sample_data(physeq_rare)[sample_data(physeq_rare)$Method == "blank", ]

#Meta2<-meta(physeq_rare)
#Meta2$logdepth <- log(Meta2$LibrarySize + 1)
#log<- as.factor(Meta2$logdepth)
#log






################
################
##Works for method but bars are not ordered for abundance so not
bp_ab <- taxa_barplot(physeq_rare,
                      target_glom = "Order",              # Taxonomic level to glom
                      treatment_variable = "Method", # Variable to use on x-axis
                      abundance_type = "relative",         # Plot absolute abundance
                      x_angle = 0,                       # Rotate x-axis labels by 90 degrees
                      fill_variable = "Order",            # Fill bars by Order
                      facet_variable = "GPTno",            # Facet the plot by Diet
                      top_n_taxa = 100,                    # Show the top 20 taxa
                      palette = color_palette$extended_palette)                    # Use custom color palette (MG)

print(bp_ab$barplot)   



bp_ab <- taxa_barplot(physeq,
                      target_glom = "Order",              
                      treatment_variable = "Method", 
                      abundance_type = "relative",        
                      x_angle = 0,                       
                      fill_variable = "Order",            
                      facet_variable = "GPTno",           
                      top_n_taxa = 120,                   
                      palette = color_palette$extended_palette)

# Modify the ggplot object directly from the returned list
bp_ab$barplot <- bp_ab$barplot + 
  theme(legend.position = "bottom",                    
        legend.box = "horizontal",                    
        legend.spacing.x = unit(0.5, "cm")) +         
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))  

# Print the adjusted plot
print(bp_ab$barplot)



