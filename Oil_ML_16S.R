install.packages("janitor")
library(tidyverse)
library(janitor)
library(phyloseq)
setwd("~/Documents/Experimental/Oil_ML")

################ IMPORT MTU DATA####################

MTU_taxa<-read_csv("mtuTaxtab2.csv")
MTU_taxa<-column_to_rownames(MTU_taxa,var="...1")
MTU_taxa2 <- as.matrix(MTU_taxa)      #convert to matrix format

MTU_table<-read_csv("mtuFrame2.csv")
SampleName<-as.data.frame(MTU_table$SampleName)
colnames(SampleName)<-"SampleNames"
MTU_table<-MTU_table%>%
  select(-c("rowname", "Run","SampleName","...1"))
counts2 <- as.matrix(MTU_table)            #convert into matrix format
counts2<-as.data.frame(counts2)
counts2[1:5,1:5]
MTU_table2<-cbind(SampleName,counts2)
MTU_table2<-column_to_rownames(MTU_table2,var="SampleNames")
MTU_table2<-t(MTU_table2)

MTU_meta<-SampleName%>%
  separate(SampleNames, c("Site","Time","Condition","Replicate"), sep="[.]", remove=FALSE)%>%
  mutate(Oiled=ifelse(Condition == "CONTROL", "Non-oiled","Oiled"))%>%
  mutate(database = "MTU")%>%
  column_to_rownames(var="SampleNames")

dim(MTU_meta)
dim(MTU_table2)
dim(MTU_taxa)

samdata = sample_data(MTU_meta)
seqtab = otu_table(MTU_table2, taxa_are_rows = TRUE)
taxtab = tax_table(MTU_taxa2)
MTU_ps = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))

taxa_names(seqtab)
sample_names(samdata)

################ IMPORT SEASONS DATA####################

seasons_spring_ASV<-read_csv('spring_seasonsFrame.csv')
seasons_spring_taxa<-read_csv('spring_seasonsTaxtab.csv')
seasons_other_ASV<-read_csv('seasonsFrame.csv')
seasons_other_taxa<-read_csv('seasonsTaxtab.csv')

seasons_spring_ASV<-seasons_spring_ASV%>%
  column_to_rownames(var="...1")
seasons_other_ASV<-seasons_other_ASV%>%
  column_to_rownames(var="...1")

seasons_spring_ASVt<-as.data.frame(t(as.matrix(seasons_spring_ASV)))
seasons_spring_ASVt<-seasons_spring_ASVt%>%
  rownames_to_column(var="ASV")

seasons_other_ASVt<-as.data.frame(t(as.matrix(seasons_other_ASV)))
seasons_other_ASVt<-seasons_other_ASVt%>%
  rownames_to_column(var="ASV")

seasons_all_ASV<-full_join(seasons_other_ASVt,seasons_spring_ASVt)

seasons_all_ASV<-seasons_all_ASV%>%
  replace(is.na(.),0)%>%
  column_to_rownames(var="ASV")
seasons_all_ASV<-as.matrix(seasons_all_ASV)
colnames(seasons_all_ASV)


seasons_all_taxa<-full_join(seasons_other_taxa,seasons_spring_taxa,by="...1")%>%
  select("...1","Kingdom.x","Phylum.x","Class.x","Order.x","Family.x","Genus.x")
colnames(seasons_all_taxa)<-c("ASV","Kingdom","Phylum","Class","Order","Family","Genus")
colnames(seasons_all_taxa)

seasons_all_taxa2<-as.matrix(seasons_all_taxa%>%
  column_to_rownames(var="ASV"))

seasons_meta<-read_csv("Seasons_metadata_new.csv")
seasons_meta<-seasons_meta%>%
  mutate(database="Seasons")%>%
  column_to_rownames(var="Sample_Name")

samdata = sample_data(seasons_meta)
seqtab = otu_table(seasons_all_ASV, taxa_are_rows = TRUE)
taxtab = tax_table(seasons_all_taxa2)
seasons_ps = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))
seasons_ps<-subset_samples(seasons_ps,Oil_Type!="Blank")




################ IMPORT straits DATA####################
straits_taxa<-read_csv("straitsTaxtab.csv")
straits_taxa<-column_to_rownames(straits_taxa,var="...1")
straits_taxa2 <- as.matrix(straits_taxa)      #convert to matrix format

straits_table<-read_csv("straitsFrame.csv")
SampleName<-as.data.frame(straits_table[,1])
colnames(SampleName)<-"SampleNames"
straits_table<-straits_table%>%
  select(-c("rowname","Run","...1"))
counts2_st <- as.matrix(straits_table)            #convert into matrix format
counts2_st<-as.data.frame(counts2_st)
counts2_st[1:5,1:5]
straits_table2<-cbind(SampleName,counts2_st)
straits_table2<-column_to_rownames(straits_table2,var="SampleNames")
straits_table2<-t(straits_table2)

straits_meta<-SampleName%>%
  separate(SampleNames, c("Site","Location","Condition","Replicate"), sep="[.]", remove=FALSE)%>%
  unite("Site_Name",Site:Location,remove = FALSE)%>%
  mutate(Oiled=ifelse(Replicate == "CONTROL", "Non-oiled",ifelse(Replicate=="BAK","Oiled", ifelse(Replicate=="DBIT","Oiled","Environmental"))))%>%
  mutate(database="Straits")%>%
  column_to_rownames(var="SampleNames")

dim(straits_meta)
dim(straits_table2)
dim(straits_taxa)

samdata = sample_data(straits_meta)
seqtab = otu_table(straits_table2, taxa_are_rows = TRUE)
taxtab = tax_table(straits_taxa2)
straits_ps = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))

taxa_names(seqtab)
sample_names(seqtab)


merged<-merge_phyloseq(MTU_ps, straits_ps,seasons_ps)

samplesover1000_all <- subset_samples(merged, sample_sums(merged) > 1000)

any(taxa_sums(samplesover1000_all) == 0)

sum(taxa_sums(samplesover1000_all) == 0)

prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)

any(taxa_sums(prune_samplesover1000_all) == 0)

#for reproducible data
set.seed(81)

rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))

pps<- rarefy_samplesover1000_all

MTU_ps<-subset_samples(pps,database=="MTU")
straits_ps<-subset_samples(pps,database=="Straits")
seasons_ps<-subset_samples(pps,database=="Seasons")

MTU_pcoa <- ordinate(
  physeq = MTU_ps, 
  method = "PCoA", 
  distance = "bray"
)

straits_pcoa <- ordinate(
  physeq = straits_ps, 
  method = "PCoA", 
  distance = "bray"
)

seasons_pcoa <- ordinate(
  physeq = seasons_ps, 
  method = "PCoA", 
  distance = "bray"
)


colors <- c("purple", "orange", "blue", "red", "white")

straits_ps@sam_data
plot_ordination(
  physeq = straits_ps,                                                         #phyloseq object
  ordination = straits_pcoa, type="samples", color="Oiled",shape="Site_Name")                                                #ordination
  geom_point(aes(fill = Oiled, shape = database), size = 5) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  scale_fill_manual(values = colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command


### Straits Plots ###
genusabundance_straits <- straits_ps %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus) 
head(genusabundance_straits)
  
all_straits <- genusabundance_straits %>%
  select(Phylum, Class, Family, Genus, Sample, Abundance, Sample, Oiled, Site_Name,Replicate) %>%
  filter(Abundance != 0) %>%
  mutate(
    Phylum = as.character(Phylum),
    Class = as.character(Class),
    Family = as.character(Family),
    Genus = as.character(Genus)
  )
head(all_straits)  

phylum_straits <- all_straits %>%
  select(Oiled, Site_Name, Phylum, Abundance) %>%  #choose variables to work with
  group_by(Oiled, Site_Name) %>%                   #group by variables used to plot NOT taxonomic level
  mutate(totalSum = sum(Abundance)) %>%                                 #calculate total abundance of each Phylum
  ungroup() %>%                                                         #remove grouping variables
  group_by(Oiled, Site_Name,, Phylum) %>%           #group by same variables PLUS taxonomic level
  summarise(                                                            
    Abundance = sum(Abundance),                                         #sum abundance in each phylum for each  group
    totalSum,
    RelAb = Abundance/totalSum) %>%                                     #calculate relative abundance
  unique()

Color_22_Palette_straits <- c("black", "aquamarine", "darkgray", "gray88", "green","deepskyblue","darkred","orangered","magenta","orange","mediumvioletred","purple","white","darkolivegreen")
Color_22_Palette_seasons <- c("black", "aquamarine", "chocolate1", "darkgray", "gray88", "green","deepskyblue","darkred","orangered","magenta","yellow","orange","cyan","seagreen3","mediumvioletred","purple","white","olivedrab2","lightpink","darkolivegreen")
ggplot(phylum_straits)+
  geom_col(mapping = aes(x = Oiled, y = RelAb, fill = Phylum), color = "black", position = "stack", show.legend = TRUE)+
  facet_grid(cols = vars(Site_Name))+
  ylab("Proportion of Community") +
  xlab(NULL)+
  scale_fill_manual(values = Color_22_Palette_straits) +
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=3,byrow=TRUE))

class_straits <- all_straits %>%
  select(Oiled, Site_Name, Class, Abundance) %>%  #choose variables to work with
  group_by(Oiled, Site_Name) %>%                   #group by variables used to plot NOT taxonomic level
  mutate(totalSum = sum(Abundance)) %>%                                 #calculate total abundance of each Phylum
  ungroup() %>%                             
  group_by(Oiled, Site_Name, Class, totalSum) %>%
  summarise(
    Abundance = sum(Abundance),
    Class = ifelse(Abundance < 0.01, "< 1 %", Class)) %>%               #change Genus label to group low abundance taxa together
  group_by(Oiled, Site_Name, Class, totalSum) %>%  #now group and summarize again to group newly labeled low abundance taxa together
  summarise(
    Abundance = sum(Abundance),
    RelAb = Abundance/totalSum) %>%
  unique()

ggplot(class_straits)+
  geom_col(mapping = aes(x = Oiled, y = RelAb, fill = Class), color = "black", position = "stack", show.legend = TRUE)+
  facet_grid(cols = vars(Site_Name))+
  ylab("Proportion of Community") +
  xlab(NULL)+
  scale_fill_manual(values = Color_22_Palette_straits) +
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, angle =45, vjust = 1, hjust = 1, color = "black"),
        legend.text = element_text(size = 10),
        legend.position = "right",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 12, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 10))+
  guides(fill=guide_legend(ncol=2,byrow=TRUE))

### Sesaons Plots ###
test<-as.data.frame(seasons_ps@otu_table)
genusabundance_season <- seasons_ps %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus) 
head(genusabundance_season)

all_season <- genusabundance_season %>%
  select(Phylum, Class, Family, Genus, Sample, Abundance, Oiled, Season) %>%
  filter(Abundance != 0) %>%
  mutate(
    Phylum = as.character(Phylum),
    Class = as.character(Class),
    Family = as.character(Family),
    Genus = as.character(Genus)
  )
head(all_season)  
class_season <- all_season %>%
  dplyr::group_by(Oiled, Season, Class)%>%
  mutate(Class.1p = ifelse(Abundance < 0.01, "<1%", Class))
head(class_season)

phylum <- all_season %>%
  select(Oiled, Season, Phylum, Abundance) %>%  #choose variables to work with
  group_by(Oiled, Season) %>%                   #group by variables used to plot NOT taxonomic level
  mutate(totalSum = sum(Abundance)) %>%                                 #calculate total abundance of each Phylum
  ungroup() %>%                                                         #remove grouping variables
  group_by(Oiled, Season,, Phylum) %>%           #group by same variables PLUS taxonomic level
  summarise(                                                            
    Abundance = sum(Abundance),                                         #sum abundance in each phylum for each  group
    totalSum,
    RelAb = Abundance/totalSum) %>%                                     #calculate relative abundance
  unique()
"pink","tan"
Color_22_Palette_seasons <- c("black", "aquamarine", "chocolate1", "darkgray", "gray88", "green","deepskyblue","darkred","orangered","magenta","yellow","orange","cyan","seagreen3","mediumvioletred","purple","white","olivedrab2","lightpink","darkolivegreen")


ggplot(phylum)+
  geom_col(mapping = aes(x = Oiled, y = RelAb, fill = Phylum), color = "black", position = "stack", show.legend = TRUE)+
  facet_grid(cols = vars(Season))+
  ylab("Proportion of Community") +
  xlab(NULL)+
  scale_fill_manual(values = Color_22_Palette_seasons) +
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=3,byrow=TRUE))



class_season <- all_season%>%
  select(Oiled, Season, Class, Abundance) %>%  #choose variables to work with
  group_by(Oiled, Season) %>%                   #group by variables used to plot NOT taxonomic level
  mutate(totalSum = sum(Abundance)) %>%                                 #calculate total abundance of each Phylum
  ungroup() %>%                             
  group_by(Oiled, Season, Class, totalSum) %>%
  summarise(
    Abundance = sum(Abundance),
    Class = ifelse(Abundance < 0.01, "< 1 %", Class)) %>%               #change Genus label to group low abundance taxa together
  group_by(Oiled, Season, Class, totalSum) %>%  #now group and summarize again to group newly labeled low abundance taxa together
  summarise(
    Abundance = sum(Abundance),
    RelAb = Abundance/totalSum) %>%
  unique()

ggplot(class_season)+
  geom_col(mapping = aes(x = Oiled, y = RelAb, fill = Class), color = "black", position = "stack", show.legend = TRUE)+
  facet_grid(cols = vars(Season))+
  ylab("Proportion of Community") +
  xlab(NULL)+
  scale_fill_manual(values = Color_22_Palette_seasons) +
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, angle =45, vjust = 1, hjust = 1, color = "black"),
        legend.text = element_text(size = 10),
        legend.position = "right",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 12, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 10))+
  guides(fill=guide_legend(ncol=2,byrow=TRUE))









### MTU Plots ###
MTU_ps@sam_data
genusabundance_MTU <- MTU_ps %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus) 
head(genusabundance_MTU)

all_MTU <- genusabundance_MTU %>%
  select(Phylum, Class, Family, Genus, Sample, Abundance, Oiled, Site ,Condition) %>%
  filter(Abundance != 0) %>%
  mutate(
    Phylum = as.character(Phylum),
    Class = as.character(Class),
    Family = as.character(Family),
    Genus = as.character(Genus)
  )
head(all_MTU)  
phylum_MTU <- all_MTU %>%
  select(Oiled, Site, Phylum, Abundance) %>%  #choose variables to work with
  group_by(Oiled, Site) %>%                   #group by variables used to plot NOT taxonomic level
  mutate(totalSum = sum(Abundance)) %>%                                 #calculate total abundance of each Phylum
  ungroup() %>%                                                         #remove grouping variables
  group_by(Oiled, Site,, Phylum) %>%           #group by same variables PLUS taxonomic level
  summarise(                                                            
    Abundance = sum(Abundance),                                         #sum abundance in each phylum for each  group
    totalSum,
    RelAb = Abundance/totalSum) %>%                                     #calculate relative abundance
  unique()

Color_22_Palette_MTU <- c("black", "aquamarine", "chocolate1", "darkgray", "gray88", "green","deepskyblue","darkred","orangered","magenta","orange","cyan","blue","seagreen3","royalblue","mediumvioletred","purple","white","sienna4","darkolivegreen","darkseagreen","mediumvioletred")

ggplot(phylum_MTU)+
  geom_col(mapping = aes(x = Oiled, y = RelAb, fill = Phylum), color = "black", position = "stack", show.legend = TRUE)+
  facet_grid(cols = vars(Site))+
  ylab("Proportion of Community") +
  xlab(NULL)+
  scale_fill_manual(values = Color_22_Palette_MTU) +
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=3,byrow=TRUE))


class_MTU <- all_MTU%>%
  select(Oiled, Site, Class, Abundance) %>%  #choose variables to work with
  group_by(Oiled, Site) %>%                   #group by variables used to plot NOT taxonomic level
  mutate(totalSum = sum(Abundance)) %>%                                 #calculate total abundance of each Phylum
  ungroup() %>%                             
  group_by(Oiled, Site, Class, totalSum) %>%
  summarise(
    Abundance = sum(Abundance),
    Class = ifelse(Abundance < 0.01, "< 1 %", Class)) %>%               #change Genus label to group low abundance taxa together
  group_by(Oiled, Site, Class, totalSum) %>%  #now group and summarize again to group newly labeled low abundance taxa together
  summarise(
    Abundance = sum(Abundance),
    RelAb = Abundance/totalSum) %>%
  unique()

ggplot(class_MTU)+
  geom_col(mapping = aes(x = Oiled, y = RelAb, fill = Class), color = "black", position = "stack", show.legend = TRUE)+
  facet_grid(cols = vars(Site))+
  ylab("Proportion of Community") +
  xlab(NULL)+
  scale_fill_manual(values = Color_22_Palette_MTU) +
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, angle =45, vjust = 1, hjust = 1, color = "black"),
        legend.text = element_text(size = 10),
        legend.position = "right",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 12, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 10))+
  guides(fill=guide_legend(ncol=2,byrow=TRUE))
Color_22_Palette <- c("black", "blue", "aquamarine", "chocolate1", "darkgray", "gray88", "green","deepskyblue","darkred","orangered","magenta","orange","cyan","black","seagreen3","royalblue","red","purple","white","sienna4","darkolivegreen","darkseagreen","mediumvioletred")

