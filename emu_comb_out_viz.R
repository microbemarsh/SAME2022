###############################################################################
####################___CMI_LAB_EMU_AUTOMATED_PLOTTING____######################
###############################################################################

### Created by Austin G. Marshall and Daniel T. Fuller ###

###############################################################################

### Load libraries for the packages used in this workflow ###
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(janitor)
library(phyloseq)

###############################################################################

### Loading in combined abundance data from Emu ###
a = read.delim(file = "emu-combined-abundance-genus.tsv")

###############################################################################
############___Change_the_order_around_to_plot_in_certain_order___#############
###############___Yes_this_is_ugly_but_who_cares_if_it_works___################
###############################################################################

b = a %>%
  rename("SLR69" = "X22_08_26_SLR69") %>%
  rename("SLR68" = "X22_08_26_SLR68") %>%
  rename("SLR6" = "X22_08_26_SLR6") %>%
  rename("SLR3" = "X22_08_26_SLR3") %>%
  rename("SLR16" = "X22_08_26_SLR16") %>%
  rename("SLR14" = "X22_08_26_SLR14") %>%
  rename("SLR72" = "X22_08_26_SLR72") %>%
  rename("SLR4" = "X22_08_26_SLR4") %>%
  rename("SLR2" = "X22_09_01_SLR2") %>%
  rename("SLR29" = "X22_09_01_SLR29") %>%
  rename("SLR27" = "X22_09_01_SLR27") %>%
  rename("SLR26" = "X22_09_01_SLR26") %>%
  rename("SLR25" = "X22_09_01_SLR25") %>%
  rename("SLR22" = "X22_09_01_SLR22") %>%
  rename("SLR19" = "X22_09_01_SLR19") %>%
  rename("SLR18" = "X22_09_01_SLR18") %>%
  rename("SLR17" = "X22_09_01_SLR17") %>%
  rename("SLR73" = "X22_09_01_SLR73") %>%
  rename("SLR56" = "X22_09_12_SLR56") %>%
  rename("SLR55" = "X22_09_12_SLR55") %>%
  rename("SLR53" = "X22_09_12_SLR53") %>%
  rename("SLR52" = "X22_09_12_SLR52") %>%
  rename("SLR50" = "X22_09_12_SLR50") %>%
  rename("SLR48" = "X22_09_12_SLR48") %>%
  rename("SLR47" = "X22_09_12_SLR47") %>%
  rename("SLR45" = "X22_09_12_SLR45") %>%
  rename("SLR44" = "X22_09_12_SLR44") %>%
  rename("SLR43" = "X22_09_12_SLR43") %>%
  rename("SLR41" = "X22_09_12_SLR41") %>%
  rename("SLR35" = "X22_09_12_SLR35") %>%
  rename("SLR67" = "X22_09_22_SLR67") %>%
  rename("SLR66" = "X22_09_22_SLR66") %>%
  rename("SLR46" = "X22_09_22_SLR46") %>%
  rename("SLR36" = "X22_09_22_SLR36") %>%
  rename("SLR23" = "X22_09_22_SLR23") %>%
  rename("SLR21" = "X22_09_22_SLR21") %>%
  rename("SLR71" = "X22_09_22_SLR71") %>%
  rename("SLR42" = "X22_09_22_SLR42") %>%
  rename("SLR34" = "X22_09_22_SLR34") %>%
  rename("SLR38" = "X22_10_07_SLR38") %>%
  rename("SLR31" = "X22_10_07_SLR31") %>%
  rename("SLR8" = "X22_10_07_SLR8") %>%
  rename("SLR7" = "X22_10_07_SLR7") %>%
  rename("SLR70" = "X22_10_07_SLR70") %>%
  rename("SLR62" = "X22_10_07_SLR62") %>%
  rename("SLR60" = "X22_10_07_SLR60") %>%
  rename("SLR51" = "X22_10_07_SLR51") %>%
  rename("SLR33" = "X22_10_07_SLR33") %>%
  rename("SLR32" = "X22_10_07_SLR32") %>%
  rename("SLR28" = "X22_10_07_SLR28") %>%
  rename("SLR1" = "X22_10_07_SLR1") %>%
  rename("SLR58" = "X22_10_12_SLR58") %>%
  rename("SLR57" = "X22_10_12_SLR57") %>%
  rename("SLR5" = "X22_10_12_SLR5") %>%
  rename("SLR40" = "X22_10_12_SLR40") %>%
  rename("SLR49" = "X22_10_12_SLR49") %>%
  rename("SLR61" = "X22_11_02_SLR61") %>%
  rename("SLR39" = "X22_11_02_SLR39") %>%
  rename("SLR37" = "X22_11_02_SLR37") %>%
  rename("SLR20" = "X22_11_02_SLR20") %>%
  rename("SLR30" = "X22_11_02_SLR30") %>%
  rename("SLR13" = "X22_11_02_SLR13") %>%
  rename("SLR10" = "X22_11_02_SLR10") %>%
  rename("SLR24" = "X22_11_02_SLR24") %>%
  rename("SLR63" = "X22_11_02_SLR63") %>%
  rename("SLR59" = "X22_11_02_SLR59") %>%
  rename("SLR12" = "X22_11_16_SLR12") %>%
  rename("SLR11" = "X22_11_16_SLR11") %>%
  rename("SLR9" = "X22_11_16_SLR9") %>%
  mutate_all(~replace(., is.na(.), 0))

### Grabbing the data without the blank taxa ###
without_blanks = b %>%
  na_if("") %>%
  na.omit()

### Grabbing the blank taxa only ###
with_blanks = b %>%
  filter(genus=="") %>%
  na_if("")

### Summing the blanks together ###
blanksums = with_blanks %>%
  select(-genus) %>%
  colSums() %>%
  t() %>%
  as.data.frame()

### Setting up the naming stuff for the filtering step ###
rownames(without_blanks)<-without_blanks$genus
before_filt <- without_blanks %>%
  select(-genus)

###############################################################################
################ Applying filter and binding rows together ####################
###############################################################################

abund_filt <- apply(before_filt, 2, FUN = function(x) (x/sum(x))<=0.035)
data.cleaned <- before_filt[-which(rowSums(abund_filt)==ncol(before_filt)),]
othercat <- colSums(before_filt[which(rowSums(abund_filt)==ncol(before_filt)),])
bind1 <- rbind(data.cleaned, othercat)
bind2 <- rbind(bind1, blanksums)
final_w_name = bind2 %>%
  rownames_to_column(var = "Sample")

###############################################################################
##########__Stop_here_and_change_what_samples_need_to_be_grabbed__#############
###############################################################################

grab_sums = bind2[c('23','1'),]
grab_blanksums = filter(final_w_name, final_w_name$Sample==1)
grab_taxasums = filter(final_w_name, final_w_name$Sample==23)

### Need to change these numbers to correspond to the number of samples ###
final_other = grab_taxasums[2:71] - grab_blanksums[2:71]

### Need to change the numbers to account for changing datasets ### ### Stuck here on big one ###
bind2_minusothers <- bind2[c(1:22),]

done <- rbind(bind2_minusothers, final_other) 
rownames(done)[length(done[,1])]="Other"

done3 = t(done) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")

###############################################################################
###### Gotta change the formatting once more for ggplot2 to like it ###########
###############################################################################

tbl_reshape <- melt(done3, id.vars = "Sample" , variable.name = "Genus")

### Now We plot, change the title to whatever your experiment is ###

genusplot <- ggplot(tbl_reshape, aes(x = fct_inorder(Sample), y = value, fill = Genus)) +
  geom_bar(colour = "black", stat = "identity", position = "fill") +
  labs(x = "Sample", y = "Relative abundance", title = "SLR_16S_Code_trial") +
  theme_bw() 

###############################################################################
#################### Change color patterns using ggsci ########################
###############################################################################

genusplot_goodone = genusplot + scale_fill_ucscgb(palette = "default", alpha = 1) 

### Final plotting of the figure ###

genusplot_goodone
