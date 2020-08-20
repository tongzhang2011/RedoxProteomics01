
#### link ms-gf IDs with MASIC reporter ion intensities----------
MSGF_link_MASIC <- function(msgf, meta, masic) {
  # Dataset_ID will be used to link to MASIC data.
  meta1 <- select(meta, Job, Dataset, Dataset_ID)
  
  d <- merge(msgf, meta1, by = "Job", all.x = TRUE)
  
  d <- select(d,Job,Dataset,Dataset_ID, Scan,Peptide,Protein,MSGF_SpecProb,MSGFDB_SpecEValue,EValue,
              MH,Charge,PrecursorMZ,DelM_PPM,NTT,MSGFScore,PepQValue)
  
  # Note that 130N, Ion_130.135,   is left blank
  ri1 <- select(ri,Dataset, ScanNumber,Ion_126.128:Ion_129.138, Ion_130.141, Ion_131.138)
  
  colnames(d)[3] <- "m_Dataset_ID"
  colnames(d)[4] <- "m_scan_num"
  
  colnames(ri1)[1] <- "m_Dataset_ID"
  colnames(ri1)[2] <- "m_scan_num"
  
  dm <- merge(d, ri1, by = c("m_Dataset_ID", "m_scan_num"), all.x = TRUE)
  
  return(dm)
}
#### link ms-gf IDs with MASIC reporter ion intensities----------

### psm_filter----- START---------
psm_filter <- function(df,PepQValue_cutoff) {
  
  ppm_center <- median(df$DelM_PPM)
  ppm_left <- ppm_center - 10
  ppm_right <- ppm_center + 10
  
  df <- filter(df,DelM_PPM < ppm_right & DelM_PPM > ppm_left)
  df <- filter(df, PepQValue < PepQValue_cutoff)
  
  # df <- filter(df, !grepl("Contaminant_", Protein))
  # df <- filter(df, !grepl("XXX_", Protein))
  return(df)
}
### psm_filter----- END---------



# plot for QC------ delppm and pepQ value----------
plot_delm <- function(df) {
  hist(df$DelM_PPM,
       xlab = "Mass error (ppm)", main = "",
       font.axis = 2, font.lab = 2)
}

plot_pepQ <- function(df) {
  hist(df$PepQValue,
       xlab = "PepQValue", main = "",
       font.axis = 2, font.lab = 2)
}

# plot for QC------ delppm and pepQ value----------

#####################  plot for FDR  #####################################
plot_FDR <- function(df,PepQValue_cutoff) {
  
  ppm_center <- median(df$DelM_PPM)
  ppm_left <- ppm_center - 10
  ppm_right <- ppm_center + 10
  
  df <- filter(df,DelM_PPM < ppm_right & DelM_PPM > ppm_left)
  
  df <- filter(df, PepQValue < PepQValue_cutoff)
  
  psm <- as.numeric(table(grepl("XXX_",df$Protein)))
  
  Decoy <- psm[2]
  True <- psm[1]
  Total <- sum(psm)
  fdr <- round(100* 2 * Decoy / Total, digits = 2)
  
  title <- paste(paste("PSM IDs (",
                       paste(paste("FDR =", fdr),
                             "%", sep = ""), sep = ""),
                 ")", sep = "")
  
  num <- c(Decoy,True)
  type <- c( "Decoy", "True")
  df <- data.frame(type = type, num = num)
  
  ylim <- 1.1*True
  
  p1 <- ggplot(df, aes(type, num)) +
    geom_bar(stat = "identity", width = 0.5, 
             color = "black", fill = c("red", "green"))+
    scale_y_continuous(limits = c(0, ylim))+
    xlab("") +
    ylab("Numer of PSMs") +
    geom_text(aes(label = num), 
              position = position_dodge(width = 0.6),
              size = 4, fontface = 2,
              vjust = -0.6, angle = 0,  hjust = 0.5) +
    theme_bw() +
    theme(axis.title.y = element_text(size = 11, colour = "black", face = "bold"),
          axis.text.x = element_text(size = 11, colour = "black", face = "bold"),
          axis.text.y = element_text(size = 10, colour = "black", face = "bold"),
          legend.title = element_blank()) +
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p1)
}

#####################  plot for FDR  #####################################

#####################  enrichemtn effiency ##########################

plot_Cys_spec <- function(dc) {
  
  
  dc$clean_pep <- str_sub(dc$Peptide,3,str_length(dc$Peptide)-2)
  pn <- as.numeric(table(grepl("C", dc$clean_pep)))
  total_pep <- sum(pn)
  non_cys <- pn[1]
  cys_pep <- pn[2]
  Enrich_eff <- round(100 *cys_pep / total_pep, digits = 2)
  
  title <- paste("Enrichment Specificity:",
                 paste(Enrich_eff, "%", sep = ""),
                 sep = " ")
  
  num <- c(non_cys,cys_pep)
  type <- c( "Non-Cys", "Cys")
  df <- data.frame(type = type, num = num)
  
  x_limits <- c( "Non-Cys", "Cys")
  ylim <- 1.1*cys_pep
  
  p1 <- ggplot(df, aes(type, num)) +
    geom_bar(stat = "identity", width = 0.5, 
             color = "black", fill = c("red", "green"))+
    scale_x_discrete(limits = x_limits)+
    scale_y_continuous(limits = c(0, ylim))+
    xlab("") +
    ylab("Numer of PSMs") +
    geom_text(aes(label = num), 
              position = position_dodge(width = 0.6),
              size = 4, fontface = 2,
              vjust = -0.6, angle = 0,  hjust = 0.5) +
    theme_bw() +
    theme(axis.title.y = element_text(size = 11, colour = "black", face = "bold"),
          axis.text.x = element_text(size = 11, colour = "black", face = "bold"),
          axis.text.y = element_text(size = 10, colour = "black", face = "bold"),
          legend.title = element_blank()) +
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p1)
}

#####################   enrichemtn effiency ##########################


# aggregate at the peptide level-------- START-----------
Aggregate_to_Uniqe_Peptide <- function(dm) {
  dm2 <- select(dm, Peptide, Protein, Ion_126.128 : Ion_131.138 )
  
  dm2 <- dm2 %>% 
    group_by(Peptide) %>% 
    summarise(Protein = first(Protein),    # protein by the first protein
              TMT1 = sum(Ion_126.128),     # for PSMs belongint to the same peptide, add.
              TMT2 = sum(Ion_127.125),
              TMT3 = sum(Ion_127.131),
              TMT4 = sum(Ion_128.128),
              
              TMT5 = sum(Ion_128.134),
              TMT6 = sum(Ion_129.131),
              TMT7 = sum(Ion_129.138),
              TMT8 = sum(Ion_130.141),
              
              TMT9 = sum(Ion_131.138),
              PSM_cnt = n())             # PSM count 
  
  d_comp <- dm2[rowSums(dm2[,3:11] == 0) == 0, ]
  
  d_comp[,3:11] <- log(d_comp[,3:11],2) 
  
  return(d_comp)
}
# aggregate at the peptide level-------- END -----------

# Remove decoy hits, contamination, and non-cys peps-------- START-----------
Remove_Decoy_Contam_NonCys <- function(df) {
  df <- filter(df, !grepl("Contaminant", df$Protein))  # remove Contaminant peptide
  df <- filter(df, !grepl("XXX_", df$Protein))   # remove decoy hit
  df$clean_pep1 <- str_sub(df$Peptide,3,str_length(df$Peptide)-2)
  df <- filter(df, grepl("C", df$clean_pep1))  # remove non-Cys peptide
  return(df)
}
# Remove decoy hits, contamination, and non-cys peps-------- START-----------


# within-group normalization-------------START--------------
med_norm <- function(df) {
  norm.coeff <- apply(df, 2, median)  # calculate the median for each channel
  df1 <- sweep( df, 2, norm.coeff, '-')  # bring the median of each channel to 0.
  avg_of_median <- mean(norm.coeff )   # the new median for each channel.
  df1 <- df1 +  avg_of_median
  return(df1)
}

Norm_Within_Group <- function(df) {
  df[,3:6] <- med_norm(df[,3:6])  
  df[,7:10] <- med_norm(df[,7:10])
  return(df)
}
# within-group normalization-------------END--------------



# TMT intensity distribution plot-------------START--------------
TMT_Int_Distribution <- function(df1, df2) {
  col_set1 <- c(rep("grey", 4),
                rep("darkolivegreen3", 4),  
                "deepskyblue2")
  par(mfrow = c(1,2))
  
  med_line <- mean(round(apply(df1[,3:10], 2,median),2))
  boxplot(df1[,3:11],col = col_set1, las = 2,
          cex.axis = 1.1, font.axis = 2,
          main = "TMT Int.Distribution: Before Normalization",
          ylab = "log2 Intensity",
          cex.lab = 1.2, font.lab =2)
  abline(h = med_line, col = "Red", lty = 5, lwd = 3)
  
  med_line <- mean(round(apply(df2[,3:10], 2,median),2))
  boxplot(df2[,3:11],col = col_set1, las = 2,
          cex.axis = 1.1, font.axis = 2,
          main = "TMT Int.Distribution: After Normalization",
          ylab = "log2 Intensity",
          cex.lab = 1.2, font.lab =2)
  abline(h = med_line, col = "Red", lty = 5, lwd = 3)
  
}
# TMT intensity distribution plot-------------START--------------



# calculat p value, log2FC, occupancy and sd.-------START----------------------------------------

# a function to calculate for indivudal occ: occ_indi
occ_i <- function(x,y) {
  2^(x-y)*100
}

# a fucntion to calculate sd for each group
occ_sd <- function(x1,x2,x3,x4,y){
  occ1 <- occ_i(x1,y)
  occ2 <- occ_i(x2,y)
  occ3 <- occ_i(x3,y)
  occ4 <- occ_i(x4,y)
  
  sdx <- sd(c(occ1,occ2, occ3, occ4))
  
  return(sdx)
}

### main funciton
Stats_Calc <- function(df){
  df$p.value <- apply(df[,3:10], 1, function(x) t.test(x[1:4],x[5:8])$p.value)
  df$log2FC <- apply(df[,3:10], 1, function(x) mean(x[5:8]) - mean(x[1:4]))
  df$adjusted_p <- p.adjust(df$p.value, method = "BH")
  df$occ_1 <- apply(df[,3:11], 1, function(x) 100*2^(mean(x[1:4])- x[9]))
  df$occ_2 <- apply(df[,3:11], 1, function(x) 100*2^(mean(x[5:8])- x[9]))
  
  # add the sd info by the apply funciton
  df$sd1 <- apply(df[,c(3:6,11)], 1, function(x) occ_sd(x[1], x[2], x[3], x[4], x[5]))
  df$sd2 <- apply(df[,c(7:10,11)], 1, function(x) occ_sd(x[1], x[2], x[3], x[4], x[5]))
  return(df)
}
# calculat p value, log2FC, occupancy and sd.-------END----------------------------------------


#  volcano plot with denisty info------------------------------------------START---------------
vplot_q <- function(df,mymain) {
  x1 <- df$log2FC
  x2 <- -log(df$adjusted_p,10)
  df <- data.frame(x1,x2)
  
  ## Use densCols() output to get density at each point
  x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  
  ## Map densities to colors
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                              "#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  
  ## Plot it, reordering rows so that densest points are plotted on top
  
  par(mar = c(6,6,6,6))
  
  plot(x2~x1, data=df[order(df$dens),], 
       pch = 20, col=col, cex = 0.3,
       xlab = "log2 (Fold Change)",
       ylab = "-log10 (adj.p value)",
       cex.lab = 1.8, cex.axis = 1.5,
       main = mymain, cex.main = 2)
  
  n = dim(df)[1]
  usr <- par("usr")
  x_pos <- usr[2] * 0.7
  y_pos <- usr[4] * 0.95
  text(x_pos,y_pos,paste("n = ",n, sep = ""), cex = 1.5, font = 2)
  box(lwd = 2)
}

# with original p vlaue
vplot_p <- function(df,mymain) {
  x1 <- df$log2FC
  x2 <- -log(df$p.value,10)
  df <- data.frame(x1,x2)
  
  ## Use densCols() output to get density at each point
  x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  
  ## Map densities to colors
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                              "#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  
  ## Plot it, reordering rows so that densest points are plotted on top
  
  par(mar = c(6,6,6,6))
  
  plot(x2~x1, data=df[order(df$dens),], 
       pch = 20, col=col, cex = 0.3,
       xlab = "log2 (Fold Change)",
       ylab = "-log10 (p value)",
       cex.lab = 1.8, cex.axis = 1.5,
       main = mymain, cex.main = 2)
  
  n = dim(df)[1]
  usr <- par("usr")
  x_pos <- usr[2] * 0.7
  y_pos <- usr[4] * 0.95
  text(x_pos,y_pos,paste("n = ",n, sep = ""), cex = 1.5, font = 2)
  box(lwd = 2)
}
#  volcano plot with denisty info------------------------------------------END---------------

################ p value distribution ###############-----------START--------
plot_p_Dist <- function(df) {
  hist(df$p.value, breaks = seq(0,1,0.01),
       xlab = "p value", main = "",
       las = 1,
       col = c("blue", rep("white",99)))
}
################ p value distribution ###############-----------END--------


################ occupancy distribution ###############-----------START--------
occ_Dist <- function(df){
  s1 <- filter(df, occ_1 < 100, occ_2 < 100 )
  occ <- as.numeric(s1$occ_1)
  
  hist(occ, breaks = 100,
       main = "", xlab = "",
       cex.lab = 1, font.lab = 1, col = "darkolivegreen3",
       xlim = c(0, 110),
      
       yaxt = "n")
  abline(v = median(occ), col = "red", lty = 2, lwd =2)
  usr <- par("usr")
  x_range <- usr[2] -usr[1]
  x_pos <- usr[1] + 0.4*x_range
  
  y_range <- usr[4] -usr[2]
  y_pos <- usr[3] + 0.9*y_range
  
  txt <- paste("Median = ",
               paste(round(median(occ),2), "%", sep = ""),
               sep = "")
  text(x_pos, y_pos,txt, cex = 1, font = 1, col = "red" )
  axis(side = 2,las = 2,
       mgp = c(3, 0.75, 0))
}
################ occupancy distribution ###############-----------END--------


# occ correlation plot-------------====================###################### START 
plot_Occ_Corr <- function(df) {
  s1 <- filter(df, occ_1 < 100, occ_2 < 100 )
  colnames(s1)
  x1 <- s1$occ_1
  x2 <- s1$occ_2
  
  df <- data.frame(x1,x2)
  
  ## Use densCols() output to get density at each point
  x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  
  ## Map densities to colors
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                              "#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  
  ## Plot it, reordering rows so that densest points are plotted on top
  plot(x2~x1, data=df[order(df$dens),], 
       xlim = c(0,100),
       ylim = c(0,100),
       pch = 20, col=col, cex = 0.2)
  abline(0,1, col = "red", lty = 2, lwd = 2)
}
# occ correlation plot-------------====================######################  END


######## CYs site location ##############------------START----------------------------------------------
# the starting location of peptide on protein----fucntion 1 within Cys_Site_Location-----------------START-----------------------
# function pept_loc >>>>>>>>>>>
pept_loc <- function(pep_seq, protein_seq){
  pep_seq <- as.character(pep_seq) 
  protein_seq <- as.character(protein_seq)
  
  clean_pep1 <- str_sub(pep_seq,3,str_length(pep_seq)-2) #remove the previous aa residue and the dot
  clean_pep2 <- gsub("[^A-Z]","",clean_pep1) #remove the modification sign
  
  pept_loc1 <- str_locate(protein_seq,clean_pep2) # result is a list
  pept_loc2 <- unlist(pept_loc1[1])[1] 
  return(pept_loc2)
}
# End <<<<<<<<<  function pept_loc
# the starting location of peptide on protein----fucntion 1 within Cys_Site_Location-----------------END-----------------------


# function mod_loc---------------------fucntion 2 within Cys_Site_Location-------------START---------------------
# -------------------------------------------------------------
# modification list-----------------------
#	-	229.162933	K	S	TMT6Tag	74307
#	-	229.162933	<	T	TMT6Tag	40768
#	*	15.994915	M	D	Plus1Oxy	11502
#	#	125.047676	C	D	NEM	17866  # Note in this search, there is no NEM modification
# modification list-----------------------
# first: remove any -, *, so peptide only have # mod
# second: turn C# (NEM modififed cys) into c (little c), now the whole sequence do not have any mod # no need any more
# third: get the location of free thiol on the sequence
# here all Cys are supposed to be free, and then labeld with IAA

# count the position of modification on peptide

mod_loc <- function(x){
  x <- as.character(x) 
  
  clean_pep <- str_sub(x,3,str_length(x)-2) # remove head and tail;
  clean_pep1 <- gsub("\\*|-","",clean_pep) # remove * and -
  pep_C_NEM_small <- gsub("C#","c",clean_pep1)
  
  star_pos1 <- str_locate_all(pep_C_NEM_small ,"C") # location of free thiols
  star_pos2 <- star_pos1[[1]]
  mod_all <- as.numeric( star_pos2[,1])
  mod_all <- paste(mod_all, collapse= ",")
  return(mod_all)
}

# function mod_loc---------------------fucntion 2 within Cys_Site_Location----------------END------------------


# sum_non_NA -----------------------FUNCTION 3 within Cys_Site_Location-----------------START---------------------
# Now a function to calculate the position of modification on protein level
# function sum_non_NA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# x will be the position of peptide on protein
# y will be the modificaiton on peptide (it could be two numbers)

sum_non_NA <- function(x,y) {
  x <- as.numeric(x)
  if(is.na(y)) {pos = NA}
  
  if(!is.na(y)) {  # only one number
    if(!grepl(",",y)) {
      y <- as.numeric(y)
      pos = x + y -1
    }
    
    if(grepl(",",y)) {  # more than one number
      y <- as.numeric(unlist(strsplit(y,",")))
      pos = x + y -1
      pos <- do.call(paste, c(as.list(pos), sep = ";"))
    }}
  return(pos)
}
# sum_non_NA -----------------------FUNCTION 3 within Cys_Site_Location-----------------END---------------------

Cys_Site_Location <- function(df,db) {
  df1 <- select(df, Peptide,Protein) # only peptide and protein are needed.
  df1$id <- 1:dim(df1)[1] # need to keep df2 in the order of df1
  df2 <- merge(df1,db, 
               by.x = "Protein", by.y = "ID",
               type = "left")
  df2 <- df2[order(df2$id),] # reorder df2
  df2 <- select(df2, -id)
  
  
  # step 1: determine the start position of peptide on proteins
  df2$pep_on_pro_loc <- apply(df2[,c(2,4)], 1, function(x) pept_loc(x[1], x[2]))
  
  
  # step 2: determine the positions of free -SH (initially oxidized) on peptide--------
  Pep_mod <- NULL
  for (i in 1:dim(df2)[1]) {
    pep_mod_pos <- mod_loc(df2[i,2])
    Pep_mod <- c(Pep_mod,pep_mod_pos)
  }
  df2$Cys_on_pep_loc <- Pep_mod
  
  # step 3: Cys mod on protein level
  df2$Cys_on_pro_loc <- apply(df2[,5:6],1,function(x) sum_non_NA(x[1],x[2]))
  
  # step 4: add localization and protein description back to the origindal df
  df$pep_on_pro_loc <- df2$pep_on_pro_loc
  df$Cys_on_pep_loc <- df2$Cys_on_pep_loc
  df$Cys_on_pro_loc <- df2$Cys_on_pro_loc
  
return(df)
  
}

######## CYs site location ##############------------END------------------------------------------------

# add protein annoation info: gene symbol, subcellular loc, etc######----START-----
Protein_Annotation <- function(df1, df2){
  df3 <- merge(df1, df2, 
               by.x = "Protein", by.y = "Protein_ID", 
               all.x = TRUE)
  return(df3)
}
# add protein annoation info: gene symbol, subcellular loc, etc######----END-----



############### Aggregate to Cys site--------##################----------------------------START-----
Aggregate_to_Cys_Site <- function(df){
  df$Cys_ID <- paste(df$Protein, df$Cys_on_pro_loc, sep = "_Cys")
  
  # create dataset for merge
  dm <- select(df,
               Cys_ID,Protein,Cys_on_pep_loc,
               Entry: subcellular_localisation)
  dm <- dm[!duplicated(dm$Cys_ID),]
  
  # select Cys_ID, and TMT quan data
  d1 <- df[,c(28,3:11)]
  
  d1[,2:10] <- 2^d1[,2:10] # get back to raw intensity
  d2 <-  aggregate(.~ Cys_ID, d1,sum) # aggregation
  d2[,2:10] <- log(d2[,2:10],2) # back to log2
  
  d2$p.value <- apply(d2[,2:9], 1, function(x) t.test(x[1:4],x[5:8])$p.value)
  d2$log2FC <- apply(d2[,2:9], 1, function(x) mean(x[5:8]) - mean(x[1:4]))
  d2$adjusted_p <- p.adjust(d2$p.value, method = "BH")
  d2$occ_1 <- apply(d2[,2:10], 1, function(x) 100*2^(mean(x[1:4])- x[9]))
  d2$occ_2 <- apply(d2[,2:10], 1, function(x) 100*2^(mean(x[5:8])- x[9]))
  
  colnames(d2)
  # add the sd info by the apply funciton
  d2$sd1 <- apply(d2[,c(2:5,10)], 1, function(x) occ_sd(x[1], x[2], x[3], x[4], x[5]))
  d2$sd2 <- apply(d2[,c(6:9,10)], 1, function(x) occ_sd(x[1], x[2], x[3], x[4], x[5]))
  
  
  d3 <- merge(d2,dm, by = "Cys_ID")
  return(d3)
  
}

############### Aggregate to Cys site--------##################----------------------------END-----
