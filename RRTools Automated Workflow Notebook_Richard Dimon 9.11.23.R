# library(pacman)
# pacman::p_load(ggsn, ggspatial, colorRamps,ggtree,tanggle,png, spatstat,ICSNP,openxlsx,ggmap, ade4, adegenet, animation, ape, car, caTools, cowplot, dartR, devtools, diveRsity, dplyr, gganimate, ggh4x, ggplot2, ggrepel, grid, gridExtra, LEA, magick, mapplots, OptGenMix, oz, ozmaps, plyr, poppr, readxl, reshape2, rgl, RRtools, sfsCalcs, tidyr, png, vegan, webshot2, geosphere, HWxtest, phangorn, phytools)
load("C:/Users/dimonr/OneDrive - DPIE/R/Scripts/PSFsaved_scripts4autom.R")
maindir<- "E:/rrspecies/"
setwd(maindir)
RandRbase <- ""
spls <- read.csv(paste0(maindir,"RRSpeciesParameters.csv"),header = TRUE)
nsp <- length(spls$species) #number of species

# for (z in 15:45) {
z=47

##### Import and filter data - with or without QC stats #####

species <- spls$species[z]
dataset <- spls$dataset[z]
analysis <- spls$analysis[z]
topskip <- spls$topskip[z]
nmetavar <- spls$nmetavar[z]
max_missing <- spls$max_missing[z]
min_repro <- spls$min_repro[z]
print(paste0("Starting ", species, analysis))


source("C:/Users/dimonr/OneDrive - DPIE/R/Packages/RRtools/R/report.dart.qc.stats.R")

# create an output directory
outputloc <- paste0(maindir, species, "/output_", analysis)
outputloc2 <- paste0(maindir, species, "/output_", analysis, "/")
if (!dir.exists(outputloc)) {
  dir.create(outputloc)
}


qcdirectory <- paste0(outputloc, "/qual_stat/") #create a qual_stat file inside each output folder
if (!dir.exists(qcdirectory)) {
  dir.create(qcdirectory)
}

#with QC analyses
# d1        <- read.dart.xls.onerow(RandRbase,species,dataset,topskip, nmetavar, euchits=FALSE, altcount=TRUE)
# missing <- is.na(d1$gt)
# count_of_missing_by_locus <- rowSums(missing)
# count_of_missing_by_sample <- colSums(missing)
# qc1       <- report.dart.qc.stats(d1, outputloc, species, dataset, threshold_missing_loci = 0.8)
# d2        <- remove.poor.quality.snps(d1, min_repro=min_repro, max_missing=max_missing)
# qc2       <- report.dart.qc.stats(d2, outputloc, species, dataset)
# d3        <- sample.one.snp.per.locus.random(d2, seed=12345)
# qc3       <- report.dart.qc.stats(d3, outputloc, species, dataset)
# d3_rm <-  remove.duplicate.samples(d3, least_missing=TRUE, remove_fixed_loci=TRUE)
# d3 <- d3_rm


#without QC analyses
d1        <- read.dart.xls.onerow(RandRbase,species,dataset,topskip, nmetavar, euchits=FALSE, altcount=TRUE)
missing <- is.na(d1$gt)
count_of_missing_by_locus <- rowSums(missing)
count_of_missing_by_sample <- colSums(missing)
d2        <- remove.poor.quality.snps(d1, min_repro=min_repro, max_missing=max_missing)
d3        <- sample.one.snp.per.locus.random(d2, seed=12345)
d3_rm <-  remove.duplicate.samples(d3, least_missing=TRUE, remove_fixed_loci=TRUE)
d3 <- d3_rm

#read metadata file
metafile <- read.xlsx(paste0(maindir, species,"/meta/",species,"_",dataset,"_meta.xlsx"))
m1 <- read.meta.data(d3, RandRbase, species, dataset, fields = as.numeric(ncol(metafile)-4))
colnames(m1$analyses)
dm <- dart.meta.data.merge(d3, m1)
fields <- c(analysis)
dms <- data.by.meta.fields(dm, fields, RandRbase, species, dataset, object = analysis)

# correct inconsistencies in J's scripts for naming treatment which is for naming files
treatment <- paste0("raw_SNPFilt_1SNPperClone_Field_", analysis)
dms$treatment <- treatment
gta <- dms$gt

#qualStatDMS
# qcdms <- report.dart.qc.stats(dms, outputloc, species, dataset)
# summarise_QC(species, RandRbase, dataset) #deafult Dart QC summary figures

# 
# #assign the same shapes for both splitstree and PCA
# shapeassign <- c()
# shapenum <- c()
# shapenum <- c(7, 8, 9, 10, 12, 14, 15, 16, 17, 18, 23, 25) #my fave shapes!
# shps <- function(x, n){
#   sample(c(x, sample(shapenum, n-length(x), replace=TRUE)))
# }
# shapeassign <- shps(shapenum, length(unique(dms$meta$analyses[,analysis])))
# 
# 

##### Migrate Analysis limiting 5 samples per pop #####

gst_dir <- paste0(outputloc, "/gst_migrate/")
if (!dir.exists(gst_dir)) {
  dir.create(gst_dir)
}

source("C:/Users/dimonr/OneDrive - DPIE/R/Packages/sos_functions.R")



m1 <- read.meta.data(d3, RandRbase, species, dataset, fields = as.numeric(ncol(metafile)-4))
temp<-as.data.frame(m1$analyses)

#set a sample threshold for each pop:
samplethreshold <-5# chnage this to what number of samples you want toi have uniform across all pops



#### Below lines remove sites with less than threshold number of samples, and subsample sites that are above the threshold
for (x in unique(na.omit(temp[, analysis]))) {
  pop_samples <- which(temp[, analysis] == x)
  if (length(pop_samples) < samplethreshold) {
    temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
    print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
  }
  if (length(pop_samples) == samplethreshold) {
    print(paste0("Keeping ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
  }
  if (length(pop_samples) > samplethreshold) {
    subsamplepops <- pop_samples[sample(1:length(pop_samples), size=samplethreshold, replace = F)]
    for (b in 1:length(pop_samples)){
      if (pop_samples[b] %in% subsamplepops == FALSE){
        temp[, analysis][pop_samples[b]] <- NA
      }
      print(paste0("Reducing ",x, " to a sample size of ", samplethreshold," samples"))
    }
  } else {}
}




# ## If you want to just remove any sites that are below the threshold (but keep individuals that are above a threshold, then unhash the below section)
# for (x in unique(na.omit(temp[, analysis]))) {
#   pop_samples <- which(temp[, analysis] == x)
#   if (length(pop_samples) < samplethreshold) {
#     temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
#     print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
#   } else {}
# }

m1$analyses<-temp
dm        <- dart.meta.data.merge(d3, m1)
fields    <- c(analysis)
dms     <- data.by.meta.fields(dm, fields, RandRbase, species, dataset, object=analysis); print(dms$meta$analyses[,analysis])

#correct inconsistencies in J's scripts for naming treatment which is for naming files
treatment <- paste0("raw_SNPFilt_1SNPperClone_Field_", analysis)
dms$treatment <- treatment
gta <- dms$gt

# c <- make_genepop_file(dms_pf, maf=0.05, missing=0.2, group=paste0(species, analysis), grouping=dms_pf$meta$analyses[,analysis])

c   <- dart2genepop(dms, RandRbase, species, dataset, dms$meta$analyses[,analysis])


# Diversity was measured as mean allelic richness for every population using the resampling method in the divBasic function in the R package diveRsity (46); to convert these diversity values into the data structure of pairwise relationships between populations, we calculated genetic diversity ratios for each population pair in each direction as the ratio of allelic richness in the destination versus the origin population. Pairwise genetic similarity was measured as 1/Fst.

# diveRsity::divBasic()
# AR destination/AR origin for each population 
# pairwise genetic similarity 1/fst 


#Migration was estimated using the divMigrate method (19) with Jost’s D metric of differentiation (47), as implemented in the R package diveRsity (46). This approach uses allele frequency differences between population pairs to estimate rates of migration in each direction; note that these rates are relative to other population pairs in the same data set and cannot be compared across data sets.

v <- diveRsity::divMigrate(infile=c, outfile=NULL, stat="gst",plot_network=TRUE, filter_threshold = 0.2, boots=1000, para=TRUE)



d_mig <- v$gRelMig #v$dRelMig

# from is rows, to is columns
colnames(d_mig) <- unique(dms$meta$analyses[,analysis]) 
rownames(d_mig) <- unique(dms$meta$analyses[,analysis]) 

# Filter by significance -- significance is if there is a significant difference in the directions
# Test for overlap of the estimated 95% confidence intervals. Where there is no overlap, the directional gene flow components are said to be significantly different (asymmetric).
mig_sig <- v$gRelMigSig #v$dRelMigSig
colnames(mig_sig) <- unique(dms$meta$analyses[,analysis])
rownames(mig_sig) <- unique(dms$meta$analyses[,analysis])

d_mig[mig_sig>0.05] <- NA

# qgraph::qgraph(d_mig,legend = F, edge.labels = F,
# curve = 2.5, mar = c(2, 2, 5, 5))

long_mig <- melt(d_mig)

# meta_agg <- m2 %>%
#   group_by(pop_large_short,pop_large) %>%
#   summarize(lat = mean(lat, na.rm=TRUE),
#             long = mean(long,na.rm=TRUE),
#             .groups = 'drop')%>%
#   subset(.,pop_large_short!="Ex_situ_PF")



dms_df <- data.frame(cbind(dms$meta$analyses[,analysis],dms$meta$lat, dms$meta$long))
colnames(dms_df) <- c("ansites", "lat", "long")

dms_df$lat <- as.numeric(dms_df$lat)
dms_df$long <- as.numeric(dms_df$long)


meta_agg <- dms_df %>% 
  group_by(dms_df$ansites) %>% 
  dplyr::summarize(lat = mean(lat, na.rm=TRUE),
                   long = mean(long, na.rm=TRUE), 
                   .groups = 'drop')

meta_agg <- na.omit(meta_agg)


colnames(meta_agg) <- c("gsites", "lat", "long")

long_mig <- merge(long_mig, distinct(meta_agg[, c("gsites","lat","long")]), by.x = "Var1", by.y = "gsites", all.y = FALSE)
long_mig <- merge(long_mig, distinct(meta_agg[, c("gsites","lat","long")]), by.x = "Var2", by.y = "gsites", all.y = FALSE)
long_mig<- distinct(long_mig)  

colnames(long_mig)[1:3] <- c("from", "to", "gst")

# long_mig <- long_mig[long_mig$from!="Ex_situ_PF"&long_mig$to!="Ex_situ_PF"&!is.na(long_mig$gst),]

long_mig <- long_mig[!is.na(long_mig$gst),] # remove any NA values between the same sites



# 
# map <- ggmap(get_map(location = bound, source = "stamen",
#                      zoom = 14, scale = 4,
#                      maptype = "terrain",
#                      color = "bw"))
# ggplot()+
# gst_map <-map+coord_cartesian()+coord_fixed()+
#   theme_bw()+
#   geom_curve(data=long_mig[long_mig$gst>0.3,],#
#                aes(x=long.x, y=lat.x,
#                    xend = long.y, yend = lat.y, colour=gst), #alpha=0.7+gst,colour=gst, 
#               size = 0.9,na.rm = TRUE,curvature=0.3,#lineend = "round",
#              # curvature=-0.3,,
#                arrow=arrow(angle=20, ends="first",type="closed", length=unit(2, "mm")))+
#   scale_color_gradient(low = "white", high = "red")+ #midpoint=0.5, mid="blue",
#   # geom_text(data=long_mig[long_mig$gst>=0.27,],
#   #           aes(label = round(gst, 2), x = (long.x + long.y) / 2, y = (lat.y+lat.x) / 2),
#   #           size = 2,
#   #           color = "black")+
#   geom_point(data=meta_agg, mapping=aes(x=long, y=lat), colour="black")+
#   xlim(lims[1], lims[2])+
#   ylim(lims[3], lims[4])+   labs(x = "Longitude", y = "Latitude", colour="Gst") +guides(size = "none", alpha="none")+
#   ggrepel::geom_label_repel(data = meta_agg,aes(x = long, y = lat, label=pop_large_short),
#                             min.segment.length=0.25, color="black",fill="white",size=3, segment.colour="white",
#                             alpha=0.9, label.size=0, nudge_y=0.003)
# 
# 
# ggsave("PherFitz/outputs/plots/gst_map.png", plot = gst_map, width = 200, height = 110, dpi = 300, units = "mm")
# 










#Mapping
row_sub = apply(data.frame(dms$meta$long), 1, function(row) all(row !=0 ))
newdmslong <- data.frame(dms$meta$long)[row_sub,]
row_sub = apply(data.frame(dms$meta$lat), 1, function(row) all(row !=0 ))
newdmslat <- data.frame(dms$meta$lat)[row_sub,]

#this bit looks for range of latlong values to zoom in on and calulcates range for zoom for ggmap
divxlims <- c(min(na.omit(newdmslong))-1,max(na.omit(newdmslong))+1) #default for large distribution s around 0.6
divylims <- c(min(na.omit(newdmslat))-1,max(na.omit(newdmslat))+1) #default for large distribution s around 0.4
divxlimsrange <- abs(c(min(na.omit(newdmslong))-max(na.omit(newdmslong))))
divylimsrange <- abs(c(min(na.omit(newdmslat))-max(na.omit(newdmslat))))

if(divylimsrange > 5) {gmapzoom= 4}
if(divylimsrange > 4 && divylimsrange < 5)  {gmapzoom= 5}
if(divylimsrange > 3 && divylimsrange < 4)  {gmapzoom= 8}
if(divylimsrange > 2 && divylimsrange < 3)  {gmapzoom= 9}
if(divylimsrange > 1 && divylimsrange < 2)  {gmapzoom= 10}
if(divylimsrange > 0.2 && divylimsrange < 1) {gmapzoom= 11}
if(divylimsrange < 0.2) {gmapzoom= 12}

library(ggmap)
# register_google(key = "AIzaSyBZIYV071yUULZTtcUMUZTMwOgRJrKAcaw")
# 
# googlemapbase = get_map(location = c(left = divxlims[1], bottom = divylims[1], right = divxlims[2], top=divylims[2]), zoom = gmapzoom, scale = "auto", maptype = "satellite", source = c("google"))
# 

register_stadiamaps(key = "4fe26ca2-aa51-4e9e-8151-dca93798847e")

stadiamapbase = get_stadiamap(bbox = c(left = divxlims[1], bottom = divylims[1], right = divxlims[2], top=divylims[2]), zoom = 9, maptype = "stamen_terrain", crop=F)


#import GHA shapefiles
library(sf)

BR <- read_sf("C:/Users/dimonr/OneDrive - DPIE/R/Projects/GHA/Border Ranges GHA 10km Buffer.shp")
WP <- read_sf("C:/Users/dimonr/OneDrive - DPIE/R/Projects/GHA/Washpool GHA 10km Buffer.shx")
IL <- read_sf("C:/Users/dimonr/OneDrive - DPIE/R/Projects/GHA/Illuka GHA 10km Buffer.shp")
DR <- read_sf("C:/Users/dimonr/OneDrive - DPIE/R/Projects/GHA/Dorrigo GHA 10km Buffer.shp")
WK <- read_sf("C:/Users/dimonr/OneDrive - DPIE/R/Projects/GHA/Werrikimbe GHA 10km Buffer.shp")
BA <- read_sf("C:/Users/dimonr/OneDrive - DPIE/R/Projects/GHA/Barrington GHA 10km Buffer.shp")
DRWK <- read_sf("C:/Users/dimonr/OneDrive - DPIE/R/Projects/GHA/DorrigoANDWerrikimbe GHA 10km Buffer.shp")

BR <- as_Spatial(BR)
WP <- as_Spatial(WP)
IL <- as_Spatial(IL)
DR <- as_Spatial(DR)
WK <- as_Spatial(WK)
BA <- as_Spatial(BA)
DRWK <- as_Spatial(DRWK)


library(ggthemes)

long_mig <- long_mig[order(long_mig$gst),]


ggmap(stadiamapbase) + theme_few()+
  coord_cartesian()+ coord_fixed()+
  geom_polygon(data=BR, aes(long, lat, group=group),colour = alpha("darkred", 0.5), linewidth = 1,fill="darkred", alpha = 0.2)+
  geom_polygon(data=WP, aes(long, lat, group=group),colour = alpha("darkorange", 0.5), linewidth = 1, fill="darkorange", alpha = 0.2)+
  # geom_polygon(data=IL, aes(long, lat, group=group),colour = alpha("pink", 0.5), linewidth = 1, alpha = 0.3)+
  # geom_polygon(data=DRWK, aes(long, lat, group=group),colour = alpha("cyan", 0.5), linewidth = 1, fill="cyan", alpha = 0.2)+
  geom_polygon(data=DR, aes(long, lat, group=group),colour = alpha("cyan", 0.5), linewidth = 1, fill="cyan", alpha = 0.2)+
  geom_polygon(data=WK, aes(long, lat, group=group),colour = alpha("green", 0.5), linewidth = 1, fill="green", alpha = 0.2)+
  geom_polygon(data=BA, aes(long, lat, group=group),colour = alpha("yellow", 0.5), linewidth = 1, fill="yellow", alpha = 0.2)+
  geom_curve(data=long_mig[long_mig$gst>0,],#
             aes(x=long.x, y=lat.x,
                 xend = long.y, yend = lat.y, colour=gst), #alpha=0.7+gst,colour=gst, 
             linewidth = 3,na.rm = TRUE,curvature=0.3, alpha=1, #lineend = "round",
             arrow=arrow(angle=20, ends="first",type="closed", length=unit(5, "mm")))+
  scale_color_gradient(low = "white", high = "darkblue")+ #midpoint=0.5, mid="blue",
  # geom_text(data=long_mig[long_mig$gst>=0.27,],
  #           aes(label = round(gst, 2), x = (long.x + long.y) / 2, y = (lat.y+lat.x) / 2),
  #           size = 2,
  #           color = "black")+
  geom_point(data=meta_agg, mapping=aes(x=long, y=lat), colour="black", size=2)+
  labs(x = "Longitude", y = "Latitude", colour="GST") +guides(size = "none", alpha="none")+
  # ggsn::scalebar(data=meta_agg, dist = 25, dist_unit = "km", location = "bottomleft", 
  #              st.bottom = T, st.size = 5, st.dist = 0.01, border.size = 0.1,
  #              transform = TRUE, model = "WGS84", height = 0.015, anchor = c(x = min(meta_agg$long)-0.3, y = min(meta_agg$lat)-0.3)) +
  xlim(150.8, 153.9)+ # increase to change width of plot
  ylim(-32.5, -27.7)+
  annotation_north_arrow(location = "br", which_north = "true", pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"), style = north_arrow_fancy_orienteering)+
  ggtitle(paste0(species," ",analysis, " DivMigrate"))

# ggrepel::geom_label_repel(data = meta_agg,aes(x = long, y = lat, label=gsites),
#                           min.segment.length=0.25, color="black",fill="white",size=1, segment.colour="white",
#                           alpha=0.9, label.size=0, nudge_y=0.003)


# gst_fst <- ggarrange(fst_map, gst_no_map, nrow=2, labels=c("A","B"), align="hv")

ggsave(paste0(species,"_", analysis, "gst_Map.tiff"), path = paste0(gst_dir), width = 10, height = 12, dpi = 300, units = "in")













##### FST, IBD and distance map with limiting X samples per pop #####

m1 <- read.meta.data(d3, RandRbase, species, dataset, fields = as.numeric(ncol(metafile)-4))
temp<-as.data.frame(m1$analyses)

#set a sample threshold for each pop:
samplethreshold <-5# change this to what number of samples you want toi have uniform across all pops


#### Below lines remove sites with less than threshold number of samples, and subsample sites that are above the threshold
for (x in unique(na.omit(temp[, analysis]))) {
  pop_samples <- which(temp[, analysis] == x)
  if (length(pop_samples) < samplethreshold) {
    temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
    print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
  }
  if (length(pop_samples) == samplethreshold) {
    print(paste0("Keeping ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
  }
  if (length(pop_samples) > samplethreshold) {
    subsamplepops <- pop_samples[sample(1:length(pop_samples), size=samplethreshold, replace = F)]
    for (b in 1:length(pop_samples)){
      if (pop_samples[b] %in% subsamplepops == FALSE){
        temp[, analysis][pop_samples[b]] <- NA
      }
      print(paste0("Reducing ",x, " to a sample size of ", samplethreshold," samples"))
    }
  } else {}
}




# ## If you want to just remove any sites that are below the threshold (but keep individuals that are above a threshold, then unhash the below section)
# for (x in unique(na.omit(temp[, analysis]))) {
#   pop_samples <- which(temp[, analysis] == x)
#   if (length(pop_samples) < samplethreshold) {
#     temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
#     print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
#   } else {}
# }

m1$analyses<-temp
dm        <- dart.meta.data.merge(d3, m1)
fields    <- c(analysis)
dms     <- data.by.meta.fields(dm, fields, RandRbase, species, dataset, object=analysis); print(dms$meta$analyses[,analysis])


# Generate pairwise fst for populations
# Determine if IBD is significant - Mantel test


#check this function uses a MAF of 0.05 rather than 0.2 (20% is quite high)
pFst <- population.pw.Fst(dms, dms$meta$analyses[, analysis], RandRbase, species, dataset)
pS <- population.pw.spatial.dist(dms, dms$meta$analyses[, analysis])

fst_dir <- paste0(outputloc, "/fst_limiting_", samplethreshold, "_Samps/")
if (!dir.exists(fst_dir)) {
  dir.create(fst_dir)
}

# Plot of fst vs geographic distance
tiff(paste0(fst_dir, species, "_fstplot_limiting_", samplethreshold, "_Samps.tiff"),
     units = "in", width = 13.3, height = 7.5, res = 300)
par(mfrow = c(1, 2))
diag(pS$S) <- NA
diag(pFst$Fst) <- NA

Fst_sig <- cbind(melt(pS$S), unlist(as.list(pFst$Fst)))

colnames(Fst_sig)[3] <- "Geo_dist"
colnames(Fst_sig)[4] <- "Fst"
Fst_sig$Geo_dist2 <- Fst_sig$Geo_dist / 1000

plot(Fst_sig$Geo_dist2, Fst_sig$Fst, xlab = "distance (km)", ylab = "Fst", cex = 1,
     font = 4, cex.main = 1)
title(main = paste0(species, " pairwise fst plots"), adj = 0.001, font.main = 4)
man <- mantel(xdis = pS$S, ydis = pFst$Fst, permutations = 999, na.rm = TRUE)
legend("bottomright", bty = "n", cex = 0.8, text.col = "blue", element_text(face = "bold"),
       legend = paste("Mantel statistic r is ",
                      format(man$statistic, digits = 4),
                      " P =", format(man$signif)))

# Plot of linearised fst vs geographic distance
Fst_sig <- cbind(melt(pS$S), unlist(as.list(pFst$Fst)))
colnames(Fst_sig)[3] <- "Geo_dist"
colnames(Fst_sig)[4] <- "Fst"
Fst_sig$lin_fst <- Fst_sig$Fst / (1 - Fst_sig$Fst)
Fst_sig$Geo_dist2 <- Fst_sig$Geo_dist / 1000
Fst_sig$log10Geo_dist <- log10(Fst_sig$Geo_dist2)
write.table(Fst_sig,
            paste0(fst_dir, species, " fst_geo data.csv"))
plot(Fst_sig$log10Geo_dist, Fst_sig$lin_fst, xlab = "log10(distance)",
     ylab = "Linearised Fst", font = 4, cex.main = 1)
man2 <- mantel(xdis = pS$S, ydis = (pFst$Fst) / (1 - pFst$Fst), permutations = 999, na.rm = TRUE)
legend("bottomright", bty = "n", cex = 0.8, text.col = "blue", element_text(face = "bold"),
       legend = paste("Mantel statistic r is ",
                      format(man2$statistic, digits = 4),
                      " P =", format(man2$signif)))
dev.off()

cat(paste0("pairwise fst and linearised fst drawn!", "\n"))

# Heatmap of just geographic distance or fst
par(mfrow = c(2, 1), oma = c(0, 0, 1, 0))
geo_d <- pS$S
geo_d[upper.tri(geo_d)] <- NA
rownames(geo_d) <- colnames(pS$S)

dimnames <- list(var1 = colnames(pS$S), var2 = colnames(pS$S))
mat <- matrix(geo_d, ncol = length(colnames(geo_d)), nrow = length(colnames(geo_d)), dimnames = dimnames)
df <- as.data.frame(as.table(mat))

p1 <- ggplot(df, aes(var1, var2)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  xlab("") + ylab("") +
  theme(legend.position = "none", axis.text = element_text(size = 5), plot.title = element_text(face = "bold.italic")) +
  ggtitle(paste0(species, " ", analysis, " heatmaps", "\n", "heatmap of pairwise geographic distance limiting ", samplethreshold, " samples per pop"))

genetic_d <- pFst$Fst
genetic_d[upper.tri(genetic_d)] <- NA
rownames(genetic_d) <- colnames(pFst$Fst)

dimnames2 <- list(var1 = colnames(pFst$Fst), var2 = colnames(pFst$Fst))
mat2 <- matrix(genetic_d, ncol = length(colnames(geo_d)), nrow = length(colnames(geo_d)), dimnames = dimnames)
df2 <- as.data.frame(as.table(mat2))
df3 <- df2[complete.cases(df2$Freq),]

p2 <- ggplot(df3, aes(var1, var2)) +
  geom_tile(aes(fill = Freq), colour = "white", na.rm = TRUE) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = round(Freq, 3)), size = 5, df3) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  xlab("") + ylab("") +
  theme(legend.position = "none", axis.text = element_text(size = 8),
        plot.title = element_text(face = "bold.italic")) + ggtitle("heatmap of pairwise Fst")

plot_grid(p1, p2, ncol = 1)

ggsave(paste0(fst_dir, species, "_fstheatmaps_limiting_", samplethreshold, "_Samps.tiff"),
       width = 13.3, height = 7.5, dpi = 300, units = "in", device = 'tiff')

cat(paste0("geographic and pairwise fst heatmap drawn!", "\n"))

####save matrix
genetic_d <-pFst$Fst
write.table(genetic_d, paste0(fst_dir,species," fst matrix.csv"),sep = ",")
geo_d <-pS$S
geo_d[upper.tri(geo_d)] <-geo_d[lower.tri(geo_d)]
new <- matrix(NA, nrow = dim(geo_d)[1], ncol = dim(geo_d)[2])
new[upper.tri(new)] <- geo_d[upper.tri(geo_d)]
new[lower.tri(new)] <- genetic_d[lower.tri(genetic_d)]
colnames(new) <- colnames(geo_d)
rownames(new) <- rownames(geo_d)

write.table(new, paste0(fst_dir,species," heatmap matrix.csv"),sep = ",")


#### FST distance script 

fst_dms <- dms

# calculate FST and geodist
pop_gds_file <- dart2gds(fst_dms, RandRbase, species, dataset)
pop_pFst <- population.pw.Fst(fst_dms, fst_dms$meta$analyses[,analysis], RandRbase, species, dataset, maf_val = 0.05, miss_val = 0.2)
pop_pS <- population.pw.spatial.dist(fst_dms, fst_dms$meta$analyses[,analysis])

# meta_agg <- m2 %>%
#   group_by(pop_large, genetic_group) %>%
#   summarize(lat = mean(lat, na.rm=TRUE),
#             long = mean(long,na.rm=TRUE),
#             .groups = 'drop')



dms_df <- data.frame(cbind(fst_dms$meta$analyses[,analysis],fst_dms$meta$lat, fst_dms$meta$long))
colnames(dms_df) <- c("ansites", "lat", "long")

dms_df$lat <- as.numeric(dms_df$lat)
dms_df$long <- as.numeric(dms_df$long)


meta_agg <- dms_df %>% 
  group_by(dms_df$ansites) %>% 
  dplyr::summarize(lat = mean(lat, na.rm=TRUE),
                   long = mean(long, na.rm=TRUE), 
                   .groups = 'drop')


meta_agg <- na.omit(meta_agg)


colnames(meta_agg) <- c("gsites", "lat", "long")


####plot IBD plot

# Make self comparisons NA
diag(pop_pFst$Fst) <- NA
diag(pop_pS$S) <- NA

# Mantel test
pop_man <- mantel(xdis = pop_pS$S, ydis = pop_pFst$Fst, permutations = 10000, na.rm = TRUE)
pop_man

# mantel plot
pop_Fst_sig <- cbind(melt(pop_pS$S), unlist(as.list(pop_pFst$Fst)))
colnames(pop_Fst_sig)[3] <- "Geo_dist"
colnames(pop_Fst_sig)[4] <- "Fst"
pop_Fst_sig$Geo_dist2 <- pop_Fst_sig$Geo_dist / 1000
pop_Fst_sig <- pop_Fst_sig[!is.na(pop_Fst_sig$Geo_dist),]

# # adding metadata for pop_larges
# pop_Fst_sig2 <- merge(pop_Fst_sig, distinct(meta_agg[, c("pop_large", "genetic_group","lat","long")]), by.x = "Var1", by.y = "pop_large", all.y = FALSE)
# pop_Fst_sig2 <- merge(pop_Fst_sig2, distinct(meta_agg[, c("pop_large", "genetic_group","lat","long")]), by.x = "Var2", by.y = "pop_large", all.y = FALSE)
# pop_Fst_sig2$same_sp <- ifelse(pop_Fst_sig2$genetic_group.x == pop_Fst_sig2$genetic_group.y, "Within group", "Between group")
# pop_Fst_sig2<- distinct(pop_Fst_sig2)  



# adding metadata for pop_larges
pop_Fst_sig2 <- merge(pop_Fst_sig, distinct(meta_agg[,c("gsites","lat","long")]), by.x = "Var1", by.y = "gsites", all.y = FALSE)
pop_Fst_sig2 <- merge(pop_Fst_sig2, distinct(meta_agg[, c("gsites","lat","long")]), by.x = "Var2", by.y = "gsites", all.y = FALSE)
# pop_Fst_sig2$same_sp <- ifelse(pop_Fst_sig2$genetic_group.x == pop_Fst_sig2$genetic_group.y, "Within group", "Between group")
pop_Fst_sig2<- distinct(pop_Fst_sig2)  





#Mapping
row_sub = apply(data.frame(fst_dms$meta$long), 1, function(row) all(row !=0 ))
newdmslong <- data.frame(fst_dms$meta$long)[row_sub,]
row_sub = apply(data.frame(fst_dms$meta$lat), 1, function(row) all(row !=0 ))
newdmslat <- data.frame(fst_dms$meta$lat)[row_sub,]

#this bit looks for range of latlong values to zoom in on and calulcates range for zoom for ggmap
divxlims <- c(min(na.omit(newdmslong))-1,max(na.omit(newdmslong))+1) #default for large distribution s around 0.6
divylims <- c(min(na.omit(newdmslat))-1,max(na.omit(newdmslat))+1) #default for large distribution s around 0.4
divxlimsrange <- abs(c(min(na.omit(newdmslong))-max(na.omit(newdmslong))))
divylimsrange <- abs(c(min(na.omit(newdmslat))-max(na.omit(newdmslat))))

if(divylimsrange > 5) {gmapzoom= 4}
if(divylimsrange > 4 && divylimsrange < 5)  {gmapzoom= 5}
if(divylimsrange > 3 && divylimsrange < 4)  {gmapzoom= 6}
if(divylimsrange > 2 && divylimsrange < 3)  {gmapzoom= 7}
if(divylimsrange > 1 && divylimsrange < 2)  {gmapzoom= 8}
if(divylimsrange > 0.2 && divylimsrange < 1) {gmapzoom= 10}
if(divylimsrange < 0.2) {gmapzoom= 12}


data <- data.frame(sample=dms$sample_names,
                   site=dms$meta$site,
                   analysis=dms$meta$analyses[,analysis],
                   lat=dms$meta$lat,
                   long=dms$meta$long)
data2 <- data[order(data$site),]
data3 <- data2[order(-data2$lat),]
data3$site <- factor(data3$site, levels = unique(data3$site),ordered = TRUE)

#this arranges the sites according to latitude
rep <- ddply(data3,
             .(analysis),
             summarise,
             lat  = mean(lat),
             long = mean(long))
rep2 <- rep[order(-rep$lat),]
rep2$analysis <-factor(rep2$analysis, levels=unique(rep2$analysis))

###colours for each group
rep_an <- ddply(data3,
                .(analysis),
                summarise,
                lat  = mean(lat),
                long = mean(long))
rep_an2 <- rep_an[order(-rep_an$lat),]
rval <- as.numeric(length(unique(data3$analysis)))    ##this counts the number(length) of unique characters for the variable set as analysis at the beginning
rep_an2$bgcols <- c(rev(matlab.like2(rval)))
rep_an3 <- rep_an2[,c("analysis","bgcols")]
rep_all <- merge(rep2,rep_an3, by= "analysis")
rep_all$analysis <-factor(rep_all$analysis, levels=unique(rep_an3$analysis))
rep_all <- rep_all[order(-rep_all$lat),]

for (q in 1:length(rep_all$analysis)){
  rep_all$number[q] <- q
  rep_all$siteAndNumber[q] <- paste0(rep_all$number[q]," ",rep_all$analysis[q]," (", length(which(data3$analysis==rep_all$analysis[q])), ")")
}

#googlemaps base and plot
library(ggmap)
register_google(key = "AIzaSyBZIYV071yUULZTtcUMUZTMwOgRJrKAcaw")

#gmapzoom = 5

if (length(rep_an3$analysis) > 30){legendcolumns = 3}
if (length(rep_an3$analysis) > 20 && length(rep_an3$analysis) < 30){legendcolumns = 2}
if (length(rep_an3$analysis) < 20){legendcolumns = 1}

googlemapbase = get_map(location = c(left = divxlims[1], bottom = divylims[1], right = divxlims[2], top=divylims[2]), zoom = gmapzoom, scale = "auto", maptype = "satellite", source = c("google"))

pop_Fst_sig2 <- pop_Fst_sig2[order(-pop_Fst_sig2$Fst),]


Fstthreshvals <-c(0.1, 0.2, 0.3, 0.4, 0.5)
pop_Fst_sig2$lineWid <- NA

for (m in 1:length(pop_Fst_sig2$Fst)){
  if (pop_Fst_sig2$Fst[m]<Fstthreshvals[1]){pop_Fst_sig2$lineWid[m] <- 0.1}
  if (pop_Fst_sig2$Fst[m]>Fstthreshvals[1] && pop_Fst_sig2$Fst[m]<Fstthreshvals[2]){pop_Fst_sig2$lineWid[m] <- 0.1}
  if (pop_Fst_sig2$Fst[m]>Fstthreshvals[2] && pop_Fst_sig2$Fst[m]<Fstthreshvals[3]){pop_Fst_sig2$lineWid[m] <-  0.8}
  if (pop_Fst_sig2$Fst[m]>Fstthreshvals[3] && pop_Fst_sig2$Fst[m]<Fstthreshvals[4]){pop_Fst_sig2$lineWid[m] <- 0.05}
  if (pop_Fst_sig2$Fst[m]>Fstthreshvals[4] && pop_Fst_sig2$Fst[m]<Fstthreshvals[5]){pop_Fst_sig2$lineWid[m] <- 0.01}
  if (pop_Fst_sig2$Fst[m]>Fstthreshvals[5]){pop_Fst_sig2$lineWid[m] <- 0.005}
}



register_stadiamaps(key = "4fe26ca2-aa51-4e9e-8151-dca93798847e")

stadiamapbase = get_stadiamap(bbox = c(left = divxlims[1], bottom = divylims[1], right = divxlims[2], top=divylims[2]), zoom = 9, maptype = "stamen_terrain", crop=F)


#import GHA shapefiles
library(sf)

BR <- read_sf("C:/Users/dimonr/OneDrive - DPIE/R/Projects/GHA/Border Ranges GHA 10km Buffer.shp")
WP <- read_sf("C:/Users/dimonr/OneDrive - DPIE/R/Projects/GHA/Washpool GHA 10km Buffer.shx")
IL <- read_sf("C:/Users/dimonr/OneDrive - DPIE/R/Projects/GHA/Illuka GHA 10km Buffer.shp")
DR <- read_sf("C:/Users/dimonr/OneDrive - DPIE/R/Projects/GHA/Dorrigo GHA 10km Buffer.shp")
WK <- read_sf("C:/Users/dimonr/OneDrive - DPIE/R/Projects/GHA/Werrikimbe GHA 10km Buffer.shp")
BA <- read_sf("C:/Users/dimonr/OneDrive - DPIE/R/Projects/GHA/Barrington GHA 10km Buffer.shp")
DRWK <- read_sf("C:/Users/dimonr/OneDrive - DPIE/R/Projects/GHA/DorrigoANDWerrikimbe GHA 10km Buffer.shp")

BR <- as_Spatial(BR)
WP <- as_Spatial(WP)
IL <- as_Spatial(IL)
DR <- as_Spatial(DR)
WK <- as_Spatial(WK)
BA <- as_Spatial(BA)
DRWK <- as_Spatial(DRWK)

# google map cropped by lat and long
ggmap(stadiamapbase) + theme_few()+
  coord_cartesian()+ coord_fixed()+
  geom_polygon(data=BR, aes(long, lat, group=group),colour = alpha("darkred", 0.5), linewidth = 1,fill="darkred", alpha = 0.2)+
  geom_polygon(data=WP, aes(long, lat, group=group),colour = alpha("darkorange", 0.5), linewidth = 1, fill="darkorange", alpha = 0.2)+
  # geom_polygon(data=IL, aes(long, lat, group=group),colour = alpha("pink", 0.5), linewidth = 1, alpha = 0)+
  # geom_polygon(data=DRWK, aes(long, lat, group=group),colour = alpha("cyan", 0.5), linewidth = 1, fill="cyan", alpha = 0.2)+
  geom_polygon(data=DR, aes(long, lat, group=group),colour = alpha("cyan", 0.5), linewidth = 1, fill="cyan", alpha = 0.2)+
  geom_polygon(data=WK, aes(long, lat, group=group),colour = alpha("green", 0.5), linewidth = 1, fill="green", alpha = 0.2)+
  geom_polygon(data=BA, aes(long, lat, group=group),colour = alpha("yellow", 0.5), linewidth = 1, fill="yellow", alpha = 0.2)+
  geom_segment(data=pop_Fst_sig2[pop_Fst_sig2$Fst<=1,], linewidth=3, alpha=1,na.rm = TRUE,
               # geom_segment(data=pop_Fst_sig2,
               aes(x=long.x, y=lat.x,
                   xend = long.y, yend = lat.y,
                   colour=Fst), lineend = "round")+
  # geom_segment(data=pop_Fst_sig2[pop_Fst_sig2$Fst<=0.1,], linewidth=1, alpha=0.5,na.rm = TRUE,
  #              # geom_segment(data=pop_Fst_sig2,
  #              aes(x=long.x, y=lat.x,
  #                  xend = long.y, yend = lat.y,
  #                  colour=Fst), lineend = "round")+
  # scale_linewidth_manual(name="Fst", values = as.character(pop_Fst_sig2$Fst)) +
  # scale_linewidth_manual(pop_Fst_sig2$Fst, limits=factor(0.0001, 0.8))+
  scale_color_gradient(low = "red", high = "white")+
  # geom_text(data=pop_Fst_sig2[pop_Fst_sig2$Fst<=0.5,],
  #           aes(label = round(Fst, 2), x = (long.x + long.y) / 2, y = (lat.y+lat.x) / 2),
  #           size = 2,
  #           color = "black")+
  geom_point(data=meta_agg, mapping=aes(x=long, y=lat), col="black", size=1)+
  # xlim(lims[1], lims[2])+
  # ylim(lims[3], lims[4])+   
  labs(x = "Longitude", y = "Latitude", colour="FST") +guides(size = "none", alpha="none")+
  xlim(150.8, 153.9)+ # increase to change width of plot
  ylim(-32.5, -27.7)+
  annotation_north_arrow(location = "br", which_north = "true", pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"), style = north_arrow_fancy_orienteering)+
  ggtitle(paste0(species," ",analysis, " Fst Direction"))
# ggrepel::geom_label_repel(data = meta_agg,aes(x = long, y = lat, label=gsites),
#                           min.segment.length=0.25, color="black",fill="white",size=3, segment.colour="white", alpha=0.9, label.size=0, nudge_y=0.004)

ggsave(paste0(species,"_", analysis, "Fst_Direction_Map.tiff"), path = paste0(gst_dir), width = 10, height = 12, dpi = 300, units = "in")



















##### DIVERSITY Stats limiting X samples per pop #####

m1 <- read.meta.data(d3, RandRbase, species, dataset, fields = as.numeric(ncol(metafile)-4))
temp<-as.data.frame(m1$analyses)

#set a sample threshold for each pop:
samplethreshold <-5# chnage this to what number of samples you want toi have uniform across all pops



#### Below lines remove sites with less than threshold number of samples, and subsample sites that are above the threshold
for (x in unique(na.omit(temp[, analysis]))) {
  pop_samples <- which(temp[, analysis] == x)
  if (length(pop_samples) < samplethreshold) {
    temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
    print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
  }
  if (length(pop_samples) == samplethreshold) {
    print(paste0("Keeping ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
  }
  if (length(pop_samples) > samplethreshold) {
    subsamplepops <- pop_samples[sample(1:length(pop_samples), size=samplethreshold, replace = F)]
    for (b in 1:length(pop_samples)){
      if (pop_samples[b] %in% subsamplepops == FALSE){
        temp[, analysis][pop_samples[b]] <- NA
      }
      print(paste0("Reducing ",x, " to a sample size of ", samplethreshold," samples"))
    }
  } else {}
}




# ## If you want to just remove any sites that are below the threshold (but keep individuals that are above a threshold, then unhash the below section)
# for (x in unique(na.omit(temp[, analysis]))) {
#   pop_samples <- which(temp[, analysis] == x)
#   if (length(pop_samples) < samplethreshold) {
#     temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
#     print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
#   } else {}
# }

m1$analyses<-temp
dm        <- dart.meta.data.merge(d3, m1)
fields    <- c(analysis)
dms     <- data.by.meta.fields(dm, fields, RandRbase, species, dataset, object=analysis); print(dms$meta$analyses[,analysis])

#correct inconsistencies in J's scripts for naming treatment which is for naming files
treatment <- paste0("raw_SNPFilt_1SNPperClone_Field_", analysis)
dms$treatment <- treatment
gta <- dms$gt

gp   <- dart2genepop(dms, RandRbase, species, dataset, dms$meta$analyses[,analysis])

#calculate basic pop gen stats
bs <- basicStats(infile = gp, outfile = NULL,
                 fis_ci = FALSE, ar_ci = TRUE,
                 ar_boots = 999,
                 rarefaction = FALSE, ar_alpha = 0.05)

#create map plot and table

divrda_dir <- paste0(RandRbase,species,"/popgen/",treatment,"/genepop/bs.rda")
save(bs, file = divrda_dir)

if(!dir.exists(outputloc)){
  dir.create(outputloc)
}

tempdir <- paste0(outputloc, "/temp/")
if(!dir.exists(tempdir)){
  dir.create(tempdir)
}

Div_dir <- paste0(outputloc,"/diversity limiting ", samplethreshold, " samples/")
if(!dir.exists(Div_dir)){
  dir.create(Div_dir)
}

{
  dmsanalysis <- dms$meta$analyses[,analysis]
  meta <- data.frame(table(dmsanalysis))
  z=meta$Freq<2
  
  if (sum(z, na.rm = TRUE)>=1) {#if there are any pops with less than 2 indiv...
    meta_lesthan2 <-as.character(meta[,1][meta$Freq<2]) #the pops with less than 2 indiv are...
    myvars <- dmsanalysis%in%meta_lesthan2
    newdmsanalysis <- dmsanalysis[!myvars]
    excluded_sample_names <- dms$sample_names[myvars]
    dms_d <- exclude.samples(dms, excluded_sample_names, remove_fixed_loci=TRUE)
    cat("metadata contains pops less than 2 individual","\n")
    print(meta_lesthan2)
  } else {
    meta_lesthan2 <-as.character(meta[,1][meta$Freq<2]) #the pops with less than 2 indiv are...
    myvars <- dmsanalysis%in%meta_lesthan2
    newdmsanalysis <- dms$meta$analyses[,analysis]
    dms_d <-dms
    excluded_sample_names <- NULL
    cat("metadata is alright without removal of pops","\n")
  }
  
  #making dataframe for  metadata without pops with less than 2 indiv
  dms_meta <- cbind.data.frame(dms$meta$sample_names, dms$meta$analyses[,analysis],dms$meta$lat, dms$meta$long)
  colnames(dms_meta) <- c("sample_names","site", "lat", "long")
  #pp <- dms_meta[!myvars,]
  #newdms_meta <- pp[apply(table(pp$site)>0,1,any),]
  newdms_meta <- dms_meta[!myvars,]
  
  colnames(newdms_meta) <- c("sample_names","site", "lat", "long")
  npop <-length(unique(newdms_meta$site))
  
  #making dataframe of averages (lat/longs) for metadata
  
  newmetalatlongs <- ddply(newdms_meta,
                           c("site"),
                           summarise,
                           lat  = mean(lat),
                           long  = mean(long),
                           N=length(site))
}

npop <- length(unique(newdmsanalysis))
result <- mat.or.vec(npop, 11)
measurement_names <- rownames(bs$main_tab[[1]])
population_names  <-
  names(bs$main_tab) #ls() rearranges the names
rownames(result) <- population_names
colnames(result) <- measurement_names

for (r in 1:npop) {
  popstats <- bs$main_tab[[r]][, "overall"] ##extract from a list
  result[r, ] <- popstats
}

result <- as.data.frame(result)
result$sample_names <- rownames(result)

###getting latlongs into the diversity results
metnew2 <- merge(newdms_meta[, 1:2], newmetalatlongs[, 1:4], by = "site")

data2 <- merge(result, metnew2, by = "sample_names", all.x = TRUE) #merge data
colnames(data2)[which(names(data2) == "dms.meta.analyses...analysis.")] <-
  "dms.meta.site"
data2$species <- as.character(species)

gp   <- dart2genepop0(dms, RandRbase, species, dataset, newdmsanalysis, maf_val = 0)

gp_genind <- read.genepop(gp, ncode = 2)

newdms_meta$site <- factor(newdms_meta$site, levels = rev(unique(newdms_meta$site)))
gp_genind@other <- newdms_meta
strata(gp_genind) <- gp_genind@other
setPop(gp_genind) <- ~ site

p_allele <-
  data.frame(rowSums(private_alleles(gp_genind, locus ~ site,
                                     count.alleles = F)))
p_allele$site <- rownames(p_allele)
rownames(p_allele) <- NULL
names(p_allele)[names(p_allele) == "rowSums.private_alleles.gp_genind..locus...site..count.alleles...F.."] <-
  "n_pa"
data <- merge(data2, p_allele, by.x = "site", by.y = "site")

write.table(data,
            paste0(Div_dir, "/", species, " diveRsity stats.csv"),
            sep = ",",
            row.names = F)

###plot

cexsize <- 0.5
ptsize <-0.5

tiff(paste0(outputloc,"/temp/",species," diveRsity stats2.tiff"), units="in", width=11.7, height=4, res=300)
par(mfrow=c(1,3), ## only works if plot out and not export tiff
    mar=c(1,1,1,1),
    oma=c(0,0,0.2,0)) #margins between nsw plots

row_sub_long = apply(data.frame(newdms_meta$long), 1, function(row) all(row !=0 ))
newdmslong <- data.frame(newdms_meta$long)[row_sub_long,]
row_sub_lat = apply(data.frame(newdms_meta$lat), 1, function(row) all(row !=0 ))
newdmslat <- data.frame(newdms_meta$lat)[row_sub_lat,]

if (file.exists(paste0(outputloc,"/temp/site_number.csv"))) {
  rep2 <- read.csv(paste0(outputloc,"/temp/site_number.csv"),sep=" ")
  
  testResult <- unlist(lapply(unique(dms$meta$analyses[,analysis]), function(s) {all(s %in% rep2$site)}))
  if (!all(testResult)) {
    dms_meta <- cbind.data.frame(dms$meta$sample_names, dms$meta$analyses[,analysis],dms$meta$lat, dms$meta$long)
    colnames(dms_meta) <- c("sample_names","site", "lat", "long")
    dms_meta$lat <- as.numeric(dms_meta$lat)
    dms_meta$long <- as.numeric(dms_meta$long)
    ###
    dms_meta2 <- dms_meta[order(dms_meta$site), ]
    dms_meta2 <- dms_meta2[order(-dms_meta2$lat), ]
    ####make sure your rep has nums_n_site which is in df2!!!
    dms_meta2$site <- factor(dms_meta2$site, levels = unique(dms_meta2$site),ordered = TRUE)
    
    rep <- ddply(dms_meta2,
                 .(site),
                 summarise,
                 lat  = mean(lat),
                 long = mean(long))
    rep2 <- rep[order(-rep$lat), ]
    x <- length(unique(dms_meta2$site))
    bgcols = rev(matlab.like2(x))
    rep2$nums_n_site <- paste0(1:x, ": ", rep2$site)
    rep2$num_site <- 1:x
    write.table(rep2, paste0(outputloc,"/temp/site_number_fromdiv.csv"))
  }
  
} else {
  dms_meta <- cbind.data.frame(dms$meta$sample_names, dms$meta$analyses[,analysis],dms$meta$lat, dms$meta$long)
  colnames(dms_meta) <- c("sample_names","site", "lat", "long")
  dms_meta$lat <- as.numeric(dms_meta$lat)
  dms_meta$long <- as.numeric(dms_meta$long)
  ###
  dms_meta2 <- dms_meta[order(dms_meta$site), ]
  dms_meta2 <- dms_meta2[order(-dms_meta2$lat), ]
  ####make sure your rep has nums_n_site which is in df2!!!
  dms_meta2$site <- factor(dms_meta2$site, levels = unique(dms_meta2$site),ordered = TRUE)
  rep <- ddply(dms_meta2,
               .(site),
               summarise,
               lat  = mean(lat),
               long = mean(long))
  rep2 <- rep[order(-rep$lat), ]
  x <- length(unique(dms_meta2$site))
  bgcols = rev(matlab.like2(x))
  rep2$nums_n_site <- paste0(1:x, ": ", rep2$site)
  rep2$num_site <- 1:x
  write.table(rep2, paste0(outputloc,"/temp/site_number_fromdiv.csv"))
}


data11 <- merge(data,rep2[,c("site","nums_n_site","num_site")],by="site")
data11$lab <- data11$num_site

divxlims <- c(min(na.omit(rep2$long))-0.4,max(na.omit(rep2$long))+0.4)
divylims <- c(min(na.omit(rep2$lat))-0.4,max(na.omit(rep2$lat))+0.4)

if ((divxlims[2]-divxlims[1]) <3) {
  rad <- 0.2
} else {
  rad <- 0.5}

## allelic richness
oz(xlim=divxlims, ylim=divylims, lwd=0.3)#divxlims taken from PCA section
title(expression(italic("ar")),font=2,cex.main=2)
points(x=data11$long, y=data11$lat,pch=16,bg="black",cex=ptsize, lwd=0.3)
draw.bubble(data11$long, data11$lat, (data11$ar)^50, bg=alpha("red",0.2), pch=21, maxradius = rad, lwd=0.3)
data11 <- data11[complete.cases(data11), ]#remove latlongs that are NAs
pointLabel(data11$long, data11$lat, labels = paste("  ", data11$lab, "  ", sep=""), cex=0.5,offset=0.25)
box()

## expected heterozygosity
oz(xlim=divxlims, ylim=divylims, lwd=0.3)
points(x=data11$long, y=data11$lat,pch=16,bg="black",cex=ptsize, lwd=0.3)
draw.bubble(data11$long, data11$lat, (data11$exp_het)^15, bg=alpha("red",0.2), pch=21, maxradius = rad, lwd=0.3)
title(expression(italic("exp_het")),font=4,cex.main=2)
pointLabel(data11$long, data11$lat, labels = paste("  ", data11$lab, "  ", sep=""), cex=0.5,offset=0.25)
box()

## observed heterozygosity
oz(xlim=divxlims, ylim=divylims, lwd=0.3)
points(x=data11$long, y=data11$lat,pch=16,bg="black",cex=ptsize, lwd=0.3)
pointLabel(data11$long, data11$lat, labels = paste("  ", data11$lab, "  ", sep=""), cex=0.5,offset=0.25)

if(min(data11$obs_het)-max(data11$obs_het)==0) {
  title(expression(italic("obs_het = 0")),font=4,cex.main=0.7)
} else {
  draw.bubble(data11$long, data11$lat, (data11$obs_het)^10, bg=alpha("red",0.2),pch=21, maxradius = rad, lwd=0.3)
  title(expression(italic("obs_het")),cex.main=2)
}
box()
dev.off()


#this bit checks for zeros in lats and longs and cut out those samples
row_sub = apply(data.frame(dms$meta$long), 1, function(row) all(row !=0 ))
newdmslong <- data.frame(dms$meta$long)[row_sub,]
row_sub = apply(data.frame(dms$meta$lat), 1, function(row) all(row !=0 ))
newdmslat <- data.frame(dms$meta$lat)[row_sub,]

#this bit looks for range of latlong values to zoom in on and calulcates range for zoom for ggmap
divxlims <- c(min(na.omit(newdmslong))-1,max(na.omit(newdmslong))+1)
divylims <- c(min(na.omit(newdmslat))-1,max(na.omit(newdmslat))+1)
divxlimsrange <- abs(c(min(na.omit(newdmslong))-max(na.omit(newdmslong))))
divylimsrange <- abs(c(min(na.omit(newdmslat))-max(na.omit(newdmslat))))


coordinate_cities <- read.csv(paste0(maindir,"CityMapCordinates.csv"))
sf_oz <- ozmap_data("states")

#order and map expected het as a coloured range
data11 <- data11[order(data11$exp_het),]

ggplot(sf_oz, aes(fill="white")) + geom_sf() +
  scale_fill_manual(values = "white", guide = "none") +
  theme_void() + geom_sf() +
  xlim(divxlims) + ylim(divylims)+
  geom_point(data = data11, mapping = aes(x = long, y = lat, color = exp_het), size=5, alpha=0.9) +
  scale_color_gradient(low="yellow", high="red", name="Exp Het") +
  geom_point(data=coordinate_cities, mapping = aes(x = long, y = lat), size=1)+
  geom_text(data=coordinate_cities, mapping = aes(x = long, y = lat, label=City),hjust=0,nudge_y=0.03, size=3)


ggsave(paste0(Div_dir, species, "_Expected_Het_Heatmap.tiff"), width = 10, height = 12, dpi = 300, units = "in")




#order and map observed het as a coloured range
data11 <- data11[order(data11$obs_het),]


ggplot(sf_oz, aes(fill="white")) + geom_sf() +
  scale_fill_manual(values = "white", guide = "none") +
  theme_void() + geom_sf() +
  xlim(divxlims) + ylim(divylims)+
  geom_point(data = data11, mapping = aes(x = long, y = lat, color = obs_het), size=5, alpha=0.9) +
  scale_color_gradient(low="yellow", high="red", name="Obs Het") +
  geom_point(data=coordinate_cities, mapping = aes(x = long, y = lat), size=1)+
  geom_text(data=coordinate_cities, mapping = aes(x = long, y = lat, label=City),hjust=0,nudge_y=0.03, size=3)


ggsave(paste0(Div_dir, species, "_Observed_Het_Heatmap.tiff"), width = 10, height = 12, dpi = 300, units = "in")




#order and map allelic richness as a coloured range
data11 <- data11[order(data11$ar),]


ggplot(sf_oz, aes(fill="white")) + geom_sf() +
  scale_fill_manual(values = "white", guide = "none") +
  theme_void() + geom_sf() +
  xlim(divxlims) + ylim(divylims)+
  geom_point(data = data11, mapping = aes(x = long, y = lat, color = ar), size=5, alpha=0.9) +
  scale_color_gradient(low="yellow", high="red", name="ar") +
  geom_point(data=coordinate_cities, mapping = aes(x = long, y = lat), size=1)+
  geom_text(data=coordinate_cities, mapping = aes(x = long, y = lat, label=City),hjust=0,nudge_y=0.03, size=3)


ggsave(paste0(Div_dir, species, "_Allelic_Richness_Heatmap.tiff"), width = 10, height = 12, dpi = 300, units = "in")







cat(paste0("all diversity maps are drawn!","\n"))
################combine div maps with a table
setwd(outputloc)
###table
geodist_data <- get_geodistinfo(newdms_meta)

data11 <- merge(data,rep2[,c("site","nums_n_site","num_site")],by="site")
data11$lab <- data11$num_site

if ("n_pa" %in% names(data11)){
}else{
  data11$n_pa <- data$rowSums.private_alleles.gp_genind..locus...site..count.alleles...F..
}

data3 <- merge(data11, geodist_data, by="site", all.x=TRUE)
data2grob <- cbind.data.frame(data3$num_site,data3$site,data3$lat,data3$long, data3$ar, data3$exp_het,data3$obs_het,data3$n_pa, data3$fis, data3$Freq, data3$meanDist)
colnames(data2grob) <- c("","site","lat","long","allelic_richness (ar)"," expected_het "," observed_het ","n_private_alleles","      fis      "," n_samples ","  meanDist  ")
data2grob <- data2grob[order(-data2grob$lat),]

ftsz <- 0.2
tiff(paste0("temp/",species," diveRsity grobtable.tiff"),units="in", width=11.7, height=5, res=300)
tt <- ttheme_default(base_size =7,
                     core=list(fg_params=list(fontface=3)),
                     colhead=list(fg_params=list(col="navyblue", fontface=2L)),
                     rowhead=list(fg_params=list(col="black", fontface=2L)))
tg <- tableGrob(data2grob,theme=tt,rows=NULL)
tg$heights <- rep(unit(0.2,"null"), nrow(tg))
grid.draw(tg)
dev.off()

cat(paste0("diversity table drawn!","\n"))

###now combine!

loc_diversity<- dir(path = "temp", pattern = "stats2.tiff$", full.names = TRUE, recursive = TRUE) ###must set working directory to make readPNG work
loc_diversitytab <- dir(path = "temp", pattern = "grobtable.tiff$", full.names = TRUE, recursive = TRUE) ###must set working directory to make readPNG work

lloo <- c(loc_diversity,loc_diversitytab)

plots <- lapply(ll <- lloo,function(x){
  img <- as.raster(readTIFF(x)) ##change this to readpng if readPNG doesnt work
  rasterGrob(img, interpolate = FALSE)
})

# lay <- rbind(c(1),c(1), c(2),c(2),c(2))

lay <- rbind(c(1),c(2))

ggsave(paste0(Div_dir,species, "_divntab_limting5Samps.tiff"),
       width=11.7, height=8.3,dpi = 300, units = "in", device='tiff',
       marrangeGrob(grobs = plots, layout_matrix = lay,top=textGrob(paste0(analysis, " Diversity Outputs Limting ", samplethreshold, " Individuals per Pop (equal numbers)"), gp=gpar(fontsize=16,fontface=4))))

file.remove(lloo)

cat(paste0("diversity maps and table combined!","\n"))

setwd(maindir)



##### Splitstree #####
#create an input file to then open and run on Splitstree


library(tanggle)EA
library(phangorn)
library(ggtree)

splitsdir <- paste0(outputloc,"/Splitstree/")

if(!dir.exists(splitsdir)){
  dir.create(splitsdir)
}

m1 <- read.meta.data(d3, RandRbase, species, dataset, fields = as.numeric(ncol(metafile)-4))
temp<-as.data.frame(m1$analyses)

#set a sample threshold for each pop:
samplethreshold <-1# chnage this to what number of samples you want toi have uniform across all pops

## If you want to just remove any sites that are below the threshold (but keep individuals that are above a threshold, then unhash the below section)
for (x in unique(na.omit(temp[, analysis]))) {
  pop_samples <- which(temp[, analysis] == x)
  if (length(pop_samples) < samplethreshold) {
    temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
    print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
  } else {}
}
m1$analyses<-temp
dm        <- dart.meta.data.merge(d3, m1)
fields    <- c(analysis)
dms     <- data.by.meta.fields(dm, fields, RandRbase, species, dataset, object=analysis); print(dms$meta$analyses[,analysis])

#Note, this modified version of dart2splitstree needs outputloc/splitsreee/ dorectory to be created
source("C:/Users/dimonr/OneDrive - DPIE/R/Packages/RRtools/R/dart2splitstree.R")
snp <- dart2splitstree(dms, RandRbase, species, dataset,  dms$meta$analyses[,analysis], add_pop=TRUE, output=splitsdir)
#output nexus file will be in main directory after RandRbase
# paste(splitsdir,species,"_",analysis,"_SPLITSTREE.nex",sep="")




##### Mapping Script #####
m1 <- read.meta.data(d3, RandRbase, species, dataset, fields = as.numeric(ncol(metafile)-4))
temp<-as.data.frame(m1$analyses)

#set a sample threshold for each pop:
samplethreshold <-1# chnage this to what number of samples you want toi have uniform across all pops

## If you want to just remove any sites that are below the threshold (but keep individuals that are above a threshold, then unhash the below section)
for (x in unique(na.omit(temp[, analysis]))) {
  pop_samples <- which(temp[, analysis] == x)
  if (length(pop_samples) < samplethreshold) {
    temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
    print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
  } else {}
}
m1$analyses<-temp
dm        <- dart.meta.data.merge(d3, m1)
fields    <- c(analysis)
dms     <- data.by.meta.fields(dm, fields, RandRbase, species, dataset, object=analysis); print(dms$meta$analyses[,analysis])


if(!dir.exists(paste0(outputloc, "/Map"))){
  dir.create(paste0(outputloc, "/Map"))
}

#this bit checks for zeros in lats and longs and cut out those samples
row_sub = apply(data.frame(dms$meta$long), 1, function(row) all(row !=0 ))
newdmslong <- data.frame(dms$meta$long)[row_sub,]
row_sub = apply(data.frame(dms$meta$lat), 1, function(row) all(row !=0 ))
newdmslat <- data.frame(dms$meta$lat)[row_sub,]

#this bit looks for range of latlong values to zoom in on and calulcates range for zoom for ggmap
divxlims <- c(min(na.omit(newdmslong))-1,max(na.omit(newdmslong))+1) #default for large distribution s around 0.6
divylims <- c(min(na.omit(newdmslat))-1,max(na.omit(newdmslat))+1) #default for large distribution s around 0.4
divxlimsrange <- abs(c(min(na.omit(newdmslong))-max(na.omit(newdmslong))))
divylimsrange <- abs(c(min(na.omit(newdmslat))-max(na.omit(newdmslat))))

if(divylimsrange > 5) {gmapzoom= 4}
if(divylimsrange > 4 && divylimsrange < 5)  {gmapzoom= 5}
if(divylimsrange > 3 && divylimsrange < 4)  {gmapzoom= 6}
if(divylimsrange > 2 && divylimsrange < 3)  {gmapzoom= 7}
if(divylimsrange > 1 && divylimsrange < 2)  {gmapzoom= 8}
if(divylimsrange > 0.2 && divylimsrange < 1) {gmapzoom= 10}
if(divylimsrange < 0.2) {gmapzoom= 12}


data <- data.frame(sample=dms$sample_names,
                   site=dms$meta$site,
                   analysis=dms$meta$analyses[,analysis],
                   lat=dms$meta$lat,
                   long=dms$meta$long)
data2 <- data[order(data$site),]
data3 <- data2[order(-data2$lat),]
data3$site <- factor(data3$site, levels = unique(data3$site),ordered = TRUE)

#this arranges the sites according to latitude
rep <- ddply(data3,
             .(analysis),
             summarise,
             lat  = mean(lat),
             long = mean(long))
rep2 <- rep[order(-rep$lat),]
rep2$analysis <-factor(rep2$analysis, levels=unique(rep2$analysis))

###colours for each group
rep_an <- ddply(data3,
                .(analysis),
                summarise,
                lat  = mean(lat),
                long = mean(long))
rep_an2 <- rep_an[order(-rep_an$lat),]
rval <- as.numeric(length(unique(data3$analysis)))    ##this counts the number(length) of unique characters for the variable set as analysis at the beginning
rep_an2$bgcols <- c(rev(matlab.like2(rval)))
rep_an3 <- rep_an2[,c("analysis","bgcols")]
rep_all <- merge(rep2,rep_an3, by= "analysis")
rep_all$analysis <-factor(rep_all$analysis, levels=unique(rep_an3$analysis))
rep_all <- rep_all[order(-rep_all$lat),]

for (q in 1:length(rep_all$analysis)){
  rep_all$number[q] <- q
  rep_all$siteAndNumber[q] <- paste0(rep_all$number[q]," ",rep_all$analysis[q]," (", length(which(data3$analysis==rep_all$analysis[q])), ")")
}

#googlemaps base and plot
install_github("stadiamaps/ggmap")
library(ggmap)
register_google(key = "AIzaSyBZIYV071yUULZTtcUMUZTMwOgRJrKAcaw")
register_stadiamaps(key = "4fe26ca2-aa51-4e9e-8151-dca93798847e")


#gmapzoom = 5

if (length(rep_an3$analysis) > 30){legendcolumns = 3}
if (length(rep_an3$analysis) > 20 && length(rep_an3$analysis) < 30){legendcolumns = 2}
if (length(rep_an3$analysis) < 20){legendcolumns = 1}



googlemapbase = get_map(location = c(left = divxlims[1], bottom = divylims[1], right = divxlims[2], top=divylims[2]), zoom = gmapzoom, scale = "auto", maptype = "satellite", source = c("google"))

# google map cropped by lat and long
ggmap(googlemapbase) +
  geom_point(data = rep_all, mapping = aes(x = long, y = lat, col=analysis, shape=analysis), size=6, , stroke=1.5)+
  geom_text_repel(data=rep_all, mapping=aes(x = long, y = lat, label=number, col=analysis, fontface = "bold"), size=3, hjust=0, nudge_x = 0.01) +
  # geom_point(data=coordinate_cities, mapping = aes(x = long, y = lat), size=1) +
  # geom_text(data=coordinate_cities, mapping = aes(x = long, y = lat, label=City),hjust=0,nudge_y=0.03,fontface = "bold", size=5)+
  guides(col = guide_legend(ncol = legendcolumns)) +
  theme(legend.text=element_text(size=12)) +
  ylab(paste0("Lat"))+
  xlab(paste0("Long"))+
  scale_color_manual(labels = rep_all$siteAndNumber, values=rep_an3$bgcols, name = "Sites (samples)")+
  scale_shape_manual(labels = rep_all$siteAndNumber, values=shapeassign, name = "Sites (samples)") +
  ggtitle(paste0("Sampling Distribution and Number of Samples - Analysis: ", analysis))


ggsave(paste0(species,"_Satellite_Distribution_Map.tiff"), path = paste0(outputloc, "/Map/"), width = 13.3, height = 8, dpi = 300, units = "in")




googlemapbase = get_map(location = c(left = divxlims[1], bottom = divylims[1], right = divxlims[2], top=divylims[2]), zoom = gmapzoom, scale = "auto", maptype = "terrain", source = c("google"))

# google map cropped by lat and long
ggmap(googlemapbase) +
  geom_point(data = rep_all, mapping = aes(x = long, y = lat, col=analysis, shape=analysis), size=6, , stroke=1.5)+
  geom_text_repel(data=rep_all, mapping=aes(x = long, y = lat, label=number, col=analysis, fontface = "bold"), size=3, hjust=0, nudge_x = 0.01) +
  # geom_point(data=coordinate_cities, mapping = aes(x = long, y = lat), size=1) +
  # geom_text(data=coordinate_cities, mapping = aes(x = long, y = lat, label=City),hjust=0,nudge_y=0.03,fontface = "bold", size=5)+
  guides(col = guide_legend(ncol = legendcolumns)) +
  theme(legend.text=element_text(size=12)) +
  ylab(paste0("Lat"))+
  xlab(paste0("Long"))+
  scale_color_manual(labels = rep_all$siteAndNumber, values=rep_an3$bgcols, name = "Sites (samples)")+
  scale_shape_manual(labels = rep_all$siteAndNumber, values=shapeassign, name = "Sites (samples)") +
  ggtitle(paste0("Sampling Distribution and Number of Samples - Analysis: ", analysis))


ggsave(paste0(species,"_terrain_Distribution_Map.tiff"), path = paste0(outputloc, "/Map/"), width = 13.3, height = 8, dpi = 300, units = "in")





stadiamapbase = get_stadiamap(bbox = c(left = divxlims[1], bottom = divylims[1], right = divxlims[2], top=divylims[2]), zoom = gmapzoom, maptype = "stamen_terrain", crop=T)

# googlemapbase = get_map(location = c(left = divxlims[1], bottom = divylims[1], right = divxlims[2], top=divylims[2]), zoom = gmapzoom, scale = "auto", maptype = "terrain", source = c("stamen"))

# google map cropped by lat and long
# ggmap(googlemapbase) +
  ggmap(stadiamapbase) +
    geom_point(data = rep_all, mapping = aes(x = long, y = lat, col=analysis, shape=analysis), size=6, , stroke=1.5)+
  geom_text_repel(data=rep_all, mapping=aes(x = long, y = lat, label=number, col=analysis, fontface = "bold"), size=3, hjust=0, nudge_x = 0.01) +
  # geom_point(data=coordinate_cities, mapping = aes(x = long, y = lat), size=1) +
  # geom_text(data=coordinate_cities, mapping = aes(x = long, y = lat, label=City),hjust=0,nudge_y=0.03,fontface = "bold", size=5)+
  guides(col = guide_legend(ncol = legendcolumns)) +
  theme(legend.text=element_text(size=12)) +
  ylab(paste0("Lat"))+
  xlab(paste0("Long"))+
  scale_color_manual(labels = rep_all$siteAndNumber, values=rep_an3$bgcols, name = "Sites (samples)")+
  scale_shape_manual(labels = rep_all$siteAndNumber, values=shapeassign, name = "Sites (samples)") +
  ggtitle(paste0("Sampling Distribution and Number of Samples - Analysis: ", analysis))


ggsave(paste0(species,"_stamen_Distribution_Map.tiff"), path = paste0(outputloc, "/Map/"), width = 13.3, height = 8, dpi = 300, units = "in")





#ozmaps simfplified white background map
coordinate_cities <- read.csv(paste0(maindir,"DetailedCityMapCordinates.csv"))
sf_oz <- ozmap_data("states")
divxlims <- c(min(na.omit(newdmslong))-1,max(na.omit(newdmslong))+1) #default for large distribution s around 0.6
divylims <- c(min(na.omit(newdmslat))+1,max(na.omit(newdmslat))-1) #default for large distribution s around 0.4

ggplot(sf_oz, aes(fill="white")) + geom_sf() +
  scale_fill_manual(values = "white", guide = "none") +
  theme_void() + geom_sf() +
  xlim(divxlims) + ylim(divylims)+
  geom_point(data = rep_all, mapping = aes(x = long, y = lat, col=analysis, shape=analysis), size=6, stroke=1.5)+
  geom_text_repel(data=rep_all, mapping=aes(x = long, y = lat, label=number, col=analysis, fontface = "bold"), size=3, hjust=0, nudge_x = 0.01) +
  geom_point(data=coordinate_cities, mapping = aes(x = long, y = lat), size=1) +
  geom_text(data=coordinate_cities, mapping = aes(x = long, y = lat, label=City),hjust=0,nudge_y=0.01,fontface = "bold", size=5)+
  guides(col = guide_legend(ncol = legendcolumns)) +
  theme(legend.text=element_text(size=12)) +
  ylab(paste0("Lat"))+
  xlab(paste0("Long"))+
  scale_color_manual(labels = rep_all$siteAndNumber, values=rep_an3$bgcols, name = "Sites (samples)")+
  scale_shape_manual(labels = rep_all$siteAndNumber, values=shapeassign, name = "Sites (samples)") +
  ggtitle(paste0("Sampling Distribution and Number of Samples - Analysis: ", analysis))

ggsave(paste0(species,"ozmaps_Distribution_Map.tiff"), path = paste0(outputloc, "/Map/"), width = 13.3, height = 8, dpi = 300, units = "in")





##### LEA - Entropy, Barplot and Piechart with atleast 2 samples per pop #####

#read metadata file
metafile <- read.xlsx(paste0(maindir, species,"/meta/",species,"_",dataset,"_meta.xlsx"))
m1 <- read.meta.data(d3, RandRbase, species, dataset, fields = as.numeric(ncol(metafile)-4))
colnames(m1$analyses)
# create a temp meta dataframe to remove sites/samples that include less than X individuals
temp <- as.data.frame(m1$analyses)

samplethreshold <- 2 # need a minimum of 2 samples per pop to run LEA


##### Below lines remove sites with less than threshold number of samples, and subsample sites that are above the threshold
# for (x in unique(na.omit(temp[, analysis]))) {
#   pop_samples <- which(temp[, analysis] == x)
#   if (length(pop_samples) < samplethreshold) {
#     temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
#     print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
#   }
#   if (length(pop_samples) == samplethreshold) {
#     print(paste0("Keeping ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
#   }
#   if (length(pop_samples) > samplethreshold) {
#     subsamplepops <- pop_samples[sample(1:length(pop_samples), size=samplethreshold, replace = F)]
#     for (b in 1:length(pop_samples)){
#       if (pop_samples[b] %in% subsamplepops == FALSE){
#         temp[, analysis][pop_samples[b]] <- NA
#       }
#       print(paste0("Reducing ",x, " to a sample size of ", samplethreshold," samples"))
#     }
#   } else {}
# }

##### If you want to just remove any sites that are below the threshold (but keep individuals that are above a threshold, then unhash the below section)
for (x in unique(na.omit(temp[, analysis]))) {
  pop_samples <- which(temp[, analysis] == x)
  if (length(pop_samples) < samplethreshold) {
    temp[, analysis] <- gsub(x, replacement = NA, temp[, analysis])
    print(paste0("Removing ", x, " with ", " with a sample size of ", length(pop_samples), " samples"))
  } else {}
}

m1$analyses <- temp
dm <- dart.meta.data.merge(d3, m1)
fields <- c(analysis)
dms <- data.by.meta.fields(dm, fields, RandRbase, species, dataset, object = analysis)
# correct inconsistencies in J's scripts for naming treatment which is for naming files
treatment <- paste0("raw_SNPFilt_1SNPperClone_Field_", analysis)
dms$treatment <- treatment
gta <- dms$gt


kvalrange <- 1:10
real_kvalrange <- 2:7

LEA_dir <- paste0(outputloc,"/LEA_atleast_2_samples/")
if(!dir.exists(LEA_dir)){
  dir.create(LEA_dir)
}
tempdir <- paste0(outputloc, "/temp/")
if(!dir.exists(tempdir)){
  dir.create(tempdir)
}
##check for LEA folder and if it does not exist, run LEA and create map plot
lea_popdir <- paste0(maindir,species,"/popgen/",treatment,"/lea/")

#delete current directory if it exists
lea_popdir2 <- paste0(maindir,species,"/popgen/",treatment,"/")
if (dir.exists(lea_popdir2)) {
  unlink(lea_popdir2,recursive = TRUE)
  cat("Existing Directory has been deleted to create new popgen file")
}

if (!dir.exists(lea_popdir)) {
  nd_lea <- dart2lea(dms, RandRbase, species, dataset)
  snmf1=snmf(nd_lea, K=kvalrange, entropy = TRUE, repetitions = 1, project = "new")
  plot(snmf1, lwd = 6, col = "red", pch=1)
  
  require(cowplot)
  require(plyr)
  require(ggplot2)
  require(car)
  require(png)
  require(grid)
  require(gridExtra)
  require(tiff)
  require(oz)
  require(mapplots)
  
  lea_popdir <- paste0(maindir,species,"/popgen/",treatment,"/lea/",species,"_",dataset,".snmfProject")
  snmf_project <- load.snmfProject(lea_popdir)
  
  tiff(paste0(LEA_dir, species, " LEA Entropy.tiff"), units="in", width=12, height=6, res=300)
  par(mar=c(4,4,0,0)+0.1)
  plot(snmf_project, lwd = 6, col = "red", pch=1,)
  dev.off()
  
  ###this bit to tidy meta data - ignore pops with 1 indiv
    dmsanalysis <- dms$meta$analyses[,analysis]
    meta <- data.frame(table(dmsanalysis))
    z=meta$Freq<2
    
    if (sum(z, na.rm = TRUE)>=1) {#if there are any pops with less than 2 indiv...
      meta_lesthan2 <-as.character(meta[,1][meta$Freq<2]) #the pops with less than 2 indiv are...
      myvars <- dmsanalysis%in%meta_lesthan2
      newdmsanalysis <- dmsanalysis[!myvars]
      excluded_sample_names <- dms$sample_names[myvars]
      dms_d <- exclude.samples(dms, excluded_sample_names, remove_fixed_loci=TRUE)
      cat("metadata contains pops less than 2 individual","\n")
      print(meta_lesthan2)
    } else {
      meta_lesthan2 <-as.character(meta[,1][meta$Freq<2]) #the pops with less than 2 indiv are...
      myvars <- dmsanalysis%in%meta_lesthan2
      newdmsanalysis <- dms$meta$analyses[,analysis]
      dms_d <-dms
      excluded_sample_names <- NULL
      cat("metadata is alright without removal of pops","\n")
    }
    
    #making dataframe for  metadata without pops with less than 2 indiv
    dms_meta <- cbind.data.frame(dms$meta$sample_names, dms$meta$analyses[,analysis],dms$meta$lat, dms$meta$long)
    colnames(dms_meta) <- c("sample_names","site", "lat", "long")
    #pp <- dms_meta[!myvars,]
    #newdms_meta <- pp[apply(table(pp$site)>0,1,any),]
    newdms_meta <- dms_meta[!myvars,]
    colnames(newdms_meta) <- c("sample_names","site", "lat", "long")
    npop <-length(unique(newdms_meta$site))
    
    #making dataframe of averages (lat/longs) for metadata
    
    newmetalatlongs <- ddply(newdms_meta, c("site"),summarise,lat  = mean(lat),long  = mean(long), N=length(site))
  
  for (j in 1:length(real_kvalrange)) {
     kval <- real_kvalrange[j]
    if(kval==2) {
      kval_col <- c("red", "blue")    } 
    if (kval==3){
      kval_col <- c("red", "blue", "green")    } 
    if (kval==4){
      kval_col <- c("red", "blue", "green", "yellow")    } 
    if (kval==4){
      kval_col <- c("red", "blue", "green", "orange")    } 
    if (kval==5){
      kval_col <- c("red", "blue", "green", "orange", "purple")    } 
    if (kval==6){
      kval_col <- c("red", "blue", "green", "orange", "purple", "yellow")    } 
    if (kval==7){
      kval_col <- c("red", "blue", "green", "orange", "purple", "yellow", "pink")    } 
    if (kval>7){
      kval_col <- rev(matlab.like2(kval))   } 
    
    
    # kval_col <- PCAcols
    
    ###choose the best run
    ce           <- cross.entropy(snmf_project, K = kval)
    Rbest        <- which.min(ce) # with the lowest cross entropy criterion
    
    lea_popdirRbest <- paste0(maindir,RandRbase,species,"/popgen/",treatment,"/lea/",species,"_",dataset,".snmf/K",kval,"/run",Rbest,"/")
    longQfile <- paste0(lea_popdirRbest,species,"_",dataset,"_r",Rbest,".",kval,".Q")
    d_kval <- read.csv(longQfile, header=FALSE, sep=" ")
    
    ##### draw LEA pies on map
    d_kval_data <- data.frame(d_kval)
    dmsanalysis <- dms$meta$analyses[,analysis]
    #dmsanalysis2 <- substr(dmsanalysis,1,nchar(dmsanalysis)-11)
    meta <- data.frame(table(dmsanalysis))
    z=meta$Freq<2
    
    if (sum(z, na.rm = TRUE)>=1) { #if there are any pops with less than 2 indiv...
      d_kval_data3 <- cbind(d_kval_data[!myvars,],newdms_meta)
      d_kval2 <- d_kval[!myvars,] #original matrix without meta data for when making pies
      # cat(paste0("there are pops with less than 2 individuals","\n"))
    } else {
      d_kval_data3 <- cbind(d_kval_data,newdms_meta)
      # cat(paste0("pops are ok, more than 2 individuals","\n"))
      d_kval2 <- d_kval}
    
    #d_kval_data3 <- d_kval_data3[order(-d_kval_data3$lat),]
    d_kval_data3$site <- factor(d_kval_data3$site, levels = unique(d_kval_data3$site),ordered = TRUE)
    
    # d_kval3_dir <- paste0(maindir,LEA_dir,"/",species," LEA qmatrix k", kval,"_run",Rbest,"_plotted.csv")
    # write.table(d_kval_data3,d_kval3_dir,sep="\t",row.names = F)
    
    d_kval_data_all <- cbind(d_kval_data,dms_meta)
    d_all_dir <- paste0(LEA_dir,species," LEA qmatrix k", kval,"_run",Rbest,"_all.csv")
    write.table(d_kval_data_all,d_all_dir,sep="\t",row.names = F)
    
    cat(paste0("lea with meta info saved in ",d_all_dir ,"\n"))
    
    coord <- data.frame(d_kval_data3$long, d_kval_data3$lat)
    pop.factor <- factor(d_kval_data3$site)
    pop = as.numeric(pop.factor)
    qpop = matrix(NA, ncol = kval, nrow = npop)
    coord.pop = matrix(NA, ncol = 2, nrow = npop)
    
    #calculate aver age q values for each pop
    for (c in unique(pop)){
      qpop[c,] = apply(d_kval2[pop == c,], 2, mean)#apply(m,2,mean), the "2" means get the mean of the columns
      coord.pop[c,] = apply(coord[pop == c,], 2, mean)}
    
    ###save lea pop averaged qvalues
    leapoppie_output <- cbind(qpop,coord.pop)
    rownames(leapoppie_output)<-levels(d_kval_data3$site)
    leapiescoord_file <-paste0(lea_popdir, species," LEApiescoord k=",kval,".csv")
    write.table(leapoppie_output, leapiescoord_file)

    tiff(paste0(LEA_dir, species," LEA pies_bar k=",kval,".tiff"), units="in", width=14, height=18, res=300)
    par(mai=c(0,0,0,0))
    plot(1, type="n", axes=F, xlab="", ylab="")
    box()
    par(new = TRUE)
    par(fig=c(0,1,0.3,1)) #c(x1, x2, y1, y2)
    # par(fig=c(0.1,0.8,0.1,0.9)) #c(x1, x2, y1, y2)
    max_lat <- max(na.omit(newdms_meta$lat))
    
    row_sub_long = apply(data.frame(newdms_meta$long), 1, function(row) all(row !=0 ))
    newdmslong <- data.frame(newdms_meta$long)[row_sub_long,]
    row_sub_lat = apply(data.frame(newdms_meta$lat), 1, function(row) all(row !=0 ))
    newdmslat <- data.frame(newdms_meta$lat)[row_sub_lat,]
    divxlims <- c(min(na.omit(newdmslong))-1,max(na.omit(newdmslong))+1)
    divylims <- c(min(na.omit(newdmslat))-1,max(na.omit(newdmslat))+1)
    
    if (max_lat >= -20) {
      oz(xlim=divxlims, ylim=divylims)
      #text(min(na.omit(qmatrix_data3$long))-0.7, max(na.omit(qmatrix_data3$lat))+1.2, paste0(" k=", kval),font=c(3,2),cex=1)
    } else {
      if ( max_lat >= -28 | max_lat >= -35) { #OR
        oz(xlim=divxlims, ylim=divylims)
        # text(min(na.omit(qmatrix_data3$long))-0.1, max(na.omit(qmatrix_data3$lat))-0.1, paste0(" k=", kval),font=c(3,2),cex=1)

      } else {
        oz(xlim=divxlims, ylim=divylims)
        #text(min(na.omit(qmatrix_data3$long))-0.7, max(na.omit(qmatrix_data3$lat))+1.2, paste0(" k=", kval),font=c(3,2),cex=1)
      }
    }
    #draw pies on map
    for (n in 1:npop){
      add.pie(z = qpop[n,], x = coord.pop[n,1], y = coord.pop[n,2],
              radius=(divylims[2]-divylims[1])/40,
              labels = "",
              density = 70, # removing density fills pie slices
              col = kval_col)
      }
  
    d_kval <- read.csv(longQfile, header=FALSE, sep=" ")
    d_kval2 <-d_kval[!myvars,]
    d_kval2$lat <- newdms_meta$lat
    d_kval2 <- d_kval2[order(-d_kval2$lat),]
    d_kval2$lat <- NULL
    colnames(d_kval2) <- NULL
    rownames(d_kval2) <- NULL
    coordinate_cities <- read.csv(paste0(maindir,"CityMapCordinates.csv"))
    points(x=coordinate_cities$long, y= coordinate_cities$lat, pch=17, cex=1.5)
    text(x=coordinate_cities$long, y= coordinate_cities$lat, label=coordinate_cities$City, pos=4, font=2, cex=1)
    newdms_meta <- newdms_meta[order(-newdms_meta$lat),]
    
    #####now drawing barplot
    # par(fig=c(0.6,0.9,0.1,0.3),new = TRUE) #c(x1, x2, y1, y2)
    par(fig=c(0.1,1,0.05,0.28),new = TRUE) #c(x1, x2, y1, y2)
    barplot(t(as.matrix(d_kval2)), col=kval_col, border = T, space = 0, xlab = "Individuals", ylab = "Admixture coefficients",
            names.arg=newdms_meta$site, las=2,cex.names=0.4)
    mtext(paste0("      k= ", kval), side=3, adj=1, line=1, font=2, cex=6)
    dev.off()
    # par(resetPar())
  } 
  
  #####draw barplots, pies and save qmatrix csv files
  tiff(paste0(LEA_dir,species," samp_distrib.tiff"),
       units="in", width=6, height=10, res=200)
  #par(mar=c(0,2,0,0)) #c(bottom, left, top, right)
  # par(mai=c(0.1,0.1,0.2,0.1))
  par(mai=c(0.1, 0.1, 0.1, 0.1))
  plot(1, type="n", axes=F, xlab="", ylab="")
  par(new = TRUE)
  # par(fig=c(0.1,0.8,0.1,0.9)) #c(x1, x2, y1, y2)
  par(fig=c(0, 1, 0, 1)) #c(x1, x2, y1, y2)
  oz(xlim=divxlims, ylim=divylims, lwd=0.5)
  coord.pop <- data.frame(coord.pop)
  colnames(coord.pop) <- c("long", "lat")
  coord.pop$site <- sort(unique(d_kval_data3$site))
  
  if (file.exists(paste0(LEA_dir,"site_number.csv"))) {
    rep2 <- read.csv(paste0(LEA_dir,"site_number.csv"),sep=" ")
  } else {
    dms_meta <- cbind.data.frame(dms$meta$sample_names, dms$meta$analyses[,analysis],dms$meta$lat, dms$meta$long)
    colnames(dms_meta) <- c("sample_names","site", "lat", "long")
    dms_meta$lat <- as.numeric(dms_meta$lat)
    dms_meta$long <- as.numeric(dms_meta$long)
    ###
    dms_meta2 <- dms_meta[order(dms_meta$site), ]
    dms_meta2 <- dms_meta2[order(-dms_meta2$lat), ]
    ####make sure your rep has nums_n_site which is in df2!!!
    dms_meta2$site <- factor(dms_meta2$site, levels = unique(dms_meta2$site),ordered = TRUE)
    
    rep <- ddply(dms_meta2,
                 .(site),
                 summarise,
                 lat  = mean(lat),
                 long = mean(long))
    rep2 <- rep[order(-rep$lat), ]
    x <- length(unique(dms_meta2$site))
    bgcols = rev(matlab.like2(x))
    rep2$nums_n_site <- paste0(1:x, ": ", rep2$site)
    rep2$num_site <- 1:x
    write.table(rep2, paste0(LEA_dir,"site_number_fromLEA.csv"))
  }
  
  coord.pop11 <- merge(coord.pop,rep2[,c("site","nums_n_site","num_site")],by="site")
  coord.pop11$lab <- coord.pop11$num_site
  points(coord.pop11$long,coord.pop11$lat,pch=21,lwd=0.2,bg="blue",cex=1)
  #text(x=150.8758, y=-34.40813+0.1,labels="2 to 29",pos=3,cex=0.7,font=4)
  coord.pop11 <- coord.pop11[complete.cases(coord.pop11), ]#remove latlongs that are NAs
  pointLabel(coord.pop11$long, coord.pop11$lat, labels = paste("  ", coord.pop11$lab, "  ", sep=""), cex=0.7, font=4,offset=0.25)
  
  if(exists("excluded_sample_names")) {
    excluded1pops <- dms$meta$analyses[,analysis][(dms$meta$sample_names%in%excluded_sample_names)]
    # cat("samples from 1indiv pop excluded in map","\n")
  }else {
    excluded1pops <- NULL
  }
  
  if(sum(is.na(newdms_meta$lat))>0){
    excludedNApops <- as.character(unique(newdms_meta$site[is.na(newdms_meta$lat)])) #pops excluded from map because latlongs are NA
    # cat("samples wihtout latlongs excluded in map","\n")
  } else {
    excludedNApops <- NULL
  }
  
  if(exists("excludedNApops") | exists("excluded1pops")){
    par(fig=c(0.1, 0.9, 0.1, 0.9),new = TRUE) #c(x1, x2, y1, y2)
    # par(fig=c(0,0.5,0,0.5),new = TRUE) #c(x1, x2, y1, y2)
    all_excluded <- c(excluded1pops,excludedNApops)
    rep2$site <- as.factor(rep2$site)
    all_ex_nums <-  rep2$nums_n_site[(rep2$site%in%all_excluded)]
    all_ex_nums2 <- paste0(all_ex_nums, "\n",collapse = "")
    all_ex_nums3 <- paste0("Not included:\n",all_ex_nums2)
    mtext( all_ex_nums3,cex = 0.7,side=1,adj=0)
    dev.off()
  } else {
    dev.off()
  }

  ##################TABLE
  coord.pop2 <- cbind.data.frame(coord.pop11$lab,coord.pop11$site,coord.pop11$lat,coord.pop11$long)
  colnames(coord.pop2)<-c("","site","lat","long")
  coord.pop2 <- coord.pop2[order(-coord.pop2$lat),]
  tg = gridExtra::tableGrob(coord.pop2,rows = NULL)
  h = grid::convertHeight(sum(tg$heights), "in", TRUE)+0.5
  w = grid::convertWidth(sum(tg$widths), "in", TRUE)+0.5
  ggplot2::ggsave(paste0(LEA_dir,species," samp_distrib2.tiff"), tg, width=w, height=h)
  # cat(paste0("samp_distrib2 figure done!","\n"))
  
  ############combine with sampl_dist
  setwd(outputloc)
  loc_1<- dir(path = "LEA", pattern = "samp_distrib2.tiff$", full.names = TRUE, recursive = TRUE) ###must set working directory to make readPNG work
  loc_2 <- dir(path = "LEA", pattern = "samp_distrib.tiff$", full.names = TRUE, recursive = TRUE) ###must set working directory to make readPNG work
  loc_3 <- dir(path = "LEA", pattern = "Entropy.tiff$", full.names = TRUE, recursive = TRUE) ###must set working directory to make readPNG work
  lloo <- c(loc_3)
  plots <- lapply(ll <- lloo,function(x){
    img <- as.raster(readTIFF(x)) ##change this to readTIFF if readPNG doesnt work
    rasterGrob(img, interpolate = FALSE)
  })
  #lay <- rbind(c(2,2,1,1),c(2,2,1,1),c(2,2,1,1),c(3,3,3,3))
  lay <- rbind(c(1,1))
  ggsave(paste0("LEA/",species, " distrib_entrp.tiff"),
         width=8, height=10, units = "in", device='tiff',dpi = 300,
         marrangeGrob(grobs = plots, layout_matrix = lay,top=textGrob(paste0(analysis," LEA outputs"), gp=gpar(fontsize=20,fontface=4))))
  
  ####combine all pies
  piedatalist=list()
  pie_ls_files <- c("LEA pies_bar k=2.tiff$", "LEA pies_bar k=3.tiff$", "LEA pies_bar k=4.tiff$","LEA pies_bar k=5.tiff$","LEA pies_bar k=6.tiff$","LEA pies_bar k=7.tiff$")
  
  for (f in 1:length(pie_ls_files)) {
    ll_filenamesleapies<- dir(path = "LEA", pattern = pie_ls_files[f], full.names = TRUE, recursive = TRUE) ###must set working directory to make readPNG work
    piedatalist[[f]] <- ll_filenamesleapies
  }
  
  pie_plots <- lapply(ll <- piedatalist,function(x){
    img <- as.raster(readTIFF(x)) ##change this to readTIFF if readPNG doesnt work
    rasterGrob(img, interpolate = FALSE)
  })
  
  
  lay <- rbind(c(1,2,3),c(4,5,6)) #for kvals from 2 to 7
  ggsave(paste0("LEA/",species, " LEA pies_bar_ALL.tiff"),
         width=16, height=14,dpi = 300, units = "in", device='tiff',
         marrangeGrob(grobs = pie_plots, layout_matrix = lay,top=NULL))
  # cat(paste0("combined all lea piebar plots!","\n"))
  
  ###combine LEA figures, pies for k equals 2,3,4,5, sampling distrib, entropy into a single figure
  datalist=list()
  
  ls_files <- c("distrib_entrp.tiff$",
                "LEA pies_bar_ALL.tiff$")
  
  for (g in 1:length(ls_files)) {
    ll_filenames <- dir(path = "LEA_atleast_2_samples", pattern = ls_files[g], full.names = TRUE, recursive = TRUE) ###must set working directory to make readPNG work
    datalist[[g]] <- ll_filenames
  }
  locsLEA <- unlist(datalist)
  
  plots <- lapply(ll <- locsLEA,function(x){
    img <- as.raster(readTIFF(x)) ##change this to readTIFF if readPNG doesnt work
    rasterGrob(img, interpolate = FALSE)
  })
  
  #layoo <- rbind(c(1,2))
  layoo <- rbind(c(1,2,2,2),c(1,2,2,2),c(1,2,2,2),c(1,2,2,2))
  
  ggsave(paste0(LEA_dir,species, " lea k2345loc bar test.tiff"),
         width=12, height=8.3,dpi = 600, units = "in", device='tiff',
         marrangeGrob(grobs = plots, layout_matrix = layoo,top=NULL))
  cat(paste0("all LEA plots done!","\n"))
  
  setwd(maindir)
  
} else {
  print("remove LEA folder under popgen folder first!!!")
  #or if you have already run LEA, and just want to replot everything?
}




####draw pretty structure pies with stamen maps for a particualr kval
###this bit to tidy meta data - ignore pops with 1 indiv
dmsanalysis <- dms$meta$analyses[,analysis]
meta <- data.frame(table(dmsanalysis))
z=meta$Freq<2

if (sum(z, na.rm = TRUE)>=1) {#if there are any pops with less than 2 indiv...
  meta_lesthan2 <-as.character(meta[,1][meta$Freq<2]) #the pops with less than 2 indiv are...
  myvars <- dmsanalysis%in%meta_lesthan2
  newdmsanalysis <- dmsanalysis[!myvars]
  excluded_sample_names <- dms$sample_names[myvars]
  dms_d <- exclude.samples(dms, excluded_sample_names, remove_fixed_loci=TRUE)
  cat("metadata contains pops less than 2 individual","\n")
  print(meta_lesthan2)
} else {
  meta_lesthan2 <-as.character(meta[,1][meta$Freq<2]) #the pops with less than 2 indiv are...
  myvars <- dmsanalysis%in%meta_lesthan2
  newdmsanalysis <- dms$meta$analyses[,analysis]
  dms_d <-dms
  excluded_sample_names <- NULL
  cat("metadata is alright without removal of pops","\n")
}

#making dataframe for  metadata without pops with less than 2 indiv
dms_meta <- cbind.data.frame(dms$meta$sample_names, dms$meta$analyses[,analysis],dms$meta$lat, dms$meta$long)
colnames(dms_meta) <- c("sample_names","site", "lat", "long")
#pp <- dms_meta[!myvars,]
#newdms_meta <- pp[apply(table(pp$site)>0,1,any),]
newdms_meta <- dms_meta[!myvars,]
colnames(newdms_meta) <- c("sample_names","site", "lat", "long")
npop <-length(unique(newdms_meta$site))

#making dataframe of averages (lat/longs) for metadata

newmetalatlongs <- ddply(newdms_meta, c("site"),summarise,lat  = mean(lat),long  = mean(long), N=length(site))
kval = 2 # change to whatever kval you want to map
if(kval==2) {
  kval_col <- c("red", "blue"); finalcolname <- c("long", "lat","sites","A","B")} 
if (kval==3){
  kval_col <- c("red", "blue", "green"); finalcolname <- c("long", "lat","sites","A","B","C")    } 
if (kval==4){
  kval_col <- c("red", "blue", "green", "yellow"); finalcolname <- c("long", "lat","sites","A","B","C", "D")     } 
if (kval==4){
  kval_col <- c("red", "blue", "green", "orange"); finalcolname <- c("long", "lat","sites","A","B","C", "D", "E")     } 
if (kval==5){
  kval_col <- c("red", "blue", "green", "orange", "purple"); finalcolname <- c("long", "lat","sites","A","B","C", "D", "E", "F")     } 
if (kval==6){
  kval_col <- c("red", "blue", "green", "orange", "purple", "yellow"); finalcolname <- c("long", "lat","sites","A","B","C", "D", "E", "F", "G")     } 
if (kval==7){
  kval_col <- c("red", "blue", "green", "orange", "purple", "yellow", "pink"); finalcolname <- c("long", "lat","sites","A","B","C")     } 
if (kval>7){
  kval_col <- rev(matlab.like2(kval))   } 


###choose the best run
ce           <- cross.entropy(snmf_project, K = kval)
Rbest        <- which.min(ce) # with the lowest cross entropy criterion
lea_popdirRbest <- paste0(maindir,RandRbase,species,"/popgen/",treatment,"/lea/",species,"_",dataset,".snmf/K",kval,"/run",Rbest,"/")
longQfile <- paste0(lea_popdirRbest,species,"_",dataset,"_r",Rbest,".",kval,".Q")
d_kval <- read.csv(longQfile, header=FALSE, sep=" ")

##### draw LEA pies on map
d_kval_data <- data.frame(d_kval)
dmsanalysis <- dms$meta$analyses[,analysis]
#dmsanalysis2 <- substr(dmsanalysis,1,nchar(dmsanalysis)-11)
meta <- data.frame(table(dmsanalysis))
z=meta$Freq<2

if (sum(z, na.rm = TRUE)>=1) { #if there are any pops with less than 2 indiv...
  d_kval_data3 <- cbind(d_kval_data[!myvars,],newdms_meta)
  d_kval2 <- d_kval[!myvars,] #original matrix without meta data for when making pies
  # cat(paste0("there are pops with less than 2 individuals","\n"))
} else {
  d_kval_data3 <- cbind(d_kval_data,newdms_meta)
  # cat(paste0("pops are ok, more than 2 individuals","\n"))
  d_kval2 <- d_kval}

#d_kval_data3 <- d_kval_data3[order(-d_kval_data3$lat),]
d_kval_data3$site <- factor(d_kval_data3$site, levels = unique(d_kval_data3$site),ordered = TRUE)
d_kval_data_all <- cbind(d_kval_data,dms_meta)
d_all_dir <- paste0(LEA_dir,species," LEA qmatrix k", kval,"_run",Rbest,"_all.csv")
write.table(d_kval_data_all,d_all_dir,sep="\t",row.names = F)

cat(paste0("lea with meta info saved in ",d_all_dir ,"\n"))

coord <- data.frame(d_kval_data3$long, d_kval_data3$lat)
pop.factor <- factor(d_kval_data3$site)
pop = as.numeric(pop.factor)
qpop = matrix(NA, ncol = kval, nrow = npop)
coord.pop = matrix(NA, ncol = 2, nrow = npop)

#calculate aver age q values for each pop
for (c in unique(pop)){
  qpop[c,] = apply(d_kval2[pop == c,], 2, mean)#apply(m,2,mean), the "2" means get the mean of the columns
  coord.pop[c,] = apply(coord[pop == c,], 2, mean)
}

coord.pop <- cbind(coord.pop, unique(pop))


###save lea pop averaged qvalues
leapoppie_output <- cbind(qpop,coord.pop)
rownames(leapoppie_output)<-levels(d_kval_data3$site)
leapiescoord_file <-paste0(lea_popdir, species," LEApiescoord k=",kval,".csv")
write.table(leapoppie_output, leapiescoord_file)

max_lat <- max(na.omit(newdms_meta$lat))
row_sub_long = apply(data.frame(newdms_meta$long), 1, function(row) all(row !=0 ))
newdmslong <- data.frame(newdms_meta$long)[row_sub_long,]
row_sub_lat = apply(data.frame(newdms_meta$lat), 1, function(row) all(row !=0 ))
newdmslat <- data.frame(newdms_meta$lat)[row_sub_lat,]

divxlims <- c(min(na.omit(newdmslong))-1,max(na.omit(newdmslong))+1)
divylims <- c(min(na.omit(newdmslat))-1,max(na.omit(newdmslat))+1)
divxlimsrange <- abs(c(min(na.omit(newdmslong))-max(na.omit(newdmslong))))
divylimsrange <- abs(c(min(na.omit(newdmslat))-max(na.omit(newdmslat))))

if(divylimsrange > 5) {gmapzoom= 4}
if(divylimsrange > 4 && divylimsrange < 5)  {gmapzoom= 5}
if(divylimsrange > 3 && divylimsrange < 4)  {gmapzoom= 6}
if(divylimsrange > 2 && divylimsrange < 3)  {gmapzoom= 7}
if(divylimsrange > 1 && divylimsrange < 2)  {gmapzoom= 8}
if(divylimsrange > 0.2 && divylimsrange < 1) {gmapzoom= 10}
if(divylimsrange < 0.2) {gmapzoom= 12}

register_stadiamaps(key = "4fe26ca2-aa51-4e9e-8151-dca93798847e")

stadiamapbase = get_stadiamap(bbox = c(left = divxlims[1], bottom = divylims[1], right = divxlims[2], top=divylims[2]), zoom = 11, maptype = "stamen_terrain", crop=F)




cnames <- c()
for (b in 1:ncol(qpop)){
  cnames[b] <- b
}
qpop <- data.frame(qpop)
colnames(qpop) <- cnames
finalpie <- cbind(coord.pop, qpop)
finalpie <- data.frame(finalpie)
colnames(finalpie) <- finalcolname
finalpie$radius <- 6 * abs(rnorm(nrow(finalpie)))



tiff(paste0(LEA_dir, species," LEA Stamen Map k=",kval,".tiff"), units="in", width=12, height=18, res=300)

ggmap(stadiamapbase) + coord_fixed() + 
  geom_scatterpie(data=finalpie, cols=LETTERS[1:2], aes(x= long, y = lat, group=sites, r=(divylims[2]-divylims[1])/60), alpha=0.8) +
  theme(legend.position="none") +
  scale_fill_manual(values=kval_col) +
  ggsn::scalebar(data=finalpie, dist = 5, dist_unit = "km", location = "bottomleft", 
                 st.bottom = T, st.size = 5, st.dist = 0.01, border.size = 0.1,
                 transform = TRUE, model = "WGS84", height = 0.015, anchor = c(x = min(finalpie$long)-0.02, y = min(finalpie$lat)-0.01)) +
  xlim(min(finalpie$long)-0.02, max(finalpie$long)+0.02)+ # increase to change width of plot
  ylim(min(finalpie$lat)-0.02, max(finalpie$lat)+0.02)+
  annotation_north_arrow(location = "tr", which_north = "true", pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"), style = north_arrow_fancy_orienteering)
  # geom_point(data=finalpie, aes(x= long , y=lat)) +
  # geom_label_repel(data=finalpie, aes(x= long , y=lat, label = sites), fontface="bold", size = 2.5)


dev.off()




d_kval <- read.csv(longQfile, header=FALSE, sep=" ")
d_kval2 <-d_kval[!myvars,]
d_kval2$lat <- newdms_meta$lat
d_kval2 <- d_kval2[order(-d_kval2$lat),]
d_kval2$lat <- NULL
colnames(d_kval2) <- NULL
rownames(d_kval2) <- NULL

newdms_meta <- dms_meta[!myvars,]
newdms_meta2 <- newdms_meta[order(-newdms_meta$lat),]




tiff(paste0(LEA_dir, species," LEA Barplot Only k=",kval,".tiff"), units="in", width=14, height=8, res=300)
par(fig=c(0,1,0.2,1),new = TRUE) #c(x1, x2, y1, y2)
barplot(t(as.matrix(d_kval2)), col=kval_col, border = T, space = 0, xlab = "", ylab = "Admixture coefficients",
        names.arg=newdms_meta2$site, las=2,cex.names=0.5)
# mtext(paste0("      k= ", kval), side=3, adj=1, line=1, font=2, cex=6)
dev.off()









##### LEA - Entropy, Barplot and Piechart limiting X samples per pop #####

#read metadata file
metafile <- read.xlsx(paste0(maindir, species,"/meta/",species,"_",dataset,"_meta.xlsx"))
m1 <- read.meta.data(d3, RandRbase, species, dataset, fields = as.numeric(ncol(metafile)-4))
colnames(m1$analyses)
# create a temp meta dataframe to remove sites/samples that include less than X individuals
temp <- as.data.frame(m1$analyses)

samplethreshold <- 5 # need a minimum of 2 samples per pop to run LEA

##### Below lines remove sites with less than threshold number of samples, and subsample sites that are above the threshold
for (x in unique(na.omit(temp[, analysis]))) {
  pop_samples <- which(temp[, analysis] == x)
  if (length(pop_samples) < samplethreshold) {
    temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
    print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
  }
  if (length(pop_samples) == samplethreshold) {
    print(paste0("Keeping ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
  }
  if (length(pop_samples) > samplethreshold) {
    subsamplepops <- pop_samples[sample(1:length(pop_samples), size=samplethreshold, replace = F)]
    for (b in 1:length(pop_samples)){
      if (pop_samples[b] %in% subsamplepops == FALSE){
        temp[, analysis][pop_samples[b]] <- NA
      }
      print(paste0("Reducing ",x, " to a sample size of ", samplethreshold," samples"))
    }
  } else {}
}

##### If you want to just remove any sites that are below the threshold (but keep individuals that are above a threshold, then unhash the below section)
# for (x in unique(na.omit(temp[, analysis]))) {
#   pop_samples <- which(temp[, analysis] == x)
#   if (length(pop_samples) < samplethreshold) {
#     temp[, analysis] <- gsub(x, replacement = NA, temp[, analysis])
#     print(paste0("Removing ", x, " with ", " with a sample size of ", length(pop_samples), " samples"))
#   } else {}
# }

m1$analyses <- temp
dm <- dart.meta.data.merge(d3, m1)
fields <- c(analysis)
dms <- data.by.meta.fields(dm, fields, RandRbase, species, dataset, object = analysis)
# correct inconsistencies in J's scripts for naming treatment which is for naming files
treatment <- paste0("raw_SNPFilt_1SNPperClone_Field_", analysis)
dms$treatment <- treatment
gta <- dms$gt



kvalrange <- 1:25
real_kvalrange <- 2:7

LEA_dir <- paste0(outputloc,"/LEA_limiting_",samplethreshold, "_samples/")
if(!dir.exists(LEA_dir)){
  dir.create(LEA_dir)
}
tempdir <- paste0(outputloc, "/temp/")
if(!dir.exists(tempdir)){
  dir.create(tempdir)
}
##check for LEA folder and if it does not exist, run LEA and create map plot
lea_popdir <- paste0(maindir,species,"/popgen/",treatment,"/lea/")

#delete current directory if it exists
lea_popdir2 <- paste0(maindir,species,"/popgen/",treatment,"/")
if (dir.exists(lea_popdir2)) {
  unlink(lea_popdir2,recursive = TRUE)
  cat("Existing Directory has been deleted to create new popgen file")
}

if (!dir.exists(lea_popdir)) {
  nd_lea <- dart2lea(dms, RandRbase, species, dataset)
  snmf1=snmf(nd_lea, K=kvalrange, entropy = TRUE, repetitions = 1, project = "new")
  plot(snmf1, lwd = 6, col = "red", pch=1)
  
  require(cowplot)
  require(plyr)
  require(ggplot2)
  require(car)
  require(png)
  require(grid)
  require(gridExtra)
  require(tiff)
  require(oz)
  require(mapplots)
  
  lea_popdir <- paste0(maindir,species,"/popgen/",treatment,"/lea/",species,"_",dataset,".snmfProject")
  snmf_project <- load.snmfProject(lea_popdir)
  
  tiff(paste0(LEA_dir,species, " LEA Entropy.tiff"), units="in", width=12, height=6, res=300)
  par(mar=c(4,4,0,0)+0.1)
  plot(snmf_project, lwd = 6, col = "red", pch=1,)
  dev.off()
  
  ###this bit to tidy meta data - ignore pops with 1 indiv
    dmsanalysis <- dms$meta$analyses[,analysis]
    meta <- data.frame(table(dmsanalysis))
    z=meta$Freq<2
    
    if (sum(z, na.rm = TRUE)>=1) {#if there are any pops with less than 2 indiv...
      meta_lesthan2 <-as.character(meta[,1][meta$Freq<2]) #the pops with less than 2 indiv are...
      myvars <- dmsanalysis%in%meta_lesthan2
      newdmsanalysis <- dmsanalysis[!myvars]
      excluded_sample_names <- dms$sample_names[myvars]
      dms_d <- exclude.samples(dms, excluded_sample_names, remove_fixed_loci=TRUE)
      cat("metadata contains pops less than 2 individual","\n")
      print(meta_lesthan2)
    } else {
      meta_lesthan2 <-as.character(meta[,1][meta$Freq<2]) #the pops with less than 2 indiv are...
      myvars <- dmsanalysis%in%meta_lesthan2
      newdmsanalysis <- dms$meta$analyses[,analysis]
      dms_d <-dms
      excluded_sample_names <- NULL
      cat("metadata is alright without removal of pops","\n")
    }
    
    #making dataframe for  metadata without pops with less than 2 indiv
    dms_meta <- cbind.data.frame(dms$meta$sample_names, dms$meta$analyses[,analysis],dms$meta$lat, dms$meta$long)
    colnames(dms_meta) <- c("sample_names","site", "lat", "long")
    #pp <- dms_meta[!myvars,]
    #newdms_meta <- pp[apply(table(pp$site)>0,1,any),]
    newdms_meta <- dms_meta[!myvars,]
    colnames(newdms_meta) <- c("sample_names","site", "lat", "long")
    npop <-length(unique(newdms_meta$site))
    
    #making dataframe of averages (lat/longs) for metadata
    newmetalatlongs <- ddply(newdms_meta, c("site"),summarise,lat  = mean(lat),long  = mean(long), N=length(site))

  
  for (j in 1:length(real_kvalrange)) {
    
    kval <- real_kvalrange[j]
    if(kval==2) {
      kval_col <- c("red", "blue")    } 
    if (kval==3){
      kval_col <- c("red", "blue", "green")    } 
    if (kval==4){
      kval_col <- c("red", "blue", "green", "yellow")    } 
    if (kval==4){
      kval_col <- c("red", "blue", "green", "orange")    } 
    if (kval==5){
      kval_col <- c("red", "blue", "green", "orange", "purple")    } 
    if (kval==6){
      kval_col <- c("red", "blue", "green", "orange", "purple", "yellow")    } 
    if (kval==7){
      kval_col <- c("red", "blue", "green", "orange", "purple", "yellow", "pink")    } 
    if (kval>7){
      kval_col <- rev(matlab.like2(kval))   } 
    
    
    # kval_col <- PCAcols
    
    ###choose the best run
    ce           <- cross.entropy(snmf_project, K = kval)
    Rbest        <- which.min(ce) # with the lowest cross entropy criterion
    
    lea_popdirRbest <- paste0(maindir,RandRbase,species,"/popgen/",treatment,"/lea/",species,"_",dataset,".snmf/K",kval,"/run",Rbest,"/")
    longQfile <- paste0(lea_popdirRbest,species,"_",dataset,"_r",Rbest,".",kval,".Q")
    d_kval <- read.csv(longQfile, header=FALSE, sep=" ")
    
    ##### draw LEA pies on map
    d_kval_data <- data.frame(d_kval)
    dmsanalysis <- dms$meta$analyses[,analysis]
    #dmsanalysis2 <- substr(dmsanalysis,1,nchar(dmsanalysis)-11)
    meta <- data.frame(table(dmsanalysis))
    z=meta$Freq<2
    
    if (sum(z, na.rm = TRUE)>=1) { #if there are any pops with less than 2 indiv...
      d_kval_data3 <- cbind(d_kval_data[!myvars,],newdms_meta)
      d_kval2 <- d_kval[!myvars,] #original matrix without meta data for when making pies
      # cat(paste0("there are pops with less than 2 individuals","\n"))
    } else {
      d_kval_data3 <- cbind(d_kval_data,newdms_meta)
      # cat(paste0("pops are ok, more than 2 individuals","\n"))
      d_kval2 <- d_kval}
    
    #d_kval_data3 <- d_kval_data3[order(-d_kval_data3$lat),]
    d_kval_data3$site <- factor(d_kval_data3$site, levels = unique(d_kval_data3$site),ordered = TRUE)
    
    # d_kval3_dir <- paste0(maindir,LEA_dir,"/",species," LEA qmatrix k", kval,"_run",Rbest,"_plotted.csv")
    # write.table(d_kval_data3,d_kval3_dir,sep="\t",row.names = F)
    
    d_kval_data_all <- cbind(d_kval_data,dms_meta)
    d_all_dir <- paste0(LEA_dir,species," LEA qmatrix k", kval,"_run",Rbest,"_all.csv")
    write.table(d_kval_data_all,d_all_dir,sep="\t",row.names = F)
    
    cat(paste0("lea with meta info saved in ",d_all_dir ,"\n"))
    
    coord <- data.frame(d_kval_data3$long, d_kval_data3$lat)
    pop.factor <- factor(d_kval_data3$site)
    pop = as.numeric(pop.factor)
    qpop = matrix(NA, ncol = kval, nrow = npop)
    coord.pop = matrix(NA, ncol = 2, nrow = npop)
    
    #calculate aver age q values for each pop
    for (c in unique(pop)){
      qpop[c,] = apply(d_kval2[pop == c,], 2, mean)#apply(m,2,mean), the "2" means get the mean of the columns
      coord.pop[c,] = apply(coord[pop == c,], 2, mean)}
    
    ###save lea pop averaged qvalues
    leapoppie_output <- cbind(qpop,coord.pop)
    rownames(leapoppie_output)<-levels(d_kval_data3$site)
    leapiescoord_file <-paste0(lea_popdir, species," LEApiescoord k=",kval,".csv")
    write.table(leapoppie_output, leapiescoord_file)
    
    
    tiff(paste0(LEA_dir, species," LEA pies_bar k=",kval,".tiff"), units="in", width=14, height=18, res=300)
    par(mai=c(0,0,0,0))
    plot(1, type="n", axes=F, xlab="", ylab="")
    box()
    par(new = TRUE)
    par(fig=c(0,1,0.3,1)) #c(x1, x2, y1, y2)
    # par(fig=c(0.1,0.8,0.1,0.9)) #c(x1, x2, y1, y2)
    max_lat <- max(na.omit(newdms_meta$lat))
    
    row_sub_long = apply(data.frame(newdms_meta$long), 1, function(row) all(row !=0 ))
    newdmslong <- data.frame(newdms_meta$long)[row_sub_long,]
    row_sub_lat = apply(data.frame(newdms_meta$lat), 1, function(row) all(row !=0 ))
    newdmslat <- data.frame(newdms_meta$lat)[row_sub_lat,]
    
    divxlims <- c(min(na.omit(newdmslong))-1,max(na.omit(newdmslong))+1)
    divylims <- c(min(na.omit(newdmslat))-1,max(na.omit(newdmslat))+1)
    
    
    
    if (max_lat >= -20) {
      oz(xlim=divxlims, ylim=divylims)
      #text(min(na.omit(qmatrix_data3$long))-0.7, max(na.omit(qmatrix_data3$lat))+1.2, paste0(" k=", kval),font=c(3,2),cex=1)
    } else {
      if ( max_lat >= -28 | max_lat >= -35) { #OR
        oz(xlim=divxlims, ylim=divylims)
        # text(min(na.omit(qmatrix_data3$long))-0.1, max(na.omit(qmatrix_data3$lat))-0.1, paste0(" k=", kval),font=c(3,2),cex=1)
        
      } else {
        oz(xlim=divxlims, ylim=divylims)
        #text(min(na.omit(qmatrix_data3$long))-0.7, max(na.omit(qmatrix_data3$lat))+1.2, paste0(" k=", kval),font=c(3,2),cex=1)
      }
    }
    
    
    
    #draw pies on map
    for (n in 1:npop){
      add.pie(z = qpop[n,], x = coord.pop[n,1], y = coord.pop[n,2],
              radius=(divylims[2]-divylims[1])/40,
              labels = "",
              density = 70, # removing density fills pie slices
              col = kval_col)
    }
    
    d_kval <- read.csv(longQfile, header=FALSE, sep=" ")
    d_kval2 <-d_kval[!myvars,]
    d_kval2$lat <- newdms_meta$lat
    d_kval2 <- d_kval2[order(-d_kval2$lat),]
    d_kval2$lat <- NULL
    colnames(d_kval2) <- NULL
    rownames(d_kval2) <- NULL
    
    coordinate_cities <- read.csv(paste0(maindir,"CityMapCordinates.csv"))
    points(x=coordinate_cities$long, y= coordinate_cities$lat, pch=17, cex=1.5)
    text(x=coordinate_cities$long, y= coordinate_cities$lat, label=coordinate_cities$City, pos=4, font=2, cex=1)
    
    
    newdms_meta <- newdms_meta[order(-newdms_meta$lat),]
    
    #####now drawing barplot
    # par(fig=c(0.6,0.9,0.1,0.3),new = TRUE) #c(x1, x2, y1, y2)
    par(fig=c(0.1,1,0.05,0.28),new = TRUE) #c(x1, x2, y1, y2)
    barplot(t(as.matrix(d_kval2)), col=kval_col, border = T, space = 0, xlab = "Individuals", ylab = "Admixture coefficients",
            names.arg=newdms_meta$site, las=2,cex.names=0.4)
    mtext(paste0("      k= ", kval), side=3, adj=1, line=1, font=2, cex=6)
    dev.off()
    # par(resetPar())
  
    
    
      
    
  } 

#draw barplots, pies and save qmatrix csv files
  
  ####draw sampling distribution with numbering
  tiff(paste0(LEA_dir,species," samp_distrib.tiff"),
       units="in", width=6, height=10, res=200)
  #par(mar=c(0,2,0,0)) #c(bottom, left, top, right)
  # par(mai=c(0.1,0.1,0.2,0.1))
  par(mai=c(0.1, 0.1, 0.1, 0.1))
  plot(1, type="n", axes=F, xlab="", ylab="")
  par(new = TRUE)
  # par(fig=c(0.1,0.8,0.1,0.9)) #c(x1, x2, y1, y2)
  par(fig=c(0, 1, 0, 1)) #c(x1, x2, y1, y2)
  oz(xlim=divxlims, ylim=divylims, lwd=0.5)
  coord.pop <- data.frame(coord.pop)
  colnames(coord.pop) <- c("long", "lat")
  coord.pop$site <- sort(unique(d_kval_data3$site))
  
  if (file.exists(paste0(LEA_dir,"site_number.csv"))) {
    rep2 <- read.csv(paste0(LEA_dir,"site_number.csv"),sep=" ")
  } else {
    dms_meta <- cbind.data.frame(dms$meta$sample_names, dms$meta$analyses[,analysis],dms$meta$lat, dms$meta$long)
    colnames(dms_meta) <- c("sample_names","site", "lat", "long")
    dms_meta$lat <- as.numeric(dms_meta$lat)
    dms_meta$long <- as.numeric(dms_meta$long)
    ###
    dms_meta2 <- dms_meta[order(dms_meta$site), ]
    dms_meta2 <- dms_meta2[order(-dms_meta2$lat), ]
    ####make sure your rep has nums_n_site which is in df2!!!
    dms_meta2$site <- factor(dms_meta2$site, levels = unique(dms_meta2$site),ordered = TRUE)
    
    rep <- ddply(dms_meta2,
                 .(site),
                 summarise,
                 lat  = mean(lat),
                 long = mean(long))
    rep2 <- rep[order(-rep$lat), ]
    x <- length(unique(dms_meta2$site))
    bgcols = rev(matlab.like2(x))
    rep2$nums_n_site <- paste0(1:x, ": ", rep2$site)
    rep2$num_site <- 1:x
    write.table(rep2, paste0(LEA_dir,"site_number_fromLEA.csv"))
  }
  
  coord.pop11 <- merge(coord.pop,rep2[,c("site","nums_n_site","num_site")],by="site")
  
  #coord.pop <- coord.pop[order(-coord.pop$lat),]
  coord.pop11$lab <- coord.pop11$num_site
  
  points(coord.pop11$long,coord.pop11$lat,pch=21,lwd=0.2,bg="blue",cex=1)
  #text(x=150.8758, y=-34.40813+0.1,labels="2 to 29",pos=3,cex=0.7,font=4)
  
  coord.pop11 <- coord.pop11[complete.cases(coord.pop11), ]#remove latlongs that are NAs
  pointLabel(coord.pop11$long, coord.pop11$lat, labels = paste("  ", coord.pop11$lab, "  ", sep=""), cex=0.7, font=4,offset=0.25)
  
  if(exists("excluded_sample_names")) {
    excluded1pops <- dms$meta$analyses[,analysis][(dms$meta$sample_names%in%excluded_sample_names)]
    # cat("samples from 1indiv pop excluded in map","\n")
  }else {
    excluded1pops <- NULL
  }
  
  if(sum(is.na(newdms_meta$lat))>0){
    excludedNApops <- as.character(unique(newdms_meta$site[is.na(newdms_meta$lat)])) #pops excluded from map because latlongs are NA
    # cat("samples wihtout latlongs excluded in map","\n")
  } else {
    excludedNApops <- NULL
  }
  
  if(exists("excludedNApops") | exists("excluded1pops")){
    par(fig=c(0.1, 0.9, 0.1, 0.9),new = TRUE) #c(x1, x2, y1, y2)
    # par(fig=c(0,0.5,0,0.5),new = TRUE) #c(x1, x2, y1, y2)
    all_excluded <- c(excluded1pops,excludedNApops)
    rep2$site <- as.factor(rep2$site)
    all_ex_nums <-  rep2$nums_n_site[(rep2$site%in%all_excluded)]
    all_ex_nums2 <- paste0(all_ex_nums, "\n",collapse = "")
    all_ex_nums3 <- paste0("Not included:\n",all_ex_nums2)
    mtext( all_ex_nums3,cex = 0.7,side=1,adj=0)
    dev.off()
  } else {
    dev.off()
  }
  
  # cat(paste0("sampled distribution figure done!","\n"))
  
  ##################TABLE
  coord.pop2 <- cbind.data.frame(coord.pop11$lab,coord.pop11$site,coord.pop11$lat,coord.pop11$long)
  colnames(coord.pop2)<-c("","site","lat","long")
  coord.pop2 <- coord.pop2[order(-coord.pop2$lat),]
  tg = gridExtra::tableGrob(coord.pop2,rows = NULL)
  h = grid::convertHeight(sum(tg$heights), "in", TRUE)+0.5
  w = grid::convertWidth(sum(tg$widths), "in", TRUE)+0.5
  ggplot2::ggsave(paste0(LEA_dir,species," samp_distrib2.tiff"), tg, width=w, height=h)
  # cat(paste0("samp_distrib2 figure done!","\n"))
  
  ##################combine with sampl_dist
  setwd(outputloc)
  
  loc_1<- dir(path = "LEA", pattern = "samp_distrib2.tiff$", full.names = TRUE, recursive = TRUE) ###must set working directory to make readPNG work
  loc_2 <- dir(path = "LEA", pattern = "samp_distrib.tiff$", full.names = TRUE, recursive = TRUE) ###must set working directory to make readPNG work
  loc_3 <- dir(path = "LEA", pattern = "Entropy.tiff$", full.names = TRUE, recursive = TRUE) ###must set working directory to make readPNG work
  
  lloo <- c(loc_3)
  plots <- lapply(ll <- lloo,function(x){
    img <- as.raster(readTIFF(x)) ##change this to readTIFF if readPNG doesnt work
    rasterGrob(img, interpolate = FALSE)
  })
  #lay <- rbind(c(2,2,1,1),c(2,2,1,1),c(2,2,1,1),c(3,3,3,3))
  lay <- rbind(c(1,1))
  ggsave(paste0("LEA/",species, " distrib_entrp.tiff"),
         width=8, height=10, units = "in", device='tiff',dpi = 300,
         marrangeGrob(grobs = plots, layout_matrix = lay,top=textGrob(paste0(analysis," LEA outputs"), gp=gpar(fontsize=20,fontface=4))))
  
  # cat(paste0("combined map and site for more lea outputs!","\n"))
  
  
  ####combine all pies
  piedatalist=list()
  pie_ls_files <- c("LEA pies_bar k=2.tiff$",
                    "LEA pies_bar k=3.tiff$",
                    "LEA pies_bar k=4.tiff$",
                    "LEA pies_bar k=5.tiff$",
                    "LEA pies_bar k=6.tiff$",
                    "LEA pies_bar k=7.tiff$")
  
  for (f in 1:length(pie_ls_files)) {
    ll_filenamesleapies<- dir(path = "LEA", pattern = pie_ls_files[f], full.names = TRUE, recursive = TRUE) ###must set working directory to make readPNG work
    piedatalist[[f]] <- ll_filenamesleapies
  }
  
  pie_plots <- lapply(ll <- piedatalist,function(x){
    img <- as.raster(readTIFF(x)) ##change this to readTIFF if readPNG doesnt work
    rasterGrob(img, interpolate = FALSE)
  })
  
  
  lay <- rbind(c(1,2,3),c(4,5,6)) #for kvals from 2 to 7
  ggsave(paste0("LEA/",species, " LEA pies_bar_ALL.tiff"),
         width=16, height=14,dpi = 300, units = "in", device='tiff',
         marrangeGrob(grobs = pie_plots, layout_matrix = lay,top=NULL))
  # cat(paste0("combined all lea piebar plots!","\n"))
  
  
  ###combine LEA figures, pies for k equals 2,3,4,5, sampling distrib, entropy into a single figure
  datalist=list()
  
  ls_files <- c("distrib_entrp.tiff$",
                "LEA pies_bar_ALL.tiff$")
  
  for (g in 1:length(ls_files)) {
    ll_filenames <- dir(path = "LEA", pattern = ls_files[g], full.names = TRUE, recursive = TRUE) ###must set working directory to make readPNG work
    datalist[[g]] <- ll_filenames
  }
  locsLEA <- unlist(datalist)
  
  plots <- lapply(ll <- locsLEA,function(x){
    img <- as.raster(readTIFF(x)) ##change this to readTIFF if readPNG doesnt work
    rasterGrob(img, interpolate = FALSE)
  })
  
  #layoo <- rbind(c(1,2))
  layoo <- rbind(c(1,2,2,2),c(1,2,2,2),c(1,2,2,2),c(1,2,2,2))
  
  ggsave(paste0(LEA_dir,species, " lea k2345loc bar test.tiff"),
         width=12, height=8.3,dpi = 600, units = "in", device='tiff',
         marrangeGrob(grobs = plots, layout_matrix = layoo,top=NULL))
  cat(paste0("all LEA plots done!","\n"))
  
  # locsLEA
  #"LEA/PersHirs2 distrib_entrp.tiff"    "LEA/PersHirs2 LEA pies_bar_ALL.tiff"
  # pielocsLEA
  #[1] "temp/PersHirs2 LEA pies_bar k=2.tiff" "temp/PersHirs2 LEA pies_bar k=3.tiff"
  # [3] "temp/PersHirs2 LEA pies_bar k=4.tiff" "temp/PersHirs2 LEA pies_bar k=5.tiff"
  # [5] "temp/PersHirs2 LEA pies_bar k=6.tiff" "temp/PersHirs2 LEA pies_bar k=7.tiff"
  # [7] "temp/PersHirs2 LEA pies_bar k=8.tiff"
  # lloo
  # [1] "LEA/PersHirs2 samp_distrib2.tiff" "LEA/PersHirs2 samp_distrib.tiff"
  
  
  
  # file.remove(locsLEA)
  # pielocsLEA <- unlist(piedatalist)
  # file.remove(pielocsLEA)
  # file.remove(lloo)
  
  setwd(maindir)
  
  
} else {
  print("remove LEA folder under popgen folder first!!!")
  #or if you have already run LEA, and just want to replot everything?
}



##### PCA, latitude PCA and Centroid PCA #####
m1 <- read.meta.data(d3, RandRbase, species, dataset, fields = as.numeric(ncol(metafile)-4))
temp<-as.data.frame(m1$analyses)

#set a sample threshold for each pop:
samplethreshold <-1# chnage this to what number of samples you want toi have uniform across all pops

## If you want to just remove any sites that are below the threshold (but keep individuals that are above a threshold, then unhash the below section)
for (x in unique(na.omit(temp[, analysis]))) {
  pop_samples <- which(temp[, analysis] == x)
  if (length(pop_samples) < samplethreshold) {
    temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
    print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
  } else {}
}
m1$analyses<-temp
dm        <- dart.meta.data.merge(d3, m1)
fields    <- c(analysis)
dms     <- data.by.meta.fields(dm, fields, RandRbase, species, dataset, object=analysis); print(dms$meta$analyses[,analysis])

# PCAcols <- matlab.like2(length(unique(dms$meta$analyses[, analysis])))


PCAcols <- rev(matlab.like2(length(unique(dms$meta$analyses[, analysis]))))
# PCAcols <- c("red", "blue")
# matlab.like2(n)
# blue2red(n)
# blue2green2red(n)


nd_gl <- dart2gl(dms, RandRbase, species, dataset) # converts the cleaned data to genlight format
nd_pca <- glPca(nd_gl, nf = 5, parallel = FALSE) # nf indicates the number of principal components to be retained
scatter(nd_pca) # a quick check of PCA

PCAdirectory <- paste0(outputloc, "/PCA/")
if (!dir.exists(PCAdirectory)) {
  dir.create(PCAdirectory)
}

sample_PC1_PC2 <- cbind(dms$meta$sample_names,
                        dms$meta$analyses[, analysis],
                        dms$meta$lat, dms$meta$long,
                        nd_pca$scores[, 1], # PC1
                        nd_pca$scores[, 2], # PC2
                        nd_pca$scores[, 3]) # PC3

colnames(sample_PC1_PC2) <- c("names", "site", "lat", "long", "PC1", "PC2", "PC3")
write.table(sample_PC1_PC2,
            paste0(outputloc, "/PCA/", species, " sample_PC1_PC2_PC3.csv"),
            sep = ",", col.names = TRUE, row.names = FALSE)

names <- rownames(sample_PC1_PC2)
rownames(sample_PC1_PC2) <- NULL
data <- cbind(names, sample_PC1_PC2)

df <- data.frame(data, stringsAsFactors = FALSE)
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)
df$PC3 <- as.numeric(df$PC3)
df$lat <- as.numeric(df$lat)
df$long <- as.numeric(df$long)

df2new <- df[order(-df$lat),]
df2new$site <- factor(df2new$site, levels = unique(df2new$site), ordered = TRUE)

rep <- ddply(df2new,
             .(site),
             summarise,
             lat = mean(lat),
             long = mean(long),
             PC1 = mean(PC1),
             PC2 = mean(PC2),
             PC3 = mean(PC3))

rep2 <- rep[order(-rep$lat),]
x <- length(unique(df2new$site))
bgcols <- rev(matlab.like2(x))
rep2$nums_n_site <- paste0(1:x, ": ", rep2$site)
rep2$num_site <- 1:x
rep222 <- rep2
rep222$nums_n_site <- factor(rep222$nums_n_site, levels = unique(rep222$nums_n_site))
df3 <- merge(df, rep222[, c(1, 7, 8)], by = "site")
f <- nd_pca$eig[nd_pca$eig > sum(nd_pca$eig / length(nd_pca$eig))]
e <- round(f * 100 / sum(nd_pca$eig), 1)

percvar <- cbind(e[1], e[2], e[3])
colnames(percvar) <- c("PC1", "PC2", "PC3")
write.table(percvar,
            paste0(outputloc, "/PCA/", species, " percent_variance.csv"),
            sep = ",", col.names = TRUE, row.names = FALSE)

options(ggrepel.max.overlaps = Inf)



PC1_PC2_plot <- ggplot(data = df3, aes(x = PC1, y = PC2, z = PC3, col = nums_n_site, shape = nums_n_site)) +
  geom_point(alpha = 1, size = 3, stroke=1.5) +
  theme_bw() +
  theme(axis.line = element_line(colour = "grey50"), legend.text = element_text(size = 10),
        axis.title = element_text(size = 30, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold.italic", size = 10)) +
  ylab(paste0("PC2", " (", e[2], "%)")) +
  xlab(paste0("PC1", " (", e[1], "%)")) +
  scale_shape_manual(labels = rep222$nums_n_site,
                     values = shapeassign)+
  scale_color_manual(labels = rep222$nums_n_site,
                     values = PCAcols)+
  guides(col = guide_legend(ncol = 6, byrow = F))+
  geom_text_repel(aes(label = num_site, fontface="bold"), col = PCAcols, size = 3, data = rep222, segment.color = NA, min.segment.length = unit(0.1, "lines"))
# geom_text_repel(aes(label = num_site, fontface="bold"), col = PCAcols, size = 3, data = df3, segment.color = NA, min.segment.length = unit(0.1, "lines"))
PC1_PC2_plot
ggsave(paste0(outputloc, "/PCA/", species, "_PC1-PC2.tiff"), width = 16, height = 10, dpi = 300, units = "in")


PC1_PC2_plot <- ggplot(data = df3, aes(x = PC1, y = PC2, z = PC3, col = nums_n_site, shape = nums_n_site)) +
  geom_point(alpha = 1, size = 3, stroke=1.5) +
  theme_bw() +
  theme(axis.line = element_line(colour = "grey50"), legend.text = element_text(size = 10),
        axis.title = element_text(size = 30, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold.italic", size = 10)) +
  guides(col = guide_legend(ncol = 6)) +
  ylab(paste0("PC2", " (", e[2], "%)")) +
  xlab(paste0("PC1", " (", e[1], "%)")) +
  scale_shape_manual(labels = rep222$nums_n_site,
                     values = shapeassign)+
  scale_color_manual(labels = rep222$nums_n_site,
                     values = PCAcols)+
  theme(legend.position="none")
# geom_text_repel(aes(label = num_site, fontface="bold"), col = PCAcols, size = 3, data = rep222, segment.color = NA, min.segment.length = unit(0.1, "lines"))
# geom_text_repel(aes(label = num_site, fontface="bold"), col = PCAcols, size = 3, data = df3, segment.color = NA, min.segment.length = unit(0.1, "lines"))
PC1_PC2_plot
ggsave(paste0(outputloc, "/PCA/", species, "_PC1-PC2-NoLegend.tiff"), width = 10, height = 10, dpi = 300, units = "in")


PC1_PC3_plot <- ggplot(data = df3, aes(x = PC1, y = PC3, col = nums_n_site, shape = nums_n_site)) +
  geom_point(alpha = 1, size = 3, stroke=1.5) +
  geom_text_repel(aes(label = num_site), col = PCAcols, size = 2, data = rep222, segment.color = NA) +
  theme_bw() +
  theme(axis.line = element_line(colour = "grey50"), legend.text = element_text(size = 10),
        axis.title = element_text(size = 30, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold.italic", size = 10)) +
  guides(col = guide_legend(ncol = 6)) +
  ylab(paste0("PC3", " (", e[3], "%)")) +
  xlab(paste0("PC1", " (", e[1], "%)")) +
  scale_color_manual(labels = rep222$nums_n_site,
                     values = PCAcols)+
  scale_shape_manual(labels = rep222$nums_n_site,
                     values = shapeassign)
PC1_PC3_plot
ggsave(paste0(outputloc, "/PCA/", species, "_PC1-PC3.tiff"), width = 11, height = 10, dpi = 300, units = "in")

PC2_PC3_plot <- ggplot(data=df3, aes(x=PC2, y=PC3,col=nums_n_site, shape = nums_n_site))+
  geom_point(alpha = 1, size = 3, stroke=1.5) +
  geom_text_repel(aes(label = num_site),col=PCAcols,size=2,data =rep222,segment.color = NA)+
  theme_bw() +
  theme(axis.line = element_line(colour = "grey50"), legend.text=element_text(size=10),
        axis.title=element_text(size=30, face= "bold"),
        axis.text=element_text(size=20, face= "bold"),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold.italic", size=10))+
  guides(col = guide_legend(ncol = 6)) + #how many columns in the legend
  ylab(paste0("PC3"," (",e[3],"%)"))+
  xlab(paste0("PC2"," (",e[2],"%)"))+
  scale_color_manual(labels = rep222$nums_n_site,
                     values = PCAcols)+
  scale_shape_manual(labels = rep222$nums_n_site,
                     values = shapeassign)
PC2_PC3_plot
ggsave(paste0(outputloc, "/PCA/", species,"_PC2-PC3.tiff"), width = 11, height = 10, dpi = 300, units = "in")


#now plot all 3 PCAs together
options(ggrepel.max.overlaps = Inf)
PC1_PC2_plot <- ggplot(data=df3, aes(x=PC1, y=PC2,col=nums_n_site, shape = nums_n_site))+
  geom_point(alpha = 1, size = 2, stroke=1) +
  # geom_text_repel(aes(label = num_site),col=PCAcols,size=2,data =rep222,segment.color = NA)+
  theme_bw() +
  theme(axis.line = element_line(colour = "grey50"), legend.text=element_text(size=10),
        axis.title=element_text(size=15, face= "bold"),
        axis.text=element_text(size=10, face= "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold.italic", size=10))+
  guides(col = guide_legend(ncol = 6)) + #how many columns in the legend
  ylab(paste0("PC2"," (",e[2],"%)"))+
  xlab(paste0("PC1"," (",e[1],"%)"))+
  scale_color_manual(labels = rep222$nums_n_site,
                     values = PCAcols)+
  scale_shape_manual(labels = rep222$nums_n_site,
                     values = shapeassign)
PC1_PC3_plot <- ggplot(data=df3, aes(x=PC1, y=PC3,col=nums_n_site, shape = nums_n_site))+
  geom_point(alpha = 1, size = 2, stroke=1) +
  # geom_text_repel(aes(label = num_site),col=PCAcols,size=2,data =rep222,segment.color = NA)+
  theme_bw() +
  theme(axis.line = element_line(colour = "grey50"), legend.text=element_text(size=10),
        axis.title=element_text(size=15, face= "bold"),
        axis.text=element_text(size=10, face= "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold.italic", size=10))+
  guides(col = guide_legend(ncol = 6)) + #how many columns in the legend
  ylab(paste0("PC3"," (",e[3],"%)"))+
  xlab(paste0("PC1"," (",e[1],"%)"))+
  scale_color_manual(labels = rep222$nums_n_site,
                     values = PCAcols)+
  scale_shape_manual(labels = rep222$nums_n_site,
                     values = shapeassign)
PC2_PC3_plot <- ggplot(data=df3, aes(x=PC2, y=PC3,col=nums_n_site, shape = nums_n_site))+
  geom_point(alpha = 1, size = 2, stroke=1) +
  # geom_text_repel(aes(label = num_site),col=PCAcols,size=2,data =rep222,segment.color = NA)+
  theme_bw() +
  theme(axis.line = element_line(colour = "grey50"), legend.text=element_text(size=10),
        axis.title=element_text(size=15, face= "bold"),
        axis.text=element_text(size=10, face= "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold.italic", size=10))+
  guides(col = guide_legend(ncol = 6)) + #how many columns in the legend
  ylab(paste0("PC3"," (",e[3],"%)"))+
  xlab(paste0("PC2"," (",e[2],"%)"))+
  scale_color_manual(labels = rep222$nums_n_site,
                     values = PCAcols)+
  scale_shape_manual(labels = rep222$nums_n_site,
                     values = shapeassign)
legend_a <- get_legend(PC_plot1 + theme(legend.position="bottom",
                                        legend.text=element_text(size=8)))
prow1 <- plot_grid(PC1_PC2_plot,PC1_PC3_plot,PC2_PC3_plot, ncol=3)
p1 <- plot_grid( prow1, legend_a, ncol = 1, rel_heights = c(1, .3))
p1
ggsave(paste0(outputloc, "/PCA/", species,"_Combined PCA.tiff"), width = 16, height = 8, dpi = 300, units = "in")






# 
# PC1_PC2_plot <- ggplot(data = df3, aes(x = PC1, y = PC2, z = PC3, col = nums_n_site)) +
#   geom_point(alpha = 0.8, size = 5) +
#   geom_text_repel(aes(label = num_site), col = PCAcols, size = 2, data = rep222, segment.color = NA) +
#   theme_minimal() +
#   theme(axis.line = element_line(colour = "grey50"), legend.text = element_text(size = 10),
#         axis.title = element_text(size = 30, face = "bold"),
#         axis.text = element_text(size = 20, face = "bold"),
#         legend.title = element_blank(),
#         legend.position = "bottom",
#         plot.title = element_text(face = "bold.italic", size = 10)) +
#   guides(col = guide_legend(ncol = 6)) +
#   ylab(paste0("PC2", " (", e[2], "%)")) +
#   xlab(paste0("PC1", " (", e[1], "%)")) +
#   scale_color_manual(labels = rep222$nums_n_site,
#                      values = PCAcols)
# PC1_PC2_plot
# ggsave(paste0(outputloc, "/PCA/", species, "_PC1-.tiff"), width = 14, height = 10, dpi = 300, units = "in")
# 
# PC1_PC3_plot <- ggplot(data = df3, aes(x = PC1, y = PC3, col = nums_n_site)) +
#   geom_point(alpha = 0.8, size = 5) +
#   geom_text_repel(aes(label = num_site), col = PCAcols, size = 2, data = rep222, segment.color = NA) +
#   theme_minimal() +
#   theme(axis.line = element_line(colour = "grey50"), legend.text = element_text(size = 10),
#         axis.title = element_text(size = 30, face = "bold"),
#         axis.text = element_text(size = 20, face = "bold"),
#         legend.title = element_blank(),
#         legend.position = "bottom",
#         plot.title = element_text(face = "bold.italic", size = 10)) +
#   guides(col = guide_legend(ncol = 6)) +
#   ylab(paste0("PC3", " (", e[3], "%)")) +
#   xlab(paste0("PC1", " (", e[1], "%)")) +
#   scale_color_manual(labels = rep222$nums_n_site,
#                      values = PCAcols)
# PC1_PC3_plot
# ggsave(paste0(outputloc, "/PCA/", species, "_PC1-PC3.tiff"), width = 14, height = 10, dpi = 300, units = "in")
# 
# 
# PC2_PC3_plot <- ggplot(data=df3, aes(x=PC2, y=PC3,col=nums_n_site))+
#   geom_point(alpha=0.8, size=5) + #alpha is for transparency
#   geom_text_repel(aes(label = num_site),col=PCAcols,size=2,data =rep222,segment.color = NA)+
#   theme_minimal()+
#   theme(axis.line = element_line(colour = "grey50"), legend.text=element_text(size=10),
#         axis.title=element_text(size=30, face= "bold"),
#         axis.text=element_text(size=20, face= "bold"),
#         legend.title = element_blank(),
#         legend.position = "bottom",
#         plot.title = element_text(face = "bold.italic", size=10))+
#   guides(col = guide_legend(ncol = 6)) + #how many columns in the legend
#   ylab(paste0("PC3"," (",e[3],"%)"))+
#   xlab(paste0("PC2"," (",e[2],"%)"))+
#   scale_color_manual(labels=rep222$nums_n_site,
#                      values = PCAcols)
# PC2_PC3_plot
# ggsave(paste0(outputloc, "/PCA/", species,"_PC2-PC3.tiff"), width = 14, height = 10, dpi = 300, units = "in")
# 
# #now plot all 3 PCAs together
# options(ggrepel.max.overlaps = Inf)
# PC1_PC2_plot <- ggplot(data=df3, aes(x=PC1, y=PC2,col=nums_n_site))+
#   geom_point(alpha=0.8, size=4) + #alpha is for transparency
#   geom_text_repel(aes(label = num_site),col=PCAcols,size=2,data =rep222,segment.color = NA)+
#   theme_minimal()+
#   theme(axis.line = element_line(colour = "grey50"), legend.text=element_text(size=10),
#         axis.title=element_text(size=15, face= "bold"),
#         axis.text=element_text(size=10, face= "bold"),
#         legend.title = element_blank(),
#         legend.position = "none",
#         plot.title = element_text(face = "bold.italic", size=10))+
#   guides(col = guide_legend(ncol = 6)) + #how many columns in the legend
#   ylab(paste0("PC2"," (",e[2],"%)"))+
#   xlab(paste0("PC1"," (",e[1],"%)"))+
#   scale_color_manual(labels=rep222$nums_n_site,
#                      values = PCAcols)
# 
# PC1_PC3_plot <- ggplot(data=df3, aes(x=PC1, y=PC3,col=nums_n_site))+
#   geom_point(alpha=0.8, size=4) + #alpha is for transparency
#   geom_text_repel(aes(label = num_site),col=PCAcols,size=2,data =rep222,segment.color = NA)+
#   theme_minimal()+
#   theme(axis.line = element_line(colour = "grey50"), legend.text=element_text(size=10),
#         axis.title=element_text(size=15, face= "bold"),
#         axis.text=element_text(size=10, face= "bold"),
#         legend.title = element_blank(),
#         legend.position = "none",
#         plot.title = element_text(face = "bold.italic", size=10))+
#   guides(col = guide_legend(ncol = 6)) + #how many columns in the legend
#   ylab(paste0("PC3"," (",e[3],"%)"))+
#   xlab(paste0("PC1"," (",e[1],"%)"))+
#   scale_color_manual(labels=rep222$nums_n_site,
#                      values = PCAcols)
# 
# PC2_PC3_plot <- ggplot(data=df3, aes(x=PC2, y=PC3,col=nums_n_site))+
#   geom_point(alpha=0.8, size=4) + #alpha is for transparency
#   geom_text_repel(aes(label = num_site),col=PCAcols,size=2,data =rep222,segment.color = NA)+
#   theme_minimal()+
#   theme(axis.line = element_line(colour = "grey50"), legend.text=element_text(size=10),
#         axis.title=element_text(size=15, face= "bold"),
#         axis.text=element_text(size=10, face= "bold"),
#         legend.title = element_blank(),
#         legend.position = "none",
#         plot.title = element_text(face = "bold.italic", size=10))+
#   guides(col = guide_legend(ncol = 6)) + #how many columns in the legend
#   ylab(paste0("PC3"," (",e[3],"%)"))+
#   xlab(paste0("PC2"," (",e[2],"%)"))+
#   scale_color_manual(labels=rep222$nums_n_site,
#                      values = PCAcols)
# 
# theme_set(theme_gray())
# PC_plot1 <- ggplot(data=df3, aes(x=PC2, y=PC3,col=nums_n_site))+
#   geom_point(alpha=0.8, size=4) + #alpha is for transparency
#   geom_text_repel(aes(label = num_site),col=PCAcols,size=2,data =rep222,segment.color = NA)+
#   theme_minimal()+
#   theme(axis.line = element_line(colour = "grey50"), legend.text=element_text(size=10),
#         axis.title=element_text(size=15, face= "bold"),
#         axis.text=element_text(size=10, face= "bold"),
#         legend.title = element_blank(),
#         legend.position = "bottom",
#         plot.title = element_text(face = "bold.italic", size=10))+
#   guides(col = guide_legend(ncol = 6)) + #how many columns in the legend
#   ylab(paste0("PC3"," (",e[3],"%)"))+
#   xlab(paste0("PC2"," (",e[2],"%)"))+
#   scale_color_manual(labels=rep222$nums_n_site,
#                      #values = c("red", "green", "blue", "orange", "purple"))      #matlab.like2(length(unique(df2new$site))))
#                      values = PCAcols)
# 
# legend_a <- get_legend(PC_plot1 + theme(legend.position="bottom",
#                                         legend.text=element_text(size=8)))
# prow1 <- plot_grid(PC1_PC2_plot,PC1_PC3_plot,PC2_PC3_plot, ncol=3)
# p1 <- plot_grid( prow1, legend_a, ncol = 1, rel_heights = c(1, .3))
# p1
# 
# ggsave(paste0(outputloc, "/PCA/", species,"_Combined PCA.tiff"), width = 15, height = 8, dpi = 300, units = "in")
# 
# 
# 















## Latitude PCA 

PC1_PC2_plot <- ggplot(df2new, aes(PC1, PC2, col=lat))+ geom_point(alpha=0.5) +
  theme(legend.position="none",aspect.ratio = 1,
        plot.title = element_text(face = "bold.italic", size=9),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7))+
  ylab(paste0("PCA2"," (",e[2],"%)"))+ xlab(paste0("PCA1"," (",e[1],"%)"))+
  ggtitle(paste0(species, "\n"," PC1 vs PC2")) + scale_color_gradient(low="blue", high="red")

###PC1 vs PC3
PC1_PC3_plot <- ggplot(df2new, aes(PC1, PC3, col=lat))+ geom_point(alpha=0.5) +
  theme(legend.position="none",aspect.ratio = 1,
        plot.title = element_text(face = "bold.italic", size=9),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7))+
  ylab(paste0("PCA3"," (",e[3],"%)"))+ xlab(paste0("PCA1"," (",e[1],"%)"))+
  ggtitle(paste0("\n"," PC1 vs PC3")) + scale_color_gradient(low="blue", high="red")

###PC2 vs PC3
PC2_PC3_plot <- ggplot(df2new, aes(PC2, PC3, col=lat))+ geom_point(alpha=0.5) +
  ylab(paste0("PCA3"," (",e[3],"%)"))+ xlab(paste0("PCA2"," (",e[2],"%)"))+
  theme(legend.position="none",aspect.ratio = 1,
        plot.title = element_text(face = "bold.italic", size=9),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7)) +
  scale_color_gradient(low="blue", high="red")+
  ggtitle(paste0("\n"," PC2 vs PC3"))

PC_plot <- ggplot(df2new, aes(PC2, PC3, col=lat))+ geom_point() +
  ylab(paste0("PCA3"," (",e[3],"%)"))+ xlab(paste0("PCA2"," (",e[2],"%)"))+
  scale_color_gradient(low="blue", high="red")+
  ggtitle(paste0(" PC2 vs PC3"))

theme_set(theme_gray())

legend_b <- get_legend(PC_plot + theme(legend.position="bottom",
                                       legend.text=element_text(size=4))+ labs(col = "Latitude"))

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
prow <- plot_grid(PC1_PC2_plot,PC1_PC3_plot,PC2_PC3_plot, ncol=3)
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .3))
p

ggsave(paste0(outputloc, "/PCA/", species, " Latitude PCA.tiff"), width = 12, height = 6, dpi = 300, units = "in")


#calculate distance from PCA centroid for PC1 and PC2
meta2new <- df2new[order(-df2new$lat),] # this dataset is used in generating a summary of pops and their average latlongs for making maps/tables
meta2new$site <- factor(meta2new$site, levels = unique(meta2new$site),ordered = TRUE)
reg_names_xna <- unique(meta2new$reg)
num_reg <- 1

pcaData2 <- read.csv(paste0(outputloc,"/PCA/",species," sample_PC1_PC2_PC3.csv"),header=T)
pcaData <- pcaData2[order(-pcaData2$lat),]
pcaData$site <- factor(pcaData$site, levels = unique(pcaData$site),ordered = TRUE)

#d2 <- merge(pcaData, meta2new[,c("names")], by.x="X", by.y="names")
d <- df2new[order(-df2new$lat),]
d$site <- factor(d$site, levels = unique(d$site),ordered = TRUE)
pcaVAR <- read.csv(paste0(outputloc,"/PCA/",species," percent_variance.csv"),header=T)

# Compute median coordinates of PC1 and PC2
medCoords <- unlist(c(spatial.median(as.numeric(pcaData$PC1)), spatial.median(as.numeric(pcaData$PC2)), spatial.median(as.numeric(pcaData$PC3))))

###plot out PC123 with median included
pca_plot12 <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(shape = 16, colour = "darkgrey", size = 0.5)+ ggtitle(species)+
  labs(x=paste0("PC1"," (",pcaVAR[1],"%)"), y=paste0("PC2"," (",pcaVAR[2],"%)"))+theme_classic()+theme(plot.title = element_text(size=5))
pca_plot13 <- ggplot(pcaData, aes(x = PC1, y = PC3)) +
  geom_point(shape = 16, colour = "darkgrey", size = 0.5)+
  labs(x=paste0("PC1"," (",pcaVAR[1],"%)"),y=paste0("PC3"," (",pcaVAR[3],"%)"))+theme_classic()
pca_plot23 <- ggplot(pcaData, aes(x = PC2, y = PC3)) +
  geom_point(shape = 16, colour = "darkgrey", size = 0.5)+
  labs(x=paste0("PC2"," (",pcaVAR[2],"%)"), y=paste0("PC3"," (",pcaVAR[3],"%)"))+theme_classic()

for (i in 1:num_reg){
  pca_plot12 <- pca_plot12+
    geom_point(aes(x = medCoords[1], y = medCoords[2]), shape = "+", col = "red", size = 5)
  pca_plot13 <- pca_plot13+
    geom_point(aes(x = medCoords[1], y = medCoords[2]), shape = "+", col = "red", size = 5)
  pca_plot23 <- pca_plot23+
    geom_point(aes(x = medCoords[1], y = medCoords[2]), shape = "+", col = "red", size = 5)
  
}

pca_plot12 #check, the red dot will be the median centroid, each point represent an individual on a multivariate space,
med123 <- plot_grid(pca_plot12,pca_plot13,pca_plot23, ncol=3)
ggsave(med123, file=paste0(outputloc,"/PCA/",species,"_median_PC123.tiff"),
       width = 8, height = 3, dpi = 300, units = "in", device='tiff')

for (i in 1:length(df2new$names)){
  df2new$DistPointPC12[i] <- crossdist(df2new$PC1[i], df2new$PC2[i],medCoords[1],medCoords[2])
  
}


#this bit checks for zeros in lats and longs and cut out those samples
row_sub = apply(data.frame(dms$meta$long), 1, function(row) all(row !=0 ))
newdmslong <- data.frame(dms$meta$long)[row_sub,]
row_sub = apply(data.frame(dms$meta$lat), 1, function(row) all(row !=0 ))
newdmslat <- data.frame(dms$meta$lat)[row_sub,]

#this bit looks for range of latlong values to zoom in on and calulcates range for zoom for ggmap
divxlims <- c(min(na.omit(newdmslong))-1,max(na.omit(newdmslong))+1)
divylims <- c(min(na.omit(newdmslat))-1,max(na.omit(newdmslat))+1)
divxlimsrange <- abs(c(min(na.omit(newdmslong))-max(na.omit(newdmslong))))
divylimsrange <- abs(c(min(na.omit(newdmslat))-max(na.omit(newdmslat))))


coordinate_cities <- read.csv(paste0(maindir,"CityMapCordinates.csv"))
sf_oz <- ozmap_data("states")

#order and map expected het as a coloured range
df2new <- df2new[order(df2new$lat),]

ggplot(sf_oz, aes(fill="white")) + geom_sf() +
  scale_fill_manual(values = "white", guide = "none") +
  theme_void() + geom_sf() +
  xlim(divxlims) + ylim(divylims)+
  geom_point(data = df2new, mapping = aes(x = long, y = lat, color = DistPointPC12), size=5, alpha=0.9) +
  scale_color_gradient(low="yellow", high="red", name="Dist From Centroid") +
  geom_point(data=coordinate_cities, mapping = aes(x = long, y = lat), size=1)+
  geom_text(data=coordinate_cities, mapping = aes(x = long, y = lat, label=City),hjust=0,nudge_y=0.03, size=3)


ggsave(paste0(outputloc, "/PCA/", species, "Dist_from_PCA12_Centroid.tiff"), width = 10, height = 12, dpi = 300, units = "in")



##### FST, IBD and distance map with at least 1 sample per pop #####

m1 <- read.meta.data(d3, RandRbase, species, dataset, fields = as.numeric(ncol(metafile)-4))
temp<-as.data.frame(m1$analyses)

#set a sample threshold for each pop:
samplethreshold <-2# chnage this to what number of samples you want toi have uniform across all pops



# #### Below lines remove sites with less than threshold number of samples, and subsample sites that are above the threshold
# for (x in unique(na.omit(temp[, analysis]))) {
#   pop_samples <- which(temp[, analysis] == x)
#   if (length(pop_samples) < samplethreshold) {
#     temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
#     print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
#   }
#   if (length(pop_samples) == samplethreshold) {
#     print(paste0("Keeping ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
#   }
#   if (length(pop_samples) > samplethreshold) {
#     subsamplepops <- pop_samples[sample(1:length(pop_samples), size=samplethreshold, replace = F)]
#     for (b in 1:length(pop_samples)){
#       if (pop_samples[b] %in% subsamplepops == FALSE){
#         temp[, analysis][pop_samples[b]] <- NA
#       }
#       print(paste0("Reducing ",x, " to a sample size of ", samplethreshold," samples"))
#     }
#   } else {}
# }
# 
# 


## If you want to just remove any sites that are below the threshold (but keep individuals that are above a threshold, then unhash the below section)
for (x in unique(na.omit(temp[, analysis]))) {
  pop_samples <- which(temp[, analysis] == x)
  if (length(pop_samples) < samplethreshold) {
    temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
    print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
  } else {}
}

m1$analyses<-temp
dm        <- dart.meta.data.merge(d3, m1)
fields    <- c(analysis)
dms     <- data.by.meta.fields(dm, fields, RandRbase, species, dataset, object=analysis); print(dms$meta$analyses[,analysis])


# Generate pairwise fst for populations
# Determine if IBD is significant - Mantel test
pFst <- population.pw.Fst(dms, dms$meta$analyses[, analysis], RandRbase, species, dataset)
pS <- population.pw.spatial.dist(dms, dms$meta$analyses[, analysis])

fst_dir <- paste0(outputloc, "/fst_atleast1Samp/")
if (!dir.exists(fst_dir)) {
  dir.create(fst_dir)
}

# Plot of fst vs geographic distance
tiff(paste0(fst_dir, species, "_fstplot_atleast1Samps.tiff"),
     units = "in", width = 15, height = 8, res = 300)
par(mfrow = c(1, 2))
diag(pS$S) <- NA
diag(pFst$Fst) <- NA

Fst_sig <- cbind(melt(pS$S), unlist(as.list(pFst$Fst)))

colnames(Fst_sig)[3] <- "Geo_dist"
colnames(Fst_sig)[4] <- "Fst"
Fst_sig$Geo_dist2 <- Fst_sig$Geo_dist / 1000

plot(Fst_sig$Geo_dist2, Fst_sig$Fst, xlab = "distance (km)", ylab = "Fst", cex = 1,
     font = 4, cex.main = 1, cex.lab=1.6)
title(main = paste0(species, " ", analysis, " pairwise fst plots"), adj = 0.001, font.main = 4)
man <- mantel(xdis = pS$S, ydis = pFst$Fst, permutations = 999, na.rm = TRUE)
legend("bottomright", bty = "n", cex = 2, text.col = "blue", element_text(face = "bold"),
       legend = paste("Mantel r = ",
                      format(man$statistic, digits = 4),
                      " P =", format(man$signif)))

# Plot of linearised fst vs geographic distance
Fst_sig <- cbind(melt(pS$S), unlist(as.list(pFst$Fst)))
colnames(Fst_sig)[3] <- "Geo_dist"
colnames(Fst_sig)[4] <- "Fst"
Fst_sig$lin_fst <- Fst_sig$Fst / (1 - Fst_sig$Fst)
Fst_sig$Geo_dist2 <- Fst_sig$Geo_dist / 1000
Fst_sig$log10Geo_dist <- log10(Fst_sig$Geo_dist2)
write.table(Fst_sig,
            paste0(fst_dir, species, " fst_geo data.csv"))
plot(Fst_sig$log10Geo_dist, Fst_sig$lin_fst, xlab = "log10(distance)",
     ylab = "Linearised Fst", font = 4, cex.main = 1, cex.lab=1.6)
man2 <- mantel(xdis = pS$S, ydis = (pFst$Fst) / (1 - pFst$Fst), permutations = 999, na.rm = TRUE)
legend("topleft", bty = "n", cex = 2, text.col = "blue", element_text(face = "bold"),
       legend = paste("Mantel r = ",
                      format(man2$statistic, digits = 4),
                      " P =", format(man2$signif)))
dev.off()

cat(paste0("pairwise fst and linearised fst drawn!", "\n"))

# Heatmap of just geographic distance or fst
par(mfrow = c(2, 1), oma = c(0, 0, 1, 0))
geo_d <- pS$S
geo_d[upper.tri(geo_d)] <- NA
rownames(geo_d) <- colnames(pS$S)

dimnames <- list(var1 = colnames(pS$S), var2 = colnames(pS$S))
mat <- matrix(geo_d, ncol = length(colnames(geo_d)), nrow = length(colnames(geo_d)), dimnames = dimnames)
df <- as.data.frame(as.table(mat))

p1 <- ggplot(df, aes(var1, var2)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  xlab("") + ylab("") +
  theme(legend.position = "none", axis.text = element_text(size = 5), plot.title = element_text(face = "bold.italic")) +
  ggtitle(paste0(species, " ", analysis, " heatmaps", "\n", "heatmap of pairwise geographic distance"))

genetic_d <- pFst$Fst
genetic_d[upper.tri(genetic_d)] <- NA
rownames(genetic_d) <- colnames(pFst$Fst)

dimnames2 <- list(var1 = colnames(pFst$Fst), var2 = colnames(pFst$Fst))
mat2 <- matrix(genetic_d, ncol = length(colnames(geo_d)), nrow = length(colnames(geo_d)), dimnames = dimnames)
df2 <- as.data.frame(as.table(mat2))
df3 <- df2[complete.cases(df2$Freq),]

p2 <- ggplot(df3, aes(var1, var2)) +
  geom_tile(aes(fill = Freq), colour = "white", na.rm = TRUE) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = round(Freq, 3)), size = 2, df3) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  xlab("") + ylab("") +
  theme(legend.position = "none", axis.text = element_text(size = 5),
        plot.title = element_text(face = "bold.italic")) + ggtitle("heatmap of pairwise Fst with atleast 1 sample per pop")

plot_grid(p1, p2, ncol = 1)

ggsave(paste0(fst_dir, species, "_fstheatmaps_atleast1Samps.tiff"),
       width = 13.3, height = 7.5, dpi = 300, units = "in", device = 'tiff')

cat(paste0("geographic and pairwise fst heatmap drawn!", "\n"))

####save matrix
genetic_d <-pFst$Fst
write.table(genetic_d, paste0(fst_dir,species," fst matrix.csv"),sep = ",")
geo_d <-pS$S
geo_d[upper.tri(geo_d)] <-geo_d[lower.tri(geo_d)]
new <- matrix(NA, nrow = dim(geo_d)[1], ncol = dim(geo_d)[2])
new[upper.tri(new)] <- geo_d[upper.tri(geo_d)]
new[lower.tri(new)] <- genetic_d[lower.tri(genetic_d)]
colnames(new) <- colnames(geo_d)
rownames(new) <- rownames(geo_d)

write.table(new, paste0(fst_dir,species," heatmap matrix.csv"),sep = ",")



#### FST distance script 

fst_dms <- dms

# calculate FST and geodist
pop_gds_file <- dart2gds(fst_dms, RandRbase, species, dataset)
pop_pFst <- population.pw.Fst(fst_dms, fst_dms$meta$analyses[,analysis], RandRbase, species, dataset, maf_val = 0.05, miss_val = 0.2)
pop_pS <- population.pw.spatial.dist(fst_dms, fst_dms$meta$analyses[,analysis])

# meta_agg <- m2 %>%
#   group_by(pop_large, genetic_group) %>%
#   summarize(lat = mean(lat, na.rm=TRUE),
#             long = mean(long,na.rm=TRUE),
#             .groups = 'drop')



dms_df <- data.frame(cbind(fst_dms$meta$analyses[,analysis],fst_dms$meta$lat, fst_dms$meta$long))
colnames(dms_df) <- c("ansites", "lat", "long")

dms_df$lat <- as.numeric(dms_df$lat)
dms_df$long <- as.numeric(dms_df$long)


meta_agg <- dms_df %>% 
  group_by(dms_df$ansites) %>% 
  dplyr::summarize(lat = mean(lat, na.rm=TRUE),
                   long = mean(long, na.rm=TRUE), 
                   .groups = 'drop')


meta_agg <- na.omit(meta_agg)


colnames(meta_agg) <- c("gsites", "lat", "long")


####plot IBD plot

# Make self comparisons NA
diag(pop_pFst$Fst) <- NA
diag(pop_pS$S) <- NA

# Mantel test
pop_man <- mantel(xdis = pop_pS$S, ydis = pop_pFst$Fst, permutations = 10000, na.rm = TRUE)
pop_man

# mantel plot
pop_Fst_sig <- cbind(melt(pop_pS$S), unlist(as.list(pop_pFst$Fst)))
colnames(pop_Fst_sig)[3] <- "Geo_dist"
colnames(pop_Fst_sig)[4] <- "Fst"
pop_Fst_sig$Geo_dist2 <- pop_Fst_sig$Geo_dist / 1000
pop_Fst_sig <- pop_Fst_sig[!is.na(pop_Fst_sig$Geo_dist),]

# # adding metadata for pop_larges
# pop_Fst_sig2 <- merge(pop_Fst_sig, distinct(meta_agg[, c("pop_large", "genetic_group","lat","long")]), by.x = "Var1", by.y = "pop_large", all.y = FALSE)
# pop_Fst_sig2 <- merge(pop_Fst_sig2, distinct(meta_agg[, c("pop_large", "genetic_group","lat","long")]), by.x = "Var2", by.y = "pop_large", all.y = FALSE)
# pop_Fst_sig2$same_sp <- ifelse(pop_Fst_sig2$genetic_group.x == pop_Fst_sig2$genetic_group.y, "Within group", "Between group")
# pop_Fst_sig2<- distinct(pop_Fst_sig2)  



# adding metadata for pop_larges
pop_Fst_sig2 <- merge(pop_Fst_sig, distinct(meta_agg[,c("gsites","lat","long")]), by.x = "Var1", by.y = "gsites", all.y = FALSE)
pop_Fst_sig2 <- merge(pop_Fst_sig2, distinct(meta_agg[, c("gsites","lat","long")]), by.x = "Var2", by.y = "gsites", all.y = FALSE)
# pop_Fst_sig2$same_sp <- ifelse(pop_Fst_sig2$genetic_group.x == pop_Fst_sig2$genetic_group.y, "Within group", "Between group")
pop_Fst_sig2<- distinct(pop_Fst_sig2)  





#Mapping
row_sub = apply(data.frame(fst_dms$meta$long), 1, function(row) all(row !=0 ))
newdmslong <- data.frame(fst_dms$meta$long)[row_sub,]
row_sub = apply(data.frame(fst_dms$meta$lat), 1, function(row) all(row !=0 ))
newdmslat <- data.frame(fst_dms$meta$lat)[row_sub,]

#this bit looks for range of latlong values to zoom in on and calulcates range for zoom for ggmap
divxlims <- c(min(na.omit(newdmslong))-1,max(na.omit(newdmslong))+1) #left-right direction   default for large distribution s around 0.6
divylims <- c(min(na.omit(newdmslat))-1,max(na.omit(newdmslat))+1) #up-down direction - default for large distribution s around 0.4
divxlimsrange <- abs(c(min(na.omit(newdmslong))-max(na.omit(newdmslong))))
divylimsrange <- abs(c(min(na.omit(newdmslat))-max(na.omit(newdmslat))))

if(divylimsrange > 5) {gmapzoom= 4}
if(divylimsrange > 4 && divylimsrange < 5)  {gmapzoom= 5}
if(divylimsrange > 3 && divylimsrange < 4)  {gmapzoom= 6}
if(divylimsrange > 2 && divylimsrange < 3)  {gmapzoom= 7}
if(divylimsrange > 1 && divylimsrange < 2)  {gmapzoom= 8}
if(divylimsrange > 0.2 && divylimsrange < 1) {gmapzoom= 10}
if(divylimsrange < 0.2) {gmapzoom= 12}


data <- data.frame(sample=dms$sample_names,
                   site=dms$meta$site,
                   analysis=dms$meta$analyses[,analysis],
                   lat=dms$meta$lat,
                   long=dms$meta$long)
data2 <- data[order(data$site),]
data3 <- data2[order(-data2$lat),]
data3$site <- factor(data3$site, levels = unique(data3$site),ordered = TRUE)

#this arranges the sites according to latitude
rep <- ddply(data3,
             .(analysis),
             summarise,
             lat  = mean(lat),
             long = mean(long))
rep2 <- rep[order(-rep$lat),]
rep2$analysis <-factor(rep2$analysis, levels=unique(rep2$analysis))

###colours for each group
rep_an <- ddply(data3,
                .(analysis),
                summarise,
                lat  = mean(lat),
                long = mean(long))
rep_an2 <- rep_an[order(-rep_an$lat),]
rval <- as.numeric(length(unique(data3$analysis)))    ##this counts the number(length) of unique characters for the variable set as analysis at the beginning
rep_an2$bgcols <- c(rev(matlab.like2(rval)))
rep_an3 <- rep_an2[,c("analysis","bgcols")]
rep_all <- merge(rep2,rep_an3, by= "analysis")
rep_all$analysis <-factor(rep_all$analysis, levels=unique(rep_an3$analysis))
rep_all <- rep_all[order(-rep_all$lat),]

for (q in 1:length(rep_all$analysis)){
  rep_all$number[q] <- q
  rep_all$siteAndNumber[q] <- paste0(rep_all$number[q]," ",rep_all$analysis[q]," (", length(which(data3$analysis==rep_all$analysis[q])), ")")
}

#googlemaps base and plot
library(ggmap)
register_google(key = "AIzaSyBZIYV071yUULZTtcUMUZTMwOgRJrKAcaw")

#gmapzoom = 5

if (length(rep_an3$analysis) > 30){legendcolumns = 3}
if (length(rep_an3$analysis) > 20 && length(rep_an3$analysis) < 30){legendcolumns = 2}
if (length(rep_an3$analysis) < 20){legendcolumns = 1}

googlemapbase = get_map(location = c(left = divxlims[1], bottom = divylims[1], right = divxlims[2], top=divylims[2]), zoom = gmapzoom, scale = "auto", maptype = "satellite", source = c("google"))

pop_Fst_sig2 <- pop_Fst_sig2[order(-pop_Fst_sig2$Fst),]


Fstthreshvals <-c(0.1, 0.2, 0.3, 0.4, 0.5)
pop_Fst_sig2$lineWid <- NA

for (m in 1:length(pop_Fst_sig2$Fst)){
  if (pop_Fst_sig2$Fst[m]<Fstthreshvals[1]){pop_Fst_sig2$lineWid[m] <- 0.1}
  if (pop_Fst_sig2$Fst[m]>Fstthreshvals[1] && pop_Fst_sig2$Fst[m]<Fstthreshvals[2]){pop_Fst_sig2$lineWid[m] <- 0.1}
  if (pop_Fst_sig2$Fst[m]>Fstthreshvals[2] && pop_Fst_sig2$Fst[m]<Fstthreshvals[3]){pop_Fst_sig2$lineWid[m] <-  0.8}
  if (pop_Fst_sig2$Fst[m]>Fstthreshvals[3] && pop_Fst_sig2$Fst[m]<Fstthreshvals[4]){pop_Fst_sig2$lineWid[m] <- 0.05}
  if (pop_Fst_sig2$Fst[m]>Fstthreshvals[4] && pop_Fst_sig2$Fst[m]<Fstthreshvals[5]){pop_Fst_sig2$lineWid[m] <- 0.01}
  if (pop_Fst_sig2$Fst[m]>Fstthreshvals[5]){pop_Fst_sig2$lineWid[m] <- 0.005}
}



# google map cropped by lat and long
ggmap(googlemapbase) + theme_bw()+
  geom_segment(data=pop_Fst_sig2[which(pop_Fst_sig2$Fst>=0.4&pop_Fst_sig2$Fst<1),], linewidth=1, alpha=0.5,na.rm = TRUE,
               # geom_segment(data=pop_Fst_sig2,
               aes(x=long.x, y=lat.x,xend = long.y, yend = lat.y,colour=Fst), lineend = "round")+
  geom_segment(data=pop_Fst_sig2[which(pop_Fst_sig2$Fst>=0.1&pop_Fst_sig2$Fst<0.5),], linewidth=1, alpha=0.5,na.rm = TRUE,
               # geom_segment(data=pop_Fst_sig2,
               aes(x=long.x, y=lat.x, xend = long.y, yend = lat.y, colour=Fst), lineend = "round")+
  geom_segment(data=pop_Fst_sig2[which(pop_Fst_sig2$Fst>=0&pop_Fst_sig2$Fst<0.1),], linewidth=1, alpha=1,na.rm = TRUE,
               # geom_segment(data=pop_Fst_sig2,
               aes(x=long.x, y=lat.x, xend = long.y, yend = lat.y,colour=Fst), lineend = "round")+
  # scale_linewidth_manual(name="Fst", values = as.character(pop_Fst_sig2$Fst)) +
  # scale_linewidth_manual(pop_Fst_sig2$Fst, limits=factor(0.0001, 0.8))+
  scale_color_gradient(low = "red", high = "white")+
  # geom_text(data=pop_Fst_sig2[pop_Fst_sig2$Fst<=0.5,],
  #           aes(label = round(Fst, 2), x = (long.x + long.y) / 2, y = (lat.y+lat.x) / 2),
  #           size = 2,
  #           color = "black")+
  geom_point(data=meta_agg, mapping=aes(x=long, y=lat), size=1, col="blue")+
  # xlim(lims[1], lims[2])+
  # ylim(lims[3], lims[4])+   
  labs(x = "Longitude", y = "Latitude", colour="FST") +guides(size = "none", alpha="none")
# ggrepel::geom_label_repel(data = meta_agg,aes(x = long, y = lat, label=gsites),
#                           min.segment.length=0.25, color="black",fill="white",size=3, segment.colour="white", alpha=0.9, label.size=0, nudge_y=0.004)

ggsave(paste0(species,"_", analysis, "Fst_Direction_Map.tiff"), path = paste0(fst_dir), width = 8, height = 12, dpi = 300, units = "in")



















##### Kinship analysis and UPGMA #####

if(!dir.exists(paste0(outputloc,"/Kinship/"))){
  dir.create(paste0(outputloc,"/Kinship/"))
}

m1 <- read.meta.data(d3, RandRbase, species, dataset, fields = as.numeric(ncol(metafile)-4))
temp<-as.data.frame(m1$analyses)

#set a sample threshold for each pop:
samplethreshold <-1# chnage this to what number of samples you want toi have uniform across all pops

## If you want to just remove any sites that are below the threshold (but keep individuals that are above a threshold, then unhash the below section)
for (x in unique(na.omit(temp[, analysis]))) {
  pop_samples <- which(temp[, analysis] == x)
  if (length(pop_samples) < samplethreshold) {
    temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
    print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
  } else {}
}
m1$analyses<-temp
dm        <- dart.meta.data.merge(d3, m1)
fields    <- c(analysis)
dms     <- data.by.meta.fields(dm, fields, RandRbase, species, dataset, object=analysis); print(dms$meta$analyses[,analysis])






dC <- dms

#UPGMA
library(ape)
library(phangorn)



rownames(dC$gt) <- paste( rownames(dms$gt),dms$meta$analyses[,analysis],sep="_")
SNAPPC <- dist(dC$gt)
tCUNNINGHAMII <- upgma(SNAPPC)
tCUNNINGHAMII <- ladderize(tCUNNINGHAMII, right = FALSE)
is_tip <- tCUNNINGHAMII$edge[,2] <= length(tCUNNINGHAMII$tip.label)
ordered_tips <- tCUNNINGHAMII$edge[is_tip, 2]
tip_order <- tCUNNINGHAMII$tip.label[ordered_tips]

par(mar=c(0,0,0,0))
plot(tCUNNINGHAMII, cex=0.1)
dev.off()

#kin
#note: if your populations are highly differentiated, it is best to run kinship on individual populations not the entire species.
#expect that a clonal population might have an excess of heterozygotes because if the genet has 50% heterozygotes,
#and the genet has multiple ramets, heterozygosity estimate of the population is more likely higher
#as compared to a non-clonal population of genets with different number of heterozygotes
#e.g. (50% x 3 vs 50%,30%,10%)


iIBD      <- individual.pw.IBD(dC,RandRbase,species,dataset)
kin       <-  iIBD$kinship #note if there is no kinship matrix,
#it could be because there are no SNPs to compare
#or the snpgdsIBDMoM() in individual.pw.IBD script doesnt have: kinship=TRUE

tiff(file = paste0(outputloc,"/Kinship/", analysis,"_Kinship_heatmap.tiff"),
     width = 11.7, height = 8.3, res = 300, units = "in")

#adds the row and column names corresponding to the tree to the spacial matrix
rownames(kin) <- paste0(dC$meta$sample_names,"_",dC$meta$analyses[,analysis])
colnames(kin) <- dC$meta$sample_names
# order kinship matrix using the tree
ik <- match(tip_order, rownames(kin))
ko <- kin[ik,ik]



#kinship heatmap
library(phytools)

HeatmapPlot <- phylo.heatmap(tCUNNINGHAMII,as.matrix(ko),fsize=0.5,ylim=c(0,1.1),split=c(0.2,0.80),grid=F)
HeatmapPlot

dev.off()

ko[upper.tri(ko, diag=T)] <- ""

ko[ko==0] <- ""

write.table(as.matrix(ko),file = paste0(outputloc, "/Kinship/", analysis," kinship.csv"),row.names = T,col.names = NA, sep = ",")



# order kinship matrix using the tree
ik <- match(tip_order, rownames(kin))
ko <- kin[ik,ik]

ko[upper.tri(ko, diag=T)] <- ""

ko[ko==0] <- ""

clones <- which(ko>0.45)

cloneslist <- c()

for (x in 1:length(clones)){
  k <- arrayInd(clones[x], dim(ko))
  n1 <- unlist(strsplit(rownames(ko)[k[,2]], "_"))[1]
  cloneslist$Ramet1[x] <- n1
  wn1 <- which(n1==dms$meta$sample_names)
  cloneslist$R1Site[x] <- paste0(dms$meta$analyses[,analysis][wn1])
  cloneslist$R1GPS[x] <- paste0(dms$meta$lat[wn1], ", ", dms$meta$long[wn1])
  
  cloneslist$Ramet2[x] <- colnames(ko)[k[,1]]
  n2 <- cloneslist$Ramet2[x]
  wn2 <- which(n2==dms$meta$sample_names)
  cloneslist$R2Site[x] <- paste0(dms$meta$analyses[,analysis][wn2])
  cloneslist$R2GPS[x] <- paste0(dms$meta$lat[wn2], ", ", dms$meta$long[wn2])
  
  cloneslist$Dist_Between_Points_Metres[x] <- distm(c(dms$meta$long[wn1], dms$meta$lat[wn1]), c(dms$meta$long[wn2], dms$meta$lat[wn2]))
  
}



write.table(cloneslist,file = paste0(outputloc, "/Kinship/", analysis," CloneList.csv"),row.names = F ,col.names = T, sep = ",")



##### DIVERSITY Stats with at least 5 samples per pop #####

m1 <- read.meta.data(d3, RandRbase, species, dataset, fields = as.numeric(ncol(metafile)-4))
temp<-as.data.frame(m1$analyses)

#set a sample threshold for each pop:
samplethreshold <-5# chnage this to what number of samples you want toi have uniform across all pops



# #### Below lines remove sites with less than threshold number of samples, and subsample sites that are above the threshold
# for (x in unique(na.omit(temp[, analysis]))) {
#   pop_samples <- which(temp[, analysis] == x)
#   if (length(pop_samples) < samplethreshold) {
#     temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
#     print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
#   }
#   if (length(pop_samples) == samplethreshold) {
#     print(paste0("Keeping ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
#   }
#   if (length(pop_samples) > samplethreshold) {
#     subsamplepops <- pop_samples[sample(1:length(pop_samples), size=samplethreshold, replace = F)]
#     for (b in 1:length(pop_samples)){
#       if (pop_samples[b] %in% subsamplepops == FALSE){
#         temp[, analysis][pop_samples[b]] <- NA
#       }
#       print(paste0("Reducing ",x, " to a sample size of ", samplethreshold," samples"))
#     }
#   } else {}
# }
# 
# 


## If you want to just remove any sites that are below the threshold (but keep individuals that are above a threshold, then unhash the below section)
for (x in unique(na.omit(temp[, analysis]))) {
  pop_samples <- which(temp[, analysis] == x)
  if (length(pop_samples) < samplethreshold) {
    temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
    print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
  } else {}
}

m1$analyses<-temp
dm        <- dart.meta.data.merge(d3, m1)
fields    <- c(analysis)
dms     <- data.by.meta.fields(dm, fields, RandRbase, species, dataset, object=analysis); print(dms$meta$analyses[,analysis])

#correct inconsistencies in J's scripts for naming treatment which is for naming files
treatment <- paste0("raw_SNPFilt_1SNPperClone_Field_", analysis)
dms$treatment <- treatment
gta <- dms$gt

gp   <- dart2genepop(dms, RandRbase, species, dataset, dms$meta$analyses[,analysis])

#calculate basic pop gen stats
bs <- basicStats(infile = gp, outfile = NULL,
                 fis_ci = FALSE, ar_ci = TRUE,
                 ar_boots = 999,
                 rarefaction = FALSE, ar_alpha = 0.05)

#create map plot and table

divrda_dir <- paste0(RandRbase,species,"/popgen/",treatment,"/genepop/bs.rda")
save(bs, file = divrda_dir)

if(!dir.exists(outputloc)){
  dir.create(outputloc)
}

tempdir <- paste0(outputloc, "/temp/")
if(!dir.exists(tempdir)){
  dir.create(tempdir)
}

Div_dir <- paste0(outputloc,"/diversity atleast 5 samples/")
if(!dir.exists(Div_dir)){
  dir.create(Div_dir)
}

{
  dmsanalysis <- dms$meta$analyses[,analysis]
  meta <- data.frame(table(dmsanalysis))
  z=meta$Freq<2
  
  if (sum(z, na.rm = TRUE)>=1) {#if there are any pops with less than 2 indiv...
    meta_lesthan2 <-as.character(meta[,1][meta$Freq<2]) #the pops with less than 2 indiv are...
    myvars <- dmsanalysis%in%meta_lesthan2
    newdmsanalysis <- dmsanalysis[!myvars]
    excluded_sample_names <- dms$sample_names[myvars]
    dms_d <- exclude.samples(dms, excluded_sample_names, remove_fixed_loci=TRUE)
    cat("metadata contains pops less than 2 individual","\n")
    print(meta_lesthan2)
  } else {
    meta_lesthan2 <-as.character(meta[,1][meta$Freq<2]) #the pops with less than 2 indiv are...
    myvars <- dmsanalysis%in%meta_lesthan2
    newdmsanalysis <- dms$meta$analyses[,analysis]
    dms_d <-dms
    excluded_sample_names <- NULL
    cat("metadata is alright without removal of pops","\n")
  }
  
  #making dataframe for  metadata without pops with less than 2 indiv
  dms_meta <- cbind.data.frame(dms$meta$sample_names, dms$meta$analyses[,analysis],dms$meta$lat, dms$meta$long)
  colnames(dms_meta) <- c("sample_names","site", "lat", "long")
  #pp <- dms_meta[!myvars,]
  #newdms_meta <- pp[apply(table(pp$site)>0,1,any),]
  newdms_meta <- dms_meta[!myvars,]
  
  colnames(newdms_meta) <- c("sample_names","site", "lat", "long")
  npop <-length(unique(newdms_meta$site))
  
  #making dataframe of averages (lat/longs) for metadata
  
  newmetalatlongs <- ddply(newdms_meta,
                           c("site"),
                           summarise,
                           lat  = mean(lat),
                           long  = mean(long),
                           N=length(site))
}

npop <- length(unique(newdmsanalysis))
result <- mat.or.vec(npop, 11)
measurement_names <- rownames(bs$main_tab[[1]])
population_names  <-
  names(bs$main_tab) #ls() rearranges the names
rownames(result) <- population_names
colnames(result) <- measurement_names

for (r in 1:npop) {
  popstats <- bs$main_tab[[r]][, "overall"] ##extract from a list
  result[r, ] <- popstats
}

result <- as.data.frame(result)
result$sample_names <- rownames(result)

###getting latlongs into the diversity results
metnew2 <- merge(newdms_meta[, 1:2], newmetalatlongs[, 1:4], by = "site")

data2 <- merge(result, metnew2, by = "sample_names", all.x = TRUE) #merge data
colnames(data2)[which(names(data2) == "dms.meta.analyses...analysis.")] <-
  "dms.meta.site"
data2$species <- as.character(species)

gp   <- dart2genepop0(dms, RandRbase, species, dataset, newdmsanalysis, maf_val = 0)

gp_genind <- read.genepop(gp, ncode = 2)

newdms_meta$site <- factor(newdms_meta$site, levels = rev(unique(newdms_meta$site)))
gp_genind@other <- newdms_meta
strata(gp_genind) <- gp_genind@other
setPop(gp_genind) <- ~ site

p_allele <-
  data.frame(rowSums(private_alleles(gp_genind, locus ~ site,
                                     count.alleles = F)))
p_allele$site <- rownames(p_allele)
rownames(p_allele) <- NULL
names(p_allele)[names(p_allele) == "rowSums.private_alleles.gp_genind..locus...site..count.alleles...F.."] <-
  "n_pa"
data <- merge(data2, p_allele, by.x = "site", by.y = "site")

write.table(data,
            paste0(Div_dir, "/", species, " diveRsity stats.csv"),
            sep = ",",
            row.names = F)

###plot

cexsize <- 0.5
ptsize <-0.5

tiff(paste0(outputloc,"/temp/",species," diveRsity stats2.tiff"), units="in", width=11.7, height=4, res=300)
par(mfrow=c(1,3), ## only works if plot out and not export tiff
    mar=c(1,1,1,1),
    oma=c(0,0,0.2,0)) #margins between nsw plots

row_sub_long = apply(data.frame(newdms_meta$long), 1, function(row) all(row !=0 ))
newdmslong <- data.frame(newdms_meta$long)[row_sub_long,]
row_sub_lat = apply(data.frame(newdms_meta$lat), 1, function(row) all(row !=0 ))
newdmslat <- data.frame(newdms_meta$lat)[row_sub_lat,]

if (file.exists(paste0(outputloc,"/temp/site_number.csv"))) {
  rep2 <- read.csv(paste0(outputloc,"/temp/site_number.csv"),sep=" ")
  
  testResult <- unlist(lapply(unique(dms$meta$analyses[,analysis]), function(s) {all(s %in% rep2$site)}))
  if (!all(testResult)) {
    dms_meta <- cbind.data.frame(dms$meta$sample_names, dms$meta$analyses[,analysis],dms$meta$lat, dms$meta$long)
    colnames(dms_meta) <- c("sample_names","site", "lat", "long")
    dms_meta$lat <- as.numeric(dms_meta$lat)
    dms_meta$long <- as.numeric(dms_meta$long)
    ###
    dms_meta2 <- dms_meta[order(dms_meta$site), ]
    dms_meta2 <- dms_meta2[order(-dms_meta2$lat), ]
    ####make sure your rep has nums_n_site which is in df2!!!
    dms_meta2$site <- factor(dms_meta2$site, levels = unique(dms_meta2$site),ordered = TRUE)
    
    rep <- ddply(dms_meta2,
                 .(site),
                 summarise,
                 lat  = mean(lat),
                 long = mean(long))
    rep2 <- rep[order(-rep$lat), ]
    x <- length(unique(dms_meta2$site))
    bgcols = rev(matlab.like2(x))
    rep2$nums_n_site <- paste0(1:x, ": ", rep2$site)
    rep2$num_site <- 1:x
    write.table(rep2, paste0(outputloc,"/temp/site_number_fromdiv.csv"))
  }
  
} else {
  dms_meta <- cbind.data.frame(dms$meta$sample_names, dms$meta$analyses[,analysis],dms$meta$lat, dms$meta$long)
  colnames(dms_meta) <- c("sample_names","site", "lat", "long")
  dms_meta$lat <- as.numeric(dms_meta$lat)
  dms_meta$long <- as.numeric(dms_meta$long)
  ###
  dms_meta2 <- dms_meta[order(dms_meta$site), ]
  dms_meta2 <- dms_meta2[order(-dms_meta2$lat), ]
  ####make sure your rep has nums_n_site which is in df2!!!
  dms_meta2$site <- factor(dms_meta2$site, levels = unique(dms_meta2$site),ordered = TRUE)
  rep <- ddply(dms_meta2,
               .(site),
               summarise,
               lat  = mean(lat),
               long = mean(long))
  rep2 <- rep[order(-rep$lat), ]
  x <- length(unique(dms_meta2$site))
  bgcols = rev(matlab.like2(x))
  rep2$nums_n_site <- paste0(1:x, ": ", rep2$site)
  rep2$num_site <- 1:x
  write.table(rep2, paste0(outputloc,"/temp/site_number_fromdiv.csv"))
}


data11 <- merge(data,rep2[,c("site","nums_n_site","num_site")],by="site")
data11$lab <- data11$num_site

divxlims <- c(min(na.omit(rep2$long))-0.4,max(na.omit(rep2$long))+0.4)
divylims <- c(min(na.omit(rep2$lat))-0.4,max(na.omit(rep2$lat))+0.4)

if ((divxlims[2]-divxlims[1]) <3) {
  rad <- 0.2
} else {
  rad <- 0.5}

## allelic richness
oz(xlim=divxlims, ylim=divylims, lwd=0.3)#divxlims taken from PCA section
title(expression(italic("ar")),font=2,cex.main=2)
points(x=data11$long, y=data11$lat,pch=16,bg="black",cex=ptsize, lwd=0.3)
draw.bubble(data11$long, data11$lat, (data11$ar)^50, bg=alpha("red",0.2), pch=21, maxradius = rad, lwd=0.3)
data11 <- data11[complete.cases(data11), ]#remove latlongs that are NAs
pointLabel(data11$long, data11$lat, labels = paste("  ", data11$lab, "  ", sep=""), cex=0.5,offset=0.25)
box()

## expected heterozygosity
oz(xlim=divxlims, ylim=divylims, lwd=0.3)
points(x=data11$long, y=data11$lat,pch=16,bg="black",cex=ptsize, lwd=0.3)
draw.bubble(data11$long, data11$lat, (data11$exp_het)^15, bg=alpha("red",0.2), pch=21, maxradius = rad, lwd=0.3)
title(expression(italic("exp_het")),font=4,cex.main=2)
pointLabel(data11$long, data11$lat, labels = paste("  ", data11$lab, "  ", sep=""), cex=0.5,offset=0.25)
box()

## observed heterozygosity
oz(xlim=divxlims, ylim=divylims, lwd=0.3)
points(x=data11$long, y=data11$lat,pch=16,bg="black",cex=ptsize, lwd=0.3)
pointLabel(data11$long, data11$lat, labels = paste("  ", data11$lab, "  ", sep=""), cex=0.5,offset=0.25)

if(min(data11$obs_het)-max(data11$obs_het)==0) {
  title(expression(italic("obs_het = 0")),font=4,cex.main=0.7)
} else {
  draw.bubble(data11$long, data11$lat, (data11$obs_het)^10, bg=alpha("red",0.2),pch=21, maxradius = rad, lwd=0.3)
  title(expression(italic("obs_het")),cex.main=2)
}
box()
dev.off()


#this bit checks for zeros in lats and longs and cut out those samples
row_sub = apply(data.frame(dms$meta$long), 1, function(row) all(row !=0 ))
newdmslong <- data.frame(dms$meta$long)[row_sub,]
row_sub = apply(data.frame(dms$meta$lat), 1, function(row) all(row !=0 ))
newdmslat <- data.frame(dms$meta$lat)[row_sub,]

#this bit looks for range of latlong values to zoom in on and calulcates range for zoom for ggmap
divxlims <- c(min(na.omit(newdmslong))-1,max(na.omit(newdmslong))+1)
divylims <- c(min(na.omit(newdmslat))-1,max(na.omit(newdmslat))+1)
divxlimsrange <- abs(c(min(na.omit(newdmslong))-max(na.omit(newdmslong))))
divylimsrange <- abs(c(min(na.omit(newdmslat))-max(na.omit(newdmslat))))


coordinate_cities <- read.csv(paste0(maindir,"CityMapCordinates.csv"))
sf_oz <- ozmap_data("states")

#order and map expected het as a coloured range
data11 <- data11[order(data11$exp_het),]

ggplot(sf_oz, aes(fill="white")) + geom_sf() +
  scale_fill_manual(values = "white", guide = "none") +
  theme_void() + geom_sf() +
  xlim(divxlims) + ylim(divylims)+
  geom_point(data = data11, mapping = aes(x = long, y = lat, color = exp_het), size=5, alpha=0.9) +
  scale_color_gradient(low="yellow", high="red", name="Exp Het") +
  geom_point(data=coordinate_cities, mapping = aes(x = long, y = lat), size=1)+
  geom_text(data=coordinate_cities, mapping = aes(x = long, y = lat, label=City),hjust=0,nudge_y=0.03, size=3)


ggsave(paste0(Div_dir, species, "_Expected_Het_Heatmap.tiff"), width = 10, height = 12, dpi = 300, units = "in")




#order and map observed het as a coloured range
data11 <- data11[order(data11$obs_het),]


ggplot(sf_oz, aes(fill="white")) + geom_sf() +
  scale_fill_manual(values = "white", guide = "none") +
  theme_void() + geom_sf() +
  xlim(divxlims) + ylim(divylims)+
  geom_point(data = data11, mapping = aes(x = long, y = lat, color = obs_het), size=5, alpha=0.9) +
  scale_color_gradient(low="yellow", high="red", name="Obs Het") +
  geom_point(data=coordinate_cities, mapping = aes(x = long, y = lat), size=1)+
  geom_text(data=coordinate_cities, mapping = aes(x = long, y = lat, label=City),hjust=0,nudge_y=0.03, size=3)


ggsave(paste0(Div_dir, species, "_Observed_Het_Heatmap.tiff"), width = 10, height = 12, dpi = 300, units = "in")




#order and map allelic richness as a coloured range
data11 <- data11[order(data11$ar),]


ggplot(sf_oz, aes(fill="white")) + geom_sf() +
  scale_fill_manual(values = "white", guide = "none") +
  theme_void() + geom_sf() +
  xlim(divxlims) + ylim(divylims)+
  geom_point(data = data11, mapping = aes(x = long, y = lat, color = ar), size=5, alpha=0.9) +
  scale_color_gradient(low="yellow", high="red", name="ar") +
  geom_point(data=coordinate_cities, mapping = aes(x = long, y = lat), size=1)+
  geom_text(data=coordinate_cities, mapping = aes(x = long, y = lat, label=City),hjust=0,nudge_y=0.03, size=3)


ggsave(paste0(Div_dir, species, "_Allelic_Richness_Heatmap.tiff"), width = 10, height = 12, dpi = 300, units = "in")







cat(paste0("all diversity maps are drawn!","\n"))
################combine div maps with a table
setwd(outputloc)
###table
geodist_data <- get_geodistinfo(newdms_meta)

data11 <- merge(data,rep2[,c("site","nums_n_site","num_site")],by="site")
data11$lab <- data11$num_site

if ("n_pa" %in% names(data11)){
}else{
  data11$n_pa <- data$rowSums.private_alleles.gp_genind..locus...site..count.alleles...F..
}

data3 <- merge(data11, geodist_data, by="site", all.x=TRUE)
data2grob <- cbind.data.frame(data3$num_site,data3$site,data3$lat,data3$long, data3$ar, data3$exp_het,data3$obs_het,data3$n_pa, data3$fis, data3$Freq, data3$meanDist)
colnames(data2grob) <- c("","site","lat","long","allelic_richness (ar)"," expected_het "," observed_het ","n_private_alleles","      fis      "," n_samples ","  meanDist  ")
data2grob <- data2grob[order(-data2grob$lat),]

ftsz <- 0.2
tiff(paste0("temp/",species," diveRsity grobtable.tiff"),units="in", width=11.7, height=5, res=300)
tt <- ttheme_default(base_size =7,
                     core=list(fg_params=list(fontface=3)),
                     colhead=list(fg_params=list(col="navyblue", fontface=2L)),
                     rowhead=list(fg_params=list(col="black", fontface=2L)))
tg <- tableGrob(data2grob,theme=tt,rows=NULL)
tg$heights <- rep(unit(0.2,"null"), nrow(tg))
grid.draw(tg)
dev.off()

cat(paste0("diversity table drawn!","\n"))

###now combine!

loc_diversity<- dir(path = "temp", pattern = "stats2.tiff$", full.names = TRUE, recursive = TRUE) ###must set working directory to make readPNG work
loc_diversitytab <- dir(path = "temp", pattern = "grobtable.tiff$", full.names = TRUE, recursive = TRUE) ###must set working directory to make readPNG work

lloo <- c(loc_diversity,loc_diversitytab)

plots <- lapply(ll <- lloo,function(x){
  img <- as.raster(readTIFF(x)) ##change this to readpng if readPNG doesnt work
  rasterGrob(img, interpolate = FALSE)
})

# lay <- rbind(c(1),c(1), c(2),c(2),c(2))

lay <- rbind(c(1),c(2))

ggsave(paste0(Div_dir,species, "_divntab_Atleast5Samps.tiff"),
       width=11.7, height=8.3,dpi = 300, units = "in", device='tiff',
       marrangeGrob(grobs = plots, layout_matrix = lay,top=textGrob(paste0(analysis, " Diversity Outputs With Atleast 5 Individuals per Pop or greater"), gp=gpar(fontsize=16,fontface=4))))

file.remove(lloo)

cat(paste0("diversity maps and table combined!","\n"))

setwd(maindir)



##### return a structure files from DArT genotype matrix to run using Structure #####
# dart2struct(dms, RandRbase, species, dataset, use_pops = FALSE)
# print(paste0(species," ", analysis," Number of Individuals:", length(dms$sample_names)))
# print(paste0(species," ", analysis," Number of Loci:", length(dms$locus_names)))
# #note: you need to run the dos2unix command to convert the .struct file if you are generating this file on r studio using a windows PC and then wanting to run it on the linux server





##### Neighbour Joining tree #####
m1 <- read.meta.data(d3, RandRbase, species, dataset, fields = as.numeric(ncol(metafile)-4))
temp<-as.data.frame(m1$analyses)

#set a sample threshold for each pop:
samplethreshold <-1# chnage this to what number of samples you want toi have uniform across all pops

## If you want to just remove any sites that are below the threshold (but keep individuals that are above a threshold, then unhash the below section)
for (x in unique(na.omit(temp[, analysis]))) {
  pop_samples <- which(temp[, analysis] == x)
  if (length(pop_samples) < samplethreshold) {
    temp[, analysis]<-gsub(x, replacement=NA, temp[, analysis])
    print(paste0("Removing ",x, " with ", " with a sample size of ",length(pop_samples)," samples"))
  } else {}
}
m1$analyses<-temp
dm        <- dart.meta.data.merge(d3, m1)
fields    <- c(analysis)
dms     <- data.by.meta.fields(dm, fields, RandRbase, species, dataset, object=analysis); print(dms$meta$analyses[,analysis])


if(!dir.exists(paste0(outputloc,"/NJTree/"))){
  dir.create(paste0(outputloc,"/NJTree/"))
}

gl_gm  <- dart2gl(dms, RandRbase, species, dataset)
tre <- nj(dist(as.matrix(gl_gm))) #generates  neighbor-joining tree estimate from euclidean matrix of the genotype data (genlight format)
# plot(tre, type="fan", cex=0.5)
plot(tre, type="phylogram", cex=0.5)

##########colour code according to PCA values
#PCA
nd_gl  <- dart2gl(dms, RandRbase, species, dataset)
nd_pca <- glPca(nd_gl, nf=5, parallel=FALSE)
myCol <- colorplot(nd_pca$scores,nd_pca$scores, transp=FALSE, cex=2)


tiff(paste0(outputloc,"/NJTree/",species," NJ_Tree_population.tiff"), units="in", width=11.7, height=8.3, res=300)

# scatter(nd_pca, xax=2,yax=3, posi="topleft", cex=0.1)

par(mfrow=c(1,2))
plot(nd_pca$scores[,1], nd_pca$scores[,2], xlab="PC1", ylab="PC2",col=myCol, pch=16)
text(nd_pca$scores[,1], nd_pca$scores[,2], labels=dms$meta$analyses[,analysis], pos=3, cex=0.5)

#NJ tree
plot(tre, typ="phylogram", show.tip=TRUE, no.margin=TRUE, cex=0.7)
tiplabels(pch=20, col=myCol, cex=2)


dev.off()



tiff(paste0(outputloc,"/NJTree/",species,"_NJ_Tree_PCA.tiff"), units="in", width=11.7, height=7, res=300)

##########colour code according to latitude
dms_meta <- cbind.data.frame(dms$meta$sample_names, dms$meta$analyses[,analysis],dms$meta$lat, dms$meta$long)
colnames(dms_meta) <- c("sample_names","site", "lat", "long")

rbPal <- colorRampPalette(c('blue','red'))
dms_meta$col <- rbPal(50)[as.numeric(cut(dms_meta$lat,breaks = 50))]

#NJ tree
par(mfrow=c(1,2), ## only works if plot out and not export tiff
    mar = c(5, 5, 5, 2),
    oma=c(0, 0, 0, 0)) #margins between plot

PCAVals <- data.frame(cbind(nd_pca$scores[,1], nd_pca$scores[,2]))
colnames(PCAVals) <- c("PC1", "PC2")

plot(PCAVals, xlab="PC1", ylab="PC2",col=dms_meta$col, pch=16)
pointLabel(PCAVals, labels=dms$sample_name, offset=0.1, cex=0.2)
mtext(paste0(analysis), side=3, line=2, cex=1)
mtext(paste0("NJ Tree with Latitude PCA"), side=3, line=3, cex=1.5)


par(mar = c(0, 0, 0, 0))
plot(tre, typ="phylogram", show.tip=TRUE, no.margin=TRUE, cex=0.3)
tiplabels(pch=20, col=dms_meta$col, cex=1)


dev.off()













##### Create summary PDF of analyses for species ##### 
maindir<- "E:/rrspecies/"
setwd(maindir)

# combine_plots(species)

####combine figures to get species summary

loc_QC <-dir(path = outputloc, pattern = " QCdataALL.tiff$", full.names = TRUE, recursive = TRUE) ###must set working directory to make readtiff work
loc_diversityAtleast5Samps <- dir(path = outputloc, pattern = "_divntab_Atleast5Samps.tiff$", full.names = TRUE, recursive = TRUE) 
loc_diversityLimiting5Samps <- dir(path = outputloc, pattern = "_divntab_limting5Samps.tiff$", full.names = TRUE, recursive = TRUE) 
loc_PCA_Lat <- dir(path = paste0(outputloc,"/PCA/"), pattern = "_PC1_PC2_PC3.tiff$", full.names = TRUE, recursive = TRUE) 
loc_PCA_wMap <- dir(path = paste0(outputloc,"/PCA/"), pattern = "_PCA_pops.tiff$", full.names = TRUE, recursive = TRUE)
loc_PCA_1_2 <- dir(path = paste0(outputloc,"/PCA/"), pattern = "_PC1-PC2.tiff$", full.names = TRUE, recursive = TRUE)
loc_PCA_Combined <- dir(path = paste0(outputloc,"/PCA/"), pattern = "_Combined PCA.tiff$", full.names = TRUE, recursive = TRUE) 
loc_PCA_Centroidmap <- dir(path = paste0(outputloc,"/PCA/"), pattern = "Dist_from_PCA12_Centroid.tiff$", full.names = TRUE, recursive = TRUE)
loc_PCA_Centroid <- dir(path = paste0(outputloc,"/PCA/"), pattern = "_median_PC123.tiff$", full.names = TRUE, recursive = TRUE)
loc_DistributionMap <- dir(path = paste0(outputloc,"/Map/"), pattern = "_Distribution_Map.tiff$", full.names = TRUE, recursive = TRUE)
loc_NJTree <- dir(path = paste0(outputloc,"/NJTree/"), pattern = "_NJ_Tree_PCA.tiff$", full.names = TRUE, recursive = TRUE)
loc_Kinship <- dir(path = paste0(outputloc, "/Kinship/"), pattern = "_kinship_heatmap.tiff$", full.names = TRUE, recursive = TRUE)
loc_LEAbar <-dir(path = outputloc, pattern = "\\lea k2345loc bar test.tiff$", full.names = TRUE, recursive = TRUE) 
loc_Fst_atleast1Samps <- dir(path = outputloc, pattern = "_fstplot_atleast1Samps.tiff$", full.names = TRUE, recursive = TRUE)
loc_Fst_limting5Samps <- dir(path = outputloc, pattern = "_fstplot_limiting5Samps.tiff$", full.names = TRUE, recursive = TRUE)
loc_heatmapsFst_atleast1Samps <- dir(path = outputloc, pattern = "_fstheatmaps_atleast1Samps.tiff$", full.names = TRUE, recursive = TRUE)
loc_heatmapsFst_limting5Samps <- dir(path = outputloc, pattern = "_fstheatmaps_limiting5Samps.tiff$", full.names = TRUE, recursive = TRUE)




locs <- c(loc_QC,
          loc_DistributionMap,
          loc_PCA_wMap,
          loc_PCA_Lat,
          loc_PCA_1_2,
          loc_PCA_Combined, 
          loc_PCA_Centroid,
          loc_PCA_Centroidmap,
          loc_LEAbar, 
          loc_diversityAtleast5Samps,
          loc_diversityLimiting5Samps, 
          loc_Fst_atleast1Samps,
          loc_Fst_limting5Samps,
          loc_heatmapsFst_atleast1Samps,
          loc_heatmapsFst_limting5Samps,
          loc_NJTree, 
          loc_Kinship)


for (i in 1:length(locs)) {
  if (!file.exists(locs[i])){
    cat(paste0(locs[i]," doesn't exist","\n"))
  }
}

plots <- lapply(ll <- locs,function(x){
  img <- as.raster(readTIFF(x))
  rasterGrob(img, interpolate = FALSE)
})
ggsave(paste0(outputloc2,species, "_", analysis, "_outputs.pdf"),width=7, height=5,
       marrangeGrob(grobs = plots, nrow=1, ncol=1,top=NULL))#width=7, height=3
print(locs)
cat(paste0("dart summary generated in ", outputloc,"/",species, "_", analysis, "_outputs.pdf!"))




##### Migrate Analysis with at least 5 samples per pop #####

gst_dir <- paste0(outputloc, "/gst_migrate/")
if (!dir.exists(gst_dir)) {
  dir.create(gst_dir)
}

source("C:/Users/dimonr/OneDrive - DPIE/R/Packages/sos_functions.R")


# c <- make_genepop_file(dms_pf, maf=0.05, missing=0.2, group=paste0(species, analysis), grouping=dms_pf$meta$analyses[,analysis])



c   <- dart2genepop(dms, RandRbase, species, dataset, dms$meta$analyses[,analysis])


# Diversity was measured as mean allelic richness for every population using the resampling method in the divBasic function in the R package diveRsity (46); to convert these diversity values into the data structure of pairwise relationships between populations, we calculated genetic diversity ratios for each population pair in each direction as the ratio of allelic richness in the destination versus the origin population. Pairwise genetic similarity was measured as 1/Fst.

# diveRsity::divBasic()
# AR destination/AR origin for each population 
# pairwise genetic similarity 1/fst 


#Migration was estimated using the divMigrate method (19) with Jost’s D metric of differentiation (47), as implemented in the R package diveRsity (46). This approach uses allele frequency differences between population pairs to estimate rates of migration in each direction; note that these rates are relative to other population pairs in the same data set and cannot be compared across data sets.

v <- diveRsity::divMigrate(infile=c, outfile=NULL, stat="gst",plot_network=TRUE, filter_threshold = 0.2, boots=1000, para=TRUE)



d_mig <- v$gRelMig #v$dRelMig

# from is rows, to is columns
colnames(d_mig) <- unique(dms$meta$analyses[,analysis]) 
rownames(d_mig) <- unique(dms$meta$analyses[,analysis]) 

# Filter by significance -- significance is if there is a significant difference in the directions
# Test for overlap of the estimated 95% confidence intervals. Where there is no overlap, the directional gene flow components are said to be significantly different (asymmetric).
mig_sig <- v$gRelMigSig #v$dRelMigSig
colnames(mig_sig) <- unique(dms$meta$analyses[,analysis])
rownames(mig_sig) <- unique(dms$meta$analyses[,analysis])

d_mig[mig_sig>0.05] <- NA

# qgraph::qgraph(d_mig,legend = F, edge.labels = F,
# curve = 2.5, mar = c(2, 2, 5, 5))

long_mig <- melt(d_mig)

# meta_agg <- m2 %>%
#   group_by(pop_large_short,pop_large) %>%
#   summarize(lat = mean(lat, na.rm=TRUE),
#             long = mean(long,na.rm=TRUE),
#             .groups = 'drop')%>%
#   subset(.,pop_large_short!="Ex_situ_PF")



dms_df <- data.frame(cbind(dms$meta$analyses[,analysis],dms$meta$lat, dms$meta$long))
colnames(dms_df) <- c("ansites", "lat", "long")

dms_df$lat <- as.numeric(dms_df$lat)
dms_df$long <- as.numeric(dms_df$long)


meta_agg <- dms_df %>% 
  group_by(dms_df$ansites) %>% 
  dplyr::summarize(lat = mean(lat, na.rm=TRUE),
                   long = mean(long, na.rm=TRUE), 
                   .groups = 'drop')

meta_agg <- na.omit(meta_agg)


colnames(meta_agg) <- c("gsites", "lat", "long")

long_mig <- merge(long_mig, distinct(meta_agg[, c("gsites","lat","long")]), by.x = "Var1", by.y = "gsites", all.y = FALSE)
long_mig <- merge(long_mig, distinct(meta_agg[, c("gsites","lat","long")]), by.x = "Var2", by.y = "gsites", all.y = FALSE)
long_mig<- distinct(long_mig)  

colnames(long_mig)[1:3] <- c("from", "to", "gst")

# long_mig <- long_mig[long_mig$from!="Ex_situ_PF"&long_mig$to!="Ex_situ_PF"&!is.na(long_mig$gst),]

long_mig <- long_mig[!is.na(long_mig$gst),] # remove any NA values between the same sites



# 
# map <- ggmap(get_map(location = bound, source = "stamen",
#                      zoom = 14, scale = 4,
#                      maptype = "terrain",
#                      color = "bw"))
# ggplot()+
# gst_map <-map+coord_cartesian()+coord_fixed()+
#   theme_bw()+
#   geom_curve(data=long_mig[long_mig$gst>0.3,],#
#                aes(x=long.x, y=lat.x,
#                    xend = long.y, yend = lat.y, colour=gst), #alpha=0.7+gst,colour=gst, 
#               size = 0.9,na.rm = TRUE,curvature=0.3,#lineend = "round",
#              # curvature=-0.3,,
#                arrow=arrow(angle=20, ends="first",type="closed", length=unit(2, "mm")))+
#   scale_color_gradient(low = "white", high = "red")+ #midpoint=0.5, mid="blue",
#   # geom_text(data=long_mig[long_mig$gst>=0.27,],
#   #           aes(label = round(gst, 2), x = (long.x + long.y) / 2, y = (lat.y+lat.x) / 2),
#   #           size = 2,
#   #           color = "black")+
#   geom_point(data=meta_agg, mapping=aes(x=long, y=lat), colour="black")+
#   xlim(lims[1], lims[2])+
#   ylim(lims[3], lims[4])+   labs(x = "Longitude", y = "Latitude", colour="Gst") +guides(size = "none", alpha="none")+
#   ggrepel::geom_label_repel(data = meta_agg,aes(x = long, y = lat, label=pop_large_short),
#                             min.segment.length=0.25, color="black",fill="white",size=3, segment.colour="white",
#                             alpha=0.9, label.size=0, nudge_y=0.003)
# 
# 
# ggsave("PherFitz/outputs/plots/gst_map.png", plot = gst_map, width = 200, height = 110, dpi = 300, units = "mm")
# 










#Mapping
row_sub = apply(data.frame(dms$meta$long), 1, function(row) all(row !=0 ))
newdmslong <- data.frame(dms$meta$long)[row_sub,]
row_sub = apply(data.frame(dms$meta$lat), 1, function(row) all(row !=0 ))
newdmslat <- data.frame(dms$meta$lat)[row_sub,]

#this bit looks for range of latlong values to zoom in on and calulcates range for zoom for ggmap
divxlims <- c(min(na.omit(newdmslong))-1,max(na.omit(newdmslong))+1) #default for large distribution s around 0.6
divylims <- c(min(na.omit(newdmslat))-1,max(na.omit(newdmslat))+1) #default for large distribution s around 0.4
divxlimsrange <- abs(c(min(na.omit(newdmslong))-max(na.omit(newdmslong))))
divylimsrange <- abs(c(min(na.omit(newdmslat))-max(na.omit(newdmslat))))

if(divylimsrange > 5) {gmapzoom= 4}
if(divylimsrange > 4 && divylimsrange < 5)  {gmapzoom= 5}
if(divylimsrange > 3 && divylimsrange < 4)  {gmapzoom= 6}
if(divylimsrange > 2 && divylimsrange < 3)  {gmapzoom= 7}
if(divylimsrange > 1 && divylimsrange < 2)  {gmapzoom= 8}
if(divylimsrange > 0.2 && divylimsrange < 1) {gmapzoom= 10}
if(divylimsrange < 0.2) {gmapzoom= 12}

library(ggmap)
register_google(key = "AIzaSyBZIYV071yUULZTtcUMUZTMwOgRJrKAcaw")

googlemapbase = get_map(location = c(left = divxlims[1], bottom = divylims[1], right = divxlims[2], top=divylims[2]), zoom = gmapzoom, scale = "auto", maptype = "satellite", source = c("google"))
# 
# pop_Fst_sig2 <- pop_Fst_sig2[order(-pop_Fst_sig2$Fst),]



# 
# x <- d_mig
# x[x<0.2] <- NA
# # x <- x[!colnames(x)=="Ex_situ_PF", !rownames(x)=="Ex_situ_PF"]
# qgraph::qgraph(x,legend = TRUE, edge.labels = TRUE,
# curve = 2.5, mar = c(2, 2, 5, 5))

library(ggthemes)


long_mig <- long_mig[order(long_mig$gst),]


ggmap(googlemapbase) + theme_few()+
  coord_cartesian()+ coord_fixed()+
  geom_curve(data=long_mig[long_mig$gst>0.2,],#
             aes(x=long.x, y=lat.x,
                 xend = long.y, yend = lat.y, colour=gst), #alpha=0.7+gst,colour=gst, 
             linewidth = 0.5,na.rm = TRUE,curvature=0.2, alpha=0.8, #lineend = "round",
             # curvature=-0.3,,
             arrow=arrow(angle=20, ends="first",type="closed", length=unit(2, "mm")))+
  scale_color_gradient(low = "blue", high = "red")+ #midpoint=0.5, mid="blue",
  # geom_text(data=long_mig[long_mig$gst>=0.27,],
  #           aes(label = round(gst, 2), x = (long.x + long.y) / 2, y = (lat.y+lat.x) / 2),
  #           size = 2,
  #           color = "black")+
  geom_point(data=meta_agg, mapping=aes(x=long, y=lat), colour="black", size=1)+
  # xlim(lims[1], lims[2])+
  # ylim(lims[3], lims[4])+   
  labs(x = "Longitude", y = "Latitude", colour="GST") +guides(size = "none", alpha="none")
# ggrepel::geom_label_repel(data = meta_agg,aes(x = long, y = lat, label=gsites),
#                           min.segment.length=0.25, color="black",fill="white",size=1, segment.colour="white",
#                           alpha=0.9, label.size=0, nudge_y=0.003)


# gst_fst <- ggarrange(fst_map, gst_no_map, nrow=2, labels=c("A","B"), align="hv")

ggsave(paste0(species,"_", analysis, "gst_Map.tiff"), path = paste0(gst_dir), width = 8, height = 12, dpi = 300, units = "in")








##### read Splitstree nexus and automate colour assignment ################################################################################
splitsdir <- paste0(outputloc,"/Splitstree/")
if(!dir.exists(splitsdir)){
  print("Nexus file does not exist. Have you created a splitstree file for this analysis?")
}
###to add some colour to your splitstree, save the nexus file in splitstree and then come back to R
library(phangorn)
source("C:/Users/dimonr/OneDrive - DPIE/R/Packages/RRtools/R/ggnetworx.R")
library(ggplot2)
library(ggtree)
library(ape)
library(patchwork)
#to open nexus here, first load the nexus in splitstree software and save it
# Nnet <- read.nexus.networx(paste0("E:/rrspecies/SplitstreeOutputs/",species,"_",analysis,"_SPLITSTREE.nex"))
Nnet <- read.nexus.networx(paste0(splitsdir,species,"_",analysis,"_SPLITSTREE.nex"))
nsw_num_o <- do.call(rbind,strsplit(Nnet$tip.label, "_NSW")) #get NSW number in splitstree nodes
usites <- unique(nsw_num_o[,1])
dmssplit <- dms
dmssplit$meta$analyses[, analysis]
dmscondensed <- gsub(' ','',dmssplit$meta$analyses[, analysis])
usiteslat <- c()
osites <- c()
datalist <- list()
for (i in 1:length(usites)){
  datalist[[i]] <- which(nsw_num_o == usites[i])
  names(datalist)[i] <- usites[i]; print(names(datalist)[i])
  allsitelats <- dmssplit$meta$lat[which(dmscondensed == usites[i])]
  usiteslat[i] <-  mean(allsitelats)
  osites$lat[i] <- usiteslat[i]
  osites$site[i] <- usites[i]
}

osites2 <- cbind(osites$lat,osites$site)
colnames(osites2) <- c("lat", "site")
dfsplits <- data.frame(osites2, stringsAsFactors = FALSE)
dfsplits$lat <- as.numeric(dfsplits$lat)
dfsplits$site <- as.character(dfsplits$site)
dfsplits <- dfsplits[order(-dfsplits$lat),]
for (x in 1:length(usites)){
  # dfsplits$colour[x] <- rainbow(length(usites))[x]
  dfsplits$colour[x] <- rev(matlab.like2(length(usites)))[x]
  dfsplits$shapes[x] <- shapeassign[x]
}
finalcols <- c()
finalshapes <- c()
for (x in 1:length(usites)){
  finalcols[x] <- dfsplits$colour[which(dfsplits$site==usites[x])]
  finalshapes[x] <- dfsplits$shapes[which(dfsplits$site==usites[x])]
  # dfsplits$colour[x] <- rainbow(length(usites))[x]
}
# usites <- as.character(dfsplits$site)
length(usites)
fs <- 4 ###font size of points
str <- 2 ###stroke size of points
ggsplitnet(Nnet) + geom_tiplab2(size=2)
x <- data.frame(x=Nnet$.plot$vertices[,1], y=Nnet$.plot$vertices[,2],
                sample=rep(NA, nrow(Nnet$.plot$vertices)), sample_names=rep(NA, nrow(Nnet$.plot$vertices)))
x[Nnet$translate$node,"sample"] <- Nnet$translate$label
x[Nnet$translate$node,"sample_names"] <- Nnet$translate$label
x$sample_names <- sub(x = x$sample_names,"_.*","")
x <- merge(x, m1, by="sample_names", all.x=TRUE, all.y=FALSE)
x<- x[!is.na(x$x),]
net_x_axis <- max(x$x)-min(x$x)
net_y_axis <- max(x$y)-min(x$y)
cf <- 1:length(usites)
cf <- order(as.character(cf), decreasing = FALSE)
colourgroup <- c("ggsplitnet(Nnet) + ")
for (i in 1:(length(usites))){
  text <- paste0("geom_tippoint(aes(subset = label %in% Nnet$translate$label[datalist[[", i, "]]], color =  'group", cf[i],"', shape = 'group", cf[i],"'), size=fs, stroke=str) +")
  colourgroup <- paste0(colourgroup, text)
  if (i==length(usites)){
    colourgroup <- paste0(colourgroup, "geom_tippoint(aes(subset = label %in% Nnet$translate$label[datalist[[", i , "]]], color = 'group", cf[i],"',shape = 'group", cf[i],"'), size=fs, stroke=str)")
  }
}
library(plyr)
A <- colourgroup
B <- function(x){}
body(B) <- as.quoted(A)[[1]]
B(1)
png(file = paste0(splitsdir, species, analysis," Splitstree_labelled.tiff"),
    width = 14, height = 14, res = 600, units = "in")
# network <- ggsplitnet(Nnet) + as.quoted(colourgroup)
B(1) +geom_splitnet(layout = "slanted", linewidth=0.2)+
  geom_tiplab2(size=1.5, hjust=-0.2)+
  scale_color_manual("Sites", values = finalcols, labels= usites)+
  scale_shape_manual("Sites", values = finalshapes, labels= usites)+
  expand_limits(x=c(min(x$x)-0.16*net_x_axis, max(x$x)+0.16*net_x_axis),
                y=c(min(x$y)-0.16*net_y_axis, max(x$y)+0.16*net_y_axis))+
  theme_void()+
  theme(legend.position="none")+
  # theme(legend.text = element_text(face="italic"), legend.position = "bottom") +
  coord_fixed()
dev.off()
png(file = paste0(splitsdir, species, analysis," Splitstree_unlabelled.tiff"),
    width = 14, height = 14, res = 600, units = "in")
# network <- ggsplitnet(Nnet) + as.quoted(colourgroup)
B(1) +geom_splitnet(layout = "slanted", linewidth=0.2)+
  # geom_tiplab2(size=1.5, hjust=-0.2)+
  scale_color_manual("Sites", values = finalcols, labels= usites)+
  scale_shape_manual("Sites", values = finalshapes, labels= usites)+
  expand_limits(x=c(min(x$x)-0*net_x_axis, max(x$x)+0*net_x_axis),
                y=c(min(x$y)-0*net_y_axis, max(x$y)+0*net_y_axis))+
  theme_void()+
  theme(legend.position="none")+
  # theme(legend.text = element_text(face="italic"), legend.position = "bottom") +
  coord_fixed()
dev.off()
#Splitstree with Site names only
Nnet <- read.nexus.networx(paste0(splitsdir,species,"_",analysis,"_SPLITSTREE.nex"))
Nnet$translate$label <- sub(x = Nnet$translate$label,".*_","") #removes NSW numbers
colourgroup <- c("ggsplitnet(Nnet) + ")
for (i in 1:(length(usites))){
  text <- paste0("geom_tippoint(aes(subset = label %in% Nnet$translate$label[datalist[[", i, "]]], color =  'group", cf[i],"', shape = 'group", cf[i],"'), size=fs, stroke=str) +")
  colourgroup <- paste0(colourgroup, text)
  if (i==length(usites)){
    colourgroup <- paste0(colourgroup, "geom_tippoint(aes(subset = label %in% Nnet$translate$label[datalist[[", i , "]]], color = 'group", cf[i],"',shape = 'group", cf[i],"'), size=fs, stroke=str)")
  }
}
library(plyr)
A <- colourgroup
B <- function(x){}
body(B) <- as.quoted(A)[[1]]
B(1)


png(file = paste0(splitsdir, species, analysis," Splitstree_NSWOnlyLabelled.tiff"),
    width = 14, height = 14, res = 600, units = "in")
# network <- ggsplitnet(Nnet) + as.quoted(colourgroup)
B(1) +geom_splitnet(layout = "slanted", linewidth=0.2)+
  geom_tiplab2(size=2, hjust=-0.3)+
  scale_color_manual("Sites", values = finalcols, labels= usites)+
  scale_shape_manual("Sites", values = finalshapes, labels= usites)+
  expand_limits(x=c(min(x$x)-0.16*net_x_axis, max(x$x)+0.16*net_x_axis),
                y=c(min(x$y)-0.16*net_y_axis, max(x$y)+0.16*net_y_axis))+
  theme_void()+
  theme(legend.position="none")+
  # theme(legend.text = element_text(face="italic"), legend.position = "bottom") +
  coord_fixed()
dev.off()
#Splitstree with Site names only
Nnet <- read.nexus.networx(paste0(splitsdir,species,"_",analysis,"_SPLITSTREE.nex"))
Nnet$translate$label <- sub(x = Nnet$translate$label,"_.*","") #removes NSW numbers
colourgroup <- c("ggsplitnet(Nnet) + ")
for (i in 1:(length(usites))){
  text <- paste0("geom_tippoint(aes(subset = label %in% Nnet$translate$label[datalist[[", i, "]]], color =  'group", cf[i],"', shape = 'group", cf[i],"'), size=fs, stroke=str) +")
  colourgroup <- paste0(colourgroup, text)
  if (i==length(usites)){
    colourgroup <- paste0(colourgroup, "geom_tippoint(aes(subset = label %in% Nnet$translate$label[datalist[[", i , "]]], color = 'group", cf[i],"',shape = 'group", cf[i],"'), size=fs, stroke=str)")
  }
}
library(plyr)
A <- colourgroup
B <- function(x){}
body(B) <- as.quoted(A)[[1]]
B(1)
png(file = paste0(splitsdir, species, analysis," Splitstree_SiteOnlyLabelled.tiff"),
    width = 14, height = 14, res = 600, units = "in")
# network <- ggsplitnet(Nnet) + as.quoted(colourgroup)
B(1) +geom_splitnet(layout = "slanted", linewidth=0.2)+
  geom_tiplab2(size=2, hjust=-0.3)+
  scale_color_manual("Sites", values = finalcols, labels= usites)+
  scale_shape_manual("Sites", values = finalshapes, labels= usites)+
  expand_limits(x=c(min(x$x)-0.16*net_x_axis, max(x$x)+0.16*net_x_axis),
                y=c(min(x$y)-0.16*net_y_axis, max(x$y)+0.16*net_y_axis))+
  theme_void()+
  theme(legend.position="none")+
  # theme(legend.text = element_text(face="italic"), legend.position = "bottom") +
  coord_fixed()
dev.off()
#Splitstree with numebred tips corresponding to sites ordered by latitude
Nnet <- read.nexus.networx(paste0(splitsdir,species,"_",analysis,"_SPLITSTREE.nex"))
Nnet$translate$label <- sub(x = Nnet$translate$label,"_.*","") #removes NSW numbers
for (d in 1:length(dfsplits$site)){
  Nnet$translate$label[which(Nnet$translate$label == dfsplits$site[d])] <- d
}
colourgroup <- c("ggsplitnet(Nnet) + ")
for (i in 1:(length(usites))){
  text <- paste0("geom_tippoint(aes(subset = label %in% Nnet$translate$label[datalist[[", i, "]]], color =  'group", cf[i],"', shape = 'group", cf[i],"'), size=fs, stroke=str) +")
  colourgroup <- paste0(colourgroup, text)
  if (i==length(usites)){
    colourgroup <- paste0(colourgroup, "geom_tippoint(aes(subset = label %in% Nnet$translate$label[datalist[[", i , "]]], color = 'group", cf[i],"',shape = 'group", cf[i],"'), size=fs, stroke=str)")
  }
}
library(plyr)
A <- colourgroup
B <- function(x){}
body(B) <- as.quoted(A)[[1]]
B(1)
png(file = paste0(splitsdir, species, analysis," Splitstree_Numbered.tiff"),
    width = 14, height = 14, res = 600, units = "in")
# network <- ggsplitnet(Nnet) + as.quoted(colourgroup)
B(1) + geom_splitnet(layout = "slanted", linewidth=0.2)+
  geom_tiplab2(size=5, hjust=-1, vjust=0)+
  scale_color_manual("Sites", values = finalcols, labels= usites)+
  scale_shape_manual("Sites", values = finalshapes, labels= usites)+
  # expand_limits(x=c(min(x$x)-0*net_x_axis, max(x$x)+0*net_x_axis),
  #               y=c(min(x$y)-0*net_y_axis, max(x$y)+0*net_y_axis))+
  theme_void()+
  theme(legend.position="none")+
  # theme(legend.text = element_text(face="italic"), legend.position = "bottom") +
  coord_fixed()
dev.off()