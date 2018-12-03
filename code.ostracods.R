
# load packages and data
library(ape); library(geiger); library(geomorph); library(phytools); library(dichromat)

#setwd("/Users/Bryan/Dropbox/cods")
setwd("")
tree = read.nexus("dat.tree.nex")
alldat = read.csv("dat.csv")
attach(alldat)


# PGLS --------------------------------------------------------------------

## length ~ depth 
keep = complete.cases(length, depth)
dat = data.frame(length = length[keep], depth = depth[keep], pelagic = pelagic[keep])
rownames(dat) = taxon[keep]
dat = dat[which(dat[, 3] != 2),] #clean data, remove aphotic zone species

treedat = treedata(tree, dat)  
gdf = geomorph.data.frame(length = log(treedat$data[, 1]), depth = log(treedat$data[, 2]), pelagic = as.factor(treedat$data[, 3]))
pgls.1a = procD.pgls(length ~ depth * pelagic, data = gdf, iter = 99999, phy = treedat$phy) #test for interaction

treedat = treedata(tree, dat[which(dat[, 3] == 0), c(1, 2)])
gdf = geomorph.data.frame(length = log(treedat$data[, 1]), depth = log(treedat$data[, 2]))
pgls.1b = procD.pgls(length ~ depth, data = gdf, iter = 99999, phy = treedat$phy) #interpret length~depth within levels of pelagic

treedat = treedata(tree, dat[which(dat[, 3] == 1), c(1, 2)])
gdf = geomorph.data.frame(length = log(treedat$data[, 1]), depth = log(treedat$data[, 2]))
pgls.1c = procD.pgls(length ~ depth, data = gdf, iter = 99999, phy = treedat$phy) #interpret length~depth within levels of pelagic

### photic zone eyed/eyeless partitions
keep1 = complete.cases(length, depth, count)
dat = data.frame(length[keep1], depth[keep1], pelagic[keep1], count[keep1])
keep2 = which(dat[, 3] == 0)
dat = data.frame(length = dat[, 1][keep2], depth = dat[, 2][keep2], count = dat[, 4][keep2])
rownames(dat) = alldat[, 1][keep1][keep2]
dat[which(dat[, 3] > 0), 3] = 1
colnames(dat)[3] = "eyed" #clean data

treedat = treedata(tree, dat)  
gdf = geomorph.data.frame(length = log(treedat$data[, 1]), depth = log(treedat$data[, 2]), eyed = as.factor(treedat$data[, 3]))
pgls.1d = procD.pgls(length ~ depth * eyed, data = gdf, iter = 99999, treedat$phy) #test for interaction

treedat = treedata(tree, dat[which(dat[, 3] > 0),])  
gdf = geomorph.data.frame(length = log(treedat$data[, 1]), depth = log(treedat$data[, 2]))
pgls.1e = procD.pgls(length ~ depth, data = gdf, iter = 99999, treedat$phy) #examine length~depth within levels of eyed

treedat = treedata(tree, dat[which(dat[, 3] == 0),])  
gdf = geomorph.data.frame(length = log(treedat$data[, 1]), depth = log(treedat$data[, 2]))
pgls.1f = procD.pgls(length ~ depth, data = gdf, iter = 99999, treedat$phy) #examine length~depth within levels of eyed

### summarize results
pgls.1 = list(pgls.1a, pgls.1b, pgls.1c, pgls.1d, pgls.1e, pgls.1f)
pgls.1.aov = lapply(pgls.1, function(x) {for(i in 1:length(x)) {return(x$aov.table)}})
pgls.1.aov = rbind(pgls.1.aov[[1]], "", pgls.1.aov[[2]], "", pgls.1.aov[[3]], "", pgls.1.aov[[4]], "", pgls.1.aov[[5]], "", pgls.1.aov[[6]])
pgls.1.coeff = sapply(pgls.1, function(x) {for(i in 1:length(x)) {return(x$pgls.coefficients)}})

### plot length ~ depth
pdf("length.pdf", width = 4, height = 5)
par(mfrow = c(2, 1))
keep = complete.cases(length, depth)
dat = data.frame(length = length[keep], depth = depth[keep], pelagic = pelagic[keep])
rownames(dat) = taxon[keep]
dat = dat[which(dat[, 3] != 2),]
treedat = treedata(tree, dat) #clean data, remove aphotic zone species

par(mar = c(5, 4, 4, 2) + 0.1 + c(-0.5, 0.5, -3, 0)) 
plot(log(treedat$data[, 2]), log(treedat$data[, 1]), 
	xaxt = "n", yaxt = "n", 
	pch = c(16, 1)[as.factor(treedat$data[, 3])], 
	xlab = "", ylab = "") #plot data
abline(pgls.1.coeff[[2]], lwd = 2) #plot slopes
abline(pgls.1.coeff[[3]], lwd = 2, lty = 2)
par(new = T)
plot(treedat$data[, 2], treedat$data[, 1], 
	log = "xy", type = "n", 
	xlab = "Depth (m)", ylab = "Carapace Length (mm)", cex.lab = 1) #plot log-scale axes
legend("bottomright", legend = c("Photic", "Dysphotic"), pch = c(16, 1), lty = c(1, 2), cex = 0.5, box.col = "gray")

### plot eyed/eyeless
keep1 = complete.cases(length, depth, count)
dat = cbind(length[keep1], depth[keep1], pelagic[keep1], count[keep1])
keep2 = which(dat[, 3] == 0)
dat = cbind(dat[, 1][keep2], dat[, 2][keep2], dat[, 4][keep2])
rownames(dat) = alldat[, 1][keep1][keep2]
colnames(dat) = c("length", "depth", "count")
dat[which(dat[, 3] > 0), 3] = 1
colnames(dat)[3] = "eyed"
treedat = treedata(tree, dat) #clean data

par(mar = (c(5, 4, 4, 2) + 0.1 + c(-0.5, 0.5, -3, 0)) / 1)
plot(log(treedat$data[, 2]), log(treedat$data[, 1]), 
	xaxt = "n", yaxt = "n", 
	pch = c(1, 16)[as.factor(treedat$data[, 3])], 
	xlab = "", ylab = "") #plot data
abline(pgls.1.coeff[[5]], lwd = 2)#plot slopes
abline(pgls.1.coeff[[6]], lwd = 2, lty = 2)
par(new = T)
plot(treedat$data[, 2], treedat$data[, 1], 
	log = "xy", type = "n", 
	xlab = "Depth (m)", ylab = "Carapace Length (mm)", cex.lab = 1) #plot log-scale axes
legend("topleft", legend = c("Eyed species", "Eyeless species"), pch = c(16, 1), lty = c(1, 2), cex = 0.5, box.col = "gray")
dev.off()

## relative eye length ~ depth
keep1 = complete.cases(length, eyelength, depth)
dat = cbind(length[keep1], eyelength[keep1], depth[keep1])
keep2 = apply(dat, 1, function(row) {all(row != 0)})
dat = dat[keep2,]
dat = data.frame(length = dat[, 1], eyelength = dat[, 2], depth = dat[, 3], pelagic = pelagic[keep1][keep2])
rownames(dat) = taxon[keep1][keep2]
dat = dat[which(dat[, 4] != 2),] #clean data, remove aphotic zone species

treedat = treedata(tree, dat)  
gdf = geomorph.data.frame(rel.length = log(treedat$data[, 2] / treedat$data[, 1] / 1000), depth = log(treedat$data[, 3]), pelagic = as.factor(treedat$data[, 4]))
pgls.2a = procD.pgls(rel.length ~ depth * pelagic, data = gdf, iter = 99999, treedat$phy) #test for interaction

treedat = treedata(tree, dat[which(dat[, 4] == 0), c(1, 2, 3)])
gdf = geomorph.data.frame(rel.length = log(treedat$data[, 2] / treedat$data[, 1] / 1000), depth = log(treedat$data[, 3]))
pgls.2b = procD.pgls(rel.length ~ depth, data = gdf, iter = 99999, treedat$phy) #examine relative eye length~depth within levels of pelagic

treedat = treedata(tree, dat[which(dat[, 4] == 1), c(1, 2, 3)])
gdf = geomorph.data.frame(rel.length = log(treedat$data[, 2] / treedat$data[, 1] / 1000), depth = log(treedat$data[, 3]))
pgls.2c = procD.pgls(rel.length ~ depth, data = gdf, iter = 99999, treedat$phy) #examine relative eye length~depth within levels of pelagic

### summarize results
pgls.2 = list(pgls.2a, pgls.2b, pgls.2c)
pgls.2.aov = lapply(pgls.2, function(x) {for(i in 1:length(x)) {return(x$aov.table)}})
pgls.2.aov = rbind(pgls.2.aov[[1]], "", pgls.2.aov[[2]], "", pgls.2.aov[[3]])
pgls.2.coeff = sapply(pgls.2, function(x) {for(i in 1:length(x)) {return(x$pgls.coefficients)}})

### plot relative eye length
keep1 = complete.cases(length, eyelength, depth)
dat = cbind(length[keep1], eyelength[keep1], depth[keep1])
keep2 = apply(dat, 1, function(row) {all(row != 0)})
dat = dat[keep2,]
dat = data.frame(length = dat[, 1], eyelength = dat[, 2], depth = dat[, 3], pelagic = pelagic[keep1][keep2])
rownames(dat) = taxon[keep1][keep2]
dat = dat[which(dat[, 4] != 2),] #clean data, remove aphotic zone species
treedat = treedata(tree, dat)  

## eye length ~ depth
keep1 = complete.cases(eyelength, depth)
dat = cbind(eyelength[keep1], depth[keep1])
keep2 = apply(dat, 1, function(row) {all(row != 0)})
dat = dat[keep2,]
dat = data.frame(eyelength = dat[, 1], depth = dat[, 2], pelagic = pelagic[keep1][keep2])
rownames(dat) = taxon[keep1][keep2]
dat = dat[which(dat[, 3] != 2),] #clean data, remove aphotic zone species

treedat = treedata(tree, dat)
gdf = geomorph.data.frame(eyelength = log(treedat$data[, 1]), depth = log(treedat$data[, 2]), pelagic = as.factor(treedat$data[, 3]))
pgls.3a = procD.pgls(eyelength ~ depth * pelagic, data = gdf, iter = 99999, phy = treedat$phy) #test for interaction

treedat = treedata(tree, dat[which(dat[, 3] == 0), c(1, 2)])
gdf = geomorph.data.frame(eyelength = log(treedat$data[, 1]), depth = log(treedat$data[, 2]))
pgls.3b = procD.pgls(eyelength ~ depth, data = gdf, iter = 99999, treedat$phy) #examine eye length~depth within levels of pelagic

treedat = treedata(tree, dat[which(dat[, 3] == 1), c(1, 2, 3)])
gdf = geomorph.data.frame(eyelength = log(treedat$data[, 1]), depth = log(treedat$data[, 2]))
pgls.3c = procD.pgls(eyelength ~ depth, data = gdf, iter = 99999, treedat$phy) #examine eye length~depth within levels of pelagic

### summarize results
pgls.3 = list(pgls.3a, pgls.3b, pgls.3c)
pgls.3.aov = lapply(pgls.3, function(x) {for(i in 1:length(x)) {return(x$aov.table)}})
pgls.3.aov = rbind(pgls.3.aov[[1]], "", pgls.3.aov[[2]], "", pgls.3.aov[[3]])
pgls.3.coeff = sapply(pgls.3, function(x) {for(i in 1:length(x)) {return(x$pgls.coefficients)}})

## diam ~ depth ##
keep1 = complete.cases(diam, depth)
dat = cbind(diam[keep1], depth[keep1])
keep2 =  apply(dat, 1, function(row) {all(row != 0 )})
dat = dat[keep2,]
dat = data.frame(diam = dat[, 1], depth = dat[, 2], pelagic = pelagic[keep1][keep2])
rownames(dat) = taxon[keep1][keep2]
dat = dat[which(dat[, 3] != 2),] #clean data, remove aphotic zone species

treedat = treedata(tree, dat)  
gdf = geomorph.data.frame(diam = log(treedat$data[, 1]), depth = log(treedat$data[, 2]), pelagic = as.factor(treedat$data[, 3]))
pgls.4a = procD.pgls(diam ~ pelagic * depth, data = gdf, iter = 99999, phy = treedat$phy)

treedat = treedata(tree, dat[which(dat[, 3] == 0), c(1, 2)])
gdf = geomorph.data.frame(diam = log(treedat$data[, 1]), depth = log(treedat$data[, 2]))
pgls.4b = procD.pgls(diam ~ depth, data = gdf, iter = 99999, treedat$phy) #examine diam~depth within levels of pelagic

treedat = treedata(tree, dat[which(dat[, 3] == 1), c(1, 2, 3)])
gdf = geomorph.data.frame(diam = log(treedat$data[, 1]), depth = log(treedat$data[, 2]))
pgls.4c = procD.pgls(diam ~ depth, data = gdf, iter = 99999, treedat$phy) #examine diam~depth within levels of pelagic

### summarize results
pgls.4 = list(pgls.4a, pgls.4b, pgls.4c)
pgls.4.aov = lapply(pgls.4, function(x) {for(i in 1:length(x)) {return(x$aov.table)}})
pgls.4.aov = rbind(pgls.4.aov[[1]], "", pgls.4.aov[[2]], "", pgls.4.aov[[3]])
pgls.4.coeff = sapply(pgls.4, function(x) {for(i in 1:length(x)) {return(x$pgls.coefficients)}})

### plot
keep1 = complete.cases(diam, depth)
dat = cbind(diam[keep1], depth[keep1])
keep2 =  apply(dat, 1, function(row) {all(row != 0 )})
dat = dat[keep2,]
dat = data.frame(diam = dat[, 1], depth = dat[, 2], pelagic = pelagic[keep1][keep2])
rownames(dat) = taxon[keep1][keep2]
dat = dat[which(dat[, 3] != 2),] #clean data, remove aphotic zone species
treedat = treedata(tree, dat) 

ggfactor <- treedat$data[, 3]
ggfactor[ggfactor == 0] <- "Photic"
ggfactor[ggfactor == 1] <- "Dysphotic"
ggdat <- data.frame(diam = treedat$data[, 1], Zone = as.factor(ggfactor))
levels(ggdat$Zone) <- c("Photic", "Dysphotic")

linesdat <- data.frame(Zone = c("Photic", "Dysphotic"), medians = tapply(ggdat$diam, ggdat$Zone, median))

pdf("diam.pdf", width = 4, height = 3)

ggplot(ggdat, aes(x = diam, col = Zone)) + 
    geom_density() +
    geom_vline(data = linesdat, aes(xintercept = medians, color = Zone), linetype = "dashed") + 
    labs(x = expression(paste("Ommatidium Diameter (", mu, "m)"))) +
    scale_color_grey() + 
    theme_classic()

dev.off()

## count ~ depth
dat = data.frame(count, depth)[complete.cases(count, depth),]
keep = apply(dat, 1, function(row) {all(row != 0)})
dat = dat[keep,]
dat = data.frame(count = dat[, 1], depth = dat[, 2], pelagic = pelagic[complete.cases(count, depth)][keep])
rownames(dat) = taxon[complete.cases(count, depth)][keep]
dat = dat[which(dat[, 3] != 2),] #clean data, remove aphotic zone species

treedat = treedata(tree, dat)
gdf = geomorph.data.frame(count = log(treedat$data[ , 1]), depth = log(treedat$data[ , 2]), pelagic = as.factor(treedat$data[ , 3]))
pgls.5a = procD.pgls(count ~ depth * pelagic, data = gdf, iter = 99999, phy = treedat$phy) #test for interaction

treedat = treedata(tree, dat[which(dat[, 3] == 0), c(1, 2)])
gdf = geomorph.data.frame(count = log(treedat$data[, 1]), depth = log(treedat$data[, 2]))
pgls.5b = procD.pgls(count ~ depth, data = gdf, iter = 99999, phy = treedat$phy) #examine count~depth within levels of pelagic

treedat = treedata(tree, dat[which(dat[, 3] == 1), c(1, 2)])
gdf = geomorph.data.frame(count = log(treedat$data[, 1]), depth = log(treedat$data[, 2]))
pgls.5c = procD.pgls(count ~ depth, data = gdf, iter = 99999, phy = treedat$phy) #examine count~depth within levels of pelagic

### summarize results
pgls.5 = list(pgls.5a, pgls.5b, pgls.5c)
pgls.5.aov = lapply(pgls.5, function(x) {for(i in 1:length(x)) {return(x$aov.table)}})
pgls.5.aov = rbind(pgls.5.aov[[1]], "", pgls.5.aov[[2]], "", pgls.5.aov[[3]])
pgls.5.coeff = sapply(pgls.5, function(x) {for(i in 1:length(x)) {return(x$pgls.coefficients)}})

### plot 
dat = data.frame(count, depth)[complete.cases(count, depth),]
keep = apply(dat, 1, function(row) {all(row != 0)})
dat = dat[keep,]
dat = data.frame(count = dat[, 1], depth = dat[, 2], pelagic = pelagic[complete.cases(count, depth)][keep])
rownames(dat) = taxon[complete.cases(count, depth)][keep]
dat = dat[which(dat[, 3] != 2),] #clean data, remove aphotic zone species
treedat = treedata(tree, dat) 

pdf("count.pdf", width = 5, height = 5)

par(mar = c(5, 4, 4, 2) + 0.1 + c(-0.5, 0.5, 0, 0))
plot(log(treedat$data[, 2]), log(treedat$data[, 1]),
	xaxt = "n", yaxt = "n", 
	pch = c(16, 1)[as.factor(treedat$data[, 3])], 
	xlab = "", ylab = "") #plot data 
abline(pgls.5.coeff[[2]], lwd = 2) #plot slopes
abline(pgls.5.coeff[[3]], lwd = 2, lty = 2)
par(new = T)
plot(treedat$data[, 2], treedat$data[, 1], 
	log = "xy", type = "n", 
	xlab = "Depth (m)", ylab = "Ommatidia Count", cex.lab = 1.25) #plot log-scale axes
legend("topleft", legend = c("Photic", "Dysphotic"), pch = c(16, 1), lty = c(1, 2), cex = 0.9, box.col = "gray")

dev.off()

# Phylomorphospace --------------------------------------------------------

## count ~ depth
dat = data.frame(count, depth)[complete.cases(count, depth),]
keep = apply(dat, 1, function(row) {all(row != 0)})
dat = dat[keep,]
dat = data.frame(count = dat[, 1], depth = dat[, 2], pelagic = pelagic[complete.cases(count, depth)][keep])
rownames(dat) = taxon[complete.cases(count, depth)][keep]
treedat = treedata(tree, dat)

pdf("phylomorphospace.pdf", width = 6, height = 6)

par(mar = c(5, 4, 4, 2) + 0.1 + c(-0.5, 0.5, 0, 0)) #set plotting parameters
plot(log(dat[, 2]), log(dat[, 1]), type = "n", 
	xaxt = "n", yaxt = "n", xlab = "", ylab = "", cex.lab = 1.25) #plot space
source("code.phylomorph.R")
phylomorphospace.cods(treedat$phy, log(treedat$data[, c(2, 1)]), xlab = "", ylab = "", label = "off", node.size = c(0.5, .75), add = T) #plot phylomorphospace
par(new = T)
plot(dat[, 2], dat[, 1], 
	log = "xy", type = "n", 
	xlab = "Depth (m)", ylab = "Ommatidia Count", cex.lab = 1.25) #plot log-scale axes
points(103, 11.5, pch = 17, col = "gray62", cex = 2)
abline(v = 200, lty = 2)

dev.off()

# Stochastic Trait Mapping ------------------------------------------------
colfunc <- colorRampPalette(c("black", "white", "red"))

## count ~ depth
dat = data.frame(count, depth)[complete.cases(count, depth),]
keep = apply(dat, 1, function(row) {all(row != 0)})
dat = dat[keep,]
dat = data.frame(count = dat[, 1], depth = dat[, 2], pelagic = pelagic[complete.cases(count, depth)][keep])
rownames(dat) = taxon[complete.cases(count, depth)][keep]
treedat = treedata(tree, dat)

pdf("ace_color.pdf", width = 4, height = 4)

layout(matrix(1:3, 1, 3), widths = c(0.4, 0.2, 0.4)) #set plotting parameters
par(cex = 0.3)
obj1 = contMap(treedat$phy, treedat$data[, 1], legend = 0.1 * max(nodeHeights(treedat$phy)), plot = F, lims = c(1, 23)) #map trait 1
obj1 = setMap(obj1, colors = c(colfunc(8)))
plot(obj1, ftype = "off", fsize = c(0.1, 1.5), lwd = c(1.2, 1.75), leg.txt = "Ommatidia Count", cex = 5, legend = 200) #plot trait 1
drop = treedat$phy

ylim = c(1 - 0.12 * (length(drop$tip.label) - 1), length(drop$tip.label))
plot.new(); plot.window(xlim = c(-0.5, 0.5), ylim = ylim) #add new plot and set scale
is.tip = drop$edge[, 2] <= length(drop$tip.label)
ordered.tips = drop$edge[is.tip, 2]
spp = gsub("_", " ", drop$tip.label[ordered.tips])
spp = gsub(" J57062", "", spp)
text(rep(0, length(drop$tip.label)), 1:length(drop$tip.label), spp, cex = 0.85, font = c(rep(3, 43), rep(4, 13), rep(3, 35)))

obj2 = contMap(treedat$phy, treedat$data[, 2], legend = 0.7 * max(nodeHeights(treedat$phy)), plot = F, lims = c(1, 670)) #map trait 2
obj2 = setMap(obj2, colors = c(colorschemes$BluetoGray.8))
plot(obj2, direction = "leftwards", ftype = "off", fsize = c(0.1, 1.5), lwd = c(1.2, 1.75), leg.txt = "Depth (m)", legend = 200) #plot trait 2

dev.off()

## keep only focal clade
pdf("ace2.pdf", width = 4, height = 4)

layout(matrix(1:3, 1, 3), widths = c(0.4, 0.2, 0.4)) #set plotting parameters
par(cex = 0.3)
obj1 = contMap(treedat$phy, treedat$data[, 1], legend = 0.7 * max(nodeHeights(treedat$phy)), plot = F, lims = c(1, 23)) #map trait 1
obj1 = drop.tip.contMap(obj1, setdiff(treedat$phy$tip.label, extract.clade(treedat$phy, 125)$tip.label)) #drop tips
obj1 = setMap(obj1, colors = c("black", "white"))
plot(obj1, ftype = "off", fsize = c(0.1, 1.25), lwd = c(5, 1.75), leg.txt = "Ommatidia Count") #plot trait 1

ylim = c(1 - 0.12 * (13 - 1), 13)
plot.new(); plot.window(xlim = c(-0.1, 0.1), ylim = ylim) #add new plot and set scale
spp = gsub("_", " ", treedat$phy$tip.label[ordered.tips][which(treedat$phy$tip.label[ordered.tips] %in% extract.clade(treedat$phy, 125)$tip.label)])
text(rep(0, 13), 1:13, spp, cex = 1.2, font = 2)

obj2 = contMap(treedat$phy, treedat$data[, 2], legend = 0.7 * max(nodeHeights(treedat$phy)), plot = F, lims = c(1, 670)) #map trait 2
obj2 = drop.tip.contMap(obj2, setdiff(treedat$phy$tip.label, extract.clade(treedat$phy, 125)$tip.label)) #drop tips
obj2 = setMap(obj2, colors = c("white", "black"))
plot(obj2, direction = "leftwards", ftype = "off", fsize = c(0.1, 1.25), lwd = c(5, 1.75), leg.txt = "Depth (m)") #plot trait 2

dev.off()

# eyeless:eyed species ratio vs depth -------------------------------------

dat = alldat[complete.cases(count, pelagic), c("count", "pelagic")] #remove NAs
attach(dat)

rats = c(length(count[(count == 0 & pelagic == 0)]) / length(which(pelagic == 0)),
         length(count[(count == 0 & pelagic == 1)]) / length(which(pelagic == 1)), 
         length(count[(count == 0 & pelagic == 2)]) / length(which(pelagic == 2))) #obtain ratios

reps = 99999 #set number of replicates
sample.n = c(dim(dat[pelagic == 0,])[[1]], dim(dat[pelagic == 1,])[[1]], dim(dat[pelagic == 2,])[[1]]) #obtain number of species in each pelagic zone
n = sum(sample.n)
resample = array(dim = c(reps, 3)) #set up empty array to hold randomly drawn samples

for(i in 0:2){
	for(j in 1:reps){
		boot = sample(n, sample.n[[i + 1]]) #sample row numbers from entire dataset equal to the number of species living in each pelagix zone
		resample[j, i + 1] = with(dat[boot,], length(which(count == 0)) / dim(dat[boot,])[[1]]) #fill array with null distribution of ratios
	}
}

upper = apply(resample, 2, function(x) {quantile(x, 0.95)}) #calculate 95% CIs
lower = apply(resample, 2, function(x) {quantile(x, 0.05)})

1 - length(which(rats[[1]] < resample[, 1])) / (reps + 1) #calculate p-values
1 - length(which(rats[[2]] > resample[, 2])) / (reps + 1)
1 - length(which(rats[[3]] > resample[, 3])) / (reps + 1)

detach(dat)

##plot

resample2 <- data.frame(rbind(
    cbind(resample[,1], rep(1, 99999)), 
    cbind(resample[,2], rep(2, 99999)), 
    cbind(resample[,3], rep(3, 99999))
    ))

resample3 <- data.frame(ratios = resample2[,1], pelagic = as.factor(resample2[,2]))
levels(resample3$pelagic) <- c("Photic", "Dysphotic", "Aphotic")

vline.data <- data.frame(z = c(rats), pelagic = c("Photic", "Dysphotic", "Aphotic"))
dat_text <- data.frame(
  label = c("*", "", "*"),
  pelagic   = c("Photic", "Dysphotic", "Aphotic"),
  x     = c(0.3, 0, 0.3),
  y     = c(6.5, 0, 4)
)

pdf("permutation_test.pdf", width = 6, height = 6)

ggplot(resample3, aes(x = ratios)) + geom_density(fill = "gray", alpha = 0.5, bw = 0.05) + facet_grid(pelagic~.) + geom_vline(aes(xintercept = z), vline.data, lwd = 1.2) + labs(x="Ratio of Eyeless:Eyed Species") + theme_classic() + geom_text(data = dat_text, mapping = aes(x = x, y = y, label = label), size = 10)

dev.off()


