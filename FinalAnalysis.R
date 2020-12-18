##### Packages ####
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(resample)

#### Load basic data ####
subjects <- 
  read.csv("~/R/StructNetwNMDA/data/clincaldata_minimal.csv")
node_names <- read.csv("~/R/StructNetwNMDA/data/node_names.csv", sep="")

## excluded cases (LG1 antibodies and MS double diagnosis)
# NMDA 007/016/027/150/151/154; LE_HC 091/165/170; LE_GK_011; RS_GK 008/020

#### Global parameters: AUC data ####

### Small-worldness ###

## HC

hcSigmaAUC <- 
  read.delim("~/R/StructNetwNMDA/data/gretna/ThesisExclMS/SmallWorld/Group1/Sigma_All_Thres.txt", 
             colClasses = c(rep('numeric',10), 'NULL'),
             header=FALSE)

# mean, sd upper and lower bound
hcSigmaAUC_uppermean <- mean(hcSigmaAUC$V1)
hcSigmaAUC_lowermean <- mean(hcSigmaAUC$V10)
hcSigmaAUC_uppersd <- sd(hcSigmaAUC$V1)
hcSigmaAUC_lowersd <- sd(hcSigmaAUC$V10)

## NMDA

nmdaSigmaAUC <- 
  read.delim("~/R/StructNetwNMDA/data/gretna/ThesisExclMS/SmallWorld/Group2/Sigma_All_Thres.txt", 
             colClasses = c(rep('numeric',10), 'NULL'),
             header=FALSE)

# mean, sd upper and lower bound
nmdaSigmaAUC_uppermean <- mean(nmdaSigmaAUC$V1)
nmdaSigmaAUC_lowermean <- mean(nmdaSigmaAUC$V10)
nmdaSigmaAUC_uppersd <- sd(nmdaSigmaAUC$V1)
nmdaSigmaAUC_lowersd <- sd(nmdaSigmaAUC$V10)

### Global efficiency ###

## HC
hcGeffAUC <- 
  read.delim("~/R/StructNetwNMDA/data/gretna/ThesisExclMS/NetworkEfficiency/Group1/Eg_All_Thres.txt", 
             colClasses = c(rep('numeric',10), 'NULL'),
             header=FALSE)

# mean, sd upper and lower bound
hcGeffAUC_mean <- mean(unlist(hcGeffAUC))
hcGeffAUC_sd <- sd(unlist(hcGeffAUC))

## NMDA

nmdaGeffAUC <- 
  read.delim("~/R/StructNetwNMDA/data/gretna/ThesisExclMS/NetworkEfficiency/Group2/Eg_All_Thres.txt", 
             colClasses = c(rep('numeric',10), 'NULL'),
             header=FALSE)

# mean, sd upper and lower bound
nmdaGeffAUC_mean <- mean(unlist(nmdaGeffAUC))
nmdaGeffAUC_sd <- sd(unlist(nmdaGeffAUC))

## Global clustering coefficient (mean over all nodes over all thresholds over all subjects)




#### Node Strength #### ####

### read all matrices and convert to 3D array

files_hc = list.files(path = "./data/conn_matrices/HC", pattern = ".csv", full.names = TRUE, recursive = FALSE)
files_nmda = list.files(path = "./data/conn_matrices/NMDA", pattern = ".csv", full.names = TRUE, recursive = FALSE)

hc_mat = list()
nmda_mat = list()

# for some reason, the .csv files of HC had different internal separators 

i = 1
for(file_mat in files_hc) {
  
  current = read.table(file = file_mat, header = FALSE)
  
  if(length(current) == 1) {
    current = read.table(file = file_mat, header = FALSE, sep = ",")
  } 
  
  if(length(current) == 1) {
    current = read.table(file = file_mat, header = FALSE, sep = ";")
  } 
  
  hc_mat[[i]] = current
  i = i + 1
}

# convert list of matrices into array for calculations
hc_array <- array(unlist(hc_mat), c(84,84,67))

# NMDA Data
i = 1
for(file_mat in files_nmda) {
  
  current = read.table(file = file_mat)
  nmda_mat[[i]] = current
  i = i + 1
}

nmda_array <- array(unlist(nmda_mat), c(84,84,67))

### Strength analysis

## AUC data

hc_strAUC <- read.csv("~/R/StructNetwNMDA/data/conn_matrices/hc_strAUC_excl.csv", header=FALSE)
nmda_strAUC <- read.csv("~/R/StructNetwNMDA/data/conn_matrices/nmda_strAUC_excl.csv", header=FALSE)
str_AUC <- rbind(hc_strAUC, nmda_strAUC)
colnames(str_AUC) <- node_names$x

## get short version of df with only MTL, DMN nodes
str_AUC_s <- str_AUC[,c(2, 3, 5, 7, 13, 15, 22, 24, 25, 26, 40, 47, 51, 52, 54, 56, 62, 64, 71, 73, 74, 75)]
str_AUC_s$Group[1:61] <- "HC"
str_AUC_s$Group[62:122] <- "NMDA"

# prepare long version for boxplot
str_AUC_l <- str_AUC_s
str_AUC_l <- gather(str_AUC_l, colnames(str_AUC_l)[1:22], key = 'node', value = 'meanstr')

p_str <-
  ggplot(str_AUC_l, aes(x=node, y=meanstr, fill=Group)) + 
  geom_boxplot(position=position_dodge(1)) + 
  coord_flip()
print(p_str)

# Permutation testing
nodes <- colnames(str_AUC_s[-23])  

i <- 1
for(x in nodes) {
  
  perm <- permutationTest2(str_AUC_s, mean(eval(parse(text = x))), treatment = Group, seed = 0)
  
  pvals[i] <- perm$stats$PValue
  
  i = i + 1
}

#### Betweenness Centrality ####

### AUC Data

hcBC_AUC <- read.delim("./data/gretna/ThesisExclMS/BetweennessCentrality/Group1/aBc.txt", header=FALSE)
nmdaBC_AUC <- read.delim("./data/gretna/ThesisExclMS/BetweennessCentrality/Group2/aBc.txt", header=FALSE)
BC_AUC <- rbind(hcBC_AUC, nmdaBC_AUC)
# get rid of empty 85th column
BC_AUC <- BC_AUC[1:84]
colnames(BC_AUC) <- node_names$x

## get short version of df with only MTL, DMN nodes
BC_AUC_s <- BC_AUC[,c(2, 3, 5, 7, 13, 15, 22, 24, 25, 26, 40, 47, 51, 52, 54, 56, 62, 64, 71, 73, 74, 75)]
BC_AUC_s$Group[1:61] <- "HC"
BC_AUC_s$Group[62:122] <- "NMDA"

# prepare long version for boxplot
BC_AUC_l <- BC_AUC_s
BC_AUC_l <- gather(BC_AUC_l, colnames(BC_AUC_l)[1:22], key = 'node', value = 'meanBC')

p_BC <-
  ggplot(BC_AUC_l, aes(x=node, y=meanBC, fill=Group)) + 
  geom_boxplot(position=position_dodge(1)) + 
  coord_flip()
print(p_BC)

### Raw data
hcBC_raw <- read.delim("./data/gretna/NoThreshold/BetweennessCentrality/Group1/Bc.txt", header=FALSE)
nmdaBC_raw <- read.delim("./data/gretna/NoThreshold/BetweennessCentrality/Group2/Bc.txt", header=FALSE)
BC_raw <- rbind(hcBC_raw, nmdaBC_raw)
# get rid of empty 85th column
BC_raw <- BC_raw[1:84]
colnames(BC_raw) <- node_names$x

## get short version of df with only MTL, DMN nodes
BC_raw_s <- BC_raw[,c(2, 3, 5, 7, 13, 15, 22, 24, 25, 26, 40, 47, 51, 52, 54, 56, 62, 64, 71, 73, 74, 75)]
BC_raw_s$Group[1:61] <- "HC"
BC_raw_s$Group[61:122] <- "NMDA"

# prepare long version for boxplot
BC_raw_l <- BC_raw_s
BC_raw_l <- gather(BC_raw_l, colnames(BC_raw_l)[1:22], key = 'node', value = 'meanBC')

p_BC_raw <-
  ggplot(BC_raw_l, aes(x=node, y=meanBC, fill=Group)) + 
  geom_boxplot(position=position_dodge(1)) + 
  coord_flip()
print(p_BC_raw)

## Run significance tests
testnode <- select(BC_AUC_s, R.HI, Group)
wilcox.test(BC_AUC_s$L.HI ~ BC_AUC_s$Group)

lapply(BC_raw_s[-23], function(node) wilcox.test(node ~ BC_raw_s$Group))

node_plot <- ggplot(testnode, aes(y=R.HI, x=Group)) +
  geom_point(position=position_jitter(height=0, width=0.1)) +
  stat_summary(fun.y=mean, geom="point", color="red", size=3)
print(node_plot)

# Permutation testing
nodes <- colnames(BC_AUC_s[-23])  
pvals <- rep(0,22)

i <- 1
for(x in nodes) {
  
  perm <- permutationTest2(BC_AUC_s, mean(eval(parse(text = x))), treatment = Group, seed = 0)
  
  pvals[i] <- perm$stats$PValue
  
  i = i + 1
}


#### Clustering Coefficient #### 
hcCC_AUC <- read.delim("./data/gretna/ThesisExclMS/NodalClustCoeff/Group1/aNCp.txt", header=FALSE)
nmdaCC_AUC <- read.delim("./data/gretna/ThesisExclMS/NodalClustCoeff/Group2/aNCp.txt", header=FALSE)
CC_AUC <- rbind(hcCC_AUC, nmdaCC_AUC)
# get rid of empty 85th column
CC_AUC <- CC_AUC[1:84]
colnames(CC_AUC) <- node_names$x

## get short version of df with only MTL, DMN nodes
CC_AUC_s <- CC_AUC[,c(2, 3, 5, 7, 13, 15, 22, 24, 25, 26, 40, 47, 51, 52, 54, 56, 62, 64, 71, 73, 74, 75)]
CC_AUC_s$Group[1:61] <- "HC"
CC_AUC_s$Group[62:122] <- "NMDA"

# Permutation testing
nodes <- colnames(CC_AUC_s[-23])  

pvalsCC <- vector(length = 22)
i <- 1
for(x in nodes) {
  
  perm <- permutationTest2(CC_AUC_s, mean(eval(parse(text = x))), treatment = Group, seed = 0)
  
  pvalsCC[i] <- perm$stats$PValue
  
  i = i + 1
}


#### Modular Analysis ####

### AUC data

## Node Strength
hc_strAUC <- read.csv("~/R/StructNetwNMDA/data/conn_matrices/hc_strAUC_excl.csv", header=FALSE)
nmda_strAUC <- read.csv("~/R/StructNetwNMDA/data/conn_matrices/nmda_strAUC_excl.csv", header=FALSE)
str_AUC <- rbind(hc_strAUC, nmda_strAUC)
colnames(str_AUC) <- node_names$x

hcstr_MTL_AUC <- rowSums(hc_strAUC[,c(5,15,40,47,54,64)])
hcstr_DMN_AUC <- rowSums(hc_strAUC[,c(2, 3, 7, 13, 22, 24, 25, 26, 51, 52, 56, 62, 71, 73, 74, 75)])
nmdastr_MTL_AUC <- rowSums(nmda_strAUC[,c(5,15,40,47,54,64)])
nmdastr_DMN_AUC <- rowSums(nmda_strAUC[,c(2, 3, 7, 13, 22, 24, 25, 26, 51, 52, 56, 62, 71, 73, 74, 75)])

str_MTL_AUC <- as.data.frame(c(hcstr_MTL_AUC, nmdastr_MTL_AUC))
str_MTL_AUC$Group[1:61] <- 'HC'
str_MTL_AUC$Group[62:122] <- 'NMDA'
colnames(str_MTL_AUC) <- c('sumstr', 'Group')
str_DMN_AUC <- as.data.frame(c(hcstr_DMN_AUC, nmdastr_DMN_AUC))
colnames(str_DMN_AUC) <- c('sumstr', 'Group')
str_DMN_AUC$Group[1:61] <- 'HC'
str_DMN_AUC$Group[62:122] <- 'NMDA'

## Testing
permutationTest2(str_DMN_AUC, mean(sumstr), treatment = Group, seed = 0)


## Shortest paths
hcAPL_AUC <- read.delim("./data/gretna/ThesisExclMS/NodalShortestPath/Group1/aNLp.txt", header=FALSE)
nmdaAPL_AUC <- read.delim("./data/gretna/ThesisExclMS/NodalShortestPath/Group2/aNLp.txt", header=FALSE)

hcAPL_MTL_AUC <- rowMeans(hcAPL_AUC[,c(5,15,40,47,54,64)])
hcAPL_DMN_AUC <- rowMeans(hcAPL_AUC[,c(2, 3, 7, 13, 22, 24, 25, 26, 51, 52, 56, 62, 71, 73, 74, 75)])
nmdaAPL_MTL_AUC <- rowMeans(nmdaAPL_AUC[,c(5,15,40,47,54,64)])
nmdaAPL_DMN_AUC <- rowMeans(nmdaAPL_AUC[,c(2, 3, 7, 13, 22, 24, 25, 26, 51, 52, 56, 62, 71, 73, 74, 75)])

APL_MTL_AUC <- as.data.frame(c(hcAPL_MTL_AUC, nmdaAPL_MTL_AUC))
APL_MTL_AUC$Group[1:61] <- 'HC'
APL_MTL_AUC$Group[62:122] <- 'NMDA'
colnames(APL_MTL_AUC) <- c('mAPL', 'Group')
APL_DMN_AUC <- as.data.frame(c(hcAPL_DMN_AUC, nmdaAPL_DMN_AUC))
colnames(APL_DMN_AUC) <- c('mAPL', 'Group')
APL_DMN_AUC$Group[1:61] <- 'HC'
APL_DMN_AUC$Group[62:122] <- 'NMDA'

## Testing
permutationTest2(APL_MTL_AUC, mean(mAPL), treatment = Group, seed = 0)

## Plotting
mod_plot <- ggplot(apl, aes(y=mAPL, x=Group)) +
  geom_point(position=position_jitter(height=0, width=0.1)) +
  stat_summary(fun.y=mean, geom="point", color="red", size=3)
print(mod_plot)


#### Hub analysis ####

### Top-rank NS
hcstr_ranked <- as.data.frame(colMeans(hc_strAUC))
hcstr_ranked$node <- node_names$x
colnames(hcstr_ranked) <- c('meanStr', 'node')

nmdastr_ranked <- as.data.frame(colMeans(nmda_strAUC))
nmdastr_ranked$node <- node_names$x
colnames(nmdastr_ranked) <- c('meanStr', 'node')

### Top-rank BC
hcBC_ranked <- as.data.frame(colMeans(BC_AUC[1:61,]))
hcBC_ranked$node <- node_names$x
colnames(hcBC_ranked) <- c('meanBC', 'node')

hcBC_ranked$node <- factor(hcBC_ranked$node, levels = hcBC_ranked$node[order(hcBC_ranked$meanBC)])
p_hcBC_ranked <- ggplot(hcBC_ranked, aes(x=node, y=meanBC)) + geom_bar(stat = 'identity') + coord_flip()

nmdaBC_ranked <- as.data.frame(colMeans(BC_AUC[62:122,]))
nmdaBC_ranked$node <- node_names$x
colnames(nmdaBC_ranked) <- c('meanBC', 'node')

nmdaBC_ranked$node <- factor(nmdaBC_ranked$node, levels = nmdaBC_ranked$node[order(nmdaBC_ranked$meanBC)])
p_nmdaBC_ranked <- ggplot(nmdaBC_ranked, aes(x=node, y=meanBC)) + geom_bar(stat = 'identity') + coord_flip()

### Top-rank APL
# Get full data
hcAPL_AUC <- read.delim("./data/gretna/ThesisExclMS/NodalShortestPath/Group1/aNLp.txt", header=FALSE)
nmdaAPL_AUC <- read.delim("./data/gretna/ThesisExclMS/NodalShortestPath/Group2/aNLp.txt", header=FALSE)
APL_AUC <- rbind(hcAPL_AUC, nmdaAPL_AUC)
# get rid of empty 85th column
APL_AUC <- APL_AUC[1:84]
colnames(APL_AUC) <- node_names$x

# rank
hcAPL_ranked <- as.data.frame(colMeans(APL_AUC[1:61,]))
hcAPL_ranked$node <- node_names$x
colnames(hcAPL_ranked) <- c('meanAPL', 'node')

hcAPL_ranked$node <- factor(hcAPL_ranked$node, 
                            levels = hcAPL_ranked$node[order(hcAPL_ranked$meanAPL, decreasing = TRUE)])
p_hcAPL_ranked <- ggplot(drop_na(hcAPL_ranked), aes(x=node, y=meanAPL)) + 
  geom_bar(stat = 'identity') + coord_flip()

nmdaAPL_ranked <- as.data.frame(colMeans(APL_AUC[62:122,]))
nmdaAPL_ranked$node <- node_names$x
colnames(nmdaAPL_ranked) <- c('meanAPL', 'node')

nmdaAPL_ranked$node <- factor(nmdaAPL_ranked$node, 
                             levels = nmdaAPL_ranked$node[order(nmdaAPL_ranked$meanAPL, decreasing = TRUE)])
p_nmdaAPL_ranked <- ggplot(drop_na(nmdaAPL_ranked), aes(x=node, y=meanAPL)) + 
  geom_bar(stat = 'identity') + coord_flip()

### Top-rank PC
# Get full data
hcPC_AUC <- read.delim("./data/gretna/ThesisExclMS/CommunityIndex/Group1/DataDrivenPc_Thres010.txt", header=FALSE)
nmdaPC_AUC <- read.delim("./data/gretna/ThesisExclMS/CommunityIndex/Group2/DataDrivenPc_Thres010.txt", header=FALSE)
PC_AUC <- rbind(hcPC_AUC, nmdaPC_AUC)
# get rid of empty 85th column
PC_AUC <- PC_AUC[1:84]
colnames(PC_AUC) <- node_names$x

# rank
hcPC_ranked <- as.data.frame(colMeans(PC_AUC[1:61,]))
hcPC_ranked$node <- node_names$x
colnames(hcPC_ranked) <- c('meanPC', 'node')

hcPC_ranked$node <- factor(hcPC_ranked$node, 
                            levels = hcPC_ranked$node[order(hcPC_ranked$meanPC)])
p_hcPC_ranked <- ggplot(drop_na(hcPC_ranked), aes(x=node, y=meanPC)) + 
  geom_bar(stat = 'identity') + coord_flip()

nmdaPC_ranked <- as.data.frame(colMeans(PC_AUC[62:122,]))
nmdaPC_ranked$node <- node_names$x
colnames(nmdaPC_ranked) <- c('meanPC', 'node')

nmdaPC_ranked$node <- factor(nmdaPC_ranked$node, 
                              levels = nmdaPC_ranked$node[order(nmdaPC_ranked$meanPC)])
p_nmdaPC_ranked <- ggplot(drop_na(nmdaPC_ranked), aes(x=node, y=meanPC)) + 
  geom_bar(stat = 'identity') + coord_flip()

### Get hub score; set top percentage using slice(1:max)
hub_hcStr <- hcstr_ranked %>% arrange(desc(meanStr)) %>% slice(1:21)
hub_hcBC <- hcBC_ranked %>% arrange(desc(meanBC)) %>% slice(1:21)
hub_hcAPL <- hcAPL_ranked %>% arrange(meanAPL) %>% slice(1:21)
hub_hcPC <- hcPC_ranked %>% arrange(desc(meanPC)) %>% slice(1:21)

hub_hcall <- bind_cols(hub_hcStr, hub_hcBC, hub_hcAPL, hub_hcPC)

hub_hcscore <- node_names
hub_hcscore$hubscore <- rep(0,84)
colnames(hub_hcscore) <- c('node', 'hubscore')

for(node in 1:84) {
  
  score = 0
  # Str ranking
  if(any(as.character(node_names[node,1]) %in% hub_hcStr$node)) {
    score = score + 1
  }
  # BC ranking
  if(any(as.character(node_names[node,1]) %in% hub_hcBC$node)) {
    score = score + 1
  }
  # APL ranking
  if(any(as.character(node_names[node,1]) %in% hub_hcAPL$node)) {
    score = score + 1
  }
  # PC ranking
  if(any(as.character(node_names[node,1]) %in% hub_hcPC$node)) {
    score = score + 1
  }

  hub_hcscore$hubscore[node] <- score
}

print(hub_hcscore %>% arrange(desc(hubscore)))

# NMDA
hub_nmdaStr <- nmdastr_ranked %>% arrange(desc(meanStr)) %>% slice(1:21)
hub_nmdaBC <- nmdaBC_ranked %>% arrange(desc(meanBC)) %>% slice(1:21)
hub_nmdaAPL <- nmdaAPL_ranked %>% arrange(meanAPL) %>% slice(1:21)
hub_nmdaPC <- nmdaPC_ranked %>% arrange(desc(meanPC)) %>% slice(1:21)

hub_nmdaall <- bind_cols(hub_nmdaStr, hub_nmdaBC, hub_nmdaAPL, hub_nmdaPC)

hub_nmdascore <- node_names
hub_nmdascore$hubscore <- rep(0,84)
colnames(hub_nmdascore) <- c('node', 'hubscore')

for(node in 1:84) {
  
  score = 0
  # Str ranking
  if(any(as.character(node_names[node,1]) %in% hub_nmdaStr$node)) {
    score = score + 1
  }
  # BC ranking
  if(any(as.character(node_names[node,1]) %in% hub_nmdaBC$node)) {
    score = score + 1
  }
  # APL ranking
  if(any(as.character(node_names[node,1]) %in% hub_nmdaAPL$node)) {
    score = score + 1
  }
  # PC ranking
  if(any(as.character(node_names[node,1]) %in% hub_nmdaPC$node)) {
    score = score + 1
  }
  
  hub_nmdascore$hubscore[node] <- score
}

print(hub_nmdascore %>% arrange(desc(hubscore)))

# Rich club, feeder and peripheral connections
hubcon <- read.csv("~/R/StructNetwNMDA/data/hubanalysis_groups.csv")

permutationTest2(hubcon, mean(periph), treatment = group, seed = 0)
