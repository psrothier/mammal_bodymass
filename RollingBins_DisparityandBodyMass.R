#BPH check 8-6-24

##==========================================================================================================================================================================## 
##
## BODY MASS AND FORELIMB DISPARITY IN MAMMALS  ## -------------------------------------------------------------------------------------------------------------------------##
##
##==========================================================================================================================================================================##


### Question: is body mass driving trait and ecological disparity in mammals?

### Script organization:
  # Part I.  Data preparation: Size correction using phylogenetic regression. We will use 100 phylo trees from Upham 2019 to incorporate phylogenetic uncertainty.
  # Part II. Subseting data into body mass rolling bins
  # Part III. Overall forelimb disparity.
  # Part IV. Ecological and phylogeny diversity
        # evaluation of influence of body mass on morphological disparity and ecological and phylogenetic diversity
  # Part V. Subclade forelimb disparity
  # Part VI. Bone disparity



library(ape); library(geiger)
library(evobiR); library(phytools); library(reshape2)
library(ggplot2); library(dplyr); library(cowplot)
library(dispRity); library(mvMORPH); library(psych)
library(pals)



##-----------------------------------------------------------------------------------------------------------------##
## loading data ##
##-----------------------------------------------------------------------------------------------------------------##

setwd("C:/Users/Lab/OneDrive/Academia/Rothier et al/Body mass and disparity/final dataset")
#setwd("~/Desktop/Current Projects/Priscila-BodyMass/RCode-Priscila/RCodeJuly2024")

load("pruned_trees2024.RData") # loading 100 pruned trees my datlist 
load("mydat.list2024_3.RData") # loading our dataset, containing trait distances, body mass and ecological classification. rownames matching pruned trees tip.label order
#Note: Each dataset has a different order since each pruned tree has a different order. This allows incorporation of phylogenetic uncertainty

mydat.list <- mydat.list3

##==========================================================================================================================================================================## 
##
#### PART I - Data preparation ---- ####
##
##==========================================================================================================================================================================##


##-----------------------------------------------------------------------------------------------------------------##
## calculating geometric means ##
##-----------------------------------------------------------------------------------------------------------------##

# we need to correct our morphological distances for size, so then we can compare weather forelimb "shape" is correlated with mass
# we will use the geometric means of forelimb traits as a proxy for body size: (5 traits per bone)
lgm <- lapply(mydat.list, function(x) log10(apply(x[,c(9:28)], 1, geometric.mean))) 

# creating an object for all 100 datasets so that only includes forelimb measures
limb <- lapply(mydat.list, function(x) as.matrix(log10(x[,c(9:28)]))) 

#Further subset to remove digit length and mid-shaft measurements, maintaining only aspect ratio
limb <- lapply(limb, function(x) x[,c(1,2,4,6,7,9,11,12,14,16,17,19)]) 


## finally we can run a phylogenetic regression of size correction

data <- list()
fit_BM <- list() #because we are using 100 trees, perform 100 regressions of limb and geometric means and retain their residual values for the next analyses
for(i in 1:100){
  data[[i]] = list(limb=limb[[i]], lgm=lgm[[i]])
  fit_BM[[i]] <- mvgls(limb~lgm, data=data[[i]], pruned_trees[[i]], model="BM", penalty="LASSO", method="LL")
}

resid.list <- lapply(fit_BM, function(x) as.data.frame(x$resid)) #list of 100 residual objects



##-----------------------------------------------------------------------------------------------------------------##
## plot residuals x bodymass ##
##-----------------------------------------------------------------------------------------------------------------##

# let's visualize the relationship between size residual values of each trait against body mass to evaluate allometric effects

resid.sample <- resid.list[[1]] # I'm just using one residual example here, no need to repeat across the 100 regressions
rownames(resid.sample) == rownames(mydat.list[[1]])
resid.sample$Body_mass <- log10(mydat.list[[1]]$Body_Mass)
resid.sample$Group <- mydat.list[[1]]$Group
resid.sample$Group <- recode_factor(resid.sample$Group, Cetartiodactyla_Non_Cet="Cetartiodactyla")
unique(resid.sample$Group)
colsg <- setNames(c("#0000FF", "#FF0000", "#00FF00", "#000033", "#FF00B6", 
                    "#005300","lightgoldenrod1", "#009FFF","darkolivegreen", "#9A4D42",
                    "#00FFBE" ,"darkorchid4" ,"#1F9698", "#FFACFD", "#B1CC71","darkgoldenrod3", 
                    "#FE8F42", "mediumpurple1", "gray73"),
                  c("Monotremata","Eulipotyphla","Cetartiodactyla", "Carnivora","Perissodactyla",         
                    "Primates","Rodentia","Scandentia","Dermoptera","Lagomorpha",            
                    "Cingulata","Pilosa","Proboscidea","Hyracoidea","Tubulidentata","Pholidota",          
                    "Macroscelidea","Afrosoricida","Marsupialia"))


plots <- list()

for (i in 1:12) {
   
      plot <- eval(substitute(
      ggplot(resid.sample, aes(x = Body_mass, y = resid.sample[, i], color = Group))+
      geom_point(aes(color=Group), size=1.5, alpha=0.3) +
      scale_color_manual(values=colsg)+
      theme_classic()+
      theme(legend.position = "none")+
      ylab(colnames(resid.sample[i]))+
      xlab("Body mass"), list(i = i)))
  
      plots[[i]] <- plot
   
}

pg1 <- plot_grid(plotlist = plots, ncol=3, nrow=4)

leg <- ggplot(resid.sample, aes(x=Body_mass, y=resid.sample[,1], color=Group))+
  geom_point(aes(color=Group), size=3, alpha=0.5) +
  scale_color_manual(values=colsg)+
  theme_classic()+
  ylab(colnames(resid.sample[1]))+
  xlab("Body mass")

p_leg <- get_legend(leg)

plot_grid(pg1, p_leg, rel_widths=c(6,1))


#BPH Generally, I see a lot of phylogenetic structure in this example, but not a lot of allometry




##==========================================================================================================================================================================## 
##
#### PART II - Subseting data into body mass rolling bins ---- ####
##
##==========================================================================================================================================================================##


##=========================================================================================================================================##
#### attributing size bins ###
##=========================================================================================================================================##

## In order to calculate rolling mass bins, we will first separate the data into 20 FIXED BODY MASS bins


##-----------------------------------------------------------------------------------------------------------------##
#### a) Creating 20 *FIXED* mass bins  ###
##-----------------------------------------------------------------------------------------------------------------##

## attributing 20 bins based on crescent values of body masses
myascdat <- mydat.list[[1]] %>% arrange(Body_Mass) # rearranging df to be sorted by body mass
rownames(myascdat) <- myascdat$Species #Make rownames species

myascdat$bins20 <- as.numeric(cut_number(myascdat$Body_Mass, 20)) #arranging sample into 20 bins
myascdat %>% group_by(bins20) %>% summarise(length(bins20)) # species counts per bin - is varying from 32 (bin 18) to 34 specimens
 
## To make the species count more homogeneous, let's set the sample size variation to 33/34 species per bin
## To do so let's place the heaviest individual from bin 17 (now with 34 individuals) into bin 18, so both of then will have 33 individuals:

tail(subset(myascdat, bins20==17)) # let's move Vombatus ursinus (heaviest species of bin 17) to bin 18
myascdat[(rownames(myascdat)== "Vombatus_ursinus"),]$bins20 <- 18
myascdat %>% group_by(bins20) %>% summarise(length(bins20))




##-----------------------------------------------------------------------------------------------------------------##
#### b) Creating 20 *ROLLING* bins ###
##-----------------------------------------------------------------------------------------------------------------##

# now that we have our fixed 20 bins, let's make them overlap, creating the rolling bins
# We will create a list, in which each element within this list corresponds to a different bin
bin.list <- list() # first I'll create a list to each fixed bin
for(j in 1:20){
  bin.list[[j]] <- subset(myascdat, bins20== j)
  bin.list[[j]] <- data.frame(Species=bin.list[[j]]$Species, Group=bin.list[[j]]$Group, Body_Mass=bin.list[[j]]$Body_Mass, 
                              bins20=bin.list[[j]]$bins20, Locomotion= bin.list[[j]]$Locomotion, Family=bin.list[[j]]$Family)
}

# assignning species names as rownames
for(j in 1:20){
  rownames(bin.list[[j]]) <- bin.list[[j]]$Species
}
str(bin.list)

## now, I'll create the overlapping subsets:
# except for the first and last bins, at each intermediate bin, we will combine the species originally assigned to it plus 5/6 from prev bin and 5/6 from the next one, 
# creating an overlap of species between adjacent body mass bins, and a total of 44 species per rolling bin
# for the first bin (that lacks previous adjacent bin), we will take 10 individuals from bin 2 (the 10 heaviest), and for the last bin, we will take 10 heaviest individuals from the ante penultimate bin. 


## because the counts of each fixed bin are different, I'll custom bin per bin to make sure they all have the same number of individuals (44). 
overl.bins <- list()
overl.bins[[1]] = rbind(bin.list[[1]], bin.list[[2]][1:10,]) # 10 from bin 2
overl.bins[[2]] <-  rbind(bin.list[[2-1]][29:nrow(bin.list[[2-1]]),], bin.list[[2]], bin.list[[2+1]][1:5,])
overl.bins[[3]] <-  rbind(bin.list[[3-1]][28:nrow(bin.list[[3-1]]),], bin.list[[3]], bin.list[[3+1]][1:5,]) #16 from bin 1
overl.bins[[4]] <-  rbind(bin.list[[4-1]][29:nrow(bin.list[[4-1]]),], bin.list[[4]], bin.list[[4+1]][1:5,]) #11 from bin 2
overl.bins[[5]] <-  rbind(bin.list[[5-1]][29:nrow(bin.list[[5-1]]),], bin.list[[5]], bin.list[[5+1]][1:5,]) #10 from bin 3
overl.bins[[6]] <-  rbind(bin.list[[6-1]][28:nrow(bin.list[[6-1]]),], bin.list[[6]], bin.list[[6+1]][1:5,]) #11 from bin 4
overl.bins[[7]] <-  rbind(bin.list[[7-1]][28:nrow(bin.list[[7-1]]),], bin.list[[7]], bin.list[[7+1]][1:5,])
overl.bins[[8]] <-  rbind(bin.list[[8-1]][29:nrow(bin.list[[8-1]]),], bin.list[[8]], bin.list[[8+1]][1:5,])
overl.bins[[9]] <-  rbind(bin.list[[9-1]][29:nrow(bin.list[[9-1]]),], bin.list[[9]], bin.list[[9+1]][1:5,])
overl.bins[[10]] <-  rbind(bin.list[[10-1]][28:nrow(bin.list[[10-1]]),], bin.list[[10]], bin.list[[10+1]][1:5,])
overl.bins[[11]] <-  rbind(bin.list[[11-1]][28:nrow(bin.list[[11-1]]),], bin.list[[11]], bin.list[[11+1]][1:5,])
overl.bins[[12]] <-  rbind(bin.list[[12-1]][29:nrow(bin.list[[12-1]]),], bin.list[[12]], bin.list[[12+1]][1:5,])
overl.bins[[13]] <-  rbind(bin.list[[13-1]][29:nrow(bin.list[[13-1]]),], bin.list[[13]], bin.list[[13+1]][1:5,])
overl.bins[[14]] <-  rbind(bin.list[[14-1]][28:nrow(bin.list[[14-1]]),], bin.list[[14]], bin.list[[14+1]][1:5,])
overl.bins[[15]] <-  rbind(bin.list[[15-1]][28:nrow(bin.list[[15-1]]),], bin.list[[15]], bin.list[[15+1]][1:5,])
overl.bins[[16]] <-  rbind(bin.list[[16-1]][29:nrow(bin.list[[16-1]]),], bin.list[[16]], bin.list[[16+1]][1:5,])
overl.bins[[17]] <-  rbind(bin.list[[17-1]][29:nrow(bin.list[[17-1]]),], bin.list[[17]], bin.list[[17+1]][1:5,])
overl.bins[[18]] <-  rbind(bin.list[[18-1]][28:nrow(bin.list[[18-1]]),], bin.list[[18]], bin.list[[18+1]][1:5,])
overl.bins[[19]] <-  rbind(bin.list[[19-1]][28:nrow(bin.list[[19-1]]),], bin.list[[19]], bin.list[[19+1]][1:5,])
overl.bins[[20]] = rbind(bin.list[[19]][24:nrow(bin.list[[19]]),],bin.list[[20]])

str(overl.bins)


for(k in 1:20){
  overl.bins[[k]]$roll.bins <- k
    }

dim(overl.bins[[19]])
overl.bins[[19]] %>% group_by(bins20)%>%summarise(length(bins20)) #checking counts of a few rolling bins
overl.bins[[10]] %>% group_by(bins20)%>%summarise(length(bins20))
overl.bins[[20]] %>% group_by(bins20)%>%summarise(length(bins20))

str(overl.bins) # all looking good


## next - create a list of lists: outer is each fitted model across 100 trees, inside is their corresponding residuals,organized into a list of rolling bins
# inside = match residual rownames with overl.bins bin row names

resid.overl.bins <- list()
for (i in 1:length(resid.list)) {
  resid.overl.bins[[i]] <- list()  # Initialize a sublist for each dataframe in resid.list
  
  # Iterate over each dataframe in overl.bins
  for (j in 1:length(overl.bins)) {
    resid.overl.bins[[i]][[j]] <- as.matrix(resid.list[[i]][rownames(overl.bins[[j]]), ])
  }
}

#checking random subsets to see if name sequence is really matching 
row.names(overl.bins[[6]])== row.names(resid.overl.bins[[23]][[6]]) # ok





##==========================================================================================================================================================================## 
##
#### PART III -  Overall forelimb disparity ---- ####
##
##==========================================================================================================================================================================##

## all good, now we can run disparity within our rolling bins, repeating the analyses 100 times to cover the residual regressions from 100 trees
# we will calculate disparity using dispRity package, which allows data to be bootstrapped to improve reliability of results

library(dispRity)


# 1) bootstrap data
boot_resid <- list()
for(i in 1:100){
boot_resid[[i]] <- lapply(resid.overl.bins[[i]], function(x) boot.matrix(x)) # default is 100 replications (*100 trees = 10,000 replications per bin)
}

# 2) calculate disparity from bootstrapped data
disp_resid <- list()
for(i in 1:100){
disp_resid[[i]] <- lapply(boot_resid[[i]], function(x) dispRity(x, metric = c(sum, variances), verbose=TRUE))
} # results of disparity are good to be visualized


# 3) extracting results to plot
mat.disp <- matrix(ncol=20, nrow=100) # 20 bins, 100 replications (per tree)
list.disp <- list() # creating a list of 100 matrices, corresponding to the disparity values obtained from the limb shape residuals  over the 100 trees
for(k in 1:100){
  list.disp[[k]] <- list()
  for(i in 1:20){
    mat.disp[,i] <- disp_resid[[k]][[i]]$disparity[[1]][[2]] # taking the disparity from each bootstrap rep
    colnames(mat.disp) <- 1:20
    list.disp[[k]] <-   mat.disp
  }
}

disp.df <- do.call(rbind, list.disp) # transforming disparity list in a dataframe (columns are the rolling bins (n = 20), rows are the bootstrap rep (n = 10000))
str(disp.df)
disp.df[701:800,7]==list.disp[[8]][,7] # comparing if dataframe rows correspond to the list values - in this case, disparity of rolling bin 7 from phyl residuals using tree 8

disp.df.plot <- melt(disp.df) # using melt function to reorganize dataframe to ggplot
disp.df.plot %>% group_by(Var2) %>% count(Var2) # each rolling bin has 10,000 disparity values, which is simply 100 bootstrap reps for 100 trees!
disp.df.plot$value[1:100]==list.disp[[1]][,1] # checking if first 100 elements of our new dataframe correspond to the disparity of rolling bin 1 from tree 1, and it does
disp.df.plot$value[101:200]==list.disp[[2]][,1] # checking a random result again, and all is good!!

disp.df.plot <- disp.df.plot[,-1] #removing useless row
colnames(disp.df.plot) <- c("rolling.bin", "disparity")


# 4) plots
#plot 1 of rolling bins showing limb shape disparity

ggplot(disp.df.plot, aes(x=rolling.bin, y=disparity, group=1)) +
  theme_minimal_vgrid()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "steelblue", fill="steelblue", alpha=0.3)+
  xlab("Rolling bins") + ylab("Forelimb shape disparity")+
  theme(axis.title.x = element_text(size=14, vjust = -1),
        axis.title.y = element_text(color = "steelblue", size=14,  vjust = 1, hjust = 1),
        axis.line.y = element_line(color = "steelblue"),
        axis.ticks.y = element_line(color = "steelblue"),
        axis.text.y = element_text(color = "steelblue"),
        axis.text = element_text(size=12))


# transforming rolling bin axis into log10 transformed body mass mean values

med.rb1 <- sapply(overl.bins, function(x) median((x$Body_Mass)))
med.rb <- format(round(med.rb1/1000, 2)) # round values
logmed.rb <- sapply(overl.bins, function(x) median(log10(x$Body_Mass)))
logmed.rb <- format(round(logmed.rb,2), nsmall=2)


logmed.rb.vec <- list() # combining these values into the plot dataframe
for(i in 1:20){
  logmed.rb.vec[[i]] <- rep(logmed.rb [i], 10000)
}
logmed.rb.vec <- unlist(logmed.rb.vec)
disp.df.plot$log.med.rb <- as.numeric(logmed.rb.vec)


kg.med.rb <- list() # same thing
for(i in 1:20){
  kg.med.rb[[i]] <- rep(med.rb [i], 10000)
}
kg.med.rb <- unlist(kg.med.rb)
disp.df.plot$kg.med.rb <- kg.med.rb


# plot scaling rolling bins x axis to their median body mass (in log)
ggplot(disp.df.plot, aes(x=log.med.rb, y=disparity, group=1)) +
  theme_bw()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "steelblue", fill="steelblue", alpha=0.3)+
  xlab("log10 body mass (rolling bin median)") + ylab("Forelimb shape disparity")+
  theme(axis.title.x = element_text(size=14, vjust = -1),
        axis.title.y = element_text(color = "steelblue", size=14,  vjust = 1, hjust = 1),
        axis.line.y = element_line(color = "steelblue"),
        axis.ticks.y = element_line(color = "steelblue"),
        axis.text.y = element_text(color = "steelblue"),
        axis.text = element_text(size=12))


#  plot scaling rolling bins x axis to their median body mass (raw gram values)
ggplot(disp.df.plot, aes(x=kg.med.rb, y=disparity, group=1)) +
  theme_minimal_vgrid()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "steelblue", fill="steelblue", alpha=0.3)+
  xlab("body mass (kg)") + ylab("Forelimb shape disparity")+
  theme(axis.title.x = element_text(size=14, vjust = -1),
        axis.title.y = element_text(color = "steelblue", size=14,  vjust = 1, hjust = 1),
        axis.line.y = element_line(color = "steelblue"),
        axis.ticks.y = element_line(color = "steelblue"),
        axis.text.y = element_text(color = "steelblue"),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size=12))






##==========================================================================================================================================================================## 
##
#### Part IV - Ecological and phylogeny diversity ---- ####
##
##==========================================================================================================================================================================##



##=========================================================================================================================================##
#### 1. Phylogenetic disparity ####
##=========================================================================================================================================##

# now we will evaluate if phylogenetic diversity changes along the body mass axis
# we will calculate the phylogenetic disparity across our rolling bins using  our 100 trees, using Faith's PD (from EcoPhyloMapper package)
library(epm)

#Generating Faith's PD for all 20 of our bins for 100 trees
faithPD20 <- list()
for(i in 1:100){
  faithPD20[[i]] <- lapply(resid.overl.bins[[i]], function(x) faithPD(pruned_trees[[i]], rownames(x)))
}

faithPD20[[1]] # Faith PD at each rolling bin, using tree 1

#Transforem this list into a dataframe
faithPD20.df <- melt(faithPD20)
colnames(faithPD20.df) <- c("faith.dist", "bin", "tree")

faithPD20.df$mass <- as.numeric(rep(logmed.rb, n=100))
namebins <- 1:20 #naming the rolling bins
seq <- faithPD20.df %>% dplyr::filter(bin %in% namebins) #putting factor in sequence so the plot is not out of order
faithPD20.df$bin <- factor(seq$bin, levels= namebins)

faithPD20.df$kg.mass <- as.factor(rep(med.rb, n=100))# including body mass values in the df so we can plot it later

#Plot by bin in x axis
ggplot(faithPD20.df, aes(x=bin, y=faith.dist, group=1)) +
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE)+
               theme_classic() +
               xlab("size bins") + ylab("Phylogenetic Diversity")

#Plot by log10 mass in x axis
ggplot(faithPD20.df, aes(x=kg.mass, y=faith.dist, group=1)) +
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE)+
               theme_classic() +
               xlab("log10 (Body Mass)") + ylab("Phylogenetic Diversity")





##=========================================================================================================================================##
#### 2. Ecological disparity ####
##=========================================================================================================================================##

# now we will calculate ecological disparity using Gower's index. Here, we are obtaining just one single value of Gower per rolling bin, derived from the dissimilarity of ecological categories between species
# Therefore, these results are not dependent on the phylogenetic history, so they do not need to be calculated 100 times
# Note that each taxon is given a single category

library(cluster) # Gower diversity will be calculated with the function daisy from the cluster package



###-----------------------------------------------------------------------------------------------------------------##
#### 2a) binary matrix of ecological categories ####
##------------------------------------------------------------------------------------------------------------------##


# to make the function work, we need to transform our ecological variables (now as factors) into a binary matrix
df <- data.frame(matrix(ncol=8, nrow=666)) # 8 ecological categories, 666 species
loc <- unique(mydat.list[[1]]$Locomotion) #taking our locomotor habits

loc_list <- list()
for(i in 1:length(loc)){
  df[,i] <- as.factor(ifelse(mydat.list[[1]]$Locomotion==loc[i], 1, 0)) # creating a matrix of 0 and 1, assigning the locomotor modes (1) to each species 
}

colnames(df) <- loc
rownames(df) <- rownames(mydat.list[[1]])


# organize ecological classification into a list of mass bins 
list.rb <- list()
for(i in 1:20){
  list.rb[[i]] <- df[rownames(overl.bins[[i]]),]
}



###-----------------------------------------------------------------------------------------------------------------##
#### 2b) gower dissimilarity ####
##------------------------------------------------------------------------------------------------------------------##

# calculate the pairwise gower dissimilarity between species at each bin 
gower.dist <- lapply(list.rb, function(x) daisy(x[,1:8], metric="gower"))

# creating a df with the average value of Gower dissimilarity
gower.mean <- lapply(gower.dist, function(x) summary(x)) 
gower.means <- unlist(lapply(gower.mean, function(x) x$summ[4]))
df.gower.means <- as.data.frame(cbind(bins=1:20, gower.means=gower.means, logmed.rb=as.numeric(logmed.rb)))
df.gower.means$kg.mass <- as.factor(med.rb)
df.gower.means$kg.mass <- factor(df.gower.means$kg.mass, levels=df.gower.means$kg.mass[order(df.gower.means$bins)])

ggplot(df.gower.means, aes(x=bins, y=gower.means)) +
  geom_line() +
  theme_classic()

ggplot(df.gower.means, aes(x=kg.mass, y=gower.means, group=1)) +
  geom_line() +
  theme_classic()


##=========================================================================================================================================##
#### 3. Plotting combined results ####
##=========================================================================================================================================##

#### ** FIGURE 1 ** ####

## 1 - all traits + phyl diversoty and gower dist
pp1.rb <- ggplot(disp.df.plot, aes(x=kg.med.rb, y=disparity*10, group=1)) +
  theme_minimal_vgrid()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "indianred1", fill="indianred1", alpha=0.3)+
  xlab("") + ylab("Forelimb disparity")+
  scale_y_continuous(limits = c(0.2, 2.3), breaks=seq(0.5, 2.0, 0.5))+
  #scale_x_continuous(limits = c(1, 20), breaks=seq(1,20,1))+
  theme(axis.title.x = element_text(size=11, vjust = -1,),
        axis.title.y = element_text(color = "indianred1", size=11,  vjust = 1, hjust = 0),
        axis.line.y = element_line(color = "indianred1"),
        axis.ticks.y = element_line(color = "indianred1"),
        axis.text.y = element_text(color = "gray33"),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size=8))



pp2.rb <- ggplot(df.gower.means, aes(x=kg.mass, y=gower.means*10, group=1)) +
  theme_minimal_vgrid()+ 
  geom_line(linetype=2, size=1, color="gray33")+
  xlab("") + ylab("Ecological diversity")+
  #scale_x_continuous(breaks=seq(1, 20, 1))+
  scale_y_continuous(limits = c(0.2, 2.3), breaks=seq(0.5, 2.0, 0.5))+
  #scale_x_continuous(limits = c(1, 20), breaks=seq(1,20,1))+
  theme(axis.title.x = element_text(size=11, vjust = -1),
        axis.title.y = element_text(color = "gray33", size=11,  vjust = 1, hjust = 1),
        axis.line.y = element_line(color = "indianred1"),
        axis.ticks.y = element_line(color = "indianred1"),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(color = "gray33"),
        axis.text = element_text(size=8))


faithPD20.df$bin <- as.numeric(faithPD20.df$bin)

pp3.rb <- ggplot(faithPD20.df, aes(x=kg.mass, y=faith.dist, group=1)) +
  theme_minimal_vgrid()+ 
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "slategray3", fill="slategray3", alpha=0.3)+
  xlab("") + ylab("Faith's Phylogenetic Diversity")+
  scale_y_continuous(position = "right") +
 # scale_x_continuous(limits = c(1, 20), breaks=seq(1,20,1))+
  theme(axis.title.x = element_text(size=11, vjust = -1),
        axis.title.y = element_text(color = "slategray3", size=11,  hjust = 1, vjust=1),
        axis.line.y = element_line(color = "slategray3"),
        axis.ticks.y = element_line(color = "slategray3"),
        axis.text.y = element_text(color = "slategray3"),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size=8))




aligned_plots.rb <- align_plots(pp1.rb, pp2.rb, pp3.rb, align="hv", axis="tblr") 


main_plot <- ggdraw(aligned_plots.rb[[2]]) + draw_plot(aligned_plots.rb[[3]]) + draw_plot(aligned_plots.rb[[1]])



## 2 - bar plots


##  bar clades --------------- 

# Make a bar plot for each clade to show what clades happen in those bins

library(viridis)
length(unique(myascdat$Group))
vir <- (viridis_pal()(19))
#replacing the yellow color for a darker tone, for better visualization
vir <- recode(vir, '#FDE725FF'="#FFD400")

clades <- c("Monotremata", "Marsupialia", "Pilosa", "Cingulata", "Proboscidea","Hyracoidea", "Tubulidentata",
            "Macroscelidea","Afrosoricida" , "Eulipotyphla", "Pholidota", "Carnivora","Perissodactyla",
            "Cetartiodactyla_Non_Cet", "Dermoptera", "Scandentia", "Primates", "Lagomorpha","Rodentia")


colsg <- setNames(vir, clades)

clade_bin_list <- lapply(overl.bins, function(x) x$Group)
df_clade_bin_list <- melt(clade_bin_list)

p_clade_bar <- ggplot(df_clade_bin_list, aes(x=L1)) + 
  geom_bar(aes(fill=value), width = 0.8)+
  scale_fill_manual(values=colsg)+
  theme_classic()+
  ylab("Clades")+
  xlab("")+
  scale_y_continuous(limits = c(0, 44), breaks=seq(0,40,20))+
  theme(axis.title.y = element_text(size=11), 
        axis.text = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y =element_blank(),
        legend.position = "none")+
  scale_x_continuous(breaks=seq(1, 20, 1))+
  easy_remove_axes("x")

leg_clade <- get_legend(ggplot(df_clade_bin_list, aes(x=L1)) + 
                          geom_bar(aes(fill=value), width = 0.8)+
                          scale_fill_manual(values=colsg)+
                          theme_classic()+
                          ylab("Clade ")+
                          xlab("")+
                          theme(axis.title.x = element_text(size=11, vjust = -1), 
                                axis.title.y = element_text(size=11), 
                                axis.text = element_text(size=8),
                                legend.text=element_text(size=10),
                                legend.title =element_text(size=11))+
                          scale_x_continuous(breaks=seq(1, 20, 1))+
                          easy_remove_axes("x"))





## bar plot ecology --------------- 

cols.loc<-setNames(c("forestgreen", "darkred", "black","olivedrab1", "skyblue", "red", "violet", "#FFC055"),
                   c("arboreal","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad"))

eco_bin_list <- lapply(overl.bins, function(x) x$Locomotion)
df_eco_bin_list <- melt(eco_bin_list) # taking the ecology diversity at each bin

#Make a bar plot for each bin to show what ecological modes happen in those bins
library(ggeasy)


p_eco_bar <- ggplot(df_eco_bin_list, aes(x=L1)) + 
  geom_bar(aes(fill=value),width=0.8)+
  scale_fill_manual(values=cols.loc)+
  theme_classic()+
  ylab("Locomotion")+
  xlab("Mass rolling bin")+
  scale_y_continuous(limits = c(0, 44), breaks=seq(0,40,20))+
  theme(axis.title.x = element_text(size=11, vjust = -1), 
        axis.title.y = element_text(size=11),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y =element_blank(),
        legend.position="none")+
  scale_x_continuous(breaks=seq(1, 20, 1))




leg_eco <- get_legend(ggplot(df_eco_bin_list, aes(x=L1)) + 
                        geom_bar(aes(fill=value),width=0.8)+
                        scale_fill_manual(values=cols.loc)+
                        theme_classic()+
                        ylab("Ecology ")+
                        xlab("Mass roling bins")+
                        theme(axis.title.x = element_text(size=11, vjust = -1), 
                              axis.title.y = element_text(size=11), 
                              axis.text = element_text(size=11),
                              legend.text=element_text(size=10),
                              legend.title =element_text(size=8))+
                        scale_x_continuous(breaks=seq(1, 20, 1)))






## Combining into ultimate Fig 1 plot --------------
bar_plots <- plot_grid( p_clade_bar, p_eco_bar, ncol=1,  align="hv")
fig1 <- plot_grid(main_plot, bar_plots, ncol=1, rel_heights = c(1.5,1),  align="hv")
leg_plot <- plot_grid(leg_clade, leg_eco, ncol=1, align="hv", scale=0.8)
fig1




##=========================================================================================================================================##
#### 4. Determining correlation ####
##=========================================================================================================================================##

# lets evaluate the correlation between disparity, body mass, phyl diversity and ecol dissimilarity
# this will be determined using parametric and non parametric tests

library(ggcorrplot)
library(nortest)

## to make the datasets comparable, we will use the mean values of the calculated variables per means 
#(except for body mass, cause median seems more appropriate to represent diversity)

# disparity
disp.means <- as.data.frame(disp.df.plot %>% group_by(rolling.bin) %>% summarise(mean(disparity)))

# faith distance
faith.means <- as.data.frame(faithPD20.df %>% group_by(bin) %>% summarise(mean(faith.dist)))

# body mass
median.bodymass <- as.numeric(logmed.rb)

# gathering all into one df:
df.mlr <- data.frame(rolling=c(1:20), median.bodymass=median.bodymass, disparity=disp.means$`mean(disparity)`, faith=faith.means$`mean(faith.dist)`, gower= df.gower.means$gower.means)



###-----------------------------------------------------------------------------------------------------------------##
#### 4a)  multiple linear model - parametric test ####
##------------------------------------------------------------------------------------------------------------------##

# checking normality to run linear model
disp.df.plot
lillie.test(df.mlr$median.bodymass)
lillie.test(df.mlr$disparity)
lillie.test(df.mlr$faith)
lillie.test(df.mlr$gower) # gower is not normal.. so its good to run a non-parametric test as well

qqnorm(df.mlr$median.bodymass)
qqline(df.mlr$median.bodymass)
qqnorm(df.mlr$disparity)
qqline(df.mlr$disparity)
qqnorm(df.mlr$faith)
qqline(df.mlr$faith)
qqnorm(df.mlr$gower)
qqline(df.mlr$gower)

# running the linear regression
fit.mlr <- lm(disparity~ median.bodymass + faith + gower, data=df.mlr)
summary(fit.mlr)
# > conclusion: disparity is strongly correlated with body mass, bot not with other variables!



## with all disparity results now
fit.disp.bm <- lm(disparity ~ log.med.rb, data=disp.df.plot)
summary(fit.disp.bm)

fit.disp.bm1 <- lm(disparity ~ log.med.rb, data=disp.df.plot[1:100000,])
summary(fit.disp.bm1)

fit.disp.bm2 <- lm(disparity ~ log.med.rb, data=disp.df.plot[100001:200000,])
summary(fit.disp.bm2)

###-----------------------------------------------------------------------------------------------------------------##
#### 4b)  Kendall test - non-parametric test ####
##------------------------------------------------------------------------------------------------------------------##

# we will run non-parametric correlation tests first using (1) all bins and then separating the (2) small from (3) large mammals 
# (so we can understand if  body mass has a different effect on disparty and ecologi in small vs larger taxa)


# exploratory overview: checking on correlation across variables:
corr_matrix = round(cor(df.mlr[,2:5]), 2)


## 4b.1 All bins ------------------------------------------

#Run tests comparing across all comparisons
k.AllBins.mbm <-cor.test(df.mlr$disparity, df.mlr$median.bodymass, method="kendall") # correlation using all bins: morphological disparity and body mass
k.AllBins.me <-cor.test(df.mlr$disparity, df.mlr$gower, method="kendall") # disparity and ecology
k.AllBins.mp <-cor.test(df.mlr$disparity,df.mlr$faith, method="kendall") # disparity and phylo div
k.AllBins.mbe <-cor.test(df.mlr$median.bodymass, df.mlr$gower, method="kendall") # mass and ecology
k.AllBins.mbp <-cor.test(df.mlr$median.bodymass,df.mlr$faith, method="kendall") # mass and phylo div
k.AllBins.pe <-cor.test(df.mlr$faith, df.mlr$gower, method="kendall") # phylo div and ecology

#make a list
list_kendall_AllBins <- list(k.AllBins.mbm, k.AllBins.me, k.AllBins.mp, 
                             k.AllBins.mbe , k.AllBins.mbp, k.AllBins.pe)

# extracting tau and p values to plot
list_kendall_AllBins
tau <- sapply(list_kendall_AllBins, function(x) x$estimate) 
pvalue <- sapply(list_kendall_AllBins, function(x) x$p.value)

#creating a matrix of tau  between parameters
#Row/column 1 = disparity, row/col 2 = body mass, row/col 3= ecology, row/col 4 = phylo
mat.tau <- matrix(ncol=4, nrow=4, 
                  c(1,           tau[1],     tau[2],     tau[3],   # disparity
                    tau[1],       1,         tau[4],     tau[5], # bodymass  
                    tau[2],      tau[4],       1,        tau[6], # ecology
                    tau[3],      tau[5],      tau[6],      1)) # faith dist

nam.mat <- c("Disparity", "Body mass", "Ecological diversity", "Faith's distance")
colnames(mat.tau) <- nam.mat
rownames(mat.tau) <- nam.mat


# creating a matrix of p-values  between parameters, so we will only plot significant taus
mat.p <- matrix(ncol=4, nrow=4,
                  c(1,           pvalue[1],     pvalue[2],     pvalue[3],   # disparity
                    pvalue[1],       1,         pvalue[4],     pvalue[5], # bodymass  
                    pvalue[2],      pvalue[4],       1,        pvalue[6], # ecology
                    pvalue[3],      pvalue[5],      pvalue[6],      1))# faith dist

colnames(mat.p) <- nam.mat
rownames(mat.p) <- nam.mat

#Visualize what is significant with boxes colored by correlation
pcor1 <- ggcorrplot(mat.tau, p.mat=mat.p, hc.order = TRUE, lab=TRUE, type = "lower", insig="blank", legend.title = "Correlation (p < 0.05)", 
           title="All rolling bins", lab_size=3, tl.cex = 10 )



## now let's separate bins into lower (first 10) and higher (last 10) masses to see if we can detect different patterns across smaller or larger mammals




## 4b.2 Lower bins (small to medium mammals) ------------------------------------------


# Parametric analysis
fit.mlr_1_10 <- lm(disparity~ median.bodymass + faith + gower, data=df.mlr[1:10,])
summary(fit.mlr_1_10 )

# Nonparametric analysis
k.1_10.mbm <-cor.test(df.mlr[1:10,]$disparity, df.mlr[1:10,]$median.bodymass, method="kendall") 
k.1_10.me <-cor.test(df.mlr[1:10,]$disparity, df.mlr[1:10,]$gower, method="kendall") 
k.1_10.mp <-cor.test(df.mlr[1:10,]$disparity,df.mlr[1:10,]$faith, method="kendall") 
k.1_10.mbe <-cor.test(df.mlr[1:10,]$median.bodymass, df.mlr[1:10,]$gower, method="kendall") 
k.1_10.mbp <-cor.test(df.mlr[1:10,]$median.bodymass,df.mlr[1:10,]$faith, method="kendall") 
k.1_10.pe <-cor.test(df.mlr[1:10,]$faith, df.mlr[1:10,]$gower, method="kendall") 

list_kendall_1_10 <- list(k.1_10.mbm, k.1_10.me, k.1_10.mp, 
                          k.1_10.mbe , k.1_10.mbp, k.1_10.pe)

#Generate a matrix as above and plot
tau_1_10 <- sapply(list_kendall_1_10, function(x) x$estimate)
pvalue_1_10 <- sapply(list_kendall_1_10, function(x) x$p.value)

mat.tau_1_10 <- matrix(ncol=4, nrow=4,
                  c(1,                tau_1_10[1],      tau_1_10[2],     tau_1_10[3],   # disparity
                    tau_1_10[1],       1,              tau_1_10[4],     tau_1_10[5], # bodymass  
                    tau_1_10[2],      tau_1_10[4],       1,             tau_1_10[6], # ecology
                    tau_1_10[3],      tau_1_10[5],      tau_1_10[6],      1))# faith dist

colnames(mat.tau_1_10) <- nam.mat
rownames(mat.tau_1_10) <- nam.mat

# plotting only significant taus
mat.p_1_10 <- matrix(ncol=4, nrow=4,
                c(1,                   pvalue_1_10[1],     pvalue_1_10[2],     pvalue_1_10[3],   # disparity
                  pvalue_1_10[1],       1,                  pvalue_1_10[4],   pvalue_1_10[5], # bodymass  
                  pvalue_1_10[2],      pvalue_1_10[4],       1,                pvalue_1_10[6], # ecology
                  pvalue_1_10[3],      pvalue_1_10[5],      pvalue_1_10[6],      1))# faith dist

colnames(mat.p_1_10) <- nam.mat
rownames(mat.p_1_10) <- nam.mat

pcor2 <- ggcorrplot(mat.tau_1_10, p.mat=mat.p_1_10, lab=TRUE, type = "lower", insig="blank", legend.title = "Correlation (p < 0.05)", 
           title="Small mammals (1-10 RB)", lab_size=3, tl.cex = 10 )



## 4b.3 Higher bins (medium to large mammals) ------------------------------------------

# Parametric analysis
fit.mlr_11_20 <- lm(disparity~ median.bodymass + faith + gower, data=df.mlr[11:20,])
summary(fit.mlr_11_20 )

# Nonparametric analysis
k.11_20.mbm <-cor.test(df.mlr[11:20,]$disparity, df.mlr[11:20,]$median.bodymass, method="kendall") 
k.11_20.me <-cor.test(df.mlr[11:20,]$disparity, df.mlr[11:20,]$gower, method="kendall") 
k.11_20.mp <-cor.test(df.mlr[11:20,]$disparity,df.mlr[11:20,]$faith, method="kendall") 
k.11_20.mbe <-cor.test(df.mlr[11:20,]$median.bodymass, df.mlr[11:20,]$gower, method="kendall") 
k.11_20.mbp <-cor.test(df.mlr[11:20,]$median.bodymass,df.mlr[11:20,]$faith, method="kendall") 
k.11_20.pe <-cor.test(df.mlr[11:20,]$faith, df.mlr[11:20,]$gower, method="kendall") 

list_kendall_11_20 <- list(k.11_20.mbm, k.11_20.me, k.11_20.mp, 
                          k.11_20.mbe , k.11_20.mbp, k.11_20.pe)

tau_11_20 <- sapply(list_kendall_11_20, function(x) x$estimate)
pvalue_11_20 <- sapply(list_kendall_11_20, function(x) x$p.value)

mat.tau_11_20 <- matrix(ncol=4, nrow=4,
                       c(1,                tau_11_20[1],      tau_11_20[2],     tau_11_20[3],   # disparity
                         tau_11_20[1],       1,              tau_11_20[4],     tau_11_20[5], # bodymass  
                         tau_11_20[2],      tau_11_20[4],       1,             tau_11_20[6], # ecology
                         tau_11_20[3],      tau_11_20[5],      tau_11_20[6],      1))# faith dist

colnames(mat.tau_11_20) <- nam.mat
rownames(mat.tau_11_20) <- nam.mat

mat.p_11_20 <- matrix(ncol=4, nrow=4,
                     c(1,                   pvalue_11_20[1],     pvalue_11_20[2],     pvalue_11_20[3],   # disparity
                       pvalue_11_20[1],       1,                  pvalue_11_20[4],   pvalue_11_20[5], # bodymass  
                       pvalue_11_20[2],      pvalue_11_20[4],       1,                pvalue_11_20[6], # ecology
                       pvalue_11_20[3],      pvalue_11_20[5],      pvalue_11_20[6],      1))# faith dist

colnames(mat.p_11_20) <- nam.mat
rownames(mat.p_11_20) <- nam.mat

pcor3 <- ggcorrplot(mat.tau_11_20, p.mat=mat.p_11_20, lab=TRUE, type = "lower", insig="blank", legend.title = "Correlation (p < 0.05)", 
           title="Large mammals (11-20 RB)", lab_size=3, tl.cex = 10 )



## 4b.4 Plotting results ------------------------------------------

plot_grid(pcor1,pcor2,pcor3,nrow=1)



##==========================================================================================================================================================================## 
##
#### PART V - Subclade forelimb disparity ----  ####
##
##==========================================================================================================================================================================##

# the idea here is to repeat the analyses previously described, but across mammal major clades
# after selecting our clades, we will prune the clade trees and calculate their size residuals. 
# Once the residuals are calculated, we will arrange the species into rolling bins and calculate the bin disparities


##------------------------------------------------------------------------------------------------------------------##
## pruning trees and size bins ##
##------------------------------------------------------------------------------------------------------------------##

clade.list <- lapply(mydat.list, function(x) recode(x$Group, Cetartiodactyla_Non_Cet="Cetartiodactyla")) # first, let's give a better name to Cetartiodactyls

for(i in 1:(length(mydat.list))){
  mydat.list[[i]]$Group <- clade.list[[i]]
}

# and let's count how many spcies we have for each major mammal clade. In order to run some robust analyses, we will ideally target a minimum sample size of at least 35 species
mydat.list[[1]] %>% group_by(Group) %>% summarise(length(Group)) # our candidate subclades are : Carnivora 89; Cetartiodactyla 86, Eulipotyphla 39, Marsupialia 63, Rodentia 254, Primates 66

clades <- c("Marsupialia","Rodentia", "Primates", "Carnivora", "Cetartiodactyla", "Eulipotyphla" )


## NOTE: the criteria we will use now is to create fixed bins with at least 10 individuals in each bin, creating a rolling bin of size 16. only exception is eulipothyphlans (small sample size) and rodents (large sample size)
        # Therefore, while the number of species within each rolling bin will be the same across our subclades (16, but, again, except for Rodents and Eulipothyphlans), 
        # the number of rolling bins per clade will be variable according to the total sample size of each group


##=========================================================================================================================================##
#### a) Marsupialia ####
##=========================================================================================================================================##

# subset data to only include marsupials
mydat.Marsupialia <- list()
for(i in 1:(length(mydat.list))){
  mydat.Marsupialia[[i]] <- subset(mydat.list[[i]], Group=="Marsupialia")
}

# prune tree to only include marsupials
pruned.trees.Marsupialia <- list()
for( i in 1: length(pruned_trees)){
  pruned.trees.Marsupialia[[i]] <- drop.tip(pruned_trees[[i]], (name.check(pruned_trees[[i]], mydat.Marsupialia[[i]])$tree_not_data))
  }

for(i in 1:(length(mydat.Marsupialia))){
  mydat.Marsupialia[[i]] <- mydat.Marsupialia[[i]][pruned.trees.Marsupialia[[i]]$tip.label,]
}

# body size correction
lgm.Marsupialia <- lapply(mydat.Marsupialia, function(x) log10(apply(x[,c(9:28)], 1, geometric.mean))) # using limb traits only to calculate geometric means
limb.Marsupialia <- lapply(mydat.Marsupialia, function(x) as.matrix(log10(x[,c(9:28)])))
limb.Marsupialia <- lapply(limb.Marsupialia, function(x) x[,c(1,2,4,6,7,9,11,12,14,16,17,19)]) #subsampling for aspect ratio only parameters - length and prox-dist width

## phylogenetic regression - size correction across 100 trees
data.Marsupialia <- list()
fit_BM.Marsupialia <- list()
for(i in 1:100){
  data.Marsupialia[[i]] = list(limb=limb.Marsupialia[[i]], lgm=lgm.Marsupialia[[i]])
  fit_BM.Marsupialia[[i]] <- mvgls(limb~lgm, data=data.Marsupialia[[i]], pruned.trees.Marsupialia[[i]], model="BM", penalty="LASSO", method="LL")
}

resid.list.Marsupialia <- lapply(fit_BM.Marsupialia, function(x) as.data.frame(x$resid))


##------------------------------------------------------------------------------------------------------------------##
## creating rolling bins
##------------------------------------------------------------------------------------------------------------------##

## attributing rolling bins, depending on sample size


myascdat.Marsupialia <- mydat.Marsupialia[[1]] %>% arrange(Body_Mass)
rownames(myascdat.Marsupialia) <- myascdat.Marsupialia$Species

myascdat.Marsupialia$bins <- as.numeric(cut_number(myascdat.Marsupialia$Body_Mass, 6))
myascdat.Marsupialia %>% group_by(bins) %>% summarise(length(bins)) # 6 rolling bins, species counts per bin is 10/11

# now, let's make this bins overlap
bin.list.Marsupialia <- list() # creating a list to each fixed bin
for(j in 1:6){
  bin.list.Marsupialia[[j]] <- subset(myascdat.Marsupialia, bins== j)
  bin.list.Marsupialia[[j]] <- data.frame(Species=bin.list.Marsupialia[[j]]$Species, Group=bin.list.Marsupialia[[j]]$Group, 
                                          Body_Mass=bin.list.Marsupialia[[j]]$Body_Mass, bins=bin.list.Marsupialia[[j]]$bins)
}
str(bin.list.Marsupialia)

for(j in 1:6){
  rownames(bin.list.Marsupialia[[j]]) <- bin.list.Marsupialia[[j]]$Species
}

## creating the overlapping subsets:
## counts are different, so I'll custom bin per bin to make sure they all have the same number of individuals (16). 
overl.bins.Marsupialia <- list()
overl.bins.Marsupialia[[1]] = rbind(bin.list.Marsupialia[[1]], bin.list.Marsupialia[[2]][1:5,]) # 5 from bin 2
overl.bins.Marsupialia[[2]] <-  rbind(bin.list.Marsupialia[[2-1]][9:nrow(bin.list.Marsupialia[[2-1]]),], bin.list.Marsupialia[[2]], bin.list.Marsupialia[[2+1]][1:3,]) # 3 from each
overl.bins.Marsupialia[[3]] <-  rbind(bin.list.Marsupialia[[3-1]][9:nrow(bin.list.Marsupialia[[3-1]]),], bin.list.Marsupialia[[3]], bin.list.Marsupialia[[3+1]][1:3,]) # 3, 2
overl.bins.Marsupialia[[4]] <-  rbind(bin.list.Marsupialia[[4-1]][9:nrow(bin.list.Marsupialia[[4-1]]),], bin.list.Marsupialia[[4]], bin.list.Marsupialia[[4+1]][1:3,])
overl.bins.Marsupialia[[5]] <-  rbind(bin.list.Marsupialia[[5-1]][8:nrow(bin.list.Marsupialia[[5-1]]),], bin.list.Marsupialia[[5]], bin.list.Marsupialia[[5+1]][1:3,])
overl.bins.Marsupialia[[6]] = rbind(bin.list.Marsupialia[[5]][6:nrow(bin.list.Marsupialia[[5]]),],bin.list.Marsupialia[[6]]) # 5 from bin 4

str(overl.bins.Marsupialia[[6]])

for(k in 1:6){
  overl.bins.Marsupialia[[k]]$roll.bins <- k
}

str(overl.bins.Marsupialia) #Now each bin contains 18 taxa with 5-6 taxa overlapping


## next - create a list of lists: outer is each fitted model across 100 trees, inside is corresponding residuals organized into rolling bin lists
# inside = match residual rownames with overl.bins bin row names

resid.overl.bins.Marsupialia <- list()
for (i in 1:length(resid.list.Marsupialia)) {
  resid.overl.bins.Marsupialia[[i]] <- list()  # Initialize a sublist for each dataframe in resid.list
  
  # Iterate over each dataframe in overl.bins
  for (j in 1:length(overl.bins.Marsupialia)) {
    resid.overl.bins.Marsupialia[[i]][[j]] <- as.matrix(resid.list.Marsupialia[[i]][rownames(overl.bins.Marsupialia[[j]]), ])
  }
}

##------------------------------------------------------------------------------------------------------------------##
## disparity
##------------------------------------------------------------------------------------------------------------------##

# bootstrap data
boot_resid.Marsupialia <- list()
for(i in 1:100){
  boot_resid.Marsupialia[[i]] <- lapply(resid.overl.bins.Marsupialia[[i]], function(x) boot.matrix(x)) # default is 100 replications (*100 trees = 10,000 replications per bin)
}

# calculate disparity from bootstrapped data
disp_resid.Marsupialia <- list()
for(i in 1:100){
  disp_resid.Marsupialia[[i]] <- lapply(boot_resid.Marsupialia[[i]], function(x) dispRity(x, metric = c(sum, variances), verbose=TRUE))
}

# extracting results to plot
mat.disp.Marsupialia <- matrix(ncol=6, nrow=100)
list.disp.Marsupialia <- list()
for(k in 1:100){
  list.disp.Marsupialia[[k]] <- list()
  for(i in 1:6){
    mat.disp.Marsupialia[,i] <- disp_resid.Marsupialia[[k]][[i]]$disparity[[1]][[2]]
    colnames(mat.disp.Marsupialia) <- 1:6
    list.disp.Marsupialia[[k]] <-   mat.disp.Marsupialia
  }
}

disp.df.Marsupialia <- do.call(rbind, list.disp.Marsupialia)
str(disp.df.Marsupialia)

disp.df.plot.Marsupialia <- melt(disp.df.Marsupialia)

disp.df.plot.Marsupialia <- disp.df.plot.Marsupialia[,-1]
colnames(disp.df.plot.Marsupialia) <- c("rolling.bin", "disparity")


# transforming rolling bin axis into body mass mean values
# log10 value
logmed.rb.Marsupialia <- sapply(overl.bins.Marsupialia, function(x) median(log10(x$Body_Mass)))
logmed.rb.Marsupialia <- format(round(logmed.rb.Marsupialia,2), nsmall=2)

logmed.rb.vec.Marsupialia <- list()
for(i in 1:6){
  logmed.rb.vec.Marsupialia[[i]] <- rep(logmed.rb.Marsupialia [i], 10000)
}
logmed.rb.vec.Marsupialia <- unlist(logmed.rb.vec.Marsupialia)
disp.df.plot.Marsupialia$log.med.rb.Marsupialia <- as.numeric(logmed.rb.vec.Marsupialia)



##=========================================================================================================================================##
#### b) Cetartiodactyla ####
##=========================================================================================================================================##

# subset data
mydat.Cetartiodactyla <- list()
for(i in 1:(length(mydat.list))){
  mydat.Cetartiodactyla[[i]] <- subset(mydat.list[[i]], Group=="Cetartiodactyla")
}

# prune tree
pruned.trees.Cetartiodactyla <- list()
for( i in 1: length(pruned_trees)){
  pruned.trees.Cetartiodactyla[[i]] <- drop.tip(pruned_trees[[i]], (name.check(pruned_trees[[i]], mydat.Cetartiodactyla[[i]])$tree_not_data))
}

for(i in 1:(length(mydat.Cetartiodactyla))){
  mydat.Cetartiodactyla[[i]] <- mydat.Cetartiodactyla[[i]][pruned.trees.Cetartiodactyla[[i]]$tip.label,]
}

# body size correction
lgm.Cetartiodactyla <- lapply(mydat.Cetartiodactyla, function(x) log10(apply(x[,c(9:28)], 1, geometric.mean))) # using limb traits only to calculate geometric means
limb.Cetartiodactyla <- lapply(mydat.Cetartiodactyla, function(x) as.matrix(log10(x[,c(9:28)])))
limb.Cetartiodactyla <- lapply(limb.Cetartiodactyla, function(x) x[,c(1,2,4,6,7,9,11,12,14,16,17,19)]) #subsampling for aspect ratio only parameters - length and mid-shaft

## phylogenetic regression - size correction
data.Cetartiodactyla <- list()
fit_BM.Cetartiodactyla <- list()
for(i in 1:100){
  data.Cetartiodactyla[[i]] = list(limb=limb.Cetartiodactyla[[i]], lgm=lgm.Cetartiodactyla[[i]])
  fit_BM.Cetartiodactyla[[i]] <- mvgls(limb~lgm, data=data.Cetartiodactyla[[i]], pruned.trees.Cetartiodactyla[[i]], model="BM", penalty="LASSO", method="LL")
}

resid.list.Cetartiodactyla <- lapply(fit_BM.Cetartiodactyla, function(x) as.data.frame(x$resid))



##------------------------------------------------------------------------------------------------------------------##
## creating rolling bins ##
##------------------------------------------------------------------------------------------------------------------##

## attributing 5 bins 
myascdat.Cetartiodactyla <- mydat.Cetartiodactyla[[1]] %>% arrange(Body_Mass)
rownames(myascdat.Cetartiodactyla) <- myascdat.Cetartiodactyla$Species

myascdat.Cetartiodactyla$bins <- as.numeric(cut_number(myascdat.Cetartiodactyla$Body_Mass, 8)) # 8 bins
myascdat.Cetartiodactyla %>% group_by(bins) %>% summarise(length(bins)) #10-11 species counts per bin

# now, let's make this bins overlap,
bin.list.Cetartiodactyla <- list() # creating a list to each fixed bin
for(j in 1:8){
  bin.list.Cetartiodactyla[[j]] <- subset(myascdat.Cetartiodactyla, bins== j)
  bin.list.Cetartiodactyla[[j]] <- data.frame(Species=bin.list.Cetartiodactyla[[j]]$Species, Group=bin.list.Cetartiodactyla[[j]]$Group, 
                                          Body_Mass=bin.list.Cetartiodactyla[[j]]$Body_Mass, bins=bin.list.Cetartiodactyla[[j]]$bins)
}
str(bin.list.Cetartiodactyla)

for(j in 1:8){
  rownames(bin.list.Cetartiodactyla[[j]]) <- bin.list.Cetartiodactyla[[j]]$Species
}


## creating the overlapping subsets:
## counts are different, so I'll custom bin per bin to make sure they all have the same number of individuals (16). 
overl.bins.Cetartiodactyla <- list()
overl.bins.Cetartiodactyla[[1]] = rbind(bin.list.Cetartiodactyla[[1]], bin.list.Cetartiodactyla[[2]][1:5,]) 
overl.bins.Cetartiodactyla[[2]] <-  rbind(bin.list.Cetartiodactyla[[2-1]][10:nrow(bin.list.Cetartiodactyla[[2-1]]),], bin.list.Cetartiodactyla[[2]], bin.list.Cetartiodactyla[[2+1]][1:3,]) 
overl.bins.Cetartiodactyla[[3]] <-  rbind(bin.list.Cetartiodactyla[[3-1]][9:nrow(bin.list.Cetartiodactyla[[3-1]]),], bin.list.Cetartiodactyla[[3]], bin.list.Cetartiodactyla[[3+1]][1:3,]) 
overl.bins.Cetartiodactyla[[4]] <-  rbind(bin.list.Cetartiodactyla[[4-1]][9:nrow(bin.list.Cetartiodactyla[[4-1]]),], bin.list.Cetartiodactyla[[4]], bin.list.Cetartiodactyla[[4+1]][1:3,])
overl.bins.Cetartiodactyla[[5]] <-  rbind(bin.list.Cetartiodactyla[[5-1]][10:nrow(bin.list.Cetartiodactyla[[5-1]]),], bin.list.Cetartiodactyla[[5]], bin.list.Cetartiodactyla[[5+1]][1:3,])
overl.bins.Cetartiodactyla[[6]] <-  rbind(bin.list.Cetartiodactyla[[6-1]][9:nrow(bin.list.Cetartiodactyla[[6-1]]),], bin.list.Cetartiodactyla[[6]], bin.list.Cetartiodactyla[[6+1]][1:3,])
overl.bins.Cetartiodactyla[[7]] <-  rbind(bin.list.Cetartiodactyla[[7-1]][9:nrow(bin.list.Cetartiodactyla[[7-1]]),], bin.list.Cetartiodactyla[[7]], bin.list.Cetartiodactyla[[7+1]][1:3,])
overl.bins.Cetartiodactyla[[8]] = rbind(bin.list.Cetartiodactyla[[7]][7:nrow(bin.list.Cetartiodactyla[[7]]),],bin.list.Cetartiodactyla[[8]]) 

str(overl.bins.Cetartiodactyla[[8]])
overl.bins.Cetartiodactyla[[8]] %>% group_by(bins)%>%summarise(length(bins))
str(overl.bins.Cetartiodactyla)

for(k in 1:8){
  overl.bins.Cetartiodactyla[[k]]$roll.bins <- k
}

## next - create a list of lists: outer is each fitted model across 100 trees, inside is corresponding residuals organized into rolling bin lists
# inside = match residual rownames with overl.bins bin row names
resid.overl.bins.Cetartiodactyla <- list()
for (i in 1:length(resid.list.Cetartiodactyla)) {
  resid.overl.bins.Cetartiodactyla[[i]] <- list()  # Initialize a sublist for each dataframe in resid.list
  
  # Iterate over each dataframe in overl.bins
  for (j in 1:length(overl.bins.Cetartiodactyla)) {
    resid.overl.bins.Cetartiodactyla[[i]][[j]] <- as.matrix(resid.list.Cetartiodactyla[[i]][rownames(overl.bins.Cetartiodactyla[[j]]), ])
  }
}

##------------------------------------------------------------------------------------------------------------------##
## disparity
##------------------------------------------------------------------------------------------------------------------##

# bootstrap data
boot_resid.Cetartiodactyla <- list()
for(i in 1:100){
  boot_resid.Cetartiodactyla[[i]] <- lapply(resid.overl.bins.Cetartiodactyla[[i]], function(x) boot.matrix(x)) # default is 100 replications (*100 trees = 10,000 replications per bin)
}

# calculate disparity from bootstraped data
disp_resid.Cetartiodactyla <- list()
for(i in 1:100){
  disp_resid.Cetartiodactyla[[i]] <- lapply(boot_resid.Cetartiodactyla[[i]], function(x) dispRity(x, metric = c(sum, variances), verbose=TRUE))
}

# extracting results to plot
mat.disp.Cetartiodactyla <- matrix(ncol=8, nrow=100)
list.disp.Cetartiodactyla <- list()
for(k in 1:100){
  list.disp.Cetartiodactyla[[k]] <- list()
  for(i in 1:8){
    mat.disp.Cetartiodactyla[,i] <- disp_resid.Cetartiodactyla[[k]][[i]]$disparity[[1]][[2]]
    colnames(mat.disp.Cetartiodactyla) <- 1:8
    list.disp.Cetartiodactyla[[k]] <-   mat.disp.Cetartiodactyla
  }
}

disp.df.Cetartiodactyla <- do.call(rbind, list.disp.Cetartiodactyla)
str(disp.df.Cetartiodactyla)

disp.df.plot.Cetartiodactyla <- melt(disp.df.Cetartiodactyla)
disp.df.plot.Cetartiodactyla %>% group_by(Var2) %>% count(Var2)

disp.df.plot.Cetartiodactyla <- disp.df.plot.Cetartiodactyla[,-1]
colnames(disp.df.plot.Cetartiodactyla) <- c("rolling.bin", "disparity")

# transforming rolling bin axis into body mass mean values
# log10 value
logmed.rb.Cetartiodactyla <- sapply(overl.bins.Cetartiodactyla, function(x) median(log10(x$Body_Mass)))
logmed.rb.Cetartiodactyla <- format(round(logmed.rb.Cetartiodactyla,2), nsmall=2)

logmed.rb.vec.Cetartiodactyla <- list()
for(i in 1:8){
  logmed.rb.vec.Cetartiodactyla[[i]] <- rep(logmed.rb.Cetartiodactyla [i], 10000)
}
logmed.rb.vec.Cetartiodactyla <- unlist(logmed.rb.vec.Cetartiodactyla)
disp.df.plot.Cetartiodactyla$log.med.rb.Cetartiodactyla <- as.numeric(logmed.rb.vec.Cetartiodactyla)


##=========================================================================================================================================##
#### c) Primates ####
##=========================================================================================================================================##

# subset data
mydat.Primates <- list()
for(i in 1:(length(mydat.list))){
  mydat.Primates[[i]] <- subset(mydat.list[[i]], Group=="Primates")
}

# prune tree
pruned.trees.Primates <- list()
for( i in 1: length(pruned_trees)){
  pruned.trees.Primates[[i]] <- drop.tip(pruned_trees[[i]], (name.check(pruned_trees[[i]], mydat.Primates[[i]])$tree_not_data))
}

for(i in 1:(length(mydat.Primates))){
  mydat.Primates[[i]] <- mydat.Primates[[i]][pruned.trees.Primates[[i]]$tip.label,]
}

# body size correction
lgm.Primates <- lapply(mydat.Primates, function(x) log10(apply(x[,c(9:28)], 1, geometric.mean))) # using limb traits only to calculate geometric means
limb.Primates <- lapply(mydat.Primates, function(x) as.matrix(log10(x[,c(9:28)])))
limb.Primates <- lapply(limb.Primates, function(x) x[,c(1,2,4,6,7,9,11,12,14,16,17,19)]) #subsampling for aspect ratio only parameters - length and mid-shaft

## phylogenetic regression - size correction
data.Primates <- list()
fit_BM.Primates <- list()
for(i in 1:100){
  data.Primates[[i]] = list(limb=limb.Primates[[i]], lgm=lgm.Primates[[i]])
  fit_BM.Primates[[i]] <- mvgls(limb~lgm, data=data.Primates[[i]], pruned.trees.Primates[[i]], model="BM", penalty="LASSO", method="LL")
}

resid.list.Primates <- lapply(fit_BM.Primates, function(x) as.data.frame(x$resid))

##------------------------------------------------------------------------------------------------------------------##
## creating rolling bins ##
##------------------------------------------------------------------------------------------------------------------##

## attributing fixed bins of minimum n 10
myascdat.Primates <- mydat.Primates[[1]] %>% arrange(Body_Mass)
rownames(myascdat.Primates) <- myascdat.Primates$Species

myascdat.Primates$bins <- as.numeric(cut_number(myascdat.Primates$Body_Mass, 6)) # 6 rolling bins
myascdat.Primates %>% group_by(bins) %>% summarise(length(bins)) #11 species per bin (7 bins would have resulted in 9-10 sp per bin)

# now, let's make this bins overlap, creating the rolling bins, 13 taxa per bin
bin.list.Primates <- list() # creating a list to each fixed bin
for(j in 1:6){
  bin.list.Primates[[j]] <- subset(myascdat.Primates, bins== j)
  bin.list.Primates[[j]] <- data.frame(Species=bin.list.Primates[[j]]$Species, Group=bin.list.Primates[[j]]$Group, 
                                              Body_Mass=bin.list.Primates[[j]]$Body_Mass, bins=bin.list.Primates[[j]]$bins)
}
str(bin.list.Primates)

for(j in 1:6){
  rownames(bin.list.Primates[[j]]) <- bin.list.Primates[[j]]$Species
}

## creating the overlapping subsets:

## counts are different, so I'll custom bin per bin to make sure they all have the same number of individuals (16). 
overl.bins.Primates <- list()
overl.bins.Primates[[1]] = rbind(bin.list.Primates[[1]], bin.list.Primates[[2]][1:5,])
overl.bins.Primates[[2]] <-  rbind(bin.list.Primates[[2-1]][9:nrow(bin.list.Primates[[2-1]]),], bin.list.Primates[[2]], bin.list.Primates[[2+1]][1:2,]) 
overl.bins.Primates[[3]] <-  rbind(bin.list.Primates[[3-1]][9:nrow(bin.list.Primates[[3-1]]),], bin.list.Primates[[3]], bin.list.Primates[[3+1]][1:2,]) 
overl.bins.Primates[[4]] <-  rbind(bin.list.Primates[[4-1]][9:nrow(bin.list.Primates[[4-1]]),], bin.list.Primates[[4]], bin.list.Primates[[4+1]][1:2,]) 
overl.bins.Primates[[5]] <-  rbind(bin.list.Primates[[5-1]][9:nrow(bin.list.Primates[[5-1]]),], bin.list.Primates[[5]], bin.list.Primates[[5+1]][1:2,]) 
overl.bins.Primates[[6]] = rbind(bin.list.Primates[[5]][7:nrow(bin.list.Primates[[5]]),],bin.list.Primates[[6]]) 


str(overl.bins.Primates[[6]])
overl.bins.Primates[[6]] %>% group_by(bins)%>%summarise(length(bins))
str(overl.bins.Primates)


for(k in 1:6){
  overl.bins.Primates[[k]]$roll.bins <- k
}


str(overl.bins.Primates)

## next - create a list of lists: outer is each fitted model across 100 trees, inside is corresponding residuals organized into rolling bin lists
# inside = match residual rownames with overl.bins bin row names
resid.overl.bins.Primates <- list()
for (i in 1:length(resid.list.Primates)) {
  resid.overl.bins.Primates[[i]] <- list()  # Initialize a sublist for each dataframe in resid.list
  
  # Iterate over each dataframe in overl.bins
  for (j in 1:length(overl.bins.Primates)) {
    resid.overl.bins.Primates[[i]][[j]] <- as.matrix(resid.list.Primates[[i]][rownames(overl.bins.Primates[[j]]), ])
  }
}

##------------------------------------------------------------------------------------------------------------------##
## disparity
##------------------------------------------------------------------------------------------------------------------##

# bootstrap data
boot_resid.Primates <- list()
for(i in 1:100){
  boot_resid.Primates[[i]] <- lapply(resid.overl.bins.Primates[[i]], function(x) boot.matrix(x)) # default is 100 replications (*100 trees = 10,000 replications per bin)
}

# calculate disparity from bootstraped data
disp_resid.Primates <- list()
for(i in 1:100){
  disp_resid.Primates[[i]] <- lapply(boot_resid.Primates[[i]], function(x) dispRity(x, metric = c(sum, variances), verbose=TRUE))
}

# extracting results to plot
mat.disp.Primates <- matrix(ncol=6, nrow=100)
list.disp.Primates <- list()
for(k in 1:100){
  list.disp.Primates[[k]] <- list()
  for(i in 1:6){
    mat.disp.Primates[,i] <- disp_resid.Primates[[k]][[i]]$disparity[[1]][[2]]
    colnames(mat.disp.Primates) <- 1:6
    list.disp.Primates[[k]] <-   mat.disp.Primates
  }
}

disp.df.Primates <- do.call(rbind, list.disp.Primates)
str(disp.df.Primates)

disp.df.plot.Primates <- melt(disp.df.Primates)

disp.df.plot.Primates <- disp.df.plot.Primates[,-1]
colnames(disp.df.plot.Primates) <- c("rolling.bin", "disparity")

# transforming rolling bin axis into body mass mean values
# log10 value
logmed.rb.Primates <- sapply(overl.bins.Primates, function(x) median(log10(x$Body_Mass)))
logmed.rb.Primates <- format(round(logmed.rb.Primates,2), nsmall=2)

logmed.rb.vec.Primates <- list()
for(i in 1:6){
  logmed.rb.vec.Primates[[i]] <- rep(logmed.rb.Primates [i], 10000)
}
logmed.rb.vec.Primates <- unlist(logmed.rb.vec.Primates)
disp.df.plot.Primates$log.med.rb.Primates <- as.numeric(logmed.rb.vec.Primates)


##=========================================================================================================================================##
#### d) Rodentia ####
##=========================================================================================================================================##


# subset data
mydat.Rodentia <- list()
for(i in 1:(length(mydat.list))){
  mydat.Rodentia[[i]] <- subset(mydat.list[[i]], Group=="Rodentia")
}

# prune tree
pruned.trees.Rodentia <- list()
for( i in 1: length(pruned_trees)){
  pruned.trees.Rodentia[[i]] <- drop.tip(pruned_trees[[i]], (name.check(pruned_trees[[i]], mydat.Rodentia[[i]])$tree_not_data))
}

for(i in 1:(length(mydat.Rodentia))){
  mydat.Rodentia[[i]] <- mydat.Rodentia[[i]][pruned.trees.Rodentia[[i]]$tip.label,]
}

# body size correction
lgm.Rodentia <- lapply(mydat.Rodentia, function(x) log10(apply(x[,c(9:28)], 1, geometric.mean))) # using limb traits only to calculate geometric means
limb.Rodentia <- lapply(mydat.Rodentia, function(x) as.matrix(log10(x[,c(9:28)])))
limb.Rodentia <- lapply(limb.Rodentia, function(x) x[,c(1,2,4,6,7,9,11,12,14,16,17,19)]) #subsampling for aspect ratio only parameters - length and mid-shaft

## phylogenetic regression - size correction
data.Rodentia <- list()
fit_BM.Rodentia <- list()
for(i in 1:100){
  data.Rodentia[[i]] = list(limb=limb.Rodentia[[i]], lgm=lgm.Rodentia[[i]])
  fit_BM.Rodentia[[i]] <- mvgls(limb~lgm, data=data.Rodentia[[i]], pruned.trees.Rodentia[[i]], model="BM", penalty="LASSO", method="LL")
}

resid.list.Rodentia <- lapply(fit_BM.Rodentia, function(x) as.data.frame(x$resid))

##------------------------------------------------------------------------------------------------------------------##
## creating rolling bins ##
##------------------------------------------------------------------------------------------------------------------##

## attributing 5 bins 
myascdat.Rodentia <- mydat.Rodentia[[1]] %>% arrange(Body_Mass)
rownames(myascdat.Rodentia) <- myascdat.Rodentia$Species

myascdat.Rodentia$bins <- as.numeric(cut_number(myascdat.Rodentia$Body_Mass, 10))
myascdat.Rodentia %>% group_by(bins) %>% summarise(length(bins)) # 25/26 species counts per bin 

# now, let's make this bins overlap, creating the rolling bins, 32 taxa per bin
bin.list.Rodentia <- list() # creating a list to each fixed bin
for(j in 1:10){
  bin.list.Rodentia[[j]] <- subset(myascdat.Rodentia, bins== j)
  bin.list.Rodentia[[j]] <- data.frame(Species=bin.list.Rodentia[[j]]$Species, Group=bin.list.Rodentia[[j]]$Group, 
                                       Body_Mass=bin.list.Rodentia[[j]]$Body_Mass, bins=bin.list.Rodentia[[j]]$bins)
}
str(bin.list.Rodentia)

for(j in 1:10){
  rownames(bin.list.Rodentia[[j]]) <- bin.list.Rodentia[[j]]$Species
}

## creating the overlapping subsets:

## counts are different, so I'll custom bin per bin to make sure they all have the same number of individuals (55). 
overl.bins.Rodentia <- list()
overl.bins.Rodentia[[1]] = rbind(bin.list.Rodentia[[1]], bin.list.Rodentia[[2]][1:6,]) 
overl.bins.Rodentia[[2]] <-  rbind(bin.list.Rodentia[[2-1]][23:nrow(bin.list.Rodentia[[2-1]]),], bin.list.Rodentia[[2]], bin.list.Rodentia[[2+1]][1:3,]) 
overl.bins.Rodentia[[3]] <-  rbind(bin.list.Rodentia[[3-1]][22:nrow(bin.list.Rodentia[[3-1]]),], bin.list.Rodentia[[3]], bin.list.Rodentia[[3+1]][1:3,]) 
overl.bins.Rodentia[[4]] <-  rbind(bin.list.Rodentia[[4-1]][23:nrow(bin.list.Rodentia[[4-1]]),], bin.list.Rodentia[[4]], bin.list.Rodentia[[4+1]][1:3,]) 
overl.bins.Rodentia[[5]] <-  rbind(bin.list.Rodentia[[5-1]][23:nrow(bin.list.Rodentia[[5-1]]),], bin.list.Rodentia[[5]], bin.list.Rodentia[[5+1]][1:3,]) 
overl.bins.Rodentia[[6]] <-  rbind(bin.list.Rodentia[[6-1]][22:nrow(bin.list.Rodentia[[6-1]]),], bin.list.Rodentia[[6]], bin.list.Rodentia[[6+1]][1:3,]) 
overl.bins.Rodentia[[7]] <-  rbind(bin.list.Rodentia[[7-1]][23:nrow(bin.list.Rodentia[[7-1]]),], bin.list.Rodentia[[7]], bin.list.Rodentia[[7+1]][1:3,]) 
overl.bins.Rodentia[[8]] <-  rbind(bin.list.Rodentia[[8-1]][23:nrow(bin.list.Rodentia[[8-1]]),], bin.list.Rodentia[[8]], bin.list.Rodentia[[8+1]][1:3,]) 
overl.bins.Rodentia[[9]] <-  rbind(bin.list.Rodentia[[9-1]][22:nrow(bin.list.Rodentia[[9-1]]),], bin.list.Rodentia[[9]], bin.list.Rodentia[[9+1]][1:3,]) 
overl.bins.Rodentia[[10]] = rbind(bin.list.Rodentia[[9]][20:nrow(bin.list.Rodentia[[9]]),],bin.list.Rodentia[[10]]) # 5 from bin 9

str(overl.bins.Rodentia)

for(k in 1:10){
  overl.bins.Rodentia[[k]]$roll.bins <- k
}


## next - create a list of lists: outer is each fitted model across 100 trees, inside is corresponding residuals organized into rolling bin lists
# inside = match residual rownames with overl.bins bin row names
resid.overl.bins.Rodentia <- list()
for (i in 1:length(resid.list.Rodentia)) {
  resid.overl.bins.Rodentia[[i]] <- list()  # Initialize a sublist for each dataframe in resid.list
  
  # Iterate over each dataframe in overl.bins
  for (j in 1:length(overl.bins.Rodentia)) {
    resid.overl.bins.Rodentia[[i]][[j]] <- as.matrix(resid.list.Rodentia[[i]][rownames(overl.bins.Rodentia[[j]]), ])
  }
}

##------------------------------------------------------------------------------------------------------------------##
## disparity
##------------------------------------------------------------------------------------------------------------------##

# bootstrap data
boot_resid.Rodentia <- list()
for(i in 1:100){
  boot_resid.Rodentia[[i]] <- lapply(resid.overl.bins.Rodentia[[i]], function(x) boot.matrix(x)) # default is 100 replications (*100 trees = 10,000 replications per bin)
}

# calculate disparity from bootstrapped data
disp_resid.Rodentia <- list()
for(i in 1:100){
  disp_resid.Rodentia[[i]] <- lapply(boot_resid.Rodentia[[i]], function(x) dispRity(x, metric = c(sum, variances), verbose=TRUE))
}

# extracting results to plot
mat.disp.Rodentia <- matrix(ncol=10, nrow=100)
list.disp.Rodentia <- list()
for(k in 1:100){
  list.disp.Rodentia[[k]] <- list()
  for(i in 1:10){
    mat.disp.Rodentia[,i] <- disp_resid.Rodentia[[k]][[i]]$disparity[[1]][[2]]
    colnames(mat.disp.Rodentia) <- 1:10
    list.disp.Rodentia[[k]] <-   mat.disp.Rodentia
  }
}

disp.df.Rodentia <- do.call(rbind, list.disp.Rodentia)
str(disp.df.Rodentia)

disp.df.plot.Rodentia <- melt(disp.df.Rodentia)

disp.df.plot.Rodentia <- disp.df.plot.Rodentia[,-1]
colnames(disp.df.plot.Rodentia) <- c("rolling.bin", "disparity")

# transforming rolling bin axis into body mass mean values
# log10 value

logmed.rb.Rodentia <- sapply(overl.bins.Rodentia, function(x) median(log10(x$Body_Mass)))
logmed.rb.Rodentia <- format(round(logmed.rb.Rodentia,2), nsmall=2)

logmed.rb.vec.Rodentia <- list()
for(i in 1:10){
  logmed.rb.vec.Rodentia[[i]] <- rep(logmed.rb.Rodentia [i], 10000)
}
logmed.rb.vec.Rodentia <- unlist(logmed.rb.vec.Rodentia)
disp.df.plot.Rodentia$log.med.rb.Rodentia <- as.numeric(logmed.rb.vec.Rodentia)


##=========================================================================================================================================##
#### e) Carnivora ####
##=========================================================================================================================================##

# subset data
mydat.Carnivora <- list()
for(i in 1:(length(mydat.list))){
  mydat.Carnivora[[i]] <- subset(mydat.list[[i]], Group=="Carnivora")
}

# prune tree
pruned.trees.Carnivora <- list()
for( i in 1: length(pruned_trees)){
  pruned.trees.Carnivora[[i]] <- drop.tip(pruned_trees[[i]], (name.check(pruned_trees[[i]], mydat.Carnivora[[i]])$tree_not_data))
}

for(i in 1:(length(mydat.Carnivora))){
  mydat.Carnivora[[i]] <- mydat.Carnivora[[i]][pruned.trees.Carnivora[[i]]$tip.label,]
}

# body size correction
lgm.Carnivora <- lapply(mydat.Carnivora, function(x) log10(apply(x[,c(9:28)], 1, geometric.mean))) # using limb traits only to calculate geometric means
limb.Carnivora <- lapply(mydat.Carnivora, function(x) as.matrix(log10(x[,c(9:28)])))
limb.Carnivora <- lapply(limb.Carnivora, function(x) x[,c(1,2,4,6,7,9,11,12,14,16,17,19)]) #subsampling for aspect ratio only parameters - length and mid-shaft

## phylogenetic regression - size correction
data.Carnivora <- list()
fit_BM.Carnivora <- list()
for(i in 1:100){
  data.Carnivora[[i]] = list(limb=limb.Carnivora[[i]], lgm=lgm.Carnivora[[i]])
  fit_BM.Carnivora[[i]] <- mvgls(limb~lgm, data=data.Carnivora[[i]], pruned.trees.Carnivora[[i]], model="BM", penalty="LASSO", method="LL")
}

resid.list.Carnivora <- lapply(fit_BM.Carnivora, function(x) as.data.frame(x$resid))


##------------------------------------------------------------------------------------------------------------------##
## creating rolling bins ##
##------------------------------------------------------------------------------------------------------------------##

## attributing  bins 
myascdat.Carnivora <- mydat.Carnivora[[1]] %>% arrange(Body_Mass)
rownames(myascdat.Carnivora) <- myascdat.Carnivora$Species

myascdat.Carnivora$bins <- as.numeric(cut_number(myascdat.Carnivora$Body_Mass, 8)) # 8 bins
myascdat.Carnivora %>% group_by(bins) %>% summarise(length(bins)) #11-12 species counts per bin

# now, let's make this bins overlap
bin.list.Carnivora <- list() # creating a list to each fixed bin
for(j in 1:8){
  bin.list.Carnivora[[j]] <- subset(myascdat.Carnivora, bins== j)
  bin.list.Carnivora[[j]] <- data.frame(Species=bin.list.Carnivora[[j]]$Species, Group=bin.list.Carnivora[[j]]$Group, 
                                       Body_Mass=bin.list.Carnivora[[j]]$Body_Mass, bins=bin.list.Carnivora[[j]]$bins)
}
str(bin.list.Carnivora)

for(j in 1:8){
  rownames(bin.list.Carnivora[[j]]) <- bin.list.Carnivora[[j]]$Species
}

## creating the overlapping subsets:

## counts are different, so I'll custom bin per bin to make sure they all have the same number of individuals (22). 
overl.bins.Carnivora <- list()
overl.bins.Carnivora[[1]] = rbind(bin.list.Carnivora[[1]], bin.list.Carnivora[[2]][1:4,]) 
overl.bins.Carnivora[[2]] <-  rbind(bin.list.Carnivora[[2-1]][10:nrow(bin.list.Carnivora[[2-1]]),], bin.list.Carnivora[[2]], bin.list.Carnivora[[2+1]][1:2,]) 
overl.bins.Carnivora[[3]] <-  rbind(bin.list.Carnivora[[3-1]][9:nrow(bin.list.Carnivora[[3-1]]),], bin.list.Carnivora[[3]], bin.list.Carnivora[[3+1]][1:2,]) 
overl.bins.Carnivora[[4]] <-  rbind(bin.list.Carnivora[[4-1]][9:nrow(bin.list.Carnivora[[4-1]]),], bin.list.Carnivora[[4]], bin.list.Carnivora[[4+1]][1:2,]) 
overl.bins.Carnivora[[5]] <-  rbind(bin.list.Carnivora[[5-1]][9:nrow(bin.list.Carnivora[[5-1]]),], bin.list.Carnivora[[5]], bin.list.Carnivora[[5+1]][1:2,])
overl.bins.Carnivora[[6]] <-  rbind(bin.list.Carnivora[[6-1]][9:nrow(bin.list.Carnivora[[6-1]]),], bin.list.Carnivora[[6]], bin.list.Carnivora[[6+1]][1:2,])
overl.bins.Carnivora[[7]] <-  rbind(bin.list.Carnivora[[7-1]][9:nrow(bin.list.Carnivora[[7-1]]),], bin.list.Carnivora[[7]], bin.list.Carnivora[[7+1]][1:2,])
overl.bins.Carnivora[[8]] = rbind(bin.list.Carnivora[[7]][7:nrow(bin.list.Carnivora[[7]]),],bin.list.Carnivora[[8]]) 

str(overl.bins.Carnivora[[8]])
str(overl.bins.Carnivora)

for(k in 1:8){
  overl.bins.Carnivora[[k]]$roll.bins <- k
}


## next - create a list of lists: outer is each fitted model across 100 trees, inside is corresponding residuals organized into rolling bin lists
# inside = match residual rownames with overl.bins bin row names

resid.overl.bins.Carnivora <- list()
for (i in 1:length(resid.list.Carnivora)) {
  resid.overl.bins.Carnivora[[i]] <- list()  # Initialize a sublist for each dataframe in resid.list
  
  # Iterate over each dataframe in overl.bins
  for (j in 1:length(overl.bins.Carnivora)) {
    resid.overl.bins.Carnivora[[i]][[j]] <- as.matrix(resid.list.Carnivora[[i]][rownames(overl.bins.Carnivora[[j]]), ])
  }
}


##------------------------------------------------------------------------------------------------------------------##
## disparity
##------------------------------------------------------------------------------------------------------------------##

# bootstrap data
boot_resid.Carnivora <- list()
for(i in 1:100){
  boot_resid.Carnivora[[i]] <- lapply(resid.overl.bins.Carnivora[[i]], function(x) boot.matrix(x)) # default is 100 replications (*100 trees = 10,000 replications per bin)
}

# calculate disparity from bootstraped data
disp_resid.Carnivora <- list()
for(i in 1:100){
  disp_resid.Carnivora[[i]] <- lapply(boot_resid.Carnivora[[i]], function(x) dispRity(x, metric = c(sum, variances), verbose=TRUE))
}

# extracting results to plot
mat.disp.Carnivora <- matrix(ncol=8, nrow=100)
list.disp.Carnivora <- list()
for(k in 1:100){
  list.disp.Carnivora[[k]] <- list()
  for(i in 1:8){
    mat.disp.Carnivora[,i] <- disp_resid.Carnivora[[k]][[i]]$disparity[[1]][[2]]
    colnames(mat.disp.Carnivora) <- 1:8
    list.disp.Carnivora[[k]] <-   mat.disp.Carnivora
  }
}

disp.df.Carnivora <- do.call(rbind, list.disp.Carnivora)
str(disp.df.Carnivora)

disp.df.plot.Carnivora <- melt(disp.df.Carnivora)

disp.df.plot.Carnivora <- disp.df.plot.Carnivora[,-1]
colnames(disp.df.plot.Carnivora) <- c("rolling.bin", "disparity")

# transforming rolling bin axis into body mass mean values
# log value

logmed.rb.Carnivora <- sapply(overl.bins.Carnivora, function(x) median(log10(x$Body_Mass)))
logmed.rb.Carnivora <- format(round(logmed.rb.Carnivora,2), nsmall=2)

logmed.rb.vec.Carnivora <- list()
for(i in 1:8){
  logmed.rb.vec.Carnivora[[i]] <- rep(logmed.rb.Carnivora [i], 10000)
}
logmed.rb.vec.Carnivora <- unlist(logmed.rb.vec.Carnivora)
disp.df.plot.Carnivora$log.med.rb.Carnivora <- as.numeric(logmed.rb.vec.Carnivora)


##=========================================================================================================================================##
#### f) Eulipotyphla ####
##=========================================================================================================================================##

# subset data
mydat.Eulipotyphla <- list()
for(i in 1:(length(mydat.list))){
  mydat.Eulipotyphla[[i]] <- subset(mydat.list[[i]], Group=="Eulipotyphla")
}

# prune tree
pruned.trees.Eulipotyphla <- list()
for( i in 1: length(pruned_trees)){
  pruned.trees.Eulipotyphla[[i]] <- drop.tip(pruned_trees[[i]], (name.check(pruned_trees[[i]], mydat.Eulipotyphla[[i]])$tree_not_data))
}

for(i in 1:(length(mydat.Eulipotyphla))){
  mydat.Eulipotyphla[[i]] <- mydat.Eulipotyphla[[i]][pruned.trees.Eulipotyphla[[i]]$tip.label,]
}

# body size correction
lgm.Eulipotyphla <- lapply(mydat.Eulipotyphla, function(x) log10(apply(x[,c(9:28)], 1, geometric.mean))) # using limb traits only to calculate geometric means
limb.Eulipotyphla <- lapply(mydat.Eulipotyphla, function(x) as.matrix(log10(x[,c(9:28)])))
limb.Eulipotyphla <- lapply(limb.Eulipotyphla, function(x) x[,c(1,2,4,6,7,9,11,12,14,16,17,19)]) #subsampling for aspect ratio only parameters - length and mid-shaft

## phylogenetic regression - size correction
data.Eulipotyphla <- list()
fit_BM.Eulipotyphla <- list()
for(i in 1:100){
  data.Eulipotyphla[[i]] = list(limb=limb.Eulipotyphla[[i]], lgm=lgm.Eulipotyphla[[i]])
  fit_BM.Eulipotyphla[[i]] <- mvgls(limb~lgm, data=data.Eulipotyphla[[i]], pruned.trees.Eulipotyphla[[i]], model="BM", penalty="LASSO", method="LL")
}

resid.list.Eulipotyphla <- lapply(fit_BM.Eulipotyphla, function(x) as.data.frame(x$resid))


##------------------------------------------------------------------------------------------------------------------##
## creating rolling bins
##------------------------------------------------------------------------------------------------------------------##

## attributing 5 bins 
myascdat.Eulipotyphla <- mydat.Eulipotyphla[[1]] %>% arrange(Body_Mass)
rownames(myascdat.Eulipotyphla) <- myascdat.Eulipotyphla$Species

myascdat.Eulipotyphla$bins <- as.numeric(cut_number(myascdat.Eulipotyphla$Body_Mass, 5))
myascdat.Eulipotyphla %>% group_by(bins) %>% summarise(length(bins)) #8-7 species counts per bin

# now, let's make this bins overlap, creating the rolling bins
bin.list.Eulipotyphla <- list() # creating a list to each fixed bin
for(j in 1:5){
  bin.list.Eulipotyphla[[j]] <- subset(myascdat.Eulipotyphla, bins== j)
  bin.list.Eulipotyphla[[j]] <- data.frame(Species=bin.list.Eulipotyphla[[j]]$Species, Group=bin.list.Eulipotyphla[[j]]$Group, 
                                        Body_Mass=bin.list.Eulipotyphla[[j]]$Body_Mass, bins=bin.list.Eulipotyphla[[j]]$bins)
}
str(bin.list.Eulipotyphla)

for(j in 1:5){
  rownames(bin.list.Eulipotyphla[[j]]) <- bin.list.Eulipotyphla[[j]]$Species
}


## creating the overlapping subsets:

## counts are different, so I'll custom bin per bin to make sure they all have the same number of individuals (12). 
overl.bins.Eulipotyphla <- list()
overl.bins.Eulipotyphla[[1]] = rbind(bin.list.Eulipotyphla[[1]], bin.list.Eulipotyphla[[2]][1:4,]) # 4 from bin 2
overl.bins.Eulipotyphla[[2]] <-  rbind(bin.list.Eulipotyphla[[2-1]][7:nrow(bin.list.Eulipotyphla[[2-1]]),], bin.list.Eulipotyphla[[2]], bin.list.Eulipotyphla[[2+1]][1:2,]) # 2 from each
overl.bins.Eulipotyphla[[3]] <-  rbind(bin.list.Eulipotyphla[[3-1]][6:nrow(bin.list.Eulipotyphla[[3-1]]),], bin.list.Eulipotyphla[[3]], bin.list.Eulipotyphla[[3+1]][1:2,]) # 3,2
overl.bins.Eulipotyphla[[4]] <-  rbind(bin.list.Eulipotyphla[[4-1]][6:nrow(bin.list.Eulipotyphla[[4-1]]),], bin.list.Eulipotyphla[[4]], bin.list.Eulipotyphla[[4+1]][1:2,]) # 2,2
overl.bins.Eulipotyphla[[5]] = rbind(bin.list.Eulipotyphla[[4]][5:nrow(bin.list.Eulipotyphla[[4]]),],bin.list.Eulipotyphla[[5]]) # 4 from bin 4

str(overl.bins.Eulipotyphla)


for(k in 1:5){
  overl.bins.Eulipotyphla[[k]]$roll.bins <- k
}
str(bin.list.Eulipotyphla)
str(overl.bins.Eulipotyphla)


## next - create a list of lists: outer is each fitted model across 100 trees, inside is corresponding residuals organized into rolling bin lists
# inside = match residual rownames with overl.bins bin row names

resid.overl.bins.Eulipotyphla <- list()
for (i in 1:length(resid.list.Eulipotyphla)) {
  resid.overl.bins.Eulipotyphla[[i]] <- list()  # Initialize a sublist for each dataframe in resid.list
  
  # Iterate over each dataframe in overl.bins
  for (j in 1:length(overl.bins.Eulipotyphla)) {
    resid.overl.bins.Eulipotyphla[[i]][[j]] <- as.matrix(resid.list.Eulipotyphla[[i]][rownames(overl.bins.Eulipotyphla[[j]]), ])
  }
}


##------------------------------------------------------------------------------------------------------------------##
## disparity
##------------------------------------------------------------------------------------------------------------------##

# bootstrap data
boot_resid.Eulipotyphla <- list()
for(i in 1:100){
  boot_resid.Eulipotyphla[[i]] <- lapply(resid.overl.bins.Eulipotyphla[[i]], function(x) boot.matrix(x)) # default is 100 replications (*100 trees = 10,000 replications per bin)
}

# calculate disparity from bootstraped data
disp_resid.Eulipotyphla <- list()
for(i in 1:100){
  disp_resid.Eulipotyphla[[i]] <- lapply(boot_resid.Eulipotyphla[[i]], function(x) dispRity(x, metric = c(sum, variances), verbose=TRUE))
}

# extracting results to plot
mat.disp.Eulipotyphla <- matrix(ncol=5, nrow=100)
list.disp.Eulipotyphla <- list()
for(k in 1:100){
  list.disp.Eulipotyphla[[k]] <- list()
  for(i in 1:5){
    mat.disp.Eulipotyphla[,i] <- disp_resid.Eulipotyphla[[k]][[i]]$disparity[[1]][[2]]
    colnames(mat.disp.Eulipotyphla) <- 1:5
    list.disp.Eulipotyphla[[k]] <-   mat.disp.Eulipotyphla
  }
}

disp.df.Eulipotyphla <- do.call(rbind, list.disp.Eulipotyphla)
str(disp.df.Eulipotyphla)

disp.df.plot.Eulipotyphla <- melt(disp.df.Eulipotyphla)

disp.df.plot.Eulipotyphla <- disp.df.plot.Eulipotyphla[,-1]
colnames(disp.df.plot.Eulipotyphla) <- c("rolling.bin", "disparity")

# transforming rolling bin axis into body mass mean values
# log value

logmed.rb.Eulipotyphla <- sapply(overl.bins.Eulipotyphla, function(x) median(log10(x$Body_Mass)))
logmed.rb.Eulipotyphla <- format(round(logmed.rb.Eulipotyphla,2), nsmall=2)

logmed.rb.vec.Eulipotyphla <- list()
for(i in 1:5){
  logmed.rb.vec.Eulipotyphla[[i]] <- rep(logmed.rb.Eulipotyphla [i], 10000)
}
logmed.rb.vec.Eulipotyphla <- unlist(logmed.rb.vec.Eulipotyphla)
disp.df.plot.Eulipotyphla$log.med.rb.Eulipotyphla <- as.numeric(logmed.rb.vec.Eulipotyphla)




##=========================================================================================================================================##
#### g) Plot ####
##=========================================================================================================================================##

#### ** FIGURE 2 ** ####


## combining dataframes to plot
combined.df.mammals <- cbind(disp.df.plot[,2:3], taxon="All mammals")
combined.df.Marsupialia <- data.frame(disparity=disp.df.plot.Marsupialia[,2], log.med.rb=disp.df.plot.Marsupialia[,3], taxon="Marsupialia")
combined.df.Cetartiodactyla <- data.frame(disparity=disp.df.plot.Cetartiodactyla[,2], log.med.rb=disp.df.plot.Cetartiodactyla[,3], taxon="Cetartiodactyla")
combined.df.Primates <- data.frame(disparity=disp.df.plot.Primates[,2], log.med.rb=disp.df.plot.Primates[,3], taxon="Primates")
combined.df.Rodentia <- data.frame(disparity=disp.df.plot.Rodentia[,2], log.med.rb=disp.df.plot.Rodentia[,3], taxon="Rodentia")
combined.df.Carnivora <- data.frame(disparity=disp.df.plot.Carnivora[,2], log.med.rb=disp.df.plot.Carnivora[,3], taxon="Carnivora")
combined.df.Eulipotyphla <- data.frame(disparity=disp.df.plot.Eulipotyphla[,2], log.med.rb=disp.df.plot.Eulipotyphla[,3], taxon="Eulipotyphla")

pp.all <- ggplot(combined.df.mammals, aes(x=log.med.rb, y=disparity*10, group=1)) +
  theme_minimal_vgrid()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "gray70", fill="gray70", alpha=0.3, lwd=0.75)+
  xlab("") + ylab("")+
  scale_y_continuous(limits = c(0, 3), breaks=seq(0.5, 2.5, 1.0))+
  scale_x_continuous(limits = c(0.8, 6), breaks=seq(1, 5, 1))+
  theme(axis.text = element_text(size=7),
        axis.line.y = element_line(size=0.5),
        axis.ticks.y = element_line(size=0.5),
        panel.grid.major.x  = element_line(size=0.5))

pp.Marsupialia <- ggplot(combined.df.Marsupialia, aes(x=log.med.rb, y=disparity*10, group=1)) +
  theme_minimal_vgrid()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "#481668FF", fill="#481668FF", alpha=0.3, lwd=0.75)+
  xlab("") + ylab("")+
  scale_y_continuous(limits = c(0, 3), breaks=seq(0.5, 2.5, 1.0))+
  scale_x_continuous(limits = c(0.8, 6), breaks=seq(1, 5, 1))+
  theme(axis.text = element_text(size=7),
       axis.line.y = element_line(size=0.5),
       axis.ticks.y = element_line(size=0.5),
       panel.grid.major.x  = element_line(size=0.5))
      


pp.Cetartiodactyla <- ggplot(combined.df.Cetartiodactyla, aes(x=log.med.rb, y=disparity*10, group=1)) +
  theme_minimal_vgrid()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "#4EC36BFF", fill="#4EC36BFF", alpha=0.3, lwd=0.75)+
  xlab("") + ylab("")+
  scale_y_continuous(limits = c(0, 3), breaks=seq(0.5, 2.5, 1.0))+
  scale_x_continuous(limits = c(0.8, 6), breaks=seq(1, 5, 1))+
  theme(axis.text = element_text(size=7),
        axis.line.y = element_line(size=0.5),
        axis.ticks.y = element_line(size=0.5),
        panel.grid.major.x  = element_line(size=0.5))

pp.Primates <- ggplot(combined.df.Primates, aes(x=log.med.rb, y=disparity*10, group=1)) +
  theme_minimal_vgrid()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "#B4DE2CFF" , fill="#B4DE2CFF" , alpha=0.3, lwd=0.75)+
  xlab("") + ylab("")+
  scale_y_continuous(limits = c(0, 3), breaks=seq(0.5, 2.5, 1.0))+
  scale_x_continuous(limits = c(0.8, 6), breaks=seq(1, 5, 1))+
  theme(axis.text = element_text(size=7),
        axis.line.y = element_line(size=0.5),
        axis.ticks.y = element_line(size=0.5),
        panel.grid.major.x  = element_line(size=0.5))

pp.Rodentia <- ggplot(combined.df.Rodentia, aes(x=log.med.rb, y=disparity*10, group=1)) +
  theme_minimal_vgrid()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "#FFD400" , fill="#FFD400" , alpha=0.3, lwd=0.75)+
  xlab("") + ylab("")+
  scale_y_continuous(limits = c(0, 3), breaks=seq(0.5, 2.5, 1.0))+
  scale_x_continuous(limits = c(0.8, 6), breaks=seq(1, 5, 1))+
  theme(axis.text = element_text(size=7),
        axis.line.y = element_line(size=0.5),
        axis.ticks.y = element_line(size=0.5),
        panel.grid.major.x  = element_line(size=0.5))

pp.Carnivora <- ggplot(combined.df.Carnivora, aes(x=log.med.rb, y=disparity*10, group=1)) +
  theme_minimal_vgrid()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color =  "#25AB82FF", fill= "#25AB82FF", alpha=0.3, lwd=0.75)+
  xlab("") + ylab("")+
  scale_y_continuous(limits = c(0, 3), breaks=seq(0.5, 2.5, 1.0))+
  scale_x_continuous(limits = c(0.8, 6), breaks=seq(1, 5, 1))+
  theme(axis.text = element_text(size=7),
        axis.line.y = element_line(size=0.5),
        axis.ticks.y = element_line(size=0.5),
        panel.grid.major.x  = element_line(size=0.5))

pp.Eulipotyphla <- ggplot(combined.df.Eulipotyphla, aes(x=log.med.rb, y=disparity*10, group=1)) +
  theme_minimal_vgrid()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color ="#21908CFF", fill="#21908CFF", alpha=0.3, lwd=0.75)+
  xlab("") + ylab("")+
  scale_y_continuous(limits = c(0, 3), breaks=seq(0.5, 2.5, 1.0))+
  scale_x_continuous(limits = c(0.8, 6), breaks=seq(1, 5, 1))+
  theme(axis.text = element_text(size=7),
        axis.line.y = element_line(size=0.5),
        axis.ticks.y = element_line(size=0.5),
        panel.grid.major.x  = element_line(size=0.5))



## create separated plots
aligned_plots.Marsupialia <- align_plots(pp.all, pp.Marsupialia, align="hv", axis="tblr")
p.disp.Marsupialia <- (ggdraw(aligned_plots.Marsupialia[[1]]) + draw_plot(aligned_plots.Marsupialia[[2]]))

aligned_plots.Cetartiodactyla <- align_plots(pp.all, pp.Cetartiodactyla, align="hv", axis="tblr")
p.disp.Cetartiodactyla <- (ggdraw(aligned_plots.Cetartiodactyla[[1]]) + draw_plot(aligned_plots.Cetartiodactyla[[2]]))

aligned_plots.Primates <- align_plots(pp.all, pp.Primates, align="hv", axis="tblr")
p.disp.Primates <- (ggdraw(aligned_plots.Primates[[1]]) + draw_plot(aligned_plots.Primates[[2]]))

aligned_plots.Rodentia <- align_plots(pp.all, pp.Rodentia, align="hv", axis="tblr")
p.disp.Rodentia <- (ggdraw(aligned_plots.Rodentia[[1]]) + draw_plot(aligned_plots.Rodentia[[2]]))

aligned_plots.Carnivora <- align_plots(pp.all, pp.Carnivora, align="hv", axis="tblr")
p.disp.Carnivora <- (ggdraw(aligned_plots.Carnivora[[1]]) + draw_plot(aligned_plots.Carnivora[[2]]))

aligned_plots.Eulipotyphla <- align_plots(pp.all, pp.Eulipotyphla, align="hv", axis="tblr")
p.disp.Eulipotyphla <- (ggdraw(aligned_plots.Eulipotyphla[[1]]) + draw_plot(aligned_plots.Eulipotyphla[[2]]))


plot_subgroups <- plot_grid(p.disp.Marsupialia, p.disp.Cetartiodactyla, p.disp.Primates, p.disp.Rodentia, p.disp.Carnivora, p.disp.Eulipotyphla, 
                            labels=c("Marsupialia", "Cetartiodactyla", "Primates", "Rodentia", "Carnivora", "Eulipotyphla"),  label_size=8)

library(grid)
library(gridExtra)
library(ggeasy)
y.lab <- textGrob("Forelimb disparity",  gp=gpar(col="black", fontsize=11), rot=90)
x.lab <- textGrob("Log10 body mass (g)",  gp=gpar(col="black", fontsize=11))
fig2 <- grid.arrange(arrangeGrob(plot_subgroups, left = y.lab, bottom = x.lab))

fig2



##=========================================================================================================================================##
#### h) correlation tests ####
##=========================================================================================================================================##
#create a list with mean disparity values and median masses

group.list.res <- list(Marsupialia=disp.df.plot.Marsupialia,
                      Rodentia=disp.df.plot.Rodentia, 
                      Primates=disp.df.plot.Primates, 
                      Carnivora=disp.df.plot.Carnivora, 
                      Cetartiodactyla=disp.df.plot.Cetartiodactyla, 
                      Eulipotyphla=disp.df.plot.Eulipotyphla)

to.test.group.list <- list()
for (i in 1: length(group.list.res)){
  colnames(group.list.res[[i]]) <- c("rolling.bin", "disparity", "mass")  
  to.test.group.list[[i]] <- cbind(as.data.frame(group.list.res[[i]] %>% group_by(rolling.bin) %>% summarise(mean(disparity))),
                                (as.data.frame(group.list.res[[i]] %>% group_by(rolling.bin) %>% summarise(mean(mass))))[,2])
  colnames(to.test.group.list[[i]]) <- c("rolling.bin", "disparity", "mass")   
}
names(to.test.group.list) <- clades


#run correlation tests
cor.groups <- lapply(to.test.group.list, function(x) cor.test(x$disparity, x$mass, method="kendall"))
fit.groups <- lapply(to.test.group.list, function(x) lm(x$disparity~x$mass))
lapply(fit.groups, function(x) summary(x))


# now with all disparity values

group.list.res.sim <- list(Marsupialia=combined.df.Marsupialia,
                       Rodentia=combined.df.Rodentia, 
                       Primates=combined.df.Primates, 
                       Carnivora=combined.df.Carnivora, 
                       Cetartiodactyla=combined.df.Cetartiodactyla, 
                       Eulipotyphla=combined.df.Eulipotyphla)


fit.groups2 <- lapply(group.list.res.sim, function(x) lm(x$disparity~x$log.med.rb))
lapply(fit.groups2, function(x) summary(x)$coefficients[8])
lapply(fit.groups2, function(x) summary(x)$adj.r.squared)



##==========================================================================================================================================================================## 
##
#### PART VI - Bone disparity ---- ####
##
##==========================================================================================================================================================================## 

# we checked the influence of body mass on the overall limb disparity. but does it impact each limb segment differently?
#  separate each bone trait into a different object
hum <- lapply(limb, function(x) x[,c(1:3)])
rad <- lapply(limb, function(x) x[,c(4:6)])
met <- lapply(limb, function(x) x[,c(7:9)])
phal <- lapply(limb, function(x) x[,c(10:12)])


# and we will now repeat the residual calculation per bone.
# we are repeating it because to isolate the multivariate nature of correlated trait evolution between all bones. 



##------------------------------------------------------------------------------------------------------------------##
## calculating size residuals
##------------------------------------------------------------------------------------------------------------------##

# Humerus
data_hum <- list()
fit_BM_hum <- list() #because we are using 100 trees, perform 100 regressions of limb and geometric means and retain their residual values for the next analyses
for(i in 1:100){
  data_hum[[i]] = list(hum=hum, lgm=lgm[[i]])
  fit_BM_hum[[i]] <- mvgls(hum[[i]]~lgm, data=data_hum[[i]], pruned_trees[[i]], model="BM", penalty="LASSO", method="LL")
}

resid.list_hum <- lapply(fit_BM_hum, function(x) as.data.frame(x$resid)) #list of 100 residual objects

# Radius
data_rad <- list()
fit_BM_rad <- list() 
for(i in 1:100){
  data_rad[[i]] = list(rad=rad, lgm=lgm[[i]])
  fit_BM_rad[[i]] <- mvgls(rad[[i]]~lgm, data=data_rad[[i]], pruned_trees[[i]], model="BM", penalty="LASSO", method="LL")
}

resid.list_rad <- lapply(fit_BM_rad, function(x) as.data.frame(x$resid)) #list of 100 residual objects


# Metacarpal
data_met <- list()
fit_BM_met <- list() 
for(i in 1:100){
  data_met[[i]] = list(met=met, lgm=lgm[[i]])
  fit_BM_met[[i]] <- mvgls(met[[i]]~lgm, data=data_met[[i]], pruned_trees[[i]], model="BM", penalty="LASSO", method="LL")
}

resid.list_met <- lapply(fit_BM_met, function(x) as.data.frame(x$resid)) #list of 100 residual objects



# Phalanx
data_phal <- list()
fit_BM_phal <- list() 
for(i in 1:100){
  data_phal[[i]] = list(phal=phal, lgm=lgm[[i]])
  fit_BM_phal[[i]] <- mvgls(phal[[i]]~lgm, data=data_phal[[i]], pruned_trees[[i]], model="BM", penalty="LASSO", method="LL")
}

resid.list_phal <- lapply(fit_BM_phal, function(x) as.data.frame(x$resid)) #list of 100 residual objects



##------------------------------------------------------------------------------------------------------------------##
## organizing list
##------------------------------------------------------------------------------------------------------------------##


## next - create a list of lists: outer is each fitted model across 100 trees, inside is their corresponding residuals organized into a list of rolling bins
# inside = match residual rownames with overl.bins bin row names


# humerus
resid.overl.bins_hum <- list()
for (i in 1:length(resid.list_hum)) {
  resid.overl.bins_hum[[i]] <- list()  # Initialize a sublist for each dataframe in resid.list
  
  # Iterate over each dataframe in overl.bins
  for (j in 1:length(overl.bins)) {
    resid.overl.bins_hum[[i]][[j]] <- as.matrix(resid.list_hum[[i]][rownames(overl.bins[[j]]), ])
  }
}


# radius
resid.overl.bins_rad <- list()
for (i in 1:length(resid.list_rad)) {
  resid.overl.bins_rad[[i]] <- list()  # Initialize a sublist for each dataframe in resid.list
  
  # Iterate over each dataframe in overl.bins
  for (j in 1:length(overl.bins)) {
    resid.overl.bins_rad[[i]][[j]] <- as.matrix(resid.list_rad[[i]][rownames(overl.bins[[j]]), ])
  }
}


# metacarpal
resid.overl.bins_met <- list()
for (i in 1:length(resid.list_met)) {
  resid.overl.bins_met[[i]] <- list()  # Initialize a sublist for each dataframe in resid.list
  
  # Iterate over each dataframe in overl.bins
  for (j in 1:length(overl.bins)) {
    resid.overl.bins_met[[i]][[j]] <- as.matrix(resid.list_met[[i]][rownames(overl.bins[[j]]), ])
  }
}


# phalanx
resid.overl.bins_phal <- list()
for (i in 1:length(resid.list_phal)) {
  resid.overl.bins_phal[[i]] <- list()  # Initialize a sublist for each dataframe in resid.list
  
  # Iterate over each dataframe in overl.bins
  for (j in 1:length(overl.bins)) {
    resid.overl.bins_phal[[i]][[j]] <- as.matrix(resid.list_phal[[i]][rownames(overl.bins[[j]]), ])
  }
}



##=========================================================================================================================================##
#### a) humerus disparity ####
##=========================================================================================================================================##


# 1) bootstrap data
boot_resid_hum <- list()
for(i in 1:100){
  boot_resid_hum[[i]] <- lapply(resid.overl.bins_hum[[i]], function(x) boot.matrix(x)) # default is 100 replications (*100 trees = 10,000 replications per bin)
}

# 2) calculate disparity from bootstrapped data
disp_resid_hum <- list()
for(i in 1:100){
  disp_resid_hum[[i]] <- lapply(boot_resid_hum[[i]], function(x) dispRity(x, metric = c(sum, variances), verbose=TRUE))
} 

# 3) extracting results to plot
mat.disp_hum <- matrix(ncol=20, nrow=100) # 20 bins, 100 replications (per tree)
list.disp_hum <- list() # creating a list of 100 matrices, corresponding to the disparity values obtained from the limb shape residuals  over the 100 trees
for(k in 1:100){
  list.disp_hum[[k]] <- list()
  for(i in 1:20){
    mat.disp_hum[,i] <- disp_resid_hum[[k]][[i]]$disparity[[1]][[2]] # taking the disparity from each bootstrap rep
    colnames(mat.disp_hum) <- 1:20
    list.disp_hum[[k]] <- mat.disp_hum
  }
}

disp.df_hum <- do.call(rbind, list.disp_hum) # transforming disparity list in a dataframe (columns are the rolling bins (n = 20), rows are the bootstrap rep (n = 10000))
str(disp.df_hum)

disp.df.plot_hum <- melt(disp.df_hum) #checking if results are making sense
disp.df.plot_hum %>% group_by(Var2) %>% count(Var2) # each rolling bin has 10,000 disparity values, which is simply 100 bootstrap reps for 100 trees!
disp.df.plot_hum$value[1:100]==list.disp_hum[[1]][,1] # checking if first 100 elements of our new dataframe correspond to the disparity of rolling bin 1 from tree 1, and it does
disp.df.plot_hum$value[101:200]==list.disp_hum[[2]][,1] # all is good!!

disp.df.plot_hum <- disp.df.plot_hum[,-1] #removing useless row
colnames(disp.df.plot_hum) <- c("rolling.bin", "disparity")


# plot 1 of rolling bins showing humerus disparity

ggplot(disp.df.plot_hum, aes(x=rolling.bin, y=disparity, group=1)) +
  theme_minimal_vgrid()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "steelblue", fill="steelblue", alpha=0.3)+
  xlab("Rolling bins") + ylab("Humerus shape disparity")+
  theme(axis.title.x = element_text(size=14, vjust = -1),
        axis.title.y = element_text(color = "steelblue", size=14,  vjust = 1, hjust = 1),
        axis.line.y = element_line(color = "steelblue"),
        axis.ticks.y = element_line(color = "steelblue"),
        axis.text.y = element_text(color = "steelblue"),
        axis.text = element_text(size=12))




logmed.rb.vec_hum <- list() # combining these values intio the plot dataframe
for(i in 1:20){
  logmed.rb.vec_hum[[i]] <- rep(logmed.rb [i], 10000)
}
logmed.rb.vec_hum <- unlist(logmed.rb.vec_hum)
disp.df.plot_hum$log.med.rb <- as.numeric(logmed.rb.vec_hum)


kg.med.rb_hum <- list() # same thing
for(i in 1:20){
  kg.med.rb_hum[[i]] <- rep(med.rb [i], 10000)
}
kg.med.rb_hum <- unlist(kg.med.rb_hum)
disp.df.plot_hum$kg.med.rb_hum <- kg.med.rb_hum


# plot scaling rolling bins x axis to their median body mass (in log)
ggplot(disp.df.plot_hum, aes(x=log.med.rb, y=disparity, group=1)) +
  theme_bw()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "steelblue", fill="steelblue", alpha=0.3)+
  xlab("log10 body mass (rolling bin median)") + ylab("Humerus shape disparity")+
  theme(axis.title.x = element_text(size=14, vjust = -1),
        axis.title.y = element_text(color = "steelblue", size=14,  vjust = 1, hjust = 1),
        axis.line.y = element_line(color = "steelblue"),
        axis.ticks.y = element_line(color = "steelblue"),
        axis.text.y = element_text(color = "steelblue"),
        axis.text = element_text(size=12))




##=========================================================================================================================================##
#### b) radius disparity ####
##=========================================================================================================================================##


# 1) bootstrap data
boot_resid_rad <- list()
for(i in 1:100){
  boot_resid_rad[[i]] <- lapply(resid.overl.bins_rad[[i]], function(x) boot.matrix(x)) # default is 100 replications (*100 trees = 10,000 replications per bin)
}

# 2) calculate disparity from bootstrapped data
disp_resid_rad <- list()
for(i in 1:100){
  disp_resid_rad[[i]] <- lapply(boot_resid_rad[[i]], function(x) dispRity(x, metric = c(sum, variances), verbose=TRUE))
} 


# 3) extracting results to plot
mat.disp_rad <- matrix(ncol=20, nrow=100) # 20 bins, 100 replications (per tree)
list.disp_rad <- list() # creating a list of 100 matrices, corresponding to the disparity values obtained from the limb shape residuals  over the 100 trees
for(k in 1:100){
  list.disp_rad[[k]] <- list()
  for(i in 1:20){
    mat.disp_rad[,i] <- disp_resid_rad[[k]][[i]]$disparity[[1]][[2]] # taking the disparity from each bootstrap rep
    colnames(mat.disp_rad) <- 1:20
    list.disp_rad[[k]] <- mat.disp_rad
  }
}

disp.df_rad <- do.call(rbind, list.disp_rad) # transforming disparity list in a dataframe (columns are the rolling bins (n = 20), rows are the bootstrap rep (n = 10000))
str(disp.df_rad)

disp.df.plot_rad <- melt(disp.df_rad) #checking if results are making sense
disp.df.plot_rad %>% group_by(Var2) %>% count(Var2) # each rolling bin has 10,000 disparity values, which is simply 100 bootstrap reps for 100 trees!

disp.df.plot_rad <- disp.df.plot_rad[,-1] #removing useless row
colnames(disp.df.plot_rad) <- c("rolling.bin", "disparity")


# plot 1 of rolling bins showing radius disparity

ggplot(disp.df.plot_rad, aes(x=rolling.bin, y=disparity, group=1)) +
  theme_minimal_vgrid()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "steelblue", fill="steelblue", alpha=0.3)+
  xlab("Rolling bins") + ylab("Radius shape disparity")+
  theme(axis.title.x = element_text(size=14, vjust = -1),
        axis.title.y = element_text(color = "steelblue", size=14,  vjust = 1, hjust = 1),
        axis.line.y = element_line(color = "steelblue"),
        axis.ticks.y = element_line(color = "steelblue"),
        axis.text.y = element_text(color = "steelblue"),
        axis.text = element_text(size=12))




logmed.rb.vec_rad <- list() # combining these values into the plot dataframe
for(i in 1:20){
  logmed.rb.vec_rad[[i]] <- rep(logmed.rb [i], 10000)
}
logmed.rb.vec_rad <- unlist(logmed.rb.vec_rad)
disp.df.plot_rad$log.med.rb <- as.numeric(logmed.rb.vec_rad)


kg.med.rb_rad <- list() # same thing
for(i in 1:20){
  kg.med.rb_rad[[i]] <- rep(med.rb [i], 10000)
}
kg.med.rb_rad <- unlist(kg.med.rb_rad)
disp.df.plot_rad$kg.med.rb_rad <- kg.med.rb_rad


# plot scaling rolling bins x axis to their median body mass (in log)
ggplot(disp.df.plot_rad, aes(x=log.med.rb, y=disparity, group=1)) +
  theme_bw()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "steelblue", fill="steelblue", alpha=0.3)+
  xlab("log10 body mass (rolling bin median)") + ylab("Radius shape disparity")+
  theme(axis.title.x = element_text(size=14, vjust = -1),
        axis.title.y = element_text(color = "steelblue", size=14,  vjust = 1, hjust = 1),
        axis.line.y = element_line(color = "steelblue"),
        axis.ticks.y = element_line(color = "steelblue"),
        axis.text.y = element_text(color = "steelblue"),
        axis.text = element_text(size=12))



##=========================================================================================================================================##
#### c) metacarpal disparity ####
##=========================================================================================================================================##


# 1) bootstrap data
boot_resid_met <- list()
for(i in 1:100){
  boot_resid_met[[i]] <- lapply(resid.overl.bins_met[[i]], function(x) boot.matrix(x)) # default is 100 replications (*100 trees = 10,000 replications per bin)
}

# 2) calculate disparity from bootstrapped data
disp_resid_met <- list()
for(i in 1:100){
  disp_resid_met[[i]] <- lapply(boot_resid_met[[i]], function(x) dispRity(x, metric = c(sum, variances), verbose=TRUE))
} 


# 3) extracting results to plot
mat.disp_met <- matrix(ncol=20, nrow=100) # 20 bins, 100 replications (per tree)
list.disp_met <- list() # creating a list of 100 matrices, corresponding to the disparity values obtained from the limb shape residuals  over the 100 trees
for(k in 1:100){
  list.disp_met[[k]] <- list()
  for(i in 1:20){
    mat.disp_met[,i] <- disp_resid_met[[k]][[i]]$disparity[[1]][[2]] # taking the disparity from each bootstrap rep
    colnames(mat.disp_met) <- 1:20
    list.disp_met[[k]] <- mat.disp_met
  }
}

disp.df_met <- do.call(rbind, list.disp_met) # transforming disparity list in a dataframe (columns are the rolling bins (n = 20), rows are the bootstrap rep (n = 10000))
str(disp.df_met)

disp.df.plot_met <- melt(disp.df_met) #checking if results are making sense
disp.df.plot_met %>% group_by(Var2) %>% count(Var2) # each rolling bin has 10,000 disparity values, which is simply 100 bootstrap reps for 100 trees!

disp.df.plot_met <- disp.df.plot_met[,-1] #removing useless row
colnames(disp.df.plot_met) <- c("rolling.bin", "disparity")


# plot 1 of rolling bins showing metacarpal disparity

ggplot(disp.df.plot_met, aes(x=rolling.bin, y=disparity, group=1)) +
  theme_minimal_vgrid()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "steelblue", fill="steelblue", alpha=0.3)+
  xlab("Rolling bins") + ylab("Metacarpal shape disparity")+
  theme(axis.title.x = element_text(size=14, vjust = -1),
        axis.title.y = element_text(color = "steelblue", size=14,  vjust = 1, hjust = 1),
        axis.line.y = element_line(color = "steelblue"),
        axis.ticks.y = element_line(color = "steelblue"),
        axis.text.y = element_text(color = "steelblue"),
        axis.text = element_text(size=12))


logmed.rb.vec_met <- list() # combining these values into the plot dataframe
for(i in 1:20){
  logmed.rb.vec_met[[i]] <- rep(logmed.rb [i], 10000)
}
logmed.rb.vec_met <- unlist(logmed.rb.vec_met)
disp.df.plot_met$log.med.rb <- as.numeric(logmed.rb.vec_met)


kg.med.rb_met <- list() # same thing
for(i in 1:20){
  kg.med.rb_met[[i]] <- rep(med.rb [i], 10000)
}
kg.med.rb_met <- unlist(kg.med.rb_met)
disp.df.plot_met$kg.med.rb_met <- kg.med.rb_met


# plot scaling rolling bins x axis to their median body mass (in log)
ggplot(disp.df.plot_met, aes(x=log.med.rb, y=disparity, group=1)) +
  theme_bw()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "steelblue", fill="steelblue", alpha=0.3)+
  xlab("log10 body mass (rolling bin median)") + ylab("Metacarpal shape disparity")+
  theme(axis.title.x = element_text(size=14, vjust = -1),
        axis.title.y = element_text(color = "steelblue", size=14,  vjust = 1, hjust = 1),
        axis.line.y = element_line(color = "steelblue"),
        axis.ticks.y = element_line(color = "steelblue"),
        axis.text.y = element_text(color = "steelblue"),
        axis.text = element_text(size=12))



##=========================================================================================================================================##
#### d) phalanx disparity ####
##=========================================================================================================================================##


# 1) bootstrap data
boot_resid_phal <- list()
for(i in 1:100){
  boot_resid_phal[[i]] <- lapply(resid.overl.bins_phal[[i]], function(x) boot.matrix(x)) # default is 100 replications (*100 trees = 10,000 replications per bin)
}

# 2) calculate disparity from bootstrapped data
disp_resid_phal <- list()
for(i in 1:100){
  disp_resid_phal[[i]] <- lapply(boot_resid_phal[[i]], function(x) dispRity(x, metric = c(sum, variances), verbose=TRUE))
} 


# 3) extracting results to plot
mat.disp_phal <- matrix(ncol=20, nrow=100) # 20 bins, 100 replications (per tree)
list.disp_phal <- list() # creating a list of 100 matrices, corresponding to the disparity values obtained from the limb shape residuals  over the 100 trees
for(k in 1:100){
  list.disp_phal[[k]] <- list()
  for(i in 1:20){
    mat.disp_phal[,i] <- disp_resid_phal[[k]][[i]]$disparity[[1]][[2]] # taking the disparity from each bootstrap rep
    colnames(mat.disp_phal) <- 1:20
    list.disp_phal[[k]] <- mat.disp_phal
  }
}

disp.df_phal <- do.call(rbind, list.disp_phal) # transforming disparity list in a dataframe (columns are the rolling bins (n = 20), rows are the bootstrap rep (n = 10000))
str(disp.df_phal)

disp.df.plot_phal <- melt(disp.df_phal) #checking if results are making sense
disp.df.plot_phal %>% group_by(Var2) %>% count(Var2) # each rolling bin has 10,000 disparity values, which is simply 100 bootstrap reps for 100 trees!

disp.df.plot_phal <- disp.df.plot_phal[,-1] #removing useless row
colnames(disp.df.plot_phal) <- c("rolling.bin", "disparity")


# plot 1 of rolling bins showing phalanx disparity

ggplot(disp.df.plot_phal, aes(x=rolling.bin, y=disparity, group=1)) +
  theme_minimal_vgrid()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "steelblue", fill="steelblue", alpha=0.3)+
  xlab("Rolling bins") + ylab("Phalanx shape disparity")+
  theme(axis.title.x = element_text(size=14, vjust = -1),
        axis.title.y = element_text(color = "steelblue", size=14,  vjust = 1, hjust = 1),
        axis.line.y = element_line(color = "steelblue"),
        axis.ticks.y = element_line(color = "steelblue"),
        axis.text.y = element_text(color = "steelblue"),
        axis.text = element_text(size=12))




logmed.rb.vec_phal <- list() # combining these values intio the plot dataframe
for(i in 1:20){
  logmed.rb.vec_phal[[i]] <- rep(logmed.rb [i], 10000)
}
logmed.rb.vec_phal <- unlist(logmed.rb.vec_phal)
disp.df.plot_phal$log.med.rb <- as.numeric(logmed.rb.vec_phal)


kg.med.rb_phal <- list() # same thing
for(i in 1:20){
  kg.med.rb_phal[[i]] <- rep(med.rb [i], 10000)
}
kg.med.rb_phal <- unlist(kg.med.rb_phal)
disp.df.plot_phal$kg.med.rb_phal <- kg.med.rb_phal


# plot scaling rolling bins x axis to their median body mass (in log)
ggplot(disp.df.plot_phal, aes(x=log.med.rb, y=disparity, group=1)) +
  theme_bw()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "steelblue", fill="steelblue", alpha=0.3)+
  xlab("log10 body mass (rolling bin median)") + ylab("Phalanx shape disparity")+
  theme(axis.title.x = element_text(size=14, vjust = -1),
        axis.title.y = element_text(color = "steelblue", size=14,  vjust = 1, hjust = 1),
        axis.line.y = element_line(color = "steelblue"),
        axis.ticks.y = element_line(color = "steelblue"),
        axis.text.y = element_text(color = "steelblue"),
        axis.text = element_text(size=12))



##=========================================================================================================================================##
#### e) Plot ####
##=========================================================================================================================================##

#### ** FIGURE 3 ** ####

# to compare the disparity between bones, let's set scale them all in the sam y axis

phum <- ggplot(disp.df.plot_hum, aes(x=log.med.rb, y=disparity*10, group=1)) + # multiplying disparity by 10 to facilitate visualization
  theme_classic()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "indianred1", fill="indianred1", alpha=0.3)+
  xlab("") + ylab("")+
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0.3, 1.2, 0.3))+
  theme(axis.text.y = element_text(color = "gray33"),
        axis.text = element_text(size=8))


prad <- ggplot(disp.df.plot_rad, aes(x=log.med.rb, y=disparity*10, group=1)) +
  theme_classic()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "indianred1", fill="indianred1", alpha=0.3)+
  xlab("") + ylab("")+
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0.3, 1.2, 0.3))+
  theme(axis.text.y = element_text(color = "gray33"),
        axis.text = element_text(size=8))



pmet <- ggplot(disp.df.plot_met, aes(x=log.med.rb, y=disparity*10, group=1)) +
  theme_classic()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "indianred1", fill="indianred1", alpha=0.3)+
  xlab("") + ylab("")+
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0.3, 1.2, 0.3))+
  theme(axis.text.y = element_text(color = "gray33"),
        axis.text = element_text(size=8))



pphal <- ggplot(disp.df.plot_phal, aes(x=log.med.rb, y=disparity*10, group=1)) +
  theme_classic()+
  stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1),
               geom = 'smooth', se = TRUE, color = "indianred1", fill="indianred1", alpha=0.3)+
  xlab("") + ylab("")+
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0.3, 1.2, 0.3))+
  theme(axis.text.y = element_text(color = "gray33"),
        axis.text = element_text(size=8))


library(grid)
library(gridExtra)
library(ggeasy)

pgrid <- plot_grid(phum, prad, pmet, pphal, ncol=2, labels=c("humerus", "radius", "metacarpal", "phalanx"), align = "hv", label_x = 0.15, label_y = 1,label_size = 10, scale=1)
y.lab<- textGrob("Bone shape disparity",  gp=gpar(col="black", fontsize=11), rot=90)
x.lab<- textGrob("log10 mass (g)",  gp=gpar(col="black", fontsize=11))
fig3 <- grid.arrange(arrangeGrob(pgrid, left = y.lab, bottom = x.lab))

fig3


##=========================================================================================================================================##
#### f) Linear models ####
##=========================================================================================================================================##


# to make the datasets comparable, we will use the mean values of the calculated variables per means 
#(except for body mass, that we will use median)

# disparity
disp.means_hum <- as.data.frame(disp.df.plot_hum %>% group_by(rolling.bin) %>% summarise(mean(disparity)))
disp.means_rad <- as.data.frame(disp.df.plot_rad %>% group_by(rolling.bin) %>% summarise(mean(disparity)))
disp.means_met <- as.data.frame(disp.df.plot_met %>% group_by(rolling.bin) %>% summarise(mean(disparity)))
disp.means_phal <- as.data.frame(disp.df.plot_phal %>% group_by(rolling.bin) %>% summarise(mean(disparity)))


# faith distance
faith.means <- as.data.frame(faithPD20.df %>% group_by(bin) %>% summarise(mean(faith.dist)))

# body mass
median.bodymass <- as.numeric(logmed.rb)

# gathering data into separated dataframes:
df.mlr.hum <- data.frame(rolling=c(1:20), mass=median.bodymass, 
                     disparity=disp.means_hum$`mean(disparity)`, 
                     faithPD=faith.means$`mean(faith.dist)`, gower= df.gower.means$gower.means)

df.mlr.rad <- data.frame(rolling=c(1:20), mass=median.bodymass, 
                         disparity=disp.means_rad$`mean(disparity)`, 
                         faithPD=faith.means$`mean(faith.dist)`, gower= df.gower.means$gower.means)

df.mlr.met <- data.frame(rolling=c(1:20), mass=median.bodymass, 
                         disparity=disp.means_met$`mean(disparity)`, 
                         faithPD=faith.means$`mean(faith.dist)`, gower= df.gower.means$gower.means)

df.mlr.phal <- data.frame(rolling=c(1:20), mass=median.bodymass, 
                         disparity=disp.means_phal$`mean(disparity)`, 
                         faithPD=faith.means$`mean(faith.dist)`, gower= df.gower.means$gower.means)



fit.mlr_hum <- lm(disparity~ mass + faithPD + gower, data=df.mlr.hum)
fit.mlr_rad <- lm(disparity~ mass + faithPD + gower, data=df.mlr.rad)
fit.mlr_met <- lm(disparity~ mass + faithPD + gower, data=df.mlr.met)
fit.mlr_phal <- lm(disparity~ mass + faithPD + gower, data=df.mlr.phal)

summary(fit.mlr_hum)
summary(fit.mlr_rad)
summary(fit.mlr_met)
summary(fit.mlr_phal)





##=========================================================================================================================================##
#### g) Kendall tests ####
##=========================================================================================================================================##

# we will rerun the kendall tests for each bone

param <- c("disparity", "mass", "gower", "faithPD") # names matching the col names of the data frame
pair.param <- t(combn(param,2))

nam.mat <- c("Disparity", "Body mass", "Ecological diversity", "Faith's distance") # better looking names for the plots later on


## humerus
k.hum <-list()
for (i in 1: nrow(pair.param)){
list.hum <- as.list(df.mlr.hum)
k.hum[[i]] <- cor.test(list.hum[[pair.param[i,1]]], list.hum[[pair.param[i,2]]], method="kendall") 
}

tau.hum <- sapply(k.hum, function(x) x$estimate) 
pvalue.hum <- sapply(k.hum, function(x) x$p.value)

  #tau matrix
mat.tau.hum <- matrix(ncol=4, nrow=4, 
                  c(1,           tau.hum[1],     tau.hum[2],     tau.hum[3],   # disparity
                    tau.hum[1],       1,         tau.hum[4],     tau.hum[5], # bodymass  
                    tau.hum[2],      tau.hum[4],       1,        tau.hum[6], # ecology
                    tau.hum[3],      tau.hum[5],      tau.hum[6],      1)) # faith dist

colnames(mat.tau.hum) <- nam.mat
rownames(mat.tau.hum) <- nam.mat


  # p value matrix
mat.p.hum <- matrix(ncol=4, nrow=4,
                c(1,           pvalue.hum[1],     pvalue.hum[2],     pvalue.hum[3],   # disparity
                  pvalue.hum[1],       1,         pvalue.hum[4],     pvalue.hum[5], # bodymass  
                  pvalue.hum[2],      pvalue.hum[4],       1,        pvalue.hum[6], # ecology
                  pvalue.hum[3],      pvalue.hum[5],      pvalue.hum[6],      1))# faith dist

colnames(mat.p.hum) <- nam.mat



## radius
k.rad <-list()
for (i in 1: nrow(pair.param)){
  list.rad <- as.list(df.mlr.rad)
  k.rad[[i]] <- cor.test(list.rad[[pair.param[i,1]]], list.rad[[pair.param[i,2]]], method="kendall") 
}

tau.rad <- sapply(k.rad, function(x) x$estimate) 
pvalue.rad <- sapply(k.rad, function(x) x$p.value)


  #tau matrix
mat.tau.rad <- matrix(ncol=4, nrow=4, 
                      c(1,           tau.rad[1],     tau.rad[2],     tau.rad[3],   # disparity
                        tau.rad[1],       1,         tau.rad[4],     tau.rad[5], # bodymass  
                        tau.rad[2],      tau.rad[4],       1,        tau.rad[6], # ecology
                        tau.rad[3],      tau.rad[5],      tau.rad[6],      1)) # faith dist

colnames(mat.tau.rad) <- nam.mat
rownames(mat.tau.rad) <- nam.mat


  # p value matrix
mat.p.rad <- matrix(ncol=4, nrow=4,
                    c(1,           pvalue.rad[1],     pvalue.rad[2],     pvalue.rad[3],   # disparity
                      pvalue.rad[1],       1,         pvalue.rad[4],     pvalue.rad[5], # bodymass  
                      pvalue.rad[2],      pvalue.rad[4],       1,        pvalue.rad[6], # ecology
                      pvalue.rad[3],      pvalue.rad[5],      pvalue.rad[6],      1))# faith dist

colnames(mat.p.rad) <- nam.mat



## metacarpal
k.met <-list()
for (i in 1: nrow(pair.param)){
  list.met <- as.list(df.mlr.met)
  k.met[[i]] <- cor.test(list.met[[pair.param[i,1]]], list.met[[pair.param[i,2]]], method="kendall") 
}

tau.met <- sapply(k.met, function(x) x$estimate) 
pvalue.met <- sapply(k.met, function(x) x$p.value)


  #tau matrix
mat.tau.met <- matrix(ncol=4, nrow=4, 
                      c(1,           tau.met[1],     tau.met[2],     tau.met[3],   # disparity
                        tau.met[1],       1,         tau.met[4],     tau.met[5], # bodymass  
                        tau.met[2],      tau.met[4],       1,        tau.met[6], # ecology
                        tau.met[3],      tau.met[5],      tau.met[6],      1)) # faith dist

colnames(mat.tau.met) <- nam.mat
rownames(mat.tau.met) <- nam.mat


  # p value matrix
mat.p.met <- matrix(ncol=4, nrow=4,
                    c(1,           pvalue.met[1],     pvalue.met[2],     pvalue.met[3],   # disparity
                      pvalue.met[1],       1,         pvalue.met[4],     pvalue.met[5], # bodymass  
                      pvalue.met[2],      pvalue.met[4],       1,        pvalue.met[6], # ecology
                      pvalue.met[3],      pvalue.met[5],      pvalue.met[6],      1))# faith dist

colnames(mat.p.met) <- nam.mat



## phalanx
k.phal <-list()
for (i in 1: nrow(pair.param)){
  list.phal <- as.list(df.mlr.phal)
  k.phal[[i]] <- cor.test(list.phal[[pair.param[i,1]]], list.phal[[pair.param[i,2]]], method="kendall") 
}

tau.phal <- sapply(k.phal, function(x) x$estimate) 
pvalue.phal <- sapply(k.phal, function(x) x$p.value)


  #tau matrix
mat.tau.phal <- matrix(ncol=4, nrow=4, 
                      c(1,           tau.phal[1],     tau.phal[2],     tau.phal[3],   # disparity
                        tau.phal[1],       1,         tau.phal[4],     tau.phal[5], # bodymass  
                        tau.phal[2],      tau.phal[4],       1,        tau.phal[6], # ecology
                        tau.phal[3],      tau.phal[5],      tau.phal[6],      1)) # faith dist

colnames(mat.tau.phal) <- nam.mat
rownames(mat.tau.phal) <- nam.mat


  # p value matrix
mat.p.phal <- matrix(ncol=4, nrow=4,
                    c(1,           pvalue.phal[1],     pvalue.phal[2],     pvalue.phal[3],   # disparity
                      pvalue.phal[1],       1,         pvalue.phal[4],     pvalue.phal[5], # bodymass  
                      pvalue.phal[2],      pvalue.phal[4],       1,        pvalue.phal[6], # ecology
                      pvalue.phal[3],      pvalue.phal[5],      pvalue.phal[6],      1))# faith dist

colnames(mat.p.phal) <- nam.mat



## plots

pcorhum <- ggcorrplot(mat.tau.hum, p.mat=mat.p.hum, lab=TRUE, type = "lower", insig="blank", legend.title = "Correlation (p < 0.05)", 
                    title="humerus", lab_size=3, tl.cex = 10 )


pcorrad <- ggcorrplot(mat.tau.rad, p.mat=mat.p.rad, lab=TRUE, type = "lower", insig="blank", legend.title = "Correlation (p < 0.05)", 
                      title="radius", lab_size=3, tl.cex = 10 )


pcormet <- ggcorrplot(mat.tau.met, p.mat=mat.p.met, lab=TRUE, type = "lower", insig="blank", legend.title = "Correlation (p < 0.05)", 
                      title="metacarpal", lab_size=3, tl.cex = 10 )


pcorphal <- ggcorrplot(mat.tau.phal, p.mat=mat.p.phal, lab=TRUE, type = "lower", insig="blank", legend.title = "Correlation (p < 0.05)", 
                      title="phalanx", lab_size=3, tl.cex = 10 )


plot_grid(pcorhum, pcorrad, pcormet, pcorphal)


#### ** Plots for the paper summarized ** ####

fig1 <- plot_grid(main_plot, bar_plots, ncol=1, rel_heights = c(1.5,1),  align="hv")
fig2 <- grid.arrange(arrangeGrob(plot_subgroups, left = y.lab, bottom = x.lab))
fig3 <- grid.arrange(arrangeGrob(pgrid, left = y.lab, bottom = x.lab))




##=========================================================================================================================================##
#### SUPPORTING INFORMATION ####
##=========================================================================================================================================##

#### ** FIGURE SI ** ####

pca <- prcomp(resid.sample[1:12])
summary(pca)
scores <- as.data.frame(pca$x)
scores$Group <- mydat.list[[1]]$Group
scores$Group <- recode_factor(scores$Group, Cetartiodactyla_Non_Cet="Cetartiodactyla")
names(colsg) <- as.character(recode_factor(names(colsg), Cetartiodactyla_Non_Cet="Cetartiodactyla"))
scores$Locomotion <- as.character(mydat.list[[1]]$Locomotion)

pca_clade <- ggplot(scores, aes(x=PC1, y=PC2, color=Group))+
  geom_point(aes(color=Group), size=2.5, alpha=0.5) +
  scale_color_manual(values=colsg)+
  theme_classic()+
  ylab("")+
  xlab("")+
 # theme(legend.position = "none")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed")

pca_eco <- ggplot(scores, aes(x=PC1, y=PC2, color=Locomotion))+
    geom_point(size=2.5, alpha=0.5) +
    theme_classic()+
    scale_color_manual(values= c("forestgreen", "darkred", "black", "olivedrab1", "skyblue", "red", "violet", "#FFC055"),
                       name= "Locomotion", breaks=c("arboreal", "fossorial", "gliding", "scansorial", "semiaquatic", "semifossorial", "terrestrialbip", "terrestrialquad"), 
                       labels = c("Arboreal", "Fossorial", "Gliding", "Scansorial", "Semiaquatic", "Semifossorial", "Terr. Bip", "Terr. Quad"))+
    ylab("")+
    xlab("")+
 # theme(legend.position = "none")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed")


pca.plot <- plot_grid(pca_clade, pca_eco)
y.lab.p1 <- textGrob("PC2 (19.5 %)",  gp=gpar(col="black", fontsize=12, type="bold"), rot=90)
x.lab.p1 <- textGrob("PC1 (58.0 %)",  gp=gpar(col="black", fontsize=12, type="bold"))
fig.extra <- grid.arrange(arrangeGrob(pca.plot , left = y.lab.p1, bottom = x.lab.p1))





# create one list 20 dataframes, that are repeating the scores from pca. each df will correspond to a different rolling bin, that will be specified in the "Bin" column
# add two columns. a) $bin (two factors = "TRUE" (if species corresponds to bin x), "FALSE") and b) $locomotion (with "na" to species that are not in the bin)
# increase dot size of corresponding bin. color $locomotion 

list.scores.bins <- list()

for(i in 1:20){
  list.scores.bins[[i]] <- scores
  list.scores.bins[[i]]$Bin <- rownames(list.scores.bins[[i]]) %in% rownames(overl.bins[[i]])
  list.scores.bins[[i]]$Locomotion[list.scores.bins[[i]]$Bin==FALSE] <- "NA"
  
}
test <- list.scores.bins[[2]]
cols.loc2 <- setNames(c("forestgreen", "darkred", "black","olivedrab1", "skyblue", "red", "violet", "#FFC055", "#E8E8E8"),
                      c("arboreal","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad", "NA"))

# compact function:

for (i in 1:20){ 
  plot <- eval(substitute(
    ggplot(list.scores.bins[[i]], aes(x = PC1, y = PC2, size=Bin))+
      geom_point(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == FALSE,], aes (color=Locomotion)) +
      # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
      geom_point(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
      scale_color_manual(values=cols.loc2)+
      scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
      theme_classic() +
      xlab("") +
      ylab("")+
      theme(legend.position = "none"),
    list(i = i)))
  
  plots[[i]] <- plot
}

plot_grid(plotlist = plots, ncol=4, nrow=5)


# plot by plot to custom presence/absence of axis element:

{
  pb1 <- ggplot(list.scores.bins[[1]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[1]][list.scores.bins[[1]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[1]][list.scores.bins[[1]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[1]  2.3 - 16.9 g")+
    theme(axis.text.x= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  
  pb2 <- ggplot(list.scores.bins[[2]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[2]][list.scores.bins[[2]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[2]][list.scores.bins[[2]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[2]  11.3 - 24 g")+
    theme(axis.text.x= element_blank(),
          axis.text.y= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  pb3 <- ggplot(list.scores.bins[[3]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[3]][list.scores.bins[[3]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[3]][list.scores.bins[[3]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[3]  21.4 - 38.8 g")+
    theme(axis.text.x= element_blank(),
          axis.text.y= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  pb4 <- ggplot(list.scores.bins[[4]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[4]][list.scores.bins[[4]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[4]][list.scores.bins[[4]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[4]  33.5 - 53.3 g")+
    theme(axis.text.x= element_blank(),
          axis.text.y= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  
  pb5 <- ggplot(list.scores.bins[[5]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[5]][list.scores.bins[[5]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[5]][list.scores.bins[[5]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[5]  49 - 69.8 g")+
    theme(axis.text.x= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  
  pb6 <- ggplot(list.scores.bins[[6]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[6]][list.scores.bins[[6]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[6]][list.scores.bins[[6]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[6]  63.2 - 100 g")+
    theme(axis.text.x= element_blank(),
          axis.text.y= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  pb7 <- ggplot(list.scores.bins[[7]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[7]][list.scores.bins[[7]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[7]][list.scores.bins[[7]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[7]  89.8 - 152.2 g")+
    theme(axis.text.x= element_blank(),
          axis.text.y= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  pb8 <- ggplot(list.scores.bins[[8]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[8]][list.scores.bins[[8]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[8]][list.scores.bins[[8]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[8]  138.2 - 217 g")+
    theme(axis.text.x= element_blank(),
          axis.text.y= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  
  
  pb9 <- ggplot(list.scores.bins[[9]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[9]][list.scores.bins[[9]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[9]][list.scores.bins[[9]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[9]  200 - 326.7 g")+
    theme(axis.text.x= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  
  pb10 <- ggplot(list.scores.bins[[10]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[10]][list.scores.bins[[10]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[10]][list.scores.bins[[10]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[10]  284.5 - 558 g")+
    theme(axis.text.x= element_blank(),
          axis.text.y= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  pb11 <- ggplot(list.scores.bins[[11]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[11]][list.scores.bins[[11]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[11]][list.scores.bins[[11]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[11]  445 - 916 g")+
    theme(axis.text.x= element_blank(),
          axis.text.y= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  pb12 <- ggplot(list.scores.bins[[12]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[12]][list.scores.bins[[12]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[12]][list.scores.bins[[12]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[12]  0.8 - 1.7 kg")+
    theme(axis.text.x= element_blank(),
          axis.text.y= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  
  pb13 <- ggplot(list.scores.bins[[13]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[13]][list.scores.bins[[13]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[13]][list.scores.bins[[13]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[13]  1.4 - 2.9 kg")+
    theme(axis.text.x= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  
  pb14 <- ggplot(list.scores.bins[[14]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[14]][list.scores.bins[[14]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[14]][list.scores.bins[[14]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[14]  2.6 - 4.5 kg")+
    theme(axis.text.x= element_blank(),
          axis.text.y= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  pb15 <- ggplot(list.scores.bins[[15]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[15]][list.scores.bins[[15]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[15]][list.scores.bins[[15]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[15]  4.1 - 7.7 kg")+
    theme(axis.text.x= element_blank(),
          axis.text.y= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  pb16 <- ggplot(list.scores.bins[[16]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[16]][list.scores.bins[[16]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[16]][list.scores.bins[[16]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[16]  6.5 - 12.7 kg")+
    theme(axis.text.x= element_blank(),
          axis.text.y= element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  
  pb17 <- ggplot(list.scores.bins[[17]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[17]][list.scores.bins[[17]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[17]][list.scores.bins[[17]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[17]  10.9 - 29.5 kg")+
    theme(plot.title = element_text(size=10),
          panel.border = element_rect(fill = NA, color = "black"), 
          legend.position="none")
  
  
  pb18 <- ggplot(list.scores.bins[[18]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[18]][list.scores.bins[[18]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[18]][list.scores.bins[[18]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[18]  22 - 69.5 kg")+
    theme(#axis.text.x= element_blank(),
      axis.text.y= element_blank(),
      plot.title = element_text(size=10),
      panel.border = element_rect(fill = NA, color = "black"), 
      legend.position="none")
  
  pb19 <- ggplot(list.scores.bins[[19]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[19]][list.scores.bins[[19]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[19]][list.scores.bins[[19]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[19]  52.6 - 206.1 kg")+
    theme(#axis.text.x= element_blank(),
      axis.text.y= element_blank(),
      plot.title = element_text(size=10),
      panel.border = element_rect(fill = NA, color = "black"), 
      legend.position="none")
  
  pb20 <- ggplot(list.scores.bins[[20]], aes(x = PC1, y = PC2, size=Bin))+
    geom_point(data=list.scores.bins[[20]][list.scores.bins[[20]]$Bin == FALSE,], aes (color=Locomotion)) +
    # geom_polygon(data=list.scores.bins[[i]][list.scores.bins[[i]]$Bin == TRUE,], stat = "ellipse", size=0.5, fill="black", alpha=0.2)+
    geom_point(data=list.scores.bins[[20]][list.scores.bins[[20]]$Bin == TRUE,], aes (color=Locomotion), alpha=0.6) +
    scale_color_manual(values=cols.loc2)+
    scale_size_manual(values = setNames(c(1.3,2.2),c("FALSE", "TRUE")))+
    theme_classic() +
    xlab("") +
    ylab("")+
    ggtitle("[20]  123.2 - 3824.5 kg")+
    theme(#axis.text.x= element_blank(),
      axis.text.y= element_blank(),
      plot.title = element_text(size=10),
      panel.border = element_rect(fill = NA, color = "black"), 
      legend.position="none")
  
}

plots.bins <- plot_grid(pb1, pb2, pb3, pb4, 
                        pb5, pb6, pb7, pb8, 
                        pb9, pb10, pb11, pb12,
                        pb13, pb14, pb15, pb16,
                        pb17, pb18, pb19, pb20,
                        nrow=5, ncol=4, align="hv", scale=1)

summary(pca)

library(grid)
library(gridExtra)


y.lab.p1 <- textGrob("PC2 (19.5 %)",  gp=gpar(col="black", fontsize=12, type="bold"), rot=90)
x.lab.p1 <- textGrob("PC1 (58.0 %)",  gp=gpar(col="black", fontsize=12, type="bold"))
fig.S3 <- grid.arrange(arrangeGrob(plots.bins, left = y.lab.p1, bottom = x.lab.p1))

leg_pc_plot <- get_legend(
  ggplot(scores, aes(x=PC1, y=PC2, color=Locomotion))+
  geom_point(size=3) +
  theme_classic()+
  scale_color_manual(values= c("forestgreen", "darkred", "black", "olivedrab1", "skyblue", "red", "violet", "#FFC055"),
      name= "", breaks=c("arboreal", "fossorial", "gliding", "scansorial", "semiaquatic", "semifossorial", "terrestrialbip", "terrestrialquad"), 
      labels = c("Arboreal", "Fossorial", "Gliding", "Scansorial", "Semi-aquatic", "Semi-fossorial", "Terr. Bip", "Terr. Quad"))+
  theme(legend.position = "bottom")
)

plot_grid(fig.S3, leg_pc_plot, ncol=1, rel_heights = c(5,2))
