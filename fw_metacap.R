# required packages
library(ggplot2)
library(dplyr)
library(ggforce)
library(readxl)
library(ggpubr)
library(ggrepel)
library(metacapa)
library(RColorBrewer)
library(tidyr)
library(DescTools)

# specify site (Site1 or Site2)
site = "Site1" 

# metapopulation parameters
alpha = 4      # dispersal measure from occupancy model
top = 12       # number of top sites to look at restoring
hec = 10000     # scaling for m2 area into hectares 
ex = 0.1       # scaling extinction with patch area (smaller value = less difference between sites)
prop = seq(0.1, 1, 0.1) # proportion of area restored in each site

E = 0.1 # median extinction rate from Martin et al. 2025
C = 0.011 # median colonization rate from Martin et al. 2025

Ndf = read_xlsx(paste(site, "_wetland_attributes.xlsx", sep = ""))
Ndm = read_xlsx(paste(site, "_Distance_Matrix.xlsx", sep = ""))

# format data frame and extract site IDs
Ndm = as.data.frame(Ndm)
rownames(Ndm) = as.character(Ndm[,1])
Ndm = Ndm[,-1]

# split sites with and without habitat
Ndf.extant = Ndf[which(Ndf$Suitable_Veg_Area_m2 >= 100),]
Ndf.potential = Ndf[which(Ndf$Suitable_Veg_Area_m2 < 100),]

###############
# habitat stats

quantile(Ndf$Area_m2/hec, probs = c(0, 0.5, 1))
sum(Ndf$Area_m2/hec)

quantile(Ndf.extant$Suitable_Veg_Area_m2/hec, probs = c(0, 1))
quantile(Ndf.extant$Suitable_Veg_Area_m2/Ndf.extant$Area_m2, probs = c(0, 1))
sum(Ndf.extant$Suitable_Veg_Area_m2/hec)

###############

# scale distances to km
Ndm = Ndm/1000

##################
## distance stats

# all sites
ds = unlist(Ndm)
ds = ds[ds != 0] # remove self connections
ds = ds[ds < 0.2] # remove distant connections
ds = ds[!duplicated(ds)] 

# current
d.ex = Ndm[as.character(Ndf.extant$Wetland_ID), as.character(Ndf.extant$Wetland_ID)]
ds = unlist(d.ex)
ds = ds[ds != 0]
ds = ds[ds < 0.2] # remove distant connections
ds = ds[!duplicated(ds)] 

###########

# potential landscape capacity, current landscape capacity, new sites capacity 
sol.potential1 <- meta_capacity(as.dist(d.ex), Ndf.extant$Area_m2/hec, f = dispersal_negexp(alpha), patch_mc = T, ex = ex)
sol.potential2 <- meta_capacity(as.dist(Ndm), Ndf$Area_m2/hec, f = dispersal_negexp(alpha), patch_mc = T, ex = ex)
sol.current <- meta_capacity(as.dist(d.ex), Ndf.extant$Suitable_Veg_Area_m2/hec, f = dispersal_negexp(alpha), patch_mc = T, ex = ex)

tops1 = Ndf[order(sol.potential1$patch_mc), "RowID"]$RowID
tops2 = Ndf[order(sol.potential2$patch_mc), "RowID"]$RowID

top.existing = Ndf[tops1[tops1 %in% Ndf.extant$RowID][1:top],"RowID"]
top.potential = Ndf[tops2[tops2 %in% Ndf.potential$RowID][1:top],"RowID"]

# threshold for metapopulation viability
cap.threshold = 100 * (E/C - sol.current$capacity) / sol.current$capacity


#############
# simulations
#############

# simulate restoration scenarios at existing patches
df = tibble()

count=1

for(j in 1:length(top.potential$RowID)) {
  
  s = combn(top.potential$RowID, j)
  
  for(i in 1:length(s[1,])){
    
    new.sites = Ndf[which(Ndf$RowID %in% s[,i]),]
    
    for(t in 1:length(prop)){
      
      df[count,1] = list(list(new.sites[,"Wetland_ID"]$Wetland_ID))
      
      df.temp = Ndf
      rs = prop[t] * new.sites$Area_m2
      df.temp[which(df.temp$RowID %in% s[,i]), "Suitable_Veg_Area_m2"] = rs
      
      As = df.temp[which(df.temp$Suitable_Veg_Area_m2 >= 100),]
      
      d = Ndm
      d = d[as.character(As$Wetland_ID), as.character(As$Wetland_ID)]
      
      sol <- meta_capacity(as.dist(d), As$Suitable_Veg_Area_m2/hec, f = dispersal_negexp(alpha), patch_mc = T, ex = ex)
      
      df[count,2] = (sol$capacity - sol.current$capacity) / sol.current$capacity
      df[count,3] = sum(rs)/hec
      df[count,4] = j
      df[count,5] = prop[t]
      df[count,6] = sol$capacity
      
      count = count + 1
      
    }
  }
}


colnames(df) = c("sites", "dPC", "area", "n.sites", "prop", "MC")


# simulate restoration scenarios at existing patches
df2 = tibble()

count=1

for(j in 1:length(top.existing$RowID)) {
  
  s = combn(top.existing$RowID, j)
  
  for(i in 1:length(s[1,])){
    
    new.sites = Ndf.extant[which(Ndf.extant$RowID %in% s[,i]),]
    
    for(t in 1:length(prop)){
      
      df2[count,1] = list(list(new.sites[,"Wetland_ID"]$Wetland_ID))
      
      df.temp = Ndf.extant
      rs = new.sites$Suitable_Veg_Area_m2 + (prop[t] * (new.sites$Area_m2-new.sites$Suitable_Veg_Area_m2))
      df.temp[which(df.temp$RowID %in% s[,i]), "Suitable_Veg_Area_m2"] = rs
      
      sol <- meta_capacity(as.dist(d.ex), df.temp$Suitable_Veg_Area_m2/hec, f = dispersal_negexp(alpha), patch_mc = T, ex = ex)
      
      df2[count,2] = (sol$capacity - sol.current$capacity) / sol.current$capacity
      df2[count,3] = sum(rs)/hec
      df2[count,4] = j
      df2[count,5] = prop[t]
      df2[count,6] = sol$capacity
      
      count = count + 1
      
    }
  }
}

colnames(df2) = c("sites", "dPC", "area", "n.sites", "prop", "MC")


# add status column and merge results for existing and potential patch restoration scenarios
df$status = "Restore potential patches"
df2$status = "Increase habitat in existing patches"


##############################################
# extract priority sites for combined strategy

# best scenarios by total effort (total hectares restored)
df$cat.cost = cut(df$area, breaks = seq(0, max(df$area)+1, 1))
tops = df %>% group_by(cat.cost) %>% slice_max(order_by = MC, n = 1, with_ties = FALSE)
tops$bpc = tops$dPC/tops$area
tops$MC / tops$area # benefit per unit area restored

# identify 6 best potential sites
n.pull = round((0.01 * nrow(df)) / length(unique(df$cat.cost)))
top6 = df %>% group_by(cat.cost) %>% slice_max(order_by = dPC/area, n = n.pull, with_ties = FALSE)
top6 = as.data.frame(table(unlist(top6$sites)))
colnames(top6) = c("Wetland_ID", "freq")
top6 = top6[rev(order(top6$freq)),]
top6 = top6[1:6,"Wetland_ID"]


# best scenarios by total effort (total hectares restored)
df2$cat.cost = cut(df2$area, breaks = seq(0, max(df2$area)+1, 1))
tops2 = df2 %>% group_by(cat.cost) %>% slice_max(order_by = MC, n = 1, with_ties = FALSE)
tops2$bpc = tops2$dPC/tops2$area
tops2$MC / tops2$area # benefit per unit area restored

# identify 6 best existing sites
n.pull = round((0.01 * nrow(df2)) / length(unique(df2$cat.cost)))
top26 = df2 %>% group_by(cat.cost) %>% slice_max(order_by = dPC/area, n = n.pull, with_ties = FALSE)
top26 = as.data.frame(table(unlist(top26$sites)))
colnames(top26) = c("Wetland_ID", "freq")
top26 = top26[rev(order(top26$freq)),]
top26 = top26[1:6,"Wetland_ID"]

# combine top 6 potential and top 6 existing patches
c.top = data.frame(Wetland_ID = as.numeric(as.character(c(top6, top26))), status = rep(c("potential", "existing"), each = 6))
combs = Ndf %>% filter(Wetland_ID %in% c.top$Wetland_ID) %>% select(RowID, Wetland_ID) %>%
  left_join(c.top, by = "Wetland_ID") %>%
  select(-Wetland_ID)


###################################################################################
# simulate restoration scenarios at a combination of existing and potential patches

df3 = tibble()

count=1

for(j in 1:length(combs$RowID)) {
  
  s = combn(combs$RowID, j)
  
  for(i in 1:length(s[1,])){

    new.sites = Ndf.potential[which(Ndf.potential$RowID %in% s[,i]),]
    existing.sites = Ndf.extant[which(Ndf.extant$RowID %in% s[,i]),]
    
    
    for(t in 1:length(prop)){
      
      df3[count,1] = list(list(c(new.sites[,"Wetland_ID"]$Wetland_ID, existing.sites[,"Wetland_ID"]$Wetland_ID)))
      
      df.temp = Ndf
      rs = prop[t] * new.sites$Area_m2
      df.temp[which(df.temp$RowID %in% new.sites$RowID), "Suitable_Veg_Area_m2"] = rs
      
      rs2 = existing.sites$Suitable_Veg_Area_m2 + (prop[t] * (existing.sites$Area_m2-existing.sites$Suitable_Veg_Area_m2))
      df.temp[which(df.temp$RowID %in% existing.sites$RowID), "Suitable_Veg_Area_m2"] = rs2
      
      As = df.temp[which(df.temp$Suitable_Veg_Area_m2 >= 100),]
      
      d = Ndm
      d = d[as.character(As$Wetland_ID), as.character(As$Wetland_ID)]
      
      sol <- meta_capacity(as.dist(d), As$Suitable_Veg_Area_m2/hec, f = dispersal_negexp(alpha), patch_mc = T, ex = ex)
      
      df3[count,2] = (sol$capacity - sol.current$capacity) / sol.current$capacity
      df3[count,3] = sum(rs, rs2)/hec
      df3[count,4] = j
      df3[count,5] = prop[t]
      df3[count,6] = sol$capacity
      
      count = count + 1
      
    }
  }
}


colnames(df3) = c("sites", "dPC", "area", "n.sites", "prop", "MC")

# add status column and merge results for existing and potential patch restoration scenarios
df3$status = "Combined strategy"
df3$cat.cost = cut(df3$area, breaks = seq(0, max(df3$area)+1, 1))
dfa = rbind(df, df2, df3)


###########
# plotting
###########

source("pub_theme.R")

# order status for plotting
dfa$status = factor(dfa$status, levels = c("Restore potential patches", "Increase habitat in existing patches", "Combined strategy"))

# % increase in metapopulation capacity as a function of hectares resotred and # sites targeted for restoration
ggplot(dfa, aes(area, 100*dPC, color = n.sites)) +
  geom_point() + theme_Publication() +
  xlab("Restored hectares") + 
  ylab("% Increase in metapopulation capacity") +
  theme(legend.position = "none") +
  facet_wrap(~status) + 
  theme(legend.position = "top") + guides(color = guide_colorbar(title.position = "top")) +
  scale_colour_gradientn(colors = viridis::viridis(11))+
  labs(color = "Number of sites") + 
  geom_hline(yintercept = cap.threshold, linetype = "dashed", linewidth = 1)


# % increase in metapopulation capacity per hectare restored as a function of # sites targeted for restoration 
ggplot(dfa, aes(factor(n.sites), 100*dPC/area, color = n.sites)) +
  geom_sina() + theme_Publication() + 
  xlab("Number of sites restored") + 
  ylab("Benefit per hectare (% increase in MC)") + 
  theme(legend.position = "none") +
  facet_wrap(~status) +
  scale_colour_gradientn(colors = viridis::viridis(11)) 


#####################
# combined strategies

# extract best strategies for each level of effort (total area restored)
dfa$cat.cost = cut(dfa$area, breaks = seq(0, max(dfa$area)+1, 1))
tops.all = dfa %>% group_by(cat.cost, status) %>% slice_max(order_by = dPC, n = 1, with_ties = FALSE)

ggplot(tops.all, aes(area, dPC*100, color = factor(status))) +
  geom_smooth(se = F, linewidth = 2) + theme_Publication() + 
  xlab("Restored hectares") + 
  ylab("% increase in metapopulation capacity") +
  theme(legend.title=element_blank())+
  scale_colour_Publication() + 
  geom_hline(yintercept = cap.threshold, linetype = "dashed", linewidth = 1)

ggplot(tops.all, aes(area, 100*dPC/area, color = factor(status))) +
  geom_smooth(se = F, linewidth = 2) + theme_Publication() + 
  xlab("Restored hectares") + 
  ylab("Benefit per hectare (% increase in MC)") +
  theme(legend.title=element_blank()) +
  scale_colour_Publication()


# filter out scenarios that resulted in a viable metapopulation and extract the top 10% most cost-effective
dfa.all = dfa %>% filter(MC > cap.threshold/100 * sol.current$capacity)
n.pull = round((0.1 * nrow(dfa.all)) / length(unique(dfa.all$cat.cost)))
tops.all = dfa.all %>% group_by(cat.cost) %>% slice_max(order_by = dPC/area, n = n.pull, with_ties = FALSE)

# calculate the frequency that each patch appears in the top 10% of scenarios 
values = as.data.frame(table(unlist(tops.all$sites)))
colnames(values) = c("Wetland_ID", "freq")


# Note that due to the sensitive nature of site locations
# the GPS coordinates are not included in the Zenodo repository
# as such, the following code on lines 342-377 to produce Figure 5 will not run

# read in lat lons
Nll = read_xlsx(paste(site, "_gps.xlsx", sep = ""))

# merge lat lons with main wetland dataframe
dat = merge(Ndf, Nll, by.x = "Wetland_ID", by.y = "Wet_ID")
dat$status = ifelse(dat$Suitable_Veg_Area_m2 > 100, "existing", "potential")

# merge with frequecny of occurance in top 10% of scenarios
dat = merge(dat, values, by = "Wetland_ID", all.x = T)
dat[is.na(dat$freq), "freq"] = 1

# convert frequency to proportion
dat$p = dat$freq / nrow(tops.all)


# priority patches across landscape
ggplot(dat, aes(x = Cent_X, y = Cent_Y, col = status)) + 
  geom_point(aes(size = p)) +
  geom_label_repel(aes(label = Wetland_ID), size = 3) +
  theme_Publication() + #theme(legend.position = "none")+
  scale_colour_Publication()+ 
  theme_void() +
  theme(legend.title = element_blank()) +
  theme(plot.margin = margin(t = 10, r = 10, b = 20, l = 10))


# extract patches included in top scenarios
dat = dat[dat$freq > 1, ]

# size of priority patches
ggplot(dat, aes(x = (Area_m2.x - Suitable_Veg_Area_m2) / hec, y = p, col = status)) + 
  geom_point(aes(size = p)) +
  geom_label_repel(aes(label = Wetland_ID), size = 3) +
  theme_Publication() + 
  theme(legend.position = "none")+
  scale_colour_Publication() +
  xlab("Patch area (hectares)") + ylab("Proportion of priority scenarios")


# summarizing characteristics of best scenarios
counts = data.frame(n = unlist(tops.all$n.sites), x = "priority scenarios")
counts2 = data.frame(n = unlist(dfa.all$n.sites), x = "all scenarios")
props = data.frame(n = unlist(tops.all$prop), x = "priority scenarios")
props2 = data.frame(n = unlist(dfa.all$prop), x = "all scenarios")

xy = rbind(counts, counts2)
xy2 = rbind(props, props2)


# histogram of # patches in top scenarios
plot_multi_histogram <- function(df, feature, label_column, means) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(aes(y = after_stat(ncount)), alpha=0.7, binwidth = 0.5, center = 0, position = 'dodge') +
    geom_vline(xintercept=means, color=brewer.pal(n = 3, name = "Dark2")[2:1], linetype="dashed", linewidth=1) +
    labs(x="Number of sites", y = "Frequency density")
  plt + guides(fill=guide_legend(title=NULL))
}


plot_multi_histogram(xy, 'n', 'x', c(Mode(counts$n)[1], Mode(counts2$n)[1])) + 
  theme_Publication() + 
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.position = "top") +
  scale_x_continuous(limits = c(1,12), breaks = seq(0, 12, 1))


# histogram of % habitat restored in top scenarios
plot_multi_histogram <- function(df, feature, label_column, means) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(aes(y = after_stat(ncount)), alpha=0.7, binwidth = 0.05, center = 0, position = 'dodge') +
    geom_vline(xintercept=means, color=brewer.pal(n = 3, name = "Dark2")[2:1], linetype="dashed", linewidth=1) +
    labs(x="Number of sites", y = "Frequency density")
  plt + guides(fill=guide_legend(title=NULL))
}

plot_multi_histogram(xy2, 'n', 'x', c(mean(props$n)[1], mean(props2$n)[1])) + 
  theme_Publication() + 
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.position = "top") +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  xlab("Proportion of patch restored") + xlim(0,1.05)
