# required libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readxl)
library(metacapa)

# specify site (Site1 or Site2)
site = "Site1" 

# set model parameters
alpha = 2      # dispersal measure from Martin et al. 2025
top = 12       # number of top sites to look at restoring
hec = 10000     # scaling for m2 area into hectares
ex = 0.1       # scaling extinction with patch area
prop = seq(0.1, 1, 0.1) # proportion of area restored in each site


# read in wetland attributes
Ndf = read_xlsx(paste(site, "_wetland_attributes.xlsx", sep = ""))

# read in distance matrix
Ndm = read_xlsx(paste(site, "_Distance_Matrix.xlsx", sep = ""))
Ndm = as.data.frame(Ndm)
rownames(Ndm) = as.character(Ndm[,1])
Ndm = Ndm[,-1]
Ndm = Ndm/1000 # scale distances to km


# split sites with and without habitat
Ndf.extant = Ndf[which(Ndf$Suitable_Veg_Area_m2 >= 100),]
Ndf.potential = Ndf[which(Ndf$Suitable_Veg_Area_m2 < 100),]

# current distance matrix
d.ex = Ndm[as.character(Ndf.extant$Wetland_ID), as.character(Ndf.extant$Wetland_ID)]


###########
# potential landscape capacity, current landscape capacity, new sites capacity 
sol.potential1 <- meta_capacity(as.dist(d.ex), Ndf.extant$Area_m2/hec, f = dispersal_negexp(alpha), patch_mc = T, ex = ex)
sol.potential2 <- meta_capacity(as.dist(Ndm), Ndf$Area_m2/hec, f = dispersal_negexp(alpha), patch_mc = T, ex = ex)
sol.current <- meta_capacity(as.dist(d.ex), Ndf.extant$Suitable_Veg_Area_m2/hec, f = dispersal_negexp(alpha), patch_mc = T, ex = ex)

# extract top existing sites for habitat improvement
tops1 = Ndf[order(sol.potential1$patch_mc), "RowID"]$RowID
top.existing = Ndf[tops1[tops1 %in% Ndf.extant$RowID][1:top],"RowID"]

# extract top potential sites for restoration
tops2 = Ndf[order(sol.potential2$patch_mc), "RowID"]$RowID
top.potential = Ndf[tops2[tops2 %in% Ndf.potential$RowID][1:top],"RowID"]


# restoration scenarios at potential sites
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
df$status = "Restore potential patches"


# removal experiment to estimate optimal restoration sites (partial restoration based on fixed area)
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
df2$status = "Increase habitat in existing patches"


###


# merge scenarios for restoration of potential and existing sites
dfa = rbind(df, df2)

# calculate improvement per hectare restored
dfa$cat.cost = cut(dfa$area, breaks = seq(0, max(dfa$area), 1))

# extract optimal scenarios for each cost category and strategy
tops.all = dfa %>% group_by(cat.cost, status) %>% top_n(1, dPC)


###########
# plotting
###########

theme_Publication <- function(base_size=12, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"),
            panel.border = element_rect(color = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90, margin = margin(r = 10)),
            axis.title.x = element_text(margin = margin(t = 10)),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="gray",fill="gray"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}


gg1 = ggplot(dfa, aes(area, 100*dPC, color = factor(n.sites))) +
  geom_point() + theme_Publication() + 
  xlab("Restored hectares") + 
  ylab("Improvement in metapopulation capacity") +
  theme(legend.position = "none") +
  facet_wrap(~status)

gg2 = ggplot(dfa, aes(factor(n.sites), 100*dPC/area, color = factor(n.sites))) +
  geom_sina() + theme_Publication() + 
  xlab("Number of sites restored") + 
  ylab("Benefit per hectare (% increase in MC)") + 
  theme(legend.position = "none") +
  facet_wrap(~status)

gg3 = ggplot(tops.all, aes(area, dPC*100, color = factor(status))) +
  geom_smooth(se = F, linewidth = 2) + theme_Publication() + 
  geom_smooth(data = tops.combine, se = F, linewidth = 2) +
  xlab("Restored hectares") + 
  ylab("% Improvement in metapopulation capacity") +
  theme(legend.title=element_blank())+
  scale_colour_Publication()

gg4 = ggplot(tops.all, aes(area, 100*dPC/area, color = factor(status))) +
  geom_smooth(se = F, linewidth = 2) + theme_Publication() + 
  geom_smooth(data = tops.combine, se = F, linewidth = 2) +
  xlab("Restored hectares") + 
  ylab("Benefit per hectare (% increase in MC)") +
  theme(legend.title=element_blank()) +
  scale_colour_Publication()