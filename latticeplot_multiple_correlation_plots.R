library(lattice)
library(dplyr)
library(reshape2)
# first we have to convert the MEs dataframe to long form
# add the x value to the Mes dataframe
MEs_cetuximab_trait <- as.data.frame(cbind(cetuximab_trait,MEs))
measure_vars <- names(MEs)
# melt the MEs values and convert to long
MEs_cetuximab_trait <- MEs_cetuximab_trait%>%melt(id.vars="cetuximab_trait",measure.vars=measure_vars)
names(MEs_cetuximab_trait) <- c("cetuximab_trait","MEs","correlation")
# plot the data
xyplot(correlation~cetuximab_trait|MEs,data = MEs_cetuximab_trait,cex=1.5,
       # size of axis scales
       scales=list(cex=c(1.4,1.4)),
       # x label and y label and title
       xlab=list(label="cetuximab trait", fontsize=20),
       ylab=list(label="correlation", fontsize=20),
       main=list(label="correlation of eigengenes with Cetuximab trait",fontsize=16),
       # column for distinguishing colors
       groups = cetuximab_trait,
       # legend
       key=list(space="right",
       points=list(col=c("blue","magenta"), lty=c(3,2), lwd=6),
       text=list(c("non-cetuximab treated"," cetuximab treated"))))

