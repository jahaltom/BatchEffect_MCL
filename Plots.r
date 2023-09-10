
library(ggplot2)



#Read in quant. Modified gene name will be rownames, see README.txt.
rand = read.table("FILE_BestPValue.Random.txt",header=FALSE,sep = '\t',quote="",check.names=FALSE)
exp = read.table("FILE_BestPValue.txt",header=FALSE,sep = '\t',quote="",check.names=FALSE)




png("FILE.png",width=10,height=5,units="in",res=1000)
ggplot(rand, aes(x=V1)) + 
  geom_histogram(color="black", fill="white", binwidth=0.0009) + xlim(0, 0.03)  + ylim(0, 15)  + xlab("Best P-Value") + ylab("Counts") + geom_segment(aes(x = exp$V1,
                   y = 12,
                   xend = exp$V1,
                   yend = 0.1,lineend="butt",linejoin="mitre"),size = 1.2,
               arrow = arrow(length = unit(0.3, "cm"),type = 'closed'),color="red") +
    annotate(
      geom = "text",
      x = exp$V+0.0015, y = 11,
      label = toString(round(exp[1],3)),
      color = "red" ,fontface = "bold",size = 5) +
    annotate(
      geom = "text",
      x = exp$V+0.00275, y = 11.9,
      label = "Experimental",
      color = "red",fontface = "bold",size = 5
    )
  
dev.off() 
