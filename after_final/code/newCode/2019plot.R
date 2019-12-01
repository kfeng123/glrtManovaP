library(ggplot2)
library(latex2exp)

jiujiuPlot <- function( myFile = 't1l', myLegendPosition = c(0.9,0.4)){
    theOut <- read.csv(paste0(myFile,'.csv'))
    
    tSC <- data.frame(
        SNR=theOut$SNR,
        power=theOut$SC,
        method="Sc"
    )
    tCX <- data.frame(
        SNR=theOut$SNR,
        power=theOut$CX,
        method="CX"
    )
    tHBWW <- data.frame(
        SNR=theOut$SNR,
        power=theOut$HBWW,
        method="HBWW"
        
    )
    tZGZ <- data.frame(
        SNR=theOut$SNR,
        power=theOut$ZGZ,
        method="ZGZ"
    )
    tLFD <- data.frame(
        SNR=theOut$SNR,
        power=theOut$Asy,
        method="LFD"
    )
    df=rbind(tSC,tCX,tHBWW,tZGZ,tLFD)
    
    
    cols <- c("Sc"="red","LFD"="magenta","CX"="orange1","HBWW"="green3","ZGZ"="lightskyblue2")
    lineT <- c("Sc"=2,"LFD"=1,"CX"=3,"HBWW"=4,"ZGZ"=5)
    Labels <- c("Sc"="Sc","LFD"="LFD","CX"="CX","HBWW"="HBWW","ZGZ"="ZGZ")
    g <- guide_legend("Test methods",keywidth=unit(2.25,"cm"),override.aes=list(shape=NA))
    myPlot <- ggplot(df,aes(SNR,power,colour=method,linetype=method))+
        geom_point(size=2)+
        geom_line(size=1.2)+
        scale_linetype_manual(name="Methods",values=lineT,labels=Labels)+
        scale_colour_manual(name="Methods",values=cols,labels=Labels)+
        geom_hline(yintercept = 0.05,linetype=2)+
        ylab("Empirical power")+
        xlab("SNR")+
        scale_y_continuous(breaks=c(seq(0,1,0.2),0.05), minor_breaks=seq(0,1,0.1))+
        theme_bw()+
        guides(colour=g,linetype=g)+
        theme(
            axis.text.y = element_text(size=rel(1.7)),
            axis.text.x = element_text(size=rel(1.7)),
            axis.title.y = element_text(size=rel(1.9)),
            axis.title.x = element_text(size=rel(1.9)),
            legend.position= myLegendPosition,
            legend.margin = margin(3,3,3,3),
            legend.box.background = element_rect(),
            legend.text = element_text(size=rel(1.6)),
            legend.title = element_text(size=rel(1.8)),
            legend.spacing.y = unit(0.5,"line"),
            legend.key.height = unit(1.8,"line"),
            #panel.grid.minor = element_blank(),
            aspect.ratio =0.9
        )
    
    ggsave(paste0(myFile,".eps"),myPlot)
}

jiujiuPlot("t1l",myLegendPosition=c(0.157,0.8))
jiujiuPlot("t1lb",myLegendPosition="none")
jiujiuPlot("t1r",myLegendPosition="none")
jiujiuPlot("t1rb",myLegendPosition="none")


jiujiuPlot("t2l",myLegendPosition = c(0.8,0.4))
jiujiuPlot("t2lb",myLegendPosition="none")
jiujiuPlot("t2r",myLegendPosition="none")
jiujiuPlot("t2rb",myLegendPosition="none")

jiujiuPlot("t3l",myLegendPosition=c(0.157,0.8))
jiujiuPlot("t3lb",myLegendPosition="none")
jiujiuPlot("t3r",myLegendPosition="none")
jiujiuPlot("t3rb",myLegendPosition="none")


jiujiuPlot("t4l",myLegendPosition = c(0.8,0.4))
jiujiuPlot("t4lb",myLegendPosition="none")
jiujiuPlot("t4r",myLegendPosition="none")
jiujiuPlot("t4rb",myLegendPosition="none")