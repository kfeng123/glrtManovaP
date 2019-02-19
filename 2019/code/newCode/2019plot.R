library(ggplot2)
library(latex2exp)

jiujiuPlot <- function( myFile = 't1l'){
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
    g <- guide_legend("Test methods",keywidth=unit(1,"cm"),override.aes=list(shape=NA))
    library(ggplot2)
    library(latex2exp)
    myPlot <- ggplot(df,aes(SNR,power,colour=method,linetype=method))+
        geom_point()+
        geom_line()+
        scale_linetype_manual(name="Methods",values=lineT,labels=Labels)+
        scale_colour_manual(name="Methods",values=cols,labels=Labels)+
        ylab("Empirical power")+
        xlab("SNR")+
        #scale_x_continuous(breaks=seq(0,0.9,0.1))+
        theme_bw()+
        guides(colour=g,linetype=g)+
        theme(
            axis.title.y = element_text(size=rel(1.2)),
            axis.title.x = element_text(size=rel(1.2)),
            legend.position= c(0.9,0.4),
            legend.margin = margin(3,3,3,3),
            legend.box.background = element_rect(),
            aspect.ratio =0.8
        )
    
    ggsave(paste0(myFile,".eps"),myPlot)
}

jiujiuPlot("t1l")
jiujiuPlot("t1lb")
jiujiuPlot("t1r")
jiujiuPlot("t1rb")


jiujiuPlot("t2l")
jiujiuPlot("t2lb")
jiujiuPlot("t2r")
jiujiuPlot("t2rb")

jiujiuPlot("t3l")
jiujiuPlot("t3lb")
jiujiuPlot("t3r")
jiujiuPlot("t3rb")


jiujiuPlot("t4l")
jiujiuPlot("t4lb")
jiujiuPlot("t4r")
jiujiuPlot("t4rb")