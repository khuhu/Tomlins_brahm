ascat_cn <- function(df, chromTextSpec = NULL, prefix = NULL){
  
  if(is.null(chromTextSpec)){
    chromTextdf <- read.table("/home/kevhu/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "mm10"){
    chromTextdf <- read.table("/home/kevhu/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "hg19") {
    chromTextdf <- read.table("/br_z1/kevin_storage/misc/20210801hg19_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  }
  
  chromBreak <- c(0, chromTextdf$chromBreaksPos)
  for (i in unique(df$sample)) {
    if (is.null(prefix)) {
      main <- i
    } else{
      main  <- paste("gam", prefix, i, sep = "_")
    }
    
    print(paste0("main ", main))
    tmpDf <- df[which(df$sample == i),]
    tmpDf_major <- tmpDf[,1:5]
    tmpDf_major$col <- "red"
    tmpDf_minor <- tmpDf[,c(1:4, 6)]
    tmpDf_minor$col <- "blue"
    colnames(tmpDf_major)[5] <- "cn"
    colnames(tmpDf_minor)[5] <- "cn"
    tmpDf_major$cn <- ifelse(tmpDf_major$cn == 0, tmpDf_major$cn, tmpDf_major$cn + 0.1)
    tmpDf_minor$cn <- ifelse(tmpDf_minor$cn == 0, tmpDf_minor$cn + 0.1, tmpDf_minor$cn - 0.4)
    tmpDf2 <- rbind(tmpDf_major, tmpDf_minor)
    tmpDf2 <- tmpDf2[which(tmpDf2$chr %in% c(1:19)),]
    tmpDf2 <- tmpDf2[order(as.numeric(tmpDf2$chr)),]
    colnames(tmpDf2) <- c("ID", "chrom", "loc.start", "loc.end", "seg.mean", "col")
    
    tmpDf2$loc.start <- tmpDf2$loc.start/1e6
    tmpDf2$loc.end <- tmpDf2$loc.end/1e6
    
    for (j in unique(tmpDf2$chrom)) {
      tmpDf2$loc.start[which(tmpDf2$chrom == j)] <- tmpDf2$loc.start[which(tmpDf2$chrom == j)] +
        chromTextdf$graphingStart[which(chromTextdf$chrom == j)]
      tmpDf2$loc.end[which(tmpDf2$chrom == j)] <- tmpDf2$loc.end[which(tmpDf2$chrom == j)] +
        chromTextdf$graphingStart[which(chromTextdf$chrom == j)]
    }
    
    tmpDf2 <- tmpDf2[,c("chrom", "loc.start", "loc.end", "seg.mean", "col")]
    colnames(tmpDf2) <- c("chrom", "start", "end","cn", "col")
    
    cnplot <- ggplot(tmpDf2) + geom_hline(yintercept = c(0, 1, 2, 3,4, 5), color = "#000000") + 
      geom_rect(aes(xmin = start, xmax = end, ymin = cn, ymax = cn + 0.3), fill = tmpDf2$col) + 
      geom_vline(xintercept=chromBreak) + theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 5, by = 1),
                                                             limits = c(0, 5.2)) +
      geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 4.3, label = chrom), size = 2.5) + ggtitle(main) +
      theme(plot.title = element_text(hjust = 0.5))
    
    cnplot_grob <- ggplotGrob(cnplot)
    
    dev.off()
    # pdf(file = paste0("/br_z1/kevin_storage/misc/20211101noTcAmpAnalysis/20211122_filtComps_20_", i, ".pdf"), width = 15, height = 7)
    png(filename = paste0("combinedGraphs_", main, ".png"), width = 1800, height = 500)
    grid::grid.newpage()
    grid::grid.draw(rbind(cnplot_grob, tmpBafGrob))
    dev.off()
  }
}

ascat_cn2 <- function(df, chromTextSpec = NULL, prefix = NULL){
  
  if(is.null(chromTextSpec)){
    chromTextdf <- read.table("/home/kevhu/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "mm10"){
    chromTextdf <- read.table("/home/kevhu/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "hg19") {
    chromTextdf <- read.table("/br_z1/kevin_storage/misc/20210801hg19_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  }
  
  chromBreak <- c(0, chromTextdf$chromBreaksPos)
  for (i in unique(df$sample)) {
    if (is.null(prefix)) {
      main <- i
    } else{
      main  <- paste("gam", prefix, i, sep = "_")
    }
    
    print(paste0("main ", main))
    tmpDf <- df[which(df$sample == i),]
    tmpDf$cn <- tmpDf$nMajor + tmpDf$nMinor
    tmpDf$col <- "#013220"
    tmpDf <- tmpDf[which(tmpDf$chr %in% c(1:19)),]
    tmpDf <- tmpDf[, c("sample", "chr", "startpos", "endpos", "cn", "col")]
    colnames(tmpDf) <- c("ID", "chrom", "loc.start", "loc.end", "seg.mean", "col")
    tmpDf2 <- tmpDf[order(as.numeric(tmpDf$chr)),]
    
    
    tmpDf2$loc.start <- tmpDf2$loc.start/1e6
    tmpDf2$loc.end <- tmpDf2$loc.end/1e6
    
    for (j in unique(tmpDf2$chrom)) {
      tmpDf2$loc.start[which(tmpDf2$chrom == j)] <- tmpDf2$loc.start[which(tmpDf2$chrom == j)] +
        chromTextdf$graphingStart[which(chromTextdf$chrom == j)]
      tmpDf2$loc.end[which(tmpDf2$chrom == j)] <- tmpDf2$loc.end[which(tmpDf2$chrom == j)] +
        chromTextdf$graphingStart[which(chromTextdf$chrom == j)]
    }
    
    tmpDf2 <- tmpDf2[,c("chrom", "loc.start", "loc.end", "seg.mean", "col")]
    colnames(tmpDf2) <- c("chrom", "start", "end","cn", "col")
    
    cnplot <- ggplot(tmpDf2) + geom_hline(yintercept = c(0, 1, 2, 3,4, 5,6,7), color = "#000000") + 
      geom_rect(aes(xmin = start, xmax = end, ymin = cn, ymax = cn + 0.3), fill = tmpDf2$col) + 
      geom_vline(xintercept=chromBreak) + theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 7, by = 1),
                                                             limits = c(0, 7.2)) +
      geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 6.3, label = chrom), size = 2.5) + ggtitle(main) +
      theme(plot.title = element_text(hjust = 0.5))
    
    cnplot_grob <- ggplotGrob(cnplot)
    
    dev.off()
    # pdf(file = paste0("/br_z1/kevin_storage/misc/20211101noTcAmpAnalysis/20211122_filtComps_20_", i, ".pdf"), width = 15, height = 7)
    png(filename = paste0("combinedGraphs_", main, ".png"), width = 1800, height = 500)
    grid::grid.newpage()
    grid::grid.draw(rbind(cnplot_grob, tmpBafGrob))
    dev.off()
  }
}



freqPlot_baf <- function(df, main = "no title", chromTextSpec = NULL){
  require(ggplot2)
  
  if(is.null(chromTextSpec)){
    chromTextdf <- read.table("/home/kevhu/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "mm10"){
    chromTextdf <- read.table("/home/kevhu/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "hg19") {
    chromTextdf <- read.table("/br_z1/kevin_storage/misc/20210801hg19_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  }
  
  
  
  chromBreak <- c(0, chromTextdf$chromBreaksPos)
  
  ### divide positions by megabase locations
  
  df$CHROM <- as.numeric(str_remove(df$CHROM, "chr"))
  df$point <- df$POS/1e6
  
  for (i in unique(df$CHROM)) {
    df$point[which(df$CHROM == i)] <- df$point[which(df$CHROM == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  
  df <- df[,c("CHROM", "point", "AF")]
  colnames(df) <- c("chrom", "pos", "af")
  
  
  if (nrow(df) == 0) {
    ggplot(df, aes(x = pos, y = af)) + geom_hline(yintercept = c(1,0.75, 0.5, 0.25, 0), color = "#D4D4D4")+
      geom_point(size = 0.25, color = df$color) + geom_vline(xintercept=chromBreak) + theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 1, by = 0.10),
                                                             limits = c(0,1.2)) +
      geom_text(data = chromTextdf, aes(x = xpos, y = ypos + .2, label = chrom), size = 2.5) 
  } else{
    df$color <- "#000000"
    df$color[which(df$af < 0.3)] <- "#FF0000"
    df$color[which(df$af > 0.7)] <- "#FF0000"
    
    ggplot(df, aes(x = pos, y = af)) + geom_hline(yintercept = c(1,0.75, 0.5, 0.25, 0), color = "#D4D4D4")+
      geom_point(size = 0.25, color = df$color) + geom_vline(xintercept=chromBreak) + theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 1, by = 0.10),
                                                             limits = c(0,1.2)) +
      geom_text(data = chromTextdf, aes(x = xpos, y = ypos + .2, label = chrom), size = 2.5) 
  }
}