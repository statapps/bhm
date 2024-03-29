#' Creates a Kaplan-Meier plot with at risk tables below
#' @param sfit: a survfit object
#' @param table: logical: Create a table graphic below the K-M plot, indicating at-risk numbers?
#' @param xlabs: x-axis label
#' @param ylabs: y-axis label
#' @param xlims: numeric: list of min and max for x-axis. Default = c(0,max(sfit$time))
#' @param ylims: numeric: list of min and max for y-axis. Default = c(0,1)
#' @param ystratalabs: character list. A list of names for each strata. Default = names(sfit$strata)
#' @param ystrataname: The legend name. Default = "Strata"
#' @param timeby numeric: control the granularity along the time-axis; defaults to 7 time-points. Default = signif(max(sfit$time)/7, 1)
#' @param main plot title
#' @param pval logical: add the pvalue to the plot?
#' @param marks logical: should censoring marks be added?
#' @param shape: what shape should the censoring marks be, default is a vertical line
#' @param legend: logical. should a legend be added to the plot?
#' @param legendposition: numeric. x, y position of the legend if plotted. Default=c(0.85,0.8)
#' @param ci: logical. Should confidence intervals be plotted. Default = FALSE
#' @param subs = NULL,
#' @param linecols: Character. Colour brewer pallettes too colour lines. Default ="Set1",
#' @param dashed: logical. Should a variety of linetypes be used to identify lines. Default = FALSE
#' @author Michael Way, but heavily modified version of a script created by Abhijit Dasgupta with contributions by Gil Tomas.
#' \url{http://statbandit.wordpress.com/2011/03/08/an-enhanced-kaplan-meier-plot/}
#' I have packaged this function, added functions to namespace and included a range of new parameters.
#' @examples
#'  library(survival)
#'  data(colon)
#'  fit <- survfit(Surv(time,status)~rx, data=colon)
#'  ggkm(fit, timeby=500)
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_step
#' @importFrom ggplot2 scale_linetype_manual
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_blank
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 scale_colour_brewer
#' @importFrom ggplot2 geom_ribbon
#' @importFrom grid unit
#' @importFrom gridExtra grid.arrange
#' @importFrom plyr rbind.fill
#' @export

.rbindx <- function(..., dfs=list(...)) {
  ns <- unique(unlist(sapply(dfs, names)))
  do.call(rbind, lapply(dfs, function(x) {
    for(n in ns[! ns %in% names(x)]) {x[[n]] <- NA}; x }))
}

ggkm <- function(sfit,
                 table = FALSE,
                 xlabs = "Time-to-event",
                 ylabs = "Survival (%)",
                 xlims = c(0,max(sfit$time)),
                 ylims = c(0,1),
                 ystratalabs = names(sfit$strata),
                 ystrataname = "Strata",
                 timeby = signif(max(sfit$time)/7, 1),
                 main = "",
                 pval = FALSE,
                 pvposition=c(0.30,0.60),
                 marks = TRUE,
                 shape = 3,
                 legend = TRUE,
                 legendposition=c(0.85,0.8),
                 ci = FALSE,
                 subs = NULL,
                 linecols="Set1",
                 dashed= FALSE,
		 aspectRatio = 0.7143, 
		 black = FALSE,
		 data = NULL,
		 HR = FALSE, 
		 incid = FALSE, 
		 pvaltxt = NULL, 
		 hrtxt = NULL, 
                 ...) {
  

  #################################
  # sorting the use of subsetting #
  #################################
  
  times <- seq(0, max(sfit$time), by = timeby)
  
  if(is.null(subs)){
    if(length(levels(summary(sfit)$strata)) == 0) {
      subs1 <- 1
      subs2 <- 1:length(summary(sfit,censored=T)$time)
      subs3 <- 1:length(summary(sfit,times = times,extend = TRUE)$time)
    } else {
      subs1 <- 1:length(levels(summary(sfit)$strata))
      subs2 <- 1:length(summary(sfit,censored=T)$strata)
      subs3 <- 1:length(summary(sfit,times = times,extend = TRUE)$strata)
    }
  } else{
    for(i in 1:length(subs)){
      if(i==1){
        ssvar <- paste("(?=.*\\b=",subs[i],sep="")
      }
      if(i==length(subs)){
        ssvar <- paste(ssvar,"\\b)(?=.*\\b=",subs[i],"\\b)",sep="")
      }
      if(!i %in% c(1, length(subs))){
        ssvar <- paste(ssvar,"\\b)(?=.*\\b=",subs[i],sep="")
      }
      if(i==1 & i==length(subs)){
        ssvar <- paste("(?=.*\\b=",subs[i],"\\b)",sep="")
      }
    }
    subs1 <- which(regexpr(ssvar,levels(summary(sfit)$strata), perl=T)!=-1)
    subs2 <- which(regexpr(ssvar,summary(sfit,censored=T)$strata, perl=T)!=-1)
    subs3 <- which(regexpr(ssvar,summary(sfit,times = times,extend = TRUE)$strata, perl=T)!=-1)
  }
  
  if(!is.null(subs)) pval <- FALSE
  
  ##################################
  # data manipulation pre-plotting #
  ##################################
  
  
  
  if(length(levels(summary(sfit)$strata)) == 0) {
    #[subs1]
    if(is.null(ystratalabs)) ystratalabs <- as.character(sub("group=*","","All"))
  } else {
    #[subs1]
    if(is.null(ystratalabs)) ystratalabs <- as.character(sub("group=*","",names(sfit$strata)))
  }
  
  if(is.null(ystrataname)) ystrataname <- "Strata"
  m <- max(nchar(ystratalabs))
  times <- seq(0, max(sfit$time), by = timeby)
  
  if(length(levels(summary(sfit)$strata)) == 0) {
    Factor <- factor(rep("All",length(subs2)))
  } else {
    Factor <- factor(summary(sfit, censored = T)$strata[subs2])
  }
  
  #Data to be used in the survival plot
  surv = sfit$surv[subs2]
  if (ylims[2] == 100) {
    surv = surv*100
    s0 = 100
  } else s0 = 1

  if(incid == TRUE) surv = 1 - surv

  df <- data.frame(
    time = sfit$time[subs2],
    n.risk = sfit$n.risk[subs2],
    n.event = sfit$n.event[subs2],
    n.censor = sfit$n.censor[subs2],
    #surv = sfit$surv[subs2],
    surv = surv, 
    strata = Factor,
    upper = sfit$upper[subs2],
    lower = sfit$lower[subs2]
  )
  
  #Final changes to data for survival plot
  levels(df$strata) <- ystratalabs
  zeros <- data.frame(time = 0, surv = s0,
                      strata = factor(ystratalabs, levels=levels(df$strata)),
                      upper = 1, lower = 1)

  #print(zeros)
  # df <- rbind.fill(zeros, df)
  df = .rbindx(zeros, df)
  #print(head(df))
  
  d  <- length(levels(df$strata))
  
  ###################################
  # specifying axis parameteres etc #
  ###################################
  
  if(dashed == TRUE){
    linetype=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678")
  } else {
    linetype=c("solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid")
  }
  if(black == TRUE)linetype=c("solid", "twodash", "dashed", "dotted", "dotdash", "solid", "solid", "solid", "solid", "solid", "solid")

  p <- ggplot(df, aes(x=time, y=surv, colour=strata, linetype=strata)) +
    ggtitle(main)
  
  #Set up theme elements
  p <- p + theme_bw() +
    theme(axis.title.x = element_text(vjust = 0.7),
	  aspect.ratio=aspectRatio,
          panel.grid.minor = element_blank(),
          axis.line = element_line(size =0.5, colour = "black"),
          legend.position = legendposition,
          legend.background = element_rect(fill = NULL),
          legend.key = element_rect(colour = NA),
          panel.border = element_blank(),
          plot.margin = unit(c(0, 1, .5,ifelse(m < 10, 1.5, 2.5)),"lines"),
          panel.grid.major = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")) +
    scale_x_continuous(xlabs, breaks = times, limits = xlims) +
    scale_y_continuous(ylabs, limits = ylims)
  
  
  #Add 95% CI to plot
  if(ci == TRUE)
    p <- p +  geom_ribbon(data=df, aes(ymin = df$lower, ymax = df$upper), fill = "grey", alpha=0.25, colour=NA)
  
  #Removes the legend:
  if(legend == FALSE)
    p <- p + theme(legend.position="none")
  
  #Add lines to plot
  if(black == TRUE) {
    p = p + geom_step(size = 0.75, colour = "black") + 
      scale_linetype_manual(name = ystrataname, values=linetype) + 
      scale_colour_brewer(name = ystrataname)
  } else {
    p = p + geom_step(size = 0.7) + 
      scale_linetype_manual(name = ystrataname, values=linetype) +
      scale_colour_brewer(name = ystrataname, palette=linecols)
  }
  
  #Add censoring marks to the line:
  if(marks == TRUE)
    p <- p + geom_point(data = subset(df, df$n.censor >= 1), aes(x = time, y = surv), shape = shape, colour = "black")
  
  ## Create a blank plot for place-holding
  blank.pic <- ggplot(df, aes(time, surv)) +
    ####### remove visible lines for the placeholder;
    geom_blank() + theme_bw(base_line_size = 1/222) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.title.x = element_blank(),axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),panel.border = element_blank())
  
  #####################
  # p-value placement #
  #####################
  
  if(length(levels(summary(sfit)$strata)) == 0) pval <- FALSE
  pvposition[2] = pvposition[2]*ylims[2]
  #print(pvposition)
  #print(as.integer(max(sfit$time)*(pvposition[1])))

  if(pval == TRUE) {
    if(is.null(data))stop("Need data from the input for the p-value.")
    #sdiff <- survdiff(eval(sfit$call$formula), data = eval(sfit$call$data))
    sdiff <- survdiff(eval(sfit$call$formula), data = data)
    pvalue <- pchisq(sdiff$chisq,length(sdiff$n) - 1,lower.tail = FALSE)
    if(is.null(pvaltxt)) pvaltxt = ifelse(pvalue < 0.0001,"Log-Rank test, p < 0.0001",
            paste("Log-rank test, p = ", signif(pvalue, 2), sep=""))
    # MOVE P-VALUE LEGEND HERE BELOW [set x and y]
    p <- p + annotate("text",x = as.integer(max(sfit$time)*(pvposition[1])), 
                             y = pvposition[2], label = pvaltxt)
  }

  ##########################
  # Hazard Ratio placement #
  ##########################

  if(length(levels(summary(sfit)$strata)) == 0) HR = FALSE
 
  if(HR == TRUE) {
    if(is.null(data))stop("Need data from the input for the Hazard Ratio.")
    phfit  <- coxph(eval(sfit$call$formula), data = data)
    hrd = signif(summary(phfit)$conf.int, 3)
    if(is.null(hrtxt)) hrtxt = paste("HR = ", hrd[1], " (95%CI, ", hrd[3], "-", hrd[4], ")", sep="")
    # MOVE HAZARD RATIO LEGEND HERE BELOW [set x and y]
    p <- p + annotate("text",x = as.integer(max(sfit$time)*(pvposition[1])), 
             y = pvposition[2]-0.1*(ylims[2]-ylims[1]), label = hrtxt)
  }

  
  ###################################################
  # Create table graphic to include at-risk numbers #
  ###################################################
  
  if(length(levels(summary(sfit)$strata)) == 0) {
    Factor <- factor(rep("All",length(subs3)))
  } else {
    Factor <- factor(summary(sfit,times = times,extend = TRUE)$strata[subs3])
  }
  
  if(table == TRUE) {
    n.risk = NULL
    risk.data <- data.frame(
      strata = Factor,
      time = summary(sfit,times = times,extend = TRUE)$time[subs3],
      n.risk = summary(sfit,times = times,extend = TRUE)$n.risk[subs3]
    )
    
    risk.data$strata <- factor(risk.data$strata, levels=rev(levels(risk.data$strata)))
    
    data.table <- ggplot(risk.data,aes(x = time, y = strata, label = format(n.risk, nsmall = 0))) +
      geom_text(size = 3.5) + theme_bw() +
      scale_y_discrete(breaks = as.character(levels(risk.data$strata)),
                       labels = rev(ystratalabs)) +
      scale_x_continuous("Numbers at risk", limits = xlims) +
      theme(axis.title.x = element_text(size = 10, vjust = 1),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),axis.text.x = element_blank(),
            axis.ticks = element_blank(),axis.text.y = element_text(face = "bold",hjust = 1))
    
    data.table <- data.table +
      theme(legend.position = "none") + xlab(NULL) + ylab(NULL)
    
    # ADJUST POSITION OF TABLE FOR AT RISK
    data.table <- data.table +
      theme(plot.margin = unit(c(-1.5, 1, 0.1, ifelse(m < 10, 2.5, 3.5) - 0.15 * m), "lines"))
  }
  
  
  #######################
  # Plotting the graphs #
  #######################
  
  if(table == TRUE){
    grid.arrange(p, blank.pic, data.table, clip = FALSE, nrow = 3,
                 ncol = 1, heights = unit(c(2, .1, .25), c("null", "null", "null")))
  }else {
    p
  }
}
