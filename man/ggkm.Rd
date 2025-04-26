\name{ggkm}
\alias{ggkm}
\title{Creates a Kaplan-Meier plot with at risk tables below}
\usage{
ggkm(sfit, table = FALSE, xlabs = "Time-to-event", ylabs = "Survival (\%)", 
     xlims = c(0, max(sfit$time)), ylims = c(0, 1), ystratalabs = names(sfit$strata), 
     ystrataname = "Strata", timeby = signif(max(sfit$time)/7, 1), main = "", 
     pval = FALSE, pvposition = c(0.3, 0.6), marks = TRUE, shape = 3, 
     legend = TRUE, legendposition = c(0.85, 0.8), ci = FALSE, subs = NULL, 
     linecols = "Set1", dashed = FALSE, aspectRatio = 0.7143, black = FALSE, 
     data = NULL, HR = FALSE, incid = FALSE, pvaltxt = NULL, hrtxt = NULL, ...)
}
\arguments{
\item{sfit}{a survfit object}
\item{timeby}{numeric: control the granularity along the time-axis; defaults to 7 time-points. Default = signif(max(sfit$time)/7, 1)}

\item{main}{plot title}

\item{pval}{logical: add the pvalue to the plot?}
\item{marks}{logical: should censoring marks be added?}

\item{subs}{= NULL,}


\item{table}{logical: Create a table graphic below the K-M plot, indicating at-risk numbers?}

\item{xlabs}{x-axis label}

\item{ylabs}{y-axis label}

\item{xlims}{numeric: list of min and max for x-axis. Default = c(0,max(sfit$time))}

\item{ylims}{numeric: list of min and max for y-axis. Default = c(0,1)}

\item{ystratalabs}{character list. A list of names for each strata. Default = names(sfit$strata)}

\item{ystrataname}{The legend name. Default = "Strata"}

\item{shape}{what shape should the censoring marks be, default is a vertical line}

\item{legend}{logical. should a legend be added to the plot?}

\item{legendposition}{numeric. x, y position of the legend if plotted. Default=c(0.85,0.8)}

\item{ci}{logical. Should confidence intervals be plotted. Default = FALSE}

\item{linecols}{Character. Colour brewer pallettes too colour lines. Default ="Set1",}

\item{dashed}{logical. Should a variety of linetypes be used to identify lines. Default = FALSE}
\item{pvposition}{position for the p-value, default = c(0.3, 0.6).}
\item{pvaltxt}{text for the p-value, default is NULL.}
\item{aspectRatio}{add aspect ratio of the plot, default is 0.7134.}
\item{black}{black and white plot, default is FALSE.}
\item{data}{data set for the plot, could be useful when one need HR on the plot.} 
\item{incid}{plot incidence curve instead of KM curve, default is FALSE.}
\item{HR}{add hazards ratio to the plot, default is FALSE.}
\item{hrtxt}{text for the hazards ratio, default is NULL.}
\item{\dots}{additional arguments to be passed to the ggkm function.}
}
\description{
Creates a Kaplan-Meier plot with at risk tables below.

Based on the original packge by Michael Way(https://github.com/michaelway/ggkm), 
I added some minor changes that allow users to:

1) add hazards ratio (HR) and CI for HR to the plot.

2) do the incidence plot instead of KM plot.

3) plot with black and white, sometimes required by a journal.

4) change the aspect ratio of the plot.
}
\examples{
 library(survival)
 data(colon)
 fit <- survfit(Surv(time,status)~rx, data=colon)
# ggkm(fit, timeby=500)
}
\author{
Michael Way(https://github.com/michaelway/ggkm), heavily modified version of a script created by Abhijit Dasgupta with contributions by Gil Tomas.
\url{https://statbandit.wordpress.com/2011/03/08/an-enhanced-kaplan-meier-plot/}
Michael have packaged this function, added functions to namespace and included a range of new parameters. 
Bingshu Chen added more options to the plot (see above) and made minor corrections of the R code and the R docments files to pass the R CMD check.
}
