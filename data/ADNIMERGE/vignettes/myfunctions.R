longsums <- function(Data, oc = "ADAS13")
{
  cat("\n\\clearpage\n")
  cat(paste("### ", label(Data[, oc]), "change from baseline"))
  cat("\n\n")
  dd1 <- merge(subset(dd, VISCODE == "bl", c("RID", "DX")), 
    Data[,c("RID", "VISCODE", oc)], by = "RID", all.y = TRUE)
  dd1$Month <- ifelse(dd1$VISCODE %in% c("sc", "bl"), 0, as.numeric(gsub("m", "", dd1$VISCODE)))
  dd1$Y <- dd1[, oc]
  dd1 <- merge(dd1, subset(dd1, Month == 0, c("RID", "Y")), by = "RID", 
    all.x = TRUE, suffixes = c("", ".0"))
  dd1$Y.change <- with(dd1, Y - Y.0)
  dd1 <- filter(dd1, !is.na(DX))
  dd1 <- filter(dd1, !is.na(Y.change))
  dd1 <- dd1 %>%
    group_by(Month, DX) %>%
    filter(n() > 4) %>%
    ungroup()
  
  s <- with(dd1, Hmisc::summarize(Y.change,
    by = llist(Month, DX),
    stat.name = 'Change',
    function(x) c(smean.cl.normal(x), N = sum(!is.na(x)))))
  s <- s[order(s$DX, s$Month), ]
  s <- subset(s, !is.na(DX) & N > 4)
  m <- length(unique(s$Month))
  n <- length(unique(s$DX))
  Palette <- brewer.pal(n, "Set1")
  blank.pic <- ggplot(s, aes(Month, Change)) +
      geom_blank() + theme_bw() +
      theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
           axis.title.x = element_blank(),axis.title.y = element_blank(),
           axis.ticks = element_blank(),
           panel.grid.major = element_blank(),panel.border = element_blank())
  p <- ggplot(s, aes(Month, Change, group = DX, color = DX)) +
    geom_line() +
    geom_smooth(aes(ymax = Upper, ymin = Lower), 
      size = 1.5, stat = 'identity') +
    ylab(paste(label(Data[, oc]), "change")) +
    scale_x_continuous(breaks = unique(dd1$Month)) +
    scale_colour_manual(values=Palette) +
    ggtitle("Mean and 95% CI") + theme(legend.position = "none")
  
  data.table <- ggplot(s, aes(x = Month, y = DX,
                              label = format(N, nsmall = 0), color = DX)) +
    geom_text(size = 3) + theme_bw() +
    scale_colour_manual(values=Palette) +
    scale_y_discrete(breaks = levels(s$DX),
      labels = levels(s$DX)) +
    theme(axis.title.x = element_text(vjust = 1),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_blank(),axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(hjust = 1, face = "bold", color = Palette)) +
    theme(legend.position = "none") + xlab("") + ylab(NULL) +
    theme(plot.margin = unit(c(-1.5, 1, 0.1, ifelse(m < 10, 1.5, 2.5) - 0.28 * m), "lines")) 
    # ADJUST POSITION OF TABLE FOR AT RISK

  print(p + data.table + plot_layout(ncol = 1))

  s <- with(dd1, Hmisc::summarize(Y.change,
    by = llist(Month, DX),
    stat.name = 'N',
    function(x) summary(x)[c('N', 'Min.', 'Q.1', 'Median', 'Mean', 'Q.3', 'Max.', 'SD')]))
  # browser()
  N_tab <- dd1 %>%
    group_by(Month, DX) %>%
    dplyr::summarise(N = sum(!is.na(get(oc))), 
                     num_NA= sum(is.na(get(oc)))) 
    
  s$N <- N_tab$N
  s$NA. <- N_tab$num_NA
  tab <- s[order(s$DX, s$Month), ]
  tab <- subset(tab, !is.na(DX) & Month != 0)
  rcmd <- rep("", nrow(tab))
  rcmd[seq(1, length(rcmd), by = 2)] <- "rowcolor[gray]{.9}"
  tab[, 2] <- as.character(tab[, 2])
  tab[, -c(2)] %>%
    kable(format = knitr::opts_knit$get("rmarkdown.pandoc.to"), longtable = TRUE, row.names = FALSE) %>%
    pack_rows(index = table(tab[,2]))
  
}

make_first_to_final_trans <- function(RID, stage, Year){
  dd <- data.frame(RID = RID, stage = stage, Year = Year)
  ddf <- c()
  for(rid in unique(dd$RID)){
    subdd <- subset(dd, RID == rid)
    nr <- nrow(subdd)
    ddf <- rbind(ddf, c(id = rid, start = subdd$Year[1], stop = subdd$Year[nr], 
      start.stage = subdd$stage[1], 
      end.stage = subdd$stage[nr]))  
  }
  data.frame(ddf)
}

make_all_intermediate_trans <- function(RID, stage, Year){
  dd <- data.frame(RID = RID, stage = stage, Year = Year)
  dd2 <- c()
  # id  start   stop start.stage end.stage
  for(rid in unique(dd$RID)){
    subdd <- subset(dd, RID == rid)
    if(nrow(subdd) > 1){
      start <- subdd$Year[1]
      start.stage <- subdd$stage[1]
      for(j in 1:(nrow(subdd)-1)){
        if(subdd$stage[j] != start.stage){
          dd2 <- rbind(dd2, c(id = rid, start = start, stop = subdd$Year[j], 
            start.stage = start.stage, end.stage = subdd$stage[j]))
          start.stage <- subdd$stage[j]
        }
        start <- subdd$Year[j]
      }
      dd2 <- rbind(dd2, c(id = rid, start = start, stop = subdd$Year[j+1], 
        start.stage = start.stage, 
        end.stage = ifelse(start.stage == subdd$stage[j+1], 0, subdd$stage[j+1])))  
    }
  }
  data.frame(dd2)
}

obj_memory <- function(x){
  EDU <- x[1]
  LM <- x[2]
  if(is.na(EDU) | is.na(LM)) return(NA)
  if(EDU >= 16){
    if(LM <= 8) return("LMCI")
    if(LM <= 11) return("EMCI")
    if(LM > 11) return("NL")    
  }
  if(EDU <= 15 & EDU >= 8){
    if(LM <= 4) return("LMCI")
    if(LM <= 9) return("EMCI")
    if(LM > 9) return("NL")    
  }
  if(EDU <= 7){
    if(LM <= 2) return("LMCI")
    if(LM <= 6) return("EMCI")
    if(LM > 6) return("NL")    
  }
}

####################################################################
##               My Plot Method                                   ##
####################################################################


my_msSurv_plot <- function (x, states="ALL", trans="ALL", plot.type="stateocc",
                    CI=TRUE, ci.level=0.95, ci.trans="linear", 
                    state.names=NULL, timelab="Event Times", 
                    legend.title = "Stage", ...) {


              require(ggplot2)
              plot.type <- match.arg(plot.type, c("stateocc", "transprob","entry.sub",
                                                  "entry.norm","exit.sub","exit.norm"))


              ####################################################################
              ## State occupation probabilities plot
              ####################################################################

              if (plot.type=="stateocc") {

                  CIs <- MSM.CIs(x, ci.level=0.95) ## Calling CIs
                  if (states[1]=="ALL") states <- nodes(tree(x))

                  if(is.null(state.names)){
                    labels <- unique(sort(f.st))
                  }else{
                    labels <- state.names
                  }
                  
                  f.st <- factor(states, labels = labels)
                  
                  ls <- length(states)
                  sl <- which(nodes(tree(x))%in%as.numeric(states)) ## location of states in the matrix

                  if (CI==TRUE & !is.null(var.sop(x))) {
                      rd <- CIs$CI.p
                      dimnames(rd)$state <- gsub("p", "State", dimnames(rd)$state)
                      y <- as.vector(rd[,1,sl])
                      y2 <- as.vector(rd[,2,sl]) ## lower limit
                      y3 <- as.vector(rd[,3,sl]) ## upper limit
                      xvals <- rep(et(x), length(states))
                      f.st <- factor(rep(dimnames(rd)$state[sl], each=dim(rd)[1]), labels = labels)
                      # st.plot <- xyplot(y + y2 + y3 ~ xvals | f.st, allow.multiple=TRUE,
                      #                   type="s", lty=c(1,2,2), col=c(1,2,2), ...)
                      # st.plot <- update(st.plot, main="Plot of State Occupation Probabilites",
                      #                   xlab="Event Times", ylab="State Occupation Probabilities",
                      #                   key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)),
                      #                   text=list(c("Est", "Lower CI", "Upper CI")),
                      #                   columns=3))
                      st.plot <- ggplot(data.frame(xvals, y, y2, y3, f.st),
                        aes(xvals, y, color = f.st))+
                        geom_line() +
                        geom_smooth(aes(ymax = y3, ymin = y2), stat = 'identity')
                        
                      st.plot <- st.plot +
                      ggtitle("Plot of State Occupation Probabilites") +
                        xlab(timelab) + ylab("State Occupation Probabilities") +
                        labs(color = legend.title)
                      print(st.plot)
                  } else {
                      if (CI==TRUE)
                          cat("Warning: 'var.sop'  is NULL and therefore CIs not plotted. \n")
                      rd <- CIs$CI.p
                      dimnames(rd)$state <- gsub("p", "State", dimnames(rd)$state)
                      y <- as.vector(rd[,1,sl])
                      xvals <- rep(et(x), length(states))
                      f.st <- as.factor(rep(dimnames(rd)$state[sl], each=dim(rd)[1]))
                      st.plot <- xyplot(y ~ xvals | f.st, type="s",col=1, ...)
                      st.plot <- update(st.plot, main="Plot of State Occupation Probabilites",
                                        xlab="Event Times", ylab="State Occupation Probabilities",
                                        key = list(lines=list(col=c(1), lty=c(1)), text=list(c("Est"))))
                      print(st.plot)
                  } ## end of no CIs

              } ## end of state occ plot

              ####################################################################
              ## Transition probabilities plot
              ####################################################################

              if (plot.type=="transprob") {

                  CIs <- MSM.CIs(x, ci.level, ci.trans, sop=FALSE) ## Calling CIs
                  all.trans <- pos.trans(x)
                  if (trans[1] =="ALL") trans <- all.trans

                  rd <- CIs$CI.trans
                  dimnames(rd)$trans <- gsub(" ", " to ", dimnames(rd)$trans)
                  if(!is.null(state.names)){
                    for(st in length(state.names):1){
                      dimnames(rd)$trans <- gsub(st, state.names[st], dimnames(rd)$trans)
                    }
                  }
                  tr <- which(all.trans%in%trans) ## location of states in the matrix

                  if (CI==TRUE & !is.null(cov.AJs(x))) {

                      y <- as.vector(rd[,1,tr])
                      y2 <- as.vector(rd[,2,tr]) ## lower limit
                      y3 <- as.vector(rd[,3,tr]) ## upper limit
                      xvals <- rep(et(x), length(tr))
                      f.tp <- as.factor(rep(dimnames(rd)$trans[tr], each=dim(rd)[1]))
                      # tr.plot <- xyplot(y + y2 + y3 ~ xvals | f.tp, allow.multiple=TRUE,
                      #                   type="s", lty=c(1,2,2), col=c(1,2,2),...)
                      # tr.plot <- update(tr.plot, main="Plot of Transition Probabilites",
                      #                   xlab="Event Times", ylab="Transition Probabilites",
                      #                   key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)),
                      #                   text=list(c("Est", "Lower CI", "Upper CI")),
                      #                   columns=3))
                      tr.plot <- ggplot(data.frame(xvals, y, y2, y3, f.tp),
                        aes(xvals, y, color = f.tp))+
                        geom_line() +
                        geom_smooth(aes(ymax = y3, ymin = y2), stat = 'identity')
                        
                      tr.plot <- tr.plot +
                      ggtitle("Plot of Transition Probabilites") +
                        xlab(timelab) + ylab("Transition Probabilites") +
                        labs(color = legend.title)
                      print(tr.plot)

                  } else {
                      if (CI==TRUE)
                          cat("Warning: 'cov.AJs'  is NULL and therefore CIs not plotted. \n")
                      y <- as.vector(rd[,1,tr])
                      xvals <- rep(et(x), length(tr))
                      f.tp <- as.factor(rep(dimnames(rd)$trans[tr], each=dim(rd)[1]))
                      tr.plot <- xyplot(y ~ xvals | f.tp, type="s", lty=1, col=1, ...)
                      tr.plot <- update(tr.plot, main="Plot of Transition Probabilites",
                                        xlab="Event Times", ylab="Transition Probabilities",
                                        key = list(lines=list(col=1, lty=1), text=list("Est"), columns=1))
                      print(tr.plot)

                  }

              } ## end of 'transprob' plot


              ####################################################################
              ## Entry subdistribution function
              ####################################################################

              if (plot.type=="entry.sub") {

                  enter <- names(which(!(sapply(inEdges(tree(x)), function(x) length(x) == 0))))
                  if (states[1]=="ALL") states <- enter

                  f.st <- factor(states)
                  ls <- length(states)
                  sl <- which(nodes(tree(x))%in%states) ## location of states in the matrix


                  if (CI==TRUE & !is.null(Fsub.var(x))) {

                      CIs <- Dist.CIs(x, ci.level, ci.trans, norm=FALSE) ## Calling CIs for subdistribution
                      rd <- CIs$CI.Fs
                      dimnames(rd)$state=gsub("F", "State", dimnames(rd)$state)
                      y <- as.vector(rd[,1,sl])
                      y2 <- as.vector(rd[,2,sl]) ## lower limit
                      y3 <- as.vector(rd[,3,sl]) ## upper limit
                      xvals <- rep(et(x), length(states))
                      f.st <- as.factor(rep(dimnames(rd)$state[sl], each=dim(rd)[1]))
                      ent.plot <- xyplot(y + y2 + y3 ~ xvals | f.st, allow.multiple=TRUE,
                                         type="s", lty=c(1,2,2), col=c(1,2,2), ...)
                      ent.plot <- update(ent.plot, main="Plot of State Entry Time Subdistributions",
                                         xlab="Event Times", ylab="State Entry Time Subdistributions",
                                         key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)),
                                         text=list(c("Est", "Lower CI", "Upper CI")),
                                         columns=3))
                      print(ent.plot)
                  }  else {
                          if (CI==TRUE)
                              cat("Warning: 'Fsub.var'  is NULL and therefore CIs not plotted. \n")
                          rd <- Fsub(x)
                          dimnames(rd)[[2]]=gsub("F", "State", dimnames(rd)[[2]])
                          y <- as.vector(rd[,sl])
                          xvals <- rep(et(x), length(states))
                          f.st <- as.factor(rep(dimnames(rd)[[2]][sl], each=dim(rd)[1]))
                          ent.plot <- xyplot(y ~ xvals | f.st, type="s", col=1, ...)
                          ent.plot <- update(ent.plot, main="Plot of State Entry Time Subdistributions",
                                             xlab="Event Times", ylab="State Entry Time Subdistributions",
                                             key = list(lines=list(col=c(1), lty=c(1)), text=list(c("Est"))))
                          print(ent.plot)
                      }
              } ## end of entry subdistribution plot


              ####################################################################
              ## State entry distribution (normalized)
              ####################################################################

              if (plot.type=="entry.norm") {

                  enter <- names(which(!(sapply(inEdges(tree(x)), function(x) length(x) == 0))))
                  if (states[1]=="ALL") states <- enter

                  f.st <- factor(states)
                  ls <- length(states)
                  sl <- which(nodes(tree(x))%in%states) ## location of states in the matrix


                  if (CI==TRUE & !is.null(Fnorm.var(x))) {

                      CIs <- Dist.CIs(x,ci.level,ci.trans,norm=TRUE) ## Calling CIs for normalized distribution
                      rd <- CIs$CI.Fs
                      dimnames(rd)$state=gsub("F", "State", dimnames(rd)$state)
                      y <- as.vector(rd[,1,sl])
                      y2 <- as.vector(rd[,2,sl]) ## lower limit
                      y3 <- as.vector(rd[,3,sl]) ## upper limit
                      xvals <- rep(et(x), length(states))
                      f.st <- as.factor(rep(dimnames(rd)$state[sl], each=dim(rd)[1]))
                      ent.plot <- xyplot(y + y2 + y3 ~ xvals | f.st, allow.multiple=TRUE,
                                         type="s",lty=c(1,2,2),col=c(1,2,2), ...)
                      ent.plot <- update(ent.plot, main="Plot of Normalized State Entry Time Distributions",
                                         xlab="Event Times", ylab="Normalized State Entry Time Distributions",
                                         key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)),
                                         text=list(c("Est", "Lower CI", "Upper CI")),
                                         columns=3))
                      print(ent.plot)
                      } else {
                          if (CI==TRUE)
                              cat("Warning: 'Fnorm.var' is NULL and therefore CIs not plotted. \n")
                          rd <- Fnorm(x)
                          dimnames(rd)[[2]]=gsub("F", "State", dimnames(rd)[[2]])
                          y <- as.vector(rd[,sl])
                          xvals <- rep(et(x), length(states))
                          f.st <- as.factor(rep(dimnames(rd)[[2]][sl], each=dim(rd)[1]))
                          ent.plot <- xyplot(y ~ xvals | f.st, type="s", col=1, ...)
                          ent.plot <- update(ent.plot, main="Plot of Normalized State Entry Time Distributions",
                                             xlab="Event Times", ylab="Normalized State Entry Time Distributions",
                                             key = list(lines=list(col=c(1), lty=c(1)), text=list(c("Est"))))
                          print(ent.plot)
                      }

              } ## end of normalized entry distribution plot

              ####################################################################
              ## State exit subdistribution
              ####################################################################

              if (plot.type=="exit.sub") {

                  transient <- names(which(sapply(edges(tree(x)), function(x) length(x) > 0)))
                  if (states[1]=="ALL") states <- transient

                  f.st <- factor(states)
                  ls <- length(states)
                  sl <- which(nodes(tree(x))%in%as.numeric(states)) ## location of states in the matrix


                  if (CI==TRUE & !is.null(Gsub.var(x))) {

                      CIs <- Dist.CIs(x,ci.level,ci.trans,norm=FALSE)
                      ## Calling CIs, adding norm arg for subdistribution CIs
                      rd <- CIs$CI.Gs
                      dimnames(rd)$state=gsub("G", "State", dimnames(rd)$state)
                      y <- as.vector(rd[,1,sl])
                      y2 <- as.vector(rd[,2,sl]) ## lower limit
                      y3 <- as.vector(rd[,3,sl]) ## upper limit
                      xvals <- rep(et(x), length(states))
                      f.st <- as.factor(rep(dimnames(rd)$state[sl], each=dim(rd)[1]))
                      exit.plot <- xyplot(y + y2 + y3 ~ xvals | f.st, allow.multiple=TRUE,
                                          type="s", lty=c(1,2,2), col=c(1,2,2),...)
                      exit.plot <- update(exit.plot, main="Plot of State Exit Time Distributions",
                                          xlab="Event Times", ylab="State Exit Time Distributions",
                                          key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)),
                                          text=list(c("Est", "Lower CI", "Upper CI")),
                                          columns=3))
                      print(exit.plot)
                      }  else {
                          if (CI==TRUE)
                              cat("Warning: 'Gsub.var' is NULL and therefore CIs not plotted. \n")
                          rd <- Gsub(x)
                          dimnames(rd)[[2]]=gsub("G", "State", dimnames(rd)[[2]])
                          y <- as.vector(rd[,sl])
                          xvals <- rep(et(x), length(states))
                          f.st <- as.factor(rep(dimnames(rd)[[2]][sl], each=dim(rd)[1]))
                          exit.plot <- xyplot(y ~ xvals | f.st, type="s", col=1)
                          exit.plot <- update(exit.plot, main="Plot of State Exit Time Distributions",
                                              xlab="Event Times",ylab="State Exit Time Distributions",
                                              key = list(lines=list(col=c(1), lty=c(1)), text=list(c("Est"))))
                          print(exit.plot)
                      }
              } ## end of exit.sub


              ####################################################################
              ## State exit distribution (normalized)
              ####################################################################

              if (plot.type=="exit.norm") {

                  transient <- names(which(sapply(edges(tree(x)), function(x) length(x) > 0)))
                  if (states[1]=="ALL") states<-transient

                  f.st <- factor(states)
                  ls <- length(states)
                  sl <- which(nodes(tree(x))%in%as.numeric(states)) ## location of states in the matrix


                  if (CI==TRUE & !is.null(Gnorm.var(x))) {

                      CIs <- Dist.CIs(x,ci.level,ci.trans,norm=TRUE) ## Calling CIs
                      rd <- CIs$CI.Gs
                      dimnames(rd)$state=gsub("G", "State", dimnames(rd)$state)
                      y <- as.vector(rd[,1,sl])
                      y2 <- as.vector(rd[,2,sl]) ## lower limit
                      y3 <- as.vector(rd[,3,sl]) ## upper limit
                      xvals <- rep(et(x), length(states))
                      f.st <- as.factor(rep(dimnames(rd)$state[sl], each=dim(rd)[1]))
                      exit.plot <- xyplot(y + y2 + y3 ~ xvals | f.st, allow.multiple=TRUE,
                                          type="s",lty=c(1,2,2), col=c(1,2,2), ...)
                      exit.plot <- update(exit.plot, main="Plot of State Exit Time Distributions",
                                          xlab="Event Times", ylab="State Exit Time Distributions",
                                          key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)),
                                          text=list(c("Est", "Lower CI", "Upper CI")),
                                          columns=3))
                      print(exit.plot)
                      } else {
                          if (CI==TRUE)
                              cat("Warning: 'Gnorm.var'  is NULL and therefore CIs not plotted. \n")
                          rd <- Gnorm(x)
                          dimnames(rd)[[2]]=gsub("G", "State", dimnames(rd)[[2]])
                          y <- as.vector(rd[,sl])
                          xvals <- rep(et(x), length(states))
                          f.st <- as.factor(rep(dimnames(rd)[[2]][sl], each=dim(rd)[1]))
                          exit.plot <- xyplot(y ~ xvals | f.st, type="s",col=1, ...)
                          exit.plot <- update(exit.plot, main="Plot of State Exit Time Distributions",
                                              xlab="Event Times", ylab="State Exit Time Distributions",
                                              key = list(lines=list(col=c(1), lty=c(1)), text=list(c("Est"))))
                          print(exit.plot)
                      }
              } ## end of exit.norm

          }

############################################################
##                 Adding Start Times                     ##
############################################################

## Adds starting time to Data
Add.start <- function(Data) {

    Data$start <- 0
    idx <- which(table(Data$id)>1)

    for (i in names(idx)) {
        ab <- Data[which(Data$id==i), ]
        ab <- with(ab, ab[order(ab$stop), ])
        ab2 <- which(Data$id==i) ## row numbers in Data
        start2 <- vector(length=length(ab2))
        start2[1] <- 0
        start2[2:length(ab2)] <- ab$stop[1:length(ab2)-1]
        Data$start[ab2] <- start2
    } ## end of for loop

    new.data <- data.frame(id=Data$id, start=Data$start, stop=Data$stop,
                           start.stage=Data$start.stage, end.stage=Data$end.stage)
    res <- new.data
}


####################################################################
##   Adding 'Dummy' States for Censoring and Left Truncation
####################################################################

Add.States <- function(tree, LT) {

    ## Adding censoring state to Nodes & Edges
    Nodes <- c(nodes(tree), "0")
    Edges <- edges(tree)
    Edges[["0"]] <- character(0)

    nt.states <- names(which(sapply(Edges, function(x) length(x)>0))) ## nontermewinal states

    for (stage in nt.states) {
        Edges[[stage]] <- c("0", Edges[[stage]])
    }

    ## tree for censored data
    tree0 <- new("graphNEL", nodes=Nodes, edgeL=Edges, edgemode="directed")

    ## Adding "Left Truncated" State
    if (LT) {
        Nodes <- c(nodes(tree0), "LT")
        Edges[["LT"]] <- nt.states
        nt.states.LT <- names(which(sapply(Edges, function(x) length(x)>0))) ## nonterminal states
        treeLT <- new("graphNEL", nodes=Nodes, edgeL=Edges, edgemode="directed")
    }

    if (LT) {
        return(list(tree0=tree0, nt.states=nt.states, nt.states.LT=nt.states.LT, treeLT=treeLT))
    } else {
        return(list(tree0=tree0, nt.states=nt.states))
    }

}

############################################################
##          Adding Dummy "LT" obs to Data set             ##
############################################################

LT.Data <- function(Data) {

    Data <- Data[order(Data$id), ]  ## make sure id's line up below
    ids <- unique(Data$id)
    stop.time <- with(Data, tapply(start, id, min))
    enter.st <- by(Data, Data$id, function(x) x$start.stage[which.min(x$start)])
    ## dummy initial stage
    dummy <- data.frame(id = ids, start = -1, stop = stop.time, start.stage="LT",
                        end.stage=as.character(enter.st))

    Data <- rbind(Data, dummy)
    Data <- with(Data, Data[order(id, stop), ])

    return(Data=Data)

}



############################################################
##              Counting Process & At Risk                ##
############################################################

## Assuming stage 0 is censored and stage 1 is the initial state
## state 'LT' for left truncated data

CP <- function(tree, tree0, Data, nt.states) {
    ## tree0=tree for uncens, tree0=tree0 for cens, tree0=treeLT for LT
    ## NOTE - nt.states includes 'LT' for LT data

    ## unique stop times in entire data set
    ## Exclude stop times of zero (for LT data)
    times <- sort(unique(Data$stop[!Data$stop==0]))

    ## names for dNs transitions
    lng <- sapply(edges(tree0)[nodes(tree0)%in%nt.states], length)
    ds <- paste("dN", rep(nodes(tree0)[nodes(tree0)%in%nt.states], lng),
                unlist(edges(tree0)[nodes(tree0)%in%nt.states]))

    ## names for at-risk calcs
    nt.states2 <- nt.states[!nt.states=="LT"]
    ys <- paste("y", nt.states2)  ## used for Ys

    ##  matrix of # of transitions, initialize to zeros
    dNs <- matrix(0, nrow=length(times), ncol=length(ds))

    ##  matrix of total # of transitions from a state, initialize to zeros
    sum_dNs <- matrix(0, nrow=length(times), ncol=length(nt.states))

    ##  matrix of at-risk sets for each stage at each time
    Ys <- matrix(NA, nrow=length(times), ncol=length(ys))

    ## names of rows/columns for vectors/matrices
    rownames(dNs) <- rownames(sum_dNs) <- rownames(Ys) <- times
    colnames(dNs) <- ds
    colnames(Ys) <- ys
    colnames(sum_dNs) <- paste("dN", nt.states, ".")

    ## Calculations for transitions 'dNs'
    ## Outer loop = all nodes w/transitions into them
    nodes.into <- nodes(tree0)[sapply(inEdges(tree0), function(x) length(x) > 0)]
    for (i in nodes.into) {
        ## Inner loop = all nodes which transition into node i
        nodes.from <- inEdges(tree0)[[i]]
        for (j in nodes.from) {
            nam2 <- paste("dN", j, i)
            idx <- which(Data$end.stage==i & Data$start.stage==j)
            tmp.tab <- table(Data$stop[idx][!Data$stop[idx]==0])  ##  no. trans at each trans time
            dNs[names(tmp.tab), nam2] <- tmp.tab
        }
    }


    ## start counts at time == 0 for below
    start.stages <- Data$start.stage[Data$start==0]
    start.stages <- factor(start.stages, levels=nt.states2, labels=nt.states2)
    start.cnts <- table(start.stages)

    ## Calculations for at-risk sets 'Ys'
    ## only need for non-terminal nodes - use 'nt.states2' to exclude 'LT' state
    for (i in nt.states2) {

        n <- start.cnts[i]
        nam <-paste("y", i)

        if (length(inEdges(tree0)[[i]]) > 0)
            into.node <- paste("dN", inEdges(tree0)[[i]], i) else into.node <- NULL
        if (length(edges(tree0)[[i]])>0)
            from.node <- paste("dN", i, edges(tree0)[[i]]) else from.node <- NULL

        Ys[, nam] <- c(n, n + cumsum(rowSums(dNs[, into.node, drop=FALSE])) -
                      cumsum(rowSums(dNs[, from.node, drop=FALSE])))[-(nrow(Ys)+1)]

    } ## end of loop for Ys

    ## Counting transitions from different states (ie: state sums)
    sum_dNs <- matrix(nrow=nrow(dNs), ncol=length(nt.states))
    rownames(sum_dNs) <- rownames(dNs) ##
    colnames(sum_dNs) <- paste("dN", nt.states, ".")
    a <- strsplit(colnames(sum_dNs), " ")
    a2 <- strsplit(colnames(dNs), " ")
    uni <- unique(sapply(a, function(x) x[2]))##  gives the unique states exiting

    for (i in uni) { ## calculating the dNi.s
        b <- which(sapply(a, function(x) x[2]==i))
        b2 <- which(sapply(a2, function(x) x[2]==i))
        sum_dNs[, b] <- rowSums(dNs[, b2, drop=FALSE])
    } ## end of for loop for calculating dNi.s

    list(dNs=dNs, Ys=Ys, sum_dNs=sum_dNs)

} ## end of function


############################################################
##            Datta-Satten Estimation                     ##
############################################################

DS <- function(nt.states, dNs, sum_dNs, Ys, Cens="0", cens.type) {
    ## Calculating dNs, sum_dNs, and Y from D-S(2001) paper
    ## Dividing dNs*, sum_dNs*, & Y* by K to get dNs^, sum_dNs^, & Ys^
    ## Make sure nt.states is from the non-LT

    res <- strsplit(colnames(dNs), " ") ## string splits names
    res2 <- strsplit(colnames(Ys), " ")  ## string split names of Ys
    res3 <- strsplit(colnames(sum_dNs), " ") ## string splits names of dNs

    ##  looks at censored columns, needed for D-S est
    DS.col.idx <- which(sapply(res, function(x) x[3]==Cens))
    ##  looks at censored columns, needed for D-S est
    DS2.col.idx <- which(sapply(res2, function(x) x[2]%in%nt.states))
    ##  looks at censored columns, needed for D-S est
    DS3.col.idx <- which(sapply(res3, function(x) x[2]%in%nt.states))

    ## for INDEPENDENT censoring
    if (cens.type=="ind") {

	K <- vector(length=nrow(dNs))
	dN0 <- rowSums(dNs[, DS.col.idx, drop=FALSE])
	Y0 <- rowSums(Ys[, DS2.col.idx, drop=FALSE]) ## those at risk of being censored
	N.Y <- ifelse(dN0/Y0=="NaN", 0, dN0/Y0)
	colnames(N.Y) <- NULL
	H.t <- cumsum(N.Y) ## calculating the hazard
	k <- exp(-H.t)
	K <- c(1, k[-length(k)])

	dNs.K <- dNs/K  ## D-S dNs
	Ys.K <- Ys/K  ## D-S Ys
        sum_dNs.K <- sum_dNs/K
    } ## end of independent censoring


    ## for DEPENDENT censoring
    if (cens.type=="dep") {

        dN0 <- dNs[, DS.col.idx]
        Y0 <- Ys[, DS2.col.idx] ## those at risk of being censored

        N.Y <- ifelse(dN0/Y0=="NaN", 0, dN0/Y0)
	colnames(N.Y) <- paste(colnames(dN0), "/", colnames(Y0))

        H.t <- apply(N.Y, 2, function(x) cumsum(x))
        K <- exp(-H.t)
        ## K <- apply(k, 2, function(x) c(1, x[-length(x)]))  ## maybe don't need

	dNs.K <- dNs; Ys.K <- Ys; sum_dNs.K <- sum_dNs
	for (i in nt.states) {
            K.idx <- which(sapply(strsplit(colnames(N.Y), " "), function(x) x[2]==i))
            dN.idx <- which(sapply(res, function(x) x[2]==i))
            sum_dNs.idx <- which(sapply(res3, function(x) x[2]==i))
            Ys.idx <- which(sapply(res2, function(x) x[2]==i))
            dNs.K[, dN.idx] <- dNs[, dN.idx]/K[, K.idx]
            sum_dNs.K[, sum_dNs.idx] <- sum_dNs[, sum_dNs.idx]/K[, K.idx]
            Ys.K[, Ys.idx] <- Ys[, Ys.idx]/K[, K.idx]
	}
    } ## end of dependent censoring

    res <- list(dNs.K=dNs.K, Ys.K=Ys.K, sum_dNs.K=sum_dNs.K)
    return(res)

} ## end of D-S function

############################################################
##           Reducing dNs & Ys to event times             ##
############################################################

Red <- function(tree, dNs, Ys, sum_dNs, dNs.K, Ys.K, sum_dNs.K) {

    ## tree is original tree currently inputted by user
    ## dNs, sum.dNS, & Ys come from CP function
    ## K comes from DS

    ## reducing dNs & Ys to just event times & noncens/nontruncated states
    res <- strsplit(colnames(dNs), " ") ## string splits names
    ##  looks at noncensored columns
    col.idx <- which(sapply(res, function(x) x[2]%in%nodes(tree) & x[3]%in%nodes(tree)))
    ## identifies times where transitions occur
    row.idx <- which(apply(dNs[, col.idx, drop=FALSE], 1, function(x) any(x>0)))
    dNs.et <- dNs[row.idx, col.idx, drop=FALSE] ## reduces dNs

    res2 <- strsplit(colnames(Ys), " ") ## string split names of Ys
    nt.states.f <- names(which(sapply(edges(tree), function(x) length(x)>0))) ## nonterminal states
    col2.idx <- which(sapply(res2, function(x) x[2]%in%nt.states.f)) ## ids nonterminal columns
    Ys.et <- Ys[row.idx, col2.idx, drop=FALSE] ## reduces Ys

    col3.idx <- which(sapply(strsplit(colnames(sum_dNs), " "), function(x) x[2]%in%nodes(tree)))
    sum_dNs.et <- sum_dNs[row.idx, col3.idx, drop=FALSE]

    dNs.K.et <- dNs.K[row.idx, col.idx, drop=FALSE]
    Ys.K.et <- Ys.K[row.idx, col2.idx, drop=FALSE]
    sum_dNs.K.et <- sum_dNs.K[row.idx, col3.idx, drop=FALSE]

    ans <- list(dNs=dNs.et, Ys=Ys.et, sum_dNs=sum_dNs.et, dNs.K=dNs.K.et, Ys.K=Ys.K.et, sum_dNs.K=sum_dNs.K.et)
    return(ans)

}


############################################################
##     AJ Estimates and State Occupation Probabilities    ##
############################################################

AJ.estimator <- function(ns, tree, dNs.et, Ys.et, start.probs) {

    ## currently ns is defined in main function
    ## tree needs to be uncensored tree
    cum.tm <- diag(ns)
    colnames(cum.tm) <- rownames(cum.tm) <- nodes(tree)

    ps <- matrix(NA, nrow=nrow(dNs.et), ncol=length(nodes(tree)))
    rownames(ps) <- rownames(dNs.et); colnames(ps) <- paste("p", nodes(tree))
    all.dA <- all.I.dA <- all.AJs <- array(dim=c(ns, ns, nrow(dNs.et)),
                                           dimnames=list(rows=nodes(tree),
                                           cols=nodes(tree), time=rownames(dNs.et)))

    for (i in 1:nrow(dNs.et)) { ##loop through times

        I.dA <- diag(ns) ## creates trans matrix for current time
        dA <- matrix(0, nrow=ns, ncol=ns)
        colnames(I.dA) <- rownames(I.dA) <- colnames(dA) <- rownames(dA) <- nodes(tree)

        idx <- which(dNs.et[i, , drop=FALSE]>0)  ## transition time i
        t.nam <- colnames(dNs.et)[idx] ## gets names of transitions (ie:  dN##)
        tmp <- strsplit(t.nam, " ")    ## splits title of dN##
        start <- sapply(tmp, function(x) x[2])
        end <- sapply(tmp, function(x) x[3])  ## pulls start & stop states as character strings
        idxs <- matrix(as.character(c(start, end)), ncol=2)
        idxs2 <- matrix(as.character(c(start, start)), ncol=2)

        dA[idxs] <- dNs.et[i, idx]/Ys.et[i, paste("y", start)]
        if (length(idx)==1) {
            dA[start, start] <- -dNs.et[i, idx]/Ys.et[i, paste("y", start)]
        } else {
            dA[idxs2] <- -rowSums(dA[start, ])
        }

        I.dA <- I.dA + dA ## I+dA (transition) matrix

        all.dA[, , i] <- dA     ## stores all dA matrices
        all.I.dA[, , i] <- I.dA ## array for storing all tran matrices

        cum.tm <- cum.tm %*% I.dA  ## Multiply current.tm and cum.tm to get matrix for current time
        all.AJs[, , i] <- cum.tm   ## A-J estimates, stored in array

        ## multiply by start.probs to allow for starting states other than '1'
        ps[i, ] <- start.probs%*%all.AJs[, , i] ## state occupation probabilities

    } ## end of loop

    list(ps=ps, AJs=all.AJs, I.dA=all.I.dA)
} ## end of function


############################################################
##           State Entry/Exit Distributions               ##
############################################################

Dist <- function(ps, ns, tree) {
    ## ps from AJ.estimator function
    ## tree needs to be uncensored tree

    initial <- which(sapply(inEdges(tree), function(x) !length(x)>0)) ## initial states, no Fs
    terminal <- which(sapply(edges(tree), function(x) !length(x)>0)) ## terminal states, no Gs

    Fnorm <- Fsub <- matrix(0, nrow=nrow(ps), ncol=ns) ## entry distn
    rownames(Fnorm) <- rownames(Fsub) <- rownames(ps)
    colnames(Fnorm) <- colnames(Fsub) <- paste("F", nodes(tree))

    Gsub <- Gnorm <- matrix(0, nrow=nrow(ps), ncol=ns) ## exit distn
    rownames(Gnorm) <- rownames(Gsub) <- rownames(ps)
    colnames(Gnorm) <- colnames(Gsub) <- paste("G", nodes(tree))

    ## looping through nodes
    for (i in 1:ns) {
        node <- nodes(tree)[i]
        later.stages <- names(acc(tree, node)[[1]])
        stages <- c(node, later.stages)

        Fsub[, i] <- f.numer <- rowSums(ps[, paste("p", stages), drop=FALSE])
        Fnorm[, i] <- f.numer/f.numer[length(f.numer)]

        if (length(stages)==1) next

        Gsub[, i] <- g.numer <- rowSums(ps[, paste("p", later.stages), drop=FALSE])
        Gnorm[, i] <- g.numer/f.numer[length(f.numer)]


    } ## end of for loop

    Fr <- strsplit(colnames(Fnorm), " ")
    Fs.idx <- which(sapply(Fr, function(x) x[2]%in%names(initial)))
    Fnorm[, Fs.idx] <- Fsub[, Fs.idx] <- NA

    Gr <- strsplit(colnames(Gnorm), " ")
    Gs.idx <- which(sapply(Gr, function(x) x[2]%in%names(terminal)))
    Gnorm[, Gs.idx]<- Gsub[, Gs.idx] <- NA

    list(Fnorm=Fnorm, Gsub=Gsub, Fsub=Fsub, Gnorm=Gnorm)
} ## end of function

############################################################
###         Variance of AJ estimates                    ####
############################################################

var.fn <- function(tree, ns, nt.states, dNs.et, Ys.et, sum_dNs, AJs, I.dA, ps) {

    ## covariance matrices
    state.names <- nodes(tree)
    trans.names <- paste(rep(state.names, ns), rep(state.names, each=ns))
    cov.dA <- cov.AJs <- array(0, dim = c(ns^2, ns^2, nrow(dNs.et)))
    dimnames(cov.dA) <- dimnames(cov.AJs) <- list(trans.names, trans.names, rownames(dNs.et))
    res.array <- array(0, dim=c(ns^2, ns^2))
    colnames(res.array) <- rownames(res.array) <- trans.names


    ## Ident matrix for Kronecker products
    bl.Id <- diag(1, (ns)^2)
    Id <- diag(1, ns)

    ## matrix for var matrix of state occup prob
    var.sop <- matrix(0, nrow=nrow(dNs.et), ncol=ns)
    colnames(var.sop) <- paste("Var", "p", state.names)
    rownames(var.sop) <- rownames(ps)

    ## NOTE - see note below for use of v.p ...
    v.p <- matrix(0, ns, ns)

    for (i in 1:nrow(dNs.et)) { ## loop through times


        ## VARIANCE OF A-J ( TRANS PROB MATRIX P(0, t) )
        ## Equation 4.4.20 in Anderson 1993

        ## loop on the blocks (g) - only needed for non-terminal states
	for (outer in nt.states) {

            ## find positioning of outer in state.names
            idx.outer <- which(state.names==outer)
            ## tmp matrix for calculating variance of transitions in each block (g)
            tm <- matrix(0, nrow=ns, ncol=ns)
            colnames(tm) <- rownames(tm) <- state.names

            ## loop in the blocks
            for (j in 1:ns) {

                ## this just fills in upper diagonal matrix
                ## use symmetry to fill in rest
                for (k in j:ns) {

                    statej <- state.names[j]
                    statek <- state.names[k]
                    outer.Ys <- paste("y", outer)
                    outer.sum_dNs <- paste("dN", outer, ".")

                    if (Ys.et[i, outer.Ys]==0) {  ## if Y_g = 0 then covariance = 0
        	  	tm[j, k] <- 0
	        	next
                    }

                    if (statej == outer & statek == outer) {  ## 3rd formula
			tm[j, k] <- (Ys.et[i, outer.Ys] - sum_dNs[i, outer.sum_dNs])*sum_dNs[i, outer.sum_dNs] / Ys.et[i, outer.Ys]^3

                    }  else if (statej == outer & statek != outer) {  ## 2nd formula
                        name <- paste("dN", outer, statek)
			if (!name%in%colnames(dNs.et)) next
			tm[j, k] <- -(Ys.et[i, outer.Ys] - sum_dNs[i, outer.sum_dNs])*dNs.et[i, name] / Ys.et[i, outer.Ys]^3
                    } else if (statej != outer & statek == outer) {  ## 2nd formula pt 2, for recurrent
                        name <- paste("dN", outer, statej)
			if (!name%in%colnames(dNs.et)) next
			tm[j, k] <- -(Ys.et[i, outer.Ys] - sum_dNs[i, outer.sum_dNs])*dNs.et[i, name]/Ys.et[i, outer.Ys]^3
                    } else { ## 1st formula
			namek <- paste("dN", outer, statek)
			namej <- paste("dN", outer, statej)
			if (!(namej%in%colnames(dNs.et) & namek%in%colnames(dNs.et))) next
			tm[j, k] <- (ifelse(j==k, 1, 0)*Ys.et[i, outer.Ys] - dNs.et[i, namej])*dNs.et[i, namek]/Ys.et[i, outer.Ys]^3
                    } ## end of if/else statements
                } ## end of k loop
            } ## end of j loop

            tm[lower.tri(tm)] <- t(tm)[lower.tri(tm)]

            res.array[(seq(1, ns*(ns-1)+1, by=ns) + idx.outer - 1), (seq(1, ns*(ns-1)+1, by=ns) + idx.outer - 1)] <- tm

	}## end of outer loop

        ## array holding var-cov matrix for I+dA matrix at each time (differential of NA estim)
	cov.dA[, , i] <- res.array

	if (i==1) {
            cov.AJs[, , i] <- bl.Id%*% cov.dA[, , i] %*% bl.Id
        } else {
            cov.AJs[, , i] <- (t(I.dA[, , i]) %x% Id) %*% cov.AJs[, , i-1] %*%((I.dA[, , i]) %x% Id) +
            (Id %x% AJs[, , i-1]) %*% cov.dA[, , i]  %*% (Id%x% t(AJs[, , i-1]))
        }
        ## note:  AJs[, , i-1] corresponds to P(0, t-)
        ## cov.AJs is the var/cov est of P(0, t)

        ## calculating the variance of state occupation prob
	for (j in state.names) { ## loop through states

            name <- paste(state.names[1], j)
            part1 <- var.pkj0t <- cov.AJs[name, name, i]
            res <- strsplit(colnames(ps), " ")
            col.idx <- which(sapply(res, function(x) x[2]== j)) ## looks at state transitioned to
            b.t <- AJs[, col.idx, i] ## creating vector of col for current state from trans prob

            ## NOTE - This is ZERO when all indiv start from single state ...
            part2 <- t(b.t)%*%v.p%*%b.t ## should be 0 when P(0, t)
            ## right now forced to be 0 by the way v.p defined outside of time loop
            res.varp <- part1 + part2      ## calculating var/cov matrix for time i, state j
            var.sop[i, col.idx] <- res.varp ## storing var/cov calc for time i,n state j

	} ## closes states loop
    } ## end of time loop

    list(cov.AJs=cov.AJs, cov.dA=cov.dA, var.sop=var.sop)

}## end of function


#########################################################
##                BS Variance                          ##
#########################################################

## Needed For:
## 1. Dependent Censoring
## 2. Entry / Exit functions
## 3. SOPs when > 1 starting state

BS.var <- function(Data, tree, ns, et, cens.type, B, LT) {

    n <- length(unique(Data$id)) ## sample size
    ids <- unique(Data$id)

    ## storage for bootstrap estimates of transition probability matrices
    bs.est <- array(dim=c(length(nodes(tree)), length(nodes(tree)), length(et), B),
                    dimnames=list(rows=nodes(tree), cols=nodes(tree), time=et))
    bs.ps <- array(dim=c(length(et), ns, B))
    rownames(bs.ps) <- et
    colnames(bs.ps) <- paste("p", nodes(tree))

    ## For entry / exit distributions
    bs.Fnorm <- bs.Fsub <- bs.ps; bs.Gnorm <- bs.Gsub <- bs.ps ## storage for BS entry/exit
    colnames(bs.Fnorm) <- colnames(bs.Fsub) <- paste("F", nodes(tree))
    colnames(bs.Gnorm) <- colnames(bs.Gsub) <- paste("G", nodes(tree))
    initial <- which(sapply(inEdges(tree), function(x) !length(x) > 0)) ## initial states, no Fs (entry)
    terminal <- which(sapply(edges(tree), function(x) !length(x) > 0)) ## terminal states, no Gs (exit)

    bs.var.sop <- matrix(0, nrow=length(et), ncol=ns) ## matrix for var matrix of state occup prob
    colnames(bs.var.sop) <- paste("Var", "p", nodes(tree))
    rownames(bs.var.sop) <- et

    res.array <- array(0, dim=c(ns^2, ns^2, length(et)),
                       dimnames=list(rows=paste(rep(nodes(tree), ns), sort(rep(nodes(tree), ns))),
                       cols=paste(rep(nodes(tree), ns), sort(rep(nodes(tree), ns))), time=et))

    for (b in 1:B) { ## randomly selects the indices

        ## Find the bootstrap sample, pull bs sample info from original data & put into data set
        bs <- sample(ids, n, replace=TRUE)
        bs <- factor(bs, levels=ids)
        bs.tab <- data.frame(table(bs)) ### table with the frequencies
        Data.bs <- merge(Data, bs.tab, by.x="id", by.y="bs") ## merging original data with bs freq table
        bs.id <- as.vector(unlist(apply(Data.bs[Data.bs$Freq>0, , drop=FALSE], 1, function(x) paste(x["id"], 1:x["Freq"], sep=".")))) ## creating bs id
        idx <- rep(1:nrow(Data.bs), Data.bs$Freq) ## indexing the bs sample
        Data.bs <- Data.bs[idx, ] ## creating a bs dataset
        Data.bs.originalID <- Data.bs$id
        Data.bs$id <- bs.id ## changing id column to bs.id to use functions
        Data.bs <- Data.bs[order(Data.bs$stop), ] ## ordered bs dataset

        ## Calling functions using bs dataset
        Cens <- Add.States(tree, LT)

        ## Here calculate start probs for BS data ...
        ## Data may be LT so need to account for that
        ## Put check here for whether have LT or not ...
        if (LT) {
            min.tran <- min(Data.bs$stop[Data.bs$end.stage!=0 & Data.bs$start.stage!="LT"])
            idx1 <- which(Data.bs$start.stage=="LT" & Data.bs$stop < min.tran)
            idx2 <- which(Data.bs$start.stage!="LT" & Data.bs$start < min.tran)
            start.stages <- by(Data[c(idx1, idx2), ], Data.bs$id[c(idx1, idx2)], function(x)
                               ifelse(min(x$start)<0, x$end.stage[which.min(x$start)],
                                      x$start.stage[which.min(x$start)]))
            start.stages <- factor(start.stages, levels=nodes(tree), labels=nodes(tree))
            start.cnts  <- table(start.stages)
            start.probs <- prop.table(start.cnts)
        } else {
                idx <- which(Data.bs$start < min(Data.bs$stop[Data.bs$end.stage!=0]))
                start.cnts  <- table(factor(Data.bs$start.stage[idx], levels=nodes(tree), labels=nodes(tree)))
                start.probs <- prop.table(start.cnts)
            }

        if (LT) {
            cp <- CP(tree, Cens$treeLT, Data.bs, Cens$nt.states.LT)
        }

        if (!LT) {
            cp <- CP(tree, Cens$tree0, Data.bs, Cens$nt.states)
        }

        ds.est <- DS(Cens$nt.states, cp$dNs, cp$sum_dNs,
                     cp$Ys, Cens="0", cens.type)
        cp.red <- Red(tree, cp$dNs, cp$Ys, cp$sum_dNs, ds.est$dNs.K,
                      ds.est$Ys.K, ds.est$sum_dNs.K)
        AJest <- AJ.estimator(ns, tree, cp.red$dNs.K, cp.red$Ys.K, start.probs)

        idx <- which(dimnames(bs.est)[[3]] %in% dimnames(AJest$I.dA)[[3]])
        idx2 <- which(!(dimnames(bs.est)[[3]] %in% dimnames(AJest$I.dA)[[3]]))
        bs.IA <- bs.est
        bs.IA[, , idx, b] <- AJest$I.dA
        bs.IA[, , idx2, b] <- diag(ns)

        bs.est[, , 1, b] <- bs.IA[, , 1, b]
        bs.ps[1, , b] <- start.probs%*%bs.est[, , 1, b]

        for (j in 2:length(et)) {
            bs.est[, , j, b] <- bs.est[, , j-1, b] %*% bs.IA[, , j, b]
            bs.ps[j, , b] <- start.probs%*%bs.est[, , j, b]
        } ## end of j for loop
        ## looks ok

        ## Entry / Exit variance as well
        for (i in 1:ns) {## looping through nodes
            node <- nodes(tree)[i]
            later.stages <- names(acc(tree, node)[[1]])
            stages <- c(node, later.stages)
            bs.f.numer <- rowSums(bs.ps[, paste("p", stages), b, drop=FALSE])
            ## adding bs.Fsub
            if (sum(bs.f.numer)==0)  bs.Fnorm[, i, b] <- bs.Fsub[, i, b] <- 0  else {
                bs.Fsub[, i, b] <- bs.f.numer
                bs.Fnorm[, i, b] <- bs.f.numer/bs.f.numer[length(bs.f.numer)]}
            if (length(stages)==1) next
            bs.g.numer <- rowSums(bs.ps[, paste("p", later.stages), b, drop=FALSE])
            if (sum(bs.g.numer)==0)  bs.Gnorm[, i, b] <- bs.Gsub[, i, b] <- 0 else {
                bs.Gsub[, i, b] <- bs.g.numer
                bs.Gnorm[, i, b] <- bs.g.numer/bs.g.numer[length(bs.g.numer)]}
        } ## end of for loop

    } ## end of bs loop

    ## Normalized entry / exit
    Fnorm.var <- apply(bs.Fnorm, c(1, 2), var)
    Fnorm.var[, initial] <- NA ## setting the initial state variances = NA
    Gnorm.var <- apply(bs.Gnorm, c(1, 2), var)
    Gnorm.var[, terminal] <- NA ## setting distn for terminal states to NA since don't exist

    ## Subdistribution entry / exit
    Fsub.var <- apply(bs.Fsub, c(1, 2), var)
    Fsub.var[, initial] <- NA
    Gsub.var <- apply(bs.Gsub, c(1, 2), var)
    Gsub.var[, terminal] <- NA

    ## SOP's
    bs.var.sop <- apply(bs.ps, c(1, 2), var)
    colnames(bs.var.sop) <- paste("Var", "p", nodes(tree))
    rownames(bs.var.sop) <- et

    state.names <- nodes(tree)
    trans.names <- paste(rep(state.names, ns), rep(state.names, each=ns))
    bs.cov <- array(dim=c(ns^2, ns^2, length(et)),
                    dimnames=list(rows=trans.names, cols=trans.names, time=et))
    bs.IdA.cov <- array(dim=c(ns^2, ns^2, length(et)),
                        dimnames=list(rows=trans.names, cols=trans.names, time=et))

    for (i in 1:length(et)) {
        bs.est.t <- matrix(bs.est[, , i, ], nrow=B, ncol=ns^2, byrow=TRUE)
        bs.IdA.t <- matrix(bs.IA[, , i, ], nrow=B, ncol=ns^2, byrow=TRUE)
        bs.cov[, , i] <- cov(bs.est.t)
        bs.IdA.cov[, , i] <- cov(bs.IdA.t)
    }  ## this for loop creates a B x (# of states)^2 x (# of event times)

    list(cov.AJs=bs.cov, var.sop=bs.var.sop, cov.dA=bs.IdA.cov,
         Fnorm.var=Fnorm.var, Gnorm.var=Gnorm.var, Gsub.var=Gsub.var, Fsub.var=Fsub.var)

} ## end of function




####################################################################
##          CONFIDENCE INTERVALS for p(t) & P(s, t)               ##
####################################################################

## x = msSurv object
MSM.CIs <- function(x, ci.level=0.95, ci.trans="linear", trans=TRUE, sop=TRUE) {

    if (ci.level < 0 || ci.level > 1)
        stop("confidence level must be between 0 and 1")

    z.alpha <- qnorm(ci.level + (1 - ci.level) / 2)

    ci.trans <- match.arg(ci.trans, c("linear", "log", "cloglog", "log-log"))

    ## CIs on state occup probability
    if (sop==TRUE) {
        CI.p <- array(0, dim=c(nrow(dNs(x)), 3, length(nodes(tree(x)))),
                      dimnames=list(rows=et(x),
                      cols=c("est", "lower limit", "upper limit"),
                      state=paste("p", nodes(tree(x)))))

        CI.p[ , 1, ] <- ps(x)
        switch(ci.trans[1],
               "linear" = {
                   CI.p[ , 2, ] <- ps(x) - z.alpha * sqrt(var.sop(x))
                   CI.p[ , 3, ] <- ps(x) + z.alpha * sqrt(var.sop(x))},
               "log" = {
                   CI.p[ , 2, ] <- exp(log(ps(x)) - z.alpha * sqrt(var.sop(x)) / ps(x))
                   CI.p[ , 3, ] <- exp(log(ps(x)) + z.alpha * sqrt(var.sop(x)) / ps(x))},
               "cloglog" = {
                   CI.p[ , 2, ] <- 1 - (1 - ps(x))^(exp(z.alpha * (sqrt(var.sop(x)) /
                                                                 ((1 - ps(x)) * log(1 - ps(x))))))
                   CI.p[ , 3, ] <- 1 - (1 - ps(x))^(exp(-z.alpha * (sqrt(var.sop(x)) /
                                                                  ((1 - ps(x)) * log(1 - ps(x))))))},
               "log-log" = {
                   CI.p[ , 2, ] <- ps(x)^(exp(-z.alpha * (sqrt(var.sop(x)) / (ps(x) * log(ps(x))))))
                   CI.p[ , 3, ] <- ps(x)^(exp(z.alpha * (sqrt(var.sop(x)) / (ps(x) * log(ps(x))))))})

        ## Need to loop through columns to apply pmax / pmin
        for (j in 1:length(nodes(tree(x)))) {
            CI.p[ , 2, j] <- pmax(CI.p[ , 2, j], 0)
            CI.p[ , 3, j] <- pmin(CI.p[ , 3, j], 1)
        } ## end states loop
    }

    ## CIs on transition probability matrices ##
    if (trans==TRUE) {
        CI.trans <- array(0, dim=c(nrow(dNs(x)), 4, length(pos.trans(x))),
                          dimnames = list(rows=et(x),
                          cols=c("est", "lower limit", "upper limit", "var.tp"),
                          trans=pos.trans(x)))

        for (j in 1:length(pos.trans(x))) { ## loop through possible transitions

            idx <- unlist(strsplit(pos.trans(x)[j], " "))
            CI.trans[ , 1, j] <- PE <- AJs(x)[idx[1], idx[2] , ]
            CI.trans[ , 4, j] <- var <- cov.AJs(x)[pos.trans(x)[j], pos.trans(x)[j], ]


            switch(ci.trans[1],
                   "linear" = {
                       CI.trans[ , 2, j] <- PE - z.alpha * sqrt(var)
                       CI.trans[ , 3, j] <- PE + z.alpha * sqrt(var)},
                   "log" = {
                       CI.trans[ , 2, j] <- exp(log(PE) - z.alpha * sqrt(var) / PE)
                       CI.trans[ , 3, j] <- exp(log(PE) + z.alpha * sqrt(var) / PE)},
                   "cloglog" = {
                       CI.trans[ , 2, j] <- 1 - (1 - PE)^(exp(z.alpha * (sqrt(var) /
                                                                       ((1 - PE) * log(1 - PE)))))
                       CI.trans[ , 3, j] <- 1 - (1 - PE)^(exp(-z.alpha * (sqrt(var) /
                                                                        ((1 - PE) * log(1 - PE)))))},
                   "log-log" = {
                       CI.trans[ , 2, j] <- PE^(exp(-z.alpha * (sqrt(var) / (PE * log(PE)))))
                       CI.trans[ , 3, j] <- PE^(exp(z.alpha * (sqrt(var) / (PE * log(PE)))))})

            CI.trans[ , 2, j] <- pmax(CI.trans[ , 2, j], 0)
            CI.trans[ , 3, j] <- pmin(CI.trans[ , 3, j], 1)

        } ## end j loop
    }

    if (trans==TRUE & sop==TRUE) {
        return(list(CI.p=CI.p, CI.trans=CI.trans))
    } else if (trans==TRUE) {
        return(list(CI.trans=CI.trans))
    } else if (sop==TRUE) {
        return(list(CI.p=CI.p))
    } else {
        cat("\nNothing to return \n\n")
    }

} ## end of function


####################################################################
##              CIs for Entry / Exit Functions                    ##
####################################################################


Dist.CIs <- function(x, ci.level=0.95, ci.trans="linear", norm=TRUE) {

    z.alpha <- qnorm(ci.level + (1 - ci.level) / 2)

    ci.trans <- match.arg(ci.trans, c("linear", "log", "cloglog", "log-log"))

    CI.Fs <- array(0, dim=c(length(et(x)), 3, length(nodes(tree(x)))),
                   dimnames=list(rows=et(x), cols=c("est", "lower limit", "upper limit"),
                   state=paste("F", nodes(tree(x)))))
    CI.Gs <- array(0, dim=c(length(et(x)), 3, length(nodes(tree(x)))),
                   dimnames=list(rows=et(x), cols=c("est", "lower limit", "upper limit"),
                   state=paste("G", nodes(tree(x)))))

    if (norm) {
	CI.Fs[,1,] <- Fnorm(x); CI.Gs[,1,] <- Gnorm(x)
        Fs.var <- Fnorm.var(x); Gs.var <- Gnorm.var(x)
    } else {
        CI.Fs[,1,] <- Fsub(x); CI.Gs[,1,] <- Gsub(x)
        Fs.var <- Fsub.var(x); Gs.var <- Gsub.var(x)
    }

    switch(ci.trans[1],
           "linear" = {
               CI.Fs[ , 2, ] <- CI.Fs[,1,] - z.alpha * sqrt(Fs.var)
               CI.Fs[ , 3, ] <- CI.Fs[,1,] + z.alpha * sqrt(Fs.var)
               CI.Gs[ , 2, ] <- CI.Gs[,1,] - z.alpha * sqrt(Gs.var)
               CI.Gs[ , 3, ] <- CI.Gs[,1,] + z.alpha * sqrt(Gs.var)},
           "log" = {
               CI.Fs[ , 2, ] <- exp(log(CI.Fs[,1,]) - z.alpha * sqrt(Fs.var) / CI.Fs[,1,])
               CI.Fs[ , 3, ] <- exp(log(CI.Fs[,1,]) + z.alpha * sqrt(Fs.var) / CI.Fs[,1,])
               CI.Gs[ , 2, ] <- exp(log(CI.Gs[,1,]) - z.alpha * sqrt(Gs.var) / CI.Gs[,1,])
               CI.Gs[ , 3, ] <- exp(log(CI.Gs[,1,]) + z.alpha * sqrt(Gs.var) / CI.Gs[,1,])},
           "cloglog" = {
               CI.Fs[ , 2, ] <- 1 - (1 - CI.Fs[,1,])^(exp(z.alpha * (sqrt(Fs.var) / ((1 - CI.Fs[,1,]) * log(1 - CI.Fs[,1,])))))
               CI.Fs[ , 3, ] <- 1 - (1 - CI.Fs[,1,])^(exp(-z.alpha * (sqrt(Fs.var) / ((1 - CI.Fs[,1,]) * log(1 - CI.Fs[,1,])))))
               CI.Gs[ , 2, ] <- 1 - (1 - CI.Gs[,1,])^(exp(z.alpha * (sqrt(Gs.var) / ((1 - CI.Gs[,1,]) * log(1 - CI.Gs[,1,])))))
               CI.Gs[ , 3, ] <- 1 - (1 - CI.Gs[,1,])^(exp(-z.alpha * (sqrt(Gs.var) / ((1 - CI.Gs[,1,]) * log(1 - CI.Gs[,1,])))))},
           "log-log" = {
               CI.Fs[ , 2, ] <- CI.Fs[,1,]^(exp(-z.alpha * (sqrt(Fs.var) / (CI.Fs[,1,] * log(CI.Fs[,1,])))))
               CI.Fs[ , 3, ] <- CI.Fs[,1,]^(exp(z.alpha * (sqrt(Fs.var) / (CI.Fs[,1,] * log(CI.Fs[,1,])))))
               CI.Gs[ , 2, ] <- CI.Gs[,1,]^(exp(-z.alpha * (sqrt(Gs.var) / (CI.Gs[,1,] * log(CI.Gs[,1,])))))
               CI.Gs[ , 3, ] <- CI.Gs[,1,]^(exp(z.alpha * (sqrt(Gs.var) / (CI.Gs[,1,] * log(CI.Gs[,1,])))))})

    for (j in 1:length(nodes(tree(x)))) {
        CI.Fs[ , 2, j] <- pmax(CI.Fs[ , 2, j], 0)
        CI.Fs[ , 3, j] <- pmin(CI.Fs[ , 3, j], 1)

        CI.Gs[ , 2, j] <- pmax(CI.Gs[ , 2, j], 0)
        CI.Gs[ , 3, j] <- pmin(CI.Gs[ , 3, j], 1)
    }

list(CI.Fs=CI.Fs, CI.Gs=CI.Gs)

} ## end of function
