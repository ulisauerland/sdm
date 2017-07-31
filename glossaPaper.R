
##
##  Glossa paper
##


# relies on functions in syncretism-analysis.R
setwd("/Volumes/Data/Users/us1/Documents/Projects/Syncretism")
source("syncretism-analysis.R",local=T)


# color selections from colorbrewer2.org
# Four colors:
col4 <- c('#d7191c','#fdae61','#abdda4','#2b83ba')

# Four colors, colorblind safe, print friendly, photoopy safe
col4 <- c('#e66101','#fdb863','#b2abd2','#5e3c99')

# eight colors
col8 <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666')

#  functions for vertical displays of features and paradigms

draw.paradigm.vertical <- function (p,x,y,w=.5)
{   pcol=col4
	size <- length(p)
	for (i in 1:size)
	{	rect(x,y+(i-1)*w,x+w,y+i*w,border=NA,col=pcol[p[size-i+1]])
	}
	rect(x,y,x+w,y+size*w,col=NA,border="black")
}

draw.paradigms.vertical <- function (P,x,y)
{   
	if (is.matrix(P))
	{	
	    for (i in 1:ncol(P))
		   draw.paradigm.vertical(P[,i],x-.6+.6*i,y,.5)
	}
	else
	{	    
		   draw.paradigm.vertical(P,x-.6+.6*1,y,.5)
	}	
}



draw.feature.vertical <- function (f,x,y,w,color="grey")
{   pcol=c("white",color)
	size <- length(f)
	for (i in 1:size)
	{	rect(x,y+(i-1)*w,x+w,y+i*w,border=NA,col=pcol[f[size-i+1]+1])
	}
	rect(x,y,x+w,y+size*w,col=NA,border="black")
}

draw.features.vertical <- function (F,x,y,colors=NULL)
{   
	if (is.matrix(F))
	{	if (is.null(colors))
	    for (i in 1:ncol(F))
		   draw.feature.vertical(F[,i],x-.6+.6*i,y,.5)
	else 	    for (i in 1:ncol(F))
		   draw.feature.vertical(F[,i],x-.6+.6*i,y,.5,colors[i])
	}
	else
	{	if (is.null(colors))
		   draw.feature.vertical(F,x-.6+.6*1,y,.5)
	else 	    
		   draw.feature.vertical(F,x-.6+.6*1,y,.5,colors[1])
	}
	
}



## TWO CELL CASE

cells<-2
features <- size_sort(all_features(cells))
paradigms <- all_paradigms(cells)
completeParadigm <- paradigms[,ncol(paradigms)]

fGood <- NULL
for (i in 1:(ncol(features)-1)) {     print(c(i)); for (j in (i+1):(ncol(features))) 
{
	if (!is.na(lexemesEO(features[,c(i,j)],completeParadigm)))
	      { fGood <- rbind(fGood,c(i,j)); print(c(i,j))
	}}
}

# result: 3 feature systems

m <- apply(paradigms,2,function(p)length(table(p)))
accHomph <- apply(fGood,1,
   function (r) 
   { sum( apply(paradigms,2,function(p) lexemesIO(features[,r],p)) -m) })

min(accHomph)


# plots for 2 cell case

setwd("~/Documents/Projects/Syncretism/GlossaGraphs")

w1 <- .5
w2 <- 1.15

pdf("features2-1.pdf", width=w1, height=.5)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,w2),c(-0.05,1.05))

draw.features.vertical(features[,fGood[1,]],0,0)
dev.off()

pdf("features2-2.pdf", width=w1, height=.5)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,w2),c(-0.05,1.05))

draw.features.vertical(features[,fGood[2,]],0,0)
dev.off()


pdf("features2-3.pdf", width=w1, height=.5)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,w2),c(-0.05,1.05))

draw.features.vertical(features[,fGood[3,]],0,0)
dev.off()


pdf("realization2-1-1.pdf", width=w1, height=.5)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,w2),c(-0.05,1.05))

draw.features.vertical(features[,fGood[1,]],0,0,col4)
dev.off()

pdf("realization2-2-1.pdf", width=w1, height=.5)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,w2),c(-0.05,1.05))

draw.features.vertical(features[,fGood[2,]],0,0,col4)
dev.off()

pdf("realization2-2-2.pdf", width=w1, height=.5)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,w2),c(-0.05,1.05))

draw.features.vertical(features[,fGood[2,2]],0,0,col4)
dev.off()


pdf("realization2-3-1.pdf", width=w1, height=.5)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,w2),c(-0.05,1.05))

draw.features.vertical(features[,fGood[3,]],0,0,col4)
dev.off()

pdf("realization2-3-2.pdf", width=w1, height=.5)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,w2),c(-0.05,1.05))

draw.features.vertical(features[,fGood[2,2]],0,0,col4)
dev.off()

pdf("paradigm2-1.pdf",width=.5*w1, height=.5)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,.55),c(-0.05,1.05))

draw.paradigms.vertical(paradigms[,1],0,0)
dev.off()

pdf("paradigm2-2.pdf",width=.5*w1, height=.5)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,.55),c(-0.05,1.05))

draw.paradigms.vertical(paradigms[,2],0,0)
dev.off()



######################################################################################


cells<-3
features <- size_sort(all_features(cells))
paradigms <- all_paradigms(cells)
completeParadigm <- paradigms[,ncol(paradigms)]
m <- apply(paradigms,2,function(p)length(table(p)))


setwd("~/Documents/Projects/Syncretism/GlossaGraphs")
w1 <- .85
w2 <- 1.75

pdf("features3.pdf", width=w1, height=.75)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,w2),c(-0.05,1.55))

draw.features.vertical(features[,c(4,6,2)],0,0,col4)
dev.off()

pdf("paradigm3-1.pdf",width=w1/3, , height=.75)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,.55),c(-0.05,1.55))

draw.paradigms.vertical(paradigms[,5],0,0)
dev.off()

pdf("realizations3-1.pdf", width=w1, height=.75)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,w2),c(-0.05,1.55))

draw.features.vertical(features[,c(2,4,6)],0,0,col4)
dev.off()

pdf("paradigm3-2.pdf",width=w1/3, height=.75)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,.55),c(-0.05,1.55))

draw.paradigms.vertical(paradigms[,2],0,0)
dev.off()

pdf("realizations3-2.pdf", width=w1, height=.75)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,w2),c(-0.05,1.55))

draw.features.vertical(features[,c(4,6)],0,0,col4)
dev.off()

pdf("paradigm3-3.pdf",width=w1/3, height=.75)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,.55),c(-0.05,1.55))

draw.paradigms.vertical(paradigms[,4],0,0)
dev.off()

pdf("realizations3-3.pdf", width=w1, height=.75)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(c(-0.05,w2),c(-0.05,1.55))

draw.features.vertical(features[,c(6,4)],0,0,col4)
dev.off()


# with only implicit ordering
# checking the two feature sets
fGood <- NULL
accHomph <- NULL
accHomphLI <- NULL
lex <- NULL
lexLI <- NULL

system.time(for (i in 1:(ncol(features)-1)) {     print(c(i)); for (j in (i+1):(ncol(features)))
{
	if (!is.na(lexemesEO(features[,c(i,j)],completeParadigm)))
	      { fGood <- rbind(fGood,c(i,j))
	      	print(c(i,j))
	      	le <- apply(paradigms,2,function(p) lexemesEO(features [,c(i,j)],p))
	      	accHomph <- rbind(accHomph, sum(le-m))
	      	lex <- rbind(lex,le)
	      	li <- apply(paradigms,2,function(p) lexemesIO(features [,c(i,j)],p))
	        accHomphLI <- rbind(accHomphLI, sum(li-m))
			lexLI <- rbind(lexLI,li)
	}}
})

# 0.3 seconds
# 3 feature systems



duplet <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) for (j in (i+1):nrow(fGood)) if (all(lex[i,]==lex[j,])) duplet[j] <- TRUE

# extrinsic order: no duplets
# intrinsic order: all duplets, only 1 genuine system

dupletLI <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) for (j in (i+1):nrow(fGood)) if (all(lexLI[i,]==lexLI[j,])) dupletLI[j] <- TRUE




#There are three possible systems having two features.  The all have the same theoretical rate of accidental homophony, namely 5.


# generate PDFs of all 16 valid partition spaces

setwd("~/Documents/Projects/Syncretism/GlossaGraphs")
for (i in 0:15)
{
pdf(paste("space",i,".pdf",sep=""),width=3, height=1)
	par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
	plot.new()
plot.window(c(-0.05,2.95),c(-0.05,1.55))

draw.paradigms.vertical(paradigms[,bindigit(16+i,1:5)==1],0,0)
dev.off()
}



# 
# checking the three feature sets with three cells
#

fGood <- NULL
accHomph <- NULL
accHomphLI <- NULL
lex <- NULL
lexLI <- NULL

for (i in 1:(ncol(features)-2)) {     print(c(i)); for (j in (i+1):(ncol(features)-1)) for (k in (j+1):ncol(features))
{
	if (!is.na(lexemesEO(features[,c(i,j,k)],completeParadigm)))
	      { fGood <- rbind(fGood,c(i,j,k))
	      	print(c(i,j,k))
	      	le <- apply(paradigms,2,function(p) lexemesEO(features [,c(i,j,k)],p))
	      	accHomph <- rbind(accHomph, sum(le-m))
	      	lex <- rbind(lex,le)
	      	li <- apply(paradigms,2,function(p) lexemesIO(features [,c(i,j,k)],p))
	        accHomphLI <- rbind(accHomphLI, sum(li-m))
			lexLI <- rbind(lexLI,li)
	}}
}

# 29 valid systems (out of 35 sets of 3 features)
nrow(fGood)

spacecounts <- rep(0,16)
for (i in 1:nrow(fGood))
{
 	spc <- sum(2^(0:4)[lex[i,]==m]) - 15
 	spacecounts[spc] <- spacecounts[spc]+1
}
spacecounts

spacecounts <- rep(0,16)
for (i in 1:nrow(fGood))
{
 	spc <- sum(2^(0:4)[lexLI[i,]==m]) - 15
 	spacecounts[spc] <- spacecounts[spc]+1
}
spacecounts




duplet <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) if (!duplet[i]) for (j in (i+1):nrow(fGood)) if (all(lex[i,]==lex[j,])) duplet[j] <- i

# extrinsic order: no duplets
# intrinsic order: all duplets, only 1 genuine system

dupletLI <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) if (!dupletLI[i]) for (j in (i+1):nrow(fGood)) if (all(lexLI[i,]==lexLI[j,])) dupletLI[j] <- i


plot.new()
plot.window(c(0,50),c(0,58))
title("Partition Spaces for Three Cells with Order")
number <- 1
crossindex <- rep(NA,nrow(fGood))
alternatives <- rep(0,nrow(fGood))

for (i in 1:nrow(fGood))
{
	 if (!duplet[i])
	 {
		draw.features.vertical(features[,fGood[i,]],9,(47-2*number))
		lpar <- (1:5)[lex[i,]==m]
		draw.paradigms.vertical(paradigms[,lpar],1,(47-2*number))
	 	crossindex[i] <- number
	 	number <- number+1
	 }
	 else
	 {
	 	draw.features.vertical(features[,fGood[i,]],13+4*alternatives[duplet[i]],(47-2*crossindex[duplet[i]]))
	 	alternatives[duplet[i]] <- alternatives[duplet[i]] + 1
	 }
}
for (i in 1:sum(!duplet))
  text(0,(47-2*i+.5),i,cex=.7)



plot.new()
plot.window(c(0,50),c(0,58))
title("Partition Spaces for Three Cells with Panini")
number <- 1
crossindex <- rep(NA,nrow(fGood))
alternatives <- rep(0,nrow(fGood))

for (i in 1:nrow(fGood))
{
	 if (!dupletLI[i])
	 {
		draw.features.vertical(features[,fGood[i,]],9,(47-2*number))
		lpar <- (1:5)[lexLI[i,]==m]
		draw.paradigms.vertical(paradigms[,lpar],1,(47-2*number))
	 	crossindex[i] <- number
	 	number <- number+1
	 }
	 else
	 {
	 	draw.features.vertical(features[,fGood[i,]],13+4*alternatives[dupletLI[i]],(47-2*crossindex[dupletLI[i]]))
	 	alternatives[dupletLI[i]] <- alternatives[dupletLI[i]] + 1
	 }
}
for (i in 1:sum(!dupletLI))
  text(0,(47-2*i+.5),i,cex=.7)



# 
# checking the four feature sets with three cells
#

fGood <- NULL
accHomph <- NULL
accHomphLI <- NULL
lex <- NULL
lexLI <- NULL

for (i in 1:(ncol(features)-3)) {     print(c(i)); for (j in (i+1):(ncol(features)-2)) for (k in (j+1):(ncol(features)-1))  for (l in (k+1):ncol(features))
{
	if (!is.na(lexemesEO(features[,c(i,j,k,l)],completeParadigm)))
	      { fGood <- rbind(fGood,c(i,j,k,l))
	      	print(c(i,j,k,l))
	      	le <- apply(paradigms,2,function(p) lexemesEO(features [,c(i,j,k,l)],p))
	      	accHomph <- rbind(accHomph, sum(le-m))
	      	lex <- rbind(lex,le)
	      	li <- apply(paradigms,2,function(p) lexemesIO(features [,c(i,j,k,l)],p))
	        accHomphLI <- rbind(accHomphLI, sum(li-m))
			lexLI <- rbind(lexLI,li)
	}}
}

nrow(fGood)
# any 4 feature choice is valid system



dupletLI <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) if (!dupletLI[i]) for (j in (i+1):nrow(fGood)) if (all(lexLI[i,]==lexLI[j,])) dupletLI[j] <- i

duplet <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) if (!duplet[i]) for (j in (i+1):nrow(fGood)) if (all(lex[i,]==lex[j,])) duplet[j] <- i

plot.new()
plot.window(c(0,66),c(0,58))
title("Partition Spaces for Three Cells, Four features with Order")
number <- 1
crossindex <- rep(NA,nrow(fGood))
alternatives <- rep(0,nrow(fGood))

for (i in 1:nrow(fGood))
{
	 if (!duplet[i])
	 {
		draw.features.vertical(features[,fGood[i,]],9,(47-2*number))
		lpar <- (1:5)[lex[i,]==m]
		draw.paradigms.vertical(paradigms[,lpar],1,(47-2*number))
	 	crossindex[i] <- number
	 	number <- number+1
	 }
	 else
	 {
	 	draw.features.vertical(features[,fGood[i,]],13+4*alternatives[duplet[i]],(47-2*crossindex[duplet[i]]))
	 	alternatives[duplet[i]] <- alternatives[duplet[i]] + 1
	 }
}
for (i in 1:sum(!duplet))
  text(0,(47-2*i+.5),i,cex=.7)

plot.new()
plot.window(c(0,66),c(0,58))
title("Partition Spaces for Three Cells, four features with Panini")
number <- 1
crossindex <- rep(NA,nrow(fGood))
alternatives <- rep(0,nrow(fGood))

for (i in 1:nrow(fGood))
{
	 if (!dupletLI[i])
	 {
		draw.features.vertical(features[,fGood[i,]],9,(47-2*number))
		lpar <- (1:5)[lex[i,]==m]
		draw.paradigms.vertical(paradigms[,lpar],1,(47-2*number))
	 	crossindex[i] <- number
	 	number <- number+1
	 }
	 else
	 {
	 	draw.features.vertical(features[,fGood[i,]],13+4*alternatives[dupletLI[i]],(47-2*crossindex[dupletLI[i]]))
	 	alternatives[dupletLI[i]] <- alternatives[dupletLI[i]] + 1
	 }
}
for (i in 1:sum(!dupletLI))
  text(0,(47-2*i+.5),i,cex=.7)


# 
# checking the five feature sets with three cells
#

fGood <- NULL
accHomph <- NULL
accHomphLI <- NULL
lex <- NULL
lexLI <- NULL

for (i in 1:(ncol(features)-4)) {     print(c(i)); for (j in (i+1):(ncol(features)-3)) for (k in (j+1):(ncol(features)-2))  for (l in (k+1):(ncol(features)-1)) for (ll in (l+1):(ncol(features)))
{
	if (!is.na(lexemesEO(features[,c(i,j,k,l,ll)],completeParadigm)))
	      { fGood <- rbind(fGood,c(i,j,k,l,ll))
	      	print(c(i,j,k,l,ll))
	      	le <- apply(paradigms,2,function(p) lexemesEO(features [,c(i,j,k,l,ll)],p))
	      	accHomph <- rbind(accHomph, sum(le-m))
	      	lex <- rbind(lex,le)
	      	li <- apply(paradigms,2,function(p) lexemesIO(features [,c(i,j,k,l,ll)],p))
	        accHomphLI <- rbind(accHomphLI, sum(li-m))
			lexLI <- rbind(lexLI,li)
	}}
}

nrow(fGood)

duplet <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) if (!duplet[i]) for (j in (i+1):nrow(fGood)) if (all(lex[i,]==lex[j,])) duplet[j] <- i

dupletLI <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) if (!dupletLI[i]) for (j in (i+1):nrow(fGood)) if (all(lexLI[i,]==lexLI[j,])) dupletLI[j] <- i

sum(duplet==0)
sum(dupletLI==0)


plot.new()
plot.window(c(0,90),c(0,58))
title("Partition Spaces for Three Cells, Five features with Order")
number <- 1
crossindex <- rep(NA,nrow(fGood))
alternatives <- rep(0,nrow(fGood))

for (i in 1:nrow(fGood))
{
	 if (!duplet[i])
	 {
		draw.features.vertical(features[,fGood[i,]],9,(47-2*number))
		lpar <- (1:5)[lex[i,]==m]
		draw.paradigms.vertical(paradigms[,lpar],1,(47-2*number))
	 	crossindex[i] <- number
	 	number <- number+1
	 }
	 else
	 {
	 	draw.features.vertical(features[,fGood[i,]],13+4*alternatives[duplet[i]],(47-2*crossindex[duplet[i]]))
	 	alternatives[duplet[i]] <- alternatives[duplet[i]] + 1
	 }
}
for (i in 1:sum(!duplet))
  text(0,(47-2*i+.5),i,cex=.7)

plot.new()
plot.window(c(0,66),c(0,58))
title("Partition Spaces for Three Cells, five features with Panini")
number <- 1
crossindex <- rep(NA,nrow(fGood))
alternatives <- rep(0,nrow(fGood))

for (i in 1:nrow(fGood))
{
	 if (!dupletLI[i])
	 {
		draw.features.vertical(features[,fGood[i,]],9,(47-2*number))
		lpar <- (1:5)[lex[i,]==m]
		draw.paradigms.vertical(paradigms[,lpar],1,(47-2*number))
	 	crossindex[i] <- number
	 	number <- number+1
	 }
	 else
	 {
	 	draw.features.vertical(features[,fGood[i,]],13+4*alternatives[dupletLI[i]],(47-2*crossindex[dupletLI[i]]))
	 	alternatives[dupletLI[i]] <- alternatives[dupletLI[i]] + 1
	 }
}
for (i in 1:sum(!dupletLI))
  text(0,(47-2*i+.5),i,cex=.7)


# 
# checking the six feature sets with three cells
#

fGood <- NULL
accHomph <- NULL
accHomphLI <- NULL
lex <- NULL
lexLI <- NULL

for (i in 1:(ncol(features))) 
{
	if (!is.na(lexemesEO(features[,-i],completeParadigm)))
	      { fGood <- rbind(fGood,-i)
	      	print(-i)
	      	le <- apply(paradigms,2,function(p) lexemesEO(features [,-i],p))
	      	accHomph <- rbind(accHomph, sum(le-m))
	      	lex <- rbind(lex,le)
	      	li <- apply(paradigms,2,function(p) lexemesIO(features [,-i],p))
	        accHomphLI <- rbind(accHomphLI, sum(li-m))
			lexLI <- rbind(lexLI,li)
	}
}

nrow(fGood)

duplet <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) if (!duplet[i]) for (j in (i+1):nrow(fGood)) if (all(lex[i,]==lex[j,])) duplet[j] <- i


dupletLI <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) if (!dupletLI[i]) for (j in (i+1):nrow(fGood)) if (all(lexLI[i,]==lexLI[j,])) dupletLI[j] <- i

sum(duplet==0)
sum(dupletLI==0)


plot.new()
plot.window(c(0,90),c(0,58))
title("Partition Spaces for Three Cells, Six features with Order")
number <- 1
crossindex <- rep(NA,nrow(fGood))
alternatives <- rep(0,nrow(fGood))

for (i in 1:nrow(fGood))
{
	 if (!duplet[i])
	 {
		draw.features.vertical(features[,fGood[i,]],9,(47-2*number))
		lpar <- (1:5)[lex[i,]==m]
		draw.paradigms.vertical(paradigms[,lpar],1,(47-2*number))
	 	crossindex[i] <- number
	 	number <- number+1
	 }
	 else
	 {
	 	draw.features.vertical(features[,fGood[i,]],13+4*alternatives[duplet[i]],(47-2*crossindex[duplet[i]]))
	 	alternatives[duplet[i]] <- alternatives[duplet[i]] + 1
	 }
}
for (i in 1:sum(!duplet))
  text(0,(47-2*i+.5),i,cex=.7)

plot.new()
plot.window(c(0,66),c(0,58))
title("Partition Spaces for Three Cells, Six features with Panini")
number <- 1
crossindex <- rep(NA,nrow(fGood))
alternatives <- rep(0,nrow(fGood))

for (i in 1:nrow(fGood))
{
	 if (!dupletLI[i])
	 {
		draw.features.vertical(features[,fGood[i,]],9,(47-2*number))
		lpar <- (1:5)[lex[i,]==m]
		draw.paradigms.vertical(paradigms[,lpar],1,(47-2*number))
	 	crossindex[i] <- number
	 	number <- number+1
	 }
	 else
	 {
	 	draw.features.vertical(features[,fGood[i,]],13+4*alternatives[dupletLI[i]],(47-2*crossindex[dupletLI[i]]))
	 	alternatives[dupletLI[i]] <- alternatives[dupletLI[i]] + 1
	 }
}
for (i in 1:sum(!dupletLI))
  text(0,(47-2*i+.5),i,cex=.7)







##################################################################################

cells<-4
features <- size_sort(all_features(cells))
paradigms <- all_paradigms(cells)
completeParadigm <- paradigms[,ncol(paradigms)]
m <- apply(paradigms,2,function(p)length(table(p)))


# 
# checking the three feature sets with four cells
#




#
#  rather than going through each n for having n many features manually
#  searchFeatures(n) finds the good features set with n features 
#

searchFeatures <- function(FLevel,LastFIndex=1,Findices=NULL)
{
	if (FLevel > -1)
	{	
		Coll <- NULL
		for (i in LastFIndex:(ncol(features)-FLevel))
		{
			Coll <- rbind(Coll,searchFeatures(FLevel-1,i+1,c(Findices,i)))
		}
		return(Coll)
	}
	else
	{
		if (!is.na(lexemesEO(features[,Findices],completeParadigm)))
	      { 
	      	print(Findices)
	      	le <- apply(paradigms,2,function(p) lexemesEO(features [,Findices],p))
	      	li <- apply(paradigms,2,function(p) lexemesIO(features [,Findices],p))
	      	return(data.frame(fGood=I(list(Findices)),accHomph=sum(le-m),
	      		accHomphLI=sum(li-m),lex=I(list(le)),lexLI=I(list(li))))
		 }	
	}		
}



NF <- 3
system.time(goodFeatures <- searchFeatures(2))

#the first column is a list containing the vector of lenght numberOfFeaturesPerSet
#
#features[goodFeatures[[1,1]]]
#features
#
# 



numberOfFeaturesPerSet <- 3

(goodFeatureCount<-nrow(goodFeatures))
# 140 valid systems out of 2730

duplet <- rep(FALSE,nrow(goodFeatures))
for (i in 1:(nrow(goodFeatures)-1)) if (!duplet[i]) for (j in (i+1):nrow(goodFeatures)) if (all(goodFeatures$lex[[i]]==goodFeatures$lex[[j]])) duplet[j] <- i
(differentPartitions<-sum(duplet==0))
#116 different partition sets

dupletLI <- rep(FALSE,nrow(goodFeatures))
for (i in 1:(nrow(goodFeatures)-1)) if (!dupletLI[i]) for (j in (i+1):nrow(goodFeatures)) if (all(goodFeatures$lexLI[[i]]==goodFeatures$lexLI[[j]])) dupletLI[j] <- i
(differentPartitionsLI<-sum(dupletLI==0))
#47 different partition sets



#  recode the lists in goodFeatures into arrays for older code to display the three feature
#  solutions for the four cell paradigms
fGood <- NULL
lex <- NULL
lexLI <- NULL
for (i in 1:nrow(goodFeatures)) 
{
	fGood <- rbind(fGood,as.vector(goodFeatures$fGood[[i]]))
	lex <- rbind(lex,as.vector(goodFeatures$lex[[i]]))
	lexLI <- rbind(lexLI,as.vector(goodFeatures$lexLI[[i]]))
}




plot.new()
plot.window(c(0,66),c(0,58))
title("Partition Spaces for Three Cells, Four features with Order")
number <- 1
crossindex <- rep(NA,nrow(fGood))
alternatives <- rep(0,nrow(fGood))

for (i in 1:nrow(fGood))
{
	 if (!duplet[i])
	 {
		draw.features.vertical(features[,fGood[i,]],9,(47-2*number))
		lpar <- (1:15)[lex[i,]==m]
		draw.paradigms.vertical(paradigms[,lpar],1,(47-2*number))
	 	crossindex[i] <- number
	 	number <- number+1
	 }
	 else
	 {
	 	draw.features.vertical(features[,fGood[i,]],13+4*alternatives[duplet[i]],(47-2*crossindex[duplet[i]]))
	 	alternatives[duplet[i]] <- alternatives[duplet[i]] + 1
	 }
}
for (i in 1:sum(!duplet))
  text(0,(47-2*i+.5),i,cex=.7)

plot.new()
plot.window(c(0,66),c(0,140))
title("Partition Spaces for Three Cells, four features with Panini")
number <- 1
crossindex <- rep(NA,nrow(fGood))
alternatives <- rep(0,nrow(fGood))

for (i in 1:nrow(fGood))
{
	 if (!dupletLI[i])
	 {
		draw.features.vertical(features[,fGood[i,]],9,(140-3*number))
		lpar <- (1:15)[lex[i,]==m]
		draw.paradigms.vertical(paradigms[,lpar],1,(140-3*number))
	 	crossindex[i] <- number
	 	number <- number+1
	 }
	 else
	 {
	 	draw.features.vertical(features[,fGood[i,]],13+4*alternatives[dupletLI[i]],(140-3*crossindex[dupletLI[i]]))
	 	alternatives[dupletLI[i]] <- alternatives[dupletLI[i]] + 1
	 }
}
for (i in 1:sum(!dupletLI))
  text(0,(140-3*i+.5),i,cex=.7)






# 
# checking the four feature sets with four cells
#

fGood <- NULL
accHomph <- NULL
accHomphLI <- NULL
lex <- NULL
lexLI <- NULL

for (i in 1:(ncol(features)-3)) {     print(c(i)); 
  for (j in (i+1):(ncol(features)-2)) 
    for (k in (j+1):(ncol(features)-1))  
      for (l in (k+1):ncol(features))
{
	if (!is.na(lexemesEO(features[,c(i,j,k,l)],completeParadigm)))
	      { fGood <- rbind(fGood,c(i,j,k,l))
	      	print(c(i,j,k,l))
	      	le <- apply(paradigms,2,function(p) lexemesEO(features [,c(i,j,k,l)],p))
	      	accHomph <- rbind(accHomph, sum(le-m))
	      	lex <- rbind(lex,le)
	      	li <- apply(paradigms,2,function(p) lexemesIO(features [,c(i,j,k,l)],p))
	        accHomphLI <- rbind(accHomphLI, sum(li-m))
			lexLI <- rbind(lexLI,li)
	}}
}

nrow(fGood)
# 1015 out of 15*14*13*12=32760 any 4 feature choice is valid system

duplet <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) if (!duplet[i]) for (j in (i+1):nrow(fGood)) if (all(lex[i,]==lex[j,])) duplet[j] <- i

sum(duplet==0)
#317 different partition sets

dupletLI <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) if (!dupletLI[i]) for (j in (i+1):nrow(fGood)) if (all(lexLI[i,]==lexLI[j,])) dupletLI[j] <- i

sum(dupletLI==0)
#239 different partition sets

# 
# checking the five feature sets with four cells
#

fGood <- NULL
accHomph <- NULL
accHomphLI <- NULL
lex <- NULL
lexLI <- NULL

for (ii in 1:(ncol(features)-4)) for (i in (ii+1):(ncol(features)-3)) {     print(c(ii,i)); for (j in (i+1):(ncol(features)-2)) for (k in (j+1):(ncol(features)-1))  for (l in (k+1):ncol(features))
{
	if (!is.na(lexemesEO(features[,c(ii,i,j,k,l)],completeParadigm)))
	      { fGood <- rbind(fGood,c(ii,i,j,k,l))
	      	print(c(ii,i,j,k,l))
	      	le <- apply(paradigms,2,function(p) lexemesEO(features [,c(ii,i,j,k,l)],p))
	      	accHomph <- rbind(accHomph, sum(le-m))
	      	lex <- rbind(lex,le)
	      	li <- apply(paradigms,2,function(p) lexemesIO(features [,c(ii,i,j,k,l)],p))
	        accHomphLI <- rbind(accHomphLI, sum(li-m))
			lexLI <- rbind(lexLI,li)
	}}
}

nrow(fGood)
#  2793 out of 15*14*13*12*11=360360 feature choices is a valid system

duplet <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) if (!duplet[i]) for (j in (i+1):nrow(fGood)) if (all(lex[i,]==lex[j,])) duplet[j] <- i

sum(duplet==0)
# 347 different partition sets

dupletLI <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) if (!dupletLI[i]) for (j in (i+1):nrow(fGood)) if (all(lexLI[i,]==lexLI[j,])) dupletLI[j] <- i

sum(dupletLI==0)
#402 different partition sets


# 
# checking the six feature sets with four cells
#

fGood <- NULL
accHomph <- NULL
accHomphLI <- NULL
lex <- NULL
lexLI <- NULL

for (iii in 1:(ncol(features)-5)) for (ii in (iii+1):(ncol(features)-4)) for (i in (ii+1):(ncol(features)-3)) {     print(c(ii,i)); for (j in (i+1):(ncol(features)-2)) for (k in (j+1):(ncol(features)-1))  for (l in (k+1):ncol(features))
{
	if (!is.na(lexemesEO(features[,c(iii,ii,i,j,k,l)],completeParadigm)))
	      { fGood <- rbind(fGood,c(iii,ii,i,j,k,l))
	      	print(c(iii,ii,i,j,k,l))
	      	le <- apply(paradigms,2,function(p) lexemesEO(features [,c(iii,ii,i,j,k,l)],p))
	      	accHomph <- rbind(accHomph, sum(le-m))
	      	lex <- rbind(lex,le)
	      	li <- apply(paradigms,2,function(p) lexemesIO(features [,c(iii,ii,i,j,k,l)],p))
	        accHomphLI <- rbind(accHomphLI, sum(li-m))
			lexLI <- rbind(lexLI,li)
	}}
}

nrow(fGood)
#  4935 out of 5005 feature choices is a valid system

duplet <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) if (!duplet[i]) for (j in (i+1):nrow(fGood)) if (all(lex[i,]==lex[j,])) duplet[j] <- i

sum(duplet==0)
# 310 different partition sets

dupletLI <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) if (!dupletLI[i]) for (j in (i+1):nrow(fGood)) if (all(lexLI[i,]==lexLI[j,])) dupletLI[j] <- i

sum(dupletLI==0)
#420 different partition sets




# GENERAL COMPUTATION OF TABLE FOR 3 to 15 FEATURES

cells<-4
features <- size_sort(all_features(cells))
paradigms <- all_paradigms(cells)
completeParadigm <- paradigms[,ncol(paradigms)]
m <- apply(paradigms,2,function(p)length(table(p)))


searchFeatures <- function(FLevel,LastFIndex,Findices)
{
	if (FLevel > -1)
	{	
		Coll <- NULL
		for (i in LastFIndex:(ncol(features)-FLevel))
		{
			Coll <- rbind(Coll,searchFeatures(FLevel-1,i+1,c(Findices,i)))
		}
		return(Coll)
	}
	else
	{
		if (!is.na(lexemesEO(features[,Findices],completeParadigm)))
	      { 
	      	print(Findices)
	      	le <- apply(paradigms,2,function(p) lexemesEO(features [,Findices],p))
	      	li <- apply(paradigms,2,function(p) lexemesIO(features [,Findices],p))
	      	return(data.frame(fGood=I(list(Findices)),accHomph=sum(le-m),
	      		accHomphLI=sum(li-m),lex=I(list(le)),lexLI=I(list(li))))
		 }	
	}		
}


dataFourCell <- function(NF)
{
numberOfFeaturesPerSet <- NF
goodFeatures <- searchFeatures(NF-1,1,NULL)

goodFeatureCount<-nrow(goodFeatures)

duplet <- rep(FALSE,nrow(goodFeatures))
for (i in 1:(nrow(goodFeatures)-1)) if (!duplet[i]) for (j in (i+1):nrow(goodFeatures)) if (all(goodFeatures$lex[[i]]==goodFeatures$lex[[j]])) duplet[j] <- i
differentPartitions<-sum(duplet==0)

dupletLI <- rep(FALSE,nrow(goodFeatures))
for (i in 1:(nrow(goodFeatures)-1)) if (!dupletLI[i]) for (j in (i+1):nrow(goodFeatures)) if (all(goodFeatures$lexLI[[i]]==goodFeatures$lexLI[[j]])) dupletLI[j] <- i
differentPartitionsLI<-sum(dupletLI==0)

return(list(numberOfFeaturesPerSet=numberOfFeaturesPerSet, goodFeatures=goodFeatures,goodFeatureCount=goodFeatureCount,
		differentPartitions=differentPartitions,differentPartitionsLI=differentPartitionsLI))
}

system.time(dataFourCellThree <- dataFourCell(3))
save(dataFourCellThree,"dataFourCellThree",file="dataFourCellThree",ascii=TRUE)
system.time(dataFourCellFour <- dataFourCell(4))
save(dataFourCellFour,"dataFourCellFour",file="dataFourCellFour",ascii=TRUE)
system.time(dataFourCellFive <- dataFourCell(5))
save(dataFourCellFive,"dataFourCellFive",file="dataFourCellFive",ascii=TRUE)
system.time(dataFourCellSix <- dataFourCell(6))
save(dataFourCellSix,"dataFourCellSix",file="dataFourCellSix",ascii=TRUE)
system.time(dataFourCellSeven <- dataFourCell(7))
save(dataFourCellSeven,"dataFourCellSeven",file="dataFourCellSeven",ascii=TRUE)
system.time(dataFourCellEight <- dataFourCell(8))
save(dataFourCellEight,"dataFourCellEight",file="dataFourCellEight",ascii=TRUE)
system.time(dataFourCellNine <- dataFourCell(9))
save(dataFourCellNine,"dataFourCellNine",file="dataFourCellNine",ascii=TRUE)
system.time(dataFourCellTen <- dataFourCell(10))
save(dataFourCellTen,"dataFourCellTen",file="dataFourCellTen",ascii=TRUE)
system.time(dataFourCellEleven <- dataFourCell(11))
save(dataFourCellEleven,"dataFourCellEleven",file="dataFourCellEleven",ascii=TRUE)
system.time(dataFourCellTwelve <- dataFourCell(12))
save(dataFourCellTwelve,"dataFourCellTwelve",file="dataFourCellTwelve",ascii=TRUE)
system.time(dataFourCellThirteen <- dataFourCell(13))
save(dataFourCellThirteen,"dataFourCellThirteen",file="dataFourCellThirteen",ascii=TRUE)
system.time(dataFourCellFourteen <- dataFourCell(14))
save(dataFourCellFourteen,"dataFourCellFourteen",file="dataFourCellFourteen",ascii=TRUE)
#dataFourCellFifteen <- dataFourCell(15)

load("dataFourCellThree")
load("dataFourCellFour")
load("dataFourCellFive")
load("dataFourCellSix")
load("dataFourCellSeven")
load("dataFourCellEight")
load("dataFourCellNine")
load("dataFourCellTen")
load("dataFourCellEleven")
load("dataFourCellTwelve")
load("dataFourCellThirteen")
load("dataFourCellFourteen")


dsyn <- function (n)
{
	return(switch(n,1,2,dataFourCellThree,
dataFourCellFour,
dataFourCellFive,
dataFourCellSix,
dataFourCellSeven,
dataFourCellEight,
dataFourCellNine,
dataFourCellTen,
dataFourCellEleven,
dataFourCellTwelve,
dataFourCellThirteen,
dataFourCellFourteen))
}

dsyn(14)$goodFeatureCount
dsyn(14)$differentPartitions
dsyn(14)$differentPartitionsLI

dsyn(10)$goodFeatureCount
dsyn(10)$differentPartitions
dsyn(10)$differentPartitionsLI

for (i in 3:14)
{
	print(c(i,dsyn(i)$goodFeatureCount,dsyn(i)$differentPartitions,dsyn(i)$differentPartitionsLI))
}

FC <- 1:15
for (i in 3:14) FC[i] <- dsyn(i)$goodFeatureCount
FC[1]<-0
FC[2]<-0
FC[15]<-1
sum(FC[1:15])

# Now find out sets of paradigms that aren't done by any FS

dsyn(3)$goodFeatures$lex
IX<-dsyn(3)$goodFeatures$lex[[140]]


setrep <- rep(FALSE,2^16)
repc <- rep(0,2^16)
repinJ <- rep(NA,2^16)
repinI <- rep(NA,2^16)
setrepLI <- rep(FALSE,2^16)
repcLI <- rep(0,2^16)
repinJLI <- rep(NA,2^16)
repinILI <- rep(NA,2^16)

powers <- 2^(0:14)
# check the code for 3-feature sequences
for (i in 1:nrow(dsyn(3)$goodFeatures))
{
	IX <- dsyn(3)$goodFeatures$lex[[i]]==m
	sIX <- sum(powers[IX])
	setrep[sIX]<-TRUE
	IX <- dsyn(3)$goodFeatures$lexLI[[i]]==m
	sIX <- sum(powers[IX])
	setrepLI[sIX]<-TRUE	
}

for (j in 3:14)
{
	ds <- dsyn(j)$goodFeatures
	for (i in 1:nrow(ds))
	{
	IX <- ds$lex[[i]]==m
	sIX <- sum(powers[IX])
	setrep[sIX]<-TRUE
	repc[sIX] <- repc[sIX]+1
	if (is.na(repinJ[sIX])) repinJ[sIX] <- j
	if (is.na(repinI[sIX])) repinI[sIX] <- i
	IX <- ds$lexLI[[i]]==m
	sIX <- sum(powers[IX])
	setrepLI[sIX]<-TRUE	
	repcLI[sIX] <- repcLI[sIX]+1
	if (is.na(repinJLI[sIX])) repinJLI[sIX] <- j
	if (is.na(repinILI[sIX])) repinILI[sIX] <- i
	}
}

sum(setrep)
sum(setrepLI)

sum(setrep | setrepLI)

repc[repc>0]
repcLI[repcLI>0]


(1:2^16)[setrep]
(1:2^16)[setrepLI]
as.hexmode((1:2^16)[setrep&setrepLI])
as.hexmode((1:2^16)[setrep&!setrepLI])
as.hexmode((1:2^16)[!setrep&setrepLI])

plot.new()
plot.window(xlim=c(0,2),ylim=c(0,2))
draw.paradigms.vertical(indexToParadigms(3072+16384),0,0)

indexToBoolVector <- function(IX)
{
	r<-rep(FALSE,15)
	for (i in 15:1)
	{
		if (powers[i]<=IX)
		{
			IX <- IX-powers[i];
			r[i] <- TRUE
		}
	}
	return(r)
}

indexToParadigms <- function(IX)
{
	return(paradigms[,indexToBoolVector(IX)])
}

# calculate the number of paradigms in systems generatable with/without Panini

calc <- (1:(2^16))[!setrep&!setrepLI]
calc <- calc[calc>=2^14 & calc<2^15]
ccalc <- calc
for (i in 1:length(calc)) ccalc[i] <- sum(indexToBoolVector(calc[i]))
ccalc
mean(ccalc)
calc[ccalc==3]

table(ccalc)
length(ccalc)

repcLI[calc[ccalc==2]]
as.hexmode(calc[ccalc==2])
as.hexmode(calc[ccalc==14])
as.hexmode(calc[ccalc==9])

#number of partitions in set for the systems that can be generated by only 2 feature systems
tapply((1:2^16)[repcLI==2],(1:2^16)[repcLI==2],function (x) sum(indexToBoolVector(x)))

### graphics for some examples

setwd("~/Documents/Projects/Syncretism/GlossaGraphs")



examplePdfs <- function(s)
{
	pdf(paste("fourcell",s,"Paradigms.pdf"), width=8.1/4, height=2.2/4)
	par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
	plot.new()
	plot.window(xlim=c(-0.1,8),ylim=c(-.1,2.1))
	inte <- strtoi(s,base=16)
	draw.paradigms.vertical(indexToParadigms(inte),0,0)
	dev.off()
	
	if (!is.na(repinJLI[inte]))
	{
	pdf(paste("fourcell",s,"PaniniFeatures.pdf"), width=7.6/4, height=2.2/4)
	par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
	plot.new()
	plot.window(xlim=c(-0.1,7.5),ylim=c(-.1,2.1))
	draw.features.vertical(features[,dsyn(repinJLI[inte])$goodFeatures$fGood[[repinILI[inte]]]],0,0)
	dev.off()
	}
	
	if (!is.na(repinJ[inte]))
	{
	pdf(paste("fourcell",s,"Features.pdf"), width=7.6/4, height=2.2/4)
	par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
	plot.new()
	plot.window(xlim=c(-0.1,7.5),ylim=c(-.1,2.1))
	draw.features.vertical(features[,dsyn(repinJ[inte])$goodFeatures$fGood[[repinI[inte]]]],0,0)
	dev.off()
	}

	return(c(repc[inte],repcLI[inte]))
}

examplePs <- function(s)
{
	pdf(paste("fourcell",s,"Paradigms.pdf",sep="_"), width=8.1/4, height=2.2/4)
	par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
	plot.new()
	plot.window(xlim=c(-0.1,8),ylim=c(-.1,2.1))
	inte <- strtoi(s,base=16)
	draw.paradigms.vertical(indexToParadigms(inte),0,0)
	dev.off()
}

examplePs("45db")
examplePs("7ff4")
examplePs("4c00")

examplePs("70db")
examplePs("7ddf")

examplePs("4000")
examplePs("6cba")
examplePs("7f00")

examplePs("4001")
examplePs("7770")
examplePs("7dfe")


# testing for *ABBA

# paradigms permuting ABBA: 4, 7, 9 (reversed order 0x0148,  )

examplePdfs("0148")

ib328 <- indexToBoolVector(328)
abba <- rep(FALSE,2^15)
for (i in 1:2^14)
{
	ts <- indexToBoolVector(i) & ib328
	if (ts[4]&ts[7]&ts[9])
	{
		abba[i] <- TRUE;
		abba[i+2^14] <- TRUE
	}
}

sum((setrep|setrepLI)&abba)

calc <- (1:(2^16))[(setrep|setrepLI)&abba]
calc <- calc[calc>=2^14 & calc<2^15]
ccalc <- calc
for (i in 1:length(calc)) ccalc[i] <- sum(indexToBoolVector(calc[i]))
ccalc

calc

as.hexmode((1:2^15)[setrepLI&abba][sort(ccalc,index.return=T)$ix])

as.hexmode(calc[repinJ[calc]==3])
as.hexmode(calc[repinJLI[calc]==4])

examplePdfs("7dfe")
examplePdfs("45db")


# some on seven cell case

cells<-7
features <- size_sort(all_features(cells))
paradigms <- all_paradigms(cells)
completeParadigm <- paradigms[,ncol(paradigms)]
m <- apply(paradigms,2,function(p)length(table(p)))

Findices <- c(64,73,84)

le <- apply(paradigms,2,function(p) lexemesEO(features [,Findices],p))
li <- apply(paradigms,2,function(p) lexemesIO(features [,Findices],p))

sum(le==m)

pdf(paste("sevencellParadigms.pdf"), width=8.1/4*4, height=2.2/4*4)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(xlim=c(-1,38),ylim=c(-.1,3.6))
draw.paradigms.vertical(paradigms[,(le==m)],0,0)
dev.off()

pdf(paste("sevencellFeatures.pdf"), width=1.1/4*4, height=2.2/4*4)
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot.new()
plot.window(xlim=c(-1,3),ylim=c(-.1,3.6))
draw.features.vertical(features [,Findices],0,0)
dev.off()




###########################################
## for REVISION
## Computation without and_complete
## new function lexemesIOIless in syncretism-analysis.R
##


cells<-3
features <- size_sort(all_features(cells))
paradigms <- all_paradigms(cells)
completeParadigm <- paradigms[,ncol(paradigms)]
m <- apply(paradigms,2,function(p)length(table(p)))


lexemesEOIless(features[,c(1,5,6,7)],completeParadigm)
is.na(lexemesIOIless(features[,c(5,6,7)],completeParadigm))

searchFeaturesIless <- function(FLevel,LastFIndex,Findices)
{
#  print(Findices)
  if (FLevel > -1)
  {	
    Coll <- NULL
    for (i in LastFIndex:(ncol(features)-FLevel))
    {
      Coll <- rbind(Coll,searchFeaturesIless(FLevel-1,i+1,c(Findices,i)))
    }
    return(Coll)
  }
  else
  {
    if (!is.na(lexemesIOIless(features[,Findices],completeParadigm)))
    { 
#      print(Findices)
      li <- apply(paradigms,2,function(p) lexemesIOIless(features [,Findices],p))
      return(data.frame(fGood=I(list(Findices)),
                        accHomphLIIless=sum(li-m),lexLIIless=I(list(li))))
    }	
  }		
}


searchFeaturesIless(2,1,NULL)

goodFeatures <- searchFeaturesIless(2,1,NULL)

#the first column is a list containing the vector of length numberOfFeaturesPerSet
#
#features[goodFeatures[[1,1]]]
#features
#
# 

fGood <- NULL
lex <- NULL
for (i in 1:nrow(goodFeatures))
{
  fGood <- rbind(fGood,unlist(goodFeatures$fGood[i]))
  lex <- rbind(lex,unlist(goodFeatures$lexLIIless[i]))
}

duplet <- rep(FALSE,nrow(fGood))
for (i in 1:(nrow(fGood)-1)) if (!duplet[i]) for (j in (i+1):nrow(fGood)) if (all(lex[i,]==lex[j,])) duplet[j] <- i

plot.new()
plot.window(c(0,66),c(0,58))
title("Partition Spaces for Three Cells, Four features with Order")
number <- 1
crossindex <- rep(NA,nrow(fGood))
alternatives <- rep(0,nrow(fGood))

for (i in 1:nrow(fGood))
{
  if (!duplet[i])
  {
    if (is_and_complete(features[,fGood[i,]]))
      draw.features.vertical(features[,fGood[i,]],9,(47-2*number), rep("pink",7))
    else
      draw.features.vertical(features[,fGood[i,]],9,(47-2*number), rep("grey",7))
    lpar <- (1:5)[lex[i,]==m]
    draw.paradigms.vertical(paradigms[,lpar],1,(47-2*number))
    crossindex[i] <- number
    number <- number+1
  }
  else
  {
    if (is_and_complete(features[,fGood[i,]]))
      draw.features.vertical(features[,fGood[i,]],13+4*alternatives[duplet[i]],47-2*crossindex[duplet[i]],rep("pink",7))
    else
      draw.features.vertical(features[,fGood[i,]],13+4*alternatives[duplet[i]],47-2*crossindex[duplet[i]],rep("grey",7))
  
    alternatives[duplet[i]] <- alternatives[duplet[i]] + 1
  }
}
for (i in 1:sum(!duplet))
  text(0,(47-2*i+.5),i,cex=.7)

