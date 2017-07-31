## GENERAL FUNCTIONS

bindigit <- function (n,d) (n %% 2^d) %/% 2 ^ (d-1)

all_features = function(cells)
 { r = NULL;
    f = vector(length=cells);
    for (i in 1:(2^cells-1)) { 
    	   for (j in 1:cells) f[j] = bindigit(i,j);
	   r <- cbind(r,f) }
    return(r)
    }

all_paradigms = function (n)
{ r = NULL;
  p=rep(1,n);
  if (n>=2) 
  { for (i in 1:n^(n-1))
  {  	flag = TRUE;
  	for (j in 2:n)
  	{
  	if (p[j] >= 2)
  	{
  	  if (sum(p[1:j-1]==p[j]-1) == 0) { flag = FALSE; break}
  	}
  	}
  	if (flag) r=cbind(r,p)
  	for (j in n:1)
  	   if (p[j] < n) {p[j] <- p[j]+1; break}
  	   else {p[j] <- 1}  
  } }
  return(r)
}



# lmexpected is called with three parameters
# c: the number of lexemes (cells)
# m: the number of morphemes 
# s: the rate of syncretism
# it returns the probability of a paradigm with c lexemes and m morphemes
# when the likelihood of syncretism is s
# s^{c-m} \prod_{i=1}^{m-1}\max((1-is),0)
# (there are more than one paradigm with c lexemes and m morphemes
#  except for m=c and m=1)
# if s is high, 1-s in the above may be smaller than 0, so the likelihood is 0
# an adjustment of the probabilities is necessary to make sure that
# the sum of lmexpected(c,m,s)*sum(paradigms==m) = 1

# When lmexpected is called with two arrays, c the array of lexeme counts,
# and m the array of morpheme counts, plus the syncretism level s
# lmexpected adjusts according to the above by scaling it to sum 1
# so when calling lmexpected on arrays one has to make sure to call it
# with complete arrays for one lexeme number.  I.e. there must be as
# many (c[i],m[i]) pairs as there are ways of getting an m[]

# usually c will be a constant vector when using matrices
# but the function also handles the case that c has entries for different
# lexeme numbers, then it should be complete for each lexeme number

# for example if there are three lexemes there are 3 2-morpheme possibilities
# so for s>0.5 the correct call is
# lmexpected(rep(3,5),c(1,2,2,2,3),.6)[1] = 0.333333
# not
# lmexpected(3,1,.6) = 0.36 which is unadjusted

lmexpected <- function (c,m,s)
 {  
   r <- s^(c-m);
   if (length(m)==1)
       for (i in 1:m-1)  r = r*max(0,1-i*s)
   else
       { for (j in 1:length(m)) if (m[j] > 1) for (i in 1:(m[j]-1)) { r[j] = r[j]*max(0,1-i*s) }
       	rr <- r*(r>=0)
       	sr <- sum(rr)
       	if (sr!=1) r <- rr/sr }
     return(r) 
}



# 2^n-1 poss. features 
# select up to k of those in fbase
# test whether the selection is *-neg-complete [for not two cells alls features yield the same value]: n(n-1)/2 2 cell combos, k features
# generate the *-definable set from these (up to 2^k): fplus
# orderings of fplus in forder
#     select first f_0, remainder = 1-f, next f must have sum(f * remainder) >=1
# matched paradigm: f_0, remainder_0 * f_1, ...  


# lexemesEO (F, p)
# F = array of features in columns
# p = a paradigm
# returns the number of lexemes needed by analysing p in F


lexemesEO_internal <- function (F,p,mask,av)
   # mask: vector that contains 0's for cells already used in analysis
   # av: TRUE for features still unused
{
  	  matcher <- apply(F,2,function(f) f*p*mask);
      partial_match <- apply(matcher,2,function(c) length(table(c[c!=0]))==1);
      cell_matched <- vector(mode="numeric",length=ncol(F));
      total_match <- rep(FALSE ,ncol(F))
      for (i in 1:ncol(F))
       { if (av[i] && partial_match[i])
          { cell_matched[i] <- max(matcher[,i]);
            total_match[i] <- 1== 
               all((matcher[,i] == cell_matched[i])==(p*mask==cell_matched[i]))
          } }
          # if there are total_matches:
          #      use the totally_matching cells and recurse on remainder
       if (sum(total_match)>0)
          {
          	new_mask <-mask
          	for (i in (1:nrow(F))[mask==1])
          	     new_mask[i] <- as.numeric(sum(as.matrix(F[i,total_match])) == 0)
          	if (sum(new_mask)!=0)
               	return(lexemesEO_internal(F,p,new_mask,av & ! total_match)
                        	+length(table(cell_matched[total_match])))
            else
                   return(length(table(cell_matched[total_match])))
         }
              # if there are no total matches:
              #  try out all partial matches, and recurse
         else
         {
          if (sum(partial_match)>0)
           { lex_count <- rep(NA,ncol(F))
          	 for (i in 1:ncol(F))
          	      {  if (partial_match[i] && av[i])
            	        lex_count[i] <- 
       	           lexemesEO_internal(F,p,(1-F[,i])*mask, av & 1:ncol(F)!=i) }
             if (sum(!is.na(lex_count[partial_match])) > 0)
               { return(min(lex_count[partial_match],na.rm=TRUE) + 1) }
             else return(NA)
            }
           else
             return(NA)
         }                    
}

#
# This function computes the number of lexemes required to analyze paradigm p
# using the features in F and conjunctions of those features.
# This function assumes that extrinsic ordering of lexemes is possible.
#

lexemesEO <- function(F,p)
{
    Fplus <- size_sort(and_complete(F));
    return(lexemesEO_internal(Fplus,p,rep(1,nrow(F)),rep(TRUE,ncol(Fplus))))
}



# lexemesIO (F, p)
# F = array of features in columns
# p = a paradigm
# returns the number of morphemes needed by analysing p in F

# This function computes the number of lexemes required to analyze paradigm p
# using the features in F and conjunctions of those features.
# This function assumes that only intrinsic ordering of lexemes is possible.
# It interprets unlike lexemesIOL intrinsic ordering strictly, so
# it's prohibited to use to partially overlapping lexemes even if their overlap cases
# are blocked by a more specific feature (ACTUALLY THIS only applies if one of the overlapping lexemes is total matching a part of a paradigm)                


lexemesIOS_internal <- function (F,p,mask,av)
   # mask: vector that contains 0's for cells already used in analysis
   # av: TRUE for features still unused
{
 matcher <- apply(F,2,function(f) f*p*mask);
 partial_match <- apply(matcher,2,function(c) length(table(c[c!=0]))==1);
 if (sum(partial_match)>0)
 { 
  lex_count <- rep(NA,ncol(F))
  for (i in (1:ncol(F))[partial_match & av])
  {# when trying partial match i, block features partially overlapping with i
   # even if the overlap is outside the mask
   ol <-  (F[,i]) %*% F 
   fs <- sum(F[,i])
   new_av <- av & (ol == fs | ol == 0)
   new_av[i]<-FALSE
   new_mask <- (1-F[,i])*mask
   if (sum(new_mask)>=1)
    lex_count[i] <- lexemesIOS_internal(F,p,new_mask, new_av)
   else
    lex_count[i] <- 0
  }
  if (all(is.na(lex_count[partial_match])))
   return(NA)
  else 	
  	return(min(lex_count[partial_match],na.rm=TRUE) + 1)
  }
 else
  return(NA)
} 

lexemesIOS <- function(F,p)
{
    Fplus <- size_sort(and_complete(F));
    return(lexemesIOS_internal(Fplus,p,rep(1,nrow(F)),rep(TRUE,ncol(Fplus))))
}

#
# This function computes the number of lexemes required to analyze paradigm p
# using the features in F and conjunctions of those features.
# This function assumes that only intrinsic ordering of lexemes is possible.
# It interprets unlike lexemesIOS intrinsic ordering lax, so
# it allows partially overlapping lexemes their overlap case is
# blocked by a more specific feature 
# (I think this is similar to extrinsic ordering with adding syncretism for
# overlap cases since, adding a lexeme for F*G resolves an overlap between F*G)                
# lexemesIO (F, p)
# F = array of features in columns
# p = a paradigm
# returns the number of morphemes needed by analysing p in F


lexemesIO_internal <- function (F,p,mask,av)
   # mask: vector that contains 0's for cells already used in analysis
   # av: TRUE for features still unused
{
# print(c(F,p,mask,av))
 matcher <- apply(F,2,function(f) f*p*mask);
 partial_match <- apply(matcher,2,function(c) length(table(c[c!=0]))==1);
 if (sum(partial_match)>0)
 { 
   lex_count <- rep(NA,ncol(F))
   for (i in (1:ncol(F))[partial_match & av])
   {# when trying partial match i, block features partially overlapping with i
    ol <-  (mask*F[,i]) %*% F 
    fs <- sum(F[,i]*mask)
    new_av <- av & (ol == fs | ol == 0)
    new_av[i]<-FALSE
    new_mask <- (1-F[,i])*mask
    if (sum(new_mask)>=1)
       lex_count[i] <- lexemesIO_internal(F,p,new_mask,new_av) 
    else
       lex_count[i] <- 0 
   }
   if (all(is.na(lex_count[partial_match])))
    return(NA)
   else 
    {	
    	return(min(lex_count[partial_match],na.rm=TRUE) + 1)
    }
  }
  else
   return(NA)
}

#
lexemesIO <- function(F,p)
{
    Fplus <- size_sort(and_complete(F));
    return(lexemesIO_internal(Fplus,p,rep(1,nrow(F)),rep(TRUE,ncol(Fplus))))
}

lexemesIOIless <- function(F,p)
{
  Fplus <- size_sort(F);
  return(lexemesIO_internal(Fplus,p,rep(1,nrow(F)),rep(TRUE,ncol(Fplus))))
}

lexemesEOIless <- function(F,p)
{
  Fplus <- size_sort(F);
  return(lexemesEO_internal(Fplus,p,rep(1,nrow(F)),rep(TRUE,ncol(Fplus))))
}


#generate F*, the completion of the feature set F under feature conjunction *

and_complete <- function(F)
{
new_nonsingleton <- TRUE; 
     # a flag; if there's a new non-singleton there are possibly still new combinations
previous_count <- ncol(F);
     # in later loops we only combine new features with old ones
start <- 1;
while (new_nonsingleton)
{ 
  new_nonsingleton <- FALSE;
  for (i in 1:ncol(F))
  if (sum(F[,i])>=2 && i<ncol(F))
   { for (j in max(i+1,start):ncol(F))
   	if (sum(F[,j])>=2)
   	  { c <- F[,i]*F[,j];
   	     s<-sum(c);
   	     if (s>0)
   	        { candidate_new <- TRUE;
   	        	 for (k in 1:ncol(F))
     	              { if (all(F[,k] == c)) {candidate_new <- FALSE; break} }
     	           if (candidate_new) F <- cbind(F,c);
     	        if (s>1 && candidate_new) new_nonsingleton <- TRUE;
                 } 
               } 
    }
    start <- previous_count;
    previous_count <-ncol(F);
}
return(F)
} 	


is_and_complete <- function(F)
{
  r <- TRUE ;
  for (i in 1:ncol(F))
      if (sum(F[,i])>=2 && i<ncol(F))
      { for (j in (i+1):ncol(F))
        if (sum(F[,j])>=2)
        { c <- F[,i]*F[,j];
          if (sum(c)>0)
          { 
            occurs <- FALSE
            for (k in 1:ncol(F))
              { if (all(F[,k] == c)) {occurs <- TRUE} }
            if (!occurs) {r <- FALSE}
          } 
        } 
      }
  return(r)
} 	


size_sort <- function (F)
{
	return(F[,order(apply(F,2,sum))])
	}
	




draw.paradigm <- function (p,x,y,pcol=c("red","pink","yellow","purple"))
{   
	w<-.5
	rect(x,y,x+w,y+w,border=NA,col=pcol[p[3]])
	rect(x+w,y,x+2*w,y+w,border=NA,col=pcol[p[4]])
	rect(x,y+w,x+w,y+2*w,border=NA,col=pcol[p[1]])
	rect(x+w,y+w,x+2*w,y+2*w,border=NA,col=pcol[p[2]])
	rect(x,y,x+2*w,y+2*w, col=NA,border="black")	
}

draw.feature <- function (f,x,y,w)
{   pcol=c("white","grey")
	rect(x,y,x+w,y+w,,border=NA,col=pcol[f[3]+1])
	rect(x+w,y,x+2*w,y+w,border=NA,col=pcol[f[4]+1])
	rect(x,y+w,x+w,y+2*w,border=NA,col=pcol[f[1]+1])
	rect(x+w,y+w,x+2*w,y+2*w,border=NA,col=pcol[f[2]+1])
	rect(x,y,x+2*w,y+2*w, col=NA,border="black")
}

draw.features <- function (F,x,y)
{   
	if (ncol(F)<=7)
	    for (i in 1:ncol(F))
		   draw.feature(F[,i],x-1.2+1.2*i,y,.5)
	else
	{	for (i in 1:7)
		   draw.feature(F[,i],x-1.2+1.2*i,y,.5)
		for (i in 8:ncol(F))
		    draw.feature(F[,i],x-1.2+1.2*(i-7),y+1.2,.5) }
}








