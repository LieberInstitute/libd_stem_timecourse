matchGenes_fix <- function(x,subject, type=c("any","fiveprime"),
                       promoterDist=2500,
                       skipExons=FALSE,
                       verbose=TRUE){

    if(attributes(subject)$description!="annotatedTranscripts")
        stop("subject must be the output of function annotateTranscripts or have attribute(subject)$description=\"annotatedTranscripts\".")
    
    if (is.data.frame(x)) x <- makeGRangesFromDataFrame(x)
    if(class(x)!="GRanges") stop("x must be GRanges or a data.frame with column names chr, start and end.")

    type=match.arg(type)
    
    if(type=="fiveprime"){
        map <- nearest(x,resize(subject,width=1))
    } else{
        map <- nearest(x,subject)
    }
    
    ind <- which(!is.na(map))
#    dist <- rep(NA,length(map))
#    dist[ind] <- distance(x[ind,],subject[map[ind],])
#    
#    dist[ind] <- dist[ind]*ifelse(
#        (strand(subject[map[ind],])=="+" &
#         end(subject[map[ind],]) < start(x[ind,])) |
#        (strand(subject[map[ind],])=="-" &
#         start(subject[map[ind],]) > end(x[ind,])),-1,1)
    
    type=rep("",length(x))
    subtype=rep("",length(x))
    ctype=rep("",length(x))
    dist=rep(0,length(x))	# distance to 5' end of the gene
    insidedistance<-rep(NA,length(x))
    exonnumber<-rep(NA,length(x))
    nexons <- rep(NA,length(x))
    geneL=rep(0,length(x))
    codingL=rep(0,length(x))
    subdist=rep(0,length(x))
    genenames=rep("",length(x))
    geneannotation=rep("",length(x))
    geneid=rep("",length(x))
    strands=rep("",length(x))

    
    genenames[ind]<-as.character(values(subject[map[ind],])$Gene)
    geneannotation[ind]<-as.character(values(subject[map[ind],])$Refseq)
    geneid[ind]<-as.character(values(subject[map[ind],])$Geneid)
    nexons[ind] <- values(subject[map[ind],])$Nexons
    strands[ind] <- as.character(strand(subject[map[ind],]))
    
    ## This loop could be moved to C. Note: it uses nearest
    for(j in ind){
        
        i <- map[ind][j]
        
        if(verbose & j%%100==0) cat(".")
        
        TS = start(subject)[i]
        TE = end(subject)[i]
        geneL[j] = TE-TS
        if(!is.na(subject$CSS[i])){
            CS = subject$CSS[i]
            CE = subject$CSE[i]
            codingL[j]=CE-CS
        } else {
            CS <- CE <- codingL <- NA
        }
        exons <- subject[i,]$Exons[[1]]
        Exons <- cbind(start(exons), end(exons))
        
        Strand= ifelse(strand(subject[i,])=="+",1,-1)
        S = start(x[j,])
        E = end(x[j,])
        
#        type[j]=""
        
        if(S <= TS & E >= TE){
            type[j]="covers"
            subtype[j]="covers exon(s)"
            
        } else{
            if(E < TS){
                if(Strand==1){
                    type[j]="upstream" 
                    dist[j]=TS-E
                } else{
                    type[j]="downstream"
                    dist[j]=TE-E
                }
            }
            if(S > TE){
                if(Strand==-1){
                    type[j]="upstream"
                    dist[j]=S-TE
                }  else{
                    type[j]="downstream"
                    dist[j]=S-TS
                }
            }
            ## totally within gene
            if (S >= TS & E <= TE){
                type[j]="inside"
                if(Strand==-1) dist[j]=TE-E  else dist[j]=S-TS
            }
            ## overlaps exactly one side of gene ("covers" done above)
            if(type[j]==""){
                if(S < TS & E <= TE){
                    ##OVERLAP FRONT
                    if(Strand==1) type[j]="overlaps 5'" else{
                        type[j]="overlaps 3'"
                        dist[j]=TE-E
                    }
                }
                else if (S >= TS & E > TE){
                    ##OVERLAP BACK
                    if(Strand==-1) type[j]="overlaps 5'" else{
                        type[j]="overlaps 3'"
                        dist[j]=S-TS
                    }
                }
            }
        }
        
        m1=NA;m2=NA
        if( type[j]%in%c("overlaps 5'","overlaps 3'","inside") & !skipExons ){
            
            ir=IRanges(start=c(S,E),width=1)
            map2 = nearest(ir, exons)
            dist2<-distance(ir,exons[map2])
            pos <- start(ir) ##start end the same
            dist2 <- dist2*ifelse(
                (Strand==1 & 
                 end(exons[map2,]) < pos ) |
                (Strand==-1 &
                 start(exons[map2,]) > pos),-1,1)
            
            tmp<-cbind(dist2,map2)
            m1 = tmp[1,1]
            m2 = tmp[2,1]
            exon1 = tmp[1,2]
            exon2 = tmp[2,2]
            m1m2Index=which.min(abs(c(m1,m2)))
            
            if(exon1==exon2 & m1==0 & m2==0){
                subtype[j]="inside exon"
                insidedistance[j]=0
                exonnumber[j]=c(exon1,exon2)[m1m2Index]
            }
            
            if( (sign(m1)==sign(m2) & (m1!=0 & m2!=0) & exon1==exon2) |
               (sign(m1)==-1 & sign(m2)==1 & exon2-exon1==1) ){
                subtype[j]="inside intron"
                
                insidedistance[j]=c(m1,m2)[m1m2Index]
                exonnumber[j]=c(exon1,exon2)[m1m2Index]
                
            }
            if( (exon2-exon1 > 1) |
               ( (exon2-exon1 == 1) & 
                ((sign(m1)==-1 & sign(m2)==-1) |
                 (sign(m1)==1 & sign(m2)==-1) |
                 (sign(m1)==1 & sign(m2)==1) |
                 (sign(m1)==0 & sign(m2)==-1)|
                 (sign(m1)==1 & sign(m2)==0))) |
               (exon2==exon1 & sign(m1)!=sign(m2))){
                subtype[j]="covers exon(s)"
                insidedistance[j]=0
                exonnumber[j]=c(exon1,exon2)[m1m2Index]
            }
            
            if( (exon2-exon1 <= 1 & sign(m1)==-1 & sign(m2)==0) |
               (exon2==exon1 & sign(m1)==1 & sign(m2)==0)){
                if(Strand==1) subtype[j]="overlaps exon upstream" else subtype[j]="overlaps exon downstream"
                insidedistance[j]=0
                exonnumber[j]=c(exon1,exon2)[m1m2Index]
            }
            if( (exon2-exon1 <= 1 & sign(m1)==0 & sign(m2)==1) |
               (exon2==exon1 & sign(m1)==0 & sign(m2)==-1)){
                if(Strand==-1) subtype[j]="overlaps exon upstream" else subtype[j]="overlaps exon downstream"
                insidedistance[j]=0
                exonnumber[j]=c(exon1,exon2)[m1m2Index]
            }
            if( exon2-exon1 == 1 & sign(m1)==0 & sign(m2)==0){
                subtype[j]="overlaps two exons"
                insidedistance[j]=0
                exonnumber[j]=c(exon1,exon2)[m1m2Index]
            }
            
            if(Strand!=1){
                insidedistance[j]= -insidedistance[j]
                exonnumber[j] = nrow(Exons) - exonnumber[j] + 1
            }
            
            ctype[j]="inside transcription region"
            
            if(!is.na(CS)) {
                if(S<CS & E<CS){
                    if(Strand==1) ctype[j]="5' UTR" else ctype[j]="3'UTR"
                }
            }
            
            
            if(!is.na(CE)) {
                if(S>CE & E>CE){
                    if(Strand==-1) ctype[j]="5' UTR" else ctype[j]="3'UTR"
                }
            }
            
            if(!is.na(CS) & !is.na(CE)) {
                if(S<CS & E>CE){
                    ctype[j]="covers coding region"
                }
            }
            
            if(!is.na(CS)) {
                if(S<CS & E>CS){
                    if(Strand==1) ctype[j]="overlaps 5' UTR" else ctype[j]="overlaps 3'UTR"
                }
            }
            
            if(!is.na(CE)) {
                if(S<CE & E>CE){
                    if(Strand==-1) ctype[j]="overlaps 5' UTR" else ctype[j]="overlaps 3'UTR"
                }
            }
            
            
        }
    ##    if(TE-TS<10^5){##graphical check
            
    ##         plot(0,0,ylim=c(0,0.6),xlim=range(c(start(subject[i,]),end(subject[i,]),start(x[j,]),end(x[j,]))),
    ##              xlab=paste("inside distance=",dist[j],insidedistance[j],m1,m2,exonnumber[j]))
    ##         polygon(c(TS,TE,TE,TS),c(0,0,0.5,0.5),density=0,col=2)
    ##         polygon(c(CS,CE,CE,CS),c(0.1,0.1,0.4,0.4),density=0,col=3)
    ##         abline(h=0.25,lwd=2)
    ##         apply(Exons,1,function(x) polygon(c(x[1],x[2],x[2],x[1]),
    ##                                           c(0.2,0.2,0.3,0.3),col=4))
            
    ##         polygon(c(start(x[j,]),end(x[j,]),end(x[j,]),start(x[j,])),c(0.4,0.4,0.5,0.5),col=5)
    ##         lines(c(TS,TS+1000),c(0.55,0.55),lwd=3)
    ##         title(paste(j,i,Strand,type[j],subtype[j],ctype[j],dist[j],sep=":"))
    ##     }

    }
    
    type[dist<=promoterDist & type=="upstream"] <- "promoter"
    type[dist<=promoterDist & type=="downstream"] <- "close to 3'"
    
    description=type
    tmpIndex=which(description=="inside")
    description[tmpIndex] <- subtype[tmpIndex]
    tmp <- data.frame(name=I(genenames),
                      annotation=I(geneannotation),
                      description=factor(description,levels=c("upstream","promoter","overlaps 5'","inside intron","inside exon","covers exon(s)","overlaps exon upstream","overlaps exon downstream","overlaps two exons","overlaps 3'","close to 3'","downstream","covers")),
                      region=factor(type,levels=c("upstream","promoter","overlaps 5'","inside","overlaps 3'","close to 3'","downstream","covers")),
                      distance=dist,
                      subregion=factor(subtype,levels=c("inside intron","inside exon","covers exon(s)","overlaps exon upstream","overlaps exon downstream","overlaps two exons")),
                      insideDistance=insidedistance,
                      exonnumber=exonnumber,
                      nexons=nexons,
                      UTR=factor(ctype,levels=c("inside transcription region","5' UTR","overlaps 5' UTR","3'UTR","overlaps 3'UTR","covers transcription region")),
                      strand=strands,
                      geneL=geneL,
                      codingL=codingL,
                      Geneid=geneid,
                      subjectHits=map)
    return(tmp)
}
