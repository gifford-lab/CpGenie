hashtable <- c('A','T','G','C','N')
hashtable_rev <- c('T','A','C','G','N')

translate <- function(seq){
	return(paste0(hashtable[seq+1],collapse=''))
}

translateRev <- function(seq){
	return(paste0(hashtable_rev[rev(seq+1)],collapse=''))
}

untranslate<- function(seq){
	return(sapply(1:length(seq),function(idx){which(hashtable==seq[idx])-1}))
}

readGenome<-function(genomefile,s,len){
	con <- file(genomefile,open='rb')
	seek(con,where = s)
	rbequiv=readBin(con,integer(),size=1,n=len)
	close(con)
	rbequiv
}
readRef<-function(refile,offsetfile,chr){
	offsets = as.double(readLines(offsetfile))
    chrlen=diff(offsets)
	readsize=chrlen[chr]
	con <- file(refile,open='rb')
	seek(con,where = offsets[chr])
	rbequiv=readBin(con,integer(),size=1,n=readsize)
	close(con)
	rbequiv
}

pullKmer <- function(chridx){
	thiswin = win[[chridx]]
	ref <- readRef(genomefile,offsetfile,chrs[chridx])

	return(lapply(1:length(thiswin),function(reg){
		left = min(thiswin[reg])
		right = max(thiswin[reg])
		sapply(left:(right-klen+1),function(idx){
			translate(ref[idx:(idx+klen-1)])
		})
	}))
}

pullKmerTwoDirec <- function(chridx){
	thiswin = win[[chridx]]
	ref <- readRef(genomefile,offsetfile,chrs[chridx])

	return(lapply(1:length(thiswin),function(reg){
		left = min(thiswin[reg])
		right = max(thiswin[reg])
		c( sapply(left:(right-klen+1),function(idx){
			translate(ref[idx:(idx+klen-1)])
		}),sapply((right-klen+1):left,function(idx){
		            translateRev(ref[idx:(idx+klen-1)])
		        }))
	}))
}

pullKmerVCF <- function(chridx){
	thiswin = win[[chridx]]
	thisalt = alt[[chridx]]
	ref <- readRef(genomefile,offsetfile,chrs[chridx])
	
	return(lapply(1:length(thiswin),function(reg){
		left = min(thiswin[reg])
		right = max(thiswin[reg])
		ref_part = ref[left:right]
		ref_part[(length(ref_part)/2+1)] = untranslate(thisalt[reg])
		sapply(1:(right-left + 1-klen+1),function(idx){
			translate(ref_part[idx:(idx+klen-1)])
		})
	}))
}

pullKmerVCFtwodirec <- function(chridx){
	thiswin = win[[chridx]]
	thisalt = alt[[chridx]]
	ref <- readRef(genomefile,offsetfile,chrs[chridx])
	
	return(lapply(1:length(thiswin),function(reg){
		left = min(thiswin[reg])
		right = max(thiswin[reg])
		ref_part = ref[left:right]
		ref_part[(length(ref_part)/2+1)] = untranslate(thisalt[reg])
		c(sapply(1:(right-left + 1-klen+1),function(idx){
			translate(ref_part[idx:(idx+klen-1)])
		}),sapply((right-left+1-klen+1):1,function(idx){
				translateRev(ref_part[idx:(idx+klen-1)])})
		)
	}))
}

pullKmerNonOverlap <- function(chridx){
	thiswin = win[[chridx]]
	ref <- readRef(genomefile,offsetfile,chrs[chridx])

	out =list()

	for (reg in 1:length(thiswin)){
		left = min(thiswin[reg])
		right = max(thiswin[reg])
		new = c()
		for (startpos in left:(left+klen-1)){
			kmernum = floor((right-startpos+1)/klen)
			new = c(new,sapply(1:kmernum, function(idx){
				left = startpos + (idx-1)*klen
				right = left + klen - 1
				translate(ref[left:right])
			}))
		}
		out = c(out,list(new))
	}
	return(out)
}

pullSeq <- function(loci,geomefile,offsetfile){
	if (nrow(loci)==0) return(c())
	uni_chr = unique(loci[,1])
	seq = matrix(nrow=nrow(loci),ncol=1)
	for (chridx in 1:length(uni_chr)){
		chr = uni_chr[chridx]
		print(paste0('chr',chr))
		ref <- readRef(genomefile,offsetfile,chr)
		pick = which(loci[,1]==chr)
		part = loci[pick,]
		for (row in 1:nrow(part)){
			left = part[row,2]
			right = part[row,3]
			seq[pick[row],1] = translate(ref[left:right])
		}
	}
	return(seq)
}

pullKmerNonOverlapTwoDirec <- function(chridx){
	thiswin = win[[chridx]]
	ref <- readRef(genomefile,offsetfile,chrs[chridx])

	out =list()

	for (reg in 1:length(thiswin)){
		left = min(thiswin[reg])
		right = max(thiswin[reg])
		new = c()
		for (startpos in left:(left+klen-1)){
			kmernum = floor((right-startpos+1)/klen)
			new = c(new,sapply(1:kmernum, function(idx){
				left = startpos + (idx-1)*klen
				right = left + klen - 1
				translate(ref[left:right])
			}))
		}
		for (endpos in right:(right-klen+1)){
        	kmernum = floor((endpos-left+1)/klen)
        	new = c(new,sapply(1:kmernum, function(idx){
        		right = endpos - (idx-1)*klen
        		left = right-klen + 1
        		translateRev(ref[left:right])
        	}))
        }
		out = c(out,list(new))
	}
	return(out)
}
pullSeqRev <- function(loci,geomefile,offsetfile){
	if (nrow(loci)==0) return(c())
	uni_chr = unique(loci[,1])
	seq = matrix(nrow=nrow(loci),ncol=1)
	for (chridx in 1:length(uni_chr)){
		chr = uni_chr[chridx]
		print(paste0('chr',chr))
		ref <- readRef(genomefile,offsetfile,chr)
		pick = which(loci[,1]==chr)
		part = loci[pick,]
		for (row in 1:nrow(part)){
			left = part[row,2]
			right = part[row,3]
			seq[pick[row],1] = translateRev(ref[left:right])
		}
	}
	return(seq)
}
