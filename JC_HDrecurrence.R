args=commandArgs(TRUE)
run<-as.integer(args[1])

set.seed(run)

setwd("/lustre/scratch110/sanger/dw9/Jiqiu_Mar2014")

homdel.file= "HomoD_twoevents.txt"
hemidel.file = "HemiD_oneevent_refine.txt"

homdel.data = read.table(homdel.file,sep="\t",header=T,stringsAsFactors=F)
hemidel.data = read.table(hemidel.file,sep="\t",header=T,stringsAsFactors=F)

mixtureModelInfoMajor = read.table("majorDel_MM_info.txt",sep="\t",header=T)
mixtureModelInfoMinor = read.table("minorDel_MM_info.txt",sep="\t",header=T)

#just use counts for specific tumour type
segment.counts = read.table("all_segments_data.txt",sep="\t",header=T,stringsAsFactors=F)

chr.names=c(1:22,"X")
chr.data=list(0)
chr.seg.counts = array(0,23)
for(chr in 1:23){
	chr.data[[chr]] = segment.counts[segment.counts$chr==chr.names[chr],]
	chr.seg.counts[chr] = nrow(chr.data[[chr]])
}
cum.chr.seg.counts = c(0,cumsum(chr.seg.counts))

chr.lengths = read.table("chr_lengths.txt",sep="\t",header=T)[,2]
genome.length = sum(as.numeric(chr.lengths))
#half.genome.length=round(genome.length/2)
chr.boundaries = c(0,cumsum(as.numeric(chr.lengths)))

sim.counts = vector(mode="numeric",length = nrow(segment.counts))

homdel.data$HD.length = homdel.data$HomoDend - homdel.data$HomoDstart+1

#no.perms = 10000
no.perms = 2000
#no.perms=2
for(p in 1:no.perms){
	#if(p %% 10 == 0){
		print(paste("perm ",p,sep=""))
	#}
	hemidel.data$log.length = log(hemidel.data$length*1000000)
	hemidel.data$fixed = sapply(1:nrow(hemidel.data),function(h,m,n,i){probs=c(sum(dnorm(h$log.length[i],m$mean,m$sd)*m$prop),dnorm(h$log.length[i],n$mean,n$sd));prob = probs[1]/sum(probs);runif(1)<prob},h=hemidel.data, m=mixtureModelInfoMajor, n=mixtureModelInfoMinor)
	sim.HDs = array(NA,c(0,4))
	for(sample in unique(homdel.data$sample)){
		#print(sample)
		sample.sim.HDs = NULL
		sample.data.hom = homdel.data[homdel.data$sample==sample,]
		sample.data.hemi = hemidel.data[hemidel.data$sample==sample,]

		sample.data.fixed = rbind(sample.data.hom[,3:5],sample.data.hemi[sample.data.hemi$fixed,3:5])
		temp = sample.data.hemi[!sample.data.hemi$fixed,3:5]
		names(temp) = names(sample.data.hom)[c(3,6,7)]
		sample.data.unfixed = rbind(sample.data.hom[,c(3,6,7)],temp)
		sample.data.unfixed$HD.length = sample.data.unfixed[,3] - sample.data.unfixed[,2]

		sampled.dels = array(NA,c(0,3))
		#print(paste(sample,nrow(sample.data),sep=","))
		for(i in 1:nrow(sample.data.unfixed)){
			HD.OK = F
			while(!HD.OK){
				#check that the homdel doesn't cross over 2 chromosomes
				chr=-2
				chr.end=-1
				while(chr!=chr.end)
				{
					#genome.start.pos = sample(genome.length-sample.data.unfixed$HD.length[i]+1,1)
					#genome.length is too large
					half.genome.length=round((genome.length-sample.data.unfixed$HD.length[i])/2)
					genome.start.pos = sample(half.genome.length,1)
					if(sample(2,1)==2){
						genome.start.pos = genome.start.pos + half.genome.length
					}
					genome.end.pos = genome.start.pos+sample.data.unfixed$HD.length[i]-1
					chr = sum(genome.start.pos>chr.boundaries)
					chr.end = sum(genome.end.pos>chr.boundaries)
				}
				#print(paste("chr=",chr,sep=""))
				chr.start.pos = genome.start.pos - chr.boundaries[chr]
				chr.end.pos = genome.end.pos - chr.boundaries[chr]
				#gene.index1 = findInterval(chr.start.pos,chr.data[[chr]]$Startpos)
				HD.OK=T
			}
			sampled.dels = rbind(sampled.dels,c(chr,chr.start.pos,chr.end.pos))
			#sample.sim.HDs = rbind(sample.sim.HDs,c(chr,chr.start.pos,chr.end.pos))
		}
		colnames(sampled.dels) = names(sample.data.fixed)
		all.dels = rbind(sample.data.fixed,sampled.dels)

		sample.sim.HDs = array(NA,c(0,4))
		for(chr in unique(all.dels$chr)){
			c = which(chr==chr.names)
			chr.del.data = all.dels[all.dels$chr==chr,2:3]
			names(chr.del.data) = c("startpos","endpos")
			chr.del.data$endpos = chr.del.data$endpos +1
			chr.del.data =chr.del.data[order(chr.del.data$startpos),]
			all.breakpoints = c(1,sort(unique(unlist(chr.del.data))),chr.lengths[c])
			no.breakpoints = length(all.breakpoints)-1
			chr.sample.sim.HDs = cbind(chr,all.breakpoints[1:no.breakpoints],all.breakpoints[-1],0)
			chr.sample.sim.HDs[,4] = sapply(1:nrow(chr.sample.sim.HDs),function(d,s,i){sum(d$startpos >= as.numeric(s[i,2]) & d$startpos < as.numeric(s[i,3])| d$endpos > as.numeric(s[i,2]) & d$endpos <= as.numeric(s[i,3]) | d$startpos < as.numeric(s[i,2]) & d$endpos > as.numeric(s[i,3]))},d=chr.del.data,s=chr.sample.sim.HDs)
			chr.sample.sim.HDs[1:(no.breakpoints-1),3] = as.numeric(chr.sample.sim.HDs[1:(no.breakpoints-1),3])-1
			inds = which(chr.sample.sim.HDs[,4]>1)
			if(length(inds)>1){
				sample.sim.HDs = rbind(sample.sim.HDs,cbind(sample,chr.sample.sim.HDs[inds,1:3]))
			}else if(length(inds)==1){
				sample.sim.HDs = rbind(sample.sim.HDs,c(sample,chr.sample.sim.HDs[inds,1:3]))
			}
		}
		sim.HDs = rbind(sim.HDs,sample.sim.HDs)
	}
	sim.HDs = sim.HDs[,-1]
	#segment.data = array(NA,c(0,4))
	sim.hits = vector(mode="numeric",length = nrow(segment.counts))
	for(c in 1:23){
		chr = chr.names[c]
		#print(chr)
		chr.sim.data = data.frame(startpos = as.numeric(sim.HDs[sim.HDs[,1]==c,2]),endpos = as.numeric(sim.HDs[sim.HDs[,1]==c,3]))
		all.breakpoints = c(1,sort(unique(unlist(chr.sim.data))),chr.lengths[c])
		no.breakpoints = length(all.breakpoints)-1
		chr.segment.data = cbind(chr,all.breakpoints[1:no.breakpoints],all.breakpoints[-1],0)
		chr.segment.data[,4] = sapply(1:nrow(chr.segment.data),function(d,s,i){sum(d$startpos >= as.numeric(s[i,2]) & d$startpos < as.numeric(s[i,3])| d$endpos > as.numeric(s[i,2]) & d$endpos <= as.numeric(s[i,3]) | d$startpos < as.numeric(s[i,2]) & d$endpos > as.numeric(s[i,3]))},d=chr.sim.data,s=chr.segment.data)
		chr.segment.data[1:(no.breakpoints-1),3] = as.numeric(chr.segment.data[1:(no.breakpoints-1),3])-1
		colnames(chr.segment.data)=c("chr","startpos","endpos","count")
		chr.segment.data = data.frame(chr=chr.segment.data[,1],startpos=as.numeric(chr.segment.data[,2]),endpos=as.numeric(chr.segment.data[,3]),count=as.numeric(chr.segment.data[,4]))
		#segment.data = rbind(segment.data,chr.segment.data)
		no.segments = chr.seg.counts[c]
		for(s in 1:no.segments){
			sim.hits[cum.chr.seg.counts[c]+s] = sum(sapply(1:nrow(chr.segment.data),function(c2,d,l,i){max(min(c2$endpos[s],d$endpos[i])-max(c2$startpos[s],d$startpos[i])+1,0)*d$count[i]/l},d=chr.segment.data,c2=chr.data[[c]],l=chr.data[[c]]$endpos[s] - chr.data[[c]]$startpos[s]+1))
			if(sim.hits[cum.chr.seg.counts[c]+s]>=chr.data[[c]]$count[s]){
				sim.counts[cum.chr.seg.counts[c]+s] = sim.counts[cum.chr.seg.counts[c]+s] + 1
			}
		}
	}
}
write.table(sim.counts,paste("simulations_bySegment/sim_counts_all_samples_",run,".txt",sep=""),sep="\t",quote=F)
