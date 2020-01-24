##single sims
mu=1.5e-8
recomb_rate=1e-8
Ne=10000
nBases=1e6
samplesize=200
s=0.1
fix=1
discoal_path="~/work/programs/discoal/discoal"

temp<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,sweep="hard")
G<-temp$genomes

library(jpeg)
writeJPEG(G,"~/Desktop/test.jpeg",quality=1)
  
  
  
  
