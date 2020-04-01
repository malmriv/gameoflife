#Read the data
data = read.table("results.txt")
s = data$V1

#Set size of the lattice and MC steps (check in code)
N = 128
steps = 500

#Create a position vector
x = c(1:N)
y = x

#Set (plotting) size of each node and number of frames to plot
nodesize = 0.6
layers = 500

#Create a directory to save the frames
dir.create("./frames")

for(i in 1:steps) {
  if(i%%round(steps/layers,0) == 0) {
  #Create a new .PNG file
  png(paste("./frames/",i,".png",sep=""),height=500,width=500)
  
  #Create a plot and set desired graphics
  par(mar = c(0, 0, 0, 0),col.lab="black",col.axis="black",cex.main=1.6,bg="black",xpd=NA)
  plot(0,0,type="n",xlim=c(0,N+1),ylim=c(0,N+1),asp=1)
  box(col="black")
  
  for(j in 1:N) {
    for(k in 1:N) {
      if(s[(i-1)*N^2+(j-1)*N+k] == 1) lines(x[j],y[k],type="p",pch=15,col="#fcf8e4",cex=nodesize)
      if(s[(i-1)*N^2+(j-1)*N+k] == 0) {}
    }
  }
  dev.off()
  }
}