setwd('C:\\Users\\Martin\\Desktop\\projects\\pickupanddeliveryteamorienteering\\R\\')
source('printsol.r')
source('example1.txt')
source('example2.txt')
pdf('test123.pdf', width = 25, height = 10) #paper = "a4r",
printsol(instance1,tmax,lmax,result$tour,result$intensity,type = "circle",ctourfactor=1)
dev.off()
pdf('test123.pdf', width = 25, height = 10) #paper = "a4r",
printsol(instance1,tmax,lmax,result$tour,result$intensity,type = "rectangle",ctourfactor=1)
dev.off()
save.image()
#quit(save = 'no')
