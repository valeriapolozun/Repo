setwd('C:\\Users\\User\\Documents\\Visual Studio 2012\\Projects\\Repo\\ConsoleApplication1\\Rfiles\\')
source('printsol.r')
source('example1.txt')
source('example2.txt')
source('example3.txt')
pdf('test123.pdf', width = 25, height = 10) #paper = "a4r",
printsol(instance1,tmax,lmax,result$tour,result$intensity,type = "circle",ctourfactor=1)
dev.off()
pdf('test123.pdf', width = 25, height = 10) #paper = "a4r",
printsol(instance1,tmax,lmax,result$tour,result$intensity,type = "rectangle",ctourfactor=1)
dev.off()
save.image()
#quit(save = 'no')
