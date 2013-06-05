graph.combine = function(g1,layout1,g2,layout2){
 require(igraph);
 g1 = delete.edges(g1,E(g1))
 g = graph.disjoint.union(g1, g2)
 V(g)$color = cbind(V(g1)$color,V(g2)$color);
 V(g)$label = cbind(V(g1)$label,V(g2)$label);
 V(g)$size = cbind(V(g1)$size,V(g2)$size);
 V(g)$size2 = cbind(V(g1)$size2,V(g2)$size2);
 glayout = rbind(layout1,layout2);
 list(g,glayout)
}

printsol =  function(instance1,tmax,lmax,tour,intensity1 = c(),type = "circle",ctourfactor=1)   {
# printsol(instance1,tmax,lmax,result$tour,result$intensity,type = "circle",ctourfactor=1)
	
	costsRange = range(instance1[instance1[,3]>0,4])
	profitsRange = range(instance1[instance1[,3]<0,4])
	if(costsRange[1]==costsRange[2] && profitsRange[1]==profitsRange[2]){
		 type = "circle";
		 print(type)
	}
	
	tcost = matrix(0,dim(instance1)[1],dim(instance1)[1]);
	mindist =  1e10;
	for( i in 1:dim(instance1)[1]){
	    for( j in 1:dim(instance1)[1]){
		   tcost[i,j]=ceiling(10000*sqrt((instance1[i,1]-instance1[j,1])^2+(instance1[i,2]-instance1[j,2])^2))/10000;
		   if(tcost[i,j] > 0){
				mindist = min(mindist,tcost[i,j]);
		   }
		}
	}
	
	
	require(igraph);
	pmmfrow =  par()$mfrow
	pmargins =  par()$mar

	balance = instance1[,3];
	intensity = c();
	if (length(intensity1) > 0){
	   intensity = 1 + instance1[,1] - instance1[,1];
	   intensity[tour] = intensity1;
	}
	
	g1 = graph.empty();
	g1 = add.vertices(g1, dim(instance1)[1] );
	V(g1)$color[instance1[,3]<0] = "red";
	V(g1)$color[instance1[,3]==0] = "yellow";
	V(g1)$color[instance1[,3]>0] = "grey";
	V(g1)$label = 1:length(instance1[,1]);
	
	xl = max(instance1[,1])-min(instance1[,1])
	yl = max(instance1[,2])-min(instance1[,2])
	xyl = min(xl,yl)
	if(type == "circle"){
		V(g1)$size = 5*xyl*sqrt(abs(instance1[,3]))/sqrt(median(abs(instance1[,3])));
		V(g1)$size[c(1,2)] = 5*xyl;
		V(g1)$size2 = 5*xyl*abs(instance1[,4])/median(abs(instance1[,4]));
	} else if (type == "rectangle"){
	    V(g1)$size = 15*xyl*(mindist/2)*abs(instance1[,3])/max(abs(instance1[,3]));
		V(g1)$size[c(1,2)] =  5*xyl;
		V(g1)$size2 = 15*xyl*(mindist/2)*abs(instance1[,4])/max(abs(instance1[,4]));
		V(g1)$size2[V(g1)$size2==0]=1;
		V(g1)$size2[c(1,2)] = 5*xyl;
		
		#V(g1)$size = 5*xyl*abs(instance1[,3])/median(abs(instance1[,3]));
		#V(g1)$size[c(1,2)] =  5*xyl;
		#V(g1)$size2 = 5*xyl*abs(instance1[,4])/median(abs(instance1[,4]));
		#V(g1)$size2[V(g1)$size2==0]=1
		#V(g1)$size2[c(1,2)] = 5*xyl;
	}
	
	touredges = cbind(tour[-length(tour)],tour[-1])
	g1 = add.edges(g1,as.vector(t(touredges)));
	glayout = as.matrix(instance1[,c(1,2)])
	colnames(glayout) = NULL
	#tkid = tkplot(g1,layout = glayout);
	par(mfrow = c(1,2));
	par(mar = pmargins + c(0,0,0,2));
  
	if (length(intensity) == 0){
	  loading = cumsum(instance1[tour,3]); # * tmax / lmax;
	  profits = (instance1[,3]*instance1[,4])[tour];
	} else
	{
	  loading = cumsum(instance1[tour,3]*intensity[tour]); # * tmax / lmax;
	  profits = (instance1[,3]*instance1[,4])[tour]*intensity[tour];
	}
	
  
	profitscum = cumsum(profits);
	
	tlength = 0; 
	tl=c(0);
	for(i in 1:(length(tour)-1)){
	   tlength=tlength+tcost[tour[i],tour[i+1]]; 
	   tl[i+1] = tlength;
	};
	profit = sum(profits) - ctourfactor*tlength;
	ymax = max(max(loading),tmax,lmax);
	plot(loading,xaxt = "n",col="blue",type="o",ylim=c(min(-10,min(loading)),ymax));
	axis(1, at=1:length(tour), labels=tour, cex.axis=0.4)
	
	
	lines(tl, col="blue",lty=2);
	text(length(profitscum)-1,tlength,paste("tlength:",tlength),col="blue")
	
	#lines(0.01*profitscum,col="blue",lty=3)
	
	
	lines(balance[tour], col="blue",lty=4);
	if (length(intensity) > 0){
	   lines(balance[tour]*intensity[tour], col="blue",lty=4);
	}
	
	text(length(profitscum)-1,tmax+ymax/50,"tmax")
	
	lines(matrix(tmax,length(profitscum),1));
	
	 
	text(2,lmax-ymax/50,"lmax")
	lines(matrix(lmax,length(profitscum),1));
	text(which.max(loading)+1,max(loading),paste("maxload:",max(loading)),col="blue")
	
	par(new=TRUE)
	plot(profitscum,type="l",col="red",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(floor(1.1*min(profitscum)),ceiling(1.1*max(profitscum))))
	axis(4)
	mtext("profits",side=4,line=3)
	text(length(profitscum)-1,profit,paste("profit:",profit),col="red")
	
	
	ldist = xyl/2
	
	par(mar = pmargins + c(0,2,0,0));
	if(type == "circle"){
		if (length(intensity) == 0){
			 tttt = g1;
			 plot(tttt,axes=TRUE,rescale=F,xlab = 'x',ylab = 'y', asp = 1,xlim = range(instance1[,1]),ylim = range(instance1[,2]),layout = glayout,vertex.label.family="serif",edge.label.family="Palatino",edge.arrow.size = 0.7, edge.width = 1.5, edge.color = "black",edge.arrow.width = 0.7);
		} else
		{
		  g2 = g1;
		  print(" ...  ")
		  V(g2)$size = intensity * V(g1)$size;
		  
		  tttt = graph.combine(g1,glayout,g2,glayout);
		  print(" ...  ")
		  plot(tttt[[1]],axes=TRUE,rescale=F,xlab = 'x',ylab = 'y', asp = 1,xlim = range(instance1[,1]),ylim = range(instance1[,2]),vertex.label.dist = ldist,vertex.label.degree=pi/4,layout = tttt[[2]],vertex.label.family="serif",edge.label.family="Palatino",edge.arrow.size = 0.7, edge.width = 1.5, edge.color = "black",edge.arrow.width = 0.7)		  		  
		}
	} else if (type == "rectangle"){
		if (length(intensity) == 0){
			 tttt = g1;
			 zzzz = V(tttt)$size;
			 V(tttt)$size = V(tttt)$size2;
			 V(tttt)$size2 = zzzz;
			 plot(tttt,axes=TRUE,layout = glayout,rescale=F, xlab = 'x',ylab = 'y',asp = 1,xlim = range(instance1[,1]),ylim = range(instance1[,2]),vertex.label.family="serif",edge.label.family="Palatino",vertex.shape="rectangle",edge.arrow.size = 0.7, edge.width = 1.5, edge.color = "black",edge.arrow.width = 0.7)

		} else
		{
		  print("?????????")
		  g2 = g1;
		  V(g2)$size = sqrt(intensity) * V(g1)$size;
		  V(g2)$size2 = sqrt(intensity) * V(g1)$size2;
		  
		  tttt = graph.combine(g1,glayout,g2,glayout);
		  print("?????????")
		  zzzz = V(tttt[[1]])$size;
		  V(tttt[[1]])$size = V(tttt[[1]])$size2;
		  V(tttt[[1]])$size2 = zzzz;
		  
		  print("?????????")
		  #print( V(tttt[[1]]))
		  #print( V(tttt[[1]])$size)
		  #print( V(tttt[[1]])$size2)
		  #print(tttt[[2]])		  
		  plot(tttt[[1]],axes=TRUE,rescale=F,xlab = 'x',ylab = 'y', asp = 1,xlim = range(instance1[,1]),ylim = range(instance1[,2]),vertex.label.dist = ldist,vertex.label.degree=pi/4,layout = tttt[[2]],vertex.label.family="serif",edge.label.family="Palatino",vertex.shape="rectangle",edge.arrow.size = 0.7, edge.width = 1.5, edge.color = "black",edge.arrow.width = 0.7)
		  print("###########")
		}
		#rect(1, 1, 3, 2, col="red", border="black")
		#rect(1.2, 1.1, 2.8, 1.9, col="red", border="black")
		#text(2,1.5,'profit')
		#par(srt = 90)
		#text(0.5,1.5,'profit/unit',cex = 0.8)
		#par(srt = 0)
		#text(2,2.4,'max qty',cex = 0.8)

	}
	par(mfrow = pmmfrow);  
	par(mar = pmargins); 
	list(g = tttt[[1]], glayout = glayout)
}