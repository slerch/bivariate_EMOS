
correlation_fun <- function(x1, x2, y1, y2, radius, ...){
	sect.2 <- sect.3 <- sect.4 <- sect.5 <- sect.6 <- sect.7 <- sect.8 <- sect.9 <- NULL
	
	# define sections according to wind direction and sort
	for (i in 1:length(x1)){
	  if (x1[i]^2 + y1[i]^2 <= radius^2) next
	  else phi <- atan2(-x1[i], -y1[i])/(2*pi)*360
	  if (phi < 0) phi <- 360 + phi
	  if ((phi >= 0) & (phi < 45)) sect.6 <- rbind(sect.6, cbind(x2[i], y2[i]))
	  if ((phi >= 45) & (phi < 90)) sect.7 <- rbind(sect.7, cbind(x2[i], y2[i]))
	  if ((phi >= 90) & (phi < 135)) sect.8 <- rbind(sect.8, cbind(x2[i], y2[i]))
	  if ((phi >= 135) & (phi < 180)) sect.9 <- rbind(sect.9, cbind(x2[i], y2[i]))
	  if ((phi >= 180) & (phi < 225)) sect.2 <- rbind(sect.2, cbind(x2[i], y2[i]))
	  if ((phi >= 225) & (phi < 270)) sect.3 <- rbind(sect.3, cbind(x2[i], y2[i]))
	  if ((phi >= 270) & (phi < 315)) sect.4 <- rbind(sect.4, cbind(x2[i], y2[i]))
	  if ((phi >= 315) & (phi < 360)) sect.5 <- rbind(sect.5, cbind(x2[i], y2[i]))
	}
	
	# set corr if too few or equal data 
	if (length(sect.2[,1]) %in% c(0,1)) {
	  cor.2 <- 0 
	} else if (sum(var(sect.2)==0)) {
	  cor.2 <- 1 
	} else 	{
	  cor.2 <- cor(sect.2[,1], sect.2[,2])
	}
	if (length(sect.3[,1]) %in% c(0,1)) {
	  cor.3 <- 0 
	} else if (sum(var(sect.3)==0)) {
	  cor.3 <- 1 
	} else 	{
	  cor.3 <- cor(sect.3[,1], sect.3[,2])
	}
	if (length(sect.4[,1]) %in% c(0,1)) {
	  cor.4 <- 0 
	} else if (sum(var(sect.4)==0)) {
	  cor.4 <- 1 
	} else 	{
	  cor.4 <- cor(sect.4[,1], sect.4[,2])
	}
	if (length(sect.5[,1]) %in% c(0,1)) {
	  cor.5 <- 0 
	} else if (sum(var(sect.5)==0)) {
	  cor.5 <- 1 
	} else 	{
	  cor.5 <- cor(sect.5[,1], sect.5[,2])
	}
	if (length(sect.6[,1]) %in% c(0,1)) {
	  cor.6 <- 0 
	} else if (sum(var(sect.6)==0)) {
	  cor.6 <- 1 
	} else 	{
	  cor.6 <- cor(sect.6[,1], sect.6[,2])
	}
	if (length(sect.7[,1]) %in% c(0,1)) {
	  cor.7 <- 0 
	} else if (sum(var(sect.7)==0)) {
	  cor.7 <- 1 
	} else 	{
	  cor.7 <- cor(sect.7[,1], sect.7[,2])
	}
	if (length(sect.8[,1]) %in% c(0,1)) {
	  cor.8 <- 0 
	} else if (sum(var(sect.8)==0)) {
	  cor.8 <- 1 
	} else 	{
	  cor.8 <- cor(sect.8[,1], sect.8[,2])
	}
	if (length(sect.9[,1]) %in% c(0,1)) {
	  cor.9 <- 0 
	} else if (sum(var(sect.9)==0)) {
	  cor.9 <- 1 
	} else 	{
	  cor.9 <- cor(sect.9[,1], sect.9[,2])
	}
	
	cor <- c(cor.6, cor.7, cor.8, cor.9, cor.2, cor.3, cor.4, cor.5)
	l <- c(length(sect.6[,1]), length(sect.7[,1]), length(sect.8[,1]), length(sect.9[,1]), length(sect.2[,1]), length(sect.3[,1]), length(sect.4[,1]), length(sect.5[,1]))
	y <- rbind(cor,l)
	
	# wind direction angle
	theta <- c(22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5)
	
	plot(theta, y[1,], type="b", xlab="Wind Direction in Degrees", ylab="Correlation", xlim=c(0,360), ylim=c(-1,1), ...)
	text(x = theta, y = y[1,], labels = y[2,], font=2, cex = 1.5, pos=2)
	return(y)
}



