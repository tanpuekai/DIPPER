
library(grid)
library(lattice)
library(gridExtra)
library(stringr)
##########################################
library(e1071)
SVMPLOT<-function(PCAi,Yin,percenti,DIM="12"){

	Y<-Yin
	TEXT<-""
	if(length(unique(Y))<2)return()
	if(DIM=="12"){
		inddim<-c(1,2)
	}else{
		inddim<-c(1,3)
	}
	print(length(unique(Y)))
	print(table(Y))
	print(dim(PCAi))
	dat = data.frame(PCAi$x[,inddim], y = Y)
	svmfit = svm(y ~ ., data = dat, kernel = "polynomial",coef0	=1)
	#plot(svmfit, dat)
	px1<-seq(range(PCAi$x[,inddim[1]])[1]*1.05,range(PCAi$x[,inddim[1]])[2]*1.05,len=20)
	px2<-seq(range(PCAi$x[,inddim[2]])[1]*1.05,range(PCAi$x[,inddim[2]])[2]*1.05,len=20)

	if(DIM=="12"){
		xgrid = expand.grid(PC1 = px1, PC2= px2)
	}else{
		xgrid = expand.grid(PC1 = px1, PC3= px2)
	}

	if(1){

		ygrid = predict(svmfit, xgrid)

		func = predict(svmfit, xgrid, decision.values = TRUE)
		func = attributes(func)$decision

		plot(PCAi$x[,inddim],cex=2,xlab=percenti[1],ylab=percenti[2],
			main=TEXT,
			pch=PCHanita,
			col=KOLanita)
		contour(px1, px2, matrix(func, 20, 20), level = 0, add = TRUE)
		contour(px1, px2, matrix(func, 20, 20), level = 0.5, add = TRUE, col = "red", lwd = 2,lty=2)
		contour(px1, px2, matrix(func, 20, 20), level = -0.5, add = TRUE, col = "blue", lwd = 2,lty=2)

	}
}
