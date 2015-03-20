
require(playwith)

#PNNL spektrum
spek<-"/home/pdietiker/Dokumente/CH4/CH4_25T.TXT"
# Hitran Spektrum
spek<-"/home/pdietiker/Dokumente/CH4/06_hit12.par"
# Hitran oder PNNL
t.source<-"HITRAN"
# Isotopomer
t.iso<-61
# Bereich des Spektrums
t.bereich
#breite der Gauslinie
FWHM<-1
#Gaussfunktion auf höhe normiert
t.Gaussfunktionh <- function(FWHM,Wellenzahl,Int,x0){Int*exp(-1/2*((Wellenzahl-x0)/(FWHM/(2*sqrt(2*log(2)))))^2)}
#Gaussfunktion auf fläche normiert
t.Gaussfunktionf <- function(FWHM,Wellenzahl,Int,x0){Int*1/((FWHM/(2*sqrt(2*log(2))))*sqrt(2*pi))*exp(-1/2*((Wellenzahl-x0)/(FWHM/(2*sqrt(2*log(2)))))^2)}

#Import des Spektrums
d.spektrum<-read.table(spek,sep="",skip=0)
if (t.source=="HITRAN") {
  d.spektrum<-read.fwf(spek,widths=c(3,12,10,-15))[,1:3]
  colnames(d.spektrum)<-c("Isotopomer","Wavenumber","Intensity")
  d.spektrum1<-d.spektrum[d.spektrum$Isotopomer==t.iso]
  
}
}
# colClasses=c(rep("numeric",10),rep("character",4),rep("numeric",2),"character",rep("numeric",2)),
# ,widths=c(3,1,12,10,10,5,5,10,4,8,15,15,15,15,6,12,1,7,7)
#Test stick-specktrum 
# d.spektrum<-cbind(seq(-1000,1000,by=0.1),c(rep(0,10000),1,rep(0,10000)))
t.spekplot<-as.data.frame(d.spektrum[,1:2])
t.spekplot<-t.spekplot[-c(60000,length(t.spekplot)),]
#sicherstellen dass die Länge des Spektrums durch 2 teilbar ist 
if(length(t.spekplot)%%2)
  t.spekplot<-t.spekplot[-length(t.spekplot)]

#Tabelle der Gaussfunktion 
t.span<-4*FWHM
t.stepwidth<-0.1
# t.Gausstab<-t.Gaussfunktionf(FWHM,seq(-4*FWHM,4*FWHM,by=0.1),1,0)
t.Gausstab<-t.Gaussfunktionf(FWHM,seq(-t.span/2,t.span/2,by=t.stepwidth),1,0)
#sicherstellen dass die Länge der Gausstabelle durch 2 teilbar ist 
if(length(t.Gausstab)%%2)
   t.Gausstab<-t.Gausstab[-length(t.Gausstab)]

# Convolution
t.spekplotconvolve<-convolve(d.spektrum[,2],t.Gausstab,type="filter")

# längen des Spektrums und der Convolution
ls=length(d.spektrum[,1])
lc=length(t.spekplotconvolve)

# Korrekture der Frequenzachse der Convolution
# t.spekplotconvolve1<-as.data.frame(cbind(d.spektrum[seq((ls-lc)/2+2,(lc+(ls-lc)/2+1)),1],6*FWHM/length(t.Gausstab)*t.spekplotconvolve))
t.spekplotconvolve1<-as.data.frame(cbind(d.spektrum[seq((ls-lc)/2+2,(lc+(ls-lc)/2+1)),1],t.stepwidth*t.spekplotconvolve))
# Zeichnen des Spektrums
playwith({
plot(t.spekplot[,1:2],type="l",ylim=c(0,max(t.spekplot[,2]*1.1)))
lines(t.spekplotconvolve1[,1:2],type="l",col="red",lwd=2)
},time.mode = TRUE, new=FALSE)

