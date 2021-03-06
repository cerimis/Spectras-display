#######################################################
# Imortiert Spektren aus den Datenbanken PNNL.
# Erlaubt die Convolution mit Gauss/Lorentz/Voigt und stellt die Spektren diverser Gase überlagert dar.
#######################################################



require(playwith)
# Konstanten
k <- 1.3806*10^-23
na <- 6.022*10^23
c <- 3*10^8
# rm(list = ls())
#PNNL spektrum
d.daten<- data.frame(
  Name=c("Ethylenglykol","CH4","C2H6","C3H8","H2O","C4H10","CO2"),
  ppm=c(0,1,0,0,0,0,0),
  Spektrum=c("/home/pdietiker/Dokumente/PNNL/Ethylene_glycol/ETOHOH_50T.TXT",
             "C://Daten//Laser//CH4_25T.TXT",
             "/home/pdietiker/Dokumente/PNNL/Ethane/C2H6_25T.TXT",
             "/home/pdietiker/Dokumente/PNNL/Propane/C3H8_25T.TXT",
             "C://Daten//Laser//H2O_25T.TXT",
             "/home/pdietiker/Dokumente/PNNL/n-Butane/Butane_25T.TXT",
             "/home/pdietiker/Dokumente/PNNL/Carbon_dioxide/CO2_25T.TXT"
  ),
  stringsAsFactors = FALSE
)

# Bereich des Spektrums
t.bereich<-c(2800,3200)

# Linienprofil (Gauss, Voigt, Lorenz,für keine Faltung:irgend ein anderer String)
t.profil="Gausss"
#breite der Gausslinie (Doppler)
FWHMG <-1
# Breite der Lorenzlinie
FWHML<-0.1
# x-Achse transformieren (0 für Wellenzahl, 1 für Piezospannung , 2 für Wellenlänge)
t.xtransform <- "1"
#2 Punkte der linearen Kalibrationskurve (Spannung, Wellenzahl)
t.kal <-rbind(c(0,3086),c(2,3086-220))
# Als Transmission anzeigen?
t.transmission <- TRUE
# Zellenlänge / m
t.laenge <- 36
# Schrittweite
t.schritt<-0.1
#Gaussfunktion auf fläche normiert
t.Gaussfunktionf <- function(FWHM,Wellenzahl,Int,x0){Int*1/((FWHM/(2*sqrt(2*log(2))))*sqrt(2*pi))*exp(-1/2*((Wellenzahl-x0)/(FWHM/(2*sqrt(2*log(2)))))^2)}
# Lorenzfunktion auf fläche normiert
t.Lorentzfunktionf<- function(FWHM,Wellenzahl, Int,x0){Int*FWHM/(2*pi*((x0-Wellenzahl)^2+(FWHM/2)^2))}

#Definition der Funktion zum einschalten der Summe in der Playwith Graphik
f.Sum_handler_ein<-function(widget, playState){
  
  lines(d.spektren[[dim(d.spektren)[1],"Spektrum"]][,1:2],type="l",lwd=0.5)
}
#Definition der Funktion zum einschalten der Summe ohne erste Komponente in der Playwith Graphik
f.Sumdiff_handler_ein<-function(widget, playState){
  t.sumdiff<-rep(0,dim(d.daten)[1])
  for(i in 2:(dim(d.daten)[1]-1))
  {
    t.sumdiff<- t.sumdiff+ .GlobalEnv$d.spektren[,"Spektrum"][[i]]["Absorbance"]
  }
  .GlobalEnv$t.spekplot.sumdiff <- cbind(.GlobalEnv$d.spektren[1,"Spektrum"][[1]]["Wavenumber"],t.sumdiff+playState$env$Offset)
  playState$env$Sumdiffplot<-.GlobalEnv$t.spekplot.sumdiff
  callArg(playState, "data") <- quote(Sumdiffplot)
  #   playReplot(playState)
  lines(t.spekplot.sumdiff,type="l",lwd=0.5,col="red")
}
#Definition der Funktion zum ausschalten der Summe in einer Playwith graphik
f.Sum_handler_aus<-function(widget, playState){
  playReplot(playState)
}

#Playwith objekte zum ein und ausschalten der Sume
Sumein<-list("Sumein","gtk-execute","Summe Einblenden ",callback=f.Sum_handler_ein)
Sumaus<-list("Sumaus","gtk-execute","Summe Ausblenden ",callback=f.Sum_handler_aus)
Sumdiffein<-list("Sumdiffein","gtk-execute","Differenzsumme Einblenden ",callback=f.Sumdiff_handler_ein)
#löschen der Komponenten mit ppm=0
d.daten <- subset(d.daten,ppm!=0)
# Jedes Spektrum wird importiert, anschliessend auf t.bereich eingeschränkt, auf eine gerade Anzahl 
# Elemente gekürzt, gefaltet und in der liste d.import als data.frame gespeichert.
# das Data.frame s.Spektren wird aus d.daten und den importierten Spektren zusammengesetzt.
d.import<-list()
for(i in 1:dim(d.daten)[1])
{ d.import1<-as.data.frame(read.table(d.daten$Spektrum[i],sep=""))
  colnames(d.import1)<-c("Wavenumber","Absorbance")
  
  d.import1<-d.import1[d.import1$Wavenumber>=t.bereich[1]& d.import1$Wavenumber<=t.bereich[2],]
  
  if(dim(d.import1)[1]%%2)
  {  
    d.import1<-d.import1[-dim(d.import1)[1],]
  }
  d.import[[i]]<-d.import1  
}
d.spektren<-data.frame(Name=d.daten$Name,ppm=d.daten$ppm,Spektrum=I(d.import),
                       stringsAsFactors = FALSE)
rownames(d.spektren)<-d.daten$Name


#Tabelle der Faltungsfunktion
if (t.profil=="Gauss" | t.profil=="Voigt"){
  t.span<-4*FWHMG
  t.stepwidth<-0.1
  t.lsGtab<-t.Gaussfunktionf(FWHMG,seq(-t.span/2,t.span/2,by=t.stepwidth),1,0)
  #sicherstellen dass die Länge der tabelle durch 2 teilbar ist 
  if(length(t.lsGtab)%%2)
    t.lsGtab<-t.lsGtab[-length(t.lsGtab)]
}
if (t.profil=="Lorentz"| t.profil=="Voigt"){
  t.span<-12*FWHML
  t.stepwidth<-0.1
  t.lsLtab<-t.Lorentzfunktionf(FWHML,seq(-t.span/2,t.span/2,by=t.stepwidth),1,0)
  #sicherstellen dass die Länge der tabelle durch 2 teilbar ist 
  if(length(t.lsLtab)%%2)
    t.lsLtab<-t.lsLtab[-length(t.lsLtab)]
}

#  Anpassung der intensitäten
for(i in 1:dim(d.daten)[1])
{
  d.spektren[[i,"Spektrum"]]$Absorbance=d.spektren[[i,"Spektrum"]]$Absorbance*d.spektren$ppm[i]*t.laenge
}
#Summieren der Spektren und anf?gen an d.spektren
t.sum<-rep(0,dim(d.spektren[[1,"Spektrum"]]["Absorbance"])[1])
for(i in 1:dim(d.spektren)[1])
{
  t.sum<- t.sum + d.spektren[[i,"Spektrum"]]["Absorbance"]
}
t.sumspektrum<-data.frame(d.spektren[[1,"Spektrum"]]["Wavenumber"],t.sum)
d.spektren=rbind(d.spektren,d.spektren["CH4",])
d.spektren[dim(d.spektren)[1],"Name"]<-"Summe"
d.spektren[dim(d.spektren)[1],"ppm"]<-0
d.spektren[[dim(d.spektren)[1],"Spektrum"]]<-t.sumspektrum


# Convolution
for(i in 1:(dim(d.daten)[1]+1))
{
  if (t.profil=="Gauss")
  {
    d.spektrumc<-convolve(d.spektren[[i,"Spektrum"]]$Absorbance,t.lsGtab,type="filter")
    ls=length(d.spektren[[i,"Spektrum"]][[1]])
    lc=length(d.spektrumc)
    d.spektrumc2<-as.data.frame(cbind("Wavenumber"=d.spektren[[i,"Spektrum"]][seq((ls-lc)/2+2,(lc+(ls-lc)/2+1)),1],"Absorbance"=t.stepwidth*d.spektrumc))
    d.spektren[[i,"Spektrum"]]<-d.spektrumc2
  } 
  if (t.profil=="Lorentz")
  {
    d.spektrumc<-convolve(d.spektren[[i,"Spektrum"]]$Absorbance,t.lsLtab,type="filter")
    ls=length(d.spektren[[i,"Spektrum"]][[1]])
    lc=length(d.spektrumc)
    d.spektrumc2<-as.data.frame(cbind("Wavenumber"=d.spektren[[i,"Spektrum"]][seq((ls-lc)/2+2,(lc+(ls-lc)/2+1)),1],"Absorbance"=t.stepwidth*d.spektrumc))
    d.spektren[[i,"Spektrum"]]<-d.spektrumc2
  }
  if (t.profil=="Voigt")
  {
    d.spektrumc<-convolve(convolve(d.spektren[[i,"Spektrum"]]$Absorbance,t.lsGtab,type="filter"),t.lsLtab,type="filter")
    ls=length(d.spektren[[i,"Spektrum"]][[1]])
    lc=length(d.spektrumc)
    d.spektrumc2<-as.data.frame(cbind("Wavenumber"=d.spektren[[i,"Spektrum"]][seq((ls-lc)/2+2,(lc+(ls-lc)/2+1)),1],"Absorbance"=t.stepwidth^2*d.spektrumc))
    d.spektren[[i,"Spektrum"]]<-d.spektrumc2  
  } 
}


switch(t.xtransform,
       "1" = {for(i in 1:(dim(d.daten)[1]+1))
       {
         if(t.transmission)
         {
           colnames(d.spektren[[i,"Spektrum"]]) <- c("Voltage","Transmission")
         }else
         {
           colnames(d.spektren[[i,"Spektrum"]]) <- c("Voltage","Absorbance")
         }
         test <- d.spektren[[i,"Spektrum"]]
         test$Voltage <- (test$Voltage-min(t.kal[,2]))*-diff(range(t.kal[,1]))/(diff(range(t.kal[,2])))+max(t.kal[,1])
         test <- test[test$Voltage>=0 & test$Voltage<=2,]
         if(t.transmission) test$Transmission <- 10^(-test$Transmission)
         d.spektren[[i,"Spektrum"]] <- test
       }
       },
       "2" = {for(i in 1:(dim(d.daten)[1]+1))
       {
         if(t.transmission)
         {
           colnames(d.spektren[[i,"Spektrum"]]) <- c("Wavelength","Transmission")
         }else
         {
           colnames(d.spektren[[i,"Spektrum"]]) <- c("Wavelength","Absorbance")
         }
         test <- d.spektren[[i,"Spektrum"]]
         test$Wavelength <- 1/(100*test$Wavelength)*10^6
         if(t.transmission) test$Transmission <- 10^(-test$Transmission)
         d.spektren[[i,"Spektrum"]] <- test
       }
       },
       "0" = {for(i in 1:(dim(d.daten)[1]+1))
       {
         if(t.transmission)
         {
           colnames(d.spektren[[i,"Spektrum"]]) <- c("Wavenumber","Transmission")
         }else
         {
           colnames(d.spektren[[i,"Spektrum"]]) <- c("Wavenumber","Absorbance")
         }
         test <- d.spektren[[i,"Spektrum"]]
         if(t.transmission) test$Transmission <- 10^(-test$Transmission)
         d.spektren[[i,"Spektrum"]] <- test
       }
       }
)
#Maximale Intensität aller Spektren
d.maxint=max(unlist(lapply(d.spektren[,"Spektrum"],"[[",2)))
#Maximaler ppm-Wert
d.maxppm=max(unlist(lapply(d.spektren[,"ppm"],"[[",1)))
# längster Name
d.maxwidth=d.spektren$Name[which.max(unlist(lapply(unlist(lapply(d.spektren[,"Name"],"[[",1)),nchar, type="width")))]


# Zeichnen des Spektrums
playwith(
{
  plot(d.spektren[[1,"Spektrum"]][,1:2],type="l",ylim=c(0,d.maxint*1.1),col=palette()[1])
  
  for(i in 2:dim(d.daten)[1])
  {
    lines(d.spektren[[i,"Spektrum"]][,1:2],type="l",lwd=0.5,col=palette()[i])
  }
  templegend<-legend("topright",rep(" ",dim(d.daten)[1]+1),col=palette()[seq(1,dim(d.daten)[1]+1)],lty=1,pch=".",text.width = (strwidth(paste("aaaaaa",d.maxwidth,toString(d.maxppm)))))
  text(templegend$rect$left+templegend$rect$w,templegend$text$y,mapply(paste,d.spektren$ppm,"ppm",sep=" "),pos=2)
  text(templegend$text$x,templegend$text$y,d.spektren$Name,pos=4)
}
,time.mode = TRUE,new = TRUE,
parameters = list(Offset=100),
tools=list(Sumein,Sumaus,Sumdiffein))


##############################################################################3
# graphik zum Drucken
#playwith(
#{
#  plot(d.spektren[[1,"Spektrum"]][,1:2],type="l",ylim=c(0,d.maxint*1.1),col=palette()[1],
#       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
#  par("cex")[2]
#  
#  for(i in 2:dim(d.daten)[1])
#  {
#    lines(d.spektren[[i,"Spektrum"]][,1:2],type="l",lwd=0.5,col=palette()[i])
#  }
#  templegend<-legend("topright",rep(" ",dim(d.daten)[1]),col=palette()[seq(1,dim(d.daten)[1])],lty=1,lwd=2,pch=".",text.width = 2*(strwidth(paste("aaaaaa",d.maxwidth,toString(d.maxppm)))),y.intersp=2)
#  text(templegend$rect$left+templegend$rect$w,templegend$text$y,mapply(paste,d.spektren$ppm,"ppm",sep=" "),pos=2,cex=2)
#  text(templegend$text$x,templegend$text$y,d.spektren$Name,pos=4,cex=2)
#}
#,time.mode = TRUE,new = TRUE,
#parameters = list(Offset=100),
#tools=list(Sumein,Sumaus,Sumdiffein))


# interpoliert das spektrum an points/2 punkte,invertiert das Spektrum und hängt es an das urspüngliche an
points<- 4002
d.spektruminterpol <- approx (d.spektren[[dim(d.spektren)[1],"Spektrum"]],n=points/2)
d.spektruminterpol2 <- append(d.spektruminterpol$y,rev(d.spektruminterpol$y))
write.table(format(t(d.spektruminterpol2),digits=7),file="C://Daten//Laser//CH4_simulated.txt",sep="\t",row.names=FALSE, col.names=FALSE,quote=FALSE)



# Speicher der Spektren separat für PNNL und Hitran
# if (t.source=="HITRAN") d.spektrumc2hitran<-d.spektrumc2
# if (t.source=="PNNL") {
#   d.spektrumc2pnnl<-d.spektrumc2
#   d.spektrumc2pnnl[,2]<-d.spektrumc2pnnl[,2]*8e-17
# }
# Vergleich von Hitran und PNNL
# playwith({
#   plot(d.spektrumc2pnnl[,1:2],type="l",ylim=c(0,max(d.spektrumc2pnnl[,2]*1.1)))
#   lines(d.spektrumc2hitran[,1:2],type="l",col="red",lwd=2)
# },time.mode = TRUE, new=TRUE)