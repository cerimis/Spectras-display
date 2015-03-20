#######################################################
# Imortiert Spektren aus den Datenbanken HITRAN und PNNL.
# Erlaubt die Convolution mit Gauss/Lorentz/Voigt und stellt die Spektren dar
#######################################################



require(playwith)
# Konstanten
r <- 1.3806*10^-23
na <- 6.022*10^23
# rm(list = ls())
#PNNL spektrum
spekp<-"/home/pdietiker/Dokumente/CH4/CH4_25T.TXT"
# Hitran Spektrum
spekh<-"/home/pdietiker/Dokumente/CH4/06_hit12a.par"
# Hitran oder PNNL
t.source<-"HITRAN"
# t.source<-"PNNL"
# Isotopomer
t.iso<-"All"
# Bereich des Spektrums
t.bereich<-c(2940,3000)
# kleinste Intensität
t.int<-1e-28
#Doppler verbreiterung?
t.doppler<-TRUE
#Collisional Broadening (Daten aus Hitran:"H", Weisskopf:"W", kein:"F")
# Temperatur
Temp <- 300
Druck#
press <- 1000
# Masse
M <- 18
# Linienprofil
t.profil="Voigt"
# Schrittweite
t.schritt<-0.1
#Gaussfunktion auf höhe normiert
t.Gaussfunktionh <- function(FWHM,Wellenzahl,Int,x0){Int*exp(-1/2*((Wellenzahl-x0)/(FWHM/(2*sqrt(2*log(2)))))^2)}
#Gaussfunktion auf fläche normiert
t.Gaussfunktionf <- function(FWHM,Wellenzahl,Int,x0){Int*1/((FWHM/(2*sqrt(2*log(2))))*sqrt(2*pi))*exp(-1/2*((Wellenzahl-x0)/(FWHM/(2*sqrt(2*log(2)))))^2)}
# Lorenzfunktion auf fläche normiert
t.Lorentzfunktionf<- function(FWHM,Wellenzahl, Int,x0){FWHM/(2*pi*((x0-Wellenzahl)^2+(FWHM/2)^2))}

#breite der Gausslinie (Doppler)
FWHMG <- 100*Sqrt(8*log(2)*k*T/(M/(10^3*na))]; 
# Breite der Lorenzlinie (collisional)
FWHML<-0.5

#Versuch eine Reale Lineshape-Funktion zu modelieren
t.realfunktionf <- function(FWHM,Wellenzahl,Int,x0){
  if(Wellenzahl>=(x0-FWHM/3)) Int*1/((FWHM/(2*sqrt(2*log(2))))*sqrt(2*pi))*exp(-1/2*((Wellenzahl-x0)/(FWHM/(2*sqrt(2*log(2)))))^2)
  else
    0
}
# Funtion zur bestimmung der anzahl Nachkomastellen einer Zahl
f.decimalnumcount<-function(x){stopifnot(class(x)=="character")
                             x<-gsub("(.*)(\\.)|([0]*$)","",x)
                             nchar(x)
}  

#Import des Spektrums
if (t.source=="PNNL") {
d.spektrum<-as.data.frame(read.table(spekp,sep=""))
colnames(d.spektrum)<-c("Wavenumber","Intensity")
d.spektrum2<-d.spektrum[d.spektrum$Wavenumber>=t.bereich[1]& d.spektrum$Wavenumber<=t.bereich[2],]
}
if (t.source=="HITRAN") {
  d.spektrum<-as.data.frame(read.fwf(spekh,widths=c(3,12,10,-15))[,1:3])
  colnames(d.spektrum)<-c("Isotopomer","Wavenumber","Intensity")
  # Filtern des Spektrums
  if(mode(t.iso)=="character"){
    d.spektrum1<-d.spektrum[d.spektrum$Wavenumber>=t.bereich[1]& d.spektrum$Wavenumber<=t.bereich[2] & d.spektrum$Intensity>=t.int,2:3] 
  }
  if(mode(t.iso)=="numeric"){
  d.spektrum1<-d.spektrum[is.element(d.spektrum$Isotopomer,t.iso) & d.spektrum$Wavenumber>=t.bereich[1]& d.spektrum$Wavenumber<=t.bereich[2] & d.spektrum$Intensity>=t.int,2:3]
  }
  # Runden der Position auf t.schritt
  d.spektrum1[,"Wavenumber"]=round(d.spektrum1[,"Wavenumber"],f.decimalnumcount(as.character(t.schritt)))
  # Addieren von intensitäten bei mehreren linien bei  gleicher position
  # Addiert die Intensitäten aller Linien auf einer Position. Die Summe wird bei der ersten Linie der Position gespeichert
  # Dieser Schritt ist geschwindigkeitsbestimmend
  yy<-1
  y<-1
  while(y<dim(d.spektrum1)[1]){
    yy<-y
    while(d.spektrum1$Wavenumber[yy+1]==d.spektrum1$Wavenumber[yy] &yy<dim(d.spektrum1)[1]) {
      d.spektrum1$Intensity[y]<-d.spektrum1$Intensity[y]+d.spektrum1$Intensity[yy+1]
      yy<-yy+1
    }
    y<-yy+1
  }
  # Löschen aller Linien auf einer Position ausser der ersten
  d.spektrum1<-d.spektrum1[!duplicated(d.spektrum1$Wavenumber),]
  # Erstellen eines leeren Spektrum über den Bereich t.bereich mit Schrittweite t.schritt
  d.spektrum2<-as.data.frame(cbind(seq(t.bereich[1],t.bereich[2],by=t.schritt),rep(0,length(seq(t.bereich[1],t.bereich[2],by=t.schritt)))))
  colnames(d.spektrum2)<-c("Wavenumber","Intensity")
  # Einfügen der Linien
    for (n in 1:dim(d.spektrum1)[1]){
    d.spektrum2[d.spektrum2$Wavenumber==d.spektrum1$Wavenumber[n],"Intensity"]<-d.spektrum1$Intensity[n]
  }
  #invertieren der Tabelle
  d.spektrum2<-sort(d.spektrum2)
}
# Einfügen der Lorenzlinien
#Lorenzlinien mit maximas bei jedem Punkt von t.referenz bei allen Punkten von t.freq ausgewertet und anschliessend aufsummiert.
#Falls t.referenz weniger als 1000 linien enth?lt wir dies mit apply und anschliessend rowsum gemacht
#Bei mehr als 1000 linien wir in einer schlaufe jede linie einzeln bereichnet und anschliessend zum vorhergehenden addiert
#Falls imort==TRUE wird anstelle der Faltung das Referenzspektrum importiert 
if(dim(d.spektrum)[1]<1000)
  d.spektrum2<-rowSums(apply(t.referenz,1,t.Gaussfunktion,sigma=fwhm,Int="Fluoreszenz",x0="Wellenzahl",x=t.freq))
else
{
  t.reftab<-rep(0,length(t.freq)) 
  
  for(i in 1:dim(t.referenz)[1])
  {
    t.peak<-t.Gaussfunktion2(fwhm,t.freq,t.referenz$Fluoreszenz[i],t.referenz$Wellenzahl[i])
    t.reftab<-t.reftab+t.peak
    rm(t.peak)
  }
}
#sicherstellen dass die Länge des Spektrums durch 2 teilbar ist 
if(length(d.spektrum2)%%2)
  d.spektrum2<-d.spektrum2[-length(d.spektrum2)]


#Test stick-specktrum 
# d.spektrum<-cbind(seq(-1000,1000,by=0.1),c(rep(0,10000),1,rep(0,10000)))
# t.spekplot<-as.data.frame(d.spektrum[,1:2])
# t.spekplot<-t.spekplot[-c(60000,length(t.spekplot)),]


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
# t.Gausstab<-cbind(seq(-t.span/2,t.span/2,by=t.stepwidth),0)
# for (n in 1:length(seq(-t.span/2,t.span/2,by=t.stepwidth))){
#   t.Gausstab[n,2]<-t.realfunktionf(FWHM,t.Gausstab[n,1],1,0)
# }
# t.Gausstab<-t.Gausstab[,2]


# Convolution
if (t.profil=="Gauss") d.spektrumc<-convolve(d.spektrum2[,2],t.lsGtab,type="filter")
if (t.profil=="Lorentz") d.spektrumc<-convolve(d.spektrum2[,2],t.lsLtab,type="filter")
if (t.profil=="Voigt") d.spektrumc<-convolve(convolve(d.spektrum2[,2],t.lsGtab,type="filter"),t.lsLtab,type="filter")


# längen des Spektrums und der Convolution
ls=length(d.spektrum2[,1])
lc=length(d.spektrumc)

# Korrekture der Frequenzachse der Convolution
# t.spekplotconvolve1<-as.data.frame(cbind(d.spektrum[seq((ls-lc)/2+2,(lc+(ls-lc)/2+1)),1],6*FWHM/length(t.Gausstab)*t.spekplotconvolve))
if (t.profile=="Voigt")
  d.spektrumc2<-as.data.frame(cbind(d.spektrum2[seq((ls-lc)/2+2,(lc+(ls-lc)/2+1)),1],t.stepwidth^2*d.spektrumc))
  else
  d.spektrumc2<-as.data.frame(cbind(d.spektrum2[seq((ls-lc)/2+2,(lc+(ls-lc)/2+1)),1],t.stepwidth*d.spektrumc))
# Zeichnen des Spektrums
if (t.source=="PNNL") {
playwith({
plot(d.spektrum2[,1:2],type="l",ylim=c(0,max(d.spektrum2[,2]*1.1)))
lines(d.spektrumc2[,1:2],type="l",col="red",lwd=2)
},time.mode = TRUE, new=TRUE)
}
if (t.source=="HITRAN") {
  playwith({1
    plot(d.spektrum2[,1:2],type="h",ylim=c(0,max(d.spektrum2[,2]*1.1)))
    lines(d.spektrumc2[,1:2],type="l",col="red",lwd=2)
  },time.mode = TRUE, new=TRUE)
}

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