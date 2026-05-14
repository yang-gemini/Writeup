# Daten einlesen (beachten Sie, dass die Originaldatei oft kein Header hat)
# Falls die Datei keine Überschriften hat, nutzen Sie header=FALSE
dat <- read.table("pgalpga2008.txt", header=FALSE, col.names = c("Distance", "Precision", "FM"))

# Subset erstellen: Nur weibliche Golfer (FM == 1 laut Bildbeschreibung)
datF <- subset(dat, FM = 2)

# Die Regression anpassen
fit <- lm(Precision ~ Distance, data = datF)

# Zusammenfassung anzeigen
s <- summary(fit)
print(s)
plot(fitted(fit), residuals(fit))