
load("rmse_mod_colony.Rdata")
load("rmse_raw_colony.Rdata")
load("logrmse_mod_colony.Rdata")
load("logrmse_raw_colony.Rdata")
load("rmse_mod_shell.Rdata")
load("rmse_raw_shell.Rdata")
load("logrmse_mod_shell.Rdata")
load("logrmse_raw_shell.Rdata")


p.err=c(1e-5,0.0001,0.0003,0.001,0.003,0.01,0.03,0.1,0.3)

png("RMSE_compare.png",height=2000,width=3000)
par(cex=4)
plot(p.err,rmse.raw.colony,type="b",log="xy",
  xlab="Observational error probability",ylab="RMSE",
  ylim=c(min(c(rmse.mod.colony,rmse.raw.colony,rmse.mod.shell,rmse.raw.shell)),
	 max(c(rmse.mod.colony,rmse.raw.colony,rmse.mod.shell,rmse.raw.shell))),
	 col="green",lwd=6,axes=F, main="(a)")
lines(p.err,rmse.mod.colony,type="b",col="black",lwd=6)
lines(p.err,rmse.mod.shell,type="b",col="blue",lwd=6)
lines(p.err,rmse.raw.shell,type="b",col="red",lwd=6)
axis(1, at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.3),
  label=c("0","0.0001","0.001","0.01","0.1","0.3"),outer=F,lwd=6)
axis(2, at=seq(0.02,0.1,0.02),label=seq(0.02,0.1,0.02),outer=F,lwd=6)
box(lwd=6)
dev.off()

png("logRMSE_compare.png",height=2000,width=3000)
par(cex=4)
plot(p.err,logrmse.raw.colony,type="b",log="xy",
  xlab="Observational error probability",ylab="Log-RMSE",
  ylim=c(min(c(logrmse.mod.colony,logrmse.raw.colony,logrmse.mod.shell,logrmse.raw.shell)),
	 max(c(logrmse.mod.colony,logrmse.raw.colony,logrmse.mod.shell,logrmse.raw.shell))),
	 col="green",lwd=6,axes=F, main="(b)")
lines(p.err,logrmse.mod.colony,type="b",col="black",lwd=6)
lines(p.err,logrmse.mod.shell,type="b",col="blue",lwd=6)
lines(p.err,logrmse.raw.shell,type="b",col="red",lwd=6)
axis(1, at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.3),
  label=c("0","0.0001","0.001","0.01","0.1","0.3"),outer=F,lwd=6)
axis(2, at=c(0.3,0.5,1.0,1.5,1.9),
  label=c("","0.5","1.0","1.5",""),outer=F,lwd=6)
box(lwd=6)
dev.off()

png("RMSE_compare_both.png",height=2000,width=4000)
par(mfrow=c(1,2))
par(cex=5)
plot(p.err,rmse.raw.colony,type="b",log="xy",
  xlab="Observational error probability",ylab="RMSE",
  ylim=c(min(c(rmse.mod.colony,rmse.raw.colony,rmse.mod.shell,rmse.raw.shell)),
	 max(c(rmse.mod.colony,rmse.raw.colony,rmse.mod.shell,rmse.raw.shell))),
	 col="green",lwd=6,axes=F, main="(a)")
lines(p.err,rmse.mod.colony,type="b",col="black",lwd=6)
lines(p.err,rmse.mod.shell,type="b",col="blue",lwd=6)
lines(p.err,rmse.raw.shell,type="b",col="red",lwd=6)
axis(1, at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.3),
  label=c("0","0.0001","0.001","0.01","0.1","0.3"),outer=F,lwd=6)
axis(2, at=seq(0.02,0.1,0.02),label=seq(0.02,0.1,0.02),outer=F,lwd=6)
box(lwd=6)

par(cex=5)
plot(p.err,logrmse.raw.colony,type="b",log="xy",
  xlab="Observational error probability",ylab="Log-RMSE",
  ylim=c(min(c(logrmse.mod.colony,logrmse.raw.colony,logrmse.mod.shell,logrmse.raw.shell)),
	 max(c(logrmse.mod.colony,logrmse.raw.colony,logrmse.mod.shell,logrmse.raw.shell))),
	 col="green",lwd=6,axes=F, main="(b)")
lines(p.err,logrmse.mod.colony,type="b",col="black",lwd=6)
lines(p.err,logrmse.mod.shell,type="b",col="blue",lwd=6)
lines(p.err,logrmse.raw.shell,type="b",col="red",lwd=6)
axis(1, at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.3),
  label=c("0","0.0001","0.001","0.01","0.1","0.3"),outer=F,lwd=6)
axis(2, at=c(0.3,0.5,1.0,1.5,1.9),
  label=c("","0.5","1.0","1.5",""),outer=F,lwd=6)
box(lwd=6)
dev.off()

