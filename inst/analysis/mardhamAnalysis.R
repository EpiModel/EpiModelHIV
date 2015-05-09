
if(F) {
##############################################################

mean(mardham.sim.01a$atts.list[[520]]$vl[mardham.sim.01a$atts.list[[520]]$stage%in%"AF"])
mean(mardham.sim.01a$atts.list[[520]]$vl[mardham.sim.01a$atts.list[[520]]$stage%in%"AR"])
mean(mardham.sim.01a$atts.list[[520]]$vl[mardham.sim.01a$atts.list[[520]]$stage%in%"C"])
mean(mardham.sim.01a$atts.list[[520]]$vl[mardham.sim.01a$atts.list[[520]]$stage%in%"D"])

###############################################################

times <- 1:length(mardham.sim.01b$disc.ai)

discai <- sapply(times, function(x) nrow(mardham.sim.01b$disc.ai[[x]]))
disc.uai <- sapply(times, function(x) sum(mardham.sim.01b$disc.ai[[x]]$uai==1))
mean.vl <- sapply(times, function(x) mean(mardham.sim.01b$atts.list[[x]]$vl, na.rm=T))
percent.aids <- sapply(times, function(x) mean(mardham.sim.01b$atts.list[[x]]$stage=="D", na.rm=T))
percent.chronic <- sapply(times, function(x) mean(mardham.sim.01b$atts.list[[x]]$stage=="C", na.rm=T))
percent.chronic <- sapply(times, function(x) mean(mardham.sim.01b$atts.list[[x]]$stage=="C", na.rm=T))

################################################################

window <- 200
plot(filter(
  mardham.sim.01b$summ$i - mardham.sim.01b$summ$dg/5 - mardham.sim.01b$summ$da,
  rep(1,window))/window)
lines(c(0,10000),c(0,0))

points(filter(diff(mardham.sim.01b$summ$p/mardham.sim.01b$summ$n)*10, rep(1,window)))

pdf("basics.pdf")
plot(mardham.sim.01b$summ$p, main='absolute prevalence')
plot(mardham.sim.01b$summ$p/mardham.sim.01b$summ$n, main='prevalence')
plot(mardham.sim.01b$summ$n, main='popsize')

window <- 200
plot(filter(
  mardham.sim.01b$summ$i - mardham.sim.01b$summ$dg*.23 - mardham.sim.01b$summ$da,
  rep(1,window))/window)
lines(c(0,10000),c(0,0))

dev.off()



#####

diff(mardham.sim.01c$summ$p) - 
(mardham.sim.01c$summ$i-mardham.sim.01c$summ$da-mardham.sim.01c$summ$dgp)[-1]

window<-25
plot(filter(mardham.sim.01c$summ$p, rep(1,window))/window) 
plot(filter(mardham.sim.01c$summ$p/mardham.sim.01c$summ$n, rep(1,window))/window) 
plot(filter(mardham.sim.01c$summ$i, rep(1,window))/window) 
plot(filter(mardham.sim.01c$summ$da, rep(1,window))/window) 
plot(filter(mardham.sim.01c$summ$dgp, rep(1,window))/window) 

timestep <- 1300

table(mardham.sim.01c$atts.list[[timestep]]$
        tt.traj[mardham.sim.01c$atts.list[[timestep]]$stage=="D"])
hist(mardham.sim.01c$atts.list[[timestep]]$
        inf.time[mardham.sim.01c$atts.list[[timestep]]$stage=="D"])
max(mardham.sim.01c$atts.list[[timestep]]$
        inf.time[mardham.sim.01c$atts.list[[timestep]]$stage=="D"],na.rm=T)
min(mardham.sim.01c$atts.list[[timestep]]$
       inf.time[mardham.sim.01c$atts.list[[timestep]]$stage=="D"],na.rm=T)
mean(mardham.sim.01c$atts.list[[timestep]]$
       inf.time[mardham.sim.01c$atts.list[[timestep]]$stage=="D"],na.rm=T)

aaa <- which(mardham.sim.01c$atts.list[[timestep]]$stage=="D")

mardham.sim.01c$atts.list[[timestep]]$stage[aaa]

bbb <- data.frame(
  a = mardham.sim.01c$atts.list[[timestep]]$inf.time[aaa],
  b = mardham.sim.01c$atts.list[[timestep]]$tt.traj[aaa],
  c = mardham.sim.01c$atts.list[[timestep]]$cum.time.off.tx[aaa],
  d = mardham.sim.01c$atts.list[[timestep]]$cum.time.on.tx[aaa],
  e = mardham.sim.01c$atts.list[[timestep]]$stage.time[aaa],
  f = mardham.sim.01c$atts.list[[timestep]]$uid[aaa]
  )

xxx <- 8
yyy <- mardham.get.id.from.uid(mardham.sim.01c, xxx)
zzz <- mardham.get.att.from.uid(mardham.sim.01c, 'stage.time', yyy, xxx)
plot(mardham.get.att.from.uid(mardham.sim.01c, 'cum.time.off.tx', yyy, xxx))
plot(mardham.get.att.from.uid(mardham.sim.01c, 'vl', yyy, xxx))
         
mean(mardham.sim.01c$atts.list[[1]]$vl, na.rm=T)
mean(mardham.sim.01c$atts.list[[1000]]$vl, na.rm=T)
       
##########################

plot(mardham.sim.01f$summ$p/mardham.sim.01f$summ$n, ylim=c(0,0.3))
points(mardham.sim.01c$summ$p/mardham.sim.01c$summ$n, col='red')

plot(mardham.sim.01f$summ$p.B/mardham.sim.01f$summ$n.B, ylim=c(0,0.5))
points(mardham.sim.01f$summ$p.W/mardham.sim.01f$summ$n.W, col='red')

plot(mardham.sim.01f$summ$n.B)
points(mardham.sim.01f$summ$n.W, col='red')

mardham.meanstats.01$meanstats.p[[1]]*2/network.size(mardham.basepop.01$nD.main)
mardham.meanstats.01$meanstats.i[[1]]*2/network.size(mardham.basepop.01$nD.main)

plot(mardham.sim.01$summ$md.MB, ylim=c(0,1))
points(mardham.sim.01$summ$md.MW, col='red')
lines(c(0,10000), rep(mardham.meanstats.01$meanstats.m[[1]]*2/
                        network.size(mardham.basepop.01$nD.main),2))

plot(mardham.sim.01$summ$md.PB, ylim=c(0,1))
points(mardham.sim.01$summ$md.PW, col='red')
lines(c(0,10000), rep(mardham.meanstats.01$meanstats.p[[1]]*2/
                        network.size(mardham.basepop.01$nD.main),2))

plot(mardham.sim.01$summ$md.IB, ylim=c(0,1))
points(mardham.sim.01$summ$md.IW, col='red')
lines(c(0,10000), rep(sum(mardham.meanstats.01$meanstats.i[1:12])/
                        network.size(mardham.basepop.01$nD.main),2))

######################

discai <- sapply(1:length(mardham.sim.01$disc.ai), function(x) nrow(mardham.sim.01$disc.ai[[x]]))
discai.a <- sapply(1:length(mardham.sim.01a$disc.ai), function(x) nrow(mardham.sim.01a$disc.ai[[x]]))

plot(discai)
points(discai.a,col='red')

disc.uai <- sapply(1:length(mardham.sim.01$disc.ai), function(x) sum(mardham.sim.01$disc.ai[[x]]$uai==1))
disc.uai.a <- sapply(1:length(mardham.sim.01a$disc.ai), function(x) sum(mardham.sim.01a$disc.ai[[x]]$uai==1))

plot(disc.uai)
points(disc.uai.a,col='red')

mean(mardham.sim.01a$atts.list[[1]]$vl, na.rm=T)
mean(mardham.sim.01a$atts.list[[520]]$vl, na.rm=T)

mean(mardham.sim.01a$atts.list[[1]]$circ, na.rm=T)
mean(mardham.sim.01a$atts.list[[520]]$circ, na.rm=T)

table(mardham.sim.01a$atts.list[[1]]$stage)
table(mardham.sim.01a$atts.list[[520]]$stage)

mean(mardham.sim.01a$atts.list[[1]]$vl[mardham.sim.01a$atts.list[[1]]$stage%in%"AF"])
mean(mardham.sim.01a$atts.list[[1]]$vl[mardham.sim.01a$atts.list[[1]]$stage%in%"AR"])
mean(mardham.sim.01a$atts.list[[1]]$vl[mardham.sim.01a$atts.list[[1]]$stage%in%"C"])
mean(mardham.sim.01a$atts.list[[1]]$vl[mardham.sim.01a$atts.list[[1]]$stage%in%"D"])
mean(mardham.sim.01a$atts.list[[520]]$vl[mardham.sim.01a$atts.list[[520]]$stage%in%"AF"])
mean(mardham.sim.01a$atts.list[[520]]$vl[mardham.sim.01a$atts.list[[520]]$stage%in%"AR"])
mean(mardham.sim.01a$atts.list[[520]]$vl[mardham.sim.01a$atts.list[[520]]$stage%in%"C"])
mean(mardham.sim.01a$atts.list[[520]]$vl[mardham.sim.01a$atts.list[[520]]$stage%in%"D"])

mean(mardham.sim.01a$atts.list)

plot(filter(mardham.sim.01$summ$i, rep(1,25))/25)
plot(filter(mardham.sim.01$summ$dg, rep(1,25))/25)
plot(filter(mardham.sim.01$summ$da, rep(1,25))/25)


###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################

plot(mardham.sim.01b$summ$p/mardham.sim.01b$summ$n,ylim=c(0,0.3))
plot(mardham.sim.01b$summ$n)

plot(mardham.sim.01b$summ$md.MB, ylim=c(0,1))
points(mardham.sim.01b$summ$md.MW, col='red')
lines(c(0,10000), rep(mardham.meanstats.01$meanstats.m[[1]]*2/
                        network.size(mardham.basepop.01$nD.main),2))

plot(mardham.sim.01b$summ$md.PB, ylim=c(0,1))
points(mardham.sim.01b$summ$md.PW, col='red')
lines(c(0,10000), rep(mardham.meanstats.01$meanstats.p[[1]]*2/
                        network.size(mardham.basepop.01$nD.main),2))

plot(mardham.sim.01b$summ$md.IB, ylim=c(0,1))
points(mardham.sim.01b$summ$md.IW, col='red')
lines(c(0,10000), rep(sum(mardham.meanstats.01$meanstats.i[1:12])/
                        network.size(mardham.basepop.01$nD.main),2))

######################

times <- 1:length(mardham.sim.01b$disc.ai)

discai <- sapply(times, function(x) nrow(mardham.sim.01b$disc.ai[[x]]))
plot(discai)

disc.uai <- sapply(times, 
                   function(x) sum(mardham.sim.01b$disc.ai[[x]]$uai==1))
plot(disc.uai)

mean.vl <- sapply(times, 
                  function(x) mean(mardham.sim.01b$atts.list[[x]]$vl, na.rm=T))
plot(mean.vl)
percent.aids <- sapply(times,
                       function(x) mean(mardham.sim.01b$atts.list[[x]]$stage=="D", na.rm=T))
plot(percent.aids)  

percent.chronic <- sapply(times,
                          function(x) mean(mardham.sim.01b$atts.list[[x]]$stage=="C", na.rm=T))
plot(percent.chronic)  

percent.chronic <- sapply(times,
                          function(x) mean(mardham.sim.01b$atts.list[[x]]$stage=="C", na.rm=T))
plot(percent.chronic)  

plot(filter(mardham.sim.01b$summ$i, rep(1,25))/25)
plot(filter(mardham.sim.01b$summ$dg/5, rep(1,25))/25)
plot(filter(mardham.sim.01b$summ$da, rep(1,25))/25)

window <- 200
plot(filter(
  mardham.sim.01b$summ$i - mardham.sim.01b$summ$dg/5 - mardham.sim.01b$summ$da,
  rep(1,window))/window)
lines(c(0,10000),c(0,0))

points(filter(diff(mardham.sim.01b$summ$p/mardham.sim.01b$summ$n)*10, rep(1,window)))

pdf("basics.pdf")
plot(mardham.sim.01b$summ$p, main='absolute prevalence')
plot(mardham.sim.01b$summ$p/mardham.sim.01b$summ$n, main='prevalence')
plot(mardham.sim.01b$summ$n, main='popsize')

window <- 200
plot(filter(
  mardham.sim.01b$summ$i - mardham.sim.01b$summ$dg*.23 - mardham.sim.01b$summ$da,
  rep(1,window))/window)
lines(c(0,10000),c(0,0))

dev.off()



#####

diff(mardham.sim.01c$summ$p) - 
  (mardham.sim.01c$summ$i-mardham.sim.01c$summ$da-mardham.sim.01c$summ$dgp)[-1]

window<-25
plot(filter(mardham.sim.01c$summ$p, rep(1,window))/window) 
plot(filter(mardham.sim.01c$summ$p/mardham.sim.01c$summ$n, rep(1,window))/window) 
plot(filter(mardham.sim.01c$summ$i, rep(1,window))/window) 
plot(filter(mardham.sim.01c$summ$da, rep(1,window))/window) 
plot(filter(mardham.sim.01c$summ$dgp, rep(1,window))/window) 

timestep <- 1300

table(mardham.sim.01c$atts.list[[timestep]]$
        tt.traj[mardham.sim.01c$atts.list[[timestep]]$stage=="D"])
hist(mardham.sim.01c$atts.list[[timestep]]$
       inf.time[mardham.sim.01c$atts.list[[timestep]]$stage=="D"])
max(mardham.sim.01c$atts.list[[timestep]]$
      inf.time[mardham.sim.01c$atts.list[[timestep]]$stage=="D"],na.rm=T)
min(mardham.sim.01c$atts.list[[timestep]]$
      inf.time[mardham.sim.01c$atts.list[[timestep]]$stage=="D"],na.rm=T)
mean(mardham.sim.01c$atts.list[[timestep]]$
       inf.time[mardham.sim.01c$atts.list[[timestep]]$stage=="D"],na.rm=T)

aaa <- which(mardham.sim.01c$atts.list[[timestep]]$stage=="D")

mardham.sim.01c$atts.list[[timestep]]$stage[aaa]

bbb <- data.frame(
  a = mardham.sim.01c$atts.list[[timestep]]$inf.time[aaa],
  b = mardham.sim.01c$atts.list[[timestep]]$tt.traj[aaa],
  c = mardham.sim.01c$atts.list[[timestep]]$cum.time.off.tx[aaa],
  d = mardham.sim.01c$atts.list[[timestep]]$cum.time.on.tx[aaa],
  e = mardham.sim.01c$atts.list[[timestep]]$stage.time[aaa],
  f = mardham.sim.01c$atts.list[[timestep]]$uid[aaa]
)

xxx <- 8
yyy <- mardham.get.id.from.uid(mardham.sim.01c, xxx)
zzz <- mardham.get.att.from.uid(mardham.sim.01c, 'stage.time', yyy, xxx)
plot(mardham.get.att.from.uid(mardham.sim.01c, 'cum.time.off.tx', yyy, xxx))
plot(mardham.get.att.from.uid(mardham.sim.01c, 'vl', yyy, xxx))

         
mean(mardham.sim.01c$atts.list[[1]]$vl, na.rm=T)
mean(mardham.sim.01c$atts.list[[1000]]$vl, na.rm=T)
       
       ##########################
       
plot(mardham.sim.01$summ$p/mardham.sim.01$summ$n, ylim=c(0,0.3))
points(mardham.sim.01$summ$p.B/mardham.sim.01$summ$n.B, col ='red')
points(mardham.sim.01$summ$p.W/mardham.sim.01$summ$n.W, col ='blue')

aaa <- xtabs(~floor(mardham.sim.01$atts.curr$age)+
        mardham.sim.01$atts.curr$inf.status)
plot(aaa[,2]/rowSums(aaa))

bbb <- xtabs(~floor(mardham.basepop.01$atts.curr$age)+
                 mardham.basepop.01$atts.curr$inf.status)
points(bbb[,2]/rowSums(bbb), col='red')

#######################

aaa <- xtabs(~mardham.sim.01$atts.curr$diag.status+
      mardham.sim.01$atts.curr$tx.status+
      mardham.sim.01$atts.curr$race)

bbb <- xtabs(~mardham.sim.01$atts.curr$tx.status+
       (mardham.sim.01$atts.curr$vl==1.5)+
       mardham.sim.01$atts.curr$race)

aaa[2,2,]/colSums(aaa[2,,])
bbb[2,2,]/colSums(bbb[2,,])

#######################

ccc <- xtabs(~mardham.basepop.01$atts.curr$diag.status+
               mardham.basepop.01$atts.curr$tx.status+
               mardham.basepop.01$atts.curr$race)

ddd <- xtabs(~mardham.basepop.01$atts.curr$tx.status+
               (mardham.basepop.01$atts.curr$vl==1.5)+
               mardham.basepop.01$atts.curr$race)

ccc[2,2,]/colSums(ccc[2,,])
ddd[2,2,]/colSums(ddd[2,,])

###############

uids <- 1:3000

#stage.time <- matrix(NA, length(uids), 520)
#cum.time.on.tx <- matrix(NA, length(uids), 520)
#cum.time.off.tx <- matrix(NA, length(uids), 520)
tx.status <- matrix(NA, length(uids), 520)
stage <- matrix(NA, length(uids), 520)

for (uid in uids) {
  id <- mardham.get.id.from.uid(mardham.sim.01.test.tx, uid)
  if (sum(!is.na(id))>0) {
#    stage.time[uid,] <- mardham.get.att.from.uid(mardham.sim.01.test.tx, 'stage.time', id, uid)  
#    cum.time.on.tx[uid,] <- mardham.get.att.from.uid(mardham.sim.01.test.tx, 'cum.time.on.tx', id, uid)  
#    cum.time.off.tx[uid,] <- mardham.get.att.from.uid(mardham.sim.01.test.tx, 'cum.time.off.tx', id, uid)  
    tx.status[uid,] <- mardham.get.att.from.uid(mardham.sim.01.test.tx, 'tx.status', id, uid)  
    stage[uid,] <- mardham.get.att.from.uid(mardham.sim.01.test.tx, 'stage', id, uid)  
  }
  cat(uid,'\n')
}

chronic <- matrix(stage%in%"C", nrow=length(uids))
chronic <- chronic[, -ncol(chronic)]
offoff <- sapply(1:519, function(x) (tx.status[,x]%in%0) & (tx.status[,x+1]%in%0))
offon <- sapply(1:519, function(x) (tx.status[,x]%in%0) & (tx.status[,x+1]%in%1))
onoff <- sapply(1:519, function(x) (tx.status[,x]%in%1) & (tx.status[,x+1]%in%0))
onon <- sapply(1:519, function(x) (tx.status[,x]%in%1) & (tx.status[,x+1]%in%1))

chronic.sum <- colSums(chronic)
offoff.sum <- colSums(chronic & offoff)
offon.sum <- colSums(chronic & offon)
onoff.sum <- colSums(chronic & onoff)
onon.sum <- colSums(chronic & onon)

offoff.prop <- offoff.sum/chronic.sum
offon.prop <- offon.sum/chronic.sum
onoff.prop <- onoff.sum/chronic.sum
onon.prop <- onon.sum/chronic.sum

mean(offon.prop/ (offoff.prop + offon.prop))
mean(onoff.prop/ (onoff.prop + onon.prop))

######

prop.on.tx <- sapply(1:520, function(x) 
  sum(mardham.sim.01.test.tx$atts.list[[x]]$tx.status, na.rm=T) / 
  sum(mardham.sim.01.test.tx$atts.list[[x]]$diag.status, na.rm=T) )

full.prop.on.tx <- sapply(1:520, function(x) 
  sum(mardham.sim.01.test.tx$atts.list[[x]]$tx.status==1 & 
        mardham.sim.01.test.tx$atts.list[[x]]$tt.traj=="YF", na.rm=T) / 
    sum(mardham.sim.01.test.tx$atts.list[[x]]$diag.status==1 & 
          mardham.sim.01.test.tx$atts.list[[x]]$tt.traj=="YF", na.rm=T) )

part.prop.on.tx <- sapply(1:520, function(x) 
  sum(mardham.sim.01.test.tx$atts.list[[x]]$tx.status==1 & 
        mardham.sim.01.test.tx$atts.list[[x]]$tt.traj=="YP", na.rm=T) / 
    sum(mardham.sim.01.test.tx$atts.list[[x]]$diag.status==1 & 
          mardham.sim.01.test.tx$atts.list[[x]]$tt.traj=="YP", na.rm=T) )

}