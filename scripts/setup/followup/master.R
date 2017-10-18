
## build master.sh script ##

library("EpiModelHPC")

system("scp scripts/followup/*.[Rs]* hyak:/gscratch/csde/sjenness/sti/")

# Table 1 // Figures 1 and 2
vars <- list(COV = seq(0.1, 0.9, 0.1),
             PSTIINT = 182,
             RC = seq(0, 1, 0.1))
qsub_master(simno.start = 1000,
            nsubjobs = 16,
            backfill = TRUE,
            vars = vars,
            append = FALSE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")

# Table 2 // Figure 3
vars <- list(COV = 0.4,
             PSTIINT = c(30, 91, 182, 273, 364),
             RC = 0.4,
             PROBTX = 1)
qsub_master(simno.start = 2000,
            nsubjobs = 16,
            backfill = FALSE,
            vars = vars,
            append = TRUE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")

vars <- list(COV = 0.4,
             PSTIINT = 182,
             RC = 0.4,
             PROBTX = c(0, 0.25, 0.50, 0.75, 1))
qsub_master(simno.start = "auto",
            nsubjobs = 16,
            backfill = FALSE,
            vars = vars,
            append = TRUE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")

vars <- list(COV = 0.4,
             PSTIINT = 182,
             RC = 0.4,
             PROBTX = 0,
             ASYMPTX = c(0, 0.05, 0.1, 0.15, 0.2))
qsub_master(simno.start = "auto",
            nsubjobs = 16,
            backfill = FALSE,
            vars = vars,
            append = TRUE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")


# Supp Runs for Table 2

# STI test interval by risk compensation
vars <- list(COV = 0.4,
             PSTIINT = seq(30, 365, 30),
             RC = seq(0, 1, 0.1),
             PROBTX = 1,
             ASYMPTX = 0)
qsub_master(simno.start = 3000,
            nsubjobs = 16,
            backfill = FALSE,
            vars = vars,
            append = TRUE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")


# STI treatment prob by risk compensation
vars <- list(COV = 0.4,
             PSTIINT = 182,
             RC = seq(0, 1, 0.1),
             PROBTX = seq(0, 1, 0.1),
             ASYMPTX = 0)
qsub_master(simno.start = 4000,
            nsubjobs = 16,
            backfill = FALSE,
            vars = vars,
            append = TRUE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")
