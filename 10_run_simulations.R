# Run the simulation study.
source("00_load_functions.R")

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1L) args[1L] else "preliminary"
reps <- if (length(args) >= 2L) as.integer(args[2L]) else 5L
study <- run_study(profile = profile, reps = reps)
