#!/bin/bash -x
#----------------------------------------------------
# Example SLURM job script to run MPI applications on 
# TACC's Stampede system.
#
# $Id: job.mpi 1580 2013-01-08 04:10:50Z karl $
#----------------------------------------------------

#SBATCH -J job              # Job name
#SBATCH -o job.%j.out       # Name of stdout output file (%j expands to jobId)
#SBATCH -p development        # Queue name
#SBATCH -N 1                  # Total number of nodes requested (16 cores/node)
#SBATCH -n 16                 # Total number of mpi tasks requested
#SBATCH -t 02:00:00           # Run time (hh:mm:ss) - 1.5 hours

#SBATCH -A RAD-A-texana      # <-- Allocation name to charge job against
date
