export MODULES_RUN_QUARANTINE=LD_PRELOAD
module load Anaconda3/2024.02-1
module load GCCcore/12.3.0
module load GCC/12.3.0
module load OpenBLAS/0.3.23-GCC-12.3.0 
export BLASLIB=-lopenblas
export BLASDIR=${LD_LIBRARY_PATH}
