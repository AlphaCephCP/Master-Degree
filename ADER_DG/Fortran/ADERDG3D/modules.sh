#load modules for compiling PDESol
module swap PrgEnv-cray PrgEnv-intel
module load cray-tpsl
module load numlib/intel/mkl/2018.1

#user alias
alias myjob="qsub -u xitmtave"
alias clean="rm *.plt *.txt *.dat *.log *.e* *.o* *~ "

