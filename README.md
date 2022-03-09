# GQS
The Scripts is used to calculate GQS for single locus (on single file) or on given list of multiple locus (from summary statistics). 

The script should run with Python 3.8.X

usage: GQSmain.py [-h] [--ifile IFILE] [--regs REGS] [--r2_th R2_TH] [--chrm CHRM] [--chrm_h CHRM_H] [--pval_h PVAL_H] [--snp_h SNP_H] [--pos_h POS_H] [--ld_h LD_H]
                  [--refG {genome1000-EUR,genome1000-EAS,HRC-EUR,HRC-EAS}] [--addout ADDOUT]
                 
## Running on single file (with already estimated linkage disequilibrium (---ld_h-))
```
./GQSmain.py --ifile covidGWAS_chr7_test_file1.txt --r2_th 0.01 --chrm 7 --snp_h SNP --pval_h PVAL --ld_h RSQR --addout test_run
```
