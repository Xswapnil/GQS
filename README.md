# GQS
The Scripts is used to calculate GQS for single locus (on single file) or on given list of multiple locus (from summary statistics). 

The script should run with Python 3.8.X

usage: GQSmain.py [-h] [--ifile IFILE] [--regs REGS] [--r2_th R2_TH] [--chrm CHRM] [--chrm_h CHRM_H] [--pval_h PVAL_H] [--snp_h SNP_H] [--pos_h POS_H] [--ld_h LD_H]
                  [--refG {genome1000-EUR,genome1000-EAS,HRC-EUR,HRC-EAS}] [--addout ADDOUT]
                 
## Running on single file (with already estimated linkage disequilibrium (--ld_h))
```
./GQSmain.py --ifile covidGWAS_chr7_test_file1.txt --r2_th 0.01 --chrm 7 --snp_h SNP --pval_h PVAL --ld_h RSQR --addout test_run_1
```
## Running on single file without estimated linkage disequilibrium
###### You need to set paths for following in the **"gqs_config"**. You can edit the example of congif file in the repository 
1. PLINK[^1] is need to be installed
2. Path for the reference genome *(Note: the current script looks for pattern **chr** for example if --chrm is give 7 argument the script will look for 1000KG_chr7_impute.[fam,bim,bed])*

```
../GQSmain.py --ifile covidGWAS_chr7_test_file1.2c.txt --r2_th 0.01 --chrm 7 --snp_h SNP --pval_h PVAL --refG genome1000-EUR --addout test_run_2
```

## Running on summary stats
###### This requires locus region file (--regs). The GQS will be calcualted for each region defined in the rows of this file. The file should be space/tab delimited file without a header. 
- First column with Chromosome number.
- Second column with Start postion.
- Third column with End postion.

Example should of locus region file
```
1  43979146   44508946
2  200654071  201358071
2  57892490   58879490
3  135693422  136721422
3  161578675  161889805
3  180475138  181256138
3  48493599   52084599
4  176668514  176789614
4  24194810   24330510
6  152997243  153102881
```

```
../GQSmain.py --ifile small_daner --regs regions --pos_h BP --refG genome1000-EUR --snp_h SNP --pval_h P --chrm_h CHR --addout test_run_3
```


[^1]: https://zzz.bwh.harvard.edu/plink/
