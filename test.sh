#!/bin/bash

./parseInputData test.vcf test.BB test.txt
R CMD BATCH --no-save --no-restore '--args test.txt test.pur' BayCloneC.R ./test.log
