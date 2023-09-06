#!/bin/bash
#count up matching alignments to a_sym_, b_sym_, c_sym_, d_sym_
awk 'BEGIN {OFS = "\t"
sym=0;all=0;nonSym=0}
{all++
if($3 ~ "scaffold")
	sym++
else
	nonSym++}
END {print all,nonSym,symna}'

