#!/bin/bash

cython -3 --embed bwt_based.py -o bwt_based.c

gcc 							\
    -O3 						\
    bwt_based.c						\
    -I /usr/include/python3.7 	\
    -L /home/andy/miniconda3/envs/gsa_p5/include/python3.7m	\
    -L /home/andy/miniconda3/envs/gsa_p5/lib/libpython3.7m.a 			\
    -l python3.7.a 					\
    -l pthread 						\
    -l m 						\
    -l util 						\
    -l dl 						\
    -o bwt_based_andy
