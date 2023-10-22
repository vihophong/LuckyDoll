#!/bin/bash                                                                                                                  
ls -tr $1/R$2_* | wc -l > $3
ls -tr $1/R$2_* >> $3
