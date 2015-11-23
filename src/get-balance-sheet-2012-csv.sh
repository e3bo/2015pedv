#!/bin/bash

set -e

file1=$(mktemp)
file2=$(mktemp)
outfile='state-hogBalanceSheetDec2011Dec2012.csv'
archive='/home/docker/data/usda01.library.cornell.edu_usda_current_MeatAnimPr_MeatAnimPr-04-25-2013.zip'

echo 'state,inventory2011,pigCrop,inshipments' > "${file1}" 

unzip -c "${archive}" meat_p14_t014.csv \
  | grep ,\"d\" \
  | cut -d',' -f3- \
  | grep -v ^\"\" \
  | grep -v ^\"Idaho \
  | grep -v ^\"Washington \
  | sed 's/ Washington/IdahoWashington/' >> "${file1}"

echo 'state2,marketings,slaughter,deaths,inventory2012' > "${file2}"
unzip -c "${archive}" meat_p15_t015.csv \
  | grep ,\"d\" \
  | cut -d',' -f3- \
  | grep -v ^\"\" \
  | grep -v ^\"Idaho \
  | grep -v ^\"Washington \
  | sed 's/ Washington/IdahoWashington/' >> "${file2}"

paste --delimiter="," "${file1}" "${file2}" >"${outfile}"
sed -i -e 's/\r//g' "${outfile}"
rm "${file1}" "${file2}"
