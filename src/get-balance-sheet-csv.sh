#!/bin/bash

set -e

file1=$(mktemp)
file2=$(mktemp)
outfile='state-hogBalanceSheetDec2000Dec2001.csv'
archive='/home/docker/data/usda.mannlib.cornell.edu_usda_nass_MeatAnimPr__2000s_2002_MeatAnimPr_04_26_2002.zip'

echo 'state,inventory2000,pigCrop,inshipments' > "${file1}" 

unzip -c "${archive}" meat_014.csv \
  | grep -a ,\"d\" \
  | cut -d',' -f3- \
  | grep -v ^\"\" \
  | grep -v ^\"Idaho \
  | grep -v ^\"Washington \
  | sed 's/ Washington/IdahoWashington/' >> "${file1}"

echo 'state2,marketings,slaughter,deaths,inventory2001' > "${file2}"
unzip -c "${archive}" meat_015.csv \
  | grep -a ,\"d\" \
  | cut -d',' -f3- \
  | grep -v ^\"\" \
  | grep -v ^\"Idaho \
  | grep -v ^\"Washington \
  | sed 's/ Washington/IdahoWashington/' >> "${file2}"

paste --delimiter="," "${file1}" "${file2}" >"${outfile}"
sed -i -e 's/\r//g' "${outfile}"
rm "${file1}" "${file2}"
