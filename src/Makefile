RS = '/usr/bin/Rscript'
WD = /root/work
SD = /root/src

all: $(WD)/flows-prediction.eps

$(WD)/flows-prediction.eps: $(SD)/make-regression-results-figs.r $(WD)/flows.html
	cd $(WD) && $< >$<.so 2>$<.se 

$(WD)/flows.html: $(WD)/state-hogBalanceSheetDec2000Dec2001.csv $(WD)/shipment-flows-origins-on-rows-dests-on-columns.csv $(WD)/PEDvweeklyreport-state-ts-01-08-14.csv | data
	cd $(WD) && cp $(SD)/flows.r . && $(RS) --vanilla -e "library(knitr); spin(\"flows.r\")" >flows.r.so 2>flows.r.se

$(WD)/state-hogBalanceSheetDec2000Dec2001.csv: $(SD)/get-balance-sheet-csv.sh
	cd $(WD) && $<

$(WD)/shipment-flows-origins-on-rows-dests-on-columns.csv: $(SD)/clean-shipment-matrix.r
	cd $(WD) && $<

$(WD)/PEDvweeklyreport-state-ts-01-08-14.csv: $(SD)/ped-weekly-report-scraper-py3.py
	cd $(WD) && $<

data:
	cd $(WD) && cp -p ../data/* .
