RS = '/usr/bin/Rscript'
WD = /root/work
SD = /root/src

all: $(WD)/state-hogBalanceSheetDec2000Dec2001.csv $(WD)/shipment-flows-origins-on-rows-dests-on-columns.csv $(WD)/PEDvweeklyreport-state-ts-01-08-14.csv

$(WD)/state-hogBalanceSheetDec2000Dec2001.csv: $(SD)/get-balance-sheet-csv.sh
	cd $(WD) && $<

$(WD)/shipment-flows-origins-on-rows-dests-on-columns.csv: $(SD)/clean-shipment-matrix.r
	cd $(WD) && $<

$(WD)/PEDvweeklyreport-state-ts-01-08-14.csv: $(SD)/ped-weekly-report-scraper-py3.py
	cd $(WD) && $<
