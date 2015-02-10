RS = '/usr/bin/Rscript'
WD = /root/work
SD = /root/src

all: $(WD)/flagPedvCum

$(WD)/flagPedvCum: $(SD)/pedv-cum.r $(WD)/flagWeeklyReports $(WD)/shipment-flows-origins-on-rows-dests-on-columns.csv $(WD)/state-hogBalanceSheetDec2011Dec2012.csv | data
	cd $(WD) && cp $< . && $(RS) -e "library(knitr); spin(\"pedv-cum.r\")" >$(<F).so 2>$(<F).se && touch $@

$(WD)/flagTsirFigs: $(SD)/make-regression-results-figs.r $(WD)/flows.html
	cd $(WD) && $< >$(<F).so 2>$(<F).se  && touch $@

$(WD)/flows.html: $(SD)/flows.r $(WD)/state-hogBalanceSheetDec2000Dec2001.csv $(WD)/shipment-flows-origins-on-rows-dests-on-columns.csv $(WD)/flagWeeklyReports | data
	cd $(WD) && cp $< . && $(RS) -e "library(knitr); spin(\"flows.r\")" >$(<F).so 2>$(<F).se

$(WD)/state-hogBalanceSheetDec2000Dec2001.csv: $(SD)/get-balance-sheet-csv.sh
	cd $(WD) && $<

$(WD)/state-hogBalanceSheetDec2011Dec2012.csv: $(SD)/get-balance-sheet-2012-csv.sh
	cd $(WD) && $<


$(WD)/shipment-flows-origins-on-rows-dests-on-columns.csv: $(SD)/clean-shipment-matrix.r
	cd $(WD) && $<

$(WD)/flagWeeklyReports: $(SD)/ped-weekly-report-scraper-py3.py
	cd $(WD) && $< && touch $@

data:
	cd $(WD) && cp -p ../data/* .
