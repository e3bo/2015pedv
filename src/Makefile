RS = '/usr/bin/Rscript'
WD = /root/work
SD = /root/src

all: $(WD)/sim-panel.pdf

$(WD)/sim-panel.pdf: $(SD)/cor-sim-study-plots.r $(WD)/corSims-layer1.RData
	cd $(WD) && $< >$(<F).so 2>$(<F).se

$(WD)/corSims-layer1.RData: $(SD)/corSims-layer1.r $(WD)/corSims-layer0.RData
	cd $(WD) && $< >$(<F).so 2>$(<F).se

$(WD)/flagMantel: $(SD)/mantel-testing.r $(SD)/custom-plot-functions.r $(WD)/flagData $(WD)/flagWeeklyReports $(WD)/flagPedvCum
	cd $(WD) && cp $< . && cp $(SD)/custom-plot-functions.r . && ./$(<F) >$(<F).so 2>$(<F).se && touch $@

$(WD)/flagAssocPlots: $(SD)/make-assoc-plots.r $(WD)/flagPedvCum
	cd $(WD) && cp $< . && $(RS) -e "library(knitr); spin(\"make-assoc-plots.r\")" >$(<F).so 2>$(<F).se && touch $@

$(WD)/flagPedvCum: $(SD)/pedv-cum.r $(WD)/flagWeeklyReports $(WD)/shipment-flows-origins-on-rows-dests-on-columns.csv $(WD)/state-hogBalanceSheetDec2011Dec2012.csv $(WD)/flagData
	cd $(WD) && cp $< . && $(RS) -e "library(knitr); spin(\"pedv-cum.r\")" >$(<F).so 2>$(<F).se && touch $@

$(WD)/flagTsirFigs: $(SD)/make-regression-results-figs.r $(WD)/flows.html
	cd $(WD) && cp $< . && ./$(<F) >$(<F).so 2>$(<F).se && touch $@

$(WD)/flows.html: $(SD)/flows.r $(WD)/state-hogBalanceSheetDec2000Dec2001.csv $(WD)/shipment-flows-origins-on-rows-dests-on-columns.csv $(WD)/flagWeeklyReports $(WD)/flagData
	cd $(WD) && cp $< . && $(RS) -e "library(knitr); spin(\"flows.r\")" >$(<F).so 2>$(<F).se

$(WD)/state-hogBalanceSheetDec2000Dec2001.csv: $(SD)/get-balance-sheet-csv.sh
	cd $(WD) && $<

$(WD)/state-hogBalanceSheetDec2011Dec2012.csv: $(SD)/get-balance-sheet-2012-csv.sh
	cd $(WD) && $<

$(WD)/shipment-flows-origins-on-rows-dests-on-columns.csv: $(SD)/clean-shipment-matrix.r
	cd $(WD) && $<

$(WD)/flagWeeklyReports: $(SD)/ped-weekly-report-scraper-py3.py
	cd $(WD) && $< && touch $@

$(WD)/flagData:
	cd $(WD) && cp -p ../data/* . && touch $@
