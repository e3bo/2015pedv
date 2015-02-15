RS = '/usr/bin/Rscript'
WD = /root/work
SD = /root/src
VPATH = $(SD):$(WD)

all: sim-panel.pdf flagMantel

sim-panel.pdf: cor-sim-study-plots.r corSims-layer1.RData
	cd $(WD) && $< >$(<F).so 2>$(<F).se

corSims-layer1.RData: corSims-layer1.r corSims-layer0.RData
	cd $(WD) && $< >$(<F).so 2>$(<F).se

corSims-layer0.RData: corSims-layer0.r flagSimulatorInstalled
	cd $(WD) && $< >$(<F).so 2>$(<F).se

flagSimulatorInstalled: $(shell find /root/src/duoDynamics)
	cd $(WD) && R CMD build $(SD)/duoDynamics && R CMD INSTALL duoDynamics_1.0.tar.gz && cp $(SD)/.BatchJobs.R $(WD)/ && touch $@

flagMantel: mantel-testing.r custom-plot-functions.r flagData flagWeeklyReports flagPedvCum
	cd $(WD) && cp $(SD)/custom-plot-functions.r . && $< >$(<F).so 2>$(<F).se && touch $@

flagAssocPlots: make-assoc-plots.r flagPedvCum
	cd $(WD) && cp $< . && $(RS) -e "library(knitr); spin(\"make-assoc-plots.r\")" >$(<F).so 2>$(<F).se && touch $@

flagPedvCum: pedv-cum.r flagWeeklyReports shipment-flows-origins-on-rows-dests-on-columns.csv state-hogBalanceSheetDec2011Dec2012.csv flagData
	cd $(WD) && cp $< . && $(RS) -e "library(knitr); spin(\"pedv-cum.r\")" >$(<F).so 2>$(<F).se && touch $@

flagTsirFigs: make-regression-results-figs.r flows.html
	cd $(WD) && $< >$(<F).so 2>$(<F).se && touch $@

flows.html: flows.r state-hogBalanceSheetDec2000Dec2001.csv shipment-flows-origins-on-rows-dests-on-columns.csv flagWeeklyReports flagData
	cd $(WD) cp $< . && $(RS) -e "library(knitr); spin(\"flows.r\")" >$(<F).so 2>$(<F).se

state-hogBalanceSheetDec2000Dec2001.csv: get-balance-sheet-csv.sh
	cd $(WD) && $<

state-hogBalanceSheetDec2011Dec2012.csv: get-balance-sheet-2012-csv.sh
	cd $(WD) && $<

shipment-flows-origins-on-rows-dests-on-columns.csv: clean-shipment-matrix.r
	cd $(WD) && $<

flagWeeklyReports: ped-weekly-report-scraper-py3.py
	cd $(WD) && $< && touch $@

flagData:
	cd $(WD) && cp -p ../data/* . && touch $@
