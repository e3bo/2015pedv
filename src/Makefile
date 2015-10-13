RS = '/usr/bin/Rscript'
WD = /root/work
SD = /root/src
DD = /root/data
VPATH = $(SD):$(WD)

all: accessions-premises-correlation.r.so flagMantelPlots flagAssocPlots flagTsirFigs sim-panel.pdf mapPlotSeries

mapPlotSeries: map-series.r flagData
	$< && touch $@

accessions-premises-correlation.r.so: accessions-premises-correlation.r flagData
	$< >$(<F).so 2>$(<F).se

sim-panel.pdf: cor-sim-study-plots.r flagSims1
	$< >$(<F).so 2>$(<F).se

flagSims1: corSims-layer1.r flagSims0
	$< >$(<F).so 2>$(<F).se && touch $@

flagSims0: corSims-layer0.r flagSimulatorInstalled
	$< >$(<F).so 2>$(<F).se && touch $@

flagSimulatorInstalled: $(shell find $(SD)/duoDynamics) .BatchJobs.R
	R CMD build $(SD)/duoDynamics && R CMD INSTALL duoDynamics_1.0.tar.gz && cp $(SD)/.BatchJobs.R . && touch $@

flagMantelPlots: mantel-plots.r custom-plot-functions.r flagMantel
	cp $(SD)/custom-plot-functions.r . && $< >$(<F).so 2>$(<F).se && touch $@

flagMantel: mantel-testing.r flagData flagWeeklyReports flagPedvCum
	$< >$(<F).so 2>$(<F).se && touch $@

flagAssocPlots: make-assoc-plots.r flagPedvCum
	cp $< . && $(RS) -e "library(knitr); spin(\"make-assoc-plots.r\")" >$(<F).so 2>$(<F).se && touch $@

flagPedvCum: pedv-cum.r flagWeeklyReports shipment-flows-origins-on-rows-dests-on-columns.csv state-hogBalanceSheetDec2011Dec2012.csv flagData
	cp $< . && $(RS) -e "library(knitr); spin(\"pedv-cum.r\")" >$(<F).so 2>$(<F).se && touch $@

flagTsirFigs: make-regression-results-figs.r flows.html
	$< >$(<F).so 2>$(<F).se && touch $@

flows.html: flows.r state-hogBalanceSheetDec2000Dec2001.csv shipment-flows-origins-on-rows-dests-on-columns.csv flagWeeklyReports flagData
	cp $< . && $(RS) -e "library(knitr); spin(\"flows.r\")" >$(<F).so 2>$(<F).se

state-hogBalanceSheetDec2000Dec2001.csv: get-balance-sheet-csv.sh
	$<

state-hogBalanceSheetDec2011Dec2012.csv: get-balance-sheet-2012-csv.sh
	$<

shipment-flows-origins-on-rows-dests-on-columns.csv: clean-shipment-matrix.r
	$<

flagWeeklyReports: ped-weekly-report-scraper-py3.py
	$< && touch $@

flagData:
	cp -p $(DD)/* . && ln -s $(DD)/states-shapefile . && touch $@
