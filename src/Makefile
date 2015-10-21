RS = '/usr/bin/Rscript'
WD = /root/work
SD = /root/src
DD = /root/data
VPATH = $(SD):$(WD)

all: accessions-premises-correlation.r.so flagMantelPlots flagAssocPlots flagTsirFigs sim-panel.pdf mapPlotSeries

mapPlotSeries: map-series.r flagData flagWeeklyReports
	$< && touch $@

accessions-premises-correlation.r.so: accessions-premises-correlation.r flagData
	$< >$(<F).so 2>$(<F).se

sim-panel.pdf: cor-sim-study-plots.r flagSims1
	$< >$(<F).so 2>$(<F).se

flagSims1: corSims-layer1.r flagSims0
	$< >$(<F).so 2>$(<F).se && touch $@

flagSims0: corSims-layer0.r flagPackageInstalled
	$< >$(<F).so 2>$(<F).se && touch $@

flagPackageInstalled: $(shell find $(SD)/aporkalypse) .BatchJobs.R flagWeeklyReports shipment-flows-origins-on-rows-dests-on-columns.csv state-hogBalanceSheetDec2011Dec2012.csv  state-hogBalanceSheetDec2000Dec2001.csv flagData
	cp -r $(SD)/aporkalypse . && cp aporkalypse/data-raw/data-prep.R . && $(RS) -e "library(knitr); spin(\"data-prep.R\")"  && mv data aporkalypse/data && R CMD build ./aporkalypse && R CMD INSTALL aporkalypse_0.0.0.9000.tar.gz && cp $(SD)/.BatchJobs.R . && touch $@

flagMantelPlots: mantel-plots.r custom-plot-functions.r flagMantel
	cp $(SD)/custom-plot-functions.r . && $< >$(<F).so 2>$(<F).se && touch $@

flagMantel: mantel-testing.r flagData flagWeeklyReports flagPedvCum flagPackageInstalled
	$< >$(<F).so 2>$(<F).se && touch $@

flagAssocPlots: make-assoc-plots.r flagPedvCum
	cp $< . && $(RS) -e "library(knitr); spin(\"make-assoc-plots.r\")" >$(<F).so 2>$(<F).se && touch $@

flagPedvCum: pedv-cum.r flagPackageInstalled
	cp $< . && $(RS) -e "library(knitr); spin(\"pedv-cum.r\")" >$(<F).so 2>$(<F).se && touch $@

flagTsirFigs: make-regression-results-figs.r flows.html
	$< >$(<F).so 2>$(<F).se && touch $@

flows.html: flows.r flagPackageInstalled
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
	cp -pr $(DD)/* . && touch $@
