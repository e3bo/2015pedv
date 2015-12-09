RS = '/usr/bin/Rscript'
WD = /home/docker/work
SD = /home/docker/src
DD = /home/docker/data
VPATH = $(SD):$(WD)

simflags := flagSims-r-sharedBord flagSims-r-shipment flagSims-r-gcd
metamodels := sim-study-checkpoint4-r-sharedBord.rda \
sim-study-checkpoint4-r-shipment.rda \
sim-study-checkpoint4-r-gcd.rda

all: accessions-premises-correlation.r.so flagMantelPlots flagAssocPlots \
flagTsirFigs $(simflags) mapPlotSeries

mapPlotSeries: map-series.r flagData flagWeeklyReports
	$< && touch $@

accessions-premises-correlation.r.so: accessions-premises-correlation.r \
flagData
	$< >$(<F).so 2>$(<F).se

flagSims: $(simflags)

$(simflags): flagSims-%: sim-study-layer05.R sim-study-checkpoint4-%.rda
	vars=$* $< >$(<F)-$*.so 2>$(<F)-$*.se && touch $@

$(metamodels): sim-study-checkpoint4-%.rda: sim-study-layer04.R sim-study-checkpoint3.rda
	vars=$* $< >$(<F)-$*.so 2>$(<F)-$*.se

sim-study-checkpoint3.rda: sim-study-layer03.R sim-study-checkpoint2.rda
	$< >$(<F).so 2>$(<F).se

sim-study-checkpoint2.rda: sim-study-layer02.R sim-study-checkpoint1.rda
	$< >$(<F).so 2>$(<F).se

sim-study-checkpoint1.rda: sim-study-layer01.R flagPackageInstalled
	$< >$(<F).so 2>$(<F).se

flagPackageInstalled: $(shell find $(SD)/sds) flagWeeklyReports \
shipment-flows-origins-on-rows-dests-on-columns.csv \
state-hogBalanceSheetDec2011Dec2012.csv \
state-hogBalanceSheetDec2000Dec2001.csv \
flagData
	cp -r $(SD)/sds . && cp sds/data-raw/data-prep.R . \
	&& $(RS) -e "library(knitr); spin(\"data-prep.R\")" \
	&& mv data sds/data && R CMD build ./sds \
        && R CMD check --no-manual --no-vignettes sds_*.tar.gz >check.so \
	&& R CMD INSTALL sds_0.0.0.9000.tar.gz && touch $@

flagMantelPlots: mantel-plots.r custom-plot-functions.r flagMantel
	cp $(SD)/custom-plot-functions.r . && $< >$(<F).so 2>$(<F).se \
	&& touch $@

flagMantel: mantel-testing.r flagData flagWeeklyReports flagPedvCum \
	flagPackageInstalled
	$< >$(<F).so 2>$(<F).se && touch $@

flagAssocPlots: make-assoc-plots.r flagPedvCum
	cp $< . && $(RS) -e "library(knitr); spin(\"make-assoc-plots.r\")" \
	>$(<F).so 2>$(<F).se && touch $@

flagPedvCum: pedv-cum.r flagPackageInstalled
	cp $< . && $(RS) -e "library(knitr); spin(\"pedv-cum.r\")" \
	>$(<F).so 2>$(<F).se && touch $@

flagTsirFigs: make-regression-results-figs.r flows.html
	$< >$(<F).so 2>$(<F).se && touch $@

flows.html: flows.r flagPackageInstalled
	cp $< . && $(RS) -e "library(knitr); spin(\"flows.r\")" >$(<F).so \
	2>$(<F).se

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
