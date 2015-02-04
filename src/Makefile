RS = '/usr/bin/Rscript'
WD = /root/work
SD = /root/src

$(WD)/shipment-flows-origins-on-rows-dests-on-columns.csv: $(SD)/clean-shipment-matrix.r
	cd $(WD) && $(RS) $< 
