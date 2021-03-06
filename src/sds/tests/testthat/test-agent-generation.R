context("agent generation")

data('county.hogs.pigs.02.map', package='sds')
r <- raster::raster(county.hogs.pigs.02.map, ncol=285, nrow=178)
set.seed(123)

test_that("List of neighboring raster cells is correct", {
            ncell <- raster::ncell(r)
            ncol <- raster::ncol(r)
            nrow <- raster::nrow(r)
            cell.samps <- c(1, ncol, ncell - ncol + 1, ncell, sample.int(ncell, size=100))
            for(cell in cell.samps){
              nbs <- Get8nbs(r=r, cell=cell, ncells=ncell, nc=ncol, nr=nrow)
              expect_true(cell %in% nbs)
              nadj <- nrow(raster::adjacent(r, nbs, target=cell, directions=8))
              expect_equal(nadj, length(nbs) - 1)
            }
          })

test_that("Formating of fips codes is correct", {
            expect_identical(FormatFips(1, 'state'), '01')
            expect_identical(FormatFips(33, 'state'), '33')
            expect_identical(FormatFips(1, 'county'), '001')
            expect_identical(FormatFips(200, 'county'), '200')
        })

agents <- CreateAgents(census.dilation=1)

test_that("Network consistent with flow estimates", {
            types <- unique(as.character(agents$adf$abb))
            targetM <- SubsetFlows(types)
            M <- targetM
            M[TRUE] <- 0
            agents$adf <- agents$adf[order(agents$adf$id), ]
                                        # Assumed for looking up neighbors
            for(i in 1:nrow(agents$adf)){
              type.i <- as.character(agents$adf$abb[i])
              for(j in agents$net.nbs[[i]]){
                type.j <- as.character(agents$adf$abb[j])
                M[type.i, type.j] <- M[type.i, type.j] + 1
              }
            }
            cor.M.targetM <- cor(as.numeric(M), as.numeric(targetM))
            expect_more_than(cor.M.targetM, 0.99)
          })

test_that("Cell to county mappings are correct", {
            nr <- nrow(agents$adf)
            sampled.ids <- sample.int(n=nr, size=100)
            for(id in sampled.ids){
              coord.cell <- raster::cellFromXY(r, agents$adf[id, c('x', 'y')])
              expect_equal(agents$adf$cell[id], coord.cell)
              poly <- GetCountySPDF(agents$adf$stfips[id],
                                    agents$adf$cofips[id],
                                    county.hogs.pigs.02.map)
              pt <- sp::SpatialPoints(agents$adf[id, c('x', 'y')],
                                      proj4string=sp::CRS(sp::proj4string(poly)))
              over.out <- sp::over(pt, poly)
              expect_equal(nrow(over.out), 1)
            }
          })

test_that("Farm density map looks reasonable", {
          if(FALSE){
              foo <- CreateAgents(census.dilation=1)
              bar <- raster::rasterize(foo$adf[, c('x', 'y')], y=r, fun='count')
              png('farms-per-approx-16000-km-side-cell.png', width=1000, height=800) ## For sanity check
              sp::plot(county.hogs.pigs.02.map, col='grey')
              sp::plot(bar, add=TRUE)
              dev.off()
              target.raster <- bar
              save(target.raster, file='sysdata.rda')
          }
            current.raster <- raster::rasterize(agents$adf[, c('x', 'y')],
                                              y=r, fun='count')
            v1 <- raster::values(target.raster)
            v2 <- raster::values(current.raster)
            v1 <- ifelse(is.na(v1), 0, v1)
            v2 <- ifelse(is.na(v2), 0, v2)
            rmsd <- sd(v1 - v2)
            expect_true(rmsd < 5)
        })

test_that("Simulation output is sane", {
            nsteps <- 38
            adf.out <- RunSim(adf=agents$adf, net.nbs=agents$net.nbs, sp.nbs=agents$sp.nbs,
                              tprob.net=1, tprob.sp=1, rprob=1, nsteps=nsteps)
            expect_true(all(adf.out$infection.time < adf.out$recovery.time, na.rm=TRUE))
            tsl <- GetTimeSeries(adf.out, nsteps=nsteps)
            expect_true(all(rowSums(tsl$new.cases) == rowSums(tsl$no.infected)))

            cases <- sample.int(n=nrow(agents$adf), size=100)
            adf.out <- RunSim(adf=agents$adf, net.nbs=agents$net.nbs, sp.nbs=agents$sp.nbs,
                              tprob.net=1, tprob.sp=1, rprob=1, cases=cases, nsteps=nsteps)
            expect_true(all(adf.out$infection.time < adf.out$recovery.time, na.rm=TRUE))
            expect_true(all(adf.out$infection.time[cases] == 0))
            tsl <- GetTimeSeries(adf.out, nsteps=nsteps)
            expect_true(all(rowSums(tsl$new.cases) == rowSums(tsl$no.infected)))

            cases <- sample.int(n=nrow(agents$adf), size=100)
            adf.out <- RunSim(adf=agents$adf, net.nbs=agents$net.nbs, sp.nbs=agents$sp.nbs,
                              tprob.net=.1, tprob.sp=.1, rprob=0.5, cases=cases, seasonal.amplitude=1,
                              nsteps=nsteps)
            ip <- adf.out$recovery.time - adf.out$infection.time
            tt <- t.test(ip[adf.out$infection.time == 0], mu=2)
            expect_true(tt$p.value > 1e-4)

            tt2 <- t.test(ip[adf.out$infection.time != 0],
                          ip[adf.out$infection.time == 0], alternative='greater')
            expect_true(tt$p.value > 0.05)

            abs(tt$null.value - tt$estimate)
            expect_true(all(adf.out$infection.time < adf.out$recovery.time, na.rm=TRUE))
            expect_true(all(adf.out$infection.time[cases] == 0))
            tsl <- GetTimeSeries(adf.out, nsteps=nsteps)
            expect_true(all(rowSums(tsl$new.cases) <= rowSums(tsl$no.infected) | rowSums(tsl$new.cases ==0)))
       })
