context("agent generation")

r <- raster::raster(county.hogs.pigs.02.map)
raster::res(r) <- 16000

test_that("List of neighboring raster cells is correct", {
            nc <- raster::ncell(r)
            ncol <- raster::ncol(r)
            cell.samps <- c(1, ncol, nc - ncol + 1, nc, sample.int(nc, size=100))
            for(cell in cell.samps){
              nbs <- Get8nbs(r, cell)
              expect_true(cell %in% nbs)
              nadj <- nrow(raster::adjacent(r, nbs, target=cell, directions=8))
              expect_equal(nadj, length(nbs) - 1)
            }
          })


agents <- CreateAgents(census.dilation=1)

test_that("Cell to county mappings are correct", {
            nr <- nrow(agents$adf)
            sampled.ids <- sample.int(n=nr, size=100)
            for(id in sampled.ids){
              coord.cell <- raster::cellFromXY(r, agents$adf[id, c('x', 'y')])
              expect_equal(agents$adf$cell[id], coord.cell)
              poly <- GetCountySPDF(agents$adf$stfips[id],
                                    agents$adf$cofips[id])
              pt <- SpatialPoints(agents$adf[id, c('x', 'y')],
                                  proj4string=sp::CRS(sp::proj4string(poly)))
              over.out <- over(pt, poly)
              expect_equal(nrow(over.out), 1)
            }
          })

test_that("Farm density map looks reasonable", {
            if(FALSE){
              foo <- CreateAgents(census.dilation=1)
              bar <- raster::rasterize(foo$adf[, c('x', 'y')], y=r, fun='count')
              png('farms-per-16000-km2-cell.png', width=1000, height=800) ## For sanity check
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
