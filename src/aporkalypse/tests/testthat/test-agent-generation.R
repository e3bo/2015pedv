context("agent generation")

test_that("List of neighboring raster cells is correct", {
            r <- raster::raster(county.hogs.pigs.02.map)
            raster::res(r) <- 16000
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

#test_that("Cell to county mappings are correct", {
#
#          })

test_that("Farm density map looks like other maps", {
            r <- raster::raster(county.hogs.pigs.02.map)
            raster::res(r) <- 16000
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
            foo2 <- CreateAgents(census.dilation=1)
            bar2 <- raster::rasterize(foo2$adf[, c('x', 'y')], y=r, fun='count')
            v1 <- raster::values(target.raster)
            v2 <- raster::values(bar2)
            v1 <- ifelse(is.na(v1), 0, v1)
            v2 <- ifelse(is.na(v2), 0, v2)
            rmsd <- sd(v1 - v2)
            expect_true(rmsd < 5)
        })
