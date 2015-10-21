context("agent generation")

test_that("list of neighboring raster cells is correct", {
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
