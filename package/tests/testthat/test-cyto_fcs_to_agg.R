context("bla")

test_that("fcs_to_agg works", {
  fcs_dir <- "../../test_dir/"
  if (!dir.exists(fcs_dir)) fcs_dir <- "package/test_dir/"
  fcs_names <- list.files(path = fcs_dir, pattern = ".fcs")
  seed = 1
  cTotal = 30000
  output_name = paste0(fcs_dir, "/test_agg.fcs")

  ff_agg <- FlowSOM::AggregateFlowFrames(paste0(fcs_dir, "/", fcs_names),
                                         cTotal= cTotal, writeOutput = TRUE,
                                         outputFile = output_name)

  if (file.exists(output_name)) file.remove(output_name)
  ff_agg <- flowCore::transform(ff_agg,
                                flowCore::transformList(colnames(ff_agg@exprs)[c(3,17,28:62,71)], flowCore::arcsinhTransform(b=1/5, a=0, c=0)))

  expect_is(ff_agg, "flowFrame")
  expect_true(class(ff_agg) == "flowFrame")
})
