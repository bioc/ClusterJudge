library(ClusterJudge)

test_that('error when clusters is an unamed vector or entity attribute does not have 2 columns',{
  expect_error(clusterJudge(c(1,2),c(1,2)),'clusters must have as names the entity names of the reference entity.attribute')

  expect_error(clusterJudge(c('1'=1,'2'=2),c(1,2)),'entity attribute structure must have only 2 columns')

}
)