# TODO: Add comment
# 
# Author: tobias
###############################################################################



setClass("limma",
    representation = list(MArrayLM = "MArrayLM",
        geneSymbols = "character")
)

### S4 class for output of computeLogRatio function
setClass("ExpressionSetWithComputation",
		contains = "ExpressionSet")


