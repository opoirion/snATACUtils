library("Matrix")
library("data.table")

library("argparse")

parser <- ArgumentParser()


################################ Collecting input args ##############################
parser$add_argument("-x", "--xgi", required=TRUE, help="xgi file")
parser$add_argument("-y", "--ygi", required=TRUE, help="ygi file")

parser$add_argument("-f", "--coo_matrix", required=TRUE, help="COO matrix file")

parser$add_argument("-o", "--out", required=TRUE, help="Output_file")

args <- parser$parse_args()
####################################################################################

xgi = args$xgi
ygi = args$ygi
mfile = args$coo_matrix
output_file = args$out
############################### Creating matrix ####################################
timestart <- Sys.time()
peak_ids = fread(
    ygi,
    sep="\t", header = FALSE)
cell_ids = fread(
    xgi,
    sep="\t", header = FALSE)

dat = fread(mfile, sep="\t")

mat = sparseMatrix(i = dat$V1 + 1, j =dat$V2 + 1, x=dat$V3)
print(sprintf("Matrix created in %f s", Sys.time() - timestart))
####################################################################################

################################ Saving matrix ######################################
# Saving using writeMM. Can be loaded using loadMM
timestart <- Sys.time()
writeMM(mat, file = output_file)
print(sprintf("Matrix saved in %f s", Sys.time() - timestart))
######################################################################################
