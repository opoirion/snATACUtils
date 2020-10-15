library('edgeR')
library("data.table");
library("Matrix");
library("stringr")

library("argparse")

parser <- ArgumentParser()

parser$add_argument("-x", "--xgi", required=TRUE, help="barcode ID file")
parser$add_argument("-y", "--ygi", required=TRUE, help="feature name file")
parser$add_argument("-m", "--mfile", required=TRUE, help="matrix file in COO format")
parser$add_argument("-b", "--batch_file", required=TRUE, help="barcode -> batch ID file")
parser$add_argument("-c", "--case_control_file", required=TRUE, help="barcode -> group ID file")
parser$add_argument("-xs", "--xgi_subset", required=FALSE, help="barcode subset to use", default=NULL)
parser$add_argument("-ys", "--ygi_subset", required=FALSE, help="feature subset to use", default=NULL)
parser$add_argument("-ysn", "--ygi_subset_norm", required=FALSE, help="feature subset to use after edgeR normalisation", default=NULL)
parser$add_argument("-o", "--output", required=TRUE, help="output files name")
parser$add_argument("-uf", "--useFDR", required=FALSE, help="use FDR correction", default = FALSE)
parser$add_argument("-l", "--level", required=FALSE, help="level of significance", default = 0.05)
parser$add_argument("-he", "--heatm", required=FALSE, help="create HEATMAP (false|batch|group)", default = FALSE)
parser$add_argument("-ac", "--additional_cofactors", required=FALSE, help="Additional cofactors files. The header specifies datasetID and the different cofactors written as new columns", default = NULL)
parser$add_argument("-lfc", "--log_fold_change", required=FALSE, help="log fold change CUTOFF", default = 0.0)


#### example of cofactor file ####
# multiple cofactors can be added
## datasetID       sex
## CARE181125_3B   M
## CARE181125_3C   M
## CARE181125_3D   M
## CARE181213_2A   F
## CARE181213_2B   F
## CARE190307_10A  M
## CARE190307_10B  M
## CARE190307_10C  M
## CARE190307_10D  M
####################################



################## DEBUG #######################################################################################
## xgi = "/home/opoirion/data/lung_analysis/V2/edgeR_LUNG_justin/barcodes.tsv"
## ygi = "/home/opoirion/data/lung_analysis/V2/edgeR_LUNG_justin/features.tsv"
## mfile = "/home/opoirion/data/lung_analysis/V2/edgeR_LUNG_justin/counts_table.coo.gz"
## batch_file = "/home/opoirion/data/lung_analysis/V2/edgeR_LUNG_justin/group.tsv"
## case_control_file = "/home/opoirion/data/lung_analysis/V2/edgeR_LUNG_justin/timepoints/3_year_vs_30_year/airway_smooth_muscle.xgi"
## xgi_subset = "/home/opoirion/data/lung_analysis/V2/edgeR_LUNG_justin/timepoints/3_year_vs_30_year/airway_smooth_muscle.xgi"
## ygi_subset = NULL
## outFile = "/home/oliver/data/eye/DA_results_v2/DA.Cones"
## useFDR = FALSE
## level = "1e-5"
## doHEAT = TRUE
## heatType = "both"
## logFCCutoff = 0.0
################################################################################################################

args <- parser$parse_args()

xgi = args$xgi
ygi = args$ygi
mfile = args$mfile
batch_file = args$batch_file
case_control_file = args$case_control_file
xgi_subset = args$xgi_subset
ygi_subset = args$ygi_subset
ygi_subset_norm = args$ygi_subset_norm
outFile = args$output
useFDR = args$useFDR
level = as.numeric(args$level)
heatType = args$heatm
additional_cofactors = args$additional_cofactors
logFCCutoff = as.numeric(args$log_fold_change)


print(paste("#### LOG FC: ", logFCCutoff))

if (useFDR == "TRUE" || useFDR == "true" || useFDR == "True") {
    useFDR = TRUE
}

print(paste("#### HEATTYPE: ", heatType))

if (heatType != FALSE && !is.null(heatType) && heatType != "false" && heatType != "FALSE" && heatType != "False") {
    doHEAT = TRUE
} else {
    doHEAT = FALSE
}

if (substr(outFile, 1,1) != "/") {
    outFile = paste("./", outFile, sep="")
}

options(digit=14)

print("Loading Input file...")

print("Loading ygi...")
# Loading Feature names
con = file(ygi, 'r')
features = readLines(con)
features = str_trim(features)

print("Loading batch ID...")
# Loading batch data from a file having in each line the structure: barcodeID\tbatchID
grdata = read.table(file = batch_file, sep="\t", col.names = c("barcode", "batch"),
                    colClasses=c("character", "character"))

grdata$barcode = str_trim(grdata$barcode)
grdata$batch = str_trim(grdata$batch)

print("Loading case-control...")
# Loading the case-control data from a file having in each line the structure: barcodeID\caseID
ccdata = read.table(file = case_control_file, sep="\t", col.names = c("barcode", "case"),
                    colClasses=c("character", "character"))

ccdata$barcode = str_trim(ccdata$barcode)
ccdata$case = str_trim(ccdata$case)

print("Loading xgi...")
# Loading samples
xgi.data = read.table(file = xgi, stringsAsFactors = FALSE)
xgi.data = as.vector(xgi.data$V1)
xgi.data = str_trim(xgi.data)

print("Loading matrix...")
# Loading matrix
matrix = fread(mfile, sep="\t")

matrix = sparseMatrix(i = matrix$V1 + 1, j = matrix$V2 + 1, x = matrix$V3)

# Subsetting samples if provided
if (!is.null(xgi_subset)) {
    print("Subsetting xgi...")
    xgi.subset = read.table(file = xgi_subset, stringsAsFactors = FALSE)
    xgi.subset = as.vector(xgi.subset$V1)

    index = match(xgi.subset, xgi.data)

    matrix = matrix[index, ]
    xgi.data = xgi.subset

    grdata = grdata[match(xgi.subset, grdata$barcode),]
    ccdata = ccdata[match(xgi.subset, ccdata$barcode),]
}

# Subsetting features
if (!is.null(ygi_subset)) {
    print("Subsetting ygi...")
    con = file(ygi_subset, 'r')
    features.subset = readLines(con)

    index = match(features.subset, features)
    matrix = matrix[, index]
    features = features[index]
}

cofactors.table = NULL

# Loading additional cofactors
if (!is.null(additional_cofactors)) {
    cofactors.table  = read.table(file = additional_cofactors, stringsAsFactors = FALSE, header = TRUE)
    rownames(cofactors.table) = cofactors.table$datasetID
}


# Initialize the count matrix
countTable = matrix(nrow = dim(matrix)[2], ncol=0)

group = c()
batchList = c()
caseList = c()
count = 0
cofactors = list()

print("Create feature sum table...")
# Creating a sum vector for each batch and each case / control
for (batch in unique(grdata$batch)) {

    index = which(grdata$batch == batch)
    barcodesBatch = grdata$barcode[index]

    count = count + 1

    for (case in unique(ccdata$case)) {

        index2 = which(ccdata$case == case)

        barcodesCase = ccdata$barcode[index2]

        barcodes = intersect(barcodesCase, barcodesBatch)

        if (length(barcodes) < 2) {
            next
        }

        caseList = c(caseList, case)
        batchList = c(batchList, batch)

        indexMat = match(barcodes, xgi.data)
        countTable = cbind(countTable, colSums(matrix[indexMat, ]))

        group = c(group, case)

        if (!is.null(cofactors.table)) {
            for (cofactor in colnames(cofactors.table)[2: length(cofactors.table)]) {
                cofactors[[cofactor]] = append(cofactors[[cofactor]], cofactors.table[batch, cofactor])
            }
        }

        }
}

rm(matrix)
gc()

print(paste("Dimension of the EdgeR table: ", dim(countTable)))

print("EdgeR analysis...")
# Creating edgeR object
dgelist = DGEList(countTable, group=group)

# Naming row
if (length(batchList) != count) {
    batchList = paste(batchList, caseList, sep = "_")
}

colnames(dgelist) = batchList
rownames(dgelist) = features

# Calculate normalisation factor
dgelist = calcNormFactors(dgelist)

# Filter low expressed peaks
index = rowSums(cpm(dgelist$counts) > 1) >= 2
dgelist = dgelist[index, ]

# declare design of  experiment

if (is.null(cofactors.table)) {
    design = model.matrix(~group)
    factor = colnames(design)[dim(design)[2]]
} else {
    cofactors$group = group
    design = model.matrix(~ ., cofactors)
    factor = "group"

}

rownames(design) = colnames(dgelist)

print("Fitting...")
# Estimate dispertion
dgeDisp = estimateDisp(dgelist, design)

#likelihood ratio test
fittedDge = glmFit(dgeDisp, design)
fittedDge = glmLRT(fittedDge, coef=factor)

# Subsetting features after norm
if (!is.null(ygi_subset_norm)) {
    print("Subsetting normed features...")
    con = file(ygi_subset_norm, 'r')
    features.subset = readLines(con)

    index = match(features.subset, rownames(dgelist))
    index = index[!is.na(index)]
    dgeDisp = dgeDisp[index, ]
    dgelist = dgelist[index, ]
    fittedDge = fittedDge[index, ]
}

fittedDge$counts = dgelist$counts

indexSorted = order(fittedDge$table$PValue)

fittedDge = fittedDge[indexSorted, ]

# Adjust p-value
fittedDge$table$FDR = p.adjust(fittedDge$table$PValue, method = "BH")

print("Writting outputs...")
# Prepare output table
outputTable = cbind(fittedDge$counts, fittedDge$table)
rownames(outputTable) = rownames(fittedDge$counts)

# Create output and bed files
dir.create(dirname(outFile), showWarnings = FALSE)

write.table(outputTable,
            file = paste(outFile, ".tsv", sep=""),
            sep = "\t", col.names = NA)


if (isTRUE(useFDR)) {
    caseBed = rownames(outputTable)[
        which(outputTable$logFC < logFCCutoff & outputTable$FDR<level)
    ]
    controlBed = rownames(outputTable)[
        which(outputTable$logFC > logFCCutoff & outputTable$FDR<level)
    ]
} else {
    caseBed = rownames(outputTable)[
        which(outputTable$logFC < logFCCutoff & outputTable$PValue<level)
    ]
    controlBed = rownames(outputTable)[
        which(outputTable$logFC > logFCCutoff & outputTable$PValue<level)
    ]

}

write.table(caseBed, file = paste(outFile, ".case.bed", sep=""),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
write.table(controlBed, file = paste(outFile, ".control.bed", sep=""),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

print(paste("output files :", outFile))

if (isTRUE(doHEAT)) {
    library('heatmap3')
    library('RColorBrewer')

    heatm = matrix(nrow = dim(fittedDge$counts)[1], ncol=length(unique(group)))
    rownames(heatm) = rownames(fittedDge$counts)
    colnames(heatm) = unique(group)

    index = c()

    for (gr in unique(group)) {
        heatm[, gr] = rowSums(fittedDge$counts[, which(group == gr)])
    }

    if (isTRUE(useFDR)) {
        index = which(fittedDge$table$FDR < level)
    } else {
        index = which(fittedDge$table$PValue < level)

    }
    heatm = heatm[, order(colnames(heatm))]

    index = order(fittedDge$table$logFC[index])
    heatm = heatm[index, ]

    if (isTRUE(heatType != "donor")) {

        pngName = paste(outFile, ".batch.heatmap.png", sep = "")

        png(file=pngName)
        heatmap3(log(1.0 + heatm),
                 Colv = NA,
                 Rowv = NA,
                 scale = 'none',
                 labRow = FALSE, margins = c(5,2),
                 col = colorRampPalette(brewer.pal(8, "Blues"))(25),
                 main = sprintf("Nb features: %d", dim(heatm)[1]),
                 cexCol = 1.0)
        dev.off()

        print(paste("heatmap png saved:", pngName))

    }

    if (!is.na(match(heatType, c("donor", "both")))) {
        heatm = fittedDge$counts
        pngName = paste(outFile, ".donor.heatmap.png", sep = "")

        heatm = heatm[, order(colnames(heatm))]
        heatm = heatm[index, ]

        png(file=pngName)
        heatmap3(log(1.0 + heatm),
                 Colv = NA,
                 Rowv = NA,
                 scale = 'none',
                 labRow = FALSE, margins = c(10,2),
                 col = colorRampPalette(brewer.pal(8, "Blues"))(25),
                 main = sprintf("Nb features: %d", dim(heatm)[1]),
                 cexCol = 1.0)
        dev.off()

        print(paste("heatmap png saved:", pngName))

    }
}
