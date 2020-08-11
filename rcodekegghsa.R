#args[1] - file with the data table
#args[2] - header present (YES/NO)
#args[3] - number of column with gene IDs (Gene symbols or entrez IDs - autodetect)
#args[4] - number of column with gene measure
#args[5] - Statistics type (Parametric / Non-parametric)

#args[6] - sets for the analysis (KEGG_metsig / KEGG_disease / BioCarta / GO_BP / GO_CC / GO_MF ; currently - human only)
#args[7] - change direction for the analysis (Combined - both genes with the highest and the lowest measure are considered and grouped together;
#                                             Separate -  both genes with the highest and the lowest measure are considered but analysed separately;
#                                             Down - only genes with the lowest measure are considered (useful if the measure is p-value);
#                                             Up - only genes with the highest measure are considered)

#args[8] - qvalue threshold for plotting
setwd("J:/Bioinformatics-R-Projects/KEGG-Analysis/ConInf")

args <- commandArgs(trailingOnly = TRUE)

args[1] <- "ConInf.txt"
args[2] <- "YES"
args[3] <- 2 # number of column with gene IDs (Gene symbols or entrez IDs - autodetect)
args[4] <- 4# number of column with gene measure
args[5] <- "Non-Parametric" # Statistics type (Parametric / Non-parametric)

args[6] <-  "KEGG_metsig" #sets for the analysis (KEGG_metsig / KEGG_disease / BioCarta / GO_BP / GO_CC / GO_MF ; currently - human only)
args[7] <- "Separate" # change direction for the analysis (Combined - both genes with the highest and the lowest measure are considered and grouped together;
args[8] <- "0.5" #qvalue threshold for plotting


options(echo=FALSE)
options(bitmapType='cairo')


library("gage")
library("gageData")
library("pathview")

options(stringsAsFactors = FALSE)

#parse arguments

infile = args[1]

if (args[2] == "YES") {
    header_present = TRUE
} else {
    header_present = FALSE
}

ids_col = as.numeric (args[3])
meas_col = as.numeric (args[4])

stat_t = args[5]

gset = args [6]

direction = args[7]

qval_thresh = as.numeric (args[8])

#prepare input data

tdata = read.table (infile, header = header_present)

measures = tdata [, meas_col]

if (is.numeric (tdata [1, ids_col]) == TRUE) {
   # identifiers are ENTREZ IDs
   names (measures) = tdata [,ids_col]
} else {
   # convert identifiers from gene symbols to ENTREZ ids
   data(egSymb)
   names (measures) = sym2eg (tdata [, ids_col])
}

if (gset == "KEGG_metsig") {
    kegg = TRUE
    data(kegg.gs)
    gene_sets = kegg.gs
} else if (gset == "KEGG_disease") {
    kegg = TRUE
    data(kegg.gs.dise)
    gene_sets = kegg.gs.dise
} else if (gset == "BioCarta") {
    kegg = FALSE
    data(carta.hs)
    gene_sets = carta.hs
} else if (gset == "GO_BP") {
    kegg = FALSE
    data(go.sets.hs)
    data(go.subs.hs)
    gene_sets = go.sets.hs[go.subs.hs$BP]
} else if (gset == "GO_CC") {
    kegg = FALSE
    data(go.sets.hs)
    data(go.subs.hs)
    gene_sets = go.sets.hs[go.subs.hs$CC]
} else if (gset == "GO_MF") {
    kegg = FALSE
    data(go.sets.hs)
    data(go.subs.hs)
    gene_sets = go.sets.hs[go.subs.hs$MF]
} else {  #actually we'll never be here for correct arguments
   cat ("Invalid id for sets for the analysis: ", gene_sets, "\n", sep="")
   cat ("Using all kegg pathways\n");
   kegg = TRUE
   data(kegg.sets.hs)
   gene_sets = kegg.sets.hs
}


# Apply gage
if (direction == "Combined") {
    samedir = FALSE
    writegreater = TRUE
    writeless = FALSE
} else if (direction == "Separate") {
    samedir = TRUE
    writegreater = TRUE
    writeless = TRUE
} else if (direction ==  "Down") {
    samedir = TRUE
    writegreater = FALSE
    writeless = TRUE
} else { (direction ==  "Down") - #if the argument is correct
    samedir = TRUE
    writegreater = TRUE
    writeless = FALSE

}

if (stat_t == "Parametric") {
    ranktest = FALSE
} else {
    ranktest = TRUE
}

gage_res = gage(measures, gsets = as.list(gene_sets), ref = NULL, samp = NULL, same.dir = samedir, rank.test = ranktest)

if (writegreater == TRUE) {
    cat ("\"Set\"\t", file="ResSets_statsGreater.txt", append=FALSE)
    write.table (gage_res$greater [, c ("p.val", "q.val", "set.size")], append=TRUE, file="ResSets_statsGreater.txt", quote=TRUE, sep="\t", row.names=TRUE, col.names=TRUE)
}
if (writeless == TRUE) {
    cat ("\"Set\"\t", file="ResSets_statsLess.txt", append=FALSE)
    write.table (gage_res$less [, c ("p.val", "q.val", "set.size")], append=TRUE, file="ResSets_statsLess.txt", quote=TRUE, sep="\t", row.names=TRUE, col.names=TRUE)
}


#Finally, for the KEGG case prepare png figures
if (kegg == TRUE) {
    path_g_cut = vector ()
    path_l_cut = vector ()

    if (writegreater == TRUE) {
        tab=gage_res$greater
        ind = (tab[, "q.val"] < qval_thresh & !is.na(tab[, "q.val"]))
        if (length (ind) > 0) {
            path_g <- rownames(tab)[ind]
            path_g_cut <- substr(c(path_g), 1, 8)
        }
    }

    if (writeless == TRUE) {
        tab=gage_res$less
        ind = (tab[, "q.val"] < qval_thresh & !is.na(tab[, "q.val"]))
        if (length (ind) > 0) {
            path_l <- rownames(tab)[ind]
            path_l_cut <- substr(c(path_l), 1, 8)
        }
    }

    
    path_cut = unique (c (path_g_cut, path_l_cut))
    if (length (path_cut) > 0) {
        out_list <- sapply(path_cut, function(pid) pathview(gene.data=measures, pathway.id = pid, species = "hsa", out.suffix="ResSets"))
    }
}

