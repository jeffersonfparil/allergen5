args = commandArgs(trailingOnly=TRUE)
# args = c("LOC124688403-LOC124686040.pw.aln.split.kaks", "0.001", "allergen-5a")
fname_input = args[1]
alpha = as.numeric(args[2])
gene = args[3]
cex=1
fname_output_svg = paste0(fname_input, ".svg")
fname_output_sig = paste0(fname_input, "-SIGNIFICANT_PEAKS.csv")

dat = read.delim(fname_input, header=TRUE)
dat$Sequence = gsub(">", "", dat$Sequence)
dat$Sequence = gsub(")", "", dat$Sequence)
dat$Sequence = gsub("\\(", "-", dat$Sequence)
dat$Sequence = gsub("---", "-", dat$Sequence)
X = matrix(unlist(strsplit(as.character(dat$Sequence), "-")), ncol=6, byrow=TRUE)
Species1 = X[1,1]
Gene1    = X[1,2]
Species2 = X[1,3]
Gene2    = X[1,4]
df = data.frame(
    Position_ini=as.numeric(X[,5]),
    Position_fin=as.numeric(X[,6]),
    KaKs=dat$Ka.Ks,
    p=dat$P.Value.Fisher.
)
df$KaKs[is.na(df$KaKs)] = 1.0

svg(fname_output_svg, height=5, width=8)

plot(df$Position_ini, df$KaKs, type="l",
    main=paste0(Gene1, " (", Species1, ") vs\n", Gene2, " (", Species2, ")"),
    xlab="Position (bp)", ylab="Ka/Ks")
grid()
idx = c(1:nrow(df))[!is.na(df$KaKs) & (df$p <= alpha) & (df$KaKs > 1)]
for (j in idx){
    text(x=df$Position_ini[j], y=df$KaKs[j], lab="*", col="red", cex=2)
}
abline(h=1, lty=2, col="red")
if (!exists("outdf")){
    outdf = df[idx, ]
} else {
    outdf = rbind(outdf, df[idx, ])
}
dev.off()
if (nrow(outdf) > 0) {
    write.table(outdf, file=fname_output_sig, sep=",", quote=FALSE, row.names=FALSE)
}
print(paste0("OUTPUTS: ", fname_input, ".{svg, sigout}"))

