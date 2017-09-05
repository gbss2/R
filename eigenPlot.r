#==========================================
#
# (/u0254) Copyleft 2016, by GBSS and Contributors.
#
# 
# -----------------
# eigenPlot.r
# -----------------
# GNU GPL 2016, by GBSS and Contributors.
#
# Original Author: Giordano Bruno Soares-Souza
# Contributor(s): John Baums (developer of iwanthue)
# Updated by: Giordano Bruno Soares-Souza (04/09/2017)
#
# Command line: Rscript eigenPlot(plotMode("F"or"P"),<filename.evec>,popFile)
# Dependencies: R, library gtools
# Description: Plot Eigenstrat PCA plots
# "F" = family mode
# "P" = population mode
# Requirements: For pedigrees, family ids MUST BE different from Individuals IDs
# 				For reference populations, fids MUST BE the same as individuals ids 
# POP FILE Format:
#	IID1	SUBPOP1	POP1 (last column is optional)
#	IID2	SUBPOP2	POP2 
# Future Developments:
# 	1) Add option to use external file to assign individuals to pedigrees
#	2) Fix the issue related to more than 24 ref pops in family mode
#	3) Allow users to choose the pcs to be ploted
#	4) Improve plot aestethics
#	5) Group subpopulations by 3rd column 
################################################################################

if(!require(colorspace)) { install.packages("colorspace"); require(colorspace)}

# Above is the function Iwanthue written John Baums under generic use/modify license (https://gist.github.com/johnbaums/45b49da5e260a9fc1cd7 ; 	https://github.com/johnbaums/hues/blob/master/R/iwanthue.R)

#' Generate a colour palette by k-means clustering of LAB colour space.
#' 
#' Generate a palette of distinct colours through k-means clustering of LAB 
#' colour space.
#' 
#' @param n Numeric. The number of colours to generate.
#' @param hmin Numeric, in the range [0, 360]. The lower limit of the hue range 
#'   to be clustered.
#' @param hmax Numeric, in the range [0, 360]. The upper limit of the hue range 
#'   to be clustered.
#' @param cmin Numeric, in the range [0, 180]. The lower limit of the chroma 
#'   range to be clustered.
#' @param cmax Numeric, in the range [0, 180]. The upper limit of the chroma 
#'   range to be clustered.
#' @param lmin Numeric, in the range [0, 100]. The lower limit of the luminance 
#'   range to be clustered.
#' @param lmax Numeric, in the range [0, 100]. The upper limit of the luminance 
#'   range to be clustered.
#' @param plot Logical. Should the colour swatches be plotted (using 
#'   \code{\link{swatch}})?
#' @param random Logical. If \code{TRUE}, clustering will be determined by the 
#'   existing RNG state. If \code{FALSE}, the seed will be set to \code{1} for 
#'   clustering, and on exit, the function will restore the pre-existing RNG 
#'   state.
#' @return A vector of \code{n} colours (as hexadecimal strings), representing 
#'   centers of clusters determined through k-means clustering of the LAB colour
#'   space delimited by \code{hmin}, \code{hmax}, \code{cmin}, \code{cmax}, 
#'   \code{lmin} and \code{lmax}.
#' @details Note that \code{iwanthue} currently doesn't support \code{hmin} 
#'   greater than \code{hmax} (which should be allowed, since hue is circular).
#' @references 
#' \itemize{
#'   \item \href{http://tools.medialab.sciences-po.fr/iwanthue/}{iwanthue - colors for data scientists}
#'   \item \href{https://github.com/medialab/iwanthue}{iwanthue on
#'   GitHub}   
#' }
#' @seealso \code{\link{swatch}}
#' @export
#' @importFrom colorspace LAB hex coords
#' @examples 
#' iwanthue(5)
#' iwanthue(5, plot=TRUE)
#' iwanthue(5, 0, 240, 0, 24, 0, 100, plot=TRUE) # shades
#' iwanthue(5, 0, 360, 0, 54, 67, 100, plot=TRUE) # pastel
#' iwanthue(5, 0, 360, 54, 180, 27, 67, plot=TRUE) # pimp
#' iwanthue(5, 0, 360, 36, 180, 13, 73, plot=TRUE) #intense
#' iwanthue(3, 0, 300, 60, 180, 73, 100, plot=TRUE) # fluoro
#' iwanthue(3, 220, 260, 12, 150, 0, 53, plot=TRUE) # blue ocean
iwanthue <- function(n, hmin=0, hmax=360, cmin=0, cmax=180, lmin=0, lmax=100, 
                     plot=FALSE, random=FALSE) {
  stopifnot(hmin >= 0, cmin >= 0, lmin >= 0, 
            hmax <= 360, cmax <= 180, lmax <= 100, 
            hmin <= hmax, cmin <= cmax, lmin <= lmax,
            n > 0)
  if(!random) {
    if (exists(".Random.seed", .GlobalEnv)) {
      old_seed <- .GlobalEnv$.Random.seed
      on.exit(.GlobalEnv$.Random.seed <- old_seed)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    set.seed(1)
  }
  lab <- LAB(as.matrix(expand.grid(seq(0, 100, 1), 
                                   seq(-100, 100, 5), 
                                   seq(-110, 100, 5))))
  if (any((hmin != 0 || cmin != 0 || lmin != 0 ||
           hmax != 360 || cmax != 180 || lmax != 100))) {
    hcl <- as(lab, 'polarLUV')
    hcl_coords <- coords(hcl)
    hcl <- hcl[which(hcl_coords[, 'H'] <= hmax & hcl_coords[, 'H'] >= hmin &
                       hcl_coords[, 'C'] <= cmax & hcl_coords[, 'C'] >= cmin & 
                       hcl_coords[, 'L'] <= lmax & hcl_coords[, 'L'] >= lmin), ]
    lab <- as(hcl, 'LAB')    
  }
  lab <- lab[which(!is.na(hex(lab))), ]
  clus <- kmeans(coords(lab), n, iter.max=50)
  if (isTRUE(plot)) {
    swatch(hex(LAB(clus$centers)))
  }
  hex(LAB(clus$centers))
}

# Get arguments
# plotMode = F for pedigrees and P for populations 
# filename = evec file from smartpca
# popfile (optional)
args = commandArgs(trailingOnly=TRUE)
plotMode = args[1]
filename = args[2]
popfile = args[3]

# Retrive Eigenvalues and Eigenvectors from evec file (Unix based)
system(paste0("head -n 1 ",filename," | sed -r 's/^\\s+//g' | sed -r 's/\\s+/\t/g' | sed -r 's/\\#//g' > ",filename,".head.tmp"), intern=TRUE)
system(paste0("tail -n +2 ",filename," | sed -r 's/^\\s+//g' | sed -r 's/\\s+/\t/g' | sed -r 's/\\:/\\t/g' > ",filename,".body.tmp"), intern=TRUE)
writeLines("Reading EVEC file")
eigVal = read.table(paste0(filename,".head.tmp"),head=F,as.is=T)
evec = read.table(paste0(filename,".body.tmp"),head=F,as.is=T)

# IF popfile is present
writeLines("Reading POP file")

if(!is.null(popfile)){

	filePOP = read.table(popfile,head=F,as.is=T,sep="\t")
	evec[,ncol(evec)] = filePOP[match(evec[,1],filePOP[,1]),2]
	print(dim(evec))
	pops = as.vector(unique(evec[,ncol(evec)]))
	
	if(dim(filePOP)[2]==3){
		writeLines("Hierarchical Ploting Selected")
		cont = filePOP[match(evec[,1],filePOP[,1]),3]
		conts = as.vector(unique(cont))
	}
	
	if(anyNA(evec[,ncol(evec)])){ 
	
		writeLines("Missing population ID for some individuals! See individuals_without_pop_ID.txt") 
		write.table(evec,'evec_with_pop_IDS.txt',col.names=F,row.names=F,quote=F)
		sink("individuals_without_pop_ID.txt")
		print(subset(evec,is.na(evec[,ncol(evec)]),select=c(1,ncol(evec))))
		sink()
	}
	
}


# PLOTING PCA
writeLines("Ploting PCAs")
eigVal[1] = 0
pdf(gsub(".evec",".pdf",filename))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
#print(eigVal)
#print(pops)

# PEDIGREE MODE
if(plotMode=="F"){

refPops = subset(evec,V1==V2)
fam = subset(evec,V1!=V2)

pchP = match(as.vector(refPops[,ncol(evec)]),pops)

	for (i in 4: (dim(evec)[2]-1)){
		eigX = eigVal[2]/sum(eigVal)
		eigY = eigVal[(i-1)]/sum(eigVal)
		plot(evec[,3],evec[,i],type="n",xlab="",ylab="")
		points(fam[,3],fam[,i],col="black", bg=colors()[as.integer(evec[,1])], pch=21)
		points(refPops[,3],refPops[,i],col=24,pch=pchP)
		title(ylab=paste0("PC ",i-2," - Var explained = ",round(eigY[[1]],3)),xlab=paste0("PC 1 - Var explained = ",round(eigX[[1]],3)), line=2, cex.lab=1)
		legend("topright",inset=c(-0.36,0),legend=pops,pch=c(21,2:20),col=24,bty = "n")
	}

}

# POPULATION MODE
if(plotMode=="P"){

	inds = match(as.vector(evec[,ncol(evec)]),pops)
	colorPops = iwanthue(length(pops))
	colorPlot = colorPops[inds]
	
	# IF POPFILE HAS 3 columns PCHs = 3rd column and Colors = 2nd column
	if(dim(filePOP)[2]==3){
		indsCon = match(as.vector(cont),conts)
		pchs = c(15,16,17,18,25,8,3,4,20,11,9,13,7,12,14)
		pchPops = rep(pchs,len=length(conts))
		pchPlot = pchPops[indsCon]
		pops2legend = pops
		pch2legend = as.numeric(unique(cbind(evec[,ncol(evec)],cont,pchPlot))[,3])
		sink("pchs.txt")
		print(pchPlot)
		print(colorPlot)
		sink()
	
	} else {
		# IF POPFILE HAS MORE THAN 15 POPULATIONS, USE DISTINCT PCHS TO DIFERENTIATE
		if(length(pops)>15){
			pchs = c(15,16,17,18,25,8,3,4,20,11,9,13,7,12,14)
			pchPops = rep(pchs,len=length(pops))
			pchPlot = pchPops[inds]
			pops2legend = pops
			pch2legend = pchPops
			sink("pchs.txt")
			print(pchPlot)
			sink()
		} else { 
			# DEFAULT PLOT
			pchPlot = 20
			pch2legend = 20
		}
	}

	pchPlot = as.numeric(pchPlot)
	pchPops = as.numeric(pchPops)
	
	#print(colorPops)
	print(dim(evec))
	cexAdjust=1-(length(pops2legend)/100)

		for (i in 4: (dim(evec)[2]-1)){
			eigX = eigVal[2]/sum(eigVal)
			eigY = eigVal[(i-1)]/sum(eigVal)
	#		print(eig)
	#		print(eig[[1]])
			plot(evec[,3],evec[,i],xlab="",ylab="",col=colorPlot,pch=pchPlot )
			title(ylab=paste0("PC ",i-2," - Var explained = ",round(eigY[[1]],3)),xlab=paste0("PC 1 - Var explained = ",round(eigX[[1]],3)), line=2, cex.lab=1)
			legend("topright",inset=c(-0.35,0),legend=pops2legend,pch=pch2legend,col=colorPops,bty = "n",cex=cexAdjust,pt.cex=cexAdjust)
			print(i)
		}

}

dev.off()

