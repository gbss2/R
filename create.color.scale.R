## Function to get gradient colors based on a vector of numeric values 
# x = vector of numeric values
# intervals = increments for color space generation 
# round = number of digits to round numbers
# ncolors = number of colors to interpolate from a given pallete
# palette = palette used to create the gradient of colors (eg. "YlOrRd")
# ATTENTION: intervals and round options should be adjusted according to the distribution of x values, 
# e.g. for integers, intervals should be set to 1 and round to 0, instead,
# to discriminate tenths, the best values are 0.1 for intervals and 1 to round, 
# for hundredths, 0.01 and 2, respectively, and so on.


create.color.scale <- function(x, intervals, round, ncolors, palette){
  x <- as.numeric(x)
  col.seq <- round(seq(min(x), max(x), intervals), round)
  col.seq <- append(col.seq,round(max(x),round))
  colors <- colorRampPalette(RColorBrewer::brewer.pal(ncolors, palette))(length(col.seq))
  x.map <- round(x,round)
  x.index <- match(x.map, col.seq)
  return(colors[x.index])
  
}
