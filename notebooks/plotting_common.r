library(ggplot2)
library(scales)
library(tikzDevice)
library(ggrepel)
theme_set(theme_bw() + theme(text = element_text(size=10)))

# from https://stackoverflow.com/questions/30179442/plotting-minor-breaks-on-a-log-scale-with-ggplot
log_breaks = function(maj, radix=10) {
  function(x) {
    minx         = floor(min(logb(x,radix), na.rm=T)) - 1
    maxx         = ceiling(max(logb(x,radix), na.rm=T)) + 1
    n_major      = maxx - minx + 1
    major_breaks = seq(minx, maxx, by=1)
    if (maj) {
      breaks = major_breaks
    } else {
      steps = logb(1:(radix-1),radix)
      breaks = rep(steps, times=n_major) +
               rep(major_breaks, each=radix-1)
    }
    radix^breaks
  }
}
scale_x_log_eng = function(..., radix=10) {
  scale_x_continuous(...,
                     trans=log_trans(radix),
                     breaks=log_breaks(TRUE, radix),
                     minor_breaks=log_breaks(FALSE, radix))
}
scale_y_log_eng = function(..., radix=10) {
  scale_y_continuous(...,
                     trans=log_trans(radix),
                     breaks=log_breaks(TRUE, radix),
                     minor_breaks=log_breaks(FALSE, radix))
}

tikz_file <- function(filename, plot=last_plot(), width=5, height=2.5) {
        p = plot
        # p = plot + theme(text = element_text(size=16, family="CM Roman"),
        #                  axis.title.x = element_text(face="italic"),
        #                  axis.title.y = element_text(face="bold"))
        # f = paste(tools::file_path_sans_ext(filename), ".eps", sep="")
        # postscript(f, family="CM Roman", width = width, height = height, useKerning=TRUE)
        # file = set.file.extension(filename, extension=".png")
        file = filename
        print(sprintf("Writing to %s", file))
        # png(file=file)
        tikz(file = filename, width = width, height = height)
        print(plot)
        dev.off()
}
