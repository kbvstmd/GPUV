# Replot S5 from the archived source partition and allocation tables.
if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
file_arg <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE))
root <- dirname(dirname(dirname(dirname(normalizePath(file_arg, winslash='/')))))
data <- file.path(root, 'data', 'Fig_S5')
out <- file.path(root, 'pdf', 'redrawn', 'Fig_S5')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
colors <- c('#D4A24C', '#3B7EA8', '#7466A3')
labels <- c('Current study and INPHARED', 'Current study only', 'INPHARED only')

plot_pie <- function() {
  values <- read.delim(file.path(data, 'Panel_A_pie.tsv'))
  par(mar=c(2,0,4,0), family='Arial')
  pie(values$branch_length, labels=NA, col=colors, border='white', clockwise=TRUE,
      init.angle=90, radius=.79)
  angles <- pi/2 - 2*pi*(cumsum(values$fraction_total) - values$fraction_total/2)
  text(.52*cos(angles), .52*sin(angles), sprintf('%.1f%%', values$percent_total),
       col=c('#24323B', 'white', 'white'))
  mtext('Percentage of TerL phylogenetic diversity (PD)', side=3, font=2)
  legend('bottom', legend=labels, fill=colors, bty='n', cex=.73, xpd=NA)
}

plot_allocation <- function(filename) {
  d <- read.delim(file.path(data, filename))
  par(mar=c(.5, .5, 3, .5), family='Arial', xpd=NA)
  plot.new()
  n <- nrow(d)
  plot.window(xlim=c(-1.02, 1.98), ylim=c(.4, n+.6), xaxs='i', yaxs='i')
  for (i in seq_len(n)) {
    y <- n-i+1
    fractions <- unlist(d[i, c('shared_fraction','urban_only_fraction','reference_only_fraction')])
    edges <- c(0, cumsum(fractions))
    for (j in which(fractions > 0)) {
      rect(edges[j], y-.43, edges[j+1L], y+.43, col=colors[j], border='#303030')
    }
    rect(1.07, y-.43, 1.07+d$total_fraction[i], y+.43, col='#CDD5DA', border='#303030')
    text(-.035, y, d$category[i], adj=1, cex=.86)
    text(1.095+d$total_fraction[i], y, sprintf('%.2f%%', 100*d$total_fraction[i]), adj=0, cex=.82)
  }
  mtext('Percentage of PD*', side=3, at=.5, font=2)
  mtext('Total PD*', side=3, at=1.55, font=2)
}

for (panel in c('a','b','c')) {
  cairo_pdf(file.path(out, paste0('Fig_S5', panel, '.pdf')), width=if(panel=='a') 5.3 else 6,
            height=4.75, family='Arial')
  switch(panel, a=plot_pie(), b=plot_allocation('Panel_B_family_wide.tsv'),
         c=plot_allocation('Panel_C_host_phylum_wide.tsv'))
  dev.off()
}
