# Render measured timings and worker intervals, without rerunning fits.
# Rscript scripts/plot_parallel_benchmark.R [results_dir] [figures_dir]
args <- commandArgs(trailingOnly = TRUE)
results <- if (length(args)) args[1] else 'docs/results/parallel'
figures <- if (length(args) > 1) args[2] else 'docs/figures/parallel'
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
s <- read.csv(file.path(results, 'summary.csv'))
r <- read.csv(file.path(results, 'runs.csv'))
d <- read.csv(file.path(results, 'starts.csv'))
cols <- c(ACL = '#2666A5', ALSCL = '#CE682D')
meta <- readLines(file.path(results, 'environment.txt'))
years <- sub('.*; ([0-9]+) years.*', '\\1', grep('^Data:', meta, value = TRUE))
cpu <- sub('^CPU: ', '', grep('^CPU:', meta, value = TRUE))
if (!length(cpu)) cpu <- 'shared computer'
render <- function(name, draw, width = 10, height = 5) {
  svg(file.path(figures, paste0(name, '.svg')), width = width, height = height,
      family = 'sans', bg = 'white'); draw(); dev.off()
  png(file.path(figures, paste0(name, '.png')), width = width * 160, height = height * 160,
      res = 160, type = 'cairo', bg = 'white'); draw(); dev.off()
}
render('timing', function() {
  par(mfrow = c(1,2), mar = c(4.5,4.5,3.5,1), oma = c(2,0,2,0), las = 1)
  for (model in c('ACL','ALSCL')) {
    z <- s[s$model == model, ]; raw <- r[r$model == model, ]
    plot(z$ncores, z$median_seconds, type = 'n', xlim = c(.6,4.4),
         ylim = c(0,max(z$max_seconds)*1.20), xaxt = 'n',
         xlab = 'Socket workers (ncores)', ylab = 'Total fit time (seconds)', main = model)
    axis(1, at = c(1,2,4)); abline(h = axTicks(2), col = '#E8ECF0')
    points(raw$ncores, raw$elapsed_seconds, pch = 1, col = adjustcolor(cols[model],.55))
    arrows(z$ncores,z$min_seconds,z$ncores,z$max_seconds, angle = 90,code = 3,
           length = .04,col = cols[model])
    lines(z$ncores,z$median_seconds,col = cols[model],lwd = 2)
    points(z$ncores,z$median_seconds,pch = 19,col = cols[model],cex = 1.2)
    text(z$ncores,z$median_seconds,labels = sprintf('%.2f s\n%.2fx',z$median_seconds,z$speedup),
         pos = 3,cex = .85,offset = .8)
  }
  mtext('Same four starts; two optimizer passes per start', side = 3,outer = TRUE,font = 2)
  mtext(sprintf('Median and range of %s runs | Synthetic YTF, %s years | %s, shared computer',
                paste(unique(s$repetitions), collapse = '/'), years, cpu),
        side = 1,outer = TRUE,cex = .8)
})
render('worker-timeline', function() {
  par(mfrow = c(1,2), mar = c(4.5,4.6,3.5,1), oma = c(2,0,2,0), las = 1)
  for (model in c('ACL','ALSCL')) {
    run <- r[r$model == model & r$ncores == 4 & r$repetition == 1, ]
    z <- d[d$run_id == run$run_id, ]; offset <- min(z$started_at)
    lo <- z$started_at - offset; hi <- z$finished_at - offset
    plot(NA, xlim = c(0,max(hi)*1.1), ylim = c(.5,4.5), yaxt = 'n',
         xlab = 'Seconds since first worker started', ylab = '', main = model)
    axis(2, at = 1:4,labels = paste('Start',z$start_id))
    abline(v = axTicks(1),col = '#E8ECF0')
    segments(lo,1:4,hi,1:4,lwd = 15,col = adjustcolor(cols[model],.8),lend = 1)
    text((lo+hi)/2,1:4,labels = paste('PID',z$pid),col = 'white',cex = .8)
    mtext(sprintf('Worker CPU: %.2f s | interval: %.2f s',sum(z$cpu_seconds),max(hi)),
          side = 3,line = .4,cex = .78)
  }
  mtext('Four real fitting processes overlap in time',side = 3,outer = TRUE,font = 2)
  mtext('First 4-worker trial; intervals cover TMB objective construction and optimization only',
        side = 1,outer = TRUE,cex = .8)
})
