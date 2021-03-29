suppressMessages(library(tidyverse))
suppressMessages(library(pacman))
suppressMessages(library(data.table))

wkPath <- c('./1.plot',  './2.results','./3.input', './4_download')
for(i in wkPath){
  wkPathi = i
  # wkPathi = paste0(sectionName, '/', i)
  #每一个子项目都含plot、result、input
  if (!dir.exists(wkPathi)) dir.create(wkPathi)
}
rm(list=c('i', 'wkPathi', 'wkPath'))
