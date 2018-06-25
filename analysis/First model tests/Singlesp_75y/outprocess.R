library(tidyverse);
library(stringr);

indir <- file.path("/home/luf74xx/Dokumente/model/inputs"); #CCTBserver
inid <- file.path("singlespinputj1.csv") # check in simID file in output
input <- as_tibble(read.table(file.path(indir, inid), header = TRUE, sep = ","))
outdir <- file.path("/home/luf74xx/Dokumente/model/EDoutputs");
outid <- file.path("singlesp1-75y"); #simID
outraw <- as_tibble(read.table(file.path(outdir, outid, "orgsweekly.csv"), header = TRUE, sep = "\t"));
loc <- gsub("[\\(|\\)]", "", outraw$location)
loc <- as.data.frame(matrix(unlist(str_split(loc, ",")),ncol =3,byrow = T))
names(loc) = c("xloc", "yloc", "floc")
masses = c()
for(plant in outraw$biomass){
  m <- unlist(str_extract_all(plant,"([0-9]{1,4}.[0-9]{1,7})")) # extracts both veg and repr
  if(length(m) < 2){
    m <- c(0,m)
  }
  masses <- rbind(masses,m)
}
masses <- as.data.frame(masses, stringsAsFactors = FALSE); #there is probably a better way to convert into numeric
masses <- as.data.frame(sapply(masses, as.numeric));
names(masses) <-c("repr", "veg");
outdata <- as_tibble(cbind(select(outraw,-one_of("location","biomass")),
                           loc,
                           masses))
rm(masses, loc, plant, m)
saveRDS(outdata,"/home/luf74xx/Dokumente/model/analysis/Singlesp_75y/outdata.rds")