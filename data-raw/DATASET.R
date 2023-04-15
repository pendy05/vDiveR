## code to prepare `DATASET` dataset goes here
proteins_1host<-read.csv('data-raw/proteins_1host.csv')
protein_2hosts<-read.csv('data-raw/protein_2hosts.csv')

JSON_sample <- jsonlite::fromJSON(txt = readLines(con = "data-raw/proteinA.json", encoding = "UTF-8"))
metadata<-read.csv('data-raw/metadata.csv')

usethis::use_data(proteins_1host, protein_2hosts, JSON_sample, metadata, overwrite = T)
