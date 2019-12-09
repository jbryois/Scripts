library(stringr)

d <- c("Darina Czamara, Gökçen Eraslan, Christian M Page, Jari Lahti, Marius Lahti-Pulkkinen, Esa Hämäläinen, Eero Kajantie, Hannele Laivuori, Pia M Villa, Rebecca M Reynolds, Wenche Nystad, Siri E Håberg, Stephanie J London, Kieran J O’Donnell, Elika Garg, Michael J Meaney, Sonja Entringer, Pathik D Wadhwa, Claudia Buss, Meaghan J Jones, David TS Lin, Julie L MacIsaac, Michael S Kobor, Nastassja Koen, Heather J Zar, Karestan C Koenen, Shareefa Dalvie, Dan J Stein, Ivan Kondofersky, Nikola S Müller, Fabian J Theis, Katri Räikkönen, Elisabeth B Binder")
format_names <- function(d) {
d_split <- str_split(d,",") %>% unlist() %>% str_trim()
first_letter <- substr(d_split,1,1)
second_letter <- str_extract(d_split," [A-Z]{1,2} ") %>% str_trim()
second_letter[is.na(second_letter)] <- ""
initials <- paste0(first_letter,second_letter)
last_name <- str_extract(d_split," [A-Z].+$") %>% str_trim() %>% str_replace("[A-Z]{1,2} ", "")
formated_names <- paste(last_name,initials,collapse =", ")
}

formated <- format_names(d)
write(formated,"~/Desktop/formated_names.txt")