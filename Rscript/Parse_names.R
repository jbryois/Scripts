library(stringr)

d <- c("Julien Bryois, Nathan G Skene, Thomas Folkmann Hansen, Lisette Kogelman, Hunna J Watson, Anorexia Nervosa Working Group of the Psychiatric Genomics Consortium, International Headache Genetics Consortium, 23andMe, Gerome Breen, Cynthia Bulik, Ernest Arenas, Jens Hjerling-Leffler, Patrick F Sullivan")
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