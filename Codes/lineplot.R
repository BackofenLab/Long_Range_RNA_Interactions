
library(tidyverse)
library(ggnewscale) # for separated color scales

# load rri data
rri <- 
  read_delim(
    "Results/IntaRNA_raw_output.txt",
    delim=";",
    comment = "#",
    col_names = T,
    skip = 2
  ) |> 
  filter(id1 != "id1", !is.na(start1)) |> 
  mutate(id = str_remove(id1,".5UTR")) |> 
  select(-c(id1,id2,subseqDP, hybridDP)) |> 
  mutate(across(!id, as.numeric)) |> 
  distinct() |> 
  pivot_longer(
    cols=matches("[12]$"),
    names_pattern=("(.*)(1|2)"),
    names_to=c("pos","utr")
  ) |> 
  mutate( utr = ifelse(utr==1,"5","3")) |> 
  pivot_wider(
    values_from = value,
    names_from = pos
  )

# length of CDS extension of each UTR
lenCDSext <- 100

# load genome data
genome <-
  read_delim(
    "Data/parameter_table.csv"
  ) |> 
  select(-c(seq5,seq3)) |> 
  pivot_longer(starts_with("UTR"), names_to="utr",values_to = "len") |> 
  mutate(utr = str_extract(utr,"\\d")) |> 
  mutate(id = fct_reorder2(id,id,class)) |> 
  mutate(lenCDS = ifelse(utr=="5",lenCDSext,-lenCDSext))



genome |> 
  filter(utr==5) |> 
  ggplot() +
  xlim(-350,100)+
  geom_segment( aes(x = -len,y=id,yend=id), xend=0, size=2, col="gray" ) + # UTR
  geom_segment( aes(y=id,yend=id, col=class), x=0, xend=lenCDS, size=2 ) + # CDS
  new_scale_color() +
  geom_segment(
    data= rri |> filter(utr=="5"),
    aes(x=start,xend=end,y=id,yend=id,alpha=-E),size=2, col="darkblue"
  ) + 
  theme_light()
  
