
library(tidyverse)
library(ggnewscale) # for separated color scales
library(ggpubr) # multiplot grid

# load rri data
rri <- 
  read_delim(
    "Results/IntaRNA_raw_output.txt",
    delim=";",
    comment = "#",
    col_names = T,
    skip = 2
  ) |> 
  filter(id1 != "id1", # remove additional header
         !is.na(start1)) |> # remove genome name separators
  mutate(id = str_remove(id1,".5UTR")) |> # get genome id
  select(-c(id1,id2,subseqDP, hybridDP)) |> # drop obsolete stuff
  mutate(across(!id, as.numeric)) |> # get number data type
  distinct() |> # drop duplicates
# decompose start/end information for 5' and 3' UTR in 3 steps
  pivot_longer(
    cols=matches("[12]$"),
    names_pattern=("(.*)(1|2)"),
    names_to=c("pos","utr")
  ) |> # (1) merge all 4 cols into long format
  mutate( utr = ifelse(utr==1,"5","3")) |> # get correct UTR labels
  pivot_wider(
    values_from = value,
    names_from = pos
  )  # reform start/end columns for each utr

# length of CDS extension of each UTR
lenCDSext <- 100

# load genome data
genome <-
  read_delim(
    "Data/parameter_table.csv"
  ) |> 
  select(-c(seq5,seq3)) |> # drop obsolete data
  pivot_longer(starts_with("UTR"), names_to="utr",values_to = "len") |> # separate UTR lengths into rows
  mutate(utr = str_extract(utr,"\\d")) |> # reduce UTR annotation (as in rri)
  mutate(id = fct_reorder2(id,id,class)) |> # order/group genomes by class
  mutate(lenCDS = ifelse(utr=="5",lenCDSext,-lenCDSext), # get CDS ends w.r.t. UTR type
         len = ifelse(utr=="5",-len,len)) # UTR indexing/length w.r.t. UTR type



pu5 <- 
genome |> 
  filter(utr==5) |> 
  ggplot() +
  xlim(NA,100)+ # set max x since not correct via next segment call
  geom_segment( aes(x = len,y=id,yend=id), xend=0, size=2, col="gray" ) + # UTR
  geom_segment( aes(y=id,yend=id, col=class, xend=lenCDS), x=0, size=2 ) + # CDS
  new_scale_color() + # new color scale for E-based RRI coloring
  geom_segment(
    data= rri |> filter(utr=="5"),
    aes(x=start,xend=end,y=id,yend=id,alpha=-E),size=2, col="darkblue"
  ) + # RRIs
  theme_light() +
  labs(
    x="position w.r.t. 5'UTR-CDS transition",
  ) +
  theme(
    axis.title.y = element_blank()
  )

# genome |> 
#   filter(utr==3) |> 
#   ggplot() +
#   xlim(-100,1250)+
#   geom_segment( aes(x = len,y=id,yend=id), xend=0, size=2, col="gray" ) + # UTR
#   geom_segment( aes(y=id,yend=id, col=class, xend=lenCDS), x=0, size=2 ) + # CDS
#   new_scale_color() +
#   geom_segment(
#     data= rri |> filter(utr=="3"),
#     aes(x=start,xend=end,y=id,yend=id,alpha=-E),size=2, col="darkblue"
#   ) + 
#   theme_light()
  
pu3 <- 
genome |> 
  filter(utr==3) |> 
  ggplot() +
  xlim(NA,0)+ # set max x since not correct via next segment call
  geom_segment( aes(y=id,yend=id, col=class, xend=lenCDS-len, x=-len), size=2 ) + # CDS
  geom_segment( aes(x = -len,y=id,yend=id), xend=0, size=2, col="gray" ) + # UTR
  new_scale_color() + # new col scale for E-based RRI coloring
  geom_segment(
    data= rri |> filter(utr=="3") |> left_join(genome),
    aes(x=start-len,xend=end-len,y=id,yend=id,alpha=-E),size=2, col="darkblue"
  ) + # RRIs
  theme_light() +
  labs(
    x="position w.r.t. genome end",
  ) +
  theme(
    axis.title.y = element_blank()
  )


ggarrange(
  pu5,
  pu3,
  nrow = 2, ncol=1
) 

ggsave("Results/rri-position.pdf", width=21, height=30, units="cm")  
ggsave("Results/rri-position.png", width=21, height=30, units="cm")  
