library(dplyr)

DBinchikey <- read.csv("DBinchikey.csv")
PathwayID_name <- read.csv("PathwayID_name.csv")
KEGG_pathwayID <- read.csv("KEGG_pathwayID.csv")
inchkey.test <- read.csv("inchikey.test.csv")

inchilist <-as.character(inchkey.test$Inchikey)

DBinchikey %>% 
  filter(Inchikey %in% inchilist) %>% 
  as.data.frame() -> match.result  # match with local database
write.csv(match.result,"matched.result.csv")
cat("there are",(length(inchilist)-nrow(match.result)),"Inchikeys NOT matched in KEGG database")


# match.result %>%
#   select(KEGGID) %>%
#   write.csv(.,"matched KEGGID.csv")

match.result %>% 
  select(pathway_1:pathway_6) -> pathway.number # pathways related

pathway.number %>% 
  unlist() %>% 
  as.vector() %>% 
  as.factor() %>% 
  levels() %>%
  .[2:length(.)] -> 
  levesl.number   # all realted pathways


Moval <- function(vector,factor){
  vector <- vector[-which(vector == factor)]
  return(vector)
}

Moval.pathway <-  c("map01100","map01060","map01110","map01120",
                    "map02010","map01130")

for (i in 1:length(Moval.pathway)) {
  levesl.number <- 
    Moval(levesl.number, Moval.pathway[i])
}


# counts of every pathways
Counts<-NULL
Pathwayname <- NULL
Total <- NULL
for (i in 1:length(levesl.number)) {
  Counts[i] <- nrow(pathway.number %>% 
                      filter_all(any_vars(. == levesl.number[i])))
  
  Pathwayname[i] <- PathwayID_name %>%
    filter(ID == levesl.number[i]) %>%
    select(pathway) %>%
    pull %>%
    as.character() 
  
  Total[i] <-  KEGG_pathwayID %>%
    filter_all(any_vars(. == levesl.number[i])) %>% 
    select(KEGG_ID) %>%
    nrow()
}

data.frame(pathway = levesl.number,
           enrichment = Counts,
           Pathwayname = Pathwayname,
           Total = Total) %>% 
  arrange(desc(enrichment)) -> pathway.order # pathway Impact order

pathway.order %>% slice(1:20) %>% tbl_df()

probabilities <- NULL
for (i in 1:nrow(pathway.order)) {
  m <- pathway.order %>% slice(i) %>% select(Total) %>% pull
  n <- KEGG_pathwayID[KEGG_pathwayID$pathway_1 != "",] %>% nrow()
  k <- length(inchilist)
  x <- pathway.order %>% slice(i) %>% select(enrichment) %>% pull
  probabilities[i] <- dhyper(x, m, n, k, log = FALSE)
}

pathway.order$pvalue <- probabilities
pathway.order$ajust.pvalue <- round(p.adjust(probabilities), 5)
pathway.order <- pathway.order %>% arrange(ajust.pvalue)
pathway.order <- pathway.order %>% 
  mutate(Overlap = round(enrichment/Total*100,2),
        logPjust = -log10(ajust.pvalue),
        logP = -log10(pvalue)) %>%
  arrange(desc(Overlap))
pathway.order <- pathway.order %>% mutate(Co = ifelse(ajust.pvalue < 0.05, "red", "blue"))
pathway.order <- pathway.order %>% mutate(Co2 = ifelse(pvalue < 0.05, "red", "blue"))




pathway.order %>% select(pathway) %>%
  pull(1) %>%
  as.character() %>%
  .[1]

match.result %>%
filter(pathway_1 == "map02060")















library(ggplot2,ggpubr)
ggpubr::ggscatter(pathway.order,
                  x = "Overlap",
                  y = "logPjust",
                  xlab = "Overlap (%)",
                  ylab = "-log10(ajust.pvalue)",
                  size = "enrichment",
                  fill = "Co",
                  shape = 21,
                  color = "black",
                  legend = ""
                  ) +
  scale_fill_manual(values = c("grey","red")) +
  scale_size(range = c(2,10))+
  ggrepel::geom_text_repel(
    data = pathway.order %>% slice(1:20) %>% filter(ajust.pvalue  < 0.05),
    aes(x = Overlap,
        y = logPjust,
        label = Pathwayname),
    size = 3,
    segment.color = "black", show.legend = FALSE )+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")
ggsave("baboplot.ajustP.png")


ggpubr::ggscatter(pathway.order,
                  x = "Overlap",
                  y = "logP",
                  xlab = "Overlap (%)",
                  ylab = "-log10(pvalue)",
                  size = "enrichment",
                  fill = "Co2",
                  shape = 21,
                  color = "black",
                  legend = ""
) +
  scale_fill_manual(values = c("grey","red")) +
  scale_size(range = c(2,10))+
  ggrepel::geom_text_repel(
    data = pathway.order %>% slice(1:10) %>% filter(pvalue < 0.01),
    aes(x = Overlap,
        y = logP,
        label = Pathwayname),
    size = 3,
    segment.color = "black", show.legend = FALSE )+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")
ggsave("baboplot.pvalue.png")



pathway.order2 <- pathway.order %>% arrange(desc(logPjust)) %>%
  filter(pvalue < 0.05) %>% slice(1:20)
ggpubr::ggbarplot(pathway.order2, 
                  x = "Pathwayname", y = "logPjust",
          xlab = "",
          ylab = "-log10(ajust.p)",
          fill = "Co",legend = "",
          sort.val = "asc",
          font.main = c(12, "plain", "black"),            
          font.x = c(14, "plain", "black"),                  
          font.y = c(14, "plain", "black"),  
          font.legend = c(14, "plain", "black")) +
  coord_flip() +
  scale_fill_manual(values = c("grey","red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")
  
ggsave("barplot.ajustp.png")



ggpubr::ggbarplot(pathway.order2, 
                  x = "Pathwayname", y = "logP",
                  xlab = "",
                  ylab = "-log10(p)",
                  fill = "Co2",legend = "",
                  sort.val = "asc",
                  font.main = c(12, "plain", "black"),            
                  font.x = c(14, "plain", "black"),                  
                  font.y = c(14, "plain", "black"),  
                  font.legend = c(14, "plain", "black")) +
  coord_flip() +
  scale_fill_manual(values = c("red","red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

ggsave("barplot.pvalue.png")
