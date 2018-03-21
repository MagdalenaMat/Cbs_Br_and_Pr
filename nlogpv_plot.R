m <- to_GSEA$DCIS_IDC[which(to_GSEA$DCIS_IDC[,4] > 2), ]
a <- list(
  x = m[,5],
  y = m[,4],
  text = m[,1],
  xref = "x",
  yref = "y",
  showarrow = TRUE,
  arrowhead = 7,
  ax = 20,
  ay = -40
)
print(plot_ly(data = to_GSEA$DCIS_IDC, x= ~to_GSEA$DCIS_IDC[,5], y = ~to_GSEA$DCIS_IDC[,4])) %>% add_markers() %>%
  layout(annotations = a)