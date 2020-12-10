library(minisvg)
library(purrr)
library(tibble)
# library(datapasta) # format addin
vars <- tibble::tribble(
  ~id,    ~x,  ~y, ~width, ~height,
  "titlearea",     0,   0,    297,     110,
  "page",     0,   0,    297,     420,
  "body",  33.5, 240,    175,      85,
  "info",   210, 240,     70,      85,
  "footerarea",     0, 355,    297,      65,
  "picture",  33.5,  90,    230,     115,
  "caption",  33.5, 205,    230,      10,
  "qrmore", 243.5, 362,     20,      10,
  "qrcode", 243.5, 372,     20,      20,
  "address",   105, 372,    136,      20,
  "vuwlogo",  33.5, 372,     70,      20,
  "timeline",  33.5, 345,    230,      10
)

# 594 x 841mm

doc <- SVGDocument$new(width = '841mm', 
                       height = '594mm', 
                       viewBox = '0 0 841 594')
doc$add_css("
rect {
  stroke: black;
  fill: #08306B33;
  stroke-width: 0.5;
}
")

invoke(doc$append, pmap(vars, stag$rect))

doc$save("test_auto.svg")