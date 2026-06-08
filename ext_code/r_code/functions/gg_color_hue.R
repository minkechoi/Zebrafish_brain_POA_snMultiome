# ggplot2 colors # https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette/8197703#8197703
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}