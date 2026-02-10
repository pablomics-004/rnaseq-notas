library("sessioninfo")
library("here")
library("ggplot2")

## Hello world
print("Soy Pablo")

## Directorios
dir_plots <- here::here("figuras")
dir_rdata <- here::here("processed-data")

## Crear directorio para las figuras y archivos
dir.create(dir_plots, showWarnings = FALSE)
dir.create(dir_rdata, showWarnings = FALSE)

## Hacer una imagen de ejemplo
pdf(file.path(dir_plots, "mtcars_gear_vs_mpg.pdf"), useDingbats = FALSE)
ggplot(mtcars, aes(group = gear, y = mpg)) +
    geom_boxplot()
dev.off()

## Para reproducir mi código
options(width = 120)
sessioninfo::session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.5.2 (2025-10-31)
#  os       Ubuntu 24.04.3 LTS
#  system   x86_64, linux-gnu
#  ui       Positron
#  language (EN)
#  collate  C.UTF-8
#  ctype    C.UTF-8
#  tz       America/Mexico_City
#  date     2026-02-10
#  pandoc   NA
#  quarto   1.8.27 @ /usr/share/positron/resources/app/quarto/bin/quarto

# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package      * version date (UTC) lib source
#  askpass        1.2.1   2024-10-04 [1] CRAN (R 4.5.2)
#  cli            3.6.5   2025-04-23 [1] CRAN (R 4.5.1)
#  credentials    2.0.3   2025-09-12 [1] CRAN (R 4.5.2)
#  curl           7.0.0   2025-08-19 [1] CRAN (R 4.5.2)
#  dplyr          1.2.0   2026-02-03 [1] CRAN (R 4.5.2)
#  farver         2.1.2   2024-05-13 [1] CRAN (R 4.5.1)
#  fs             1.6.6   2025-04-12 [1] CRAN (R 4.5.1)
#  generics       0.1.4   2025-05-09 [1] CRAN (R 4.5.2)
#  gert           2.3.1   2026-01-11 [1] CRAN (R 4.5.2)
#  ggplot2      * 4.0.2   2026-02-03 [1] CRAN (R 4.5.2)
#  gh             1.5.0   2025-05-26 [1] CRAN (R 4.5.2)
#  gitcreds       0.1.2   2022-09-08 [1] CRAN (R 4.5.2)
#  glue           1.8.0   2024-09-30 [1] CRAN (R 4.5.1)
#  gtable         0.3.6   2024-10-25 [1] CRAN (R 4.5.1)
#  here         * 1.0.2   2025-09-15 [1] CRAN (R 4.5.1)
#  httr2          1.2.2   2025-12-08 [1] CRAN (R 4.5.2)
#  jsonlite       2.0.0   2025-03-27 [1] CRAN (R 4.5.1)
#  labeling       0.4.3   2023-08-29 [1] CRAN (R 4.5.1)
#  lifecycle      1.0.5   2026-01-08 [1] CRAN (R 4.5.2)
#  magrittr       2.0.4   2025-09-12 [1] CRAN (R 4.5.1)
#  openssl        2.3.4   2025-09-30 [1] CRAN (R 4.5.2)
#  otel           0.2.0   2025-08-29 [1] CRAN (R 4.5.2)
#  pillar         1.11.1  2025-09-17 [1] CRAN (R 4.5.1)
#  pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.5.1)
#  purrr          1.2.1   2026-01-09 [1] CRAN (R 4.5.2)
#  R6             2.6.1   2025-02-15 [1] CRAN (R 4.5.1)
#  rappdirs       0.3.4   2026-01-17 [1] CRAN (R 4.5.2)
#  RColorBrewer   1.1-3   2022-04-03 [1] CRAN (R 4.5.1)
#  rlang          1.1.7   2026-01-09 [1] CRAN (R 4.5.2)
#  rprojroot      2.1.1   2025-08-26 [1] CRAN (R 4.5.1)
#  rstudioapi     0.18.0  2026-01-16 [1] CRAN (R 4.5.2)
#  S7             0.2.1   2025-11-14 [1] CRAN (R 4.5.2)
#  scales         1.4.0   2025-04-24 [1] CRAN (R 4.5.1)
#  sessioninfo  * 1.2.3   2025-02-05 [1] CRAN (R 4.5.2)
#  sys            3.4.3   2024-10-04 [1] CRAN (R 4.5.2)
#  tibble         3.3.1   2026-01-11 [1] CRAN (R 4.5.2)
#  tidyselect     1.2.1   2024-03-11 [1] CRAN (R 4.5.2)
#  usethis        3.2.1   2025-09-06 [1] CRAN (R 4.5.2)
#  vctrs          0.7.1   2026-01-23 [1] CRAN (R 4.5.2)
#  withr          3.0.2   2024-10-28 [1] CRAN (R 4.5.1)

#  [1] /home/pablosm/R/x86_64-pc-linux-gnu-library/4.5
#  [2] /usr/local/lib/R/site-library
#  [3] /usr/lib/R/site-library
#  [4] /usr/lib/R/library
#  * ── Packages attached to the search path.

# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
