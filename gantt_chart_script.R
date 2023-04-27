library("ganttrify")
ganttrify(project = dissertation_timeline,
          by_date = TRUE,
          size_text_relative = 1.4, 
          project_start_date = "2023-05",
          mark_quarters = F,
          alpha_wp = 0.9,
          alpha_activity = 0.6,
          month_breaks = 1,
          font_family = "Roboto Condensed")
