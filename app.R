# app.R (MODIFIED: Removed VT/AFlutter, Integrated Vignettes)

# --- Load necessary libraries ---
library(shiny)
library(shinythemes)
library(ggplot2)
library(dplyr)

# Helper function for clamping values
clamp <- function(x, min_val, max_val) {
  max(min_val, min(x, max_val))
}

# --- ECG Wave Generation Functions ---
generate_segment <- function(time_points, amplitude, duration, center, width_factor = 4) {
  if (is.na(center) || duration <= 0) return(numeric(length(time_points)))
  sigma_equiv <- duration / width_factor
  if (sigma_equiv == 0) return(ifelse(time_points == center, amplitude, 0))
  scaled_time <- (time_points - center) / sigma_equiv
  wave <- amplitude * exp(-0.5 * scaled_time^2)
  return(wave)
}

generate_p_wave <- function(time_vector, type = "normal", p_start_time = 0.1, p_duration = 0.08, fs = 250, is_pvc = FALSE, atrial_rate_for_chb = 70) {
  if (type %in% c("afib", "svt") || is_pvc) return(numeric(length(time_vector))) # Removed "aflutter", "vt"
  
  p_wave <- numeric(length(time_vector)); p_amplitude_normal <- 1.5
  
  if (type == "chb"){
    p_cycle_duration_chb <- 60 / atrial_rate_for_chb
    num_p_waves <- ceiling(max(time_vector, na.rm = TRUE) / p_cycle_duration_chb)
    if(is.infinite(num_p_waves) || is.na(num_p_waves)) num_p_waves <- 0
    
    for(k in 0:(num_p_waves-1)){
      p_center_chb <- (k * p_cycle_duration_chb) + (p_duration / 2) + 0.01
      if(!is.na(max(time_vector)) && p_center_chb - p_duration/2 < max(time_vector, na.rm=TRUE) && p_center_chb + p_duration/2 > min(time_vector, na.rm=TRUE)){
        single_p <- generate_segment(time_vector, p_amplitude_normal, p_duration, p_center_chb, width_factor = 4)
        p_wave <- p_wave + single_p
      }
    }
    return(p_wave)
  }
  
  p_effective_start_time <- p_start_time; p_effective_duration <- p_duration; p_effective_amplitude <- p_amplitude_normal
  if (is.na(p_effective_start_time) || is.na(p_effective_duration)) return(numeric(length(time_vector)))
  
  if (type == "normal") {
    p_wave_segment <- generate_segment(time_vector, p_effective_amplitude, p_effective_duration, p_effective_start_time + p_effective_duration / 2, width_factor = 4)
    p_wave <- p_wave_segment
  } else if (type == "rae") {
    p_effective_amplitude <- 3.0
    p_wave_segment <- generate_segment(time_vector, p_effective_amplitude, p_effective_duration, p_effective_start_time + p_effective_duration / 2, width_factor = 5.5)
    p_wave <- p_wave_segment
  } else if (type == "lae") {
    p_effective_amplitude <- p_amplitude_normal * 1.1
    center1 <- p_effective_start_time + p_effective_duration * 0.30; center2 <- p_effective_start_time + p_effective_duration * 0.70
    duration_component <- p_effective_duration * 0.5
    p1 <- generate_segment(time_vector, p_effective_amplitude * 0.85, duration_component, center1, width_factor = 4.5)
    p2 <- generate_segment(time_vector, p_effective_amplitude * 1.0, duration_component, center2, width_factor = 4.5)
    p_wave <- p1 + p2
  } else if (type == "bae") {
    p_effective_amplitude <- 3.0
    center1 <- p_effective_start_time + p_effective_duration * 0.30; center2 <- p_effective_start_time + p_effective_duration * 0.70
    duration_component <- p_effective_duration * 0.5
    p1 <- generate_segment(time_vector, p_effective_amplitude * 1.0, duration_component, center1, width_factor = 5)
    p2 <- generate_segment(time_vector, p_effective_amplitude * 0.85, duration_component, center2, width_factor = 4.5)
    p_wave <- p1 + p2
  }
  
  if (!type %in% c("chb")) {
    if (!is.na(p_effective_start_time) && !is.na(p_effective_duration)) {
      p_wave[time_vector < p_effective_start_time | time_vector > (p_effective_start_time + p_effective_duration)] <- 0
    } else {
      p_wave <- numeric(length(time_vector))
    }
  }
  return(p_wave)
}

generate_qrs <- function(time_vector, qrs_start_time = 0.18, qrs_duration = 0.08, fs = 250, is_pvc = FALSE, rhythm_type = "normal", pvc_amplitude_scale = 1.5) {
  if(is.na(qrs_start_time) || is.na(qrs_duration) || qrs_duration <= 0) return(numeric(length(time_vector)))
  
  q_amp <- -1.0; r_amp <- 8.0; s_amp <- -2.0
  is_wide_qrs_rhythm = FALSE # rhythm_type == "vt" was here, VT is removed
  
  q_dur_prop <- 0.25
  r_dur_prop <- 0.40
  s_dur_prop <- 0.35
  
  if(is_pvc) { # Only PVCs make QRS wider/scaled now
    q_amp <- q_amp * pvc_amplitude_scale
    r_amp <- r_amp * pvc_amplitude_scale
    s_amp <- s_amp * pvc_amplitude_scale
    q_dur_prop <- 0.20
    r_dur_prop <- 0.50
    s_dur_prop <- 0.30
  }
  
  q_dur <- qrs_duration * q_dur_prop
  r_dur <- qrs_duration * r_dur_prop
  s_dur <- qrs_duration * s_dur_prop
  
  q_peak_time <- qrs_start_time + q_dur / 2
  r_peak_time <- qrs_start_time + q_dur + (r_dur / 2)
  s_peak_time <- qrs_start_time + q_dur + r_dur + (s_dur / 2)
  
  width_factor_q = if(is_pvc) 3.0 else 4.0 # Removed is_wide_qrs_rhythm
  width_factor_r = if(is_pvc) 3.5 else 4.5
  width_factor_s = if(is_pvc) 3.0 else 4.0
  
  q_wave <- generate_segment(time_vector, q_amp, q_dur, q_peak_time, width_factor = width_factor_q)
  r_wave <- generate_segment(time_vector, r_amp, r_dur, r_peak_time, width_factor = width_factor_r)
  s_wave <- generate_segment(time_vector, s_amp, s_dur, s_peak_time, width_factor = width_factor_s)
  qrs_wave_combined <- q_wave + r_wave + s_wave
  
  qrs_wave_combined[time_vector < qrs_start_time | time_vector > (qrs_start_time + qrs_duration)] <- 0
  return(qrs_wave_combined)
}

generate_t_wave <- function(time_vector, t_start_time = 0.30, t_duration = 0.16, fs = 250, is_pvc = FALSE, rhythm_type = "normal") {
  if(is.na(t_start_time) || is.na(t_duration)) return(numeric(length(time_vector)))
  t_amplitude <- 2.0
  if(is_pvc) { t_amplitude <- -t_amplitude } # Removed rhythm_type == "vt"
  t_wave_segment <- generate_segment(time_vector, t_amplitude, t_duration, t_start_time + t_duration / 2, width_factor = 3.5)
  t_wave_segment[time_vector < t_start_time | time_vector > (t_start_time + t_duration)] <- 0
  return(t_wave_segment)
}

generate_f_waves <- function(time_vector, fs, amplitude_scale = 0.25) {
  num_f_components <- sample(3:5, 1); f_waves <- numeric(length(time_vector))
  for (i in 1:num_f_components) {
    freq <- runif(1, 5, 10); phase <- runif(1, 0, 2*pi); amp <- amplitude_scale * runif(1, 0.3, 1.0)
    f_waves <- f_waves + amp * sin(2 * pi * freq * time_vector + phase)
  }
  f_waves <- f_waves + rnorm(length(time_vector), 0, amplitude_scale * 0.1)
  return(f_waves)
}

# generate_flutter_waves function is no longer needed as Atrial Flutter is removed.
# If it were needed by another part, it would be kept, but it's specific to AFlutter.

# --- Clinical Vignettes Data ---
vignettes_data <- list( # vt1 vignette removed
  list(
    id = "svt1", title = "Case 1: Sudden Racing Heart",
    scenario = "A 28-year-old female presents to the ED with a sudden onset of a 'racing heart', lightheadedness, and mild shortness of breath that started an hour ago while she was resting. She denies chest pain. No significant past medical history.",
    ecg_params = list(p_wave_type = "svt", heart_rate = 180, qrs_duration_ms = 80, st_segment_mm = 0, pr_interval_ms = NA),
    key_findings = "Regular, narrow complex tachycardia at ~180 bpm. P-waves are not clearly discernible (may be buried in QRS or retrograde). QRS duration is normal.",
    diagnosis_management = "Supraventricular Tachycardia (SVT), likely AVNRT or AVRT. Management may include vagal maneuvers, adenosine if stable, or synchronized cardioversion if unstable."
  ),
  list(
    id = "chb1", title = "Case 3: Bradycardia and Fatigue", # Renamed from Case 3 for sequential numbering if desired, but ID kept
    scenario = "An 80-year-old female reports several weeks of increasing fatigue, lightheadedness, and exertional dyspnea. Her pulse is noted to be very slow on examination.",
    ecg_params = list(p_wave_type = "chb", heart_rate = 70, qrs_duration_ms = 100, st_segment_mm = 0, pr_interval_ms = NA),
    key_findings = "Atrial activity (P-waves) present at a regular rate (e.g., ~70 bpm). Ventricular activity (QRS complexes) also regular but at a slower, independent rate (e.g., ~30-40 bpm). There is complete AV dissociation: P-waves 'march through' QRS complexes with no consistent PR interval.",
    diagnosis_management = "Third-Degree (Complete) AV Block. This patient is symptomatic and requires pacemaker implantation."
  )
)
vignette_choices <- setNames(lapply(vignettes_data, `[[`, "id"), lapply(vignettes_data, `[[`, "title"))
vignette_choices <- c(c("Select a case..." = ""), vignette_choices)

# --- UI Definition ---
ui <- fluidPage(
  theme = shinytheme("spacelab"),
  titlePanel("Comprehensive Interactive ECG Simulator"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("Global Controls"),
      sliderInput("heart_rate", "Heart Rate / Atrial Rate (bpm):", min = 30, max = 350, value = 75, step = 5),
      radioButtons("paper_speed", "Paper Speed:", choices = c("25 mm/s" = 25, "50 mm/s" = 50), selected = 25, inline = TRUE),
      radioButtons("ecg_gain", "ECG Gain:", choices = c("5 mm/mV (0.5x)" = 0.5, "10 mm/mV (1x)" = 1.0, "20 mm/mV (2x)" = 2.0), selected = 1.0, inline = TRUE),
      
      hr(),
      h4("Rhythm & Morphology"),
      selectInput("p_wave_type", "Base Rhythm / Atrial Activity:",
                  choices = c("Normal Sinus Rhythm" = "normal", # AFlutter and VT removed
                              "Right Atrial Enlargement" = "rae",
                              "Left Atrial Enlargement" = "lae",
                              "Biatrial Enlargement" = "bae",
                              "Atrial Fibrillation" = "afib",
                              "Supraventricular Tachycardia (SVT)" = "svt",
                              "Third-Degree AV Block (CHB)" = "chb"),
                  selected = "normal"),
      
      # Conditional panel for AV Conduction (AFlutter) removed
      conditionalPanel(
        condition = "input.p_wave_type == 'normal' || input.p_wave_type == 'rae' || input.p_wave_type == 'lae' || input.p_wave_type == 'bae'",
        sliderInput("pr_interval_ms", "PR Interval (ms):", min = 100, max = 400, value = 160, step = 10)
      ),
      sliderInput("qrs_duration_ms", "QRS Duration (ms):", min = 60, max = 240, value = 80, step = 10),
      sliderInput("st_segment_mm", "ST Segment Shift (mm):", min = -4, max = 4, value = 0, step = 0.5),
      
      hr(),
      h4("Ectopy & Artifacts"),
      checkboxInput("insert_pvcs", "Insert Occasional PVCs", value = FALSE),
      checkboxInput("add_baseline_wander", "Add Baseline Wander", value = FALSE),
      
      hr(),
      h5("Simulated Intervals / Info:"),
      uiOutput("interval_info_display"),
      hr(),
      h5("Caliper Measurements:"),
      uiOutput("caliper_output"),
      actionButton("clear_calipers", "Clear Calipers", class = "btn-sm btn-warning"),
      hr(),
      p(strong("Note:"), "Simplified Lead II. For educational purposes."),
      p(paste0("Version MODIFIED | ", format(Sys.Date(), "%Y"))) # Updated version note
    ),
    mainPanel(
      width = 8,
      # TabsetPanel removed, content integrated directly
      h3(textOutput("plot_title_text")),
      plotOutput("ecg_plot", height = "450px", click = "plot_click"),
      hr(),
      uiOutput("p_wave_info"),
      hr(), # Added separator
      h4("Clinical Case Vignettes"), # Moved from tab
      selectInput("selected_vignette", "Select Clinical Case:", choices = vignette_choices),
      uiOutput("vignette_display"),
      br(),
      actionButton("show_vignette_answer", "Show Interpretation & Diagnosis"),
      uiOutput("vignette_answer_display")
    )
  )
)

# --- Server Logic ---
server <- function(input, output, session) {
  
  rv <- reactiveValues(caliper_points = list(), current_vignette_id = NULL, show_answer = FALSE)
  
  observeEvent(input$clear_calipers, { rv$caliper_points <- list() })
  observeEvent(input$p_wave_type, { rv$caliper_points <- list() }) # Clear calipers on rhythm change
  observeEvent(input$heart_rate, { rv$caliper_points <- list() }) # Clear calipers on HR change
  
  observeEvent(input$plot_click, {
    current_points <- rv$caliper_points
    if (length(current_points) >= 2) { current_points <- list() }
    new_point <- data.frame(x = input$plot_click$x, y = input$plot_click$y)
    rv$caliper_points <- c(current_points, list(new_point))
  })
  
  output$caliper_output <- renderUI({
    points <- rv$caliper_points
    if (length(points) == 2) {
      p1 <- points[[1]]; p2 <- points[[2]]; dX <- abs(p2$x - p1$x); dY <- p2$y - p1$y
      HTML(sprintf("<b>dX (Time):</b> %.3f s<br/><b>dY (Amplitude):</b> %.2f mm", dX, dY))
    } else if (length(points) == 1) { HTML("First point selected. Click a second point.") }
    else { HTML("Click on ECG to place caliper points (max 2).") }
  })
  
  observeEvent(input$selected_vignette, {
    rv$current_vignette_id <- input$selected_vignette
    rv$show_answer <- FALSE
    if (input$selected_vignette != "") {
      vignette <- Filter(function(v) v$id == input$selected_vignette, vignettes_data)[[1]]
      
      updateSelectInput(session, "p_wave_type", selected = vignette$ecg_params$p_wave_type)
      updateSliderInput(session, "heart_rate", value = vignette$ecg_params$heart_rate)
      updateSliderInput(session, "qrs_duration_ms", value = vignette$ecg_params$qrs_duration_ms)
      updateSliderInput(session, "st_segment_mm", value = vignette$ecg_params$st_segment_mm)
      
      pr_applicable_rhythms <- c("normal", "rae", "lae", "bae")
      if (vignette$ecg_params$p_wave_type %in% pr_applicable_rhythms) {
        updateSliderInput(session, "pr_interval_ms", value = ifelse(is.null(vignette$ecg_params$pr_interval_ms) || is.na(vignette$ecg_params$pr_interval_ms), 160, vignette$ecg_params$pr_interval_ms) )
      }
      # Removed av_conduction update as Atrial Flutter is removed
      
      updateCheckboxInput(session, "insert_pvcs", value = if(!is.null(vignette$ecg_params$insert_pvcs)) vignette$ecg_params$insert_pvcs else FALSE)
      updateCheckboxInput(session, "add_baseline_wander", value = if(!is.null(vignette$ecg_params$add_baseline_wander)) vignette$ecg_params$add_baseline_wander else FALSE)
    }
  })
  
  output$vignette_display <- renderUI({
    req(rv$current_vignette_id)
    if (rv$current_vignette_id != "") {
      vignette <- Filter(function(v) v$id == rv$current_vignette_id, vignettes_data)[[1]]
      tagList(h5(vignette$title), p(vignette$scenario)) # Changed h4 to h5 for visual hierarchy
    }
  })
  
  observeEvent(input$show_vignette_answer, { rv$show_answer <- TRUE })
  
  output$vignette_answer_display <- renderUI({
    req(rv$current_vignette_id, rv$show_answer)
    if (rv$current_vignette_id != "" && rv$show_answer) {
      vignette <- Filter(function(v) v$id == rv$current_vignette_id, vignettes_data)[[1]]
      tagList(h5("Key ECG Findings:"), p(vignette$key_findings), h5("Diagnosis & Management Pearls:"), p(vignette$diagnosis_management))
    }
  })
  
  ecg_base_params <- reactive({
    list(fs = 250, st_segment_duration_morph = 0.08, t_duration_morph = 0.16,
         ventricular_escape_rate_chb = 40)
  })
  
  current_rhythm_params <- reactive({
    base_p <- ecg_base_params(); hr <- input$heart_rate; p_type <- input$p_wave_type
    
    req(input$qrs_duration_ms, input$st_segment_mm)
    
    qrs_dur_morph <- input$qrs_duration_ms / 1000
    
    pr_interval_actual_ms <- NA
    pr_applicable_rhythms <- c("normal", "rae", "lae", "bae")
    if (p_type %in% pr_applicable_rhythms) {
      req(input$pr_interval_ms)
      pr_interval_actual_ms <- input$pr_interval_ms
    }
    pr_interval_setting <- if (is.na(pr_interval_actual_ms)) NA else pr_interval_actual_ms / 1000
    
    p_start_time <- 0.04
    p_duration_selected <- switch(p_type, "normal" = 0.08, "rae" = 0.08, "lae" = 0.12, "bae" = 0.12, 0) # SVT will have 0
    qrs_start_time_regular <- if(!is.na(pr_interval_setting)) p_start_time + pr_interval_setting else NA
    pr_calc_display <- if(is.na(pr_interval_setting)) NA else pr_interval_setting
    qrs_calc_display <- qrs_dur_morph
    # Removed VT specific QRS widening here, as VT is removed. QRS width is now purely by slider.
    
    t_start_regular <- if(!is.na(qrs_start_time_regular)) qrs_start_time_regular + qrs_calc_display + base_p$st_segment_duration_morph else NA
    t_end_regular <- if(!is.na(t_start_regular)) t_start_regular + base_p$t_duration_morph else NA
    qt_calc_display <- if(is.na(t_end_regular) || is.na(qrs_start_time_regular) || p_type %in% c("chb")) NA else (t_end_regular - qrs_start_time_regular) # Removed "vt"
    # av_conduction_ratio removed
    atrial_rate_for_chb <- if(p_type == "chb") hr else 75 # Simplified from atrial_rate_for_chb_and_flutter
    
    list(
      fs = base_p$fs, heart_rate = hr, p_wave_type = p_type,
      paper_speed = as.numeric(input$paper_speed), ecg_gain = as.numeric(input$ecg_gain),
      beat_duration_seconds = 60 / hr,
      p_start_time_reg = p_start_time, p_duration_reg = p_duration_selected,
      qrs_start_time_reg = qrs_start_time_regular, t_start_time_reg = t_start_regular,
      qrs_duration_morph = qrs_calc_display,
      st_segment_duration_morph = base_p$st_segment_duration_morph,
      t_duration_morph = base_p$t_duration_morph, st_segment_shift_mm = input$st_segment_mm,
      pr_interval_calc_display = pr_calc_display, qrs_duration_calc_display = qrs_calc_display, qt_interval_calc_display = qt_calc_display,
      insert_pvcs = input$insert_pvcs && !p_type %in% c("afib", "svt", "chb"), # Removed "aflutter", "vt"
      add_baseline_wander = input$add_baseline_wander,
      # av_conduction_ratio = NA, # Removed
      atrial_rate_for_chb = atrial_rate_for_chb,
      ventricular_escape_rate_chb = base_p$ventricular_escape_rate_chb
    )
  })
  
  generated_ecg_trace <- reactive({
    params <- current_rhythm_params()
    base_display_duration <- 5.0
    total_duration_to_generate <- base_display_duration * (25 / params$paper_speed)
    full_time_vector <- seq(0, total_duration_to_generate, by = 1/params$fs)
    final_ecg_signal <- numeric(length(full_time_vector))
    st_shift_overlay <- numeric(length(full_time_vector))
    
    if (params$p_wave_type == "afib") {
      # ... (afib logic remains the same) ...
      final_ecg_signal <- generate_f_waves(full_time_vector, fs=params$fs, amplitude_scale = 0.25)
      current_t <- 0; avg_rr <- 60 / params$heart_rate
      qrs_t_overlay <- numeric(length(full_time_vector))
      while(current_t < total_duration_to_generate) {
        rr <- rnorm(1, avg_rr, avg_rr * 0.20); rr <- clamp(rr, 0.33, avg_rr * 1.8)
        qrs_onset_abs <- current_t
        qrs_dur = params$qrs_duration_morph; st_dur = params$st_segment_duration_morph; t_dur = params$t_duration_morph
        single_beat_shape_duration <- qrs_dur + st_dur + t_dur + 0.05
        time_vec_shape <- seq(0, single_beat_shape_duration, by = 1/params$fs)
        qrs_shape <- generate_qrs(time_vec_shape, 0, qrs_dur, params$fs, rhythm_type = params$p_wave_type)
        t_shape_start <- qrs_dur + st_dur
        t_shape <- generate_t_wave(time_vec_shape, t_shape_start, t_dur, params$fs, rhythm_type = params$p_wave_type)
        beat_shape <- qrs_shape + t_shape
        st_start_abs = qrs_onset_abs + qrs_dur; st_end_abs = qrs_onset_abs + qrs_dur + st_dur
        st_start_idx = findInterval(st_start_abs, full_time_vector); st_end_idx = findInterval(st_end_abs, full_time_vector) -1
        if(params$st_segment_shift_mm != 0 && st_start_idx > 0 && st_end_idx >= st_start_idx && st_end_idx <= length(st_shift_overlay)){ st_shift_overlay[st_start_idx:st_end_idx] <- params$st_segment_shift_mm }
        start_idx <- findInterval(qrs_onset_abs, full_time_vector); end_idx <- start_idx + length(beat_shape) - 1
        if(start_idx > 0 && end_idx <= length(full_time_vector)) { qrs_t_overlay[start_idx:end_idx] <- qrs_t_overlay[start_idx:end_idx] + beat_shape }
        current_t <- current_t + rr
      }
      final_ecg_signal <- final_ecg_signal + qrs_t_overlay
      # else if (params$p_wave_type == "aflutter") removed
    } else if (params$p_wave_type == "chb") {
      final_ecg_signal <- generate_p_wave(full_time_vector, type = "chb", p_duration = 0.08, fs = params$fs, atrial_rate_for_chb = params$atrial_rate_for_chb) # use updated param name
      qrs_t_overlay <- numeric(length(full_time_vector))
      ventricular_rr <- 60 / params$ventricular_escape_rate_chb
      current_t <- runif(1, 0, ventricular_rr*0.5)
      while(current_t < total_duration_to_generate) {
        # ... (CHB QRS/T overlay logic remains the same) ...
        qrs_onset_abs <- current_t
        qrs_dur = params$qrs_duration_morph; st_dur = params$st_segment_duration_morph; t_dur = params$t_duration_morph
        single_beat_shape_duration <- qrs_dur + st_dur + t_dur + 0.05
        time_vec_shape <- seq(0, single_beat_shape_duration, by = 1/params$fs)
        qrs_shape <- generate_qrs(time_vec_shape, 0, qrs_dur, params$fs, rhythm_type = params$p_wave_type)
        t_shape_start <- qrs_dur + st_dur
        t_shape <- generate_t_wave(time_vec_shape, t_shape_start, t_dur, params$fs, rhythm_type = params$p_wave_type)
        beat_shape <- qrs_shape + t_shape
        st_start_abs = qrs_onset_abs + qrs_dur; st_end_abs = qrs_onset_abs + qrs_dur + st_dur
        st_start_idx = findInterval(st_start_abs, full_time_vector); st_end_idx = findInterval(st_end_abs, full_time_vector) -1
        if(params$st_segment_shift_mm != 0 && st_start_idx > 0 && st_end_idx >= st_start_idx && st_end_idx <= length(st_shift_overlay)){ st_shift_overlay[st_start_idx:st_end_idx] <- params$st_segment_shift_mm }
        start_idx <- findInterval(qrs_onset_abs, full_time_vector); end_idx <- start_idx + length(beat_shape) - 1
        if(start_idx > 0 && end_idx <= length(full_time_vector)) { qrs_t_overlay[start_idx:end_idx] <- qrs_t_overlay[start_idx:end_idx] + beat_shape }
        current_t <- current_t + ventricular_rr
      }
      final_ecg_signal <- final_ecg_signal + qrs_t_overlay
    } else { # Regular rhythms (NSR, RAE, LAE, BAE, SVT) - VT removed
      is_tachy = params$p_wave_type == "svt" # VT removed from this check
      num_beats_to_generate <- ceiling(total_duration_to_generate / params$beat_duration_seconds) + 1
      for(i in 0:(num_beats_to_generate - 1)) {
        # ... (Regular beat generation logic remains largely the same) ...
        is_pvc_beat <- FALSE
        if (params$insert_pvcs && i > 0 && runif(1) < 0.1) { is_pvc_beat <- TRUE }
        time_vector_single_beat <- seq(0, params$beat_duration_seconds, by = 1/params$fs)
        
        p_wave_sgl <- generate_p_wave(time_vector_single_beat, type = params$p_wave_type,
                                      p_start_time = params$p_start_time_reg,
                                      p_duration = params$p_duration_reg,
                                      fs = params$fs, is_pvc = is_pvc_beat)
        
        qrs_sgl_start_time <- if(is_pvc_beat) params$p_start_time_reg else (if(is_tachy && !is.na(params$p_start_time_reg)) params$p_start_time_reg else params$qrs_start_time_reg )
        if(is.na(qrs_sgl_start_time) && is_tachy) qrs_sgl_start_time <- 0.01 
        qrs_sgl_duration <- if(is_pvc_beat) 0.14 else params$qrs_duration_morph
        
        if(is.na(qrs_sgl_start_time)) next
        
        qrs_sgl <- generate_qrs(time_vector_single_beat, qrs_start_time = qrs_sgl_start_time,
                                qrs_duration = qrs_sgl_duration, fs = params$fs, is_pvc = is_pvc_beat, rhythm_type = params$p_wave_type)
        
        t_sgl_start_time <- qrs_sgl_start_time + qrs_sgl_duration + (if(is_pvc_beat || is_tachy) 0.05 else params$st_segment_duration_morph)
        t_sgl <- generate_t_wave(time_vector_single_beat, t_start_time = t_sgl_start_time,
                                 t_duration = params$t_duration_morph, fs = params$fs, is_pvc = is_pvc_beat, rhythm_type = params$p_wave_type)
        
        single_beat_ecg_signal <- p_wave_sgl + qrs_sgl + t_sgl
        
        st_abs_start_this_beat = (i * params$beat_duration_seconds) + qrs_sgl_start_time + qrs_sgl_duration
        st_abs_end_this_beat = (i * params$beat_duration_seconds) + t_sgl_start_time
        st_start_idx = findInterval(st_abs_start_this_beat, full_time_vector)
        st_end_idx = findInterval(st_abs_end_this_beat, full_time_vector) -1
        if(params$st_segment_shift_mm != 0 && st_start_idx > 0 && st_end_idx >= st_start_idx && st_end_idx <= length(st_shift_overlay) && !is.na(st_start_idx) && !is.na(st_end_idx)){
          st_shift_overlay[st_start_idx:st_end_idx] <- params$st_segment_shift_mm
        }
        
        beat_start_time_abs <- i * params$beat_duration_seconds
        start_idx <- findInterval(beat_start_time_abs, full_time_vector)
        end_idx <- start_idx + length(single_beat_ecg_signal) - 1
        if(start_idx > 0 && end_idx <= length(full_time_vector) && !is.na(start_idx) && !is.na(end_idx)) {
          final_ecg_signal[start_idx:end_idx] <- final_ecg_signal[start_idx:end_idx] + single_beat_ecg_signal
        }
      }
    }
    final_ecg_signal <- final_ecg_signal + st_shift_overlay
    if(params$add_baseline_wander){
      wander_freq = runif(1, 0.1, 0.4); wander_amp = runif(1, 0.5, 1.0)
      baseline_wander_signal = wander_amp * sin(2 * pi * wander_freq * full_time_vector + runif(1, 0, 2*pi))
      final_ecg_signal = final_ecg_signal + baseline_wander_signal
    }
    final_ecg_signal <- final_ecg_signal * params$ecg_gain
    data.frame(Time = full_time_vector, Value = final_ecg_signal) %>% filter(Time <= total_duration_to_generate)
  })
  
  output$plot_title_text <- renderText({
    params <- current_rhythm_params()
    rhythm_name_map <- c("normal" = "Normal Sinus", "rae" = "RAE", "lae" = "LAE", "bae" = "BAE",
                         "afib" = "Atrial Fibrillation", # AFlutter, VT removed
                         "svt" = "Supraventricular Tachycardia",
                         "chb" = "Third-Degree AV Block")
    rhythm_name <- rhythm_name_map[params$p_wave_type]
    hr_display <- params$heart_rate
    hr_unit_note <- "bpm"
    # AFlutter specific title logic removed
    if (params$p_wave_type == "afib") {
      hr_unit_note <- "bpm (Avg Vent. Rate)"
    } else if (params$p_wave_type == "chb") {
      hr_unit_note <- paste0("bpm (Atrial Rate), Vent. Rate ~", params$ventricular_escape_rate_chb, " bpm")
    } else if (params$p_wave_type == "svt"){ # VT removed
      hr_unit_note <- "bpm (Vent. Rate)"
    }
    title <- paste0(rhythm_name, " @ ", hr_display, " ", hr_unit_note)
    if(params$insert_pvcs) { title <- paste0(title, " with PVCs") }
    if(params$add_baseline_wander) { title <- paste0(title, " (Wandering Baseline)")}
    title
  })
  
  output$interval_info_display <- renderUI({
    params <- current_rhythm_params()
    pr_text <- if(is.na(params$pr_interval_calc_display)) "<b>PR Interval:</b> N/A" else sprintf("<b>PR Interval:</b> %.0f ms", params$pr_interval_calc_display * 1000)
    qrs_text <- sprintf("<b>QRS Duration:</b> %.0f ms", params$qrs_duration_calc_display * 1000)
    qt_text <- if(is.na(params$qt_interval_calc_display)) "<b>QT Interval:</b> Variable / N/A" else sprintf("<b>QT Interval:</b> %.0f ms (approx)", params$qt_interval_calc_display * 1000)
    rhythm_type_text <- params$p_wave_type
    rhythm_desc <- "Regular"
    hr_label <- "Heart Rate:"
    hr_value <- params$heart_rate
    if(rhythm_type_text == "afib") {
      rhythm_desc <- "Irregularly irregular"; hr_label <- "Avg Vent. Rate:"
      # else if(rhythm_type_text == "aflutter") removed
    } else if (rhythm_type_text == "chb") {
      rhythm_desc <- "AV Dissociation"; hr_label <- "Atrial Rate:"
    } else if (rhythm_type_text == "svt") { hr_label <- "Vent. Rate:" } # VT removed
    if(params$insert_pvcs && !rhythm_type_text %in% c("afib", "svt", "chb")) rhythm_desc <- paste0(rhythm_desc, " (with PVCs)") # Aflutter, VT removed
    HTML(paste(
      sprintf("<b>Rhythm:</b> %s", rhythm_desc),
      sprintf("<b>%s</b> %.0f bpm", hr_label, hr_value),
      if(rhythm_type_text == "chb") sprintf("<b>Vent. Escape Rate:</b> ~%.0f bpm", params$ventricular_escape_rate_chb) else "",
      pr_text, qrs_text, qt_text,
      sprintf("<b>ST Shift:</b> %.1f mm", params$st_segment_shift_mm),
      sep = "<br/>"))
  })
  
  output$ecg_plot <- renderPlot({
    req(input$p_wave_type, input$heart_rate, input$paper_speed, input$ecg_gain,
        input$qrs_duration_ms, input$st_segment_mm)
    pr_applicable_rhythms <- c("normal", "rae", "lae", "bae")
    if (input$p_wave_type %in% pr_applicable_rhythms) { req(input$pr_interval_ms) }
    # Removed req for input$av_conduction
    
    # ... (Plotting logic remains the same) ...
    df_display <- generated_ecg_trace()
    if (is.null(df_display) || nrow(df_display) == 0 || all(is.na(df_display$Value))) {
      return(ggplot() + theme_void() + annotate("text", x=0.5, y=0.5, label="ECG data generation error or no data.", size=5) + labs(title="Error"))
    }
    x_max_plot <- max(df_display$Time, na.rm=TRUE)
    min_val_plot <- min(df_display$Value, -4, na.rm = TRUE); max_val_plot <- max(df_display$Value, 15, na.rm = TRUE)
    plot_min_y <- floor(min_val_plot / 5) * 5 - 5; plot_max_y <- ceiling(max_val_plot / 5) * 5 + 5
    if (plot_max_y < 10 && max_val_plot < 10 && !is.na(max_val_plot)) plot_max_y <- 10;
    if (plot_min_y > -5 && min_val_plot > -5 && !is.na(min_val_plot)) plot_min_y <- -5
    
    p <- ggplot(df_display, aes(x = Time, y = Value)) +
      geom_hline(yintercept = seq(plot_min_y, plot_max_y, by = 1), color = "pink", linetype = "dotted", linewidth = 0.4) +
      geom_vline(xintercept = seq(0, x_max_plot, by = 0.04), color = "pink", linetype = "dotted", linewidth = 0.4) +
      geom_hline(yintercept = seq(plot_min_y, plot_max_y, by = 5), color = "lightcoral", linetype = "dashed", linewidth = 0.6) +
      geom_vline(xintercept = seq(0, x_max_plot, by = 0.20), color = "lightcoral", linetype = "dashed", linewidth = 0.6) +
      geom_line(color = "black", linewidth = 1.2, na.rm = TRUE) +
      labs(x = "Time (seconds)", y = paste0("Amplitude (mm) at ", input$ecg_gain, "x Gain"), title = NULL) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black"), axis.title = element_text(size = 12)) +
      coord_cartesian(ylim = c(plot_min_y, plot_max_y), xlim = c(0, x_max_plot), expand = FALSE) +
      scale_x_continuous(breaks = seq(0, x_max_plot, by = 0.20), minor_breaks = seq(0, x_max_plot, by = 0.04)) +
      scale_y_continuous(breaks = seq(plot_min_y, plot_max_y, by = 5), minor_breaks = seq(plot_min_y, plot_max_y, by = 1))
    
    cal_points_data <- rv$caliper_points
    if (length(cal_points_data) > 0) {
      df_cal_points <- do.call(rbind, lapply(cal_points_data, function(pt) data.frame(x=pt$x, y=pt$y)))
      p <- p + geom_point(data = df_cal_points, aes(x=x, y=y), color="blue", size=3, shape=4, stroke=1.5)
      if (length(cal_points_data) == 2) {
        p <- p + geom_segment(data=df_cal_points, aes(x=x[1], y=y[1], xend=x[2], yend=y[2]), color="blue", linewidth=1, linetype="dashed")
      }
    }
    print(p)
  })
  
  output$p_wave_info <- renderUI({
    type <- input$p_wave_type
    params <- current_rhythm_params()
    header_style <- "font-size: 1.1em; font-weight: bold; margin-bottom: 5px; color: #007bff;"
    list_style <- "margin-left: 20px; margin-bottom: 10px;"
    base_info <- switch(type,
                        "normal" = "<li><b>P-Wave Amplitude:</b> &lt; 2.5 mm</li><li><b>P-Wave Duration:</b> &lt; 0.12s</li><li><b>Morphology:</b> Smooth, rounded, upright P before each QRS.</li>",
                        "rae" = "<li><b>P-Wave Amplitude:</b> &gt; 2.5 mm (peaked 'P pulmonale')</li><li><b>P-Wave Duration:</b> Usually normal</li>",
                        "lae" = "<li><b>P-Wave Amplitude:</b> Often normal</li><li><b>P-Wave Duration:</b> &ge; 0.12s (prolonged)</li><li><b>Morphology:</b> Often notched P ('P mitrale').</li>",
                        "bae" = "<li><b>P-Wave Amplitude:</b> Often &gt; 2.5 mm</li><li><b>P-Wave Duration:</b> Often &ge; 0.12s</li><li><b>Morphology:</b> Features of both RAE & LAE.</li>",
                        "afib" = "<li><b>Atrial Activity:</b> No P-waves. Irregular fibrillatory (f) waves.</li><li><b>Ventricular Rhythm:</b> Irregularly irregular.</li>",
                        # "aflutter" case removed
                        "svt" = "<li><b>Atrial Activity:</b> P-waves often hidden or retrograde, not clearly visible before QRS.</li><li><b>Ventricular Rhythm:</b> Regular, rapid.</li><li><b>QRS:</b> Typically narrow.</li>",
                        # "vt" case removed
                        "chb" = "<li><b>Atrial Activity:</b> P-waves present and regular (atrial rate).</li><li><b>Ventricular Rhythm:</b> QRS complexes regular but at a slower, independent escape rate.</li><li><b>AV Relationship:</b> Complete AV dissociation (no consistent PR interval).</li>"
    )
    title_text_map <- c("normal" = "Normal Sinus Rhythm:", "rae" = "RAE:", "lae" = "LAE:", "bae" = "BAE:",
                        "afib" = "Atrial Fibrillation:", # AFlutter, VT removed
                        "svt" = "Supraventricular Tachycardia:",
                        "chb" = "Third-Degree AV Block:")
    title_text <- title_text_map[type]
    full_html <- paste0("<div style='", header_style, "'>", title_text, "</div><ul style='", list_style, "'>", base_info, "</ul>")
    if(!is.null(params$insert_pvcs) && params$insert_pvcs){ full_html <- paste0(full_html, "<div style='", header_style, "margin-top:10px;'>Occasional PVCs:</div><ul style='", list_style, "'><li><b>Morphology:</b> Wide QRS, no P-wave, discordant T-wave.</li></ul>") }
    if(!is.null(params$add_baseline_wander) && params$add_baseline_wander){ full_html <- paste0(full_html, "<div style='", header_style, "margin-top:10px;'>Baseline Wander:</div><ul style='", list_style, "'><li>Slow, undulating baseline shifts.</li></ul>") }
    if(!is.null(params$st_segment_shift_mm) && params$st_segment_shift_mm != 0) { full_html <- paste0(full_html, "<div style='", header_style, "margin-top:10px;'>ST Segment Shift:</div><ul style='", list_style, "'><li>Currently shifted by ", params$st_segment_shift_mm, " mm.</li></ul>") }
    HTML(full_html)
  })
}

# --- Run the application ---
shinyApp(ui, server)
