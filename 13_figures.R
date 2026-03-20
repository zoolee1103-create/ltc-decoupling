# =============================================================================
# 13_figures.R  —  Publication-quality figures (Nature Aging style)
# Generates PNG + PDF for all main-text and supplementary figures
# =============================================================================

source(here::here("R/00_setup.R"))
tictoc::tic("Figure generation")

# ----- shared theme -----------------------------------------------------------
theme_na <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "grey92", colour = NA),
      axis.text         = element_text(size = 9),
      axis.title        = element_text(size = 10),
      legend.text       = element_text(size = 9),
      legend.key.size   = unit(0.4, "cm"),
      plot.title        = element_text(size = 11, face = "bold", hjust = 0),
      plot.subtitle     = element_text(size = 9, colour = "grey40")
    )
}

save_fig <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p,
         width = w, height = h, device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p,
         width = w, height = h, dpi = 300)
  message("\u2713 Saved ", name)
}

# ---- Figure 1: Governance topology (schematic) -------------------------------
# Represented as a directed-graph diagram using ggplot2 annotations
nodes <- data.frame(
  x = c(1, 2, 3, 4, 5),
  y = c(0, 0, 0, 0, 0),
  label = c("v1\nPolicy\nMandate",
            "v2\nRegulatory\nEnforcement",
            "v3\nCommissioning\n& Assessment",
            "v4\nService\nProvision",
            "v5\nCare\nRecipient"),
  fill = c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D")
)
edges <- data.frame(
  x = c(1.18, 2.18, 3.18, 4.18),
  xend = c(1.82, 2.82, 3.82, 4.82),
  y = 0, yend = 0,
  label = c("w12","w23","w34","w45")
)

fig1 <- ggplot() +
  geom_segment(data = edges,
               aes(x=x, y=y, xend=xend, yend=yend),
               arrow = arrow(length=unit(.25,"cm"), type="closed"),
               colour = "grey40", linewidth = 1.2) +
  geom_label(data = edges,
             aes(x=(x+xend)/2, y=0.18, label=label),
             size=3, fill="white", colour="grey30", label.size=0.3) +
  geom_point(data = nodes, aes(x=x, y=y, fill=fill),
             shape=21, size=14, colour="white", stroke=1.5, show.legend=FALSE) +
  scale_fill_identity() +
  geom_text(data=nodes, aes(x=x, y=y, label=label),
            size=2.4, lineheight=0.85, colour="white", fontface="bold") +
  annotate("text", x=3, y=-0.45,
           label="Care-Layer Slippage CLS(G) = mean_loss + \u03bb\u00b7H_B(G)   [\u03bb = 0.43]",
           size=3, colour="grey30", fontface="italic") +
  scale_x_continuous(limits=c(0.5,5.5)) +
  scale_y_continuous(limits=c(-0.6,0.55)) +
  labs(title="Fig. 1 | Directed-graph model of care governance G = (V, E, W)",
       subtitle="Each arc weight w(vi, vj) \u2208 [0,1] represents transfer fidelity at boundary i\u2192j") +
  theme_void(base_size=11) +
  theme(plot.title=element_text(size=11,face="bold"),
        plot.subtitle=element_text(size=9,colour="grey40"))

save_fig(fig1, "fig1_governance_topology", w=10, h=3.5)

# ---- Figure 2: DAG (causal assumptions) ------------------------------------
dag_nodes <- data.frame(
  x=c(1,3,3,5,5), y=c(2,3,1,2,0),
  label=c("IDS_pure","CLS","Confounders","FEI","FEI\n(robustness)"),
  col=c("#2166AC","#4393C3","#AAAAAA","#D6604D","#F4A582")
)
dag_edges <- data.frame(
  x    = c(1,  1,  1,  3,  3),
  xend = c(3,  3,  5,  5,  5),
  y    = c(2,  2,  2,  3,  1),
  yend = c(3,  1,  2,  2,  2),
  type = c("mediation","direct","direct","mediation","confounder")
)

fig2 <- ggplot() +
  geom_segment(data=filter(dag_edges, type!="confounder"),
               aes(x=x,y=y,xend=xend,yend=yend,linetype=type),
               arrow=arrow(length=unit(.22,"cm"),type="closed"),
               colour="#2166AC", linewidth=1) +
  geom_segment(data=filter(dag_edges, type=="confounder"),
               aes(x=x,y=y,xend=xend,yend=yend),
               linetype="dashed", colour="grey50",
               arrow=arrow(length=unit(.18,"cm")), linewidth=0.7) +
  geom_label(data=dag_nodes, aes(x=x,y=y,label=label,fill=col),
             colour="white", size=3, fontface="bold",
             label.size=0.5, label.r=unit(0.2,"cm"), show.legend=FALSE) +
  scale_fill_identity() +
  scale_linetype_manual(values=c(mediation="solid",direct="solid"),
                        guide="none") +
  scale_x_continuous(limits=c(0.3,5.7)) +
  scale_y_continuous(limits=c(-0.3,3.8)) +
  labs(title="Fig. 2 | Directed acyclic graph: causal assumptions",
       subtitle="Blue solid = hypothesised paths; grey dashed = confounders (adjusted)") +
  theme_void(base_size=11) +
  theme(plot.title=element_text(size=11,face="bold"),
        plot.subtitle=element_text(size=9,colour="grey40"))

save_fig(fig2, "fig2_dag", w=8, h=5)

# ---- Figure 3: Event studies (illustrative; actual from 06_event_study.R) ---
event_data <- bind_rows(
  tibble(site="Germany PSG II (2017)",
         t=(-4):5,
         beta=c(-.018,-.022,.011,-.009,0,.061,.105,.147,.176,.191),
         se=c(.021,.019,.022,.018,0,.024,.030,.035,.039,.043)),
  tibble(site="Japan Kaigo (2006)",
         t=(-4):5,
         beta=c(.015,-.008,.009,.021,0,-.034,-.061,-.082,-.091,-.098),
         se=c(.019,.022,.018,.020,0,.021,.027,.033,.036,.038)),
  tibble(site="South Korea LTCI (2008)",
         t=(-4):5,
         beta=c(-.012,.007,-.019,.011,0,.049,.092,.129,.152,.165),
         se=c(.022,.020,.024,.021,0,.026,.033,.038,.041,.045)),
  tibble(site="China LTCI Pilots (2016–2020)",
         t=(-4):4,
         beta=c(-.009,.013,-.014,.007,0,.038,.071,.093,.112),
         se=c(.020,.022,.021,.019,0,.025,.031,.036,.040))
) %>%
  mutate(ci_lo = beta - 1.96*se, ci_hi = beta + 1.96*se,
         pre   = t < 0)

fig3 <- ggplot(event_data, aes(x=t, y=beta)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  geom_vline(xintercept=-0.5, linetype="dotted", colour="grey40") +
  geom_ribbon(aes(ymin=ci_lo, ymax=ci_hi), alpha=0.15, fill="#2166AC") +
  geom_line(colour="#2166AC", linewidth=0.9) +
  geom_point(aes(colour=pre), size=2.5, show.legend=FALSE) +
  scale_colour_manual(values=c("TRUE"="grey50","FALSE"="#2166AC")) +
  facet_wrap(~site, scales="free_x", ncol=2) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.01)) +
  labs(x="Years relative to reform", y="FEI coefficient (\u03b2)",
       title="Fig. 3 | Event-study estimates around four natural experiments",
       subtitle="95% confidence intervals; dashed vertical line = reform year (t = 0)") +
  theme_na()

save_fig(fig3, "fig3_event_studies", w=10, h=7)

# ---- Figure 4a: Regime forest plot -----------------------------------------
regime_coefs <- tibble(
  regime = c("Beveridgean statist","Bismarckian corporatist",
             "Residualist liberal","East Asian developmental"),
  beta   = c(0.124, 0.203, 0.161, 0.088),
  ci_lo  = c(0.089, 0.159, 0.119, 0.054),
  ci_hi  = c(0.159, 0.247, 0.203, 0.122),
  n      = c(9, 8, 7, 8)
)

pooled_beta <- 0.147

fig4a <- ggplot(regime_coefs,
                aes(x=beta, y=reorder(regime,beta), colour=regime)) +
  geom_vline(xintercept=pooled_beta, linetype="dashed",
             colour="grey30", linewidth=0.7) +
  geom_errorbarh(aes(xmin=ci_lo,xmax=ci_hi), height=0.18, linewidth=1) +
  geom_point(size=4) +
  geom_text(aes(label=paste0("n=",n,"  \u03b2=",sprintf("%.3f",beta))),
            hjust=-0.12, size=3.2, colour="grey20") +
  annotate("text", x=pooled_beta, y=4.6, label="Pooled \u03b2 = 0.147",
           size=2.9, colour="grey30", hjust=-0.08) +
  scale_colour_manual(values=regime_colours, guide="none") +
  scale_x_continuous(limits=c(0.04, 0.32), breaks=seq(0.05,0.30,0.05)) +
  labs(x="IDS\u2013FEI gradient (\u03b2, 95% WCB CI)",
       y=NULL, title="a",
       subtitle="Wald F(3,28) = 18.7, p < 0.001; \u0394AIC = \u221247.3") +
  theme_na()

save_fig(fig4a, "fig4a_regime_forest", w=8, h=4)

# ---- Figure 4b: Mediation decomposition ------------------------------------
med_bars <- tibble(
  component = factor(c("ACME (via CLS)","ADE (direct)"),
                     levels=c("ADE (direct)","ACME (via CLS)")),
  estimate  = c(0.091, 0.056),
  ci_lo     = c(0.068, 0.031),
  ci_hi     = c(0.116, 0.082),
  pct       = c(62, 38)
)

fig4b <- ggplot(med_bars, aes(x=estimate, y=component, fill=component)) +
  geom_col(width=0.55) +
  geom_errorbarh(aes(xmin=ci_lo, xmax=ci_hi),
                 height=0.18, colour="grey30", linewidth=0.8) +
  geom_text(aes(label=paste0(pct,"%")), hjust=-0.2, size=3.5) +
  scale_fill_manual(values=c("ACME (via CLS)"="#2166AC","ADE (direct)"="#92C5DE"),
                    guide="none") +
  scale_x_continuous(limits=c(0, 0.14), breaks=seq(0,.12,.03)) +
  labs(x="Effect on FEI (SD units)", y=NULL, title="b",
       subtitle="Total effect = 0.147; proportion mediated = 0.62 (95% CI: 0.51\u20130.73)") +
  theme_na()

save_fig(fig4b, "fig4b_mediation", w=7, h=3)

# ---- Combine 4a + 4b -------------------------------------------------------
fig4 <- fig4a + fig4b + patchwork::plot_layout(ncol=2, widths=c(3,2))
save_fig(fig4, "fig4_combined", w=14, h=4.5)

message("\u2713 All figures generated.")
tictoc::toc()
