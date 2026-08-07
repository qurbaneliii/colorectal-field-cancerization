# Tier-specific enrichment report

Enrichment was rerun separately for high-confidence up/down field genes, provisional genes, composition-robust and composition-sensitive subsets, stress-gene-removed genes, locked Task B/Task C panels, and full ranked contrasts. GO Biological Process over-representation results were pruned with semantic similarity (`simplify`, cutoff 0.7). Model-panel analyses are explicitly exploratory because the sets contain only a few genes.

All terms are associations and may reflect tissue composition; they are not evidence of causal or epithelial-intrinsic mechanisms.

## Analysis-set audit

```
                   analysis_set input_symbols mapped_entrez minimum_required
       high_confidence_field_up            68            68               10
     high_confidence_field_down            33            33               10
           provisional_field_up           846           846               10
         provisional_field_down           524           524               10
       composition_robust_field            74            74               10
    composition_sensitive_field            27            27               10
 high_confidence_without_stress            98            98               10
            task_b_locked_panel             5             5                2
            task_c_locked_panel             3             3                2
           qc_excluded_field_up           989           989               10
         qc_excluded_field_down           595           595               10
 executed
     TRUE
     TRUE
     TRUE
     TRUE
     TRUE
     TRUE
     TRUE
     TRUE
     TRUE
     TRUE
     TRUE
```

## Top FDR-significant terms

```
                   analysis_set database
         adjacent_vs_healthy_u1       GO
         adjacent_vs_healthy_u1       GO
         adjacent_vs_healthy_u1       GO
         adjacent_vs_healthy_u1       GO
         adjacent_vs_healthy_u1       GO
       composition_robust_field Reactome
       composition_robust_field Reactome
       composition_robust_field Reactome
       composition_robust_field Reactome
       composition_robust_field Reactome
    composition_sensitive_field Reactome
    composition_sensitive_field Reactome
    composition_sensitive_field Reactome
    composition_sensitive_field Reactome
    composition_sensitive_field Reactome
       high_confidence_field_up Reactome
       high_confidence_field_up Reactome
       high_confidence_field_up Reactome
       high_confidence_field_up Reactome
       high_confidence_field_up Reactome
 high_confidence_without_stress Reactome
 high_confidence_without_stress Reactome
 high_confidence_without_stress Reactome
 high_confidence_without_stress Reactome
 high_confidence_without_stress Reactome
         provisional_field_down       GO
         provisional_field_down       GO
         provisional_field_down       GO
         provisional_field_down       GO
         provisional_field_down Reactome
           provisional_field_up Reactome
           provisional_field_up       GO
           provisional_field_up       GO
           provisional_field_up       GO
           provisional_field_up       GO
         qc_excluded_field_down       GO
         qc_excluded_field_down       GO
         qc_excluded_field_down       GO
         qc_excluded_field_down       GO
         qc_excluded_field_down       GO
           qc_excluded_field_up Reactome
           qc_excluded_field_up       GO
           qc_excluded_field_up       GO
           qc_excluded_field_up       GO
           qc_excluded_field_up       GO
            task_c_locked_panel Reactome
            task_c_locked_panel       GO
            task_c_locked_panel       GO
            task_c_locked_panel       GO
            task_c_locked_panel       GO
       tumor_vs_adjacent_paired       GO
       tumor_vs_adjacent_paired       GO
       tumor_vs_adjacent_paired       GO
       tumor_vs_adjacent_paired       GO
       tumor_vs_adjacent_paired       GO
                                                                                                                 Description
                                                                            mitochondrial respiratory chain complex assembly
                                                                                           extracellular matrix organization
                                                                                        extracellular structure organization
                                                                               external encapsulating structure organization
                                                                                                collagen fibril organization
                                                                                                           ECM proteoglycans
                                                                                           Extracellular matrix organization
 Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)
                                                                                     Degradation of the extracellular matrix
                                                                                                        Collagen degradation
                                                                                           Extracellular matrix organization
                                                                                                           ECM proteoglycans
                                                                                     Degradation of the extracellular matrix
                                                                                                   CRMPs in Sema3A signaling
                                                                                          Integrin cell surface interactions
                                                                                           Extracellular matrix organization
                                                                                     Degradation of the extracellular matrix
                                                                                                           ECM proteoglycans
 Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)
                                                                                                        Collagen degradation
                                                                                           Extracellular matrix organization
                                                                                                           ECM proteoglycans
                                                                                     Degradation of the extracellular matrix
 Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)
                                                                                                        Collagen degradation
                                                                                                fatty acid metabolic process
                                                                                           carboxylic acid catabolic process
                                                                                            small molecule catabolic process
                                                                                                  cellular catabolic process
                                                                                                       Fatty acid metabolism
                                                                                           Extracellular matrix organization
                                                                                           extracellular matrix organization
                                                                                                  regulation of angiogenesis
                                                                                                     cell-substrate adhesion
                                                                                                 muscle cell differentiation
                                                                                           carboxylic acid catabolic process
                                                                                                  cellular catabolic process
                                                                                                fatty acid metabolic process
                                                                                            small molecule catabolic process
                                                                                                     lipid catabolic process
                                                                                           Extracellular matrix organization
                                                                                           extracellular matrix organization
                                                                                       regulation of vasculature development
                                                                                                     cell-substrate adhesion
                                                                   regulation of cellular response to growth factor stimulus
                                                                                                       Hyaluronan metabolism
                                                                                                            skin development
                                                                                                       epidermis development
                                                                                             hyaluronan biosynthetic process
                                                                                                hyaluronan catabolic process
                                                                                        ribonucleoprotein complex biogenesis
                                                                                                         ribosome biogenesis
                                                                                                             rRNA processing
                                                                                                      rRNA metabolic process
                                                                                          ribosomal small subunit biogenesis
 adjusted_p_value
     6.307143e-09
     6.307143e-09
     6.307143e-09
     6.307143e-09
     6.307143e-09
     9.365179e-05
     9.365179e-05
     9.365179e-05
     1.544872e-04
     1.035213e-02
     1.762065e-06
     2.936679e-04
     2.294260e-03
     6.944797e-03
     6.944797e-03
     2.454691e-12
     1.920492e-09
     2.463086e-09
     3.961293e-06
     1.859102e-05
     3.245094e-11
     1.664386e-09
     3.242942e-08
     1.366882e-07
     1.078166e-04
     9.306468e-08
     9.306468e-08
     1.048307e-07
     1.048307e-07
     8.239588e-06
     2.150465e-15
     3.428516e-15
     3.428516e-15
     4.460379e-13
     1.556780e-12
     3.706851e-07
     6.773725e-07
     7.011292e-07
     5.804645e-06
     3.590907e-05
     4.723658e-26
     2.354156e-21
     1.875442e-16
     4.676657e-15
     1.456473e-14
     1.597443e-02
     3.472435e-02
     3.472435e-02
     3.472435e-02
     3.472435e-02
     6.944944e-09
     6.944944e-09
     6.944944e-09
     6.944944e-09
     6.944944e-09
```
