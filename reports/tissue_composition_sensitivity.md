# Tissue-composition sensitivity

MCP-counter 1.2.0 was run on the independently RMA-normalized GSE44076
gene-symbol matrix. The method's official HUGO marker table was pinned to
commit `b6eac73e` and SHA-256 `408f6c5d02c8f9bd2f1599c598367881853520934f8ce398a4886a8b296922bb`. MCP-counter
was developed and validated using transcriptomic data that included Affymetrix
microarrays. Scores are comparable across samples within this cohort; they are
not literal cell fractions.

              population                                           method
                 T cells MCP-counter 1.2.0; age/sex-adjusted linear model
             CD8 T cells MCP-counter 1.2.0; age/sex-adjusted linear model
   Cytotoxic lymphocytes MCP-counter 1.2.0; age/sex-adjusted linear model
               B lineage MCP-counter 1.2.0; age/sex-adjusted linear model
                NK cells MCP-counter 1.2.0; age/sex-adjusted linear model
       Monocytic lineage MCP-counter 1.2.0; age/sex-adjusted linear model
 Myeloid dendritic cells MCP-counter 1.2.0; age/sex-adjusted linear model
             Neutrophils MCP-counter 1.2.0; age/sex-adjusted linear model
       Endothelial cells MCP-counter 1.2.0; age/sex-adjusted linear model
             Fibroblasts MCP-counter 1.2.0; age/sex-adjusted linear model
 healthy_n adjacent_normal_n healthy_mean adjacent_normal_mean
        50                98     3.302987             3.281778
        50                98     2.724349             2.595164
        50                98     3.118937             2.987986
        50                98     3.335368             3.389888
        50                98     2.989355             2.834764
        50                98     4.019373             3.879736
        50                98     3.064017             2.984917
        50                98     3.659134             3.676672
        50                98     3.628215             3.996405
        50                98     5.275091             6.688717
 adjusted_difference standard_error  statistic      p_value adjusted_p_value
        -0.005446566     0.03851102 -0.1414287 8.877289e-01     9.099707e-01
        -0.133755777     0.06027442 -2.2191134 2.804341e-02     4.673902e-02
        -0.142625128     0.04753986 -3.0001167 3.181235e-03     6.362469e-03
         0.077032125     0.11141034  0.6914271 4.904104e-01     6.130130e-01
        -0.149071116     0.03401660 -4.3823054 2.249194e-05     7.497312e-05
        -0.149260160     0.04663324 -3.2007243 1.686797e-03     4.216992e-03
        -0.059959581     0.04391358 -1.3653995 1.742567e-01     2.489381e-01
        -0.002861882     0.02526501 -0.1132745 9.099707e-01     9.099707e-01
         0.359777037     0.02749316 13.0860576 2.828078e-26     2.828078e-25
         1.376027321     0.11727115 11.7337238 9.863369e-23     4.931685e-22

PC1 and PC2 explained 33.2% and 21.9% of standardized score variance.
The composition-adjusted sensitivity model included tissue, age, sex, and only these two prespecified composition PCs.
Of 101 high-confidence field genes, 74 (73.3%) retained FDR < 0.05,
|log2FC| >= 0.5, and direction after composition adjustment.
Among 10 stromal/ECM-flagged high-confidence genes, 10 attenuated by at least 25%.

Attenuation does not show that the original association was false: composition can be part of the field microenvironment while confounding a strictly epithelial interpretation. Because score genes are derived from the same bulk transcriptome, this is an overadjustment-prone sensitivity analysis rather than a causal decomposition.
