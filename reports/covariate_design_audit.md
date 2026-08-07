# Covariate design audit

Age, sex, and tumor location were taken only from deposited GSE44076 metadata.
Stage was not used in healthy-versus-adjacent or tumor-versus-healthy models because
it is structurally undefined for cancer-free donors. Age was modeled continuously;
sex and location were modeled as factors. U2 was executed only after confirming
a full-rank design and a minimum of ten observations in every factor category.

                              model            contrast
                      U0_unadjusted adjacent_vs_healthy
                U1_age_sex_adjusted adjacent_vs_healthy
       U2_age_sex_location_adjusted adjacent_vs_healthy
                      U0_unadjusted    tumor_vs_healthy
                U1_age_sex_adjusted    tumor_vs_healthy
       U2_age_sex_location_adjusted    tumor_vs_healthy
             P_patient_fixed_effect   tumor_vs_adjacent
 U1_age_sex_adjusted_qc_sensitivity adjacent_vs_healthy
                        formula samples_available samples_complete_case
                        ~tissue               148                   148
            ~tissue + age + sex               148                   148
 ~tissue + age + sex + location               148                   148
                        ~tissue               148                   148
            ~tissue + age + sex               148                   148
 ~tissue + age + sex + location               148                   148
    ~ patient_id + tissue_class               196                   196
            ~tissue + age + sex               147                   147
 healthy_or_reference_n comparison_n age_complete sex_complete
                     50           98          148          148
                     50           98          148          148
                     50           98          148          148
                     50           98          148          148
                     50           98          148          148
                     50           98          148          148
                     98           98          196          196
                     50           97          147          147
 location_complete age_mean_reference age_mean_comparison
               148           62.50000            70.54082
               148           62.50000            70.54082
               148           62.50000            70.54082
               148           62.50000            70.54082
               148           62.50000            70.54082
               148           62.50000            70.54082
               196           70.54082            70.54082
               147           62.50000            70.55670
        sex_distribution_reference       sex_distribution_comparison
                 Female:23;Male:27                 Female:27;Male:71
                 Female:23;Male:27                 Female:27;Male:71
                 Female:23;Male:27                 Female:27;Male:71
                 Female:23;Male:27                 Female:27;Male:71
                 Female:23;Male:27                 Female:27;Male:71
                 Female:23;Male:27                 Female:27;Male:71
 absorbed_by_patient_fixed_effects absorbed_by_patient_fixed_effects
                 Female:23;Male:27                 Female:27;Male:70
   location_distribution_reference  location_distribution_comparison
                  Left:23;Right:27                  Left:60;Right:38
                  Left:23;Right:27                  Left:60;Right:38
                  Left:23;Right:27                  Left:60;Right:38
                  Left:23;Right:27                  Left:60;Right:38
                  Left:23;Right:27                  Left:60;Right:38
                  Left:23;Right:27                  Left:60;Right:38
 absorbed_by_patient_fixed_effects absorbed_by_patient_fixed_effects
                  Left:23;Right:27                  Left:59;Right:38
 design_rows design_columns design_rank full_rank minimum_factor_category_n
         148              2           2      TRUE                        50
         148              4           4      TRUE                        50
         148              5           5      TRUE                        50
         148              2           2      TRUE                        50
         148              4           4      TRUE                        50
         148              5           5      TRUE                        50
         196             99          99      TRUE                        98
         147              4           4      TRUE                        50
 condition_number           coefficient excluded_samples
         3.356962 tissueadjacent_normal
       401.364523 tissueadjacent_normal
       370.585484 tissueadjacent_normal
         3.356962           tissuetumor
       401.364523           tissuetumor
       370.585484           tissuetumor
       210.512987     tissue_classtumor
       400.009537 tissueadjacent_normal       GSM1077736

The paired tumor-versus-adjacent model retains patient fixed effects, which absorb
patient-level age, sex, and location.
