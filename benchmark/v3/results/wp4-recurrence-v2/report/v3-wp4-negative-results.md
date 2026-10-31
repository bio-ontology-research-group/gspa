# WP4 negative results (§11), selected from the scored tables

- Narrowing in the PATHWAY relation is not testable: fewer labelled multi-candidate assertions with the labelled function available in both genome sets than the pre-registered power gate requires. This is neither a positive nor a negative result.
- Label-free narrowing in the PATHWAY relation is no more frequent than under shuffles.
- Narrowing in the UNIVERSAL relation is not testable: fewer labelled multi-candidate assertions with the labelled function available in both genome sets than the pre-registered power gate requires. This is neither a positive nor a negative result.
- Narrowing in the ORGANISM relation is not testable: fewer labelled multi-candidate assertions with the labelled function available in both genome sets than the pre-registered power gate requires. This is neither a positive nor a negative result.
- Label-free narrowing in the ORGANISM relation is no more frequent than under shuffles.
- Close-relative recurrence does not improve the current gbdt models in the ORGANISM relation (R2 vs R0).
- Close-relative recurrence does not improve the current logistic models in the ORGANISM relation (R2 vs R0).
- Close-relative recurrence does not improve the current gbdt models in the PATHWAY relation (R2 vs R0).
- Close-relative recurrence does not improve the current logistic models in the PATHWAY relation (R2 vs R0).
- Close-relative recurrence does not improve the current gbdt models in the UNIVERSAL relation (R2 vs R0).
- Close-relative recurrence does not improve the current logistic models in the UNIVERSAL relation (R2 vs R0).
- Relations do not differ on paired assertions: no relation advantage.
- Recurrence does not localize genuinely missing functions (R0, ORGANISM).
- Recurrence does not localize genuinely missing functions (R0, PATHWAY).
- Recurrence does not localize genuinely missing functions (R0, UNIVERSAL).
- Recurrence does not localize genuinely missing functions (R2, ORGANISM).
- Recurrence does not localize genuinely missing functions (R2, PATHWAY).
- Recurrence does not localize genuinely missing functions (R2, UNIVERSAL).

## All §11 templates (ready to use)

- `T2_no_narrowing`: No narrowing toward labelled function from close relatives in the {ctx} relation.
- `T2_not_testable`: Narrowing in the {ctx} relation is not testable: fewer labelled multi-candidate assertions with the labelled function available in both genome sets than the pre-registered power gate requires. This is neither a positive nor a negative result.
- `T2T3_background`: In the {ctx} relation the signal cannot be told apart from density-matched background.
- `T3_none`: No evidence that complementary pieces combine in the {ctx} relation.
- `T4_shuffle`: Label-free narrowing in the {ctx} relation is no more frequent than under shuffles.
- `T6_no_gain`: Close-relative recurrence does not improve the current {kind} models in the {ctx} relation ({contrast}).
- `T6_abstain`: The current {kind} models cannot use close-relative recurrence in the {ctx} relation: they abstain on R2 features.
- `relations_equal`: Relations do not differ on paired assertions: no relation advantage.
- `T7_defect`: Engine or holdout defect: a target feature or probability changed beyond tolerance in an invariance control. Stop and fix before interpreting anything.
- `T8_no_localization`: Recurrence does not localize genuinely missing functions ({regime}, {ctx}).
- `T8_harm`: Learned reassignment of the removed function is more frequent than under SEQONLY ({regime}, {ctx}, {kind}): harm.
- `T9_draft`: The benefit does not hold on draft assemblies ({ctx}, {tier} tier).
