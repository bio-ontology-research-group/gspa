package gspa.model

/** MIMAG quality tiers for MAGs */
enum QualityTier {
    HIGH,    // >= 90% complete, < 5% contamination
    MEDIUM,  // >= 50% complete, < 10% contamination
    LOW,     // Below medium thresholds
    UNKNOWN  // No quality estimates available
}
