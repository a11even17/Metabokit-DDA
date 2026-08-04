use metabokit_app_lib::postprocess::{postprocess_preview, PostOptions};

fn main() {
    let mut args = std::env::args().skip(1);
    let output_dir = args
        .next()
        .expect("usage: postprocess_smoke <results folder> <report>");
    let report = args
        .next()
        .unwrap_or_else(|| "Report_RTseparated.csv".to_string());
    let result = postprocess_preview(
        output_dir,
        report,
        PostOptions {
            minimum_detected_percent: 50.0,
            minimum_peak_shape: 0.8,
            minimum_score: 0.5,
            minimum_sn: 5.0,
            minimum_matching_peaks: 1.0,
            maximum_cv_percent: Some(25.0),
            identified_only: true,
            msms_only: true,
            remove_isf: true,
        },
    )
    .expect("post-processing preview");
    let value = serde_json::to_value(result).expect("serialize preview");
    println!(
        "post-process: {} total, {} kept, {} removal categories, {} preview rows",
        value["totalRows"],
        value["keptRows"],
        value["removedBy"].as_array().map_or(0, Vec::len),
        value["rows"].as_array().map_or(0, Vec::len),
    );
}
