use super::{Args, Matrix, path_with_suffix, rasterize_svg, xml_escape};
use crate::candidate::{self, Fit, Sequence};
use anyhow::{Result, bail};
use std::fmt::Write;
use std::fs;
pub(super) struct Analysis {
    fits: Vec<Fit>,
    region_scores: Vec<Vec<RegionScore>>,
    sequences: Vec<Sequence>,
}
struct RegionScore {
    total_excess: f64,
    rank: usize,
}

// Keep regions in genomic order; break score ties by earlier position.
fn score_regions(fit: &Fit, counts: &[usize], sequences: &[Sequence]) -> Vec<RegionScore> {
    let Some(best) = fit.best else {
        return Vec::new();
    };
    let background = fit.starts[best].fitted.means[0];
    let mut scores: Vec<_> = candidate::regions_in_sequences(&fit.states, sequences)
        .iter()
        .map(|&(a, b)| RegionScore {
            total_excess: counts[a..b].iter().map(|&y| y as f64).sum::<f64>()
                - (b - a) as f64 * background,
            rank: 0,
        })
        .collect();
    let mut order: Vec<_> = (0..scores.len()).collect();
    order.sort_by(|&a, &b| {
        scores[b]
            .total_excess
            .total_cmp(&scores[a].total_excess)
            .then(a.cmp(&b))
    });
    for (rank, index) in order.into_iter().enumerate() {
        scores[index].rank = rank + 1;
    }
    scores
}
const LABELS: [&str; 3] = ["short", "primary", "long"];
fn save(prefix: &std::path::Path, suffix: &str, data: &str) -> Result<()> {
    fs::write(path_with_suffix(prefix, suffix), data)?;
    Ok(())
}
#[cfg(test)]
pub(super) fn analyse(args: &Args, matrix: &Matrix, counts: &[usize]) -> Result<Analysis> {
    analyse_with_sequences(args, matrix, counts, &Sequence::single(counts.len()))
}
pub(super) fn analyse_with_sequences(
    args: &Args,
    matrix: &Matrix,
    counts: &[usize],
    sequences: &[Sequence],
) -> Result<Analysis> {
    let mut next = 0;
    let mut sequence_ids = std::collections::HashSet::new();
    for seq in sequences {
        if seq.start != next
            || seq.end <= seq.start
            || seq.end > counts.len()
            || seq.id.is_empty()
            || !sequence_ids.insert(seq.id.as_str())
        {
            bail!(
                "candidate sequences must uniquely identify consecutive nonempty ranges covering every feature"
            );
        }
        next = seq.end;
    }
    if next != counts.len()
        || matrix.genes.len() != counts.len()
        || matrix.frequencies.len() != counts.len()
    {
        bail!("candidate sequence/count/feature dimensions do not match");
    }
    let singleton_segments = sequences.iter().filter(|s| s.end - s.start == 1).count();
    let transition_count: usize = sequences.iter().map(|s| s.end - s.start - 1).sum();
    let inference_mode = if args.segments.is_some() {
        "segmented"
    } else {
        "single_sequence"
    };
    if !sequences.is_empty() && transition_count == 0 {
        eprintln!(
            "Candidate regions: all segments are singletons; no within-segment transition information is available"
        );
    }
    let mut owners = vec![0; counts.len()];
    for (i, seq) in sequences.iter().enumerate() {
        owners[seq.start..seq.end].fill(i);
    }

    let bg = &args.background_duration_bounds;
    if bg.len() != 2 || !bg.iter().all(|x| x.is_finite()) || bg[0] <= 1.0 || bg[1] <= bg[0] {
        bail!("background duration bounds must be finite and satisfy 1 < MIN < MAX")
    }
    if args.candidate_simulations == 0 {
        bail!("--candidate-simulations must be positive")
    }
    let (defaults, fallback) =
        candidate::duration_defaults_sequences(counts, &matrix.frequencies, sequences);
    let requested = if args.elevated_duration_means.is_empty() {
        defaults
    } else {
        if args.elevated_duration_means.len() != 3
            || args
                .elevated_duration_means
                .iter()
                .any(|x| !x.is_finite() || *x < 1.0)
        {
            bail!("supply three finite elevated duration means >= 1")
        }
        args.elevated_duration_means.clone().try_into().unwrap()
    };
    let durations = if args.elevated_duration_means.is_empty() {
        requested.map(|x| x.max(2.0))
    } else {
        requested
    };
    let mut fits = Vec::new();
    for (index, duration) in durations.iter().enumerate() {
        eprintln!(
            "Fitting candidate regions: {} (elevated mean {duration})",
            LABELS[index]
        );
        let fit = candidate::fit_sequences(counts, sequences, *duration, (bg[0], bg[1]));
        if fit.best.is_none() {
            eprintln!("Candidate regions {}: {}", LABELS[index], fit.diagnostic)
        }
        fits.push(fit);
    }
    let region_scores: Vec<_> = fits
        .iter()
        .map(|fit| score_regions(fit, counts, sequences))
        .collect();
    let mut genes = String::from(
        "sequence\tposition\tgene\tevent_count\tprevalence\tmarginal_elevated_posterior\tviterbi_state\tregion_id\tduration_setting\televated_duration\tsegment_position\n",
    );
    let mut regions = String::from(
        "sequence\tstart\tend\tordered_gene_ids\tgene_count\tevent_count_sum\tevent_count_mean\tminimum_marginal_posterior\tmean_marginal_posterior\tsingleton\tregion_id\tduration_setting\ttotal_excess\ttotal_excess_rank\tsegment_start\tsegment_end\n",
    );
    let mut starts = String::from(
        "duration_setting\trequested_duration\teffective_duration\tstart\tinitial_mu0\tinitial_mu1\tinitial_alpha\tinitial_background_duration\tmu0\tmu1\talpha\tbackground_duration\tlog_likelihood\tconverged\tvalid\tselected\tbounds_reached\titerations\tdiagnostic\n",
    );
    for (k, fit) in fits.iter().enumerate() {
        if fit.starts.is_empty() {
            writeln!(
                starts,
                "{}\t{}\t{}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tfalse\tfalse\tfalse\tNA\t0\t{}",
                LABELS[k], requested[k], fit.duration, fit.diagnostic
            )?
        }
        for (j, s) in fit.starts.iter().enumerate() {
            writeln!(
                starts,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                LABELS[k],
                requested[k],
                fit.duration,
                j + 1,
                s.initial.means[0],
                s.initial.means[1],
                s.initial.alpha,
                s.initial.background,
                s.fitted.means[0],
                s.fitted.means[1],
                s.fitted.alpha,
                s.fitted.background,
                s.likelihood,
                s.converged,
                s.valid,
                fit.best == Some(j),
                s.bounds,
                s.iterations,
                fit.diagnostic
            )?
        }
        let spans = candidate::regions_in_sequences(&fit.states, sequences);
        let mut ids = vec![String::from("NA"); counts.len()];
        for (j, &(a, b)) in spans.iter().enumerate() {
            let id = format!("{}_{}", LABELS[k], j + 1);
            ids[a..b].fill(id.clone());
            let seq = &sequences[owners[a]];
            let sum: usize = counts[a..b].iter().sum();
            let min = fit.posterior[a..b].iter().copied().fold(1.0, f64::min);
            let mean = fit.posterior[a..b].iter().sum::<f64>() / (b - a) as f64;
            writeln!(
                regions,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                seq.id,
                a + 1,
                b,
                matrix.genes[a..b].join(","),
                b - a,
                sum,
                sum as f64 / (b - a) as f64,
                min,
                mean,
                b - a == 1,
                id,
                LABELS[k],
                region_scores[k][j].total_excess,
                region_scores[k][j].rank,
                a - seq.start + 1,
                b - seq.start
            )?;
        }
        for i in 0..counts.len() {
            let prevalence = matrix.frequencies[i]
                .map(|(a, b)| format!("{a}/{b}"))
                .unwrap_or("NA".into());
            let (posterior, state) = if fit.best.is_some() {
                (fit.posterior[i].to_string(), fit.states[i].to_string())
            } else {
                ("NA".into(), "NA".into())
            };
            writeln!(
                genes,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                sequences[owners[i]].id,
                i + 1,
                matrix.genes[i],
                counts[i],
                prevalence,
                posterior,
                state,
                ids[i],
                LABELS[k],
                fit.duration,
                i - sequences[owners[i]].start + 1
            )?;
        }
    }
    starts = starts.lines().enumerate().map(|(i,line)| {
        if i==0 {
            format!("{line}\tinference_mode\tsequence_count\tsingleton_segment_count\twithin_segment_transitions\n")
        } else {
            format!("{line}\t{inference_mode}\t{}\t{singleton_segments}\t{transition_count}\n",sequences.len())
        }
    }).collect();
    let overlap = overlap_table(&fits, sequences)?;
    save(&args.output, ".candidate_genes.tsv", &genes)?;
    save(&args.output, ".candidate_regions.tsv", &regions)?;
    save(&args.output, ".candidate_fits.tsv", &starts)?;
    save(&args.output, ".candidate_overlap.tsv", &overlap)?;
    let mut diagnostics = format!(
        "# seed={} simulations={} duration_fallback={} duration_override={}\n# inference_mode={} sequence_count={} singleton_segments={} within_segment_transitions={}\n# calibration=mean_excluding_singletons automatic_minimum_duration=2\n# calibrated_requested_means={:?} requested_means={:?} background_bounds={:?}\n# Conditional fitted-model checks, not a biological null. Counts may share tracts.\nduration_setting\tmetric\tobserved\tsimulated_mean\tsimulated_q025\tsimulated_q975\tvalid_simulations\n",
        args.candidate_seed,
        args.candidate_simulations,
        fallback,
        args.elevated_duration_means.len() == 3,
        inference_mode,
        sequences.len(),
        singleton_segments,
        transition_count,
        defaults,
        requested,
        bg
    );
    let mut svg = String::from(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="1050"><rect width="100%" height="100%" fill="white"/>"##,
    );
    for (k, fit) in fits.iter().enumerate() {
        let top = 30.0 + k as f64 * 285.0;
        writeln!(
            svg,
            r#"<text x="20" y="{top}" font-family="sans-serif" font-size="14">{}: duration {}; {}</text>"#,
            LABELS[k],
            fit.duration,
            xml_escape(&fit.diagnostic)
        )?;
        let Some(best) = fit.best else { continue };
        let maximum = counts.iter().copied().max().unwrap_or(0);
        let tail = maximum.min(20).max(1);
        let observed = metrics(counts, tail, sequences);
        let mut simulations = vec![Vec::new(); observed.len()];
        for j in 0..args.candidate_simulations {
            let seed = args
                .candidate_seed
                .wrapping_add((k as u64).wrapping_mul(1_000_003))
                .wrapping_add(j as u64);
            let simulated = candidate::simulate_sequences(sequences, fit.starts[best].fitted, seed);
            for (target, (_, v)) in simulations
                .iter_mut()
                .zip(metrics(&simulated, tail, sequences))
            {
                if v.is_finite() {
                    target.push(v)
                }
            }
        }
        let mut frequency_points = String::new();
        let mut simulated_points = String::new();
        for (i, ((name, value), values)) in observed.iter().zip(&mut simulations).enumerate() {
            values.sort_by(f64::total_cmp);
            let mean = if values.is_empty() {
                f64::NAN
            } else {
                values.iter().sum::<f64>() / values.len() as f64
            };
            let lo = quantile(values, 0.025);
            let hi = quantile(values, 0.975);
            writeln!(
                diagnostics,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                LABELS[k],
                name,
                number(*value),
                number(mean),
                number(lo),
                number(hi),
                values.len()
            )?;
            if i <= tail {
                let x = 40.0 + i as f64 * 480.0 / tail as f64;
                let y = top + 145.0 - value * 100.0;
                let sy = top + 145.0 - mean * 100.0;
                write!(frequency_points, "{x:.1},{y:.1} ")?;
                write!(simulated_points, "{x:.1},{sy:.1} ")?;
            } else {
                let y = top + 30.0 + (i - tail - 1) as f64 * 19.0;
                writeln!(
                    svg,
                    r#"<text x="550" y="{y}" font-family="sans-serif" font-size="11">{}: observed {:.3}; simulated {:.3} [{:.3}, {:.3}]</text>"#,
                    name, value, mean, lo, hi
                )?;
            }
        }
        writeln!(
            svg,
            r##"<polyline points="{frequency_points}" fill="none" stroke="#222222"/><polyline points="{simulated_points}" fill="none" stroke="#2166ac"/><text x="40" y="{}" font-family="sans-serif" font-size="11">Count frequencies: black observed; blue simulated mean; bins 0..{}, final bin ≥{}</text>"##,
            top + 163.0,
            tail,
            tail
        )?;
        let valid: Vec<_> = fit
            .starts
            .iter()
            .filter(|s| s.valid)
            .map(|s| s.likelihood)
            .collect();
        let min = valid.iter().copied().fold(f64::INFINITY, f64::min);
        let max = valid.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        writeln!(
            svg,
            r#"<text x="40" y="{}" font-family="sans-serif" font-size="11">Successful starts: {}/12; likelihood range {:.5} .. {:.5} (see fits table)</text>"#,
            top + 190.0,
            valid.len(),
            min,
            max
        )?;
        for (j, start) in fit.starts.iter().enumerate() {
            if start.likelihood.is_finite() {
                let x = 45.0 + j as f64 * 35.0;
                let y = top + 239.0 - (start.likelihood - min) / (max - min).max(1.0) * 25.0;
                writeln!(
                    svg,
                    r##"<circle cx="{x}" cy="{y}" r="3" fill="{}"><title>start {}: likelihood {}; valid {}</title></circle>"##,
                    if start.valid { "#2166ac" } else { "#999999" },
                    j + 1,
                    start.likelihood,
                    start.valid
                )?;
            }
        }
    }
    for (k, fit) in fits.iter().enumerate() {
        let y = 920.0 + k as f64 * 28.0;
        writeln!(
            svg,
            r#"<text x="20" y="{}" font-family="sans-serif" font-size="11">{}</text>"#,
            y + 10.0,
            LABELS[k]
        )?;
        for (a, b) in candidate::regions_in_sequences(&fit.states, sequences) {
            let x = 110.0 + a as f64 / counts.len() as f64 * 940.0;
            let w = (b - a) as f64 / counts.len() as f64 * 940.0;
            writeln!(
                svg,
                r##"<rect x="{x}" y="{y}" width="{w}" height="15" fill="#238b45"/>"##
            )?;
        }
    }
    svg.push_str(r#"<text x="20" y="1025" font-family="sans-serif" font-size="11">Regions along tsv_order; fitted-model simulations are not a biological null.</text></svg>"#);
    save(&args.output, ".candidate_diagnostics.tsv", &diagnostics)?;
    save(&args.output, ".candidate_diagnostics.svg", &svg)?;
    fs::write(
        path_with_suffix(&args.output, ".candidate_diagnostics.png"),
        rasterize_svg(&svg)?,
    )?;
    Ok(Analysis {
        fits,
        region_scores,
        sequences: sequences.to_vec(),
    })
}
fn overlap_table(fits: &[Fit], sequences: &[Sequence]) -> Result<String> {
    let mut overlap = String::from(
        "primary_region\tduration_setting\tmatching_region\tintersection_genes\tjaccard\tprimary_coverage\tmatching_coverage\tstatus\tsequence\n",
    );
    if fits[1].best.is_none() {
        overlap.push_str("NA\tNA\tNA\tNA\tNA\tNA\tNA\tprimary_unavailable\tNA\n")
    }
    for (j, (a, b)) in candidate::regions_in_sequences(&fits[1].states, sequences)
        .into_iter()
        .enumerate()
    {
        let seq = sequences
            .iter()
            .find(|seq| a >= seq.start && b <= seq.end)
            .unwrap();
        for k in [0, 2] {
            let mut found = false;
            for (m, (c, d)) in candidate::regions_in_sequences(&fits[k].states, sequences)
                .into_iter()
                .enumerate()
            {
                let intersection = b.min(d).saturating_sub(a.max(c));
                if intersection > 0 && c >= seq.start && d <= seq.end {
                    found = true;
                    writeln!(
                        overlap,
                        "primary_{}\t{}\t{}_{}\t{}\t{}\t{}\t{}\tmatch\t{}",
                        j + 1,
                        LABELS[k],
                        LABELS[k],
                        m + 1,
                        intersection,
                        intersection as f64 / (b - a + d - c - intersection) as f64,
                        intersection as f64 / (b - a) as f64,
                        intersection as f64 / (d - c) as f64,
                        seq.id
                    )?;
                }
            }
            if !found {
                writeln!(
                    overlap,
                    "primary_{}\t{}\tNA\t0\t0\t0\tNA\t{}\t{}",
                    j + 1,
                    LABELS[k],
                    if fits[k].best.is_none() {
                        "fit_unavailable"
                    } else {
                        "no_overlap"
                    },
                    seq.id
                )?
            }
        }
    }
    Ok(overlap)
}

fn number(x: f64) -> String {
    if x.is_finite() {
        format!("{x:.6}")
    } else {
        "NA".into()
    }
}
fn quantile(x: &[f64], q: f64) -> f64 {
    if x.is_empty() {
        f64::NAN
    } else {
        x[((x.len() - 1) as f64 * q) as usize]
    }
}
fn metrics(y: &[usize], tail: usize, sequences: &[Sequence]) -> Vec<(String, f64)> {
    let n = y.len() as f64;
    let mut m: Vec<_> = (0..=tail)
        .map(|i| {
            (
                format!(
                    "frequency_{}{}",
                    i,
                    if i == tail { "_and_above" } else { "" }
                ),
                y.iter()
                    .filter(|&&v| if i == tail { v >= i } else { v == i })
                    .count() as f64
                    / n,
            )
        })
        .collect();
    let pairs: Vec<_> = sequences
        .iter()
        .flat_map(|seq| {
            y[seq.start..seq.end]
                .windows(2)
                .map(|pair| (pair[0] as f64, pair[1] as f64))
        })
        .collect();
    let correlation = if pairs.is_empty() {
        f64::NAN
    } else {
        let am = pairs.iter().map(|p| p.0).sum::<f64>() / pairs.len() as f64;
        let bm = pairs.iter().map(|p| p.1).sum::<f64>() / pairs.len() as f64;
        let covariance = pairs.iter().map(|(a, b)| (a - am) * (b - bm)).sum::<f64>();
        let variance = pairs.iter().map(|(a, _)| (a - am).powi(2)).sum::<f64>()
            * pairs.iter().map(|(_, b)| (b - bm).powi(2)).sum::<f64>();
        if variance == 0.0 {
            f64::NAN
        } else {
            covariance / variance.sqrt()
        }
    };
    let states: Vec<_> = y.iter().map(|&v| usize::from(v > 0)).collect();
    let r = candidate::regions_in_sequences(&states, sequences);
    let positive = r.iter().map(|(a, b)| b - a).sum::<usize>();
    m.extend([
        ("adjacent_correlation".into(), correlation),
        ("positive_clusters".into(), r.len() as f64),
        (
            "mean_positive_run".into(),
            if r.is_empty() {
                0.0
            } else {
                positive as f64 / r.len() as f64
            },
        ),
        (
            "max_positive_run".into(),
            r.iter().map(|(a, b)| b - a).max().unwrap_or(0) as f64,
        ),
    ]);
    m
}
pub(super) fn draw_region_shading(
    analysis: &Analysis,
    x: f64,
    top: f64,
    column: f64,
    width: f64,
    svg: &mut String,
) {
    let fit = &analysis.fits[1];
    let _ = writeln!(
        svg,
        r#"<text x="{x}" y="{}" font-family="sans-serif" font-size="10">Primary duration {}; sequences: {}; {}</text>"#,
        top + 110.0,
        fit.duration,
        analysis.sequences.len(),
        xml_escape(&fit.diagnostic)
    );
    for edge in [x, x + width] {
        let _ = writeln!(
            svg,
            r##"<line x1="{edge}" x2="{edge}" y1="{top}" y2="{}" stroke="#777777" stroke-dasharray="2 3"/>"##,
            top + 78.0
        );
    }
    if fit.best.is_none() {
        return;
    }
    for (j, (a, b)) in candidate::regions_in_sequences(&fit.states, &analysis.sequences)
        .into_iter()
        .enumerate()
    {
        let color = if analysis.region_scores[1][j].rank <= 5 {
            "#ffd700"
        } else {
            "#238b45"
        };
        let left = x + a as f64 * column;
        let region_width = (b - a) as f64 * column;
        let _ = writeln!(
            svg,
            r##"<rect x="{left}" y="{top}" width="{region_width}" height="78" fill="{color}" fill-opacity="0.20"/>"##,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn excess_ranking_breaks_ties_by_position_and_shades_exactly_five() {
        let p = crate::candidate::Parameters {
            means: [2.0, 10.0],
            alpha: 1.0,
            background: 20.0,
            elevated: 4.0,
        };
        let fit = Fit {
            duration: 4.0,
            starts: vec![crate::candidate::Start {
                initial: p,
                fitted: p,
                likelihood: -1.0,
                converged: true,
                valid: true,
                bounds: String::new(),
                iterations: 1,
            }],
            best: Some(0),
            posterior: vec![0.9; 12],
            states: vec![1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
            diagnostic: "successful".into(),
        };
        let scores = score_regions(
            &fit,
            &[8, 0, 3, 0, 8, 0, 2, 0, 10, 0, 0, 0],
            &Sequence::single(12),
        );
        assert_eq!(
            scores.iter().map(|s| s.total_excess).collect::<Vec<_>>(),
            vec![6.0, 1.0, 6.0, 0.0, 8.0, -2.0]
        );
        assert_eq!(
            scores.iter().map(|s| s.rank).collect::<Vec<_>>(),
            vec![2, 4, 3, 5, 1, 6]
        );
        let analysis = Analysis {
            fits: vec![crate::candidate::fit(&[], 2.0, (1.001, 1e6)), fit],
            region_scores: vec![vec![], scores],
            sequences: Sequence::single(12),
        };
        let mut svg = String::new();
        draw_region_shading(&analysis, 100.0, 20.0, 2.0, 24.0, &mut svg);
        assert_eq!(svg.matches("fill=\"#ffd700\"").count(), 5);
        assert_eq!(svg.matches("fill=\"#238b45\"").count(), 1);
    }
    #[test]
    fn overlap_table_retains_split_and_merged_matches() {
        let make = |states| Fit {
            duration: 4.0,
            starts: vec![],
            best: Some(0),
            posterior: vec![],
            states,
            diagnostic: String::new(),
        };
        let fits = vec![
            make(vec![0, 1, 0, 1, 0, 0]),
            make(vec![0, 1, 1, 1, 0, 1]),
            make(vec![1, 1, 1, 1, 1, 1]),
        ];
        let table = overlap_table(&fits, &Sequence::single(6)).unwrap();
        assert!(
            table.contains(
                "primary_1\tshort\tshort_1\t1\t0.3333333333333333\t0.3333333333333333\t1"
            )
        );
        assert!(table.contains("primary_1\tshort\tshort_2"));
        assert!(table.contains("primary_1\tlong\tlong_1"));
        assert!(table.contains("primary_2\tlong\tlong_1"));
        assert!(table.contains("primary_2\tshort\tNA"));
    }
    #[test]
    fn failed_fits_emit_na_and_no_green_calls() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("failure");
        let args = Args {
            tsv: "unused".into(),
            segments: None,
            tree: Some("unused".into()),
            output: prefix.clone(),
            candidate_regions: true,
            elevated_duration_means: vec![],
            background_duration_bounds: vec![1.001, 1e6],
            candidate_seed: 42,
            candidate_simulations: 2,
        };
        let matrix = Matrix {
            genes: vec!["a".into(), "b".into()],
            isolates: vec![],
            values: vec![],
            frequencies: vec![None, None],
        };
        let analysis = analyse(&args, &matrix, &[0, 0]).unwrap();
        assert!(analysis.fits.iter().all(|fit| fit.best.is_none()));
        let rows = fs::read_to_string(path_with_suffix(&prefix, ".candidate_genes.tsv")).unwrap();
        assert!(rows.contains("\tNA\tNA\tNA\tprimary\t4"));
        let mut svg = String::new();
        draw_region_shading(&analysis, 100.0, 20.0, 2.0, 4.0, &mut svg);
        assert!(!svg.contains("#238b45"));
        let regions =
            fs::read_to_string(path_with_suffix(&prefix, ".candidate_regions.tsv")).unwrap();
        assert_eq!(regions.lines().count(), 1);
    }
}

#[cfg(test)]
mod segmented_tests {
    use super::*;
    #[test]
    fn diagnostic_runs_and_correlations_exclude_boundaries() {
        let seqs = vec![
            Sequence {
                id: "a".into(),
                start: 0,
                end: 2,
            },
            Sequence {
                id: "b".into(),
                start: 2,
                end: 4,
            },
        ];
        let m = metrics(&[1, 2, 10, 20], 20, &seqs);
        let value = |name: &str| m.iter().find(|(n, _)| n == name).unwrap().1;
        assert!((value("adjacent_correlation") - 1.0).abs() < 1e-10);
        assert_eq!(value("positive_clusters"), 2.0);
        assert_eq!(value("mean_positive_run"), 2.0);
        assert_eq!(value("max_positive_run"), 2.0);
        let singles: Vec<_> = (0..4)
            .map(|i| Sequence {
                id: i.to_string(),
                start: i,
                end: i + 1,
            })
            .collect();
        assert!(
            metrics(&[1, 2, 10, 20], 20, &singles)
                .iter()
                .find(|(n, _)| n == "adjacent_correlation")
                .unwrap()
                .1
                .is_nan()
        );
    }
    #[test]
    fn adjacent_elevated_segments_have_separate_regions_and_rankings() {
        let p = candidate::Parameters {
            means: [1.0, 10.0],
            alpha: 1.0,
            background: 20.0,
            elevated: 4.0,
        };
        let make = || Fit {
            duration: 4.0,
            starts: vec![candidate::Start {
                initial: p,
                fitted: p,
                likelihood: -1.0,
                converged: true,
                valid: true,
                bounds: String::new(),
                iterations: 1,
            }],
            best: Some(0),
            posterior: vec![0.99; 3],
            states: vec![1; 3],
            diagnostic: "successful".into(),
        };
        let sequences = vec![
            Sequence {
                id: "a".into(),
                start: 0,
                end: 2,
            },
            Sequence {
                id: "b".into(),
                start: 2,
                end: 3,
            },
        ];
        let scores = score_regions(&make(), &[5, 5, 20], &sequences);
        assert_eq!(scores.len(), 2);
        assert_eq!(scores[0].total_excess, 8.0);
        assert_eq!(scores[1].total_excess, 19.0);
        assert_eq!(scores[1].rank, 1);
        let fits = vec![make(), make(), make()];
        let overlap = overlap_table(&fits, &sequences).unwrap();
        assert!(
            overlap
                .lines()
                .skip(1)
                .all(|line| line.ends_with("\ta") || line.ends_with("\tb"))
        );
        assert!(overlap.contains("primary_2\tshort\tshort_2"));
        let analysis = Analysis {
            fits,
            region_scores: vec![
                score_regions(&make(), &[5, 5, 20], &sequences),
                score_regions(&make(), &[5, 5, 20], &sequences),
                score_regions(&make(), &[5, 5, 20], &sequences),
            ],
            sequences,
        };
        let mut svg = String::new();
        draw_region_shading(&analysis, 0.0, 0.0, 2.0, 6.0, &mut svg);
        assert_eq!(svg.matches("fill=\"#ffd700\"").count(), 2);
    }
}
