//! Geometric-duration negative-binomial HMM, with stationary initial state.
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use statrs::function::gamma::ln_gamma;

/// Independent observed sequence, expressed in global half-open TSV coordinates.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Sequence {
    pub id: String,
    pub start: usize,
    pub end: usize,
}
impl Sequence {
    pub fn single(length: usize) -> Vec<Self> {
        if length == 0 {
            Vec::new()
        } else {
            vec![Self {
                id: "tsv_order".into(),
                start: 0,
                end: length,
            }]
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Parameters {
    pub means: [f64; 2],
    pub alpha: f64,
    pub background: f64,
    pub elevated: f64,
}
impl Parameters {
    fn transitions(self) -> [[f64; 2]; 2] {
        let a = 1.0 / self.background;
        let b = 1.0 / self.elevated;
        [[(1.0 - a).ln(), a.ln()], [b.ln(), (1.0 - b).ln()]]
    }
    fn initial(self) -> [f64; 2] {
        let p = self.elevated / (self.background + self.elevated);
        [(1.0 - p).ln(), p.ln()]
    }
}
fn add(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let m = a.max(b);
    m + ((a - m).exp() + (b - m).exp()).ln()
}
pub(crate) fn emission(y: usize, mu: f64, alpha: f64) -> f64 {
    let r = 1.0 / alpha;
    let y = y as f64;
    ln_gamma(y + r) - ln_gamma(r) - ln_gamma(y + 1.0) - r * (mu / r).ln_1p()
        + y * (mu / (r + mu)).ln()
}
fn forward(y: &[usize], p: Parameters) -> (f64, Vec<[f64; 2]>) {
    if y.is_empty() {
        return (0.0, vec![]);
    }
    let t = p.transitions();
    let initial = p.initial();
    let mut f = vec![[0.0; 2]; y.len()];
    for s in 0..2 {
        f[0][s] = initial[s] + emission(y[0], p.means[s], p.alpha)
    }
    for i in 1..y.len() {
        for s in 0..2 {
            f[i][s] = add(f[i - 1][0] + t[0][s], f[i - 1][1] + t[1][s])
                + emission(y[i], p.means[s], p.alpha)
        }
    }
    (add(f[y.len() - 1][0], f[y.len() - 1][1]), f)
}
pub(crate) fn infer(y: &[usize], p: Parameters) -> (f64, Vec<f64>, Vec<usize>) {
    let (ll, f) = forward(y, p);
    if y.is_empty() {
        return (ll, vec![], vec![]);
    }
    let t = p.transitions();
    let mut b = vec![[0.0; 2]; y.len()];
    for i in (0..y.len() - 1).rev() {
        for s in 0..2 {
            b[i][s] = add(
                t[s][0] + emission(y[i + 1], p.means[0], p.alpha) + b[i + 1][0],
                t[s][1] + emission(y[i + 1], p.means[1], p.alpha) + b[i + 1][1],
            )
        }
    }
    let posterior = (0..y.len())
        .map(|i| (f[i][1] + b[i][1] - ll).exp().clamp(0.0, 1.0))
        .collect();
    let mut v = vec![[0.0; 2]; y.len()];
    let mut previous = vec![[0usize; 2]; y.len()];
    v[0] = f[0];
    for i in 1..y.len() {
        for s in 0..2 {
            let a = v[i - 1][0] + t[0][s];
            let b = v[i - 1][1] + t[1][s];
            let k = usize::from(b > a);
            previous[i][s] = k;
            v[i][s] = a.max(b) + emission(y[i], p.means[s], p.alpha)
        }
    }
    let mut states = vec![0; y.len()];
    states[y.len() - 1] = usize::from(v[y.len() - 1][1] > v[y.len() - 1][0]);
    for i in (1..y.len()).rev() {
        states[i - 1] = previous[i][states[i]]
    }
    (ll, posterior, states)
}
fn likelihood_sequences(y: &[usize], sequences: &[Sequence], p: Parameters) -> f64 {
    sequences
        .iter()
        .map(|seq| forward(&y[seq.start..seq.end], p).0)
        .sum()
}
pub(crate) fn infer_sequences(
    y: &[usize],
    sequences: &[Sequence],
    p: Parameters,
) -> (f64, Vec<f64>, Vec<usize>) {
    let mut likelihood = 0.0;
    let mut posterior = vec![0.0; y.len()];
    let mut states = vec![0; y.len()];
    for seq in sequences {
        let (ll, probability, path) = infer(&y[seq.start..seq.end], p);
        likelihood += ll;
        posterior[seq.start..seq.end].copy_from_slice(&probability);
        states[seq.start..seq.end].copy_from_slice(&path);
    }
    (likelihood, posterior, states)
}

#[derive(Debug)]
pub(crate) struct Start {
    pub initial: Parameters,
    pub fitted: Parameters,
    pub likelihood: f64,
    pub converged: bool,
    pub valid: bool,
    pub bounds: String,
    pub iterations: usize,
}
pub(crate) struct Fit {
    pub duration: f64,
    pub starts: Vec<Start>,
    pub best: Option<usize>,
    pub posterior: Vec<f64>,
    pub states: Vec<usize>,
    pub diagnostic: String,
}
fn decode(x: [f64; 4], duration: f64, limits: &[(f64, f64); 4]) -> Parameters {
    let z: Vec<_> = x
        .iter()
        .zip(limits)
        .map(|(v, (lo, hi))| lo + (hi - lo) / (1.0 + (-v).exp()))
        .collect();
    Parameters {
        means: [z[0], z[0] + (limits[0].1 - z[0]) * z[1]],
        alpha: z[2],
        background: z[3],
        elevated: duration,
    }
}
fn encode(p: Parameters, limits: &[(f64, f64); 4]) -> [f64; 4] {
    let z = [
        p.means[0],
        (p.means[1] - p.means[0]) / (limits[0].1 - p.means[0]),
        p.alpha,
        p.background,
    ];
    std::array::from_fn(|i| {
        let q = ((z[i] - limits[i].0) / (limits[i].1 - limits[i].0)).clamp(1e-7, 1.0 - 1e-7);
        (q / (1.0 - q)).ln()
    })
}
#[cfg(test)]
pub(crate) fn fit(y: &[usize], duration: f64, bg: (f64, f64)) -> Fit {
    fit_sequences(y, &Sequence::single(y.len()), duration, bg)
}
pub(crate) fn fit_sequences(
    y: &[usize],
    sequences: &[Sequence],
    duration: f64,
    bg: (f64, f64),
) -> Fit {
    let mut result = Fit {
        duration,
        starts: vec![],
        best: None,
        posterior: vec![],
        states: vec![],
        diagnostic: String::new(),
    };
    if y.len() < 2 || y.iter().all(|v| *v == y[0]) {
        result.diagnostic = "insufficient count variation; no candidate calls".into();
        return result;
    }
    let mut sorted = y.to_vec();
    sorted.sort_unstable();
    let upper = (100.0 * (sorted[y.len() - 1] as f64)).max(100.0);
    let limits = [(1e-8, upper), (0.0, 1.0), (1e-6, 1e3), bg];
    let mean = y.iter().sum::<usize>() as f64 / y.len() as f64;
    for quantile in [0.5, 0.75, 0.9] {
        for alpha in [0.05, 2.0] {
            for background in [10.0f64, 1000.0] {
                let low = sorted[((y.len() - 1) as f64 * quantile) as usize] as f64;
                let initial = Parameters {
                    means: [(low * 0.5).max(0.01), (low + mean + 1.0).min(upper * 0.9)],
                    alpha,
                    background: background
                        .clamp(bg.0 + (bg.1 - bg.0) * 1e-7, bg.1 - (bg.1 - bg.0) * 1e-7),
                    elevated: duration,
                };
                let x = encode(initial, &limits);
                let objective = |x: [f64; 4]| {
                    let ll = likelihood_sequences(y, sequences, decode(x, duration, &limits));
                    if ll.is_finite() { -ll } else { f64::INFINITY }
                };
                let mut simplex = vec![(x, objective(x))];
                for j in 0..4 {
                    let mut z = x;
                    z[j] += 0.7;
                    simplex.push((z, objective(z)))
                }
                let mut converged = false;
                let mut iterations = 0;
                for iteration in 0..1500 {
                    iterations = iteration + 1;
                    simplex.sort_by(|a, b| a.1.total_cmp(&b.1));
                    let spread = simplex[4].1 - simplex[0].1;
                    let best_point = simplex[0].0;
                    let distance = simplex
                        .iter()
                        .flat_map(|(z, _)| (0..4).map(move |j| (z[j] - best_point[j]).abs()))
                        .fold(0.0, f64::max);
                    if spread < 1e-7 * (1.0 + simplex[0].1.abs()) && distance < 1e-3 {
                        converged = true;
                        break;
                    }
                    let centroid: [f64; 4] =
                        std::array::from_fn(|j| simplex[..4].iter().map(|(z, _)| z[j] / 4.0).sum());
                    let proposal = |scale: f64| {
                        std::array::from_fn(|j| {
                            (centroid[j] + scale * (centroid[j] - simplex[4].0[j]))
                                .clamp(-25.0, 25.0)
                        })
                    };
                    let reflected = proposal(1.0);
                    let fr = objective(reflected);
                    if fr < simplex[0].1 {
                        let expanded = proposal(2.0);
                        let fe = objective(expanded);
                        simplex[4] = if fe < fr {
                            (expanded, fe)
                        } else {
                            (reflected, fr)
                        };
                    } else if fr < simplex[3].1 {
                        simplex[4] = (reflected, fr)
                    } else {
                        let contracted = proposal(if fr < simplex[4].1 { 0.5 } else { -0.5 });
                        let fc = objective(contracted);
                        if fc < fr.min(simplex[4].1) {
                            simplex[4] = (contracted, fc)
                        } else {
                            let best = simplex[0].0;
                            for k in 1..5 {
                                simplex[k].0 =
                                    std::array::from_fn(|j| (simplex[k].0[j] + best[j]) / 2.0);
                                simplex[k].1 = objective(simplex[k].0)
                            }
                        }
                    }
                }
                simplex.sort_by(|a, b| a.1.total_cmp(&b.1));
                let fitted = decode(simplex[0].0, duration, &limits);
                let likelihood = -simplex[0].1;
                let bounds = [
                    fitted.means[0],
                    (fitted.means[1] - fitted.means[0]) / (upper - fitted.means[0]),
                    fitted.alpha,
                    fitted.background,
                ]
                .iter()
                .zip(limits)
                .enumerate()
                .filter_map(|(i, (v, (lo, hi)))| {
                    let q = (v - lo) / (hi - lo);
                    (q < 1e-5 || q > 1.0 - 1e-5).then(|| {
                        [
                            "mean_background",
                            "mean_separation",
                            "dispersion",
                            "background_duration",
                        ][i]
                    })
                })
                .collect::<Vec<_>>()
                .join(",");
                let mut valid =
                    converged && likelihood.is_finite() && fitted.means[1] > 1.01 * fitted.means[0];
                if valid {
                    let (_, posterior, _) = infer_sequences(y, sequences, fitted);
                    let occupancy = posterior.iter().sum::<f64>();
                    valid = occupancy >= 0.5 && posterior.len() as f64 - occupancy >= 0.5;
                }
                result.starts.push(Start {
                    initial,
                    fitted,
                    likelihood,
                    converged,
                    valid,
                    bounds,
                    iterations,
                });
            }
        }
    }
    result.best = result
        .starts
        .iter()
        .enumerate()
        .filter(|(_, s)| s.valid)
        .max_by(|a, b| a.1.likelihood.total_cmp(&b.1.likelihood))
        .map(|(i, _)| i);
    if let Some(i) = result.best {
        let (_, posterior, states) = infer_sequences(y, sequences, result.starts[i].fitted);
        let occupied = posterior.iter().sum::<f64>();
        if occupied < 0.5 || posterior.len() as f64 - occupied < 0.5 {
            result.best = None;
            result.diagnostic =
                "degenerate fit: effectively unoccupied state; no candidate calls".into();
        } else {
            result.posterior = posterior;
            result.states = states;
            result.diagnostic = "successful".into();
        }
    } else {
        result.diagnostic =
            "all starts failed convergence or produced degenerate means; no candidate calls".into()
    }
    result
}
pub(crate) fn regions(states: &[usize]) -> Vec<(usize, usize)> {
    let mut r = vec![];
    let mut i = 0;
    while i < states.len() {
        if states[i] == 0 {
            i += 1;
            continue;
        }
        let start = i;
        while i < states.len() && states[i] == 1 {
            i += 1
        }
        r.push((start, i))
    }
    r
}
pub(crate) fn regions_in_sequences(
    states: &[usize],
    sequences: &[Sequence],
) -> Vec<(usize, usize)> {
    if states.is_empty() {
        return Vec::new();
    }
    sequences
        .iter()
        .flat_map(|seq| {
            regions(&states[seq.start..seq.end])
                .into_iter()
                .map(|(a, b)| (a + seq.start, b + seq.start))
        })
        .collect()
}
#[cfg(test)]
pub(crate) fn duration_defaults(y: &[usize], freq: &[Option<(usize, usize)>]) -> ([f64; 3], bool) {
    duration_defaults_sequences(y, freq, &Sequence::single(y.len()))
}
pub(crate) fn duration_defaults_sequences(
    y: &[usize],
    freq: &[Option<(usize, usize)>],
    sequences: &[Sequence],
) -> ([f64; 3], bool) {
    let mut runs = Vec::new();
    for seq in sequences {
        let mut n = 0;
        for (v, f) in y[seq.start..seq.end].iter().zip(&freq[seq.start..seq.end]) {
            if *v > 0 && f.is_some_and(|(a, b)| a as u128 * 100 >= b as u128 * 95) {
                n += 1;
            } else if n > 0 {
                runs.push(n);
                n = 0;
            }
        }
        if n > 0 {
            runs.push(n);
        }
    }
    runs.retain(|&length| length > 1);
    if runs.is_empty() {
        return ([2.0, 4.0, 8.0], true);
    }
    let mean = runs.iter().map(|&length| length as f64).sum::<f64>() / runs.len() as f64;
    ([mean * 0.5, mean, mean * 2.0], false)
}
pub(crate) fn simulate(n: usize, p: Parameters, seed: u64) -> Vec<usize> {
    use rand::distributions::Distribution;
    use statrs::distribution::{Gamma, Poisson};
    let mut rng = StdRng::seed_from_u64(seed);
    let mut state = usize::from(rng.r#gen::<f64>() < p.initial()[1].exp());
    let mut counts = Vec::with_capacity(n);
    for i in 0..n {
        if i > 0 && rng.r#gen::<f64>() < 1.0 / if state == 0 { p.background } else { p.elevated } {
            state = 1 - state
        }
        let lambda = Gamma::new(1.0 / p.alpha, 1.0 / (p.means[state] * p.alpha))
            .unwrap()
            .sample(&mut rng);
        counts.push(if lambda <= 0.0 {
            0
        } else {
            {
                let value: u64 = Poisson::new(lambda).unwrap().sample(&mut rng);
                value as usize
            }
        });
    }
    counts
}

pub(crate) fn simulate_sequences(sequences: &[Sequence], p: Parameters, seed: u64) -> Vec<usize> {
    sequences
        .iter()
        .enumerate()
        .flat_map(|(i, seq)| {
            simulate(
                seq.end - seq.start,
                p,
                seed.wrapping_add((i as u64).wrapping_mul(1_000_000_007)),
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn posterior_and_viterbi_match_enumeration_including_censored_ends() {
        for elevated in [1.0, 3.0] {
            let p = Parameters {
                means: [0.3, 7.0],
                alpha: 0.4,
                background: 8.0,
                elevated,
            };
            for y in [vec![6], vec![7, 0, 6], vec![0, 8, 8, 0]] {
                let (ll, posterior, states) = infer(&y, p);
                let mut total = 0.0;
                let mut marginals = vec![0.0; y.len()];
                let mut best = (f64::NEG_INFINITY, vec![]);
                for mask in 0..1usize << y.len() {
                    let path: Vec<_> = (0..y.len()).map(|i| (mask >> i) & 1).collect();
                    let mut score = p.initial()[path[0]];
                    for i in 0..y.len() {
                        score += emission(y[i], p.means[path[i]], p.alpha);
                        if i > 0 {
                            score += p.transitions()[path[i - 1]][path[i]]
                        }
                    }
                    let probability = score.exp();
                    total += probability;
                    for i in 0..y.len() {
                        if path[i] == 1 {
                            marginals[i] += probability
                        }
                    }
                    if score > best.0 {
                        best = (score, path)
                    }
                }
                assert!((ll - total.ln()).abs() < 1e-10);
                for i in 0..y.len() {
                    assert!((posterior[i] - marginals[i] / total).abs() < 1e-10)
                }
                assert_eq!(states, best.1);
            }
        }
    }
    #[test]
    fn singleton_regions_and_calibration_breaks() {
        assert_eq!(regions(&[1, 0, 1, 1, 0, 1]), vec![(0, 1), (2, 4), (5, 6)]);
        let y = [1, 1, 1, 0, 1, 1];
        let f = [
            Some((100, 100)),
            Some((100, 100)),
            Some((10, 100)),
            Some((100, 100)),
            None,
            Some((100, 100)),
        ];
        assert_eq!(duration_defaults(&y, &f), ([1.0, 2.0, 4.0], false));
        assert_eq!(
            duration_defaults(&[0, 0], &[None, None]),
            ([2.0, 4.0, 8.0], true)
        );
    }
    #[test]
    fn calibration_uses_mean_of_non_singleton_runs() {
        // Runs 1, 2, 2, 8: dropping the singleton gives mean 4, not median 2.
        let counts = [1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1];
        let frequencies = vec![Some((100, 100)); counts.len()];
        assert_eq!(
            duration_defaults(&counts, &frequencies),
            ([2.0, 4.0, 8.0], false)
        );
        assert_eq!(
            duration_defaults(&[1, 0, 1], &[Some((100, 100)); 3]),
            ([2.0, 4.0, 8.0], true)
        );
    }
    #[test]
    fn simulation_is_seeded_and_fit_recovers_separated_states() {
        let p = Parameters {
            means: [0.2, 18.0],
            alpha: 0.1,
            background: 20.0,
            elevated: 5.0,
        };
        let counts = simulate(150, p, 42);
        assert_eq!(counts, simulate(150, p, 42));
        assert_ne!(counts, simulate(150, p, 43));
        let fitted = fit(&counts, 5.0, (1.001, 1000000.0));
        let best = fitted.best.expect("separated data should fit");
        let s = &fitted.starts[best];
        assert!(s.fitted.means[1] > s.fitted.means[0] * 1.01);
        assert!(
            fitted
                .posterior
                .iter()
                .all(|v| v.is_finite() && *v >= 0.0 && *v <= 1.0)
        );
        assert_eq!(fitted.starts.len(), 12);
        let second = fit(&counts, 5.0, (1.001, 1000000.0));
        assert_eq!(fitted.states, second.states);
        assert_eq!(fitted.posterior, second.posterior);
    }
    #[test]
    fn degenerate_data_has_explicit_failure() {
        for y in [vec![], vec![0], vec![0; 20], vec![3; 20]] {
            let fit = fit(&y, 4.0, (1.001, 1e6));
            assert!(fit.best.is_none());
            assert!(!fit.diagnostic.is_empty())
        }
    }
}

#[cfg(test)]
mod segmented_tests {
    use super::*;
    fn sequences() -> Vec<Sequence> {
        vec![
            Sequence {
                id: "first".into(),
                start: 0,
                end: 2,
            },
            Sequence {
                id: "second".into(),
                start: 2,
                end: 3,
            },
        ]
    }
    #[test]
    fn multiple_sequences_match_exhaustive_enumeration_and_are_conditionally_independent() {
        for elevated in [1.0, 3.0] {
            let p = Parameters {
                means: [0.3, 7.0],
                alpha: 0.4,
                background: 8.0,
                elevated,
            };
            let y = [0, 7, 2];
            let seqs = sequences();
            let (ll, posterior, states) = infer_sequences(&y, &seqs, p);
            let mut total = 0.0;
            let mut marginal = [0.0; 3];
            let mut best = (f64::NEG_INFINITY, vec![]);
            for mask in 0..8 {
                let path: Vec<_> = (0..3).map(|i| (mask >> i) & 1).collect();
                let mut score = 0.0;
                for i in 0..3 {
                    score += emission(y[i], p.means[path[i]], p.alpha);
                    score += if i == 0 || i == 2 {
                        p.initial()[path[i]]
                    } else {
                        p.transitions()[path[i - 1]][path[i]]
                    };
                }
                total += score.exp();
                for i in 0..3 {
                    if path[i] == 1 {
                        marginal[i] += score.exp();
                    }
                }
                if score > best.0 {
                    best = (score, path);
                }
            }
            assert!((ll - total.ln()).abs() < 1e-10);
            assert!((ll - infer(&y[..2], p).0 - infer(&y[2..], p).0).abs() < 1e-10);
            for i in 0..3 {
                assert!((posterior[i] - marginal[i] / total).abs() < 1e-10);
            }
            assert_eq!(states, best.1);
            assert_eq!(
                posterior[..2],
                infer_sequences(&[0, 7, 100], &seqs, p).1[..2]
            );
        }
    }
    #[test]
    fn regions_calibration_and_simulation_stop_at_boundaries() {
        let seqs = sequences();
        assert_eq!(
            regions_in_sequences(&[1, 1, 1], &seqs),
            vec![(0, 2), (2, 3)]
        );
        let frequencies = [Some((100, 100)); 3];
        assert_eq!(
            duration_defaults_sequences(&[1, 1, 1], &frequencies, &seqs),
            ([1.0, 2.0, 4.0], false)
        );
        let singles: Vec<_> = (0..3)
            .map(|i| Sequence {
                id: i.to_string(),
                start: i,
                end: i + 1,
            })
            .collect();
        assert_eq!(
            duration_defaults_sequences(&[1, 1, 1], &frequencies, &singles),
            ([2.0, 4.0, 8.0], true)
        );
        let p = Parameters {
            means: [0.2, 10.0],
            alpha: 0.5,
            background: 20.0,
            elevated: 4.0,
        };
        let expected = [simulate(2, p, 42), simulate(1, p, 42 + 1_000_000_007)].concat();
        assert_eq!(simulate_sequences(&seqs, p, 42), expected);
        assert_eq!(
            simulate_sequences(&Sequence::single(100), p, 42),
            simulate(100, p, 42)
        );
    }
    #[test]
    fn single_sequence_retains_fit_and_constant_segments_can_be_pooled() {
        let y = [0, 0, 0, 9, 9, 9];
        let legacy = fit(&y, 4.0, (1.001, 1e6));
        let single = fit_sequences(&y, &Sequence::single(y.len()), 4.0, (1.001, 1e6));
        assert_eq!(legacy.posterior, single.posterior);
        assert_eq!(legacy.states, single.states);
        let seqs = vec![
            Sequence {
                id: "zero".into(),
                start: 0,
                end: 3,
            },
            Sequence {
                id: "high".into(),
                start: 3,
                end: 6,
            },
        ];
        let pooled = fit_sequences(&y, &seqs, 4.0, (1.001, 1e6));
        assert!(pooled.best.is_some(), "{}", pooled.diagnostic);
        assert_eq!(pooled.starts.len(), 12);
        assert_eq!(pooled.posterior.len(), 6);
    }
}
