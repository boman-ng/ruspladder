// NumPy 2.2.6 numpy/_core/src/umath/loops_utils.h.src DOUBLE_pairwise_sum.
// BSD-3-Clause; see licenses/NumPy-BSD.txt. Preserve its addition order.
pub fn sum_strided(values: &[f64]) -> f64 {
    values.iter().fold(0.0, |sum, &x| sum + x)
}

pub fn sum(values: &[f64]) -> f64 {
    if values.len() < 8 {
        return values.iter().fold(-0.0, |total, &x| total + x);
    }
    if values.len() <= 128 {
        let mut accumulators: [f64; 8] = values[..8].try_into().unwrap();
        let end = values.len() - values.len() % 8;
        for block in values[8..end].chunks_exact(8) {
            for i in 0..8 {
                accumulators[i] += block[i];
            }
        }
        let r = accumulators;
        let total = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        return values[end..].iter().fold(total, |total, &x| total + x);
    }
    let half = (values.len() / 2) / 8 * 8;
    sum(&values[..half]) + sum(&values[half..])
}
