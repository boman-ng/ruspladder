//! Compatibility for NumPy 2.2.6's default int64 argsort. Equal-key order is
//! visible in graph vertices and therefore cannot be replaced by stable sort.

#[cfg(target_arch = "x86_64")]
unsafe extern "C" {
    fn ruspladder_argsort_i64_avx2(values: *mut i64, indices: *mut usize, len: usize);
    fn ruspladder_argsort_i64_avx512(values: *mut i64, indices: *mut usize, len: usize);
    fn ruspladder_argsort_f64_avx2(values: *mut f64, indices: *mut usize, len: usize);
    fn ruspladder_argsort_f64_avx512(values: *mut f64, indices: *mut usize, len: usize);
}

/// Final test p-values have already had NaNs replaced by 1.
pub fn argsort_f64(values: &[f64]) -> Vec<usize> {
    assert!(values.iter().all(|x| !x.is_nan()));
    let mut indices: Vec<_> = (0..values.len()).collect();
    if values.len() < 2 {
        return indices;
    }
    #[cfg(target_arch = "x86_64")]
    {
        let avx512 = std::is_x86_feature_detected!("avx512f")
            && std::is_x86_feature_detected!("avx512cd")
            && std::is_x86_feature_detected!("avx512bw")
            && std::is_x86_feature_detected!("avx512dq")
            && std::is_x86_feature_detected!("avx512vl");
        if avx512 || std::is_x86_feature_detected!("avx2") {
            let mut data = values.to_vec();
            unsafe {
                if avx512 {
                    ruspladder_argsort_f64_avx512(
                        data.as_mut_ptr(),
                        indices.as_mut_ptr(),
                        values.len(),
                    );
                } else {
                    ruspladder_argsort_f64_avx2(
                        data.as_mut_ptr(),
                        indices.as_mut_ptr(),
                        values.len(),
                    );
                }
            }
            return indices;
        }
    }
    scalar_argsort(values, &mut indices);
    indices
}

pub fn argsort_i64(values: &[i64]) -> Vec<usize> {
    let mut indices: Vec<_> = (0..values.len()).collect();
    if values.len() < 2 {
        return indices;
    }
    #[cfg(target_arch = "x86_64")]
    {
        let avx512 = std::is_x86_feature_detected!("avx512f")
            && std::is_x86_feature_detected!("avx512cd")
            && std::is_x86_feature_detected!("avx512bw")
            && std::is_x86_feature_detected!("avx512dq")
            && std::is_x86_feature_detected!("avx512vl");
        if avx512 || std::is_x86_feature_detected!("avx2") {
            // The native API accepts a mutable data pointer. Own its buffer
            // rather than casting away Rust's immutable borrow.
            let mut data = values.to_vec();
            // SAFETY: both arrays contain len elements, usize is size_t on
            // x86_64, and the instructions are gated by runtime CPU detection.
            unsafe {
                if avx512 {
                    ruspladder_argsort_i64_avx512(
                        data.as_mut_ptr(),
                        indices.as_mut_ptr(),
                        values.len(),
                    );
                } else {
                    ruspladder_argsort_i64_avx2(
                        data.as_mut_ptr(),
                        indices.as_mut_ptr(),
                        values.len(),
                    );
                }
            }
            return indices;
        }
    }
    scalar_argsort(values, &mut indices);
    indices
}

// Port of NumPy 2.2.6 npysort/quicksort.cpp aquicksort_, by Charles R. Harris
// and the NumPy developers (BSD-3-Clause; see licenses/NumPy-BSD.txt).
fn scalar_argsort<T: PartialOrd + Copy>(values: &[T], indices: &mut [usize]) {
    if indices.len() < 2 {
        return;
    }
    let depth = 2 * (usize::BITS - 1 - indices.len().leading_zeros()) as i32;
    let mut stack = vec![(0usize, indices.len() - 1, depth)];
    while let Some((mut left, mut right, mut depth)) = stack.pop() {
        if depth < 0 {
            scalar_heapsort(values, &mut indices[left..=right]);
            continue;
        }
        while right.saturating_sub(left) > 15 {
            let middle = left + ((right - left) >> 1);
            if values[indices[middle]] < values[indices[left]] {
                indices.swap(middle, left);
            }
            if values[indices[right]] < values[indices[middle]] {
                indices.swap(right, middle);
            }
            if values[indices[middle]] < values[indices[left]] {
                indices.swap(middle, left);
            }
            let pivot = values[indices[middle]];
            let mut i = left;
            let mut j = right - 1;
            indices.swap(middle, j);
            loop {
                i += 1;
                while values[indices[i]] < pivot {
                    i += 1;
                }
                j -= 1;
                while pivot < values[indices[j]] {
                    j -= 1;
                }
                if i >= j {
                    break;
                }
                indices.swap(i, j);
            }
            indices.swap(i, right - 1);
            depth -= 1;
            if i - left < right - i {
                stack.push((i + 1, right, depth));
                right = i - 1;
            } else {
                stack.push((left, i - 1, depth));
                left = i + 1;
            }
        }
        for i in left + 1..=right {
            let index = indices[i];
            let mut j = i;
            while j > left && values[index] < values[indices[j - 1]] {
                indices[j] = indices[j - 1];
                j -= 1;
            }
            indices[j] = index;
        }
    }
}

// NumPy npysort/npysort_heapsort.h aheapsort_ uses a one-based heap.
fn scalar_heapsort<T: PartialOrd + Copy>(values: &[T], indices: &mut [usize]) {
    let mut n = indices.len();
    for left in (1..=n / 2).rev() {
        sift(values, indices, left, n);
    }
    while n > 1 {
        indices.swap(0, n - 1);
        n -= 1;
        sift(values, indices, 1, n);
    }
}

fn sift<T: PartialOrd + Copy>(values: &[T], indices: &mut [usize], mut i: usize, n: usize) {
    let index = indices[i - 1];
    while i <= n / 2 {
        let mut j = i * 2;
        if j < n && values[indices[j - 1]] < values[indices[j]] {
            j += 1;
        }
        if values[index] < values[indices[j - 1]] {
            indices[i - 1] = indices[j - 1];
            i = j;
        } else {
            break;
        }
    }
    indices[i - 1] = index;
}
