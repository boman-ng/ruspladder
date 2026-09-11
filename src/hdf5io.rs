//! Public HDF5 arrays, including NumPy-compatible fixed-width byte strings.
use anyhow::{Result, ensure};
use hdf5::{Dataset, Group, H5Type, types::TypeDescriptor};
use hdf5_sys::{
    h5d::{H5Dread, H5Dwrite},
    h5p::H5P_DEFAULT,
    h5s::H5S_ALL,
};

pub fn array<T: H5Type>(
    group: &Group,
    name: &str,
    shape: &[usize],
    values: &[T],
    gzip: bool,
) -> Result<Dataset> {
    let builder = group.new_dataset::<T>().shape(shape);
    let dataset = if gzip {
        builder.deflate(4).chunk_min_kb(64)
    } else {
        builder
    }
    .create(name)?;
    dataset.write_raw(values)?;
    Ok(dataset)
}

pub fn strings(
    group: &Group,
    name: &str,
    shape: &[usize],
    values: &[&str],
    gzip: bool,
) -> Result<()> {
    let width = values.iter().map(|v| v.len()).max().unwrap_or(1).max(1);
    let builder = group
        .new_dataset_builder()
        .empty_as(&TypeDescriptor::FixedAscii(width))
        .shape(shape);
    let dataset = if gzip {
        builder.deflate(4).chunk_min_kb(64)
    } else {
        builder
    }
    .create(name)?;
    write_strings(&dataset, values)
}

pub fn write_strings(dataset: &Dataset, values: &[&str]) -> Result<()> {
    ensure!(
        dataset.size() == values.len(),
        "string array shape mismatch"
    );
    // SplAdder's codeUTF8 uses NumPy astype('bytes'), which accepts ASCII only.
    ensure!(
        values.iter().all(|v| v.is_ascii()),
        "SplAdder's byte-string output requires ASCII identifiers"
    );
    let dtype = dataset.dtype()?;
    let TypeDescriptor::FixedAscii(width) = dtype.to_descriptor()? else {
        anyhow::bail!("expected fixed-width byte strings in {}", dataset.name());
    };
    let mut bytes = vec![0u8; width * values.len()];
    for (value, target) in values.iter().zip(bytes.chunks_exact_mut(width)) {
        let length = value.len().min(width);
        target[..length].copy_from_slice(&value.as_bytes()[..length]);
    }
    // The buffer contains exactly dataset.size() fixed-width elements. Using
    // the dataset's own datatype avoids a variable/fixed string conversion.
    hdf5::sync::sync(|| {
        hdf5::h5check(unsafe {
            H5Dwrite(
                dataset.id(),
                dtype.id(),
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                bytes.as_ptr().cast(),
            )
        })
    })?;
    Ok(())
}

pub fn read_strings(dataset: &Dataset) -> Result<Vec<String>> {
    let dtype = dataset.dtype()?;
    let TypeDescriptor::FixedAscii(width) = dtype.to_descriptor()? else {
        anyhow::bail!("expected fixed-width byte strings in {}", dataset.name());
    };
    let mut bytes = vec![0u8; dataset.size() * width];
    hdf5::sync::sync(|| {
        hdf5::h5check(unsafe {
            H5Dread(
                dataset.id(),
                dtype.id(),
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                bytes.as_mut_ptr().cast(),
            )
        })
    })?;
    bytes
        .chunks_exact(width)
        .map(|value| {
            let stop = value.iter().position(|&b| b == 0).unwrap_or(width);
            Ok(String::from_utf8(value[..stop].to_vec())?)
        })
        .collect()
}
