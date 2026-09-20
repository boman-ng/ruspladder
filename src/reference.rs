//! Checked, half-open FASTA queries over HTSlib's faidx API.
use anyhow::{Context, Result, ensure};
use rust_htslib::htslib;
use std::{ffi::CString, path::Path, ptr::NonNull};

pub struct Reference {
    index: NonNull<htslib::faidx_t>,
}

impl Reference {
    pub fn open(path: &Path) -> Result<Self> {
        use std::os::unix::ffi::OsStrExt;
        let name = CString::new(path.as_os_str().as_bytes())?;
        // fai_load builds the standard .fai when necessary, as pysam does.
        let index = NonNull::new(unsafe { htslib::fai_load(name.as_ptr()) })
            .with_context(|| format!("open/index FASTA {}", path.display()))?;
        Ok(Self { index })
    }

    pub fn fetch(&self, chromosome: &str, start: i64, stop: i64) -> Result<Vec<u8>> {
        ensure!(start >= 0 && stop >= start, "invalid reference interval");
        let name = CString::new(chromosome)?;
        ensure!(
            unsafe { htslib::faidx_has_seq(self.index.as_ptr(), name.as_ptr()) } != 0,
            "reference has no contig {chromosome}"
        );
        if start == stop {
            return Ok(Vec::new());
        }
        let length = unsafe { htslib::faidx_seq_len64(self.index.as_ptr(), name.as_ptr()) };
        if start >= length {
            return Ok(Vec::new());
        }
        let mut fetched = 0;
        let pointer = unsafe {
            htslib::faidx_fetch_seq64(
                self.index.as_ptr(),
                name.as_ptr(),
                start,
                stop - 1,
                &mut fetched,
            )
        };
        let pointer = NonNull::new(pointer).context("fetch reference bases")?;
        if fetched < 0 {
            unsafe { libc::free(pointer.as_ptr().cast()) };
            anyhow::bail!("reference sequence retrieval failed");
        }
        let bases =
            unsafe { std::slice::from_raw_parts(pointer.as_ptr().cast::<u8>(), fetched as usize) }
                .to_vec();
        unsafe { libc::free(pointer.as_ptr().cast()) };
        Ok(bases)
    }
}

impl Drop for Reference {
    fn drop(&mut self) {
        unsafe { htslib::fai_destroy(self.index.as_ptr()) };
    }
}

pub fn reverse_complement(sequence: &[u8]) -> Result<Vec<u8>> {
    sequence
        .iter()
        .rev()
        .map(|&base| {
            Ok(match base {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                b'a' => b't',
                b't' => b'a',
                b'c' => b'g',
                b'g' => b'c',
                b'N' | b'n' => base,
                _ => anyhow::bail!("unsupported reverse-complement base {}", base as char),
            })
        })
        .collect()
}
