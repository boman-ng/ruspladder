fn main() {
    println!("cargo:rerun-if-env-changed=RUSPLADDER_BLAS_DIR");
    let blas = std::env::var("RUSPLADDER_BLAS_DIR")
        .expect("set RUSPLADDER_BLAS_DIR to the pinned OpenBLAS library directory");
    println!("cargo:rustc-link-search=native={blas}");
    println!("cargo:rustc-link-lib=dylib=scipy_openblas64_-56d6093b");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{blas}");
    println!("cargo:rerun-if-changed=native/special.cpp");
    println!("cargo:rerun-if-changed=vendor/scipy-special");
    cc::Build::new()
        .cpp(true)
        .std("c++17")
        .opt_level(2)
        .flag("-ffp-contract=off")
        .flag("-isystem")
        .flag("vendor/scipy-special")
        .file("native/special.cpp")
        .compile("ruspladder_special");
    println!("cargo:rerun-if-changed=native/argsort.cpp");
    println!("cargo:rerun-if-changed=vendor/x86-simd-sort/src");
    if std::env::var("CARGO_CFG_TARGET_ARCH").as_deref() == Ok("x86_64") {
        for (name, flags) in [
            ("avx2", vec!["-mavx2"]),
            (
                "avx512",
                vec![
                    "-mavx512f",
                    "-mavx512cd",
                    "-mavx512bw",
                    "-mavx512dq",
                    "-mavx512vl",
                ],
            ),
        ] {
            let symbol = format!("ruspladder_argsort_i64_{name}");
            let mut build = cc::Build::new();
            build
                .cpp(true)
                .std("c++17")
                .opt_level(2)
                .include("vendor/x86-simd-sort/src")
                .file("native/argsort.cpp")
                .define("RUSPLADDER_ARGSORT", symbol.as_str());
            for flag in flags {
                build.flag(flag);
            }
            build.compile(&format!("ruspladder_sort_{name}"));
        }
    }
}
