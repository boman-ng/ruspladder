.PHONY: build test

WORK_ROOT := $(if $(RUSPLADDER_WORK_ROOT),$(RUSPLADDER_WORK_ROOT),$(HOME)/data/ruspladder)
export RUSPLADDER_WORK_ROOT := $(WORK_ROOT)
export OPENBLAS_NUM_THREADS := 1
export OMP_NUM_THREADS := 1
export NUMBA_NUM_THREADS := 4
export RAYON_NUM_THREADS := 4
TARGET_DIR = $(if $(CARGO_TARGET_DIR),$(CARGO_TARGET_DIR),$(WORK_ROOT)/target)
TEST_WORK ?= $(WORK_ROOT)/runs/public-$(shell date +%Y%m%dT%H%M%S)
BINARY ?= $(TARGET_DIR)/release/ruspladder

build:
	python3 scripts/bootstrap-native.py
	bash scripts/cargo.sh build --locked --release --bin ruspladder

# Requires scripts/bootstrap.sh (the pinned Python reference environment).
test:
	bash scripts/cargo.sh test --locked --lib
	bash scripts/cargo.sh build --locked --release --bin ruspladder --example cache_export_probe
	$(WORK_ROOT)/envs/reference/bin/python scripts/test_public_data.py \
		--binary $(BINARY) --exporter $(TARGET_DIR)/release/examples/cache_export_probe \
		--upstream $(WORK_ROOT)/upstream/spladder --work $(TEST_WORK)
