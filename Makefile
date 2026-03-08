CONDA_ENV  = samflow
TESTS_DIR  = tests
LOG_DIR    = tests/logs

# Run all tests, print results to terminal
test:
	Rscript -e "testthat::test_dir('$(TESTS_DIR)/', reporter = 'progress')"
	conda run -n $(CONDA_ENV) pytest $(TESTS_DIR)/ -v

# Run all tests and save output to timestamped log files
test-log:
	@mkdir -p $(LOG_DIR)
	Rscript -e "testthat::test_dir('$(TESTS_DIR)/', reporter = 'progress')" \
		2>&1 | tee $(LOG_DIR)/r_tests_$$(date +%Y%m%d_%H%M%S).log
	conda run -n $(CONDA_ENV) pytest $(TESTS_DIR)/ -v \
		2>&1 | tee $(LOG_DIR)/py_tests_$$(date +%Y%m%d_%H%M%S).log

# Run only Tier 1 (fast, no dependencies)
test-unit:
	Rscript -e "testthat::test_file('$(TESTS_DIR)/test_rename.R', reporter = 'progress')"
	Rscript -e "testthat::test_file('$(TESTS_DIR)/test_merge.R', reporter = 'progress')"

.PHONY: test test-log test-unit
