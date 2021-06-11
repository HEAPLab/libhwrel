ifdef CONFIG_BBQUE_LIBHWREL

# Targets provided by this project
.PHONY: libhwrel clean_libhwrel

# Add this to the "external" target
external: libhwrel
clean_external: clean_libhwrel

MODULE_EXTERNAL_LIBHWREL = external/optional/libhwrel

libhwrel:
	@echo
	@echo "==== Building Hardware Reliability Library - libhwrel ($(BUILD_TYPE)) ===="
	@[ -d $(MODULE_EXTERNAL_LIBHWREL)/build/$(BUILD_TYPE) ] || \
		mkdir -p $(MODULE_EXTERNAL_LIBHWREL)/build/$(BUILD_TYPE) || \
		exit 1
	@cd $(MODULE_EXTERNAL_LIBHWREL)/build/$(BUILD_TYPE) && \
		CC=$(CC) CFLAGS=$(TARGET_FLAGS) \
		CXX=$(CXX) CXXFLAGS=$(TARGET_FLAGS) \
		cmake $(CMAKE_COMMON_OPTIONS) ../.. || \
		exit 1
	@cd $(MODULE_EXTERNAL_LIBHWREL)/build/$(BUILD_TYPE) && \
		make -j$(CPUS) install || \
		exit 1

clean_libhwrel:
	@echo
	@echo "==== Clean-up Building Hardware Reliability Library - libhwrel ($(BUILD_TYPE)) ===="
	@rm -rf $(BUILD_DIR)/lib/libhwrel* $(BUILD_DIR)/include/libhwrel
	@rm -rf $(MODULE_EXTERNAL_LIBHWREL)/build
	@echo

else # CONFIG_BBQUE_LIBHWREL

libhwrel:
	$(warning external/libhwrel module disabled by BOSP configuration)
	$(error BOSP compilation failed)

endif # CONFIG_BBQUE_LIBHWREL
