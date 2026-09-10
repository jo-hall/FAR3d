BUILD_DIR := build

.PHONY: all clean

all: $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake .. && $(MAKE)

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

clean:
	@echo "Cleaning BOOZ_XFORM build folder"
	@rm -rf $(BUILD_DIR)
