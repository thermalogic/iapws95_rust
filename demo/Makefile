# IAPWS-95 C Interface Build Script
# 
# Usage:
#   make build    - Build the Rust library and compile the example
#   make run      - Build and run the example
#   make clean    - Clean build artifacts

# Configuration
RUST_TARGET := release
CFLAGS := -Wall -Wextra
INCLUDES := -I../src

# Platform detection
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
    LIB_NAME := libiapws95.so
    LD_FLAGS := -L../target/$(RUST_TARGET) -liapws95
    EXAMPLE_OUT := iapws95_demo
else ifeq ($(UNAME_S),Darwin)
    LIB_NAME := libiapws95.dylib
    LD_FLAGS := -L../target/$(RUST_TARGET) -liapws95
    EXAMPLE_OUT := iapws95_demo
else
    # Windows
    LIB_NAME := iapws95.lib
    MSVC_PATH := $(shell where cl 2>nul)
    EXAMPLE_OUT := iapws95_demo.exe
endif

.PHONY: all build run clean

all: build run

# Build Rust library and compile C example
build:
	@echo "Building IAPWS-95 Rust library..."
	cargo build --$(RUST_TARGET)
	@echo ""
ifeq ($(UNAME_S),Linux)
	@echo "Compiling C example (Linux)..."
	gcc $(CFLAGS) $(INCLUDES) c_example.c $(LD_FLAGS) -o $(EXAMPLE_OUT)
else ifeq ($(UNAME_S),Darwin)
	@echo "Compiling C example (macOS)..."
	gcc $(CFLAGS) $(INCLUDES) c_example.c $(LD_FLAGS) -o $(EXAMPLE_OUT)
else
	@echo "Compiling C example (Windows MSVC)..."
	@if exist "$(MSVC_PATH)" ( \
		cl /I..\src c_example.c ..\target\$(RUST_TARGET)\iapws95.lib /Fe:$(EXAMPLE_OUT) \
	) else ( \
		echo "Error: MSVC compiler not found. Please use Developer Command Prompt." \
	)
endif

# Build and run
run: build
	@echo ""
	@echo "Running example..."
ifeq ($(UNAME_S),Windows)
	.\$(EXAMPLE_OUT)
else
	./$(EXAMPLE_OUT)
endif

# Clean build artifacts
clean:
	@echo "Cleaning..."
	cargo clean
	@if exist "$(EXAMPLE_OUT)" del $(EXAMPLE_OUT)
	@echo "Done."
