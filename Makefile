# Makefile for compiling MultiDraw shared library

# Variables
CXX := g++
CXXFLAGS := -fPIC -shared `root-config --cflags --libs`
SRC_DIR := src
INC_DIR := interface
LIB_DIR := lib
TARGET := $(LIB_DIR)/libMultiDraw.so
SRC := $(SRC_DIR)/MultiDraw.cc

# Automatically detect the current working directory
CWD := $(shell pwd)

# Include paths
INCLUDES := -I$(CWD)/$(INC_DIR)

# Default target
all: $(LIB_DIR) $(TARGET)

# Compile the shared library
$(TARGET): $(SRC)
	@echo "Compiling $@"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^

# Create directories if they do not exist
$(LIB_DIR):
	@echo "Creating directory $@"
	@mkdir -p $@

# Clean up generated files
clean:
	@echo "Cleaning up..."
	@rm -f $(TARGET)

# Phony targets
.PHONY: all clean
