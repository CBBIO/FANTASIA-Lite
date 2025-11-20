#!/bin/bash

# FANTASIA Lite V0 - Setup and Test Script
# Downloads lookup bundle from Zenodo and runs validation test

set -e  # Exit on any error

echo "========================================="
echo "FANTASIA Lite V0 - Setup and Test"
echo "========================================="
echo ""

# Check if wget or curl is available
if command -v wget &> /dev/null; then
    DOWNLOAD_CMD="wget -O"
elif command -v curl &> /dev/null; then
    DOWNLOAD_CMD="curl -L -o"
else
    echo "Error: Neither wget nor curl is available. Please install one of them."
    exit 1
fi

# Create data/lookup directory if it doesn't exist
echo "[1/4] Creating data/lookup directory..."
mkdir -p data/lookup

# Download the lookup bundle
ZENODO_URL="https://zenodo.org/records/17630193/files/fantasia_lite_lookup_table.zip"
BUNDLE_FILE="fantasia_lite_lookup_table.zip"

echo "[2/4] Downloading lookup bundle from Zenodo..."
echo "URL: $ZENODO_URL"
if [ -f "$BUNDLE_FILE" ]; then
    echo "Bundle file already exists. Skipping download."
else
    $DOWNLOAD_CMD "$BUNDLE_FILE" "$ZENODO_URL"
    echo "Download complete!"
fi

# Check if files already exist
echo "[3/4] Checking for existing lookup files..."
if [ -f "data/lookup/lookup_table.npz" ] && \
   [ -f "data/lookup/annotations.json" ] && \
   [ -f "data/lookup/accessions.json" ]; then
    echo "✓ lookup_table.npz"
    echo "✓ annotations.json"
    echo "✓ accessions.json"
    echo "All required files already present! Skipping extraction."
else
    # Extract the bundle
    echo "Extracting lookup bundle..."
    
    # Check if the archive contains a nested directory structure
    echo "Inspecting archive structure..."
    ARCHIVE_CONTENTS=$(tar -tzf "$BUNDLE_FILE" | head -5)
    
    # Extract to temporary location first
    TEMP_DIR=$(mktemp -d)
    tar -xzf "$BUNDLE_FILE" -C "$TEMP_DIR"
    
    # Find the actual files and move them to the correct location
    LOOKUP_NPZ=$(find "$TEMP_DIR" -name "lookup_table.npz" -type f)
    ANNOTATIONS_JSON=$(find "$TEMP_DIR" -name "annotations.json" -type f)
    ACCESSIONS_JSON=$(find "$TEMP_DIR" -name "accessions.json" -type f)
    
    if [ -n "$LOOKUP_NPZ" ] && [ -n "$ANNOTATIONS_JSON" ] && [ -n "$ACCESSIONS_JSON" ]; then
        mv "$LOOKUP_NPZ" data/lookup/
        mv "$ANNOTATIONS_JSON" data/lookup/
        mv "$ACCESSIONS_JSON" data/lookup/
        echo "Extraction complete!"
    else
        echo "Error: Could not find required files in archive."
        rm -rf "$TEMP_DIR"
        exit 1
    fi
    
    # Clean up temporary directory
    rm -rf "$TEMP_DIR"
    
    # Verify extracted files
    echo ""
    echo "Verifying extracted files:"
    if [ -f "data/lookup/lookup_table.npz" ] && \
       [ -f "data/lookup/annotations.json" ] && \
       [ -f "data/lookup/accessions.json" ]; then
        echo "✓ lookup_table.npz"
        echo "✓ annotations.json"
        echo "✓ accessions.json"
        echo "All required files present!"
    else
        echo "Error: Missing required files after extraction."
        exit 1
    fi
fi

# Clean up the downloaded archive (optional)
echo ""
read -p "Delete downloaded archive ($BUNDLE_FILE)? [y/N] " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    rm "$BUNDLE_FILE"
    echo "Archive deleted."
fi

# Run test
echo ""
echo "[4/4] Running validation test..."
echo "Command: python3 fantasia_pipeline.py --serial-models --embed-models prot_t5 fasta_test/test.fa"
echo ""
python3 fantasia_pipeline.py --serial-models --embed-models prot_t5 fasta_test/test.fa

echo ""
echo "========================================="
echo "Setup and test completed successfully!"
echo "========================================="
