# PreQC Package

## Why PreQC Pack?
Traditional solutions require multiple software combinations to calculate indicators for quality pre-evaluation, which is costly and time-consuming. Whether it is possible to develop a low-cost and low-resource-consumption software to generate the required pre-evaluation indicators at one time? So we develop the preqc-pack.

## A collection of QC Metrics for preassessment
### File Metadata
- MD5SUM
- File Size

### FastQC
- Data Size

## Installation
Build the preqc-pack on macOS.
```
# Clone the source code
git clone https://github.com/clinico-omics/preqc-pack.git

# Change the directory
cd preqc-pack

# For Linux
cargo build --release --target=x86_64-unknown-linux-musl

# For macOS
cargo build --release
```
