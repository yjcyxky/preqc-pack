# PreQC Pack: A quality control tool for high throughput sequence data.

Check out the [User Guide]() for a list of features and installation and usage information. The User Guide also serves as a demonstration to showcase what a book looks like.

If you are interested in contributing to the development of mdBook, check out the Contribution Guide.

## Why PreQC Pack?
Traditional solutions require multiple software combinations to calculate indicators for quality pre-evaluation, which is costly and time-consuming. Whether it is possible to develop a low-cost and low-resource-consumption software to generate the required pre-evaluation indicators at one time? So we develop the preqc-pack.

## A collection of QC Metrics for NGS
### Quality Metrics
- checksum
- file size

### FastQC

### Mislabeling (NGSCheckMate)

## Performance

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

## License
All the code in this repository is released under the Mozilla Public License v2.0, for more information take a look at the LICENSE file.