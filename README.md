# Multivariate Outlier Detection using Hotelling's T² Statistic

A comprehensive R implementation of multivariate outlier detection using Hotelling's T-squared statistic for statistical process control and anomaly detection.

**Author**: Sharath Kumar Dhamodaran

## Overview

This repository implements a two-phase multivariate outlier detection (MOD) system based on Hotelling's T² statistic. The methodology is designed to identify anomalous observations in multivariate datasets by considering the correlation structure between variables, making it more powerful than univariate outlier detection methods.

## Key Features

- **Two-Phase Detection**: Separate modeling (Phase I) and monitoring (Phase II) phases
- **Multiple Bootstrap Methods**: Support for standard and modified bootstrapping for small sample sizes
- **Contribution Analysis**: Identifies which variables contribute most to out-of-control observations
- **Batch Decision Rules**: Determines if entire batches are out-of-control using statistical hypothesis testing
- **Flexible Covariance Estimation**: Standard and successive difference methods
- **Weighted Mahalanobis Distance**: Optional variable weighting
- **Comprehensive Visualization**: Hotelling T² plots and contribution charts

## Methodology

### Phase I: Model Development
- Estimates the mean vector (x̄) and covariance matrix (S) from training data
- Removes outliers using adjusted boxplot rules with Mahalanobis distance
- Calculates Upper Control Limit (UCL) via:
  - Direct calculation from training data
  - Modified bootstrapping for small samples
  - Standard bootstrapping with parallel processing support
- Outputs model parameters for Phase II scoring

### Phase II: Monitoring & Detection
- Scores new observations using Phase I parameters
- Calculates Hotelling's T² statistic for each observation
- Flags observations exceeding UCL as Out-of-Control (OOC)
- Performs contribution analysis to identify root causes
- Applies batch-level decision rules (alpha/beta risk control)

## Installation

### Prerequisites

This project requires R and the following packages:

```r
# Data manipulation
install.packages(c("tidyverse", "tidyr", "data.table", "gdata", "bit64"))

# Visualization
install.packages(c("ggplot2", "ggrepel", "gridExtra", "scales", "ggthemes",
                   "viridis", "pheatmap", "RColorBrewer", "reshape2"))

# Time series & analysis
install.packages(c("tseries", "tidyquant", "caTools"))

# Utilities
install.packages("here")
```

### Setup

1. Clone this repository:
```bash
git clone https://github.com/yourusername/Multivariate-Outlier-Detection.git
cd Multivariate-Outlier-Detection
```

2. Ensure your data files are in the `data/` directory

## Usage

### Quick Start

```r
# Load libraries and functions
source("scripts/libraries.R")
source("scripts/functions/get_max_pct.R")
source("scripts/functions/mod_p1.R")
source("scripts/functions/mod_p2.R")
source("scripts/functions/boxplot_outliers.R")
source("scripts/functions/contribution_plot.R")
source("scripts/functions/get_params_to_remove.R")
source("scripts/functions/t2_plot.R")

# Load and prepare your data
# (See multivariate_outlier_detection.R for full example)

# Phase I: Build the model
p1 <- mod_p1(in.data = train,
             in.vars = params,
             in.id_var = "CancerPatients",
             UCL_pct = 99,
             REMOVE_OUTLIER = "Y")

# Phase II: Score new data
p2 <- mod_p2(p1, in.score_data = test)

# Visualize results
t2_plot(p2)
plots <- contribution_plot(p2)
plots$plot1
```

### Function Reference

#### `mod_p1()` - Phase I Model Development
**Parameters:**
- `in.data`: Input data frame (no missing values)
- `in.vars`: Character vector of numeric variable names
- `in.id_var`: Unique identifier variable name
- `in.iter`: Bootstrap iterations (default: 1000)
- `METHOD`: 0 = no bootstrap, 1 = modified bootstrap, 2 = standard bootstrap
- `REMOVE_OUTLIER`: "Y" to remove outliers, "N" to keep all
- `cov.est.method`: "ST" for standard covariance, "SD" for successive difference
- `UCL_pct`: Percentile for UCL (e.g., 99 for 99% limit)
- `weight`: Optional weight vector for variables

**Returns:** List containing covariance matrix, inverse covariance, mean vector, UCL, and model parameters

#### `mod_p2()` - Phase II Scoring
**Parameters:**
- `p1`: Output from mod_p1()
- `in.score_data`: New data to score
- `na.rm`: Remove missing values (default: TRUE)
- `batch_alpha`: Type I error for batch decision (default: 0.01)
- `batch_beta`: Type II error for batch decision (default: 0.01)

**Returns:** List with scores, contributions, OOC flags, and contribution rankings

#### `t2_plot()` - Hotelling T² Visualization
**Parameters:**
- `p2`: Output from mod_p2()
- `show_ID_VAR`: Show ID labels on x-axis (default: FALSE)

**Returns:** ggplot2 object with Hotelling T² control chart

## Project Structure

```
Multivariate Outlier Detection/
├── data/
│   └── sample.csv                          # Sample dataset
├── scripts/
│   ├── libraries.R                         # Package dependencies
│   ├── multivariate_outlier_detection.R    # Main analysis script
│   └── functions/
│       ├── mod_p1.R                        # Phase I modeling
│       ├── mod_p2.R                        # Phase II scoring
│       ├── t2_plot.R                       # T² visualization
│       ├── contribution_plot.R             # Contribution charts
│       ├── boxplot_outliers.R              # Outlier removal
│       ├── get_params_to_remove.R          # Data preprocessing
│       ├── get_max_pct.R                   # Utility function
│       └── unknown_dec_rule.R              # Batch decision rules
└── README.md
```

## Example Workflow

The [multivariate_outlier_detection.R](scripts/multivariate_outlier_detection.R) script demonstrates a complete analysis:

1. **Data Loading**: Reads all CSV files from `data/` directory
2. **Preprocessing**: Removes problematic columns and missing values
3. **Train/Test Split**: 95% training, 5% testing
4. **Phase I**: Builds model with outlier removal and 99% UCL
5. **Phase II**: Scores test data and identifies OOC observations
6. **Visualization**: Generates contribution plots and T² charts

## Statistical Background

**Hotelling's T² Statistic:**
```
T² = (x - x̄)' S⁻¹ (x - x̄)
```

Where:
- x = observation vector
- x̄ = mean vector from Phase I
- S⁻¹ = inverse covariance matrix from Phase I

An observation is flagged as OOC if T² > UCL.

**Contribution Analysis:**
For each variable i, the contribution is calculated as:
```
D²ᵢ = T² - T²₍₋ᵢ₎
```
Where T²₍₋ᵢ₎ is calculated excluding variable i.

## Applications

- **Quality Control**: Manufacturing process monitoring
- **Anomaly Detection**: Fraud detection, system monitoring
- **Healthcare**: Patient outcome monitoring
- **Environmental**: Pollution monitoring
- **Finance**: Risk management and outlier detection

## License

This project is open source and available under the MIT License.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## References

- Hotelling, H. (1947). "Multivariate Quality Control"
- Mason, R.L., & Young, J.C. (2002). "Multivariate Statistical Process Control with Industrial Applications"
- Tracy, N.D., Young, J.C., & Mason, R.L. (1992). "Multivariate Control Charts for Individual Observations"

## Contact

For questions or issues, please open an issue on GitHub or contact the author.

---

**Note**: This implementation uses a sample cancer patient dataset for demonstration purposes. Replace with your own data for specific applications.
