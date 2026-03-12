# Soybean-Phenology-QA-KPI-Reporting-System
This project is an automated quality assurance and KPI reporting pipeline developed to support the monitoring of large-scale agronomic trial datasets during the crop season.

In multinational agronomic operations, hundreds of field trials are conducted every year, generating large volumes of heterogeneous and operationally complex data. Ensuring the quality, completeness, and consistency of this information is a major challenge—especially when agronomic teams must also manage multiple parallel responsibilities across regions, crops, and stakeholders.
To address this, I developed a Python-based reporting system that:

Ingests and standardizes agronomic trial data from multiple Excel sources
Applies automated QA checks to phenology, disease severity, and row-closing observations
Tracks within-season KPIs for data quality and operational follow-up
Produces decision-ready outputs for local and global stakeholder discussions
Converts what would normally require weeks or months of manual review into a reproducible process executed in seconds

The solution was initially built for soybean phenology trial monitoring in Brazil, but its structure is modular and can be adapted to other crops and protocols with minimal edits.

# Business Context
Agronomic trial programs operate with high data volume, regional variability, and strict timing requirements. During the season, stakeholders need to quickly understand:

whether key assessments are being collected correctly,
whether growth stage information is complete and biologically consistent,
where operational gaps exist,
which regions or teams require follow-up,
and whether datasets are reliable enough to support downstream analytics and model development.

Manual QA of this process is slow, error-prone, and difficult to scale.
This script was created to solve that problem by turning raw trial exports into a structured QA and performance monitoring system, enabling teams to detect issues early, prioritize corrections, and align next steps across local and global teams.

# 3) What the script does (technical summary)
Core functionality
The script performs an end-to-end workflow for soybean phenology trial monitoring:
1. Data ingestion and preprocessing

Loads trial and seeding data from Excel files
Standardizes column names and formats
Filters the dataset by farm/season and protocol
Normalizes phenological growth stage values
Extracts operational metadata such as:

responsible person
Brazilian state / region
observation date
seeding date
DAP (Days After Planting)



# 2. Automated QA validation
The script includes several domain-specific QA checks, such as:


Growth stage trio completeness
Verifies whether minimum, major, and maximum growth stage fields are fully populated.


Growth stage internal consistency
Ensures min ≤ maj ≤ max within each record.


Growth stage regression detection
Detects biologically impossible backward progression over time.


Development gap detection
Flags jumps greater than a defined threshold between consecutive observations.


Assessment package completeness
Confirms the presence of required disease and row-closing assessments.


Disease monotonicity checks
Validates non-decreasing disease progression where biologically expected.


Row-closing monotonicity and progress checks
Tracks closure progression and identifies open backlog trajectories.


# 3. KPI calculation
The script computes composite KPIs at different aggregation levels, such as:

by responsible person
by region/state
by field
by trajectory

The KPI framework combines:

GS trio completeness
assessment package completeness
growth stage coverage
row-closing completion
inverse QA issue burden

This creates a practical scoring system for comparing data quality and operational execution across teams and geographies.

# 4. Visualization and reporting
The pipeline automatically generates:

heatmaps
coverage charts
QA issue distributions
DAP vs. growth stage analysis
progression curves
backlog summaries
technical PDF reports
Excel workbooks segmented by responsible and region

These outputs support both technical validation and stakeholder communication.

# 4) Overall impact (business value section)
## Impact
The main value of this work is its ability to transform a highly manual and fragmented process into a fast, standardized, and scalable analytical workflow.
## Operational impact

Reduced a process that could take months of manual validation and consolidation into a workflow that runs in ~30 seconds
Enabled within-season monitoring, not just post-season correction
Improved visibility of data collection gaps before they became critical downstream problems
Created a reproducible framework that can be reused across seasons and adapted to other crops

## Analytical impact

Increased confidence in phenology and assessment data used for reporting and analytics
Improved traceability of field-level and trajectory-level issues
Supported more reliable discussions on data status, quality, and next steps

## Organizational impact

Facilitated alignment between local operational teams and global stakeholders
Provided a structured basis for leadership discussions around trial progress and data readiness
Scaled the value of agronomic data analysis from a local operational task to a global decision-support capability


# 5) Why this project matters
This project demonstrates how domain knowledge in agronomy can be combined with data engineering, QA logic, and reporting automation to solve a real business bottleneck. Rather than only analyzing results after the fact, the system helps teams monitor trial quality during the season, identify issues early, and act faster. The result is a substantial gain in efficiency, consistency, and decision quality across a global trial network.
