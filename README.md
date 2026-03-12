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

## Key Features
- Data ingestion and standardization
- Automated QA checks
- KPI computation
- Visual analytics
- Excel and PDF export

## QA Checks Implemented
- GS trio completeness
- GS triplet consistency
- GS regression
- Development gaps
- Assessment package completeness
- Disease monotonicity
- RC monotonicity and backlog tracking

## Outputs
- Consolidated KPI workbooks
- Region-level Excel reports
- Responsible-level Excel reports
- PDF technical report
- Automated plots and heatmaps

## Impact
[Use the impact section]

## Tech Stack
Python, pandas, numpy, matplotlib, seaborn, reportlab, openpyxl

## Notes
This repository contains an anonymized/adapted version of an internal operational workflow developed for large-scale agronomic trial monitoring.


# 5) Why this project matters
This project demonstrates how domain knowledge in agronomy can be combined with data engineering, QA logic, and reporting automation to solve a real business bottleneck. Rather than only analyzing results after the fact, the system helps teams monitor trial quality during the season, identify issues early, and act faster. The result is a substantial gain in efficiency, consistency, and decision quality across a global trial network.
