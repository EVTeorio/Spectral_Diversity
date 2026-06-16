# Governance Standards

## Working Directory

C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity

## R and Rtools Access

The agent is authorized to discover, read, and execute installed R and Rtools binaries, libraries, package metadata, and supporting files anywhere on this computer, even when those files are outside the designated repository directory.

This exception exists only to support reproducible R execution, package discovery, package compilation, and dependency validation. The agent must not modify, delete, reinstall, or reconfigure R, Rtools, system libraries, system paths, or system environment variables without explicit user approval for that specific action.

When R or Rtools are not available through the active shell path, the agent may search common installation locations and use absolute paths to R/Rtools executables for validation and testing.

## Repository Discovery Requirements

Before performing any work, the agent must create or update the following project inventories.

### Directory Map

**Location:**
`/reports/directory_map.md`

**Purpose:**
Provide a complete hierarchical map of the repository structure.

**Include:**
- folders
- scripts
- functions
- reports
- outputs
- tests
- documentation

**Format:**
```text
project_root/
data/
    raw/
    intermediate/
    processed/

scripts/
    01_ingest/
    02_clean/
    03_transform/

outputs/

reports/
```

Update whenever the repository structure changes.

---

### Data Dictionary

**Location:**
`/reports/data_dictionary.md`

**Purpose:**
Provide a living inventory of all datasets used in the project.

#### For every dataset include:

- Dataset Name
- Description
- Source Location
- Row Count
- Column Count
- Primary Keys
- Date Fields
- Numeric Fields
- Categorical Fields
- Missing Value Summary
- Last Updated

#### For every variable include:

| Variable | Type | Description |
|----------|------|-------------|

The agent should update the data dictionary whenever:

- new datasets are introduced
- schemas change
- fields are added or removed

---

## Coding Style & Naming Conventions

### Language

Use R for all generated code unless explicitly instructed otherwise.

Generated code should remain compatible with:

`R 4.2.3`

---

### Style

Use:

- 2-space indentation
- meaningful object names
- modular functions
- concise logic

#### Requirements:

- Keep functions focused on a single responsibility.
- Prefer maintainable solutions over clever solutions.
- Add comments where logic is not obvious.
- Avoid deeply nested logic.
- Avoid duplicated code.

---

### Naming Standards

#### Functions

`snake_case()`

Examples:

```r
clean_customer_data()
calculate_growth_rate()
```

#### Variables

`snake_case`

Examples:

```r
customer_count
sales_total
```

#### Constants

`UPPER_SNAKE_CASE`

Examples:

```r
MAX_ITERATIONS
DEFAULT_LOOKBACK_DAYS
```

#### Files

`snake_case.R`

Examples:

```text
clean_customer_data.R
customer_summary_report.R
```

#### Directories

`lowercase_snake_case`

Examples:

```text
raw_data
processed_data
analysis_results
```

---

## Execution Planning Requirements

For complex requests, create an ExecPlan before making changes.

An ExecPlan should contain:

- Objective
- Requested task
- Files To Review
- Relevant files
- Proposed Changes
- Expected modifications
- Validation Plan
- How correctness will be verified
- Risks
- Potential side effects

**Store plans under:**

`/reports/execution_plans/`

**Naming:**

`YYYYMMDD_HHMM_execplan.md`

---

## Testing Standards

Testing is mandatory for:

- new features
- bug fixes
- refactors affecting logic
- cleaning or optimizing R scripts

**Preferred framework:**

`testthat`

**Location:**

`/tests`

**Naming:**

`test_feature_name.R`

### Requirements

- use small representative subsamples of data for unit tests
- test normal conditions
- test edge cases
- test failure conditions
- test missing data scenarios
- test invalid inputs

When behavior changes:

- update existing tests
- update documentation

---

## Documentation Requirements

Documentation is considered part of task completion.

Update documentation whenever:

- workflows change
- outputs change
- assumptions change
- dependencies change

Documentation locations:

- `/README.md`
- `/reports/`
- `/docs/`

---

## Rules

The agent must not:

1. Modify system paths.
2. Modify system environment variables.
3. Modify files outside repository boundaries, except for read/execute access to R and Rtools files as described above.
4. Commit code.
5. Push code.
6. Rebase branches.
7. Delete files without approval.

The agent may:

- recommend Git commands
- generate commit messages
- generate release notes

The user remains responsible for all Git actions.

---

## Repository Organization Rules

- Autonomous experimentation: `/autoresearch`
- Architecture notes: `/reports/architecture`
- Historical summaries: `/reports/history`
- Task reports: `/reports/tasks`
- Execution plans: `/reports/execution_plans`
- Validation reports: `/reports/validation`
- Logs: `/logs`
- Unit tests: `/tests`

---

## Agent Expectations

### Before Making Changes

1. Review relevant repository files.
2. Understand current implementation.
3. Review project_state.md.
4. Review data_dictionary.md.
5. Review directory_map.md.
6. Create an ExecPlan when required.
7. Identify dependencies.
8. Confirm assumptions.

### Current Data Assumptions

1. Unnamed first columns in CSV outputs are likely index columns unless later evidence shows otherwise.
2. Spectral products that produce spectral heterogeneity values are stored in `/Indices_SHPs/Spectral_diversitySHPs`.
3. Spectra cropped to quadrat extents are stored in `/Quad_Spectra`.
4. Current analysis should use current, non-legacy data products. Files under directories named `old`, `Outdated`, or `Currently not relevant` are authorized for deletion when cleanup is requested.

### During Work

1. Make minimal targeted changes.
2. Preserve existing functionality.
3. Update tests.
4. Update documentation.
5. Record actions performed.

### After Work

1. Validate outputs.
2. Run tests.
3. Update data dictionary if needed.
4. Update directory map if needed.
5. Create task report.
6. Create validation report.
7. Summarize work completed.

---

## Project State Tracking

Maintain:

`/reports/project_state.md`

### Purpose

Provide project continuity between agent sessions.

### Required Sections

- Current Objective
- Active Research Questions
- Completed Work
- Pending Work
- Known Issues
- Technical Debt
- Important Assumptions
- Next Recommended Actions

Update after every major task.

---

## Completion Checklist

A task is not complete until all applicable items are complete.

Checklist:

- [ ] Objective completed
- [ ] Code updated
- [ ] Tests added or updated
- [ ] Validation completed
- [ ] Documentation updated
- [ ] Data dictionary updated
- [ ] Directory map updated
- [ ] Project state updated
- [ ] Task report created
- [ ] Validation report created
- [ ] Summary provided to user
