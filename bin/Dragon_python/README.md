# Dragon CLI Usage

This repository includes the main sequence design program and helper scripts for preparing input and generating formatted output.

This README only covers command-line usage. The larger demo in the main directory covers how to build `pattern.txt` input files.

## Requirements

- Python 3
- A valid `pattern.txt`

## Programs

### `trace_pattern.py`

Generates `target.txt` from `pattern.txt`.

#### Usage

    python trace_pattern.py pattern.txt > target.txt

#### Input

- `pattern.txt` — pattern blueprint file

#### Output

- Writes generated target data to standard output
- Usually redirected to `target.txt`

---

### `dragon.py`

Main sequence design / mutation program.

#### Usage

    python dragon.py [folder]

#### Input Requirements

`dragon.py` requires a `target.txt` file in the working directory.

You must generate `target.txt` first by running:

    python trace_pattern.py pattern.txt > target.txt

#### Output

- Writes design results into the folder you provide
- The output folder must exist before running `dragon.py`

#### Example

    mkdir my_folder
    python dragon.py my_folder

#### Full Example

    python trace_pattern.py pattern.txt > target.txt
    mkdir my_folder
    python dragon.py my_folder

---

### `trace_analysis.py`

Generates formatted trace / structure analysis output from a pattern file and sequence file.

#### Usage

    python trace_analysis.py pattern.txt seq.txt

#### Inputs

- `pattern.txt` — pattern blueprint file
- `seq.txt` — sequence to analyze

#### Output

- Writes formatted analysis to standard output
- Can be redirected or appended to a report file

#### Example

    python trace_analysis.py pattern.txt seq.txt >> design_report.txt

---

## Basic Workflow

### 1. Generate `target.txt`

    python trace_pattern.py pattern.txt > target.txt

### 2. Create an output folder

    mkdir my_folder

### 3. Run Dragon

    python dragon.py my_folder

### 4. Run trace analysis if needed

    python trace_analysis.py pattern.txt seq.txt

---

## Minimal Example

    python trace_pattern.py pattern.txt > target.txt
    mkdir my_folder
    python dragon.py my_folder