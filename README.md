# DNA Sequence Analysis Tool

This tool analyzes DNA sequences from FASTA or FASTQ files, calculating and plotting nucleotide composition statistics. It's designed to work with yeast genome sequences but can be used with any DNA sequence file.

## Setup Instructions

1. Make sure you have Python 3 installed on your computer. You can download it from [python.org](https://python.org).

2. Open a terminal (Command Prompt on Windows, Terminal on Mac/Linux)

3. Create a new folder for this project and copy these files into it:
   - `analyze_sequence.py`
   - `requirements.txt`

4. Navigate to your project folder in the terminal:
   ```bash
   cd path/to/your/folder
   ```

5. Install the required packages:
   ```bash
   pip install -r requirements.txt
   ```

6. Make the script executable (Mac/Linux only):
   ```bash
   chmod +x analyze_sequence.py
   ```

## How to Use

The script takes a DNA sequence file (FASTA or FASTQ format) and analyzes it using a sliding window approach.

### Basic Usage

```bash
python analyze_sequence.py -f your_sequence.fastq
```

### All Available Options

- `-f` or `--file`: Your sequence file (required)
- `-w` or `--window`: Window size (default: 100)
- `-s` or `--step`: Step size (default: 25)
- `--start`: Start position for analysis (optional)
- `--end`: End position for analysis (optional)

### Examples

1. Basic analysis with default settings:
   ```bash
   python analyze_sequence.py -f sequence.fastq
   ```

2. Change window size to 200 and step size to 50:
   ```bash
   python analyze_sequence.py -f sequence.fastq -w 200 -s 50
   ```

3. Analyze only positions 1000 to 2000:
   ```bash
   python analyze_sequence.py -f sequence.fastq --start 1000 --end 2000
   ```

## Output

The script creates two files:
1. A PNG file with the plot
2. A CSV file with the raw data

The output files will be named based on your input parameters (window size, step size, and position range).

## Understanding the Plot

The plot shows two measurements:
- Black dots: The proportion of G+T bases in each window
- Red triangles: The ratio of G to T bases in each window

## Troubleshooting

1. If you get "command not found":
   - Make sure you're in the correct directory
   - Try using `python analyze_sequence.py` instead of `./analyze_sequence.py`

2. If you get import errors:
   - Make sure you've installed the requirements: `pip install -r requirements.txt`

3. If your file isn't found:
   - Make sure you're using the correct path to your sequence file
   - Check that the file exists and has the correct permissions

## Need Help?

If you encounter any issues:
1. Check that your sequence file is in FASTA or FASTQ format
2. Verify that all command line arguments are correct
3. Make sure you've followed all setup instructions
