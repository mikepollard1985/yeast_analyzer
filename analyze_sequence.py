#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description='Analyze nucleotide composition in sliding windows')
    parser.add_argument('-f', '--file', required=True, type=str,
                      help='Path to input FASTA/FASTQ file')
    parser.add_argument('-w', '--window', type=int, default=100,
                      help='Window size for sliding window analysis (default: 100)')
    parser.add_argument('-s', '--step', type=int, default=25,
                      help='Step size for sliding window (default: 25, must not exceed window size)')
    parser.add_argument('--start', type=int,
                      help='Start position for analysis (optional, 1-based)')
    parser.add_argument('--end', type=int,
                      help='End position for analysis (optional, 1-based)')
    return parser.parse_args()

def read_sequence_file(file_path):
    """Read sequence from FASTA/FASTQ file and return the sequence string."""
    with open(file_path) as f:
        # Join all lines that don't start with markers, after stripping whitespace
        sequence = ''.join(line.strip() for line in f
                         if not line.startswith('>')
                         and not line.startswith('@')
                         and not line.startswith('+'))
    return sequence

def calculate_stats(sequence, window_size, step_size):
    """Calculate G+T/window_size and G/T ratios in sliding windows."""
    seq_length = len(sequence)
    positions = []
    gt_totals = []
    g_t_ratios = []

    for i in range(0, seq_length - window_size + 1, step_size):
        window = sequence[i:i + window_size]
        g_count = window.count('G')
        t_count = window.count('T')

        # Calculate statistics
        gt_total = (g_count + t_count) / window_size
        # Avoid division by zero
        g_t_ratio = g_count / t_count if t_count > 0 else 0

        positions.append(i + window_size // 2)  # Use midpoint of window
        gt_totals.append(gt_total)
        g_t_ratios.append(g_t_ratio)

    return positions, gt_totals, g_t_ratios

def main():
    args = parse_args()
    input_file = Path(args.file)

    if not input_file.exists():
        print(f"Error: File {input_file} does not exist", file=sys.stderr)
        sys.exit(1)

    if args.step > args.window:
        print(f"Error: Step size ({args.step}) cannot exceed window size ({args.window})", file=sys.stderr)
        sys.exit(1)

    # Read sequence directly from file
    sequence = read_sequence_file(input_file)

    # Handle start and end positions (convert from 1-based to 0-based indexing)
    start_pos = (args.start - 1) if args.start is not None else 0
    end_pos = args.end if args.end is not None else len(sequence)

    if start_pos < 0:
        print(f"Error: Start position ({args.start}) must be positive", file=sys.stderr)
        sys.exit(1)
    if end_pos > len(sequence):
        print(f"Error: End position ({args.end}) exceeds sequence length ({len(sequence)})", file=sys.stderr)
        sys.exit(1)
    if start_pos >= end_pos:
        print(f"Error: Start position ({args.start}) must be less than end position ({args.end})", file=sys.stderr)
        sys.exit(1)

    # Subset the sequence
    sequence = sequence[start_pos:end_pos]

    # Calculate statistics
    positions, gt_totals, g_t_ratios = calculate_stats(sequence, args.window, args.step)

    # Adjust positions to reflect the original coordinates
    positions = [p + start_pos for p in positions]

    # Create DataFrame and save to CSV
    df = pd.DataFrame({
        'Position': positions,
        'G+T/window': gt_totals,
        'G/T_ratio': g_t_ratios
    })
    csv_output = f'window{args.window}-step{args.step}-pos{start_pos+1}-{end_pos}.csv'
    df.to_csv(csv_output, index=False)

    # Create plot
    fig, ax1 = plt.subplots(figsize=(15, 6))

    # Plot G+T/window as scatter plot
    ax1.scatter(positions, gt_totals, color='black', s=10, label='G+T')
    ax1.set_xlabel('Position (midpoint of sliding window)')
    ax1.set_ylabel('T+G / total # bases')
    ax1.set_ylim(0, 1.8)

    # Plot G/T ratio
    ax2 = ax1.twinx()
    ax2.plot(positions, g_t_ratios, 'r^', label='G/T', markersize=3)
    ax2.set_ylabel('G/T', color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.set_ylim(0, 1.8)

    # Add legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='lower right')

    plt.title('Chromosome IX')
    plt_output = f"window{args.window}-step{args.step}-pos{start_pos+1}-{end_pos}.png"
    plt.savefig(plt_output, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Analysis complete. Results saved to:")
    print(f"  CSV: {csv_output}")
    print(f"  Plot: {plt_output}")

if __name__ == '__main__':
    main()
