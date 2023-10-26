import pandas as pd
import numpy as np
import argparse
import logging
import sys
from rich.console import Console
from rich.progress import Progress


def gc_content(seq):
    """Calculate GC content of a DNA sequence."""
    return (seq.count('G') + seq.count('C')) / len(seq)


def calculate_y_pred(original, variant, gc_weight, params):
    """Calculate y_pred based on original and variant sequences."""

    # print(original, variant)

    if original is None or variant is None or original is np.nan or variant is np.nan:
        return None
    elif original == variant:
        return None
    elif len(original) != len(variant):
        return None

    y_pred = params['intercept']

    for pos, (orig_base, var_base) in enumerate(zip(original, variant)):
        if orig_base != var_base:
            y_pred += params[f'{pos}']
            y_pred += params[f'{orig_base}{var_base}']

    y_pred += gc_weight * gc_content(original)
    return y_pred


def read_parameters(file_path):
    """Read parameters from a CSV file."""
    try:
        df = pd.read_csv(file_path)
        params = df.set_index('feature')['weight'].to_dict()
        return params
    except Exception as e:
        logging.error(f'Error reading parameters file: {e}')
        sys.exit(1)


def generate_header():
    """Generate and print header row for output."""
    header = ['original', 'variant', 'change_description', 'y_pred']
    print('\t'.join(header))


def find_closest_mismatch(score, mismatches, mismatch_list):
    """Find the closest mismatch based on the desired score."""
    closest_score = None
    closest_mismatch = None
    for mismatch, mismatch_score in mismatches:
        if closest_score is None or abs(mismatch_score - score) < abs(closest_score - score):
            if mismatch not in [x[0] for x in mismatch_list]:
                closest_score = mismatch_score
                closest_mismatch = mismatch
    return closest_mismatch, closest_score


def print_mismatches(mismatch_list, spacer):
    """Print the mismatches in a formatted manner."""
    for i, (mismatch, score) in enumerate(mismatch_list):
        if mismatch is not None:
            target_mismatch = spacer[:mismatch[0]] + mismatch[1] + spacer[mismatch[0]+1:]
            change_description = f"{spacer[mismatch[0]]}{mismatch[0]+1}{mismatch[1]}"
            row = [spacer, target_mismatch, change_description, f"{score:.4f}"]
            print('\t'.join(row))


def generate_mismatches(spacers, min_score, max_score, step, parameters, spacer_original):
    """Generate mismatches for a list of spacers based on desired scores."""
    nucleotides = ['A', 'C', 'G', 'T']

    for spacer in spacers:
        desired_scores = np.arange(min_score, max_score + step, step)

        mismatches = []
        for pos in range(len(spacer)):
            for nt in nucleotides:
                if nt == spacer[pos]:
                    continue
                target_mismatch = spacer[:pos] + nt + spacer[pos+1:]
                y_pred = calculate_y_pred(spacer, target_mismatch, parameters['GC_content'], parameters)
                mismatches.append(((pos, nt), y_pred))

        mismatch_list = []
        for score in desired_scores:
            closest_mismatch, closest_score = find_closest_mismatch(score, mismatches, mismatch_list)
            if closest_mismatch is not None:
                mismatch_list.append((closest_mismatch, closest_score))

    print_mismatches(mismatch_list, spacer_original)  # Use spacer_original instead of spacer


def main(args):
    console=Console(file=sys.stderr)
    console.log("[bold red]Initializing mismatch calculator[/bold red]")


    console.log(f'Reading parameters from {args.parameters_file}...')
    params = read_parameters(args.parameters_file)

    if args.mode == 'mismatches':
        with open(args.spacers_file, 'r') as f:
            spacers_original = [line.strip() for line in f.readlines() if line.strip()]
            spacers = [spacer.upper() for spacer in spacers_original]  # Capitalize the sequences for processing

        generate_header()
        for spacer_original, spacer in zip(spacers_original, spacers):  # Modify this line to iterate over both original and capitalized sequences
            generate_mismatches([spacer], args.min, args.max, args.step, params, spacer_original)  # Pass spacer_original as an additional argument

    elif args.mode == 'recalculate':
        try:
            data = pd.read_csv(args.existing_mismatches, sep='\t')
        except Exception as e:
            console.log(f'[bold red] Error reading input data file: {e}[/bold red]')
            sys.exit(1)

        original_aliases = {'original', 'perfect', 'target'}
        variant_aliases = {'variant', 'mismatch', 'spacer'}

        original_col = original_aliases.intersection(data.columns)
        variant_col = variant_aliases.intersection(data.columns)

        if not (len(original_col) == 1 and len(variant_col) == 1):
            console.log("[bold red] Input data file must have one of [/bold red] 'original', 'target', or 'perfect' [bold red] columns and one of [/bold red] 'variant', 'spacer', or 'mismatch' columns.")
            sys.exit(1)

        original_col = original_col.pop()  # Get the actual column name used
        variant_col = variant_col.pop()  # Get the actual column name used

        data[f'{original_col}_upper'] = data[original_col].str.upper()
        data[f'{variant_col}_upper'] = data[variant_col].str.upper()

        console.log('Calculating predicted mismatch efficacy for each row...')
        new_y_pred_column_name = 'y_pred_new' if 'y_pred' in data.columns else 'y_pred'

        def calculate_and_format(row):
            y_pred_value = calculate_y_pred(row[f'{original_col}_upper'], row[f'{variant_col}_upper'], params['GC_content'], params)
            if y_pred_value is None:
                return None
            return f"{y_pred_value:.4f}"

        data[new_y_pred_column_name] = data.apply(calculate_and_format, axis=1)

        for col in data.columns:
            if data[col].dtype == 'float64':
                # Check if all non-null values in the column are essentially integers
                if data[col].dropna().apply(lambda x: x == int(x)).all():
                    data[col] = data[col].astype('Int64')


        console.log('[bold red]Displaying results[/bold red]')
        print(data.drop([f'{original_col}_upper', f'{variant_col}_upper'], axis=1).to_csv(sep='\t', index=False, na_rep='None'))

        console.log('Done.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate mismatches for a list of spacers and/or recalculate y_pred.')
    parser.add_argument('mode', choices=['mismatches', 'recalculate'], help="Choose the mode of operation: 'mismatches' to generate mismatches, 'recalculate' to recalculate y_pred.")
    parser.add_argument('--spacers_file', help='Path to the file containing original spacers (one per line) (required for mismatches mode).')
    parser.add_argument('--existing_mismatches', help='Path to the input data file (TSV format) (required for recalculate mode).')
    parser.add_argument('--parameters_file', required=True, help='Path to the parameters file (CSV format).')
    parser.add_argument('--verbosity', choices=['debug', 'info', 'warning', 'error', 'critical'], default='info', help='Set the logging verbosity level (default: info).')
    parser.add_argument('--min', type=float, default=0, help='Minimum desired efficacy (default: 0) (required for mismatches mode).')
    parser.add_argument('--max', type=float, default=1, help='Maximum desired efficacy (default: 1) (required for mismatches mode).')
    parser.add_argument('--step', type=float, default=0.1, help='Step between desired efficacies (default: 0.1) (required for mismatches mode).')
    args = parser.parse_args()

    if args.mode == 'mismatches' and args.spacers_file is None:
        parser.error('The --spacers_file option is required for mismatches mode.')
    elif args.mode == 'recalculate' and args.existing_mismatches is None:
        parser.error('The --existing_mismatches option is required for recalculate mode.')

    main(args)
