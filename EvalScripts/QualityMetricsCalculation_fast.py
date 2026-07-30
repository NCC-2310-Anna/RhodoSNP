# -*- coding: utf-8 -*-
"""
Optimized version of QualityMetricsCalculation.py
Faster, cleaner, and more memory-efficient
"""

from argparse import ArgumentParser


def read_file(filepath):
    """Generator: reads file line by line (memory efficient)"""
    with open(filepath, "r") as f:
        for line in f:
            yield line.rstrip("\n")


def extract_truth_dict(lines, key_index):
    """Create dictionary from truth data"""
    data = {}
    for line in lines:
        if line.startswith("Chr"):
            continue
        parts = line.split("\t")
        key = parts[key_index]
        data[key] = parts[:key_index] + parts[key_index + 1:]
    return data


def extract_snp_dict(lines, key_index):
    """Create dictionary from SNP data"""
    data = {}
    for line in lines:
        if line.startswith("POS"):
            continue
        parts = line.split("\t")
        key = parts[key_index]
        data[key] = parts[:key_index] + parts[key_index + 1:]
    return data


def calculate_quality_metrics(truth_data, snp_data):
    true_positive = 0
    false_positive = 0
    false_negative = 0

    truth_keys = set(truth_data.keys())
    snp_keys = set(snp_data.keys())

    # True positives & false positives
    for pos, values in snp_data.items():
        if pos in truth_data:
            truth_values = truth_data[pos]
            if (
                values[0] == truth_values[0]
                and values[3] == "SNP"
                and truth_values[4] == "*"
            ):
                true_positive += 1
        else:
            if values[3] == "SNP":
                false_positive += 1

    # False negatives
    for pos in truth_keys:
        if pos not in snp_keys:
            false_negative += 1

    # Avoid division by zero
    precision = (
        true_positive / (true_positive + false_positive)
        if (true_positive + false_positive) > 0
        else 0
    )
    recall = (
        true_positive / (true_positive + false_negative)
        if (true_positive + false_negative) > 0
        else 0
    )
    f1_score = (
        (2 * precision * recall) / (precision + recall)
        if (precision + recall) > 0
        else 0
    )

    return [
        true_positive,
        false_positive,
        false_negative,
        precision,
        recall,
        f1_score,
    ]


def parse_metadata(filepath):
    filename = filepath.split("/")[-1]
    name = filename.replace(".tsv", "")
    caller = name.split("_")[-1]
    return [filename, caller]


def main():
    parser = ArgumentParser(description="Calculate evaluation metrics")
    parser.add_argument("-t", "--truthData", required=True)
    parser.add_argument("-f", "--snpData", required=True)
    parser.add_argument("-o", "--organism", required=True)
    args = parser.parse_args()

    truth_data = extract_truth_dict(read_file(args.truthData), key_index=1)
    snp_data = extract_snp_dict(read_file(args.snpData), key_index=0)

    metadata = parse_metadata(args.snpData)
    metrics = calculate_quality_metrics(truth_data, snp_data)

    print(args.organism,*metadata, *metrics, sep="\t")


if __name__ == "__main__":
    main()
