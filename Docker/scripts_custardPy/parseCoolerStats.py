#! /usr/bin/env python
# -*- coding: utf-8 -*-
import re
import argparse

def parse_line(line):
    match = re.match(r'(.+?)\s{2,}(\S+)\s+(\S+)%', line)
    if match:
        return match.groups()
    else:
        return None

def parse_file(file_path):
    headers = []
    numbers = []
    percentages = []
    with open(file_path, 'r') as file:
        for line in file:
            parsed = parse_line(line.strip())
            if parsed:
                header, number, percentage = parsed
                headers.append(header)
                numbers.append(number)
                percentages.append(f"{percentage}%")

    return headers, numbers, percentages


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="stats file (*/qc_report/mapping_stats.txt)", type=str)
    parser.add_argument("--header", help="output header", action='store_true')

    args = parser.parse_args()
    file_path = args.input
    headers, numbers, percentages = parse_file(file_path)

    alternating = []
    for number, percentage in zip(numbers, percentages):
        alternating.append('\t'.join([number, percentage]))

    if args.header:
        print('\t%\t'.join(headers), sep='\t')
    else:
        print('\t'.join(alternating), sep='\t')
