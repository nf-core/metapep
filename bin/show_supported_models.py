#!/usr/bin/env python
# Extract supported alleles and lengths from supported_alleles.json

import argparse
import json
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract supported models information from supported_alleles.json"
    )
    parser.add_argument(
        "-j",
        "--json",
        help="Path to supported_alleles.json",
        required=True,
    )
    
    return parser.parse_args()


# Tool-specific configurations
supported = {
    "mhcflurry": {
        "version": "2.1.4",
        "lengths": list(range(5, 16)) # 5-15
    },
    "mhcnuggets": {
        "version": "2.4.1",
        "lengths": list(range(5, 16)) # 5-15 
    },
    "mhcnuggetsii": {
        "version": "2.4.1",
        "lengths": list(range(5, 31))  # 9-30
    },
    "netmhcpan": {
        "lengths": list(range(8, 15))  # 8-14
    },
    "netmhciipan": {
        "lengths": list(range(9, 51))  # 9-50
    }
}

def load_json(filepath):
    """Load JSON file"""
    with open(filepath, 'r') as f:
        return json.load(f)

def main():
    args = parse_args() 
    # Load supported_alleles.json
    supported_alleles = load_json(args.json)

    # Process each tool found in the JSON
    for method, alleles in supported_alleles.items():
        if method not in supported:
            print(f" Skipping unknown method: {method}")
            continue
        
        config = supported[method]
        version = config.get("version")
        lengths = config["lengths"]
        
        # Write supported alleles
        if version:
            alleles_file = f"{method}.v{version}.supported_alleles.txt"
            lengths_file = f"{method}.v{version}.supported_lengths.txt"
        else:
            alleles_file = f"{method}.supported_alleles.txt"
            lengths_file = f"{method}.supported_lengths.txt"
        
        with open(alleles_file, "w") as output:
            for allele in sorted(alleles):
                output.write(allele + "\n")
        
        print(f" Created {alleles_file} ({len(alleles)} alleles)")
        
        # Write supported lengths
        with open(lengths_file, "w") as output:
            for length in lengths:
                output.write(str(length) + "\n")
        
        print(f" Created {lengths_file} ({len(lengths)} lengths)") 
    print("All model files created successfully!")

if __name__ == "__main__":
    sys.exit(main())