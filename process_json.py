import json
import argparse

def merge_json(json1, json2):
    if isinstance(json1, dict) and isinstance(json2, dict):
        merged = dict(json1)  
        for key, value in json2.items():
            if key in merged:
                if merged[key] == value:
                    continue
                merged[key] = merge_json(merged[key], value)
            else:
                merged[key] = value
        return merged

    elif isinstance(json1, list) and isinstance(json2, list):
        seen = set()
        unique = []
        for item in json1 + json2:
            item_str = json.dumps(item, sort_keys=True)
            if item_str not in seen:
                seen.add(item_str)
                unique.append(item)
        return unique
    
    elif json1 == json2:
        print("Both input files are identical, nothing to process!")
        return None 

    else:
        print("Warning: Merging could not be completed. Check your inputs.")
        return None  


if __name__ == "__main__":
    # create a parser for command line arguments 
    parser = argparse.ArgumentParser(
        description='Process WG and WT json files to create merged input for Analysis Report'
    )
    parser.add_argument(
        'wgjson',
        type=str,
        help='Path to the WGS json input'
    )
    parser.add_argument(
        'wtjson',
        type=str,
        help='Path to the WT json input'
    )
    parser.add_argument(
        'study',
        type=str,
        help='Study Title'
    )
    args = parser.parse_args()

# Load the JSON data
with open(args.wgjson, 'r') as f1, open(args.wtjson, 'r') as f2:
    json1 = json.load(f1)
    json2 = json.load(f2)

merged_json = merge_json(json1, json2)

# Save the result
with open(f'{args.study}.json', 'w') as f_out:
    json.dump(merged_json, f_out, indent=2)

