import csv
from pathlib import Path

def msweep_extract_gsv(file_path):
    num_aligned = 0
    results = []

    with open(file_path, 'r', newline='') as f:
        for line in f:
            if line.startswith('#num_aligned:'):
                num_aligned = int(line.split(':')[1].strip())
                break
        
        if num_aligned < 8:
            return None
        
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 2:
                continue

            c_id = int(parts[0])
            mean_theta = float(parts[1])

            if 1 <= c_id <= 8:
                status = 1 if mean_theta >= 0.01 else 0
                results.append(status)
        
    return results

base_dir = Path('/scratch/user/nhacker/Mh_Adaptive_Sampling/mSWEEP/mSWEEP_output')
out_path=Path('/scratch/user/nhacker/Mh_Adaptive_Sampling/mSWEEP/gsv_table')
out_path.mkdir(parents=True, exist_ok=True)
out_file = out_path / 'gsv_summary.tsv'

def main():
    with open(out_file, 'w', newline='') as fout:
        writer = csv.writer(fout, delimiter='\t')
        header = ['file_name'] + [f'GSV_{i}' for i in range(1,9)]
        writer.writerow(header)

        for file_path in base_dir.glob('*.txt'):
            row_data = msweep_extract_gsv(file_path)

            if row_data is not None:
                file_name = file_path.name
                writer.writerow([file_name] + row_data)
                print(f"[OK] Processed: {file_name}")
            else:
                print(f"[SKIP] Aligned reads < 8: {file_path.name}")

    print(f"Final table written to: {out_file}")

if __name__ == '__main__':
    main()