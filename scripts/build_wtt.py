"""
Build script: converts data/wtt_stops.csv → data/wtt_data.json

Run once from the railinsights-site folder:
    python scripts/build_wtt.py

Output is a compact indexed JSON (~30-50 MB) that the page fetches directly.
"""
import csv
import json
import sys
from pathlib import Path

# Use the full multi-book WTT if available, fall back to the local copy
_FULL = Path(r'C:\Users\LOQ\Downloads\01. May 2026 timetable - separate PDFs\wtt_stops.csv')
SRC   = _FULL if _FULL.exists() else Path('data/wtt_stops.csv')
OUT   = Path('data/wtt_data.json')

def main():
    if not SRC.exists():
        print(f'ERROR: {SRC} not found. Run from the railinsights-site folder.')
        sys.exit(1)

    print(f'Reading {SRC} ({SRC.stat().st_size / 1e6:.0f} MB)…')

    locs_index = {}
    locs_list  = []
    services   = {}

    def loc_id(name):
        if not name:
            return None
        if name not in locs_index:
            locs_index[name] = len(locs_list)
            locs_list.append(name)
        return locs_index[name]

    with open(SRC, encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i % 200_000 == 0:
                print(f'  {i:,} rows…', end='\r')

            uid = row.get('uid', '').strip()
            if not uid:
                continue

            li = loc_id(row.get('location', '').strip())

            if uid not in services:
                services[uid] = {
                    'h':  row.get('headcode', '').strip(),
                    'u':  uid,
                    'op': row.get('operator', '').strip(),
                    'to': row.get('toc', '').strip(),
                    'or': loc_id(row.get('origin', '').strip()),
                    'ot': row.get('origin_time', '').strip(),
                    'ds': loc_id(row.get('destination', '').strip()),
                    'dt': row.get('destination_time', '').strip(),
                    'tl': row.get('timing_load', '').strip(),
                    'do': row.get('dates_of_operation', '').strip(),
                    'rd': row.get('running_days', '').strip(),
                    'sc': row.get('service_code', '').strip(),
                    'bk': set(),
                    's':  [],
                }
            book = row.get('wtt_book', '').strip()
            if book:
                services[uid]['bk'].add(book)

            # Stop tuple: [loc_idx, platform, arr, dep, pass_time]
            # Omit trailing empty strings to save space
            stop = [
                li,
                row.get('platform', '').strip(),
                row.get('arr', '').strip(),
                row.get('dep', '').strip(),
                row.get('pass_time', '').strip(),
            ]
            # Trim trailing empty strings
            while len(stop) > 1 and stop[-1] == '':
                stop.pop()
            services[uid]['s'].append(stop)

    print(f'\nParsed {len(services):,} services, {len(locs_list):,} locations.')

    # Convert book sets to sorted lists for JSON serialisation
    for svc in services.values():
        svc['bk'] = sorted(svc['bk'])

    data = {'locs': locs_list, 'svcs': list(services.values())}

    print(f'Writing {OUT}…')
    with open(OUT, 'w', encoding='utf-8') as f:
        json.dump(data, f, separators=(',', ':'), ensure_ascii=False)

    size_mb = OUT.stat().st_size / 1e6
    print(f'Done. {OUT}  —  {size_mb:.1f} MB')
    if size_mb > 100:
        print('WARNING: file is over 100 MB — GitHub will block this file.')
        print('Consider using Git LFS:  git lfs track "data/wtt_data.json"')
    else:
        print('File is within GitHub\'s 100 MB limit. You can commit it normally.')

if __name__ == '__main__':
    main()
