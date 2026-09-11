"""Reproduce the Gaussian checks and paired E3-prime SMART experiments."""
from __future__ import annotations
import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent

def write_json(path: Path, obj: object) -> None:
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(obj, indent=2), encoding='utf-8')
    temp.replace(path)

def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def execute(job: tuple[list[str], Path, dict[str, str]]) -> None:
    command, log, env = job
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open('w', encoding='utf-8') as stream:
        subprocess.run(command, cwd=ROOT, env=env, stdout=stream,
                       stderr=subprocess.STDOUT, check=True)
    print(f'Completed {log.parent.name}', flush=True)

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('profile', choices=('smoke', 'paper'))
    parser.add_argument('--criterion', choices=('logrss', 'wald', 'both'), default='both')
    parser.add_argument('--replicates', type=int)
    parser.add_argument('--sizes')
    parser.add_argument('--seed', type=int, default=20260909)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--skip-gaussian', action='store_true')
    args = parser.parse_args()
    reps = args.replicates if args.replicates is not None else (60 if args.profile == 'smoke' else 5000)
    sizes = [int(x) for x in (args.sizes or ('100' if args.profile == 'smoke' else '250,1000,4000')).split(',')]
    if reps < 2 or min(sizes) < 12 or len(set(sizes)) != len(sizes) or args.workers < 1:
        parser.error('Use at least two replicates, distinct sample sizes >=12 and positive workers.')
    dest = (args.output or ROOT / 'results' / args.profile).resolve()
    dest.mkdir(parents=True, exist_ok=True)
    variants = ('logrss', 'wald') if args.criterion == 'both' else (args.criterion,)
    source_names = ('feedback_tuning.py', 'verify_smart.py', 'verify_gaussian.py')
    manifest = {'profile': args.profile, 'replicates': reps, 'sizes': sizes,
                'seed': args.seed, 'criterion_variants': list(variants),
                'source_sha256': {n:file_hash(ROOT/n) for n in source_names},
                'temperature': 2.0, 'prior': 'uniform',
                'library': 'nine pairs in {0.5,1,2}^2 and the all-wide reference',
                'distinct_training_datasets': 4 * len(sizes) * reps,
                'note': 'Criterion variants share the same configuration-specific streams.'}
    saved = dest/'run_manifest.json'
    if saved.exists() and json.loads(saved.read_text()) != manifest:
        raise SystemExit('Settings or numerical source files differ from the saved run. Use a new --output directory.')
    write_json(saved,manifest)
    env = os.environ.copy()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
        env[key]='1'
    jobs=[]
    for variant in variants:
        for case in range(4):
            for n in sizes:
                folder=dest/variant/f'case{case}_n{n}'
                done=folder/'smart_results.json'
                if done.exists():
                    existing=json.loads(done.read_text())
                    if len(existing)==1 and existing[0]['replicates']==reps and (folder/'smart_results.csv').exists() and (folder/'smart_metadata.json').exists():
                        print(f'Reusing {folder.name} ({variant})',flush=True)
                        continue
                cmd=[sys.executable,str(ROOT/'verify_smart.py'),'--replicates',str(reps),
                     '--sizes',str(n),'--seed',str(args.seed),'--criterion',variant,
                     '--case',str(case),'--output',str(folder)]
                jobs.append((cmd,folder/'run.log',env))
    start=time.monotonic()
    if jobs:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
            for task in executor.map(execute,jobs):
                pass
    for variant in variants:
        combined=[]
        for case in range(4):
            for n in sizes:
                combined.extend(json.loads((dest/variant/f'case{case}_n{n}'/'smart_results.json').read_text()))
        for row in combined:
            row['risk_score_valid_replicates']=reps-row['stabilization_count']
            row['risk_score_conditional_on_success']=bool(row['stabilization_count'])
        write_json(dest/variant/'smart_results.json',combined)
        keys=[k for k in combined[0] if k not in ('mean_alpha','score_covariance')]
        with (dest/variant/'smart_results.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=keys);writer.writeheader()
            writer.writerows({k:r[k] for k in keys} for r in combined)
        paper_keys=['delta1','delta2','n','replicates','akaike_minus_fa','paired_mcse','gaussian_limit','stabilization_count']
        with (dest/variant/'main_table_gaps.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=paper_keys);writer.writeheader()
            for r in combined:
                writer.writerow(dict(delta1=r['delta1'],delta2=r['delta2'],n=r['n'],replicates=reps,
                    akaike_minus_fa=-r['fa_minus_akaike'],paired_mcse=r['fa_minus_akaike_mcse'],
                    gaussian_limit=-r['limit_gap'],stabilization_count=r['stabilization_count']))
    if not args.skip_gaussian:
        order=128 if args.profile=='smoke' else 1024
        folder=dest/'gaussian'
        if not (folder/'gaussian_results.csv').exists():
            execute(([sys.executable,str(ROOT/'verify_gaussian.py'),'--order',str(order),
                      '--output',str(folder)],folder/'run.log',env))
    import numpy, scipy
    write_json(dest/'environment.json',{'python':sys.version,'platform':platform.platform(),
        'numpy':numpy.__version__,'scipy':scipy.__version__,
        'workers':args.workers,'elapsed_seconds_this_invocation':time.monotonic()-start})
    print(f'Results: {dest}',flush=True)

if __name__=='__main__':
    main()
