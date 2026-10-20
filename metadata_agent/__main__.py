"""CLI: python -m metadata_agent MSV000102356 --out runs/MSV000102356 --engine auto"""
import argparse, json, sys

from .core import Config
from .pipeline import Run


def main(argv=None):
    p = argparse.ArgumentParser(prog='metadata_agent',
                                description='Per-deposit metadata recovery: code calls the tools, '
                                            'the model only makes judgement calls.')
    p.add_argument('accession')
    p.add_argument('--out', default=None, help='output directory (default runs/<accession>)')
    p.add_argument('--engine', default='auto',
                   choices=['auto', 'anthropic', 'claude-cli', 'file', 'abstain'],
                   help='auto = anthropic if ANTHROPIC_API_KEY, else claude CLI, else file hand-off')
    p.add_argument('--headers-per-subtree', type=int, default=2)
    p.add_argument('--grammar-coverage', type=float, default=None)
    p.add_argument('--join-match-rate', type=float, default=None)
    p.add_argument('--summary', default=None,
                   help='shared one-row-per-dataset CSV to upsert into (fleet mode)')
    p.add_argument('--extras', action='store_true',
                   help='also write the per-file derived CSV')
    p.add_argument('--keep-log', action='store_true',
                   help='also keep the decision log JSONL (default: only the two CSVs)')
    p.add_argument('--candidates', default=None,
                   help='JSON file of publication candidates for stage 5')
    a = p.parse_args(argv)

    cfg = Config()
    if a.grammar_coverage is not None:
        cfg.grammar_coverage = a.grammar_coverage
    if a.join_match_rate is not None:
        cfg.join_match_rate = a.join_match_rate
    out = a.out or ('runs/%s' % a.accession)
    cands = json.load(open(a.candidates, encoding='utf-8')) if a.candidates else None

    run = Run(a.accession, out, cfg=cfg, engine=a.engine, summary_path=a.summary,
              extras=a.extras, keep_log=a.keep_log)
    result = run.run(header_k=a.headers_per_subtree, candidates=cands)
    print(json.dumps({'accession': a.accession, 'engine': run.engine.name, 'out': out, **result},
                     indent=1, default=str))
    return 0


if __name__ == '__main__':
    sys.exit(main())
