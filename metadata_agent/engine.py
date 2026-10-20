"""Where the model is allowed to act, and the log that records every call.

A decision point hands the model a *question plus evidence* and demands a JSON object back.
The orchestrator supplies a validator; an answer that fails it is rejected and re-asked with the
failure attached. The model never writes a cell directly - it returns a proposal that code scores.
"""
from __future__ import annotations

import json, os, re, subprocess, time
from dataclasses import dataclass, field
from typing import Callable, Optional

SYSTEM = (
    "You are the judgement component of a metadata-recovery pipeline for public mass-spectrometry "
    "deposits. Deterministic code has already gathered the evidence and will score whatever you "
    "return. Reply with ONE JSON object and nothing else.\n"
    "Rules: prefer abstaining to guessing - an empty cell is cheap, a confidently wrong one poisons "
    "every downstream reanalysis. Never invent a value that the evidence does not support. When two "
    "readings fit equally well, say so in 'unresolved' rather than picking one silently."
)


@dataclass
class Decision:
    kind: str                      # D-GRAM, D-CV, D-CONFLICT, D-ASSAY, D-LINK
    question: str
    context: dict                  # evidence bundle assembled by code
    required_keys: tuple = ('chosen', 'basis')
    options: Optional[list] = None
    attempts: int = 3


class Engine:
    name = 'base'

    def ask(self, d: Decision, validator: Optional[Callable[[dict], tuple]] = None):
        """Returns (answer_dict_or_None, transcript list)."""
        transcript = []
        for i in range(d.attempts):
            raw = self._call(d, transcript)
            if raw is None:
                transcript.append({'attempt': i + 1, 'result': 'no answer from engine'})
                return None, transcript
            try:
                ans = json.loads(raw) if isinstance(raw, str) else raw
            except Exception as e:
                transcript.append({'attempt': i + 1, 'error': 'not JSON: %s' % e, 'raw': str(raw)[:300]})
                continue
            missing = [k for k in d.required_keys if k not in ans]
            if missing:
                transcript.append({'attempt': i + 1, 'error': 'missing keys %s' % missing})
                continue
            if d.options is not None and ans.get('chosen') not in d.options:
                transcript.append({'attempt': i + 1,
                                   'error': 'chosen %r is not one of the offered options' % ans.get('chosen')})
                continue
            if validator:
                ok, why = validator(ans)
                if not ok:
                    transcript.append({'attempt': i + 1, 'error': 'validator rejected: %s' % why})
                    continue
            transcript.append({'attempt': i + 1, 'result': 'accepted'})
            return ans, transcript
        return None, transcript

    def _prompt(self, d: Decision, transcript):
        p = ['QUESTION: ' + d.question, '', 'EVIDENCE:', json.dumps(d.context, indent=1, default=str)[:12000]]
        if d.options is not None:
            p += ['', 'You must choose exactly one of these options for "chosen":',
                  json.dumps(d.options, indent=1)[:4000]]
        p += ['', 'Reply with a JSON object containing: %s' % ', '.join(d.required_keys),
              'Optional keys: "unresolved" (what evidence would settle it), '
              '"options_considered" (list of {option, note}), "confidence_tier" (A/B/C).']
        failures = [t for t in transcript if 'error' in t]
        if failures:
            p += ['', 'Your previous answer was rejected: %s' % failures[-1]['error'],
                  'Return a corrected JSON object.']
        return '\n'.join(p)

    def _call(self, d, transcript):
        raise NotImplementedError


class AbstainEngine(Engine):
    """Default when no model is reachable. Every judgement call becomes an explicit abstention."""
    name = 'abstain'

    def _call(self, d, transcript):
        return None


class AnthropicEngine(Engine):
    name = 'anthropic'

    def __init__(self, model='claude-sonnet-5', max_tokens=1500):
        import anthropic
        self.client = anthropic.Anthropic()
        self.model, self.max_tokens = model, max_tokens

    def _call(self, d, transcript):
        msg = self.client.messages.create(
            model=self.model, max_tokens=self.max_tokens, system=SYSTEM,
            messages=[{'role': 'user', 'content': self._prompt(d, transcript)}])
        text = ''.join(b.text for b in msg.content if getattr(b, 'type', '') == 'text')
        m = re.search(r'\{.*\}', text, re.S)
        return m.group(0) if m else None


class ClaudeCliEngine(Engine):
    """Uses the Claude Code CLI in print mode - no API key needed if the CLI is authenticated."""
    name = 'claude-cli'

    def __init__(self, binary='claude', model=None, timeout=180):
        self.binary, self.model, self.timeout = binary, model, timeout

    def _call(self, d, transcript):
        cmd = [self.binary, '-p', '--output-format', 'json']
        if self.model:
            cmd += ['--model', self.model]
        try:
            r = subprocess.run(cmd, input=SYSTEM + '\n\n' + self._prompt(d, transcript),
                               capture_output=True, text=True, timeout=self.timeout)
        except Exception as e:
            transcript.append({'engine_error': str(e)[:200]})
            return None
        if r.returncode != 0:
            transcript.append({'engine_error': r.stderr[:300]})
            return None
        try:
            payload = json.loads(r.stdout)
            text = payload.get('result', r.stdout)
        except Exception:
            text = r.stdout
        m = re.search(r'\{.*\}', text, re.S)
        return m.group(0) if m else None


class FileEngine(Engine):
    """Hand-off mode: queue questions to disk, read answers back on the next run.

    Lets a human - or a Claude session with no API wiring - answer the judgement calls while the
    deterministic pipeline keeps running unattended.
    """
    name = 'file'

    def __init__(self, workdir):
        self.q = os.path.join(workdir, 'pending_decisions.jsonl')
        self.a = os.path.join(workdir, 'answers.jsonl')
        os.makedirs(workdir, exist_ok=True)
        self.answers = {}
        if os.path.exists(self.a):
            for line in open(self.a, encoding='utf-8'):
                line = line.strip()
                if line:
                    r = json.loads(line)
                    self.answers[r.get('id')] = r.get('answer')

    def _call(self, d, transcript):
        key = '%s:%s' % (d.kind, d.question[:80])
        if key in self.answers:
            return self.answers[key]
        if key in self._queued():
            return None                      # already asked; do not duplicate the queue
        with open(self.q, 'a', encoding='utf-8', newline='\n') as f:
            f.write(json.dumps({'id': key, 'kind': d.kind, 'question': d.question,
                                'options': d.options, 'required_keys': list(d.required_keys),
                                'context': d.context}, default=str) + '\n')
        return None

    def _queued(self):
        if not os.path.exists(self.q):
            return set()
        out = set()
        for line in open(self.q, encoding='utf-8'):
            line = line.strip()
            if line:
                try:
                    out.add(json.loads(line).get('id'))
                except ValueError:
                    pass
        return out


def build_engine(spec, workdir):
    spec = (spec or 'auto').lower()
    if spec == 'auto':
        if os.environ.get('ANTHROPIC_API_KEY'):
            spec = 'anthropic'
        else:
            from shutil import which
            spec = 'claude-cli' if which('claude') else 'file'
    if spec == 'anthropic':
        return AnthropicEngine()
    if spec == 'claude-cli':
        return ClaudeCliEngine()
    if spec == 'file':
        return FileEngine(workdir)
    return AbstainEngine()


# ----------------------------------------------------------------- the log ----


class DecisionLog:
    """Append-only, written at the moment each call is made."""

    def __init__(self, path):
        self.path = path
        self.n = 0
        self._ids = {}
        open(path, 'w', encoding='utf-8').close()

    def _next(self, prefix):
        self._ids[prefix] = self._ids.get(prefix, 0) + 1
        return '%s-%02d' % (prefix, self._ids[prefix])

    def record(self, prefix, stage, question, chosen, basis, evidence=(), tier='A',
               affects=(), options=None, unresolved=None, decided_by='code', transcript=None):
        rec = {'id': self._next(prefix), 'stage': stage, 'question': question, 'chosen': chosen,
               'basis': basis, 'evidence': list(evidence), 'confidence_tier': tier,
               'affects': list(affects), 'decided_by': decided_by}
        if options:
            rec['options_considered'] = options
        if unresolved:
            rec['unresolved'] = unresolved
        if transcript:
            rec['engine_transcript'] = transcript
        with open(self.path, 'a', encoding='utf-8', newline='\n') as f:
            f.write(json.dumps(rec, default=str) + '\n')
        self.n += 1
        return rec['id']

    def abstain(self, column, reason, basis, stage=7, evidence=()):
        return self.record('D-ABSTAIN', stage, 'Why is %s empty?' % column,
                           'empty (%s)' % reason, basis, evidence, 'A', [column])
